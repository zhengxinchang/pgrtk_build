const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use iset::set::IntervalSet;
// use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Generate diploid VCF field from paired alnmap file from two haplotype assembly
#[derive(Parser, Debug)]
#[clap(name = "pgr-generate-diploid-vcf")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the first haplotype alnmap file
    hap0_path: String,
    /// path to the second haplotype alnmap file
    hap1_path: String,
    /// path to a ctgmap.json file
    target_len_json_path: String,
    /// the prefix of the output files
    output_prefix: String,
    /// the prefix of the output files
    #[clap(long, default_value = "Sample")]
    sample_name: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

type TargetSeqLength = Vec<(u32, String, u32)>;

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32);
type VariantRecord = (String, u32, u32, u8, String, String, String); //t_name, tc, tl, hap_type, tvs, qvs, rec_type

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut target_length_json_file = BufReader::new(
        File::open(Path::new(&args.target_len_json_path)).expect("can't open the input file"),
    );
    let mut buffer = Vec::new();
    target_length_json_file.read_to_end(&mut buffer)?;
    let mut target_length: TargetSeqLength =
        serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
            .expect("can't parse the target_len.json file");

    target_length.sort();

    let hap0_alnmap_file = BufReader::new(File::open(Path::new(&args.hap0_path)).unwrap());

    let hap1_alnmap_file = BufReader::new(File::open(Path::new(&args.hap1_path)).unwrap());

    #[allow(clippy::type_complexity)]
    let get_variant_recs = |f: BufReader<File>,
                            hap_type: u8|
     -> (
        Vec<(String, u32, u32, u8, String, String, String)>,
        FxHashMap<u32, Vec<ShimmerMatchBlock>>,
        FxHashMap<u32, Vec<ShimmerMatchBlock>>,
    ) {
        let mut variant_records = Vec::<(String, u32, u32, u8, String, String, String)>::new();
        let mut aln_blocks = FxHashMap::<u32, Vec<ShimmerMatchBlock>>::default();
        let mut unique_aln_blocks = FxHashMap::<u32, Vec<ShimmerMatchBlock>>::default();

        f.lines().for_each(|line| {
            if let Ok(line) = line {
                let fields = line.split('\t').collect::<Vec<&str>>();
                assert!(fields.len() > 3);
                let rec_type = fields[1];
                if rec_type.starts_with('V') {
                    assert!(fields.len() == 15);
                    let err_msg = format!("fail to parse on {}", line);
                    let t_name = fields[2];
                    // let ts = fields[3].parse::<u32>().expect(&err_msg);
                    // let te = fields[4].parse::<u32>().expect(&err_msg);
                    // let q_name = fields[5];
                    // let qs = fields[6].parse::<u32>().expect(&err_msg);
                    // let qe = fields[7].parse::<u32>().expect(&err_msg);
                    // let orientation = fields[8].parse::<u32>().expect(&err_msg);
                    // let td = fields[9].parse::<u32>().expect(&err_msg);
                    // let qd = fields[10].parse::<u32>().expect(&err_msg);
                    let tc = fields[11].parse::<u32>().expect(&err_msg);
                    // let tt = fields[12].chars().next().expect(&err_msg);
                    let tvs = fields[13];
                    let qvs = fields[14];
                    variant_records.push((
                        t_name.to_string(),
                        tc,
                        tvs.len() as u32,
                        hap_type,
                        tvs.to_string(),
                        qvs.to_string(),
                        rec_type.to_string(),
                    ));
                };

                if rec_type.starts_with('M') || rec_type.starts_with('V') {
                    let err_msg = format!("fail to parse on {}", line);
                    let aln_block_id = fields[0].parse::<u32>().expect(&err_msg);
                    let t_name = fields[2];
                    let ts = fields[3].parse::<u32>().expect(&err_msg);
                    let te = fields[4].parse::<u32>().expect(&err_msg);
                    let q_name = fields[5];
                    let qs = fields[6].parse::<u32>().expect(&err_msg);
                    let qe = fields[7].parse::<u32>().expect(&err_msg);
                    let orientation = fields[8].parse::<u32>().expect(&err_msg);
                    let e = aln_blocks.entry(aln_block_id).or_default();
                    e.push((
                        t_name.to_string(),
                        ts,
                        te,
                        q_name.to_string(),
                        qs,
                        qe,
                        orientation,
                    ));
                    if rec_type == "M" || rec_type == "V" {
                        let e = unique_aln_blocks.entry(aln_block_id).or_default();
                        e.push((
                            t_name.to_string(),
                            ts,
                            te,
                            q_name.to_string(),
                            qs,
                            qe,
                            orientation,
                        ));
                    }
                }
            }
        });
        (variant_records, aln_blocks, unique_aln_blocks)
    };
    let (hap0_recs, hap0_aln_blocks, hap0_unique_aln_blocks) = get_variant_recs(hap0_alnmap_file, 0);
    let (hap1_recs, hap1_aln_blocks, hap1_unique_aln_blocks) = get_variant_recs(hap1_alnmap_file, 1);

    let blocks_to_intervals = |blocks: FxHashMap<u32, Vec<ShimmerMatchBlock>>| -> FxHashMap<String, IntervalSet<u32>> {
        let mut aln_intervals = FxHashMap::<String, IntervalSet<u32>>::default();
        blocks.into_iter().for_each(|(_block_id, records)| {
            records.into_iter().for_each(|rec| {
                let t_name = rec.0;
                let bgn = rec.1;
                let end = rec.2;
                let interval_set = aln_intervals.entry(t_name).or_default();
                interval_set.insert(bgn..end);
            })
        });
        aln_intervals
    };


    let hap0_aln_intervals = blocks_to_intervals(hap0_aln_blocks);
    let hap1_aln_intervals = blocks_to_intervals(hap1_aln_blocks);
    let hap0_unique_aln_intervals = blocks_to_intervals(hap0_unique_aln_blocks);
    let hap1_unique_aln_intervals = blocks_to_intervals(hap1_unique_aln_blocks);
    

    let mut out_vcf =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("vcf")).unwrap());
    let mut out_bed =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("bed")).unwrap());
    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    target_length.into_iter().for_each(|(_, t_name, t_len)| {
        writeln!(out_vcf, r#"##contig=<ID={},length={}>"#, t_name, t_len)
            .expect("fail to write the vcf file");
    });
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=DUP,Description="duplicated alignment block">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=OVLP,Description="overlapped alignment block">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )
    .expect("fail to write the vcf file");

    writeln!(
        out_vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        args.sample_name
    )
    .expect("fail to write the vcf file");

    let convert_to_vcf_record = |records: &Vec<VariantRecord>| {
        let mut ref_bases = FxHashSet::<(u32, char)>::default();
        let mut h0alleles = Vec::<(u32, VariantRecord)>::new();
        let mut h1alleles = Vec::<(u32, VariantRecord)>::new();
        let mut al_idx_map = FxHashMap::<(u32, String, String), u32>::default();
        let mut al_idx = 0_u32;

        let ref_name = records.first().unwrap().0.clone();
        let mut rec_type = Option::<String>::None;
        records.iter().for_each(|rec| {
            let (_t_name, ts, tl, ht, vts, vqs, rt) = rec;

            if rec_type.is_none() && (rt == "V_D" || rt == "V_O") {
                rec_type = Some(rt.clone());
            }

            (0..*tl).for_each(|t_pos| {
                let vts = vts.chars().collect::<Vec<_>>();
                ref_bases.insert((*ts + t_pos, vts[t_pos as usize]));
            });

            let key = (*ts, vts.clone(), vqs.clone());

            al_idx_map.entry(key).or_insert_with(|| {
                al_idx += 1;
                al_idx
            });

            if *ht == 0 {
                h0alleles.push((al_idx, rec.clone()));
            };
            if *ht == 1 {
                h1alleles.push((al_idx, rec.clone()));
            };
        });
        let mut ref_bases = ref_bases.into_iter().collect::<Vec<_>>();
        ref_bases.sort();
        let ref_str = String::from_iter(ref_bases.iter().map(|(_, c)| *c).collect::<Vec<_>>());
        assert!(ref_str.len() == ref_bases.len()); // make sure all bases at the same t_pos are the same, if not the vectors will have different lengths
        let ts0 = ref_bases.first().unwrap().0;
        let tl0 = ref_str.len() as u32;

        let mut query_alleles = al_idx_map
            .iter()
            .map(|((ts, tvs, qvs), &al_idx)| {
                let t_prefix = ref_str[0..(ts - ts0) as usize].to_string();
                let t_suffix = ref_str[(ts + tvs.len() as u32 - ts0) as usize..].to_string();
                (al_idx, [t_prefix, qvs.clone(), t_suffix].join(""))
            })
            .collect::<Vec<_>>();

        query_alleles.sort();
        let query_alleles = query_alleles
            .iter()
            .map(|(_, qs)| qs.clone())
            .collect::<Vec<_>>()
            .join(",");

        let h0_al_idx = if let Some(i_set) = hap0_aln_intervals.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                if h0alleles.is_empty() {
                    "0".to_string()
                } else {
                    format!("{}", h0alleles.last().unwrap().0)
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let h1_al_idx = if let Some(i_set) = hap1_aln_intervals.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                if h1alleles.is_empty() {
                    "0".to_string()
                } else {
                    format!("{}", h1alleles.last().unwrap().0)
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let gt = [h0_al_idx, h1_al_idx].join("|");
        (ref_name, ts0, ref_str, query_alleles, gt, rec_type)
    };

    let mut variant_records = Vec::<VariantRecord>::new();
    variant_records.extend(hap0_recs);
    variant_records.extend(hap1_recs);

    let mut variant_records = variant_records.into_iter().collect::<Vec<_>>();
    // variant_group: represent a group of overlapped variants, ref_id, ref_start, len, REF, ALT
    let mut variant_group = Vec::<VariantRecord>::new();
    // currrent_vg_end: represent the end coordinate of the current variant group
    let mut currrent_vg_end = Option::<(String, u32)>::None;
    variant_records.sort();
    variant_records
        .into_iter()
        .for_each(|(ref_name, ts, tl, ht, vts, vqs, rec_type)| {
            if let Some(currrent_vg_end) = currrent_vg_end.clone() {
                //println!("{} {} {} {} {:?} {}", ref_name, ts, tl, ts + tl ,  currrent_vg_end, variant_group.len()  );
                if ref_name == currrent_vg_end.0 && ts < currrent_vg_end.1 {
                    variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs, rec_type));
                } else if !variant_group.is_empty() {
                    // println!("X {} {} {} {}", ref_name, ts, tl, variant_group.len());
                    let (vcf_rec_ref_name, ts0, ref_str, query_alleles, gt, g_rec_type) =
                        convert_to_vcf_record(&variant_group);
                    let rt = if let Some(g_rec_type) = g_rec_type {
                        if g_rec_type == "V_D" {
                            "DUP"
                        } else if g_rec_type == "V_O" {
                            "OVLP"
                        } else {
                            "PASS"
                        }
                    } else {
                        "PASS"
                    };
                    let qv: u32 = if rt != "PASS" { 30 } else { 40 };
                    writeln!(
                        out_vcf,
                        "{}\t{}\t.\t{}\t{}\t{}\t{}\t.\tGT\t{}",
                        vcf_rec_ref_name,
                        ts0 + 1,
                        ref_str,
                        query_alleles,
                        qv,
                        rt,
                        gt,
                    )
                    .expect("fail to write the vcf file");
                    variant_group.clear();
                    variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs, rec_type));
                }
            } else {
                variant_group.push((ref_name.clone(), ts, tl, ht, vts, vqs, rec_type));
                currrent_vg_end = Some((ref_name, ts + tl));
                return;
            }
            currrent_vg_end = Some((ref_name.clone(), ts + tl));
        });
    if !variant_group.is_empty() {
        // println!("X {} {} {} {}", ref_name, ts, tl, variant_group.len());
        let (vcf_rec_ref_name, ts0, ref_str, query_alleles, gt, g_rec_type) =
            convert_to_vcf_record(&variant_group);
        let rt = if let Some(g_rec_type) = g_rec_type {
            if g_rec_type == "V_D" {
                "DUP"
            } else if g_rec_type == "V_O" {
                "OVLP"
            } else {
                "PASS"
            }
        } else {
            "PASS"
        };
        let qv: u32 = if rt != "PASS" { 30 } else { 40 };
        writeln!(
            out_vcf,
            "{}\t{}\t.\t{}\t{}\t{}\t{}\t.\tGT\t{}",
            vcf_rec_ref_name,
            ts0 + 1,
            ref_str,
            query_alleles,
            qv,
            rt,
            gt,
        )
        .expect("fail to write the vcf file");
    };

    let merge_intervals =
        |intervals: FxHashMap<String, IntervalSet<u32>>| -> FxHashMap<String, IntervalSet<u32>> {
            intervals
                .into_iter()
                .flat_map(|(t_name, i_set)| {
                    let mut intervals = i_set
                        .unsorted_iter()
                        .map(|range| (range.start, range.end))
                        .collect::<Vec<_>>();
                    if intervals.is_empty() {
                        return None;
                    }
                    intervals.sort();
                    let mut merged_intervals = IntervalSet::<u32>::new();
                    let mut current_range = intervals.first().unwrap().to_owned();
                    intervals.into_iter().for_each(|(bgn, end)| {
                        if bgn <= current_range.1 && end > current_range.1 {
                            current_range.1 = end;
                        } else if bgn > current_range.1 {
                            merged_intervals.insert(current_range.0..current_range.1);
                            current_range = (bgn, end);
                        }
                    });
                    merged_intervals.insert(current_range.0..current_range.1);
                    Some((t_name, merged_intervals))
                })
                .collect()
        };

    let hap0_aln_merged_intervals: FxHashMap<String, IntervalSet<u32>> =
        merge_intervals(hap0_unique_aln_intervals);

    let hap1_aln_merged_intervals: FxHashMap<String, IntervalSet<u32>> =
        merge_intervals(hap1_unique_aln_intervals);

    let mut t_names = hap0_aln_merged_intervals
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    t_names.sort();

    t_names.into_iter().for_each(|t_name| {
        let hap0_aln_merged_intervals = hap0_aln_merged_intervals.get(&t_name).unwrap();
        if let Some(hap1_aln_merged_intervals) = hap1_aln_merged_intervals.get(&t_name) {
            let mut intervals = hap0_aln_merged_intervals
                .unsorted_iter()
                .map(|range| (range.start, range.end))
                .collect::<Vec<_>>();
            if intervals.is_empty() {
                return;
            }
            intervals.sort();
            let msg = "can't write the output bed file";
            intervals.into_iter().for_each(|(bgn, end)| {
                hap1_aln_merged_intervals.iter(bgn..end).for_each(|range| {
                    let (bgn1, end1) = (range.start, range.end);
                    if bgn1 < bgn && end1 < end {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn, end1).expect(msg);
                    } else if bgn1 < bgn && end <= end1 {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn, end).expect(msg);
                    } else if bgn <= bgn1 && end1 < end {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn1, end1).expect(msg);
                    } else if bgn <= bgn1 && end <= end1 {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn1, end).expect(msg);
                    }
                });
            });
        }
    });

    Ok(())
}
