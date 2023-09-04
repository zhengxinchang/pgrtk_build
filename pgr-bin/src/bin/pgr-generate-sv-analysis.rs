const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use iset::set::IntervalSet;
// use rayon::prelude::*;
use pgr_db::ext::{
    get_principal_bundle_decomposition, PrincipalBundlesWithId, SeqIndexDB, VertexToBundleIdMap,
};
use pgr_db::fasta_io::reverse_complement;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// perform structural variation principle bundle decomposition
#[derive(Parser, Debug)]
#[clap(name = "pgr-sv-decompse")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the data file containing the reference and query sequence of SV candidates
    sv_candidate_seq_path: String,
    /// the prefix of the output files
    output_prefix: String,
    /// the prefix of the output files
    #[clap(long, default_value = "Sample")]
    sample_name: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

#[derive(Debug)]
struct CandidateRecord {
    svc_type: String,
    target_name: String,
    ts: u32,
    te: u32,
    query_name: String,
    qs: u32,
    qe: u32,
    orientation: u8,
    aln_type: String,
    target_sequence: Vec<u8>,
    query_sequence: Vec<u8>,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug, Serialize, Deserialize)]
struct BundleSegment {
    bgn: u32,
    end: u32,
    bundle_id: u32,
    bundle_v_count: u32,
    bundle_dir: u32,
    bundle_v_bgn: u32,
    bundle_v_end: u32,
    is_repeat: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
enum AlnType {
    Match,
    Insertion,
    Deletion,
    Begin,
}
type AlnPathElement = (usize, usize, AlnType, u32, u32, usize, usize);
type AlnPath = Vec<AlnPathElement>;

fn align_bundles(
    q_bundles: &Vec<BundleSegment>,
    t_bundles: &Vec<BundleSegment>,
) -> (f32, usize, usize, AlnPath) {
    let q_count = q_bundles.len();
    let t_count = t_bundles.len();
    let mut s_map = FxHashMap::<(usize, usize), i64>::default();
    let mut t_map = FxHashMap::<(usize, usize), AlnType>::default();

    let mut get_aln_direction_with_best_score =
        |q_idx: usize, t_idx: usize, s_map: &FxHashMap<(usize, usize), i64>| -> (AlnType, i64) {
            let mut best = (AlnType::Match, i64::MIN);
            let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
            let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
            let min_len = if q_len > t_len { t_len } else { q_len };
            let q_b_seg = q_bundles[q_idx];
            let t_b_seg = t_bundles[t_idx];
            if q_idx == 0
                && t_idx == 0
                && (q_b_seg.bundle_id == t_b_seg.bundle_id)
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (AlnType::Match, 2 * min_len)
            } else if q_idx == 0 && t_idx == 0 {
                best = (AlnType::Begin, 0)
            };

            if q_idx > 0
                && t_idx > 0
                && q_b_seg.bundle_id == t_b_seg.bundle_id
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (
                    AlnType::Match,
                    2 * min_len + s_map.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                )
            };
            if t_idx > 0 {
                let score = -2 * q_len + s_map.get(&(q_idx, t_idx - 1)).unwrap();
                if score > best.1 {
                    best = (AlnType::Deletion, score)
                };
            };
            if q_idx > 0 {
                let score = -2 * t_len + s_map.get(&(q_idx - 1, t_idx)).unwrap();
                if score > best.1 {
                    best = (AlnType::Insertion, score)
                }
            }
            t_map.insert((q_idx, t_idx), best.0);
            best
        };

    //let mut best_score = 0;
    //let mut best_q_idx = 0;
    //let mut best_t_idx = 0;
    let mut aln_path = AlnPath::new();

    (0..t_count)
        .flat_map(|t_idx| (0..q_count).map(move |q_idx| (q_idx, t_idx)))
        .for_each(|(q_idx, t_idx)| {
            //println!("{} {}", q_idx, t_idx);
            let (_, score) = get_aln_direction_with_best_score(q_idx, t_idx, &s_map);
            s_map.insert((q_idx, t_idx), score);
            /*
            if score > best_score {
                best_score = score;
                best_q_idx = q_idx;
                best_t_idx = t_idx;
            }
            */
        });
    let mut q_idx = q_count - 1;
    let mut t_idx = t_count - 1;
    let mut diff_len = 0_usize;
    let mut max_len = 1_usize;
    while let Some(aln_type) = t_map.get(&(q_idx, t_idx)) {
        let qq_idx = q_idx;
        let tt_idx = t_idx;
        let (diff_len_delta, max_len_delta) = match aln_type {
            AlnType::Match => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                let diff_len_delta = (q_len - t_len).unsigned_abs() as usize;
                let max_len_delata = if q_len > t_len {
                    q_len as usize
                } else {
                    t_len as usize
                };
                q_idx -= 1;
                t_idx -= 1;
                (diff_len_delta, max_len_delata)
            }
            AlnType::Insertion => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                q_idx -= 1;
                (q_len as usize, q_len as usize)
            }
            AlnType::Deletion => {
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                t_idx -= 1;
                (t_len as usize, t_len as usize)
            }
            AlnType::Begin => break,
        };
        diff_len += diff_len_delta;
        max_len += max_len_delta;
        aln_path.push((
            qq_idx,
            tt_idx,
            *aln_type,
            q_bundles[qq_idx].bundle_id,
            t_bundles[tt_idx].bundle_id,
            diff_len_delta,
            max_len_delta,
        ));
    }
    aln_path.reverse();
    (
        diff_len as f32 / max_len as f32,
        diff_len,
        max_len,
        aln_path,
    )
}

#[allow(clippy::type_complexity)]
fn group_smps_by_principle_bundle_id(
    smps: &[((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)],
    bundle_length_cutoff: usize,
    bundle_merge_distance: usize,
) -> Vec<Vec<((u64, u64, u32, u32, u8), usize, u32, usize)>> {
    let mut pre_bundle_id: Option<usize> = None;
    let mut pre_direction: Option<u32> = None;
    let mut all_partitions = vec![];
    let mut new_partition = vec![];
    smps.iter().for_each(|&(smp, bundle_info)| {
        if bundle_info.is_none() {
            return;
        };
        let bundle_info = bundle_info.unwrap();
        let d = if smp.4 == bundle_info.1 { 0_u32 } else { 1_u32 };
        let bid = bundle_info.0;
        let bpos = bundle_info.2;
        if pre_bundle_id.is_none() {
            new_partition.clear();
            new_partition.push((smp, bid, d, bpos));
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
            return;
        };
        if bid != pre_bundle_id.unwrap() || d != pre_direction.unwrap() {
            let l = new_partition.len();
            if new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
                > bundle_length_cutoff
            {
                all_partitions.push(new_partition.clone());
                new_partition.clear();
            } else {
                new_partition.clear();
            };
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
        };
        new_partition.push((smp, bid, d, bpos));
    });
    let l = new_partition.len();
    if l > 0
        && new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
            > bundle_length_cutoff
    {
        all_partitions.push(new_partition);
    };

    let mut rtn_partitions = vec![];

    if all_partitions.is_empty() {
        return rtn_partitions;
    }
    let mut partition = all_partitions[0].clone();
    (1..all_partitions.len()).for_each(|idx| {
        let p = all_partitions[idx].clone();
        let p_len = partition.len();
        let p_end = partition[p_len - 1].0 .3;
        let p_bid = partition[p_len - 1].1;
        let p_d = partition[p_len - 1].2;
        let np_bgn = p[0].0 .2;
        let np_bid = p[0].1;
        let np_d = p[0].2;
        if p_bid == np_bid
            && p_d == np_d
            && (np_bgn as i64 - p_end as i64).abs() < bundle_merge_distance as i64
        {
            partition.extend(p);
        } else {
            rtn_partitions.push(partition.clone());
            partition = p;
        }
    });
    if !partition.is_empty() {
        rtn_partitions.push(partition);
    }
    rtn_partitions
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let seq_pair_file = BufReader::new(File::open(Path::new(&args.sv_candidate_seq_path)).unwrap());
    let mut paired_seq_records = Vec::new();

    seq_pair_file.lines().for_each(|line| {
        if let Ok(line) = line {
            let fields = line.split('\t').collect::<Vec<&str>>();
            let paser_err_msg = "can't parse the input file";
            assert!(fields.len() == 11);
            let svc_type = fields[0].to_string();
            let target_name = fields[1].to_string();
            let ts = fields[2].parse::<u32>().expect(paser_err_msg);
            let te = fields[3].parse::<u32>().expect(paser_err_msg);
            let query_name = fields[4].to_string();
            let qs = fields[5].parse::<u32>().expect(paser_err_msg);
            let qe = fields[6].parse::<u32>().expect(paser_err_msg);
            let orientation = fields[7].parse::<u8>().expect(paser_err_msg);
            let aln_type = fields[8].to_string();
            let target_sequence = fields[9].into();
            let query_sequence = if orientation == 0 {
                fields[10].into()
            } else {
                reverse_complement(fields[10].as_bytes())
            };

            let rec = CandidateRecord {
                svc_type,
                target_name,
                ts,
                te,
                query_name,
                qs,
                qe,
                orientation,
                aln_type,
                target_sequence,
                query_sequence,
            };
            paired_seq_records.push(rec);
        }
    });

    let mut outpu_bed_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("bed"))
            .expect("can't create the output file"),
    );

    paired_seq_records.into_iter().for_each(|rec| {
        let mut sdb = SeqIndexDB::new();
        let seq_list = vec![
            (rec.target_name.clone(), rec.target_sequence),
            (rec.query_name.clone(), rec.query_sequence),
        ];
        let (w, k, r, min_span, min_cov, min_branch_size) = (31u32, 23u32, 1u32, 13u32, 0u32, 0u32);
        sdb.load_from_seq_list(seq_list, None, w, k, r, min_span)
            .expect("can't load the sequences");
        let (principal_bundles_with_id, vertex_to_bundle_id_direction_pos) = sdb
            .get_principal_bundles_with_id(
                min_cov as usize,
                min_branch_size as usize,
                Some(vec![0, 1]),
            );
        let sid_smps = get_principal_bundle_decomposition(&vertex_to_bundle_id_direction_pos, &sdb);
        let sid_smps: FxHashMap<u32, Vec<_>> = sid_smps.into_iter().collect();
        let bid_to_size = principal_bundles_with_id
            .iter()
            .map(|v| (v.0, v.2.len()))
            .collect::<FxHashMap<usize, usize>>();
        let seq_info = sdb
            .seq_info
            .unwrap()
            .into_iter()
            .map(|(k, v)| (k, v))
            .collect::<Vec<_>>();

        let bundle_length_cutoff = 0;
        let bundle_merge_distance = 0;

        let mut repeat_count = FxHashMap::<u32, Vec<u32>>::default();
        let mut non_repeat_count = FxHashMap::<u32, Vec<u32>>::default();
        let mut sid_to_bundle_segs = FxHashMap::<u32, Vec<BundleSegment>>::default();

        seq_info.iter().for_each(|(sid, sdata)| {
            let (ctg, _src, _len) = sdata;
            let smps = sid_smps.get(sid).unwrap();
            let smp_partitions = group_smps_by_principle_bundle_id(
                smps,
                bundle_length_cutoff,
                bundle_merge_distance,
            );
            let mut ctg_bundle_count = FxHashMap::<usize, usize>::default();
            smp_partitions.iter().for_each(|p| {
                let bid = p[0].1;
                *ctg_bundle_count.entry(bid).or_insert_with(|| 0) += 1;
            });
            let mut bundle_segs = Vec::<BundleSegment>::new();
            smp_partitions.into_iter().for_each(|p| {
                let b = p[0].0 .2;
                let e = p[p.len() - 1].0 .3;
                let bid = p[0].1;
                let direction = p[0].2;
                let is_repeat = if *ctg_bundle_count.get(&bid).unwrap_or(&0) > 1 {
                    repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - k);
                    true
                } else {
                    non_repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - k);
                    false
                };
                bundle_segs.push(BundleSegment {
                    bgn: b,
                    end: e,
                    bundle_id: bid as u32,
                    bundle_v_count: bid_to_size[&bid] as u32,
                    bundle_dir: direction,
                    bundle_v_bgn: p[0].3 as u32,
                    bundle_v_end: p[p.len() - 1].3 as u32,
                    is_repeat,
                });
                // println!(
                //     "{}\t{}\t{}\t{}:{}:{}:{}:{}:{}",
                //     ctg,
                //     b,
                //     e,
                //     bid,
                //     bid_to_size[&bid],
                //     direction,
                //     p[0].3,
                //     p[p.len() - 1].3,
                //     is_repeat
                // );
            });

            sid_to_bundle_segs.insert(*sid, bundle_segs);
        });

        let target_bundles = sid_to_bundle_segs.get(&0).unwrap();
        let query_bundles = sid_to_bundle_segs.get(&1).unwrap();
        let (_dist0, _diff_len0, _max_len0, aln_path) =
            align_bundles(query_bundles, target_bundles);
        println!(
            "##\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.svc_type, rec.target_name, rec.ts, rec.te, rec.query_name, rec.qs, rec.qe, rec.orientation, rec.aln_type
        );
        aln_path.into_iter().for_each(|elm| {
            let (qb_idx, tb_idx, aln_type, _qb_bid, _tb_bid, _diff_delta, _max_diff) = elm;
            let t_seg = target_bundles[tb_idx];
            let q_seg = query_bundles[qb_idx];
            let ts = t_seg.bgn + rec.ts;
            let te = t_seg.end + rec.ts;

            let qs = if rec.orientation == 0 {
                q_seg.bgn + rec.qs
            } else {
                rec.qe - q_seg.end
            };
            let qe = if rec.orientation == 0 {
                q_seg.end + rec.qs
            } else {
                rec.qe - q_seg.bgn
            };
            let target_name = rec.target_name.clone();
            let query_name = rec.query_name.clone();
            println!(
                "{} {} {} {} {} {} {} {} {} {} {} {:?} {:?} {:?}",
                target_name,
                ts,
                te,
                query_name,
                qs,
                qe,
                rec.orientation,
                t_seg.bundle_id,
                t_seg.bundle_dir,
                q_seg.bundle_id,
                q_seg.bundle_dir,
                aln_type,
                t_seg.is_repeat,
                q_seg.is_repeat
            )
        });
    });

    Ok(())
}
