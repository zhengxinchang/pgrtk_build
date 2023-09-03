const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use iset::set::IntervalSet;
// use rayon::prelude::*;
use pgr_db::ext::{
    get_principal_bundle_decomposition, PrincipalBundlesWithId, SeqIndexDB, VertexToBundleIdMap,
};
use rustc_hash::{FxHashMap, FxHashSet};
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
struct CandidateRecord {
    target_name: String,
    ts: u32,
    te: u32,
    query_name: String,
    qs: u32,
    qe: u32,
    orientation: u8,
    target_sequence: Vec<u8>,
    query_sequence: Vec<u8>,
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
            assert!(fields.len() == 9);
            let target_name = fields[0].to_string();
            let ts = fields[1].parse::<u32>().expect(paser_err_msg);
            let te = fields[2].parse::<u32>().expect(paser_err_msg);
            let query_name = fields[3].to_string();
            let qs = fields[4].parse::<u32>().expect(paser_err_msg);
            let qe = fields[5].parse::<u32>().expect(paser_err_msg);
            let orientation = fields[6].parse::<u8>().expect(paser_err_msg);
            let target_sequence = fields[7].into();
            let query_sequence = fields[8].into();
            let rec = CandidateRecord {
                target_name,
                ts,
                te,
                query_name,
                qs,
                qe,
                orientation,
                target_sequence,
                query_sequence,
            };
            paired_seq_records.push(rec);
        }
    });

    let mut outpu_bed_file =
    BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("bed")).expect("can't create the output file"));


    paired_seq_records.into_iter().for_each(|rec| {
        let mut sdb = SeqIndexDB::new();
        let seq_list = vec![
            (rec.target_name, rec.target_sequence),
            (rec.query_name, rec.query_sequence),
        ];
        let (w, k, r, min_span, min_cov, min_branch_size) = (17u32, 31u32, 1u32, 0u32, 0u32, 0u32);
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
            smp_partitions.into_iter().for_each(|p| {
                let b = p[0].0 .2;
                let e = p[p.len() - 1].0 .3 + k;
                let bid = p[0].1;
                let direction = p[0].2;
                let is_repeat = if *ctg_bundle_count.get(&bid).unwrap_or(&0) > 1 {
                    repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - k);
                    "R"
                } else {
                    non_repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - k);
                    "U"
                };
                let _ = writeln!(
                    outpu_bed_file,
                    "{}\t{}\t{}\t{}:{}:{}:{}:{}:{}",
                    ctg,
                    b,
                    e,
                    bid,
                    bid_to_size[&bid],
                    direction,
                    p[0].3,
                    p[p.len() - 1].3,
                    is_repeat
                );
            });
        });
    });

    Ok(())
}
