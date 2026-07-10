mod helpers;

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::basic::{
    BASE_DIST_PLOT_SUFFIX, BASE_DIST_SUFFIX, BaseDistributionByCycleMetric, BasicCollector,
    MEAN_QUAL_PLOT_SUFFIX, MEAN_QUAL_SUFFIX, MeanQualityByCycleMetric, QUAL_DIST_PLOT_SUFFIX,
    QUAL_DIST_SUFFIX, QualityScoreDistributionMetric,
};
use riker_lib::test_support::{ReadBuilder, SamBuilder, drive_collector, pair, read};
use tempfile::TempDir;

/// Run the `BasicCollector` over `sam` and return the temp dir plus output prefix.
fn run_basic(sam: &SamBuilder) -> (TempDir, std::path::PathBuf) {
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    let mut collector = BasicCollector::new(bam.path(), &prefix);
    drive_collector(&mut collector, bam.path()).unwrap();
    (dir, prefix)
}

/// A single mapped read at chr1:1 with an explicit name, bases, and qualities —
/// the base/quality content is exactly what these tests assert on. Callers chain
/// flag setters (`.reverse()`, `.secondary()`, …) and `.build()`.
fn read_with(name: &str, bases: &[u8], quals: &[u8]) -> ReadBuilder {
    read().name(name).bases(bases.to_vec()).quals(quals.to_vec())
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[test]
fn test_paired_reads() {
    // R1: ACGT with quals 10,20,30,40
    // R2: TTTT with quals 20,20,20,20
    let mut sam = SamBuilder::new();
    sam.add(
        pair("r1")
            .r1(|r| r.bases(b"ACGT".to_vec()).quals(vec![10u8, 20, 30, 40]))
            .r2(|r| r.reverse().bases(b"TTTT".to_vec()).quals(vec![20u8; 4])),
    );

    let (_dir, prefix) = run_basic(&sam);

    // Check base distribution
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // R1: 4 cycles, each base appears once at its respective cycle
    let r1: Vec<_> = base_dist.iter().filter(|m| m.read_end == 1).collect();
    assert_eq!(r1.len(), 4);

    // Cycle 1: A=1.0, rest=0
    assert_close!(r1[0].frac_a, 1.0, 1e-5);
    assert_close!(r1[0].frac_c, 0.0, 1e-5);
    // Cycle 2: C=1.0
    assert_close!(r1[1].frac_c, 1.0, 1e-5);
    // Cycle 3: G=1.0
    assert_close!(r1[2].frac_g, 1.0, 1e-5);
    // Cycle 4: T=1.0
    assert_close!(r1[3].frac_t, 1.0, 1e-5);

    // R2: reverse strand, stored as TTTT. Reverse-complementing to sequencing
    // order gives AAAA (complement of T is A), so each cycle should be 100% A.
    let r2: Vec<_> = base_dist.iter().filter(|m| m.read_end == 2).collect();
    assert_eq!(r2.len(), 4);
    for m in &r2 {
        assert_close!(m.frac_a, 1.0, 1e-5);
    }

    // Check mean quality by cycle
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    // R1: 4 cycles, R2: 4 cycles offset by 4
    assert_eq!(mean_qual.len(), 8);
    // R1 cycles
    assert_close!(mean_qual[0].mean_quality, 10.0, 0.01);
    assert_close!(mean_qual[1].mean_quality, 20.0, 0.01);
    assert_close!(mean_qual[2].mean_quality, 30.0, 0.01);
    assert_close!(mean_qual[3].mean_quality, 40.0, 0.01);
    // R2 cycles (all 20, offset cycles 5-8)
    assert_eq!(mean_qual[4].cycle, 5);
    for m in &mean_qual[4..8] {
        assert_close!(m.mean_quality, 20.0, 0.01);
    }

    // Check quality score distribution
    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    // Qualities: 10(1x from R1), 20(1x R1 + 4x R2 = 5), 30(1x), 40(1x)
    let q10 = qual_dist.iter().find(|m| m.quality == 10).unwrap();
    assert_eq!(q10.count, 1);
    let q20 = qual_dist.iter().find(|m| m.quality == 20).unwrap();
    assert_eq!(q20.count, 5);
    let q30 = qual_dist.iter().find(|m| m.quality == 30).unwrap();
    assert_eq!(q30.count, 1);
    let q40 = qual_dist.iter().find(|m| m.quality == 40).unwrap();
    assert_eq!(q40.count, 1);
}

#[test]
fn test_unpaired_forward_read() {
    let mut sam = SamBuilder::new();
    // Forward read: ACG, quals 10,20,30
    sam.add(read_with("r1", b"ACG", &[10, 20, 30]));

    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // All R1, 3 cycles
    assert_eq!(base_dist.len(), 3);
    assert!(base_dist.iter().all(|m| m.read_end == 1));
    // Cycle 1=A, 2=C, 3=G
    assert_eq!(base_dist[0].cycle, 1);
    assert_close!(base_dist[0].frac_a, 1.0, 1e-5);
    assert_eq!(base_dist[1].cycle, 2);
    assert_close!(base_dist[1].frac_c, 1.0, 1e-5);
    assert_eq!(base_dist[2].cycle, 3);
    assert_close!(base_dist[2].frac_g, 1.0, 1e-5);
}

#[test]
fn test_unpaired_reverse_read() {
    let mut sam = SamBuilder::new();
    // Reverse-strand read stored as ACG (quals 10,20,30). The BAM stores it
    // reverse-complemented relative to sequencing order, so the read actually
    // sequenced was revcomp(ACG) = CGT. Both the cycle index and the base
    // must be reverse-complemented to recover it.
    sam.add(read_with("r1", b"ACG", &[10, 20, 30]).reverse());

    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    assert_eq!(base_dist.len(), 3);
    // Stored pos 0 (A) -> complement T -> cycle_idx 2 (cycle 3)
    // Stored pos 1 (C) -> complement G -> cycle_idx 1 (cycle 2)
    // Stored pos 2 (G) -> complement C -> cycle_idx 0 (cycle 1)
    assert_eq!(base_dist[0].cycle, 1);
    assert_close!(base_dist[0].frac_c, 1.0, 1e-5);
    assert_eq!(base_dist[1].cycle, 2);
    assert_close!(base_dist[1].frac_g, 1.0, 1e-5);
    assert_eq!(base_dist[2].cycle, 3);
    assert_close!(base_dist[2].frac_t, 1.0, 1e-5);

    // Quality is reversed but not complemented: pos 0 (q=10) -> cycle 3,
    // pos 1 (q=20) -> cycle 2, pos 2 (q=30) -> cycle 1
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert_eq!(mean_qual.len(), 3);
    assert_close!(mean_qual[0].mean_quality, 30.0, 0.01); // cycle 1
    assert_close!(mean_qual[1].mean_quality, 20.0, 0.01); // cycle 2
    assert_close!(mean_qual[2].mean_quality, 10.0, 0.01); // cycle 3
}

#[test]
fn test_n_bases_in_base_dist_excluded_from_qual_dist() {
    let mut sam = SamBuilder::new();
    // Read with N bases: ANT, quals 10,20,30
    sam.add(read_with("r1", b"ANT", &[10, 20, 30]));

    let (_dir, prefix) = run_basic(&sam);

    // Base distribution should include N
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert_eq!(base_dist.len(), 3);
    // Cycle 2 should have N=1.0
    assert_close!(base_dist[1].frac_n, 1.0, 1e-5);

    // Quality distribution should exclude the N base (q=20)
    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();

    // Only q=10 (A) and q=30 (T) should be present
    assert_eq!(qual_dist.len(), 2);
    let q10 = qual_dist.iter().find(|m| m.quality == 10).unwrap();
    assert_eq!(q10.count, 1);
    let q30 = qual_dist.iter().find(|m| m.quality == 30).unwrap();
    assert_eq!(q30.count, 1);
    // q=20 should be absent (the N base)
    assert!(qual_dist.iter().find(|m| m.quality == 20).is_none());
}

#[test]
fn test_secondary_supplementary_skipped() {
    let mut sam = SamBuilder::new();

    // One normal read
    sam.add(read_with("r1", b"ACGT", &[30, 30, 30, 30]));

    // Secondary read — should be skipped
    sam.add(read_with("r2", b"TTTT", &[30, 30, 30, 30]).secondary());

    // Supplementary read — should be skipped
    sam.add(read_with("r3", b"GGGG", &[30, 30, 30, 30]).supplementary());

    let (_dir, prefix) = run_basic(&sam);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();

    // Only 4 bases from the normal read
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 4);
}

#[test]
fn test_qc_fail_skipped() {
    let mut sam = SamBuilder::new();

    sam.add(read_with("r1", b"AC", &[30, 30]));
    sam.add(read_with("r2", b"GT", &[30, 30]).qc_fail());

    let (_dir, prefix) = run_basic(&sam);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 2); // only r1
}

#[test]
fn test_duplicates_included() {
    let mut sam = SamBuilder::new();

    // Normal read
    sam.add(read_with("r1", b"AC", &[30, 30]));
    // Duplicate read — should still be counted
    sam.add(read_with("r2", b"GT", &[30, 30]).duplicate());

    let (_dir, prefix) = run_basic(&sam);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 4); // both reads counted
}

#[test]
fn test_mixed_read_lengths() {
    let mut sam = SamBuilder::new();

    // 3-base read
    sam.add(read_with("r1", b"ACG", &[10, 20, 30]));
    // 5-base read
    sam.add(read_with("r2", b"TTTTT", &[40, 40, 40, 40, 40]));

    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // Should have 5 cycles (max read length)
    assert_eq!(base_dist.len(), 5);

    // Cycles 1-3: mixed (A+T, C+T, G+T at respective cycles)
    // Cycle 1: A from r1, T from r2 -> frac_a=0.5, frac_t=0.5
    assert_close!(base_dist[0].frac_a, 0.5, 1e-5);
    assert_close!(base_dist[0].frac_t, 0.5, 1e-5);

    // Cycles 4-5: only from r2 (T)
    assert_close!(base_dist[3].frac_t, 1.0, 1e-5);
    assert_close!(base_dist[4].frac_t, 1.0, 1e-5);
}

#[test]
fn test_empty_bam() {
    let sam = SamBuilder::new();
    let (_dir, prefix) = run_basic(&sam);

    // All TSV files should exist and be empty (no data rows)
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert!(base_dist.is_empty());

    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert!(mean_qual.is_empty());

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    assert!(qual_dist.is_empty());
}

#[test]
fn test_multiple_reads_accumulation() {
    let mut sam = SamBuilder::new();

    // Two forward reads of length 2: AA and CC
    sam.add(read_with("r1", b"AA", &[10, 20]));
    sam.add(read_with("r2", b"CC", &[30, 40]));

    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert_eq!(base_dist.len(), 2);

    // Each cycle: A=0.5, C=0.5
    assert_close!(base_dist[0].frac_a, 0.5, 1e-5);
    assert_close!(base_dist[0].frac_c, 0.5, 1e-5);
    assert_close!(base_dist[1].frac_a, 0.5, 1e-5);
    assert_close!(base_dist[1].frac_c, 0.5, 1e-5);

    // Mean quality: cycle 1 = (10+30)/2 = 20, cycle 2 = (20+40)/2 = 30
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert_eq!(mean_qual.len(), 2);
    assert_close!(mean_qual[0].mean_quality, 20.0, 0.01);
    assert_close!(mean_qual[1].mean_quality, 30.0, 0.01);
}

#[test]
fn test_plot_files_created() {
    let mut sam = SamBuilder::new();
    sam.add(read_with("r1", b"ACGT", &[30, 30, 30, 30]));

    let (_dir, prefix) = run_basic(&sam);

    for suffix in [BASE_DIST_PLOT_SUFFIX, MEAN_QUAL_PLOT_SUFFIX, QUAL_DIST_PLOT_SUFFIX] {
        let path = std::path::PathBuf::from(format!("{}{suffix}", prefix.to_str().unwrap()));
        assert!(path.exists(), "Missing plot file: {}", path.display());
        assert!(std::fs::metadata(&path).unwrap().len() > 0, "Empty plot file: {}", path.display());
    }
}

#[test]
fn test_truncated_quality_scores_do_not_panic() {
    // Malformed input: sequence has 4 bases but quality_scores has only 2.
    // The collector hoists the seq/quals length reconciliation out of the
    // hot loop (`n = seq.len().min(quals.len())`) and only walks the
    // common prefix; trailing bases without quality data contribute to
    // neither the base distribution nor the per-cycle quality sums.
    //
    // The pair()/read() builders assert seq.len() == quals.len() (and a BAM
    // writer would reject the mismatch too), so this deliberately-malformed
    // record is built raw and fed to the collector in-process.
    use noodles::core::Position;
    use noodles::sam::Header;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::cigar::{Op, op::Kind};
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use riker_lib::collector::Collector;
    use riker_lib::sam::riker_record::RikerRecord;

    let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
    let buf = RecordBuf::builder()
        .set_name("truncated")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![30u8, 30]))
        .build();
    let header = Header::default();
    let record = RikerRecord::from_alignment_record(&header, &buf).unwrap();

    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    let mut collector = BasicCollector::new(std::path::Path::new("none"), &prefix);
    collector.initialize(&header).unwrap();
    // The line that would have panicked pre-fix:
    collector.accept(&record, &header).unwrap();
    collector.finish().unwrap();

    // Both base distribution and mean quality reflect only the 2 cycles
    // covered by quality data.
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert_eq!(base_dist.len(), 2, "expected 2 cycles in base distribution");

    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert_eq!(mean_qual.len(), 2, "expected 2 cycles with quality data");
}

#[test]
fn test_lowercase_bases_handled_case_insensitively() {
    // The `(base & 0x1F)` indexing folds case automatically, so 'a' and
    // 'A' must produce identical metrics. BAM decode normally yields
    // uppercase, but the collector contract is case-insensitive.
    let mut upper = SamBuilder::new();
    upper.add(read_with("u", b"ACGT", &[30, 30, 30, 30]));
    let (_dir_u, prefix_u) = run_basic(&upper);

    let mut lower = SamBuilder::new();
    lower.add(read_with("l", b"acgt", &[30, 30, 30, 30]));
    let (_dir_l, prefix_l) = run_basic(&lower);

    let bd_u: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&std::path::PathBuf::from(
        format!("{}{}", prefix_u.to_str().unwrap(), BASE_DIST_SUFFIX),
    ))
    .unwrap();
    let bd_l: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&std::path::PathBuf::from(
        format!("{}{}", prefix_l.to_str().unwrap(), BASE_DIST_SUFFIX),
    ))
    .unwrap();

    assert_eq!(bd_u.len(), bd_l.len());
    for (mu, ml) in bd_u.iter().zip(bd_l.iter()) {
        assert_close!(mu.frac_a, ml.frac_a, 1e-9);
        assert_close!(mu.frac_c, ml.frac_c, 1e-9);
        assert_close!(mu.frac_g, ml.frac_g, 1e-9);
        assert_close!(mu.frac_t, ml.frac_t, 1e-9);
        assert_close!(mu.frac_n, ml.frac_n, 1e-9);
    }
}

#[test]
fn test_iupac_ambiguity_codes_excluded_from_qual_dist() {
    // BAM's 4-bit sequence encoding includes the IUPAC ambiguity codes
    // (W, S, M, K, R, Y, B, D, H, V) alongside ACGTN. Picard's
    // QualityScoreDistribution excludes any non-ACGT base from the
    // quality histogram; riker matches that behaviour via the
    // `(ACGT_BITMASK >> bi) & 1` gate. These bases still contribute to
    // the per-cycle base distribution as `frac_n` (the residual bucket).
    let mut sam = SamBuilder::new();
    // Read order: A (Q10), W (Q20 — IUPAC), C (Q30)
    sam.add(read_with("r1", b"AWC", &[10, 20, 30]));
    let (_dir, prefix) = run_basic(&sam);

    // Base distribution: cycle 2 (W) lands in frac_n.
    let bd: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&std::path::PathBuf::from(
        format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX),
    ))
    .unwrap();
    assert_eq!(bd.len(), 3);
    assert_close!(bd[1].frac_n, 1.0, 1e-9);
    assert_close!(bd[1].frac_a, 0.0, 1e-9);

    // Quality distribution: only A's Q10 and C's Q30 — W's Q20 excluded.
    let qd: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&std::path::PathBuf::from(
        format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX),
    ))
    .unwrap();
    assert_eq!(qd.len(), 2, "expected only A's and C's qualities to appear");
    assert!(qd.iter().any(|m| m.quality == 10 && m.count == 1));
    assert!(qd.iter().any(|m| m.quality == 30 && m.count == 1));
    assert!(qd.iter().all(|m| m.quality != 20), "W's Q20 should be excluded");
}

#[test]
fn test_qual_histogram_bank_merge() {
    // The qual_counts is a 4-way interleaved histogram keyed by `i & 3`.
    // A single 8-base read with the same quality at every position
    // distributes counts across all four banks (i=0..7 hits banks
    // 0,1,2,3,0,1,2,3). The merged histogram in the output TSV must
    // report the full count of 8 — an off-by-one in the bank-merge loop
    // would silently drop one or more banks.
    let mut sam = SamBuilder::new();
    sam.add(read_with("r1", b"ACGTACGT", &[35; 8]));
    let (_dir, prefix) = run_basic(&sam);

    let qd: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&std::path::PathBuf::from(
        format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX),
    ))
    .unwrap();
    let q35 = qd.iter().find(|m| m.quality == 35).expect("expected Q35 entry");
    assert_eq!(q35.count, 8, "all 8 bases should be counted across the 4 banks");
}

#[test]
fn test_reverse_strand_with_heterogeneous_qualities() {
    // Walks the `process_record::<true>` (REVERSE = true) path with
    // distinct per-position qualities so that any miscalculation of the
    // reverse cycle index would immediately show up as a per-cycle
    // mean-quality mismatch. BAM stores reverse-complemented reads in
    // the opposite orientation from sequencing, so seq[0] corresponds
    // to the LAST cycle and seq[n-1] to the FIRST cycle.
    let mut sam = SamBuilder::new();
    // Stored seq: A(Q10) C(Q20) G(Q30) T(Q40)
    // After reverse-cycle mapping the per-cycle qualities should be
    // 40, 30, 20, 10 (cycle 1 = seq[3], cycle 2 = seq[2], ...).
    sam.add(read_with("rev", b"ACGT", &[10, 20, 30, 40]).reverse());
    let (_dir, prefix) = run_basic(&sam);

    let mq: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&std::path::PathBuf::from(format!(
        "{}{}",
        prefix.to_str().unwrap(),
        MEAN_QUAL_SUFFIX
    )))
    .unwrap();
    assert_eq!(mq.len(), 4);
    assert_close!(mq[0].mean_quality, 40.0, 0.01);
    assert_close!(mq[1].mean_quality, 30.0, 0.01);
    assert_close!(mq[2].mean_quality, 20.0, 0.01);
    assert_close!(mq[3].mean_quality, 10.0, 0.01);
}

#[test]
fn test_reverse_strand_bases_are_complemented() {
    // A reverse-strand read is stored in the BAM reverse-complemented relative
    // to sequencing order, so the base distribution must complement each base
    // in addition to reversing the cycle index. Stored ACTG on the reverse
    // strand was sequenced as revcomp(ACTG) = CAGT, so cycle 1..4 = C,A,G,T.
    // Four distinct bases make every complement mapping observable.
    let mut sam = SamBuilder::new();
    sam.add(read_with("rev", b"ACTG", &[30, 30, 30, 30]).reverse());
    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    assert_eq!(base_dist.len(), 4);
    assert_close!(base_dist[0].frac_c, 1.0, 1e-5); // cycle 1
    assert_close!(base_dist[1].frac_a, 1.0, 1e-5); // cycle 2
    assert_close!(base_dist[2].frac_g, 1.0, 1e-5); // cycle 3
    assert_close!(base_dist[3].frac_t, 1.0, 1e-5); // cycle 4
}

#[test]
fn test_forward_and_reverse_of_same_read_agree() {
    // Regression for the strand bug: a molecule sequenced as AACG contributes
    // the same per-cycle base distribution whether it aligned to the forward
    // strand (stored AACG) or the reverse strand (stored revcomp(AACG) = CGTT).
    // Before the complement fix the reverse read contributed its complement,
    // smearing each cycle across two bases — the same artifact seen comparing
    // a mapped EM-seq BAM to its RevertSam-reverted counterpart.
    let mut sam = SamBuilder::new();
    sam.add(read_with("fwd", b"AACG", &[30, 30, 30, 30]));
    sam.add(read_with("rev", b"CGTT", &[30, 30, 30, 30]).reverse());
    let (_dir, prefix) = run_basic(&sam);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // Both reads land in the same read-end and cycles, each cycle 100% of the
    // sequenced base: A, A, C, G.
    assert_eq!(base_dist.len(), 4);
    assert_close!(base_dist[0].frac_a, 1.0, 1e-5); // cycle 1
    assert_close!(base_dist[1].frac_a, 1.0, 1e-5); // cycle 2
    assert_close!(base_dist[2].frac_c, 1.0, 1e-5); // cycle 3
    assert_close!(base_dist[3].frac_g, 1.0, 1e-5); // cycle 4
}
