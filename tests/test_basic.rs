mod helpers;

use helpers::{SamBuilder, read_metrics_tsv};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use riker_lib::collector::Collector;
use riker_lib::commands::basic::{
    BASE_DIST_PLOT_SUFFIX, BASE_DIST_SUFFIX, BaseDistributionByCycleMetric, BasicCollector,
    MEAN_QUAL_PLOT_SUFFIX, MEAN_QUAL_SUFFIX, MeanQualityByCycleMetric, QUAL_DIST_PLOT_SUFFIX,
    QUAL_DIST_SUFFIX, QualityScoreDistributionMetric,
};
use tempfile::TempDir;

/// Run the BasicCollector on a SamBuilder and return the output prefix.
fn run_basic(builder: &SamBuilder) -> (TempDir, std::path::PathBuf) {
    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    let header = builder.header().clone();

    let mut collector = BasicCollector::new(bam.path(), &prefix);
    collector.initialize(&header).unwrap();

    let mut reader =
        riker_lib::sam::alignment_reader::AlignmentReader::open(bam.path(), None).unwrap();
    let hdr = reader.header().clone();
    let requirements = collector.field_needs();
    for result in reader.riker_records(&requirements) {
        let record = result.unwrap();
        collector.accept(&record, &hdr).unwrap();
    }
    collector.finish().unwrap();

    (dir, prefix)
}

/// Build a RecordBuf with specific sequence and qualities.
fn make_record(name: &str, flags: Flags, seq: &[u8], quals: &[u8]) -> RecordBuf {
    let cigar: Cigar = [Op::new(Kind::Match, seq.len())].into_iter().collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq.to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build()
}

/// Build a paired record with specific sequence and qualities.
fn make_paired_record(
    name: &str,
    is_r2: bool,
    is_reverse: bool,
    seq: &[u8],
    quals: &[u8],
) -> RecordBuf {
    let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
    if is_r2 {
        flags |= Flags::LAST_SEGMENT;
    } else {
        flags |= Flags::FIRST_SEGMENT;
    }
    if is_reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }

    let cigar: Cigar = [Op::new(Kind::Match, seq.len())].into_iter().collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(100).unwrap())
        .set_template_length(if is_r2 { -200 } else { 200 })
        .set_sequence(Sequence::from(seq.to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[test]
fn test_paired_reads() {
    // R1: ACGT with quals 10,20,30,40
    // R2: TTTT with quals 20,20,20,20
    let mut builder = SamBuilder::new();
    builder.add_record(make_paired_record("r1", false, false, b"ACGT", &[10, 20, 30, 40]));
    builder.add_record(make_paired_record("r1", true, true, b"TTTT", &[20, 20, 20, 20]));

    let (_dir, prefix) = run_basic(&builder);

    // Check base distribution
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // R1: 4 cycles, each base appears once at its respective cycle
    let r1: Vec<_> = base_dist.iter().filter(|m| m.read_end == 1).collect();
    assert_eq!(r1.len(), 4);

    // Cycle 1: A=1.0, rest=0
    assert_float_eq!(r1[0].frac_a, 1.0, 1e-5);
    assert_float_eq!(r1[0].frac_c, 0.0, 1e-5);
    // Cycle 2: C=1.0
    assert_float_eq!(r1[1].frac_c, 1.0, 1e-5);
    // Cycle 3: G=1.0
    assert_float_eq!(r1[2].frac_g, 1.0, 1e-5);
    // Cycle 4: T=1.0
    assert_float_eq!(r1[3].frac_t, 1.0, 1e-5);

    // R2: reverse-complemented, so cycle numbering reversed: cycle 4,3,2,1
    // but all bases are T, so each cycle should be 100% T
    let r2: Vec<_> = base_dist.iter().filter(|m| m.read_end == 2).collect();
    assert_eq!(r2.len(), 4);
    for m in &r2 {
        assert_float_eq!(m.frac_t, 1.0, 1e-5);
    }

    // Check mean quality by cycle
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    // R1: 4 cycles, R2: 4 cycles offset by 4
    assert_eq!(mean_qual.len(), 8);
    // R1 cycles
    assert_float_eq!(mean_qual[0].mean_quality, 10.0, 0.01);
    assert_float_eq!(mean_qual[1].mean_quality, 20.0, 0.01);
    assert_float_eq!(mean_qual[2].mean_quality, 30.0, 0.01);
    assert_float_eq!(mean_qual[3].mean_quality, 40.0, 0.01);
    // R2 cycles (all 20, offset cycles 5-8)
    assert_eq!(mean_qual[4].cycle, 5);
    for m in &mean_qual[4..8] {
        assert_float_eq!(m.mean_quality, 20.0, 0.01);
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
    let mut builder = SamBuilder::new();
    // Forward read: ACG, quals 10,20,30
    let flags = Flags::empty();
    builder.add_record(make_record("r1", flags, b"ACG", &[10, 20, 30]));

    let (_dir, prefix) = run_basic(&builder);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // All R1, 3 cycles
    assert_eq!(base_dist.len(), 3);
    assert!(base_dist.iter().all(|m| m.read_end == 1));
    // Cycle 1=A, 2=C, 3=G
    assert_eq!(base_dist[0].cycle, 1);
    assert_float_eq!(base_dist[0].frac_a, 1.0, 1e-5);
    assert_eq!(base_dist[1].cycle, 2);
    assert_float_eq!(base_dist[1].frac_c, 1.0, 1e-5);
    assert_eq!(base_dist[2].cycle, 3);
    assert_float_eq!(base_dist[2].frac_g, 1.0, 1e-5);
}

#[test]
fn test_unpaired_reverse_read() {
    let mut builder = SamBuilder::new();
    // Reverse-complemented read: ACG, quals 10,20,30
    // Cycle numbering reversed: stored seq position 0 -> cycle 3, pos 1 -> cycle 2, pos 2 -> cycle 1
    let flags = Flags::REVERSE_COMPLEMENTED;
    builder.add_record(make_record("r1", flags, b"ACG", &[10, 20, 30]));

    let (_dir, prefix) = run_basic(&builder);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    assert_eq!(base_dist.len(), 3);
    // Stored pos 0 (A) -> cycle_idx 2 (cycle 3)
    // Stored pos 1 (C) -> cycle_idx 1 (cycle 2)
    // Stored pos 2 (G) -> cycle_idx 0 (cycle 1)
    assert_eq!(base_dist[0].cycle, 1);
    assert_float_eq!(base_dist[0].frac_g, 1.0, 1e-5);
    assert_eq!(base_dist[1].cycle, 2);
    assert_float_eq!(base_dist[1].frac_c, 1.0, 1e-5);
    assert_eq!(base_dist[2].cycle, 3);
    assert_float_eq!(base_dist[2].frac_a, 1.0, 1e-5);

    // Quality: pos 0 (q=10) -> cycle 3, pos 1 (q=20) -> cycle 2, pos 2 (q=30) -> cycle 1
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert_eq!(mean_qual.len(), 3);
    assert_float_eq!(mean_qual[0].mean_quality, 30.0, 0.01); // cycle 1
    assert_float_eq!(mean_qual[1].mean_quality, 20.0, 0.01); // cycle 2
    assert_float_eq!(mean_qual[2].mean_quality, 10.0, 0.01); // cycle 3
}

#[test]
fn test_n_bases_in_base_dist_excluded_from_qual_dist() {
    let mut builder = SamBuilder::new();
    // Read with N bases: ANT, quals 10,20,30
    let flags = Flags::empty();
    builder.add_record(make_record("r1", flags, b"ANT", &[10, 20, 30]));

    let (_dir, prefix) = run_basic(&builder);

    // Base distribution should include N
    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert_eq!(base_dist.len(), 3);
    // Cycle 2 should have N=1.0
    assert_float_eq!(base_dist[1].frac_n, 1.0, 1e-5);

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
    let mut builder = SamBuilder::new();

    // One normal read
    builder.add_record(make_record("r1", Flags::empty(), b"ACGT", &[30, 30, 30, 30]));

    // Secondary read — should be skipped
    builder.add_record(make_record("r2", Flags::SECONDARY, b"TTTT", &[30, 30, 30, 30]));

    // Supplementary read — should be skipped
    builder.add_record(make_record("r3", Flags::SUPPLEMENTARY, b"GGGG", &[30, 30, 30, 30]));

    let (_dir, prefix) = run_basic(&builder);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();

    // Only 4 bases from the normal read
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 4);
}

#[test]
fn test_qc_fail_skipped() {
    let mut builder = SamBuilder::new();

    builder.add_record(make_record("r1", Flags::empty(), b"AC", &[30, 30]));
    builder.add_record(make_record("r2", Flags::QC_FAIL, b"GT", &[30, 30]));

    let (_dir, prefix) = run_basic(&builder);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 2); // only r1
}

#[test]
fn test_duplicates_included() {
    let mut builder = SamBuilder::new();

    // Normal read
    builder.add_record(make_record("r1", Flags::empty(), b"AC", &[30, 30]));
    // Duplicate read — should still be counted
    builder.add_record(make_record("r2", Flags::DUPLICATE, b"GT", &[30, 30]));

    let (_dir, prefix) = run_basic(&builder);

    let qd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), QUAL_DIST_SUFFIX));
    let qual_dist: Vec<QualityScoreDistributionMetric> = read_metrics_tsv(&qd_path).unwrap();
    let total: u64 = qual_dist.iter().map(|m| m.count).sum();
    assert_eq!(total, 4); // both reads counted
}

#[test]
fn test_mixed_read_lengths() {
    let mut builder = SamBuilder::new();

    // 3-base read
    builder.add_record(make_record("r1", Flags::empty(), b"ACG", &[10, 20, 30]));
    // 5-base read
    builder.add_record(make_record("r2", Flags::empty(), b"TTTTT", &[40, 40, 40, 40, 40]));

    let (_dir, prefix) = run_basic(&builder);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();

    // Should have 5 cycles (max read length)
    assert_eq!(base_dist.len(), 5);

    // Cycles 1-3: mixed (A+T, C+T, G+T at respective cycles)
    // Cycle 1: A from r1, T from r2 -> frac_a=0.5, frac_t=0.5
    assert_float_eq!(base_dist[0].frac_a, 0.5, 1e-5);
    assert_float_eq!(base_dist[0].frac_t, 0.5, 1e-5);

    // Cycles 4-5: only from r2 (T)
    assert_float_eq!(base_dist[3].frac_t, 1.0, 1e-5);
    assert_float_eq!(base_dist[4].frac_t, 1.0, 1e-5);
}

#[test]
fn test_empty_bam() {
    let builder = SamBuilder::new();
    let (_dir, prefix) = run_basic(&builder);

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
    let mut builder = SamBuilder::new();

    // Two forward reads of length 2: AA and CC
    builder.add_record(make_record("r1", Flags::empty(), b"AA", &[10, 20]));
    builder.add_record(make_record("r2", Flags::empty(), b"CC", &[30, 40]));

    let (_dir, prefix) = run_basic(&builder);

    let bd_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), BASE_DIST_SUFFIX));
    let base_dist: Vec<BaseDistributionByCycleMetric> = read_metrics_tsv(&bd_path).unwrap();
    assert_eq!(base_dist.len(), 2);

    // Each cycle: A=0.5, C=0.5
    assert_float_eq!(base_dist[0].frac_a, 0.5, 1e-5);
    assert_float_eq!(base_dist[0].frac_c, 0.5, 1e-5);
    assert_float_eq!(base_dist[1].frac_a, 0.5, 1e-5);
    assert_float_eq!(base_dist[1].frac_c, 0.5, 1e-5);

    // Mean quality: cycle 1 = (10+30)/2 = 20, cycle 2 = (20+40)/2 = 30
    let mq_path =
        std::path::PathBuf::from(format!("{}{}", prefix.to_str().unwrap(), MEAN_QUAL_SUFFIX));
    let mean_qual: Vec<MeanQualityByCycleMetric> = read_metrics_tsv(&mq_path).unwrap();
    assert_eq!(mean_qual.len(), 2);
    assert_float_eq!(mean_qual[0].mean_quality, 20.0, 0.01);
    assert_float_eq!(mean_qual[1].mean_quality, 30.0, 0.01);
}

#[test]
fn test_plot_files_created() {
    let mut builder = SamBuilder::new();
    builder.add_record(make_record("r1", Flags::empty(), b"ACGT", &[30, 30, 30, 30]));

    let (_dir, prefix) = run_basic(&builder);

    for suffix in [BASE_DIST_PLOT_SUFFIX, MEAN_QUAL_PLOT_SUFFIX, QUAL_DIST_PLOT_SUFFIX] {
        let path = std::path::PathBuf::from(format!("{}{suffix}", prefix.to_str().unwrap()));
        assert!(path.exists(), "Missing plot file: {}", path.display());
        assert!(std::fs::metadata(&path).unwrap().len() > 0, "Empty plot file: {}", path.display());
    }
}
