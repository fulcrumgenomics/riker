mod helpers;

use helpers::{FastaBuilder, SamBuilder, SortOrder, coord_builder, read_metrics_tsv};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::gcbias::{
    DETAIL_SUFFIX, GcBias, GcBiasDetailMetric, GcBiasOptions, GcBiasSummaryMetric, SUMMARY_SUFFIX,
};
use tempfile::TempDir;

// ─── Helper ──────────────────────────────────────────────────────────────────

fn make_cmd(
    bam: &std::path::Path,
    ref_fa: &std::path::Path,
    prefix: &std::path::Path,
    exclude_duplicates: bool,
    window_size: u32,
    min_mapq: u8,
    exclude_supplementary: bool,
) -> GcBias {
    GcBias {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: ReferenceOptions { reference: ref_fa.to_path_buf() },
        options: GcBiasOptions { exclude_duplicates, window_size, min_mapq, exclude_supplementary },
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

/// Uniform 0% GC: all-A reference, reads at position 1 → all in gc=0 bin.
/// With `window_size`=10 and reference of 20 A's, there are 11 windows at gc=0.
/// 5 reads start at position 0 → all land in gc=0 bin with `normalized_coverage`=1.0.
#[test]
fn test_uniform_zero_gc() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..5 {
        bld.add_unpaired(&format!("r{i}"), 0, 1, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail.len(), 101);

    // All windows should be at gc=0
    let gc0 = &detail[0];
    assert_eq!(gc0.gc, 0);
    assert!(gc0.windows > 0, "expected windows at gc=0, got {}", gc0.windows);
    assert_eq!(gc0.read_starts, 5);
    // normalized_coverage should be 1.0 since all reads are in gc=0 bin
    assert_float_eq!(gc0.normalized_coverage, 1.0, 0.01);

    // All other bins should have zero windows and reads
    for row in &detail[1..] {
        assert_eq!(row.windows, 0, "gc={} should have 0 windows", row.gc);
        assert_eq!(row.read_starts, 0, "gc={} should have 0 reads", row.gc);
    }

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary.len(), 1);
    assert_eq!(summary[0].total_clusters, 5);
    assert_eq!(summary[0].aligned_reads, 5);
}

/// Mixed GC regions: reference has 100% GC region and 0% GC region.
/// Reads in both → correct bin separation.
#[test]
fn test_mixed_gc_regions() {
    // 10 G's followed by 10 A's = 20bp contig
    let mut ref_seq = vec![b'G'; 10];
    ref_seq.extend_from_slice(&[b'A'; 10]);
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    // Read at pos 1 (0-based 0) → window [0..10] = all G → gc=100%
    bld.add_unpaired("r1", 0, 1, 60, 10, false, false, false, None);
    // Read at pos 11 (0-based 10) → window [10..20] = all A → gc=0%
    bld.add_unpaired("r2", 0, 11, 60, 10, false, false, false, None);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();

    let gc0 = &detail[0];
    let gc100 = &detail[100];
    assert!(gc0.read_starts > 0, "should have reads at gc=0");
    assert!(gc100.read_starts > 0, "should have reads at gc=100");
}

/// N bases in the reference exclude windows.
#[test]
fn test_n_bases_exclude_windows() {
    // 10 A's + 10 N's + 10 A's = 30 bp
    let mut ref_seq = vec![b'A'; 10];
    ref_seq.extend_from_slice(&[b'N'; 10]);
    ref_seq.extend_from_slice(&[b'A'; 10]);
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let bld = coord_builder(&[("chr1", 30)]);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // window_size=10 → windows that overlap the N stretch will have >4 Ns and be excluded
    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();

    // Total windows should be less than 21 (seq.len - window_size + 1) because N windows are excluded
    let total_windows: u64 = detail.iter().map(|r| r.windows).sum();
    // With 10 Ns in the middle and window_size=10, many windows overlap the N region
    assert!(
        total_windows < 21,
        "expected some windows to be excluded due to Ns, got {total_windows}"
    );
}

/// Forward vs reverse strand: verify position calculation for both strands.
#[test]
fn test_forward_vs_reverse_strand() {
    // All G's → 100% GC
    let ref_seq = vec![b'G'; 30];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 30)]);
    // Forward read at pos 1 (0-based 0), read_len=10 → window at pos 0
    bld.add_unpaired("forward", 0, 1, 60, 10, false, false, false, None);
    // Reverse read at pos 1 (0-based 0), read_len=10, CIGAR=10M → ref_span=10, alignment_end=10
    // Position for reverse: alignment_end - window_size = 10 - 10 = 0
    bld.add_unpaired("reverse", 0, 1, 60, 10, true, false, false, None);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();

    // Both reads should land in gc=100 bin
    let gc100 = &detail[100];
    assert_eq!(gc100.read_starts, 2, "both forward and reverse should land in gc=100 bin");

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary[0].aligned_reads, 2);
}

/// Duplicate handling: duplicates included by default, excluded with flag.
#[test]
fn test_duplicate_handling() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    // 2 normal reads + 1 duplicate
    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_unpaired("r1", 0, 1, 60, 10, false, false, false, None);
    bld.add_unpaired("r2", 0, 1, 60, 10, false, false, false, None);
    bld.add_unpaired("dup", 0, 1, 60, 10, false, true, false, None);
    let bam = bld.to_temp_bam().unwrap();

    // Default (include duplicates)
    let dir1 = TempDir::new().unwrap();
    let prefix1 = dir1.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix1, false, 10, 20, false).execute().unwrap();
    let summary1: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir1.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary1[0].aligned_reads, 3);

    // Exclude duplicates
    let dir2 = TempDir::new().unwrap();
    let prefix2 = dir2.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix2, true, 10, 20, false).execute().unwrap();
    let summary2: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir2.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary2[0].aligned_reads, 2);
}

/// GC dropout: reads only in low-GC regions → `gc_dropout` > 0.
#[test]
fn test_gc_dropout() {
    // Half G's, half A's
    let mut ref_seq = vec![b'G'; 50];
    ref_seq.extend_from_slice(&[b'A'; 50]);
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 100)]);
    // Put all reads in the low-GC (all-A) region → position 41 (0-based 40) onwards
    for i in 0..10 {
        bld.add_unpaired(&format!("r{i}"), 0, 41, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // Reads are only in low-GC region, so high-GC windows get no reads → gc_dropout > 0
    assert!(summary[0].gc_dropout > 0.0, "expected gc_dropout > 0, got {}", summary[0].gc_dropout);
}

/// AT dropout: reads only in high-GC regions → `at_dropout` > 0.
#[test]
fn test_at_dropout() {
    // Half A's, half G's
    let mut ref_seq = vec![b'A'; 50];
    ref_seq.extend_from_slice(&[b'G'; 50]);
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 100)]);
    // Put all reads in the high-GC (all-G) region → position 41 (0-based 40) onwards
    for i in 0..10 {
        bld.add_unpaired(&format!("r{i}"), 0, 41, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // Reads are only in high-GC region, so low-GC windows get no reads → at_dropout > 0
    assert!(summary[0].at_dropout > 0.0, "expected at_dropout > 0, got {}", summary[0].at_dropout);
}

/// Empty BAM: no qualifying reads → zero metrics.
#[test]
fn test_empty_bam() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let bld = coord_builder(&[("chr1", 20)]);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail.len(), 101);

    // All read counts should be zero
    for row in &detail {
        assert_eq!(row.read_starts, 0);
        assert_float_eq!(row.normalized_coverage, 0.0, 1e-10);
    }

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary[0].total_clusters, 0);
    assert_eq!(summary[0].aligned_reads, 0);
    assert_float_eq!(summary[0].at_dropout, 0.0, 1e-10);
    assert_float_eq!(summary[0].gc_dropout, 0.0, 1e-10);
}

/// Multiple contigs: windows and reads aggregate correctly.
#[test]
fn test_multiple_contigs() {
    // Two contigs: chr1 all-A (gc=0), chr2 all-G (gc=100)
    let refa = FastaBuilder::new()
        .add_contig("chr1", &[b'A'; 20])
        .add_contig("chr2", &[b'G'; 20])
        .to_temp_fasta()
        .unwrap();

    let mut bld = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    bld.add_unpaired("r1", 0, 1, 60, 10, false, false, false, None);
    bld.add_unpaired("r2", 1, 1, 60, 10, false, false, false, None);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();

    // chr1 windows → gc=0 bin, chr2 windows → gc=100 bin
    let gc0 = &detail[0];
    let gc100 = &detail[100];
    assert!(gc0.windows > 0, "chr1 all-A should produce windows at gc=0");
    assert!(gc100.windows > 0, "chr2 all-G should produce windows at gc=100");
    assert_eq!(gc0.read_starts, 1, "one read on chr1");
    assert_eq!(gc100.read_starts, 1, "one read on chr2");

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary[0].aligned_reads, 2);
}

/// Read filtering: secondary, supplementary, QC-fail, low MAPQ excluded.
#[test]
fn test_read_filtering() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let contigs: Vec<(String, usize)> = vec![("chr1".to_string(), 20)];
    let mut bld = SamBuilder::with_contigs(&contigs).sort_order(SortOrder::Coordinate);

    // Good read
    bld.add_unpaired("good", 0, 1, 60, 10, false, false, false, None);

    // Secondary read — should be excluded
    let secondary = RecordBuf::builder()
        .set_name("secondary")
        .set_flags(Flags::SECONDARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(vec![b'A'; 10]))
        .set_quality_scores(QualityScores::from(vec![30u8; 10]))
        .build();
    bld.add_record(secondary);

    // Supplementary read — included by default
    let supplementary = RecordBuf::builder()
        .set_name("supplementary")
        .set_flags(Flags::SUPPLEMENTARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(vec![b'A'; 10]))
        .set_quality_scores(QualityScores::from(vec![30u8; 10]))
        .build();
    bld.add_record(supplementary);

    // QC-fail read — should be excluded
    bld.add_unpaired("qcfail", 0, 1, 60, 10, false, false, true, None);

    // Low MAPQ read — should be excluded
    bld.add_unpaired("lowmapq", 0, 1, 5, 10, false, false, false, None);

    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // Good read + supplementary read should be counted (supplementary included by default)
    assert_eq!(summary[0].aligned_reads, 2, "good + supplementary reads should be counted");
}

/// Paired reads: first-of-pair counts as cluster, both count as aligned.
#[test]
fn test_paired_reads_cluster_counting() {
    let ref_seq = vec![b'A'; 200];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 200)]);
    // Add 3 FR pairs
    for i in 0..3 {
        let pos1 = (i * 30) + 1;
        let pos2 = pos1 + 20;
        bld.add_pair(&format!("pair{i}"), 0, pos1, pos2, 30, 60, 10, false, false);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // 3 pairs → 3 clusters (first-of-pair), 6 aligned reads
    assert_eq!(summary[0].total_clusters, 3);
    assert_eq!(summary[0].aligned_reads, 6);
}

/// Quintile NC values should be reasonable for uniform coverage.
#[test]
fn test_quintile_nc_uniform() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..5 {
        bld.add_unpaired(&format!("r{i}"), 0, 1, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, false).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // All reads and windows are at gc=0 → gc_0_19_normcov should be 1.0
    assert_float_eq!(summary[0].gc_0_19_normcov, 1.0, 0.01);
    // Other quintiles should be 0 since no windows there
    assert_float_eq!(summary[0].gc_20_39_normcov, 0.0, 1e-10);
    assert_float_eq!(summary[0].gc_40_59_normcov, 0.0, 1e-10);
    assert_float_eq!(summary[0].gc_60_79_normcov, 0.0, 1e-10);
    assert_float_eq!(summary[0].gc_80_100_normcov, 0.0, 1e-10);
}

/// Supplementary reads excluded when `exclude_supplementary` is set.
#[test]
fn test_exclude_supplementary() {
    let ref_seq = vec![b'A'; 20];
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let contigs: Vec<(String, usize)> = vec![("chr1".to_string(), 20)];
    let mut bld = SamBuilder::with_contigs(&contigs).sort_order(SortOrder::Coordinate);

    // Good read
    bld.add_unpaired("good", 0, 1, 60, 10, false, false, false, None);

    // Supplementary read
    let supplementary = RecordBuf::builder()
        .set_name("supplementary")
        .set_flags(Flags::SUPPLEMENTARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(vec![b'A'; 10]))
        .set_quality_scores(QualityScores::from(vec![30u8; 10]))
        .build();
    bld.add_record(supplementary);

    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, false, 10, 20, true).execute().unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // Only the good read should be counted when supplementary reads are excluded
    assert_eq!(summary[0].aligned_reads, 1, "supplementary should be excluded");
}
