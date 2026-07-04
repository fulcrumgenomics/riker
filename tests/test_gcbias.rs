mod helpers;

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::gcbias::{
    DETAIL_SUFFIX, GcBias, GcBiasDetailMetric, GcBiasOptions, GcBiasSummaryMetric, PLOT_SUFFIX,
    SUMMARY_SUFFIX,
};
use riker_lib::test_support::{BedBuilder, FastaBuilder, coord_builder, pair, read};
use tempfile::TempDir;

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Build a `gcbias` command from a BAM, reference, output prefix, and options.
fn make_cmd(
    bam: &std::path::Path,
    ref_fa: &std::path::Path,
    prefix: &std::path::Path,
    options: GcBiasOptions,
) -> GcBias {
    GcBias {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: ReferenceOptions { reference: ref_fa.to_path_buf() },
        options,
    }
}

/// The standard test options: a small 10 bp window and a MAPQ-20 threshold.
/// Individual tests override single fields via `GcBiasOptions { .., ..opts() }`.
fn opts() -> GcBiasOptions {
    GcBiasOptions { window_size: 10, min_mapq: 20, ..GcBiasOptions::default() }
}

/// Sum the binned read starts across all GC bins (i.e. the reads actually used).
fn reads_used(detail: &[GcBiasDetailMetric]) -> u64 {
    detail.iter().map(|r| r.read_starts).sum()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

/// Uniform 0% GC: all-A reference, reads at position 1 → all in gc=0 bin.
/// With `window_size`=10 and reference of 20 A's, there are 11 windows at gc=0.
/// 5 reads start at position 0 → all land in gc=0 bin with `normalized_coverage`=1.0.
#[test]
fn test_uniform_zero_gc() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(read().name(sam.next_id()).at("chr1", 1).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail.len(), 101);

    // All windows should be at gc=0
    let gc0 = &detail[0];
    assert_eq!(gc0.gc, 0);
    assert!(gc0.windows > 0, "expected windows at gc=0, got {}", gc0.windows);
    assert_eq!(gc0.read_starts, 5);
    // normalized_coverage should be 1.0 since all reads are in gc=0 bin
    assert_close!(gc0.normalized_coverage, 1.0, 0.01);

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
    let refa = FastaBuilder::new().contig("chr1", ref_seq).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    // Read at pos 1 (0-based 0) → window [0..10] = all G → gc=100%
    sam.add(read().name("r1").at("chr1", 1).len(10));
    // Read at pos 11 (0-based 10) → window [10..20] = all A → gc=0%
    sam.add(read().name("r2").at("chr1", 11).len(10));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

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
    let refa = FastaBuilder::new().contig("chr1", ref_seq).to_temp_fasta().unwrap();

    let sam = coord_builder(&[("chr1", 30)]);
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // window_size=10 → windows that overlap the N stretch will have >4 Ns and be excluded
    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

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
    let refa = FastaBuilder::new().contig("chr1", vec![b'G'; 30]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 30)]);
    // Forward read at pos 1 (0-based 0), read_len=10 → window at pos 0
    sam.add(read().name("forward").at("chr1", 1).len(10));
    // Reverse read at pos 1 (0-based 0), read_len=10, CIGAR=10M → ref_span=10, alignment_end=10
    // Position for reverse: alignment_end - window_size = 10 - 10 = 0
    sam.add(read().name("reverse").at("chr1", 1).len(10).reverse());
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

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
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    // 2 normal reads + 1 duplicate
    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(read().name("r1").at("chr1", 1).len(10));
    sam.add(read().name("r2").at("chr1", 1).len(10));
    sam.add(read().name("dup").at("chr1", 1).len(10).duplicate());
    let bam = sam.to_temp_bam().unwrap();

    // Default (include duplicates): all three reads are used.
    let dir1 = TempDir::new().unwrap();
    let prefix1 = dir1.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix1, opts()).execute(None).unwrap();
    let detail1: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir1.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary1: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir1.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary1[0].aligned_reads, 3);
    assert_eq!(reads_used(&detail1), 3);
    assert_eq!(summary1[0].filtered_reads, 0);

    // Exclude duplicates: the duplicate is still aligned, but filtered from the bins.
    let dir2 = TempDir::new().unwrap();
    let prefix2 = dir2.path().join("out");
    make_cmd(
        bam.path(),
        refa.path(),
        &prefix2,
        GcBiasOptions { exclude_duplicates: true, ..opts() },
    )
    .execute(None)
    .unwrap();
    let detail2: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir2.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary2: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir2.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary2[0].aligned_reads, 3, "the duplicate is mapped, so still aligned");
    assert_eq!(reads_used(&detail2), 2, "the duplicate is excluded from the GC bins");
    assert_eq!(summary2[0].filtered_reads, 1);
}

/// GC dropout: reads only in low-GC regions → `gc_dropout` > 0.
#[test]
fn test_gc_dropout() {
    // Half G's, half A's
    let mut ref_seq = vec![b'G'; 50];
    ref_seq.extend_from_slice(&[b'A'; 50]);
    let refa = FastaBuilder::new().contig("chr1", ref_seq).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 100)]);
    // Put all reads in the low-GC (all-A) region → position 41 (0-based 40) onwards
    for _ in 0..10 {
        sam.add(read().name(sam.next_id()).at("chr1", 41).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

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
    let refa = FastaBuilder::new().contig("chr1", ref_seq).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 100)]);
    // Put all reads in the high-GC (all-G) region → position 41 (0-based 40) onwards
    for _ in 0..10 {
        sam.add(read().name(sam.next_id()).at("chr1", 41).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // Reads are only in high-GC region, so low-GC windows get no reads → at_dropout > 0
    assert!(summary[0].at_dropout > 0.0, "expected at_dropout > 0, got {}", summary[0].at_dropout);
}

/// Empty BAM: no qualifying reads → zero metrics.
#[test]
fn test_empty_bam() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let sam = coord_builder(&[("chr1", 20)]);
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail.len(), 101);

    // All read counts should be zero
    for row in &detail {
        assert_eq!(row.read_starts, 0);
        assert_close!(row.normalized_coverage, 0.0, 1e-10);
    }

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary[0].total_clusters, 0);
    assert_eq!(summary[0].aligned_reads, 0);
    assert_close!(summary[0].at_dropout, 0.0, 1e-10);
    assert_close!(summary[0].gc_dropout, 0.0, 1e-10);
}

/// Multiple contigs: windows and reads aggregate correctly.
#[test]
fn test_multiple_contigs() {
    // Two contigs: chr1 all-A (gc=0), chr2 all-G (gc=100)
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'G'; 20])
        .to_temp_fasta()
        .unwrap();

    let mut sam = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    sam.add(read().name("r1").at("chr1", 1).len(10));
    sam.add(read().name("r2").at("chr2", 1).len(10));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

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
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);

    // Good read
    sam.add(read().name("good").at("chr1", 1).len(10));
    // Secondary read — should be excluded
    sam.add(read().name("secondary").at("chr1", 1).len(10).secondary());
    // Supplementary read — included by default
    sam.add(read().name("supplementary").at("chr1", 1).len(10).supplementary());
    // QC-fail read — should be excluded
    sam.add(read().name("qcfail").at("chr1", 1).len(10).qc_fail());
    // Low MAPQ read — should be excluded
    sam.add(read().name("lowmapq").at("chr1", 1).len(10).mapq(5));

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();

    // Secondary and QC-fail reads are discarded entirely (counted nowhere).
    // The remaining good/supplementary/low-MAPQ reads are all mapped, so all
    // three count as aligned; only good + supplementary are actually used.
    assert_eq!(summary[0].total_reads, 3, "good + supplementary + low-MAPQ survive discard");
    assert_eq!(summary[0].aligned_reads, 3, "all three survivors are mapped");
    assert_eq!(reads_used(&detail), 2, "low-MAPQ read is filtered from the GC bins");
    assert_eq!(summary[0].filtered_reads, 1, "the low-MAPQ read is the one filtered");
    assert_close!(summary[0].frac_filtered_reads, 1.0 / 3.0, 1e-5);
}

/// Paired reads: first-of-pair counts as cluster, both count as aligned.
#[test]
fn test_paired_reads_cluster_counting() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 200]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 200)]);
    // Add 3 FR pairs
    for i in 0..3 {
        let pos1 = (i * 30) + 1;
        let pos2 = pos1 + 20;
        sam.add(pair(sam.next_id()).at("chr1", pos1, pos2).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // 3 pairs → 3 clusters (first-of-pair), 6 aligned reads
    assert_eq!(summary[0].total_clusters, 3);
    assert_eq!(summary[0].aligned_reads, 6);
}

/// Quintile NC values should be reasonable for uniform coverage.
#[test]
fn test_quintile_nc_uniform() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(read().name(sam.next_id()).at("chr1", 1).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // All reads and windows are at gc=0 → gc_0_19_normcov should be 1.0
    assert_close!(summary[0].gc_0_19_normcov, 1.0, 0.01);
    // Other quintiles should be 0 since no windows there
    assert_close!(summary[0].gc_20_39_normcov, 0.0, 1e-10);
    assert_close!(summary[0].gc_40_59_normcov, 0.0, 1e-10);
    assert_close!(summary[0].gc_60_79_normcov, 0.0, 1e-10);
    assert_close!(summary[0].gc_80_100_normcov, 0.0, 1e-10);
}

/// Supplementary reads excluded when `exclude_supplementary` is set.
#[test]
fn test_exclude_supplementary() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    // Good read
    sam.add(read().name("good").at("chr1", 1).len(10));
    // Supplementary read
    sam.add(read().name("supplementary").at("chr1", 1).len(10).supplementary());

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(
        bam.path(),
        refa.path(),
        &prefix,
        GcBiasOptions { exclude_supplementary: true, ..opts() },
    )
    .execute(None)
    .unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    // The supplementary read is mapped (so still aligned) but excluded from the
    // GC bins; only the good read is used.
    assert_eq!(summary[0].aligned_reads, 2, "both reads are mapped");
    assert_eq!(reads_used(&detail), 1, "supplementary excluded from the bins");
    assert_eq!(summary[0].filtered_reads, 1);
}

/// Excluded intervals drop both the read (numerator) and its reference window
/// (denominator), keyed on the read's computed start position.
#[test]
fn test_exclude_intervals_drops_reads_and_windows() {
    // All-G reference of length 40 → with window_size 10 there are 31 windows,
    // every one at gc=100.
    let refa = FastaBuilder::new().contig("chr1", vec![b'G'; 40]).to_temp_fasta().unwrap();

    // Read A starts at 0-based 0 (we will exclude it); read B at 0-based 10 (kept).
    let mut sam = coord_builder(&[("chr1", 40)]);
    sam.add(read().name("a_excluded").at("chr1", 1).len(10));
    sam.add(read().name("b_kept").at("chr1", 11).len(10));
    let bam = sam.to_temp_bam().unwrap();

    // Baseline without exclusion: 31 windows, both reads used.
    let dir0 = TempDir::new().unwrap();
    let prefix0 = dir0.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix0, opts()).execute(None).unwrap();
    let detail0: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir0.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail0[100].windows, 31, "31 all-G windows without exclusion");
    assert_eq!(detail0[100].read_starts, 2, "both reads used without exclusion");

    // Exclude window-start positions 0..10 via a BED file.
    let dir = TempDir::new().unwrap();
    let bed = BedBuilder::new().interval("chr1", 0, 10).to_temp_bed().unwrap();
    let prefix = dir.path().join("out");
    make_cmd(
        bam.path(),
        refa.path(),
        &prefix,
        GcBiasOptions { exclude_intervals: Some(bed.path().to_path_buf()), ..opts() },
    )
    .execute(None)
    .unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();

    // Denominator: window starts 0..9 (10 windows) dropped → 21 remain.
    assert_eq!(detail[100].windows, 21, "10 excluded windows removed from the denominator");
    // Numerator: read A (start 0) dropped; read B (start 10) kept.
    assert_eq!(detail[100].read_starts, 1, "read starting in an excluded region is dropped");
    assert_eq!(reads_used(&detail), 1);
    // Both reads are aligned; the excluded one shows up as filtered.
    assert_eq!(summary[0].aligned_reads, 2);
    assert_eq!(summary[0].filtered_reads, 1);
    assert_close!(summary[0].frac_filtered_reads, 0.5, 1e-5);
}

/// The same exclusion expressed as an IntervalList (1-based, inclusive) drops
/// the same windows as the BED form.
#[test]
fn test_exclude_intervals_accepts_interval_list() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'G'; 40]).to_temp_fasta().unwrap();
    let sam = coord_builder(&[("chr1", 40)]);
    let bam = sam.to_temp_bam().unwrap();

    // 0-based [0,10) == 1-based [1,10] inclusive.
    let dir = TempDir::new().unwrap();
    let il = dir.path().join("exclude.interval_list");
    std::fs::write(&il, "@SQ\tSN:chr1\tLN:40\nchr1\t1\t10\t+\texcl\n").unwrap();
    let prefix = dir.path().join("out");
    make_cmd(
        bam.path(),
        refa.path(),
        &prefix,
        GcBiasOptions { exclude_intervals: Some(il), ..opts() },
    )
    .execute(None)
    .unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    assert_eq!(detail[100].windows, 21, "IntervalList exclusion drops the same 10 windows");
}

/// Exclusion is keyed on a reverse-strand read's *computed* start position
/// (`alignment_end - window_size`), not on whether its span overlaps the
/// interval. This is the distinguishing behavior vs. an overlap-based filter.
#[test]
fn test_exclude_intervals_reverse_strand_uses_computed_start() {
    // All-G reference, length 40, window_size 10 → 31 windows, all gc=100.
    let refa = FastaBuilder::new().contig("chr1", vec![b'G'; 40]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 40)]);
    // Reverse read, 20M at 1-based pos 1 → spans 0-based [0,20), which OVERLAPS
    // the excluded region [0,10). But its computed start is
    // alignment_end(20) - window_size(10) = 0-based 10, which is NOT excluded →
    // an overlap filter would drop it; start-position keying KEEPS it.
    sam.add(read().name("rev_kept").at("chr1", 1).len(20).reverse());
    // Reverse read, 5M at 1-based pos 11 → spans 0-based [10,15), which does NOT
    // overlap [0,10). But its computed start is end(15) - 10 = 0-based 5, which
    // IS excluded → an overlap filter would keep it; start-position keying DROPS
    // it.
    sam.add(read().name("rev_dropped").at("chr1", 11).len(5).reverse());
    let bam = sam.to_temp_bam().unwrap();

    let dir = TempDir::new().unwrap();
    let bed = BedBuilder::new().interval("chr1", 0, 10).to_temp_bed().unwrap();
    let prefix = dir.path().join("out");
    make_cmd(
        bam.path(),
        refa.path(),
        &prefix,
        GcBiasOptions { exclude_intervals: Some(bed.path().to_path_buf()), ..opts() },
    )
    .execute(None)
    .unwrap();

    let detail: Vec<GcBiasDetailMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{DETAIL_SUFFIX}"))).unwrap();
    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();

    // Only rev_kept is binned; rev_dropped is filtered. Both are aligned.
    assert_eq!(reads_used(&detail), 1, "only the read whose computed start is unexcluded is used");
    assert_eq!(detail[100].read_starts, 1);
    assert_eq!(summary[0].aligned_reads, 2);
    assert_eq!(summary[0].filtered_reads, 1);
}

/// Unmapped reads are counted in total_reads / total_clusters but not in
/// aligned_reads, and they never inflate filtered_reads.
#[test]
fn test_unmapped_counts_in_total_not_aligned() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(read().name("mapped").at("chr1", 1).len(10));
    sam.add(read().name("unmapped").unmapped().len(10));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let summary: Vec<GcBiasSummaryMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{SUMMARY_SUFFIX}"))).unwrap();
    assert_eq!(summary[0].total_reads, 2, "both reads counted in total_reads");
    assert_eq!(summary[0].total_clusters, 2, "both reads counted as clusters");
    assert_eq!(summary[0].aligned_reads, 1, "only the mapped read is aligned");
    assert_eq!(summary[0].filtered_reads, 0, "the unmapped read never enters the aligned funnel");
}

/// The diagnostic GC-bias chart PDF is written alongside the metrics.
#[test]
fn test_writes_chart_pdf() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(read().name(sam.next_id()).at("chr1", 1).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, opts()).execute(None).unwrap();

    let plot_path = dir.path().join(format!("out{PLOT_SUFFIX}"));
    assert!(plot_path.exists(), "GC-bias chart PDF should be created");
    assert!(plot_path.metadata().unwrap().len() > 0, "GC-bias chart PDF should be non-empty");
}
