mod helpers;

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::wgs::{
    COVERAGE_SUFFIX, METRICS_SUFFIX, Wgs, WgsCoverageEntry, WgsMetrics, WgsOptions,
};
use riker_lib::test_support::{FastaBuilder, ReadBuilder, coord_builder, pair, read};
use tempfile::TempDir;

// ─── Helper ──────────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn make_cmd(
    bam: &std::path::Path,
    ref_fa: &std::path::Path,
    prefix: &std::path::Path,
    intervals: Option<std::path::PathBuf>,
    include_duplicates: bool,
    include_unpaired_reads: bool,
    min_mapq: u8,
    min_bq: u8,
    coverage_cap: u16,
) -> Wgs {
    Wgs {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: ReferenceOptions { reference: ref_fa.to_path_buf() },
        options: WgsOptions {
            intervals,
            include_duplicates,
            include_unpaired_reads,
            min_mapq,
            min_bq,
            coverage_cap,
        },
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

/// 5 FR pairs, non-overlapping: r1 covers pos 0–9, r2 covers 10–19 on 20 bp reference.
/// Each position is covered by exactly 5 reads → mean=5, `frac_at_5x`=1.
#[test]
fn test_basic_coverage() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        bld.add(pair(bld.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];

    assert_eq!(m.genome_territory, 20);
    assert_close!(m.mean_coverage, 5.0, 0.01);
    assert_close!(m.frac_bases_at_1x, 1.0, 0.001);
    assert_close!(m.frac_bases_at_5x, 1.0, 0.001);
    assert_close!(m.frac_bases_at_10x, 0.0, 0.001);
}

/// N bases in the reference should be excluded from genome territory.
#[test]
fn test_n_bases_excluded() {
    let seq: Vec<u8> = std::iter::repeat_n(b'A', 10).chain(std::iter::repeat_n(b'N', 10)).collect();
    let refa = FastaBuilder::new().contig("chr1", seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for _ in 0..3 {
        bld.add(pair(bld.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    // Only the 10 non-N bases count toward genome territory.
    assert_eq!(rows[0].genome_territory, 10);
}

/// Overlapping FR pair: reads share a name at the overlap positions, so one
/// occurrence is excluded as overlap.
#[test]
fn test_overlap_detection() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    // r1 covers 0–9, r2 covers 5–14 → 5 positions of overlap (5–9)
    bld.add(pair("overlapper").at("chr1", 1, 6).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_overlap > 0.0,
        "expected overlap exclusions, got {}",
        rows[0].frac_excluded_overlap
    );
}

/// Low-MAPQ reads should be excluded; high-MAPQ reads count.
#[test]
fn test_mapq_exclusion() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(pair("lowmq").at("chr1", 1, 11).len(10).mapq(10)); // mapq=10 < 20
    bld.add(pair("highmq").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_mapq > 0.0,
        "expected mapq exclusions, got {}",
        rows[0].frac_excluded_mapq
    );
    assert_close!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Duplicate reads are excluded when `include_duplicates`=false (the default).
#[test]
fn test_dup_exclusion() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(
        pair("dup").at("chr1", 1, 11).len(10).r1(ReadBuilder::duplicate).r2(ReadBuilder::duplicate),
    );
    bld.add(pair("normal").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_dupe > 0.0,
        "expected dupe exclusions, got {}",
        rows[0].frac_excluded_dupe
    );
    assert_close!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Unpaired reads are excluded by default (--include-unpaired-reads not set).
#[test]
fn test_unpaired_excluded_by_default() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(read().name("frag").at("chr1", 1).len(20));
    bld.add(pair("pair").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_unpaired > 0.0,
        "expected unpaired exclusions, got {}",
        rows[0].frac_excluded_unpaired
    );
    // Only the pair (depth 1 across all 20 positions) counts; the unpaired
    // fragment is dropped, so mean coverage is 1.0.
    assert_close!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Unpaired reads are counted when --include-unpaired-reads is set, and are
/// no longer tallied as unpaired exclusions.
#[test]
fn test_unpaired_included_when_flag_set() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(read().name("frag").at("chr1", 1).len(20));
    bld.add(pair("pair").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert_close!(rows[0].frac_excluded_unpaired, 0.0, 1e-9);
    // The unpaired fragment now contributes depth 1 across all 20 positions on
    // top of the pair's depth 1, so mean coverage doubles to 2.0. This proves
    // the read's bases are actually counted, not merely un-excluded.
    assert_close!(rows[0].mean_coverage, 2.0, 0.01);
}

/// Depth exceeding `coverage_cap` produces `frac_excluded_capped` > 0.
#[test]
fn test_coverage_cap() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        bld.add(pair(bld.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // cap=3 → 5 reads per position, 2 excluded per position
    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 3).execute(None).unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_capped > 0.0,
        "expected capped exclusions, got {}",
        rows[0].frac_excluded_capped
    );
}

/// The coverage histogram file should have `coverage_cap`+1 rows and sum to `genome_territory`.
#[test]
fn test_coverage_histogram_file() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(pair("r0").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let cov: Vec<WgsCoverageEntry> =
        read_metrics_tsv(&dir.path().join(format!("out{COVERAGE_SUFFIX}"))).unwrap();

    assert_eq!(cov.len(), 251, "expected 251 rows (depth 0–250)");

    let total_bases: u64 = cov.iter().map(|r| r.bases).sum();
    assert_eq!(total_bases, 20, "bases column must sum to genome_territory");

    // One pair covers both halves → all 20 positions at depth=1, depth=0 nowhere.
    assert_eq!(cov[0].depth, 0);
    assert_eq!(cov[0].bases, 0);
    assert_eq!(cov[1].depth, 1);
    assert_eq!(cov[1].bases, 20);

    // bases_at_or_above[0] = all 20 positions; [1] = 20; [2] = 0
    assert_eq!(cov[0].bases_at_or_above, 20);
    assert_eq!(cov[1].bases_at_or_above, 20);
    assert_eq!(cov[2].bases_at_or_above, 0);
}

/// Known coverage distribution → verify `frac_bases_at_Nx` values.
#[test]
fn test_frac_bases_at_nx() {
    // 20 bp reference; first 10 bp covered by 15 reads, last 10 bp by 5 reads.
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for _ in 0..15 {
        bld.add(read().name(bld.next_id()).at("chr1", 1).len(10));
    }
    for _ in 0..5 {
        bld.add(read().name(bld.next_id()).at("chr1", 11).len(10));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let m = &rows[0];

    assert_close!(m.frac_bases_at_1x, 1.0, 0.001);
    assert_close!(m.frac_bases_at_5x, 1.0, 0.001);
    assert_close!(m.frac_bases_at_10x, 0.5, 0.001); // only 10/20 positions
    assert_close!(m.frac_bases_at_15x, 0.5, 0.001);
    assert_close!(m.frac_bases_at_20x, 0.0, 0.001);
    // mean = (15*10 + 5*10) / 20 = 10
    assert_close!(m.mean_coverage, 10.0, 0.01);
}

/// Uniform coverage → `fold_80_base_penalty` should equal or approach 1.0.
#[test]
fn test_fold_penalty() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 10]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 10)]);
    for _ in 0..10 {
        bld.add(read().name(bld.next_id()).at("chr1", 1).len(10));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let m = &rows[0];

    assert_close!(m.mean_coverage, 10.0, 0.01);
    // Uniform coverage → fold_80 = mean / 20th-pct = 10/10 = 1.0
    assert!(
        m.fold_80_base_penalty >= 1.0,
        "fold_80 should be ≥ 1.0, got {}",
        m.fold_80_base_penalty
    );
    assert_close!(m.fold_80_base_penalty, 1.0, 0.01);
}

/// Interval restriction: `genome_territory` should equal only the in-interval non-N bases.
#[test]
fn test_with_intervals() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add(pair("r0").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();

    // BED interval: chr1 positions 5–14 (0-based half-open → 10 positions)
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("test.bed");
    std::fs::write(&bed_path, "chr1\t5\t15\n").unwrap();

    let prefix = dir.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix, Some(bed_path), false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].genome_territory, 10);
}

/// Bases with quality below `min_bq` should be counted as baseq-excluded.
#[test]
fn test_baseq_exclusion() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 10]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 10)]);
    // Low base quality read (qual=5 < min_bq=20).
    bld.add(read().name("lowbq").at("chr1", 1).len(10).qual(5));
    // High quality read covering the same positions.
    bld.add(read().name("highbq").at("chr1", 1).len(10));

    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_baseq > 0.0,
        "expected baseq exclusions, got {}",
        rows[0].frac_excluded_baseq
    );
    assert_close!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Validates that `sd_coverage` is computed using capped depths (not raw depths)
/// and uses sample standard deviation (n-1 denominator).
///
/// Setup: 20bp all-A reference, `coverage_cap=5`.
/// 8 read pairs with `read_len=5`: r1 covers pos 1-5, r2 covers pos 6-10.
/// Raw depth: 8 at positions 1-10, 0 at positions 11-20.
/// Capped depth: 5 at positions 1-10, 0 at positions 11-20.
/// Expected mean = (10*5 + 10*0) / 20 = 2.5
/// Expected sample variance = (10*(5-2.5)^2 + 10*(0-2.5)^2) / 19 = 125/19
/// Expected sd = sqrt(125/19) ≈ 2.5643
#[test]
fn test_sd_coverage_uses_capped_depth() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for _ in 0..8 {
        bld.add(pair(bld.next_id()).at("chr1", 1, 6).len(5));
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 5).execute(None).unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];

    assert_eq!(m.genome_territory, 20);
    assert_close!(m.mean_coverage, 2.5, 0.01);

    // sd from capped depths: sqrt(125/19) ≈ 2.5643
    let expected_sd = (125.0_f64 / 19.0).sqrt();
    assert_close!(m.sd_coverage, expected_sd, 0.01);
}

/// A second contig that receives no reads exercises the zero-coverage path in
/// `finish_metrics`: its non-N bases still count toward genome territory (all at
/// depth 0). Covers the `None` (no-intervals) branch of `count_eligible_positions`,
/// which the single-contig tests never reach.
#[test]
fn zero_coverage_contig_counts_non_n_territory() {
    // chr2: 10 bases then 10 Ns → 10 non-N positions, none covered.
    let chr2_seq: Vec<u8> =
        std::iter::repeat_n(b'A', 10).chain(std::iter::repeat_n(b'N', 10)).collect();
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", chr2_seq)
        .to_temp_fasta()
        .unwrap();

    // One pair on chr1 (depth 1 across all 20 bp); chr2 gets no reads at all.
    let mut bld = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    bld.add(pair("r0").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    // chr1: 20 non-N covered; chr2: 10 non-N uncovered → 30 total territory.
    assert_eq!(rows[0].genome_territory, 30);

    let cov: Vec<WgsCoverageEntry> =
        read_metrics_tsv(&dir.path().join(format!("out{COVERAGE_SUFFIX}"))).unwrap();
    // chr2's 10 non-N bases land in the depth-0 bin; chr1's 20 sit at depth 1.
    assert_eq!(cov[0].bases, 10);
    assert_eq!(cov[1].bases, 20);
}

/// A read-free contig that is both interval-restricted and contains N gaps
/// exercises the `Some` branch of `count_eligible_positions`: only non-N
/// positions that also fall inside the interval mask are counted (the rewritten
/// run-vs-interval intersection).
#[test]
fn zero_coverage_contig_with_intervals_counts_masked_non_n() {
    // chr2: 5 bases, 5 Ns, 5 bases → non-N runs [0,5) and [10,15).
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", b"AAAAANNNNNAAAAA".to_vec())
        .to_temp_fasta()
        .unwrap();

    // Reads only on chr1; chr2 is the read-free, interval-restricted contig.
    let mut bld = coord_builder(&[("chr1", 20), ("chr2", 15)]);
    bld.add(pair("r0").at("chr1", 1, 11).len(10));
    let bam = bld.to_temp_bam().unwrap();

    // Interval on chr2 only, [3, 13). chr1 has no interval, so it contributes 0.
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("chr2.bed");
    std::fs::write(&bed_path, "chr2\t3\t13\n").unwrap();

    let prefix = dir.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix, Some(bed_path), false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    // non-N ∩ [3,13): run [0,5)→{3,4}=2, run [10,15)→{10,11,12}=3 → 5 total.
    assert_eq!(rows[0].genome_territory, 5);
}
