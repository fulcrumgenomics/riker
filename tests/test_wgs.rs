mod helpers;

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use bigwig::BigWigReader;
use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::wgs::{
    COVERAGE_SUFFIX, METRICS_SUFFIX, PER_BASE_BIGWIG_SUFFIX, PerBaseCoverageFormat, Wgs,
    WgsCoverageEntry, WgsMetrics, WgsOptions,
};
use riker_lib::test_support::{
    BedBuilder, FastaBuilder, ReadBuilder, SamBuilder, SortOrder, coord_builder, pair, read,
};
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
            per_base_coverage: None,
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

/// A FASTA with chr1 then chr2, a BAM listing chr2 then chr1 with 5× coverage on chr2 and 1× on
/// chr1, and a BED over all of chr1, all written into `dir`.
fn swapped_contig_order_inputs(dir: &TempDir) -> (tempfile::NamedTempFile, SamBuilder, PathBuf) {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr2", 20), ("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr2", 1, 11).len(10));
    }
    sam.add(pair("on_chr1").at("chr1", 1, 11).len(10));
    let bed_path = dir.path().join("chr1.bed");
    BedBuilder::new().interval("chr1", 0, 20).write_bed(&bed_path).unwrap();
    (refa, sam, bed_path)
}

/// Intervals restrict their own contig even when the BAM lists contigs in a different order from
/// the FASTA.
#[test]
fn intervals_follow_the_bam_contig_order() {
    let dir = TempDir::new().unwrap();
    let (refa, sam, bed_path) = swapped_contig_order_inputs(&dir);
    let bam = sam.to_temp_bam().unwrap();

    let prefix = dir.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix, Some(bed_path), false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows[0].genome_territory, 20);
    assert_close!(rows[0].mean_coverage, 1.0, 1e-9);
}

// ─── Per-base coverage output ─────────────────────────────────────────────────

/// Run `wgs` with `--per-base-coverage bigwig` and default filters, writing into `dir`, and
/// open the bigWig it wrote.
fn run_per_base(
    sam: &SamBuilder,
    ref_fa: &Path,
    intervals: Option<PathBuf>,
    coverage_cap: u16,
    dir: &TempDir,
) -> BigWigReader<BufReader<File>> {
    let bam = sam.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");
    let mut cmd =
        make_cmd(bam.path(), ref_fa, &prefix, intervals, false, true, 20, 20, coverage_cap);
    cmd.options.per_base_coverage = Some(PerBaseCoverageFormat::Bigwig);
    cmd.execute(None).unwrap();
    BigWigReader::open(dir.path().join(format!("out{PER_BASE_BIGWIG_SUFFIX}"))).unwrap()
}

/// `n` copies of `value`, for comparing per-base values.
fn repeat(value: f32, n: usize) -> Vec<f32> {
    vec![value; n]
}

/// `true` if every value is NaN, which the reader returns where the file has no data.
fn all_missing(values: &[f32]) -> bool {
    values.iter().all(|value| value.is_nan())
}

/// The items stored for `chrom` as `(start, end, value)`: one per position of fixedStep data, one
/// per zero-valued bedGraph item.
fn items(reader: &mut BigWigReader<BufReader<File>>, chrom: &str) -> Vec<(u32, u32, f32)> {
    reader
        .chrom_items(chrom)
        .unwrap()
        .map(|item| {
            let item = item.unwrap();
            (item.start, item.end, item.value)
        })
        .collect()
}

/// A covered contig gets a value at every position, including zero where no read lands.
#[test]
fn per_base_coverage_writes_the_depth_at_every_position_of_a_covered_contig() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 30]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 30)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    let expected: Vec<f32> = [repeat(5.0, 20), repeat(0.0, 10)].concat();
    assert_eq!(reader.values("chr1", 0, 30).unwrap(), expected);
}

/// Where the two mates of a pair overlap, the depth counts the pair once, as the metrics do.
#[test]
fn per_base_coverage_counts_overlapping_mates_once() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    // r1 covers 0–9 and r2 covers 5–14, overlapping at 5–9.
    sam.add(pair("overlapper").at("chr1", 1, 6).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    let expected: Vec<f32> = [repeat(1.0, 15), repeat(0.0, 5)].concat();
    assert_eq!(reader.values("chr1", 0, 20).unwrap(), expected);
}

/// Duplicates, low-MAPQ reads and low base qualities add no depth, as in the metrics.
#[test]
fn per_base_coverage_excludes_what_the_metrics_exclude() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(
        pair("dup").at("chr1", 1, 11).len(10).r1(ReadBuilder::duplicate).r2(ReadBuilder::duplicate),
    );
    sam.add(pair("low_mapq").at("chr1", 1, 11).len(10).mapq(10));
    sam.add(pair("low_bq").at("chr1", 1, 11).len(10).quals(vec![5; 10]));
    sam.add(pair("kept").at("chr1", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(1.0, 20));
}

/// Per-base depth is reported in full above `--coverage-cap`, while the metrics still cap it.
#[test]
fn per_base_coverage_is_not_capped_by_coverage_cap() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 3, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(5.0, 20));
    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_close!(rows[0].mean_coverage, 3.0, 1e-9);
}

/// Per-base depth is stored in 16 bits, so a position deeper than 65,535 is reported as 65,535.
#[test]
fn per_base_coverage_reports_depth_above_65535_as_65535() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for i in 0..65_540 {
        sam.add(read().name(format!("deep{i}")).at("chr1", 1).len(5));
    }
    sam.add(read().name("shallow").at("chr1", 6).len(5));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    let expected: Vec<f32> = [repeat(65_535.0, 5), repeat(1.0, 5), repeat(0.0, 10)].concat();
    assert_eq!(reader.values("chr1", 0, 20).unwrap(), expected);
}

/// A contig that no read reaches, between two that reads do, has zero depth at every position.
#[test]
fn per_base_coverage_writes_zeros_on_a_contig_without_reads_between_covered_contigs() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .contig("chr3", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 20), ("chr2", 20), ("chr3", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    sam.add(pair("r1").at("chr3", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(1.0, 20));
    assert_eq!(reader.values("chr2", 0, 20).unwrap(), repeat(0.0, 20));
    assert_eq!(reader.values("chr3", 0, 20).unwrap(), repeat(1.0, 20));
}

/// Contigs after the last one that reads reach have zero depth at every position.
#[test]
fn per_base_coverage_writes_zeros_on_contigs_after_the_last_covered_contig() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(1.0, 20));
    assert_eq!(reader.values("chr2", 0, 20).unwrap(), repeat(0.0, 20));
}

/// A contig whose reads all fail the filters has zero depth at every position.
#[test]
fn per_base_coverage_writes_zeros_on_a_contig_whose_reads_are_all_filtered() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    sam.add(pair("low_bq").at("chr2", 1, 11).len(10).quals(vec![5; 10]));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr2", 0, 20).unwrap(), repeat(0.0, 20));
}

/// A BAM with no reads still gets zero depth at every position of every contig.
#[test]
fn per_base_coverage_writes_zeros_on_every_contig_of_a_bam_without_reads() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let sam = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(0.0, 20));
    assert_eq!(reader.values("chr2", 0, 20).unwrap(), repeat(0.0, 20));
}

/// A contig without coverage is stored as one zero-valued item, whether no read reached it or its
/// reads were all filtered.
#[test]
fn per_base_coverage_stores_a_contig_without_coverage_as_one_zero_item() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .contig("chr3", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 20), ("chr2", 20), ("chr3", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    sam.add(pair("low_bq").at("chr3", 1, 11).len(10).quals(vec![5; 10]));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(items(&mut reader, "chr2"), [(0, 20, 0.0)]);
    assert_eq!(items(&mut reader, "chr3"), [(0, 20, 0.0)]);
}

/// With `--intervals`, only interval positions get data; overlapping intervals are merged.
#[test]
fn per_base_coverage_with_intervals_writes_only_the_interval_positions() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("targets.bed");
    std::fs::write(&bed_path, "chr1\t5\t12\nchr1\t10\t15\n").unwrap();

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let values = reader.values("chr1", 0, 20).unwrap();
    assert!(all_missing(&values[..5]));
    assert_eq!(values[5..15], repeat(1.0, 10));
    assert!(all_missing(&values[15..]));
}

/// The two contig and interval layouts of the tests of intervals without coverage: a 40 bp chr1 whose
/// one pair covers positions 0–19, with intervals at 5–15 and 25–35, and a 20 bp chr2 no read
/// reaches, with an interval at 5–15.
fn intervals_without_coverage_inputs(
    dir: &TempDir,
) -> (tempfile::NamedTempFile, SamBuilder, PathBuf) {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 40])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 40), ("chr2", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    let bed_path = dir.path().join("targets.bed");
    BedBuilder::new()
        .interval("chr1", 5, 15)
        .interval("chr1", 25, 35)
        .interval("chr2", 5, 15)
        .write_bed(&bed_path)
        .unwrap();
    (refa, sam, bed_path)
}

/// With `--intervals`, an interval no read covers has zero depth at every position, on a contig
/// that reads reach elsewhere.
#[test]
fn per_base_coverage_with_intervals_writes_zeros_on_an_interval_without_coverage() {
    let dir = TempDir::new().unwrap();
    let (refa, sam, bed_path) = intervals_without_coverage_inputs(&dir);

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let values = reader.values("chr1", 0, 40).unwrap();
    assert!(all_missing(&values[..5]));
    assert_eq!(values[5..15], repeat(1.0, 10));
    assert!(all_missing(&values[15..25]));
    assert_eq!(values[25..35], repeat(0.0, 10));
    assert!(all_missing(&values[35..]));
}

/// With `--intervals`, the intervals of a contig no read reaches have zero depth at every position.
#[test]
fn per_base_coverage_with_intervals_writes_zeros_on_a_contig_without_reads() {
    let dir = TempDir::new().unwrap();
    let (refa, sam, bed_path) = intervals_without_coverage_inputs(&dir);

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let values = reader.values("chr2", 0, 20).unwrap();
    assert!(all_missing(&values[..5]));
    assert_eq!(values[5..15], repeat(0.0, 10));
    assert!(all_missing(&values[15..]));
}

/// With `--intervals`, a BAM with no reads still gets zero depth at every interval position.
#[test]
fn per_base_coverage_with_intervals_writes_zeros_on_every_interval_of_a_bam_without_reads() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let sam = coord_builder(&[("chr1", 20), ("chr2", 20)]);
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("targets.bed");
    BedBuilder::new().interval("chr1", 5, 15).interval("chr2", 0, 10).write_bed(&bed_path).unwrap();

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let chr1 = reader.values("chr1", 0, 20).unwrap();
    assert!(all_missing(&chr1[..5]));
    assert_eq!(chr1[5..15], repeat(0.0, 10));
    assert!(all_missing(&chr1[15..]));
    let chr2 = reader.values("chr2", 0, 20).unwrap();
    assert_eq!(chr2[..10], repeat(0.0, 10));
    assert!(all_missing(&chr2[10..]));
}

/// With `--intervals`, an interval without coverage is stored as one zero-valued item, and an
/// interval with coverage as per-position data.
#[test]
fn per_base_coverage_with_intervals_stores_an_interval_without_coverage_as_one_zero_item() {
    let dir = TempDir::new().unwrap();
    let (refa, sam, bed_path) = intervals_without_coverage_inputs(&dir);

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let expected: Vec<(u32, u32, f32)> =
        (5..15).map(|pos| (pos, pos + 1, 1.0)).chain([(25, 35, 0.0)]).collect();
    assert_eq!(items(&mut reader, "chr1"), expected);
    assert_eq!(items(&mut reader, "chr2"), [(5, 15, 0.0)]);
}

/// The first zoom level is 1 kb even when a long contig without reads, stored as one item, comes
/// before the covered ones.
#[test]
fn per_base_coverage_zoom_levels_start_at_1_kb_when_leading_contigs_have_no_reads() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 1_000_000])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 1_000_000), ("chr2", 20)]);
    sam.add(pair("r0").at("chr2", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.zoom_levels()[0].reduction(), 1_000);
}

/// Contigs keep their BAM header order, with `chromId` equal to the BAM `tid`, even when that
/// order is not sorted by name.
#[test]
fn per_base_coverage_keeps_contigs_in_bam_header_order() {
    let refa = FastaBuilder::new()
        .contig("chr2", vec![b'A'; 20])
        .contig("chr10", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr2", 20), ("chr10", 20)]);
    sam.add(pair("r0").at("chr2", 1, 11).len(10));
    sam.add(pair("r1").at("chr10", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    let names: Vec<&str> = reader.chroms().iter().map(|chrom| chrom.name.as_str()).collect();
    assert_eq!(names, ["chr2", "chr10"]);
    assert_eq!(reader.values("chr10", 0, 20).unwrap(), repeat(1.0, 20));
}

/// With intervals, each contig's data follows its own intervals even when the BAM lists contigs in
/// a different order from the FASTA.
#[test]
fn per_base_coverage_with_intervals_follows_the_bam_contig_order() {
    let dir = TempDir::new().unwrap();
    let (refa, sam, bed_path) = swapped_contig_order_inputs(&dir);

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(1.0, 20));
    assert!(all_missing(&reader.values("chr2", 0, 20).unwrap()));
}

/// Length of the contig in the tests of contigs longer than one push to the writer (8,192
/// positions).
const LONG_CONTIG_LEN: u32 = 20_000;

/// Pairs of 100 bp reads on a `LONG_CONTIG_LEN` chr1, straddling the push boundaries at 8,192 and
/// 16,384 and covering the contig's last base, with the per-base depth they give.
fn long_contig_pairs() -> (SamBuilder, Vec<f32>) {
    let mut sam = coord_builder(&[("chr1", LONG_CONTIG_LEN as usize)]);
    let mut expected = vec![0.0; LONG_CONTIG_LEN as usize];
    // 1-based mate starts; the last pair's mates overlap completely and count once.
    for (r1_start, r2_start) in [(8_150, 8_260), (16_340, 16_450), (19_901, 19_901)] {
        sam.add(pair(sam.next_id()).at("chr1", r1_start, r2_start).len(100));
        for start in [r1_start, r2_start] {
            expected[start - 1..start + 99].fill(1.0);
        }
    }
    (sam, expected)
}

/// A contig longer than one push to the writer is written as one continuous run.
#[test]
fn per_base_coverage_writes_contigs_longer_than_one_push() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; LONG_CONTIG_LEN as usize])
        .to_temp_fasta()
        .unwrap();
    let (sam, expected) = long_contig_pairs();
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, LONG_CONTIG_LEN).unwrap(), expected);
}

/// An interval that starts mid-contig and spans several pushes to the writer gets data at exactly
/// its positions.
#[test]
fn per_base_coverage_with_intervals_writes_intervals_longer_than_one_push() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; LONG_CONTIG_LEN as usize])
        .to_temp_fasta()
        .unwrap();
    let (sam, expected) = long_contig_pairs();
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("targets.bed");
    BedBuilder::new().interval("chr1", 5_000, 18_000).write_bed(&bed_path).unwrap();

    let mut reader = run_per_base(&sam, refa.path(), Some(bed_path), 250, &dir);

    let values = reader.values("chr1", 0, LONG_CONTIG_LEN).unwrap();
    assert!(all_missing(&values[..5_000]));
    assert_eq!(values[5_000..18_000], expected[5_000..18_000]);
    assert!(all_missing(&values[18_000..]));
}

/// Each contig's depth starts from zero, including a contig that follows a shorter one and is as
/// long as an earlier covered one.
#[test]
fn per_base_coverage_starts_each_contig_from_zero() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 40])
        .contig("chr2", vec![b'A'; 10])
        .contig("chr3", vec![b'A'; 40])
        .to_temp_fasta()
        .unwrap();
    let mut sam = coord_builder(&[("chr1", 40), ("chr2", 10), ("chr3", 40)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    sam.add(pair("r1").at("chr2", 1, 1).len(10));
    sam.add(pair("r2").at("chr3", 21, 31).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    let chr1: Vec<f32> = [repeat(1.0, 20), repeat(0.0, 20)].concat();
    let chr3: Vec<f32> = [repeat(0.0, 20), repeat(1.0, 20)].concat();
    assert_eq!(reader.values("chr1", 0, 40).unwrap(), chr1);
    assert_eq!(reader.values("chr2", 0, 10).unwrap(), repeat(1.0, 10));
    assert_eq!(reader.values("chr3", 0, 40).unwrap(), chr3);
}

/// A run that fails after the per-base file is opened removes it, leaving no partial output.
#[test]
fn per_base_coverage_is_removed_when_the_run_fails() {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    // Declared coordinate-sorted, but returns to chr1 after chr1 has been written.
    let mut sam = two_contig_builder().declare_coordinate_sorted();
    sam.add(read().name("chr1_first").at("chr1", 1).len(10));
    sam.add(read().name("chr2_between").at("chr2", 1).len(10));
    sam.add(read().name("chr1_again").at("chr1", 1).len(10));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let mut cmd =
        make_cmd(bam.path(), refa.path(), &dir.path().join("out"), None, false, true, 0, 0, 250);
    cmd.options.per_base_coverage = Some(PerBaseCoverageFormat::Bigwig);

    cmd.execute(None).expect_err("wgs accepted out-of-order input");

    assert!(!dir.path().join(format!("out{PER_BASE_BIGWIG_SUFFIX}")).exists());
}

/// A BAM whose header lists no contigs has no per-base data, so no per-base file is written.
#[test]
fn per_base_coverage_is_skipped_when_the_bam_lists_no_contigs() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let bam = coord_builder(&[]).to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let mut cmd =
        make_cmd(bam.path(), refa.path(), &dir.path().join("out"), None, false, true, 20, 20, 250);
    cmd.options.per_base_coverage = Some(PerBaseCoverageFormat::Bigwig);

    cmd.execute(None).unwrap();

    assert!(dir.path().join(format!("out{METRICS_SUFFIX}")).exists());
    assert!(!dir.path().join(format!("out{PER_BASE_BIGWIG_SUFFIX}")).exists());
}

/// The per-base file is finished, and readable, even when the run has no genome territory.
#[test]
fn per_base_coverage_is_finished_when_genome_territory_is_empty() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'N'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    let dir = TempDir::new().unwrap();

    let mut reader = run_per_base(&sam, refa.path(), None, 250, &dir);

    assert_eq!(reader.values("chr1", 0, 20).unwrap(), repeat(1.0, 20));
}

/// Turning on per-base output changes none of the other outputs.
#[test]
fn per_base_coverage_leaves_the_metrics_unchanged() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    sam.add(pair("overlapper").at("chr1", 1, 6).len(10));
    let bam = sam.to_temp_bam().unwrap();
    let without_dir = TempDir::new().unwrap();
    let with_dir = TempDir::new().unwrap();

    make_cmd(
        bam.path(),
        refa.path(),
        &without_dir.path().join("out"),
        None,
        false,
        true,
        20,
        20,
        3,
    )
    .execute(None)
    .unwrap();
    let mut with_cmd = make_cmd(
        bam.path(),
        refa.path(),
        &with_dir.path().join("out"),
        None,
        false,
        true,
        20,
        20,
        3,
    );
    with_cmd.options.per_base_coverage = Some(PerBaseCoverageFormat::Bigwig);
    with_cmd.execute(None).unwrap();

    for suffix in [METRICS_SUFFIX, COVERAGE_SUFFIX] {
        let without = std::fs::read(without_dir.path().join(format!("out{suffix}"))).unwrap();
        let with = std::fs::read(with_dir.path().join(format!("out{suffix}"))).unwrap();
        assert_eq!(with, without, "{suffix} differs with per-base output on");
    }
}

/// Without `--per-base-coverage`, no per-base file is written.
#[test]
fn per_base_coverage_is_off_by_default() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    sam.add(pair("r0").at("chr1", 1, 11).len(10));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();

    make_cmd(bam.path(), refa.path(), &dir.path().join("out"), None, false, true, 20, 20, 250)
        .execute(None)
        .unwrap();

    assert!(!dir.path().join(format!("out{PER_BASE_BIGWIG_SUFFIX}")).exists());
}

/// Every base above `--coverage-cap` is a capped exclusion: 5 pairs over 20 positions with a cap
/// of 3 tally 3 × 20 = 60 bases and exclude 2 × 20 = 40, so 40 / 100.
#[test]
fn capped_exclusions_count_every_base_above_the_cap() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();

    make_cmd(bam.path(), refa.path(), &dir.path().join("out"), None, false, true, 20, 20, 3)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_close!(rows[0].frac_excluded_capped, 0.4, 1e-5);
}

/// Bases above the cap at N positions are capped exclusions too, though N positions are outside
/// the territory: 10 non-N positions tally 3 × 10 = 30 bases, and all 20 positions exclude 2
/// each, so 40 / 70.
#[test]
fn capped_exclusions_include_bases_at_n_positions() {
    let seq: Vec<u8> = std::iter::repeat_n(b'A', 10).chain(std::iter::repeat_n(b'N', 10)).collect();
    let refa = FastaBuilder::new().contig("chr1", seq).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();

    make_cmd(bam.path(), refa.path(), &dir.path().join("out"), None, false, true, 20, 20, 3)
        .execute(None)
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_close!(rows[0].frac_excluded_capped, 40.0 / 70.0, 1e-5);
}

/// Bases above the cap outside the intervals are capped exclusions too: the 10 interval
/// positions tally 3 × 10 = 30 bases, and all 20 positions exclude 2 each, so 40 / 70.
#[test]
fn capped_exclusions_include_bases_outside_the_intervals() {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 20]).to_temp_fasta().unwrap();
    let mut sam = coord_builder(&[("chr1", 20)]);
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 1, 11).len(10));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("targets.bed");
    std::fs::write(&bed_path, "chr1\t5\t15\n").unwrap();

    make_cmd(
        bam.path(),
        refa.path(),
        &dir.path().join("out"),
        Some(bed_path),
        false,
        true,
        20,
        20,
        3,
    )
    .execute(None)
    .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows[0].genome_territory, 10);
    assert_close!(rows[0].frac_excluded_capped, 40.0 / 70.0, 1e-5);
}

// ─── Input requirements ───────────────────────────────────────────────────────

/// A builder over two 20 bp contigs, `chr1` and `chr2`, in that header order.
fn two_contig_builder() -> SamBuilder {
    SamBuilder::with_contigs(&[("chr1".to_string(), 20), ("chr2".to_string(), 20)])
}

/// Run `wgs` on input that is not in coordinate order, check that it wrote no output,
/// and return its error.
fn expect_order_rejection(sam: &SamBuilder) -> anyhow::Error {
    let refa = FastaBuilder::new()
        .contig("chr1", vec![b'A'; 20])
        .contig("chr2", vec![b'A'; 20])
        .to_temp_fasta()
        .unwrap();
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    let err = make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 0, 0, 250)
        .execute(None)
        .expect_err("wgs accepted out-of-order input");
    let metrics_path = dir.path().join(format!("out{METRICS_SUFFIX}"));
    assert!(!metrics_path.exists(), "rejected run wrote {}", metrics_path.display());
    err
}

/// wgs sweeps the reference one contig at a time, so input whose header does not
/// declare coordinate order is refused.
#[test]
fn test_unsorted_input_is_rejected() {
    let mut sam = two_contig_builder().sort_order(SortOrder::Unsorted);
    sam.add(read().name("r1").at("chr1", 1).len(10));
    let err = expect_order_rejection(&sam);
    assert!(err.to_string().contains("@HD SO:coordinate"), "unexpected error: {err}");
}

/// A header that declares `@HD SO:coordinate` over records that return to an earlier
/// contig passes the header check, so the runtime guard must refuse it.
#[test]
fn test_mislabeled_contig_revisit_is_rejected() {
    let mut sam = two_contig_builder().declare_coordinate_sorted();
    sam.add(read().name("chr1_first").at("chr1", 1).len(10));
    sam.add(read().name("chr2_between").at("chr2", 1).len(10));
    sam.add(read().name("chr1_again").at("chr1", 1).len(10));
    let err = expect_order_rejection(&sam);
    assert!(err.to_string().contains("chr1:1 after chr2:1"), "unexpected error: {err}");
}

/// Mates of an overlapping pair written right-mate first under a header that declares
/// `@HD SO:coordinate`: the runtime guard must refuse the position regression, which
/// would otherwise double-count the overlap.
#[test]
fn test_mislabeled_out_of_order_positions_are_rejected() {
    // `build()` resolves contigs against the default reference, whose `chr1` is also
    // the first contig in this header.
    let (r1, r2) = pair("p").at("chr1", 1, 6).len(10).build();
    let mut sam = two_contig_builder().declare_coordinate_sorted();
    sam.add(r2);
    sam.add(r1);
    let err = expect_order_rejection(&sam);
    assert!(err.to_string().contains("chr1:1 after chr1:6"), "unexpected error: {err}");
}

/// A position regression on a record the filters drop (here a duplicate, excluded by
/// default) must still be refused: the guard checks order before filtering.
#[test]
fn test_mislabeled_order_on_filtered_record_is_rejected() {
    let mut sam = two_contig_builder().declare_coordinate_sorted();
    sam.add(read().name("kept_first").at("chr1", 11).len(5));
    sam.add(read().name("filtered_dup").at("chr1", 1).len(5).duplicate());
    sam.add(read().name("kept_second").at("chr1", 15).len(5));
    let err = expect_order_rejection(&sam);
    assert!(err.to_string().contains("chr1:1 after chr1:11"), "unexpected error: {err}");
}
