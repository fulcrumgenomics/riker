mod helpers;

use helpers::{FastaBuilder, coord_builder, read_metrics_tsv};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::wgs::{
    COVERAGE_SUFFIX, METRICS_SUFFIX, Wgs, WgsCoverageEntry, WgsMetrics, WgsOptions,
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
    exclude_unpaired_reads: bool,
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
            exclude_unpaired_reads,
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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..5 {
        bld.add_pair(&format!("r{i}"), 0, 1, 11, 20, 60, 10, false, false);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];

    assert_eq!(m.genome_territory, 20);
    assert_float_eq!(m.mean_coverage, 5.0, 0.01);
    assert_float_eq!(m.frac_bases_at_1x, 1.0, 0.001);
    assert_float_eq!(m.frac_bases_at_5x, 1.0, 0.001);
    assert_float_eq!(m.frac_bases_at_10x, 0.0, 0.001);
}

/// N bases in the reference should be excluded from genome territory.
#[test]
fn test_n_bases_excluded() {
    let seq: Vec<u8> = std::iter::repeat_n(b'A', 10).chain(std::iter::repeat_n(b'N', 10)).collect();
    let refa = FastaBuilder::new().add_contig("chr1", &seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..3 {
        bld.add_pair(&format!("r{i}"), 0, 1, 11, 20, 60, 10, false, false);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    // r1 covers 0–9, r2 covers 5–14 → 5 positions of overlap (5–9)
    bld.add_pair("overlapper", 0, 1, 6, 15, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_pair("lowmq", 0, 1, 11, 20, 10, 10, false, false); // mapq=10 < 20
    bld.add_pair("highmq", 0, 1, 11, 20, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_mapq > 0.0,
        "expected mapq exclusions, got {}",
        rows[0].frac_excluded_mapq
    );
    assert_float_eq!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Duplicate reads are excluded when `include_duplicates`=false (the default).
#[test]
fn test_dup_exclusion() {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_pair("dup", 0, 1, 11, 20, 60, 10, true, false);
    bld.add_pair("normal", 0, 1, 11, 20, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_dupe > 0.0,
        "expected dupe exclusions, got {}",
        rows[0].frac_excluded_dupe
    );
    assert_float_eq!(rows[0].mean_coverage, 1.0, 0.01);
}

/// Unpaired reads are excluded when --exclude-unpaired-reads is set.
#[test]
fn test_exclude_unpaired() {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_unpaired("frag", 0, 1, 60, 20, false, false, false, None);
    bld.add_pair("pair", 0, 1, 11, 20, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, true, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_unpaired > 0.0,
        "expected unpaired exclusions, got {}",
        rows[0].frac_excluded_unpaired
    );
}

/// Depth exceeding `coverage_cap` produces `frac_excluded_capped` > 0.
#[test]
fn test_coverage_cap() {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..5 {
        bld.add_pair(&format!("r{i}"), 0, 1, 11, 20, 60, 10, false, false);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // cap=3 → 5 reads per position, 2 excluded per position
    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 3).execute().unwrap();

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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_pair("r0", 0, 1, 11, 20, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..15 {
        bld.add_unpaired(&format!("a{i}"), 0, 1, 60, 10, false, false, false, None);
    }
    for i in 0..5 {
        bld.add_unpaired(&format!("b{i}"), 0, 11, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let m = &rows[0];

    assert_float_eq!(m.frac_bases_at_1x, 1.0, 0.001);
    assert_float_eq!(m.frac_bases_at_5x, 1.0, 0.001);
    assert_float_eq!(m.frac_bases_at_10x, 0.5, 0.001); // only 10/20 positions
    assert_float_eq!(m.frac_bases_at_15x, 0.5, 0.001);
    assert_float_eq!(m.frac_bases_at_20x, 0.0, 0.001);
    // mean = (15*10 + 5*10) / 20 = 10
    assert_float_eq!(m.mean_coverage, 10.0, 0.01);
}

/// Uniform coverage → `fold_80_base_penalty` should equal or approach 1.0.
#[test]
fn test_fold_penalty() {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 10]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 10)]);
    for i in 0..10 {
        bld.add_unpaired(&format!("r{i}"), 0, 1, 60, 10, false, false, false, None);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let m = &rows[0];

    assert_float_eq!(m.mean_coverage, 10.0, 0.01);
    // Uniform coverage → fold_80 = mean / 20th-pct = 10/10 = 1.0
    assert!(
        m.fold_80_base_penalty >= 1.0,
        "fold_80 should be ≥ 1.0, got {}",
        m.fold_80_base_penalty
    );
    assert_float_eq!(m.fold_80_base_penalty, 1.0, 0.01);
}

/// Interval restriction: `genome_territory` should equal only the in-interval non-N bases.
#[test]
fn test_with_intervals() {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    bld.add_pair("r0", 0, 1, 11, 20, 60, 10, false, false);
    let bam = bld.to_temp_bam().unwrap();

    // BED interval: chr1 positions 5–14 (0-based half-open → 10 positions)
    let dir = TempDir::new().unwrap();
    let bed_path = dir.path().join("test.bed");
    std::fs::write(&bed_path, "chr1\t5\t15\n").unwrap();

    let prefix = dir.path().join("out");
    make_cmd(bam.path(), refa.path(), &prefix, Some(bed_path), false, false, 20, 20, 250)
        .execute()
        .unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].genome_territory, 10);
}

/// Bases with quality below `min_bq` should be counted as baseq-excluded.
#[test]
fn test_baseq_exclusion() {
    use noodles::core::Position;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 10]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 10)]);

    // Low base quality read (qual=5 < min_bq=20)
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let seq: Sequence = vec![b'A'; 10].into();
    let qual = QualityScores::from(vec![5u8; 10]);
    let record = RecordBuf::builder()
        .set_name("lowbq")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0usize)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(seq)
        .set_quality_scores(qual)
        .build();
    bld.add_record(record);

    // High quality read covering the same positions
    bld.add_unpaired("highbq", 0, 1, 60, 10, false, false, false, None);

    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 250).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].frac_excluded_baseq > 0.0,
        "expected baseq exclusions, got {}",
        rows[0].frac_excluded_baseq
    );
    assert_float_eq!(rows[0].mean_coverage, 1.0, 0.01);
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
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 20)]);
    for i in 0..8 {
        bld.add_pair(&format!("r{i}"), 0, 1, 6, 10, 60, 5, false, false);
    }
    let bam = bld.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), refa.path(), &prefix, None, false, false, 20, 20, 5).execute().unwrap();

    let rows: Vec<WgsMetrics> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];

    assert_eq!(m.genome_territory, 20);
    assert_float_eq!(m.mean_coverage, 2.5, 0.01);

    // sd from capped depths: sqrt(125/19) ≈ 2.5643
    let expected_sd = (125.0_f64 / 19.0).sqrt();
    assert_float_eq!(m.sd_coverage, expected_sd, 0.01);
}
