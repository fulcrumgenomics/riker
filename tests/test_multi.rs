mod helpers;

use anyhow::Result;
use helpers::{FastaBuilder, SamBuilder, SortOrder, read_metrics_tsv};
use riker_lib::commands::alignment::{
    AlignmentSummaryMetric, METRICS_SUFFIX as ALIGNMENT_SUFFIX, MultiAlignmentOptions,
};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{
    InputOptions, OptionalReferenceOptions, OutputOptions, ReferenceOptions,
};
use riker_lib::commands::error::MultiErrorOptions;
use riker_lib::commands::gcbias::{
    DETAIL_SUFFIX as GCBIAS_DETAIL_SUFFIX, MultiGcBiasOptions,
    SUMMARY_SUFFIX as GCBIAS_SUMMARY_SUFFIX,
};
use riker_lib::commands::hybcap::{
    HybCapMetric, METRICS_SUFFIX as HYBCAP_SUFFIX, MultiHybCapOptions,
};
use riker_lib::commands::isize::{
    InsertSizeMetric, METRICS_SUFFIX as ISIZE_SUFFIX, MultiIsizeOptions,
};
use riker_lib::commands::multi::{CollectorKind, Multi};
use riker_lib::commands::rna::{METRICS_SUFFIX as RNA_METRICS_SUFFIX, RnaOptions, RnaSeqMetric};
use riker_lib::commands::wgs::{
    COVERAGE_SUFFIX as WGS_COVERAGE_SUFFIX, METRICS_SUFFIX as WGS_SUFFIX, MultiWgsOptions, Wgs,
    WgsMetrics, WgsOptions,
};
use std::path::PathBuf;
use tempfile::TempDir;

// ─── Helper ───────────────────────────────────────────────────────────────────

/// Create default per-tool options for the multi command.
///
/// For hybcap, the required paths (baits/targets) are set to `None` so that
/// tests that don't use hybcap don't accidentally pass validation.
fn default_per_tool_opts() -> (
    MultiWgsOptions,
    MultiIsizeOptions,
    MultiHybCapOptions,
    MultiGcBiasOptions,
    MultiAlignmentOptions,
    MultiErrorOptions,
) {
    // For hybcap, we need to use From but then clear the required paths back to None,
    // since HybCapOptions::default() has empty PathBuf values that would pass validate().
    let mut hybcap_opts: MultiHybCapOptions =
        riker_lib::commands::hybcap::HybCapOptions::default().into();
    hybcap_opts.hybcap_baits = None;
    hybcap_opts.hybcap_targets = None;

    // For error, clear the required reference path back to None.
    let mut error_opts: MultiErrorOptions =
        riker_lib::commands::error::ErrorOptions::default().into();
    error_opts.error_reference = None;

    (
        WgsOptions::default().into(),
        riker_lib::commands::isize::IsizeOptions::default().into(),
        hybcap_opts,
        riker_lib::commands::gcbias::GcBiasOptions::default().into(),
        riker_lib::commands::alignment::AlignmentOptions::default().into(),
        error_opts,
    )
}

fn make_multi(
    bam_path: &std::path::Path,
    prefix: &std::path::Path,
    collectors: Vec<CollectorKind>,
) -> Multi {
    let (wgs_opts, isize_opts, hybcap_opts, gcbias_opts, alignment_opts, error_opts) =
        default_per_tool_opts();
    Multi {
        input: InputOptions { input: bam_path.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: OptionalReferenceOptions { reference: None },

        tools: collectors,
        wgs_opts,
        isize_opts,
        hybcap_opts,
        gcbias_opts,
        alignment_opts,
        error_opts,
        rna_opts: riker_lib::commands::rna::RnaOptions::default().into(),
    }
}

fn make_multi_with_ref(
    bam_path: &std::path::Path,
    prefix: &std::path::Path,
    ref_path: &std::path::Path,
    collectors: Vec<CollectorKind>,
) -> Multi {
    let (wgs_opts, isize_opts, hybcap_opts, gcbias_opts, alignment_opts, error_opts) =
        default_per_tool_opts();
    Multi {
        input: InputOptions { input: bam_path.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: OptionalReferenceOptions { reference: Some(ref_path.to_path_buf()) },

        tools: collectors,
        wgs_opts,
        isize_opts,
        hybcap_opts,
        gcbias_opts,
        alignment_opts,
        error_opts,
        rna_opts: riker_lib::commands::rna::RnaOptions::default().into(),
    }
}

fn build_test_bam() -> Result<tempfile::NamedTempFile> {
    let mut builder = SamBuilder::new();
    for i in 0..5 {
        builder.add_pair(
            &format!("read{i}"),
            0,
            100 + i * 200,
            300 + i * 200,
            200,
            60,
            100,
            false,
            false,
        );
    }
    builder.to_temp_bam()
}

/// Build a coordinate-sorted BAM with a known contig length for WGS tests.
fn build_wgs_test_bam(contig_len: usize) -> Result<tempfile::NamedTempFile> {
    let contigs = vec![("chr1".to_string(), contig_len)];
    let mut builder = SamBuilder::with_contigs(&contigs).sort_order(SortOrder::Coordinate);
    for i in 0..5 {
        builder.add_pair(&format!("r{i}"), 0, 1, 11, 20, 60, 10, false, false);
    }
    builder.to_temp_bam()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[test]
fn test_both_collectors() -> Result<()> {
    let bam = build_test_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi(bam.path(), &prefix, vec![CollectorKind::Isize, CollectorKind::Alignment])
        .execute(None)?;

    // Both output files should exist and contain data.
    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));

    assert!(isize_path.exists(), "isize metrics file should exist");
    assert!(alignment_path.exists(), "alignment metrics file should exist");

    let isize_metrics: Vec<InsertSizeMetric> = read_metrics_tsv(&isize_path)?;
    assert!(!isize_metrics.is_empty(), "isize metrics should have rows");
    assert_eq!(isize_metrics[0].pair_orientation, "FR");
    assert_eq!(isize_metrics[0].read_pairs, 5);

    let alignment_metrics: Vec<AlignmentSummaryMetric> = read_metrics_tsv(&alignment_path)?;
    assert!(!alignment_metrics.is_empty(), "alignment metrics should have rows");
    let pair = alignment_metrics.iter().find(|m| m.category == "pair").unwrap();
    assert_eq!(pair.total_reads, 10);

    Ok(())
}

#[test]
fn test_isize_only() -> Result<()> {
    let bam = build_test_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi(bam.path(), &prefix, vec![CollectorKind::Isize]).execute(None)?;

    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));

    assert!(isize_path.exists(), "isize metrics file should exist");
    assert!(!alignment_path.exists(), "alignment metrics file should NOT exist");

    let isize_metrics: Vec<InsertSizeMetric> = read_metrics_tsv(&isize_path)?;
    assert_eq!(isize_metrics[0].read_pairs, 5);

    Ok(())
}

#[test]
fn test_alignment_only() -> Result<()> {
    let bam = build_test_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi(bam.path(), &prefix, vec![CollectorKind::Alignment]).execute(None)?;

    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));

    assert!(!isize_path.exists(), "isize metrics file should NOT exist");
    assert!(alignment_path.exists(), "alignment metrics file should exist");

    let alignment_metrics: Vec<AlignmentSummaryMetric> = read_metrics_tsv(&alignment_path)?;
    let pair = alignment_metrics.iter().find(|m| m.category == "pair").unwrap();
    assert_eq!(pair.total_reads, 10);

    Ok(())
}

#[test]
fn test_rna_through_multi_writes_metrics() -> Result<()> {
    // rna is wired into multi (CollectorKind::Rna + MultiRnaOptions); exercise that path
    // end-to-end so an integration break in build_collectors/execute is caught. rna needs a
    // coordinate-sorted input and a gene model.
    let mut b = SamBuilder::with_contigs(&[("chr1".to_string(), 100_000)])
        .sort_order(SortOrder::Coordinate);
    b.add_unpaired("coding", 0, 1500, 60, 100, false, false, false, None);
    b.add_unpaired("utr", 0, 1000, 60, 100, false, false, false, None);
    b.add_unpaired("intron", 0, 2400, 60, 100, false, false, false, None);
    let bam = b.to_temp_bam()?;

    // Two-exon coding gene MYGENE on chr1 (refFlat, 0-based half-open).
    let gm = tempfile::NamedTempFile::new()?;
    std::fs::write(
        gm.path(),
        "MYGENE\tNM_1\tchr1\t+\t999\t4000\t1499\t3500\t2\t999,2999,\t2000,4000,\n",
    )?;

    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");
    let (wgs_opts, isize_opts, hybcap_opts, gcbias_opts, alignment_opts, error_opts) =
        default_per_tool_opts();
    let multi = Multi {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        reference: OptionalReferenceOptions { reference: None },
        tools: vec![CollectorKind::Rna],
        wgs_opts,
        isize_opts,
        hybcap_opts,
        gcbias_opts,
        alignment_opts,
        error_opts,
        rna_opts: RnaOptions { gene_model: gm.path().to_path_buf(), ..Default::default() }.into(),
    };
    multi.execute(None)?;

    let metrics: Vec<RnaSeqMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{RNA_METRICS_SUFFIX}", prefix.display())))?;
    assert_eq!(metrics.len(), 1, "rna via multi should write exactly one summary row");
    assert_eq!(metrics[0].mapped_reads, 3, "all three mapped reads should be counted");

    Ok(())
}

#[test]
fn test_matches_standalone_isize() -> Result<()> {
    use riker_lib::commands::isize::{InsertSize, IsizeOptions};

    let bam = build_test_bam()?;

    // Run multi with isize only.
    let dir_multi = TempDir::new()?;
    let prefix_multi = dir_multi.path().join("out");
    make_multi(bam.path(), &prefix_multi, vec![CollectorKind::Isize]).execute(None)?;

    // Run standalone isize with matching default options.
    let dir_standalone = TempDir::new()?;
    let prefix_standalone = dir_standalone.path().join("out");
    let standalone = InsertSize {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix_standalone.clone() },
        reference: OptionalReferenceOptions { reference: None },
        options: IsizeOptions::default(),
    };
    standalone.execute(None)?;

    let multi_metrics: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix_multi.display())))?;
    let standalone_metrics: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix_standalone.display())))?;

    assert_eq!(multi_metrics.len(), standalone_metrics.len());
    for (m, s) in multi_metrics.iter().zip(standalone_metrics.iter()) {
        assert_eq!(m.pair_orientation, s.pair_orientation);
        assert_eq!(m.read_pairs, s.read_pairs);
        assert!((m.mean_insert_size - s.mean_insert_size).abs() < 1e-6);
        assert!((m.median_insert_size - s.median_insert_size).abs() < 1e-6);
        assert!((m.standard_deviation - s.standard_deviation).abs() < 1e-6);
    }

    Ok(())
}

#[test]
fn test_empty_bam() -> Result<()> {
    let builder = SamBuilder::new();
    let bam = builder.to_temp_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi(bam.path(), &prefix, vec![CollectorKind::Isize, CollectorKind::Alignment])
        .execute(None)?;

    // Both files should exist (even if empty / with zero counts).
    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));

    assert!(isize_path.exists(), "isize metrics file should exist for empty BAM");
    assert!(alignment_path.exists(), "alignment metrics file should exist for empty BAM");

    Ok(())
}

// ─── WGS multi tests ─────────────────────────────────────────────────────────

#[test]
fn test_wgs_collector() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let bam = build_wgs_test_bam(20)?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi_with_ref(bam.path(), &prefix, refa.path(), vec![CollectorKind::Wgs])
        .execute(None)?;

    let wgs_path = PathBuf::from(format!("{}{WGS_SUFFIX}", prefix.display()));
    let cov_path = PathBuf::from(format!("{}{WGS_COVERAGE_SUFFIX}", prefix.display()));

    assert!(wgs_path.exists(), "WGS metrics file should exist");
    assert!(cov_path.exists(), "WGS coverage file should exist");

    let wgs_metrics: Vec<WgsMetrics> = read_metrics_tsv(&wgs_path)?;
    assert_eq!(wgs_metrics.len(), 1);
    assert_eq!(wgs_metrics[0].genome_territory, 20);

    Ok(())
}

#[test]
fn test_all_collectors() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let bam = build_wgs_test_bam(20)?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi_with_ref(
        bam.path(),
        &prefix,
        refa.path(),
        vec![CollectorKind::Alignment, CollectorKind::Isize, CollectorKind::Wgs],
    )
    .execute(None)?;

    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));
    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let wgs_path = PathBuf::from(format!("{}{WGS_SUFFIX}", prefix.display()));

    assert!(alignment_path.exists(), "alignment metrics file should exist");
    assert!(isize_path.exists(), "isize metrics file should exist");
    assert!(wgs_path.exists(), "WGS metrics file should exist");

    Ok(())
}

#[test]
fn test_matches_standalone_wgs() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let bam = build_wgs_test_bam(20)?;

    // Run multi with WGS only.
    let dir_multi = TempDir::new()?;
    let prefix_multi = dir_multi.path().join("out");
    make_multi_with_ref(bam.path(), &prefix_multi, refa.path(), vec![CollectorKind::Wgs])
        .execute(None)?;

    // Run standalone WGS with matching default options.
    let dir_standalone = TempDir::new()?;
    let prefix_standalone = dir_standalone.path().join("out");
    let standalone = Wgs {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix_standalone.clone() },
        reference: ReferenceOptions { reference: refa.path().to_path_buf() },

        options: WgsOptions::default(),
    };
    standalone.execute(None)?;

    let multi_metrics: Vec<WgsMetrics> =
        read_metrics_tsv(&PathBuf::from(format!("{}{WGS_SUFFIX}", prefix_multi.display())))?;
    let standalone_metrics: Vec<WgsMetrics> =
        read_metrics_tsv(&PathBuf::from(format!("{}{WGS_SUFFIX}", prefix_standalone.display())))?;

    assert_eq!(multi_metrics.len(), standalone_metrics.len());
    assert_eq!(multi_metrics.len(), 1);
    let m = &multi_metrics[0];
    let s = &standalone_metrics[0];
    assert_eq!(m.genome_territory, s.genome_territory);
    assert!((m.mean_coverage - s.mean_coverage).abs() < 1e-6);
    assert!((m.median_coverage - s.median_coverage).abs() < 1e-6);
    assert!((m.sd_coverage - s.sd_coverage).abs() < 1e-6);
    assert!((m.frac_bases_at_1x - s.frac_bases_at_1x).abs() < 1e-6);
    assert!((m.frac_bases_at_5x - s.frac_bases_at_5x).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_wgs_requires_reference() {
    let bam = build_wgs_test_bam(20).unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // Multi with WGS but no reference should fail.
    let result = make_multi(bam.path(), &prefix, vec![CollectorKind::Wgs]).execute(None);
    assert!(result.is_err(), "WGS without reference should error");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("reference"),
        "error message should mention 'reference', got: {err_msg}"
    );
}

// ─── GC bias through multi ───────────────────────────────────────────────────

#[test]
fn test_gcbias_collector() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let bam = build_wgs_test_bam(20)?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi_with_ref(bam.path(), &prefix, refa.path(), vec![CollectorKind::GcBias])
        .execute(None)?;

    let detail_path = PathBuf::from(format!("{}{GCBIAS_DETAIL_SUFFIX}", prefix.display()));
    let summary_path = PathBuf::from(format!("{}{GCBIAS_SUMMARY_SUFFIX}", prefix.display()));

    assert!(detail_path.exists(), "gcbias detail file should exist");
    assert!(summary_path.exists(), "gcbias summary file should exist");

    Ok(())
}

#[test]
fn test_gcbias_requires_reference() {
    let bam = build_wgs_test_bam(20).unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    let result = make_multi(bam.path(), &prefix, vec![CollectorKind::GcBias]).execute(None);
    assert!(result.is_err(), "GC bias without reference should error");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("reference"),
        "error message should mention 'reference', got: {err_msg}"
    );
}

// ─── Deduplication ───────────────────────────────────────────────────────────

#[test]
fn test_duplicate_collectors_deduplicated() -> Result<()> {
    let bam = build_test_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    // Pass alignment twice — should be deduplicated to one.
    make_multi(
        bam.path(),
        &prefix,
        vec![CollectorKind::Alignment, CollectorKind::Alignment, CollectorKind::Isize],
    )
    .execute(None)?;

    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));
    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));

    assert!(alignment_path.exists(), "alignment metrics should exist");
    assert!(isize_path.exists(), "isize metrics should exist");

    // Verify we got exactly one set of alignment metrics (not doubled).
    let metrics: Vec<AlignmentSummaryMetric> = read_metrics_tsv(&alignment_path)?;
    let pair = metrics.iter().find(|m| m.category == "pair").unwrap();
    assert_eq!(pair.total_reads, 10);

    Ok(())
}

// ─── Parallel (threaded) tests ──────────────────────────────────────────────
//
// The thread count is supplied through `.execute(Some(n))` (the toolkit-wide
// `--threads`); `make_multi`/`make_multi_with_ref` build the unthreaded config.

#[test]
fn test_parallel_matches_single_threaded_isize_alignment() -> Result<()> {
    let bam = build_test_bam()?;
    let collectors = vec![CollectorKind::Isize, CollectorKind::Alignment];

    // Single-threaded run.
    let single_dir = TempDir::new()?;
    let single_prefix = single_dir.path().join("out");
    make_multi(bam.path(), &single_prefix, collectors.clone()).execute(None)?;

    // Parallel run with 2 threads.
    let parallel_dir = TempDir::new()?;
    let parallel_prefix = parallel_dir.path().join("out");
    make_multi(bam.path(), &parallel_prefix, collectors).execute(Some(2))?;

    // Compare isize metrics.
    let single_isize: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", single_prefix.display())))?;
    let parallel_isize: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", parallel_prefix.display())))?;
    assert_eq!(single_isize.len(), parallel_isize.len());
    for (s, p) in single_isize.iter().zip(parallel_isize.iter()) {
        assert_eq!(s.pair_orientation, p.pair_orientation);
        assert_eq!(s.read_pairs, p.read_pairs);
        assert!((s.mean_insert_size - p.mean_insert_size).abs() < 1e-6);
    }

    // Compare alignment metrics.
    let single_align: Vec<AlignmentSummaryMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", single_prefix.display())))?;
    let parallel_align: Vec<AlignmentSummaryMetric> = read_metrics_tsv(&PathBuf::from(format!(
        "{}{ALIGNMENT_SUFFIX}",
        parallel_prefix.display()
    )))?;
    assert_eq!(single_align.len(), parallel_align.len());
    for (s, p) in single_align.iter().zip(parallel_align.iter()) {
        assert_eq!(s.category, p.category);
        assert_eq!(s.total_reads, p.total_reads);
    }

    Ok(())
}

/// multi output must be byte-identical across thread counts when the budget is
/// set via the top-level `--threads` (the non-deprecated path). Compares every
/// `.txt` output file against the single-threaded run; PDF plots embed
/// non-deterministic bytes and are skipped.
#[test]
fn test_output_byte_identical_across_global_thread_counts() -> Result<()> {
    let bam = build_test_bam()?;
    let collectors = vec![CollectorKind::Alignment, CollectorKind::Basic, CollectorKind::Isize];

    // Baseline: single-threaded via the top-level --threads path.
    let base_dir = TempDir::new()?;
    make_multi(bam.path(), &base_dir.path().join("out"), collectors.clone()).execute(Some(1))?;

    let mut baseline: Vec<(std::ffi::OsString, Vec<u8>)> = Vec::new();
    for entry in std::fs::read_dir(base_dir.path())? {
        let path = entry?.path();
        if path.extension().is_some_and(|e| e == "txt") {
            let name = path.file_name().expect("output file has a name").to_owned();
            baseline.push((name, std::fs::read(&path)?));
        }
    }
    assert!(!baseline.is_empty(), "expected at least one .txt output file");

    for threads in [2u8, 3, 4] {
        let dir = TempDir::new()?;
        make_multi(bam.path(), &dir.path().join("out"), collectors.clone())
            .execute(Some(threads))?;
        for (name, bytes) in &baseline {
            let got = std::fs::read(dir.path().join(name))?;
            assert_eq!(&got, bytes, "output {name:?} differs at --threads {threads}");
        }
    }
    Ok(())
}

#[test]
fn test_parallel_all_collectors() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let bam = build_wgs_test_bam(20)?;

    // Single-threaded.
    let single_dir = TempDir::new()?;
    let single_prefix = single_dir.path().join("out");
    let collectors = vec![
        CollectorKind::Alignment,
        CollectorKind::Isize,
        CollectorKind::GcBias,
        CollectorKind::Wgs,
    ];
    make_multi_with_ref(bam.path(), &single_prefix, refa.path(), collectors.clone())
        .execute(None)?;

    // Parallel with 3 threads.
    let parallel_dir = TempDir::new()?;
    let parallel_prefix = parallel_dir.path().join("out");
    make_multi_with_ref(bam.path(), &parallel_prefix, refa.path(), collectors).execute(Some(3))?;

    // Compare WGS metrics.
    let single_wgs: Vec<WgsMetrics> =
        read_metrics_tsv(&PathBuf::from(format!("{}{WGS_SUFFIX}", single_prefix.display())))?;
    let parallel_wgs: Vec<WgsMetrics> =
        read_metrics_tsv(&PathBuf::from(format!("{}{WGS_SUFFIX}", parallel_prefix.display())))?;
    assert_eq!(single_wgs.len(), parallel_wgs.len());
    assert_eq!(single_wgs[0].genome_territory, parallel_wgs[0].genome_territory);
    assert!((single_wgs[0].mean_coverage - parallel_wgs[0].mean_coverage).abs() < 1e-6);

    // Compare alignment metrics.
    let single_align: Vec<AlignmentSummaryMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", single_prefix.display())))?;
    let parallel_align: Vec<AlignmentSummaryMetric> = read_metrics_tsv(&PathBuf::from(format!(
        "{}{ALIGNMENT_SUFFIX}",
        parallel_prefix.display()
    )))?;
    assert_eq!(single_align.len(), parallel_align.len());
    for (s, p) in single_align.iter().zip(parallel_align.iter()) {
        assert_eq!(s.category, p.category);
        assert_eq!(s.total_reads, p.total_reads);
    }

    Ok(())
}

/// Run with a thread budget larger than the collector count, across enough
/// records to produce multiple batches. With a single collector only one
/// worker group is created regardless of the budget; the surplus threads go
/// to BGZF decode (per `plan_multi`). Output must still match the serial path.
#[test]
fn test_parallel_more_threads_than_collectors() -> Result<()> {
    // BATCH_SIZE is 128 inside multi.rs — build enough pairs to produce
    // several batches so the worker drains across multiple channel fills.
    let mut builder = SamBuilder::new();
    for i in 0..1000 {
        builder.add_pair(
            &format!("read{i}"),
            0,
            100 + i * 50,
            300 + i * 50,
            200,
            60,
            100,
            false,
            false,
        );
    }
    let bam = builder.to_temp_bam()?;

    // Single-threaded baseline.
    let single_dir = TempDir::new()?;
    let single_prefix = single_dir.path().join("out");
    make_multi(bam.path(), &single_prefix, vec![CollectorKind::Isize]).execute(None)?;

    // Parallel with --threads 4 but only 1 collector: partition_collectors
    // makes a single worker group, and plan_multi hands the surplus budget to
    // BGZF decode threads. Output must still match the serial path.
    let parallel_dir = TempDir::new()?;
    let parallel_prefix = parallel_dir.path().join("out");
    make_multi(bam.path(), &parallel_prefix, vec![CollectorKind::Isize]).execute(Some(4))?;

    let single: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", single_prefix.display())))?;
    let parallel: Vec<InsertSizeMetric> =
        read_metrics_tsv(&PathBuf::from(format!("{}{ISIZE_SUFFIX}", parallel_prefix.display())))?;
    assert_eq!(single.len(), parallel.len());
    assert!(!single.is_empty());
    for (s, p) in single.iter().zip(parallel.iter()) {
        assert_eq!(s.pair_orientation, p.pair_orientation);
        assert_eq!(s.read_pairs, p.read_pairs);
        assert!((s.mean_insert_size - p.mean_insert_size).abs() < 1e-6);
    }
    Ok(())
}

#[test]
fn test_parallel_empty_bam() -> Result<()> {
    let builder = SamBuilder::new();
    let bam = builder.to_temp_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    make_multi(bam.path(), &prefix, vec![CollectorKind::Isize, CollectorKind::Alignment])
        .execute(Some(2))?;

    let isize_path = PathBuf::from(format!("{}{ISIZE_SUFFIX}", prefix.display()));
    let alignment_path = PathBuf::from(format!("{}{ALIGNMENT_SUFFIX}", prefix.display()));

    assert!(isize_path.exists(), "isize metrics file should exist for empty BAM");
    assert!(alignment_path.exists(), "alignment metrics file should exist for empty BAM");

    Ok(())
}

// ─── HybCap through multi ────────────────────────────────────────────────────

#[test]
fn test_hybcap_requires_targets_and_baits() {
    let bam = build_test_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // HybCap selected but no targets/baits → error.
    let result = make_multi(bam.path(), &prefix, vec![CollectorKind::HybCap]).execute(None);
    assert!(result.is_err(), "hybcap without targets should error");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("hybcap::baits") || err_msg.contains("hybcap::targets"),
        "error should mention required hybcap paths, got: {err_msg}"
    );
}

#[test]
fn test_hybcap_targets_not_required_when_not_selected() -> Result<()> {
    let bam = build_test_bam()?;
    let dir = TempDir::new()?;
    let prefix = dir.path().join("out");

    // Only isize selected — hybcap targets/baits should not be required.
    make_multi(bam.path(), &prefix, vec![CollectorKind::Isize]).execute(None)?;
    Ok(())
}

#[test]
fn test_hybcap_via_multi() -> Result<()> {
    use helpers::coord_builder;
    use std::io::Write;

    let dir = TempDir::new()?;
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");

    // Create BED files: one 100bp target/bait on chr1:100-200.
    {
        let mut f = std::fs::File::create(&bait_path)?;
        writeln!(f, "chr1\t100\t200")?;
    }
    {
        let mut f = std::fs::File::create(&target_path)?;
        writeln!(f, "chr1\t100\t200")?;
    }

    // Build a BAM with reads overlapping the target.
    let mut bld = coord_builder(&[("chr1", 10_000)]);
    for i in 0..5 {
        bld.add_pair(&format!("r{i}"), 0, 101 + i * 10, 151 + i * 10, 50, 60, 50, false, false);
    }
    let bam = bld.to_temp_bam()?;

    let prefix = dir.path().join("out");
    let (wgs_opts, isize_opts, _, gcbias_opts, alignment_opts, error_opts) =
        default_per_tool_opts();

    // Build hybcap multi opts with the required baits/targets paths.
    let hybcap_opts: MultiHybCapOptions = riker_lib::commands::hybcap::HybCapOptions {
        baits: bait_path,
        targets: target_path,
        ..riker_lib::commands::hybcap::HybCapOptions::default()
    }
    .into();

    let multi = Multi {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        reference: OptionalReferenceOptions { reference: None },

        tools: vec![CollectorKind::HybCap],
        wgs_opts,
        isize_opts,
        hybcap_opts,
        gcbias_opts,
        alignment_opts,
        error_opts,
        rna_opts: riker_lib::commands::rna::RnaOptions::default().into(),
    };
    multi.execute(None)?;

    let metrics_path = PathBuf::from(format!("{}{HYBCAP_SUFFIX}", prefix.display()));
    assert!(metrics_path.exists(), "hybcap metrics file should exist");

    let metrics: Vec<HybCapMetric> = read_metrics_tsv(&metrics_path)?;
    assert_eq!(metrics.len(), 1);
    assert!(metrics[0].total_reads > 0, "should have counted some reads");

    Ok(())
}

// ─── From conversion tests ──────────────────────────────────────────────────

#[test]
fn test_validate_roundtrip_preserves_defaults() {
    // Verify that converting default Options → Multi → validate → back to Options
    // preserves all default values.
    let wgs_default = WgsOptions::default();
    let multi_wgs: MultiWgsOptions = wgs_default.clone().into();
    let roundtripped: WgsOptions = multi_wgs.validate().expect("validate should succeed");
    assert_eq!(roundtripped.include_duplicates, wgs_default.include_duplicates);
    assert_eq!(roundtripped.min_mapq, wgs_default.min_mapq);
    assert_eq!(roundtripped.min_bq, wgs_default.min_bq);
    assert_eq!(roundtripped.coverage_cap, wgs_default.coverage_cap);
}
