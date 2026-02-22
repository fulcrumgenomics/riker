mod helpers;

use helpers::{SamBuilder, read_metrics_tsv};
use riker_lib::commands::command::Command;
use riker_lib::commands::isize::{
    HISTOGRAM_SUFFIX, InsertSize, InsertSizeHistogramEntry, InsertSizeMetric, IsizeOptions,
    METRICS_SUFFIX, PLOT_SUFFIX,
};
use tempfile::TempDir;

// ─── Helper ──────────────────────────────────────────────────────────────────

/// Build an `InsertSize` command pointing at `bam_path` with a given output prefix.
fn make_cmd(
    bam_path: &std::path::Path,
    prefix: &std::path::Path,
    include_duplicates: bool,
    min_frac: f64,
) -> InsertSize {
    use riker_lib::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
    InsertSize {
        input: InputOptions { input: bam_path.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: OptionalReferenceOptions { reference: None },
        options: IsizeOptions { include_duplicates, min_frac, deviations: 10.0 },
    }
}

/// Assert two f64 values are approximately equal.
macro_rules! assert_f64_eq {
    ($a:expr, $b:expr) => {{
        let a: f64 = $a;
        let b: f64 = $b;
        assert!((a - b).abs() < 1e-6, "expected {b} got {a} (diff={})", (a - b).abs());
    }};
}

// ─── Integration tests ────────────────────────────────────────────────────────

#[test]
fn test_basic_fr_five_pairs() {
    let mut builder = SamBuilder::new();
    // 5 FR pairs, all with insert size 100.
    for i in 0..5 {
        let name = format!("read{i}");
        // pos1=100, pos2=200, read_len=100 → FR orientation, insert_size=100
        builder.add_pair(&name, 0, 100, 200, 100, 60, 100, false, false);
    }
    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute().unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1, "expected one orientation row");
    let m = &rows[0];
    assert_eq!(m.pair_orientation, "FR");
    assert_eq!(m.read_pairs, 5);
    assert_f64_eq!(m.median_insert_size, 100.0);
    assert_eq!(m.mode_insert_size, 100);
    assert_eq!(m.min_insert_size, 100);
    assert_eq!(m.max_insert_size, 100);
    assert_f64_eq!(m.mean_insert_size, 100.0);
    assert_f64_eq!(m.standard_deviation, 0.0);
    assert_f64_eq!(m.median_absolute_deviation, 0.0);
}

#[test]
fn test_varied_insert_sizes() {
    // 6 FR pairs with distinct insert sizes.
    let sizes: [u32; 6] = [90, 95, 100, 100, 105, 110];
    let mut builder = SamBuilder::new();
    for (i, &sz) in sizes.iter().enumerate() {
        #[allow(clippy::cast_possible_wrap)]
        builder.add_pair(
            &format!("r{i}"),
            0,
            100,
            100 + sz as usize,
            sz as i32,
            60,
            1,
            false,
            false,
        );
    }
    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute().unwrap();

    let metrics_path = dir.path().join(format!("out{METRICS_SUFFIX}"));
    let rows: Vec<InsertSizeMetric> = read_metrics_tsv(&metrics_path).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];
    assert_eq!(m.pair_orientation, "FR");
    assert_eq!(m.read_pairs, 6);
    // median of {90,95,100,100,105,110} = (100+100)/2 = 100
    assert_f64_eq!(m.median_insert_size, 100.0);
    assert_eq!(m.mode_insert_size, 100);
    assert_eq!(m.min_insert_size, 90);
    assert_eq!(m.max_insert_size, 110);
    assert_f64_eq!(m.mean_insert_size, 100.0);
}

#[test]
fn test_duplicate_filter() {
    let mut builder = SamBuilder::new();
    // 5 normal pairs
    for i in 0..5 {
        builder.add_pair(&format!("norm{i}"), 0, 100, 200, 100, 60, 50, false, false);
    }
    // 3 duplicate pairs
    for i in 0..3 {
        builder.add_pair(&format!("dup{i}"), 0, 100, 200, 100, 60, 50, true, false);
    }
    let bam = builder.to_temp_bam().unwrap();

    let dir_excl = TempDir::new().unwrap();
    let prefix_excl = dir_excl.path().join("out");
    make_cmd(bam.path(), &prefix_excl, false, 0.0).execute().unwrap();
    let rows_excl: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir_excl.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows_excl[0].read_pairs, 5, "should exclude duplicates");

    let dir_incl = TempDir::new().unwrap();
    let prefix_incl = dir_incl.path().join("out");
    make_cmd(bam.path(), &prefix_incl, true, 0.0).execute().unwrap();
    let rows_incl: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir_incl.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows_incl[0].read_pairs, 8, "should include duplicates");
}

#[test]
fn test_qc_fail_filter() {
    let mut builder = SamBuilder::new();
    // 5 PF pairs + 2 QC-fail pairs
    for i in 0..5 {
        builder.add_pair(&format!("pf{i}"), 0, 100, 200, 100, 60, 50, false, false);
    }
    for i in 0..2 {
        builder.add_pair(&format!("fail{i}"), 0, 100, 200, 100, 60, 50, false, true);
    }
    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute().unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows[0].read_pairs, 5, "QC-fail reads should be excluded");
}

#[test]
fn test_min_frac_excludes_minority_orientation() {
    let mut builder = SamBuilder::new();
    // 20 FR pairs: pos1 < pos2
    for i in 0..20 {
        builder.add_pair(&format!("fr{i}"), 0, 100, 200, 100, 60, 50, false, false);
    }
    // 1 RF pair: pos1 > pos2 (read2 at pos2=100, mate at pos1=200)
    builder.add_pair("rf0", 0, 200, 100, 100, 60, 50, false, false);

    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // RF = 1/21 ≈ 4.8% < 5% → excluded
    make_cmd(bam.path(), &prefix, false, 0.05).execute().unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert!(rows.iter().all(|r| r.pair_orientation == "FR"), "RF should be excluded by min_frac");
    assert_eq!(rows.len(), 1);
}

#[test]
fn test_min_frac_zero_includes_all_orientations() {
    let mut builder = SamBuilder::new();
    for i in 0..10 {
        builder.add_pair(&format!("fr{i}"), 0, 100, 200, 100, 60, 50, false, false);
    }
    builder.add_pair("rf0", 0, 200, 100, 100, 60, 50, false, false);

    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute().unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let orientations: Vec<&str> = rows.iter().map(|r| r.pair_orientation.as_str()).collect();
    assert!(orientations.contains(&"FR"), "FR should be present");
    assert!(orientations.contains(&"RF"), "RF should be present");
}

#[test]
fn test_histogram_output() {
    let mut builder = SamBuilder::new();
    builder.add_pair("r100a", 0, 100, 200, 100, 60, 50, false, false);
    builder.add_pair("r100b", 0, 100, 200, 100, 60, 50, false, false);
    builder.add_pair("r101", 0, 100, 201, 101, 60, 50, false, false);
    builder.add_pair("r102", 0, 100, 202, 102, 60, 50, false, false);

    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute().unwrap();

    let hist: Vec<InsertSizeHistogramEntry> =
        read_metrics_tsv(&dir.path().join(format!("out{HISTOGRAM_SUFFIX}"))).unwrap();
    assert_eq!(hist.len(), 3, "three distinct insert sizes");

    let e100 = hist.iter().find(|e| e.insert_size == 100).unwrap();
    assert_eq!(e100.fr_count, 2);
    assert_eq!(e100.rf_count, 0);
    assert_eq!(e100.tandem_count, 0);

    let e101 = hist.iter().find(|e| e.insert_size == 101).unwrap();
    assert_eq!(e101.fr_count, 1);

    let e102 = hist.iter().find(|e| e.insert_size == 102).unwrap();
    assert_eq!(e102.fr_count, 1);
}

#[test]
fn test_plot_output_created() {
    let mut builder = SamBuilder::new();
    for i in 0..5 {
        builder.add_pair(&format!("r{i}"), 0, 100, 200, 100, 60, 50, false, false);
    }
    let bam = builder.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute().unwrap();

    let plot_path = dir.path().join(format!("out{PLOT_SUFFIX}"));
    assert!(plot_path.exists(), "PDF plot file should be created");
    assert!(plot_path.metadata().unwrap().len() > 0, "PDF plot file should be non-empty");
}
