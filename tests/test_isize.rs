mod helpers;

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::isize::{
    HISTOGRAM_SUFFIX, InsertSize, InsertSizeHistogramEntry, InsertSizeMetric, IsizeOptions,
    METRICS_SUFFIX, PLOT_SUFFIX,
};
use riker_lib::test_support::{ReadBuilder, SamBuilder, pair};
use tempfile::TempDir;

// ─── Helpers ──────────────────────────────────────────────────────────────────

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

// ─── Integration tests ────────────────────────────────────────────────────────

#[test]
fn test_basic_fr_five_pairs() {
    let mut sam = SamBuilder::new();
    // 5 FR pairs, all with insert size 100.
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute(None).unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 1, "expected one orientation row");
    let m = &rows[0];
    assert_eq!(m.pair_orientation, "FR");
    assert_eq!(m.read_pairs, 5);
    assert_close!(m.median_insert_size, 100.0);
    assert_eq!(m.mode_insert_size, 100);
    assert_eq!(m.min_insert_size, 100);
    assert_eq!(m.max_insert_size, 100);
    assert_close!(m.mean_insert_size, 100.0);
    assert_close!(m.standard_deviation, 0.0);
    assert_close!(m.median_absolute_deviation, 0.0);
}

#[test]
fn test_varied_insert_sizes() {
    // 6 FR pairs with distinct insert sizes. The 90 and 95 bp inserts are shorter
    // than the 100 bp reads, so the pair builder soft-clips the adapter read-through.
    let sizes: [usize; 6] = [90, 95, 100, 100, 105, 110];
    let mut sam = SamBuilder::new();
    for &sz in &sizes {
        sam.add(pair(sam.next_id()).fr("chr1", 100..(100 + sz)));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute(None).unwrap();

    let metrics_path = dir.path().join(format!("out{METRICS_SUFFIX}"));
    let rows: Vec<InsertSizeMetric> = read_metrics_tsv(&metrics_path).unwrap();
    assert_eq!(rows.len(), 1);
    let m = &rows[0];
    assert_eq!(m.pair_orientation, "FR");
    assert_eq!(m.read_pairs, 6);
    // median of {90,95,100,100,105,110} = (100+100)/2 = 100
    assert_close!(m.median_insert_size, 100.0);
    assert_eq!(m.mode_insert_size, 100);
    assert_eq!(m.min_insert_size, 90);
    assert_eq!(m.max_insert_size, 110);
    assert_close!(m.mean_insert_size, 100.0);
}

#[test]
fn test_duplicate_filter() {
    let mut sam = SamBuilder::new();
    // 5 normal pairs
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    // 3 duplicate pairs
    for _ in 0..3 {
        sam.add(
            pair(sam.next_id())
                .fr("chr1", 100..200)
                .r1(ReadBuilder::duplicate)
                .r2(ReadBuilder::duplicate),
        );
    }
    let bam = sam.to_temp_bam().unwrap();

    let dir_excl = TempDir::new().unwrap();
    let prefix_excl = dir_excl.path().join("out");
    make_cmd(bam.path(), &prefix_excl, false, 0.0).execute(None).unwrap();
    let rows_excl: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir_excl.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows_excl[0].read_pairs, 5, "should exclude duplicates");

    let dir_incl = TempDir::new().unwrap();
    let prefix_incl = dir_incl.path().join("out");
    make_cmd(bam.path(), &prefix_incl, true, 0.0).execute(None).unwrap();
    let rows_incl: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir_incl.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows_incl[0].read_pairs, 8, "should include duplicates");
}

#[test]
fn test_qc_fail_filter() {
    let mut sam = SamBuilder::new();
    // 5 PF pairs + 2 QC-fail pairs
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    for _ in 0..2 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200).qc_fail());
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute(None).unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(rows[0].read_pairs, 5, "QC-fail reads should be excluded");
}

#[test]
fn test_min_frac_excludes_minority_orientation() {
    let mut sam = SamBuilder::new();
    // 20 FR pairs
    for _ in 0..20 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    // 1 RF pair: the reverse mate sits upstream of the forward mate ("outie").
    sam.add(pair("rf0").at("chr1", 200, 100));

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    // RF = 1/21 ≈ 4.8% < 5% → excluded
    make_cmd(bam.path(), &prefix, false, 0.05).execute(None).unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert!(rows.iter().all(|r| r.pair_orientation == "FR"), "RF should be excluded by min_frac");
    assert_eq!(rows.len(), 1);
}

#[test]
fn test_min_frac_zero_includes_all_orientations() {
    let mut sam = SamBuilder::new();
    for _ in 0..10 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    // 1 RF pair: the reverse mate sits upstream of the forward mate ("outie").
    sam.add(pair("rf0").at("chr1", 200, 100));

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute(None).unwrap();

    let rows: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let orientations: Vec<&str> = rows.iter().map(|r| r.pair_orientation.as_str()).collect();
    assert!(orientations.contains(&"FR"), "FR should be present");
    assert!(orientations.contains(&"RF"), "RF should be present");
}

#[test]
fn test_histogram_output() {
    let inserts: [(&str, usize); 4] =
        [("r100a", 100), ("r100b", 100), ("r101", 101), ("r102", 102)];
    let mut sam = SamBuilder::new();
    for (name, insert) in inserts {
        sam.add(pair(name).fr("chr1", 100..(100 + insert)));
    }

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.0).execute(None).unwrap();

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

/// Anomalous long-tail insert sizes (chimeras / mismapped pairs) should be
/// dropped from the histogram TSV — only entries within the
/// `median + deviations * MAD` window should survive.
#[test]
fn test_histogram_trims_long_tail() {
    let mut sam = SamBuilder::new();
    // 100 typical FR pairs spread over insert sizes 200..209.
    for i in 0..100usize {
        let insert = 200 + i % 10;
        sam.add(pair(sam.next_id()).fr("chr1", 100..(100 + insert)));
    }
    // 3 anomalous FR pairs with very large insert sizes that should not
    // appear in the trimmed output.
    for _ in 0..3 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..50_100));
    }

    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    make_cmd(bam.path(), &prefix, false, 0.0).execute(None).unwrap();

    let metrics: Vec<InsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    let fr = metrics.iter().find(|m| m.pair_orientation == "FR").expect("FR row");
    let trim_max = fr.median_insert_size + 10.0 * fr.median_absolute_deviation;
    // Sanity-check the test setup: tail values must actually be above the
    // computed trim threshold; otherwise the assertions below are vacuous.
    assert!(50_000.0 > trim_max, "tail (50000) should exceed trim_max ({trim_max})");

    let hist: Vec<InsertSizeHistogramEntry> =
        read_metrics_tsv(&dir.path().join(format!("out{HISTOGRAM_SUFFIX}"))).unwrap();
    assert!(!hist.is_empty(), "trimmed histogram should still contain the bulk distribution");
    // `trim_max` fits comfortably below `f64`'s 2^53 exact-integer range for
    // realistic insert-size data; the cast back to u64 just lets us compare
    // against the raw histogram keys.
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "trim_max is non-negative and within u64 range for any plausible insert size"
    )]
    let trim_cutoff = trim_max as u64;
    let beyond_trim: Vec<u64> =
        hist.iter().map(|e| e.insert_size).filter(|&s| s > trim_cutoff).collect();
    assert!(
        beyond_trim.is_empty(),
        "histogram contained entries past trim_max={trim_max}: {beyond_trim:?}"
    );
    // The bulk distribution must survive — guards against an over-aggressive
    // trim that drops everything.
    assert!(hist.iter().any(|e| e.insert_size == 200));
    assert!(hist.iter().any(|e| e.insert_size == 209));
    assert!(hist.iter().all(|e| e.insert_size != 50_000));
}

#[test]
fn test_plot_output_created() {
    let mut sam = SamBuilder::new();
    for _ in 0..5 {
        sam.add(pair(sam.next_id()).fr("chr1", 100..200));
    }
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");

    make_cmd(bam.path(), &prefix, false, 0.05).execute(None).unwrap();

    let plot_path = dir.path().join(format!("out{PLOT_SUFFIX}"));
    assert!(plot_path.exists(), "PDF plot file should be created");
    assert!(plot_path.metadata().unwrap().len() > 0, "PDF plot file should be non-empty");
}
