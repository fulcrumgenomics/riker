mod helpers;

use std::path::{Path, PathBuf};

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use riker_lib::commands::methyl::{
    BETA_SUFFIX, CONTIGS_SUFFIX, COVERAGE_SUFFIX, Chemistry, MBIAS_SUFFIX, METRICS_SUFFIX, Methyl,
    MethylBetaEntry, MethylContigMetric, MethylCoverageEntry, MethylMbiasEntry, MethylMetrics,
    MethylOptions, MethylReadConversionEntry, READ_CONVERSION_SUFFIX,
};
use riker_lib::test_support::{
    BedBuilder, FastaBuilder, Methylation, ReadBuilder, SamBuilder, coord_builder, pair, read,
};
use tempfile::TempDir;

// ─── Fixtures ────────────────────────────────────────────────────────────────

/// The 16 bp motif tiled to build `chr1`. Top-strand cytosines: a CpG at offset 2 and a CpA at
/// offset 7. Bottom-strand cytosines (reference Gs): the CpG's G at offset 3, a CpT at offset 12
/// (preceded by A) and a CpC at offset 13 (preceded by G). The motif starts with A and ends with
/// T, so tiling creates no extra CpGs.
const MOTIF: &[u8] = b"AACGTTACATTAGGTT";

/// Number of motif copies in `chr1`.
const MOTIF_COPIES: usize = 40;

/// Length of `chr1`.
const CHR1_LEN: usize = MOTIF.len() * MOTIF_COPIES;

/// Read length used throughout: four motif copies, so four CpGs per read.
const READ_LEN: usize = 64;

/// The tiled-motif reference: `chr1` (640 bp, 40 CpGs), plus an optional `lambda` control
/// contig of the same motif tiled 10 times.
fn reference(with_lambda: bool) -> FastaBuilder {
    let fasta = FastaBuilder::new().contig("chr1", MOTIF.repeat(MOTIF_COPIES));
    if with_lambda { fasta.contig("lambda", MOTIF.repeat(10)) } else { fasta }
}

/// A coordinate-sorted builder whose header matches [`reference`].
fn builder(with_lambda: bool) -> SamBuilder {
    if with_lambda {
        coord_builder(&[("chr1", CHR1_LEN), ("lambda", MOTIF.len() * 10)])
    } else {
        coord_builder(&[("chr1", CHR1_LEN)])
    }
}

/// A forward single-end read of [`READ_LEN`] bases at 1-based `start`, derived from the
/// reference.
fn fragment(name: &str, start: usize) -> ReadBuilder {
    read().name(name).at("chr1", start).len(READ_LEN).matching_ref(&reference(true))
}

/// Ten forward single-end reads tiling `chr1` end to end, each converted with `methylation`.
fn tiling_reads(methylation: &Methylation) -> Vec<ReadBuilder> {
    (0..10)
        .map(|i| fragment(&format!("r{i}"), 1 + i * READ_LEN).converted(methylation.clone()))
        .collect()
}

/// Write the reference and `reads` to temp files, run `riker methyl` with `options`, and return
/// the output directory (kept alive by the returned `TempDir`) and prefix.
fn run(with_lambda: bool, reads: Vec<ReadBuilder>, options: MethylOptions) -> (TempDir, PathBuf) {
    let fasta = reference(with_lambda).to_temp_fasta().unwrap();
    let mut sam = builder(with_lambda);
    for read in reads {
        sam.add(read);
    }
    let bam = sam.to_temp_bam().unwrap();
    run_bam(bam.path(), fasta.path(), options)
}

/// Run `riker methyl` on an existing BAM and reference.
fn run_bam(bam: &Path, fasta: &Path, options: MethylOptions) -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    Methyl {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        reference: ReferenceOptions { reference: fasta.to_path_buf() },
        options,
    }
    .execute(None)
    .unwrap();
    (dir, prefix)
}

/// The output file with `suffix` for `prefix`.
fn output(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.display()))
}

/// The single summary metrics row.
fn metrics(prefix: &Path) -> MethylMetrics {
    let mut rows: Vec<MethylMetrics> = read_metrics_tsv(&output(prefix, METRICS_SUFFIX)).unwrap();
    assert_eq!(rows.len(), 1);
    rows.remove(0)
}

// ─── Yield and methylation ───────────────────────────────────────────────────

#[test]
fn fully_cpg_methylated_reads_give_full_cpg_methylation_and_full_cph_conversion() {
    let (_dir, prefix) = run(false, tiling_reads(&Methylation::AllCpg), MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.cpg_territory, 40);
    assert_eq!(m.cpg_obs_potential, 40);
    assert_eq!(m.cpg_obs, 40);
    assert_close!(m.frac_excluded_total, 0.0);
    assert_close!(m.cpg_methylation, 1.0);
    assert_close!(m.cpg_conversion, 0.0);
    assert_close!(m.cpa_conversion, 1.0);
    assert_close!(m.cph_conversion, 1.0);
}

#[test]
fn unmethylated_reads_give_zero_cpg_methylation() {
    let (_dir, prefix) = run(false, tiling_reads(&Methylation::None), MethylOptions::default());
    let m = metrics(&prefix);
    assert_close!(m.cpg_methylation, 0.0);
    assert_close!(m.cph_conversion, 1.0);
}

#[test]
fn methylation_of_named_sites_is_reported_as_their_fraction() {
    // Sites 1-based at 3, 19, ... (offset 2 of each motif copy); methylate 10 of the 40.
    let methylated: Vec<usize> = (0..10).map(|i| 3 + i * MOTIF.len()).collect();
    let (_dir, prefix) =
        run(false, tiling_reads(&Methylation::CpgSites(methylated)), MethylOptions::default());
    assert_close!(metrics(&prefix).cpg_methylation, 0.25);
}

#[test]
fn taps_chemistry_reads_converted_cytosines_as_methylated() {
    let reads = (0..10)
        .map(|i| fragment(&format!("r{i}"), 1 + i * READ_LEN).taps_converted(Methylation::AllCpg))
        .collect();
    let options = MethylOptions { chemistry: Chemistry::Taps, ..MethylOptions::default() };
    let (_dir, prefix) = run(false, reads, options);
    let m = metrics(&prefix);
    assert_eq!(m.chemistry, "taps");
    assert_close!(m.cpg_methylation, 1.0);
    assert_close!(m.cph_conversion, 0.0);
}

#[test]
fn observations_are_normalised_per_read_and_per_mb() {
    let (_dir, prefix) = run(false, tiling_reads(&Methylation::AllCpg), MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.reads, 10);
    assert_eq!(m.bases_sequenced, 640);
    assert_close!(m.cpg_obs_per_read, 4.0);
    assert_close!(m.cpg_obs_per_mb_sequenced, 40.0 / 640.0 * 1_000_000.0, 0.01);
}

// ─── Strands and coverage ────────────────────────────────────────────────────

#[test]
fn top_and_bottom_strand_observations_merge_into_one_site() {
    // A forward and a reverse read over the same 64 bp: the forward read reports each CpG's C,
    // the reverse read its G, and both land on the same site.
    let reads = vec![
        fragment("fwd", 1).converted(Methylation::AllCpg),
        fragment("rev", 1).reverse().converted(Methylation::AllCpg),
    ];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let rows: Vec<MethylCoverageEntry> =
        read_metrics_tsv(&output(&prefix, COVERAGE_SUFFIX)).unwrap();
    assert_eq!(rows[2].cpgs, 4, "four sites at depth 2");
    assert_eq!(rows[0].cpgs, 36, "the other 36 sites at depth 0");
    assert_close!(metrics(&prefix).cpg_methylation, 1.0);
}

#[test]
fn bottom_strand_reads_measure_bottom_strand_cph_conversion() {
    let reads = (0..10)
        .map(|i| {
            fragment(&format!("r{i}"), 1 + i * READ_LEN).reverse().converted(Methylation::AllCpg)
        })
        .collect();
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_close!(m.cpt_conversion, 1.0);
    assert_close!(m.cpc_conversion, 1.0);
    assert_close!(m.cpg_methylation, 1.0);
    assert_eq!(m.cpg_obs, 40);
}

#[test]
fn coverage_statistics_cover_every_territory_cpg() {
    // Only the first half of chr1 is covered: 20 of 40 CpGs at depth 1.
    let reads = tiling_reads(&Methylation::AllCpg).into_iter().take(5).collect();
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_close!(m.frac_cpgs_at_1x, 0.5);
    assert_close!(m.mean_cpg_depth, 0.5);
}

// ─── Read-level exclusions ───────────────────────────────────────────────────

#[test]
fn read_level_exclusions_are_attributed_in_order() {
    // Five reads over the same four CpGs: one passing, one QC-fail and duplicate (counted as
    // QC-fail), one duplicate, one low-MAPQ, one QC-fail.
    let reads = vec![
        fragment("pass", 1).converted(Methylation::AllCpg),
        fragment("qc_dup", 1).converted(Methylation::AllCpg).qc_fail().duplicate(),
        fragment("dup", 1).converted(Methylation::AllCpg).duplicate(),
        fragment("mapq", 1).converted(Methylation::AllCpg).mapq(5),
        fragment("qc", 1).converted(Methylation::AllCpg).qc_fail(),
    ];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.cpg_obs_potential, 20);
    assert_eq!(m.cpg_obs, 4);
    assert_close!(m.frac_excluded_qc_fail, 0.4);
    assert_close!(m.frac_excluded_dupe, 0.2);
    assert_close!(m.frac_excluded_mapq, 0.2);
    assert_close!(m.frac_excluded_total, 0.8);
}

#[test]
fn include_duplicates_counts_duplicate_reads() {
    let reads = vec![
        fragment("pass", 1).converted(Methylation::AllCpg),
        fragment("dup", 1).converted(Methylation::AllCpg).duplicate(),
    ];
    let options = MethylOptions { include_duplicates: true, ..MethylOptions::default() };
    let (_dir, prefix) = run(false, reads, options);
    assert_eq!(metrics(&prefix).cpg_obs, 8);
}

#[test]
fn unmapped_reads_count_only_towards_reads_and_bases_sequenced() {
    let reads = vec![
        fragment("pass", 1).converted(Methylation::AllCpg),
        read().name("unmapped").len(READ_LEN).unmapped(),
    ];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.reads, 2);
    assert_eq!(m.bases_sequenced, 128);
    assert_eq!(m.cpg_obs_potential, 4);
}

#[test]
fn secondary_and_supplementary_alignments_are_ignored() {
    let reads = vec![
        fragment("pass", 1).converted(Methylation::AllCpg),
        fragment("secondary", 1).converted(Methylation::AllCpg).secondary(),
        fragment("supplementary", 65).converted(Methylation::AllCpg).supplementary(),
    ];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.reads, 1);
    assert_eq!(m.cpg_obs_potential, 4);
}

// ─── Base-level exclusions ───────────────────────────────────────────────────

#[test]
fn low_quality_bases_are_excluded() {
    let reads = vec![fragment("lowq", 1).converted(Methylation::AllCpg).qual(10)];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.cpg_obs, 0);
    assert_close!(m.frac_excluded_baseq, 1.0);
}

#[test]
fn a_cpg_base_that_is_neither_c_nor_t_is_non_informative() {
    // Offset 2 is the first CpG's C; read it as A.
    let reads = vec![fragment("snp", 1).converted(Methylation::AllCpg).sub(2, b'A')];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.cpg_obs, 3);
    assert_close!(m.frac_excluded_non_informative, 0.25);
}

#[test]
fn overlapping_mates_count_each_cpg_once() {
    // An FR pair whose mates cover the same 64 bp: both report the top strand, so read 2's
    // CpG calls duplicate read 1's.
    let reads_pair = pair("p1")
        .at("chr1", 1, 1)
        .len(READ_LEN)
        .matching_ref(&reference(false))
        .converted(Methylation::AllCpg);
    let fasta = reference(false).to_temp_fasta().unwrap();
    let mut sam = builder(false);
    sam.add(reads_pair);
    let bam = sam.to_temp_bam().unwrap();
    let (_dir, prefix) = run_bam(bam.path(), fasta.path(), MethylOptions::default());
    let m = metrics(&prefix);
    assert_eq!(m.cpg_obs_potential, 8);
    assert_eq!(m.cpg_obs, 4);
    assert_close!(m.frac_excluded_overlap, 0.5);
}

#[test]
fn end_masks_exclude_bases_by_cycle() {
    // Masking 3 cycles at the 5' end of read 1 covers offset 2, the first CpG's C.
    let reads = vec![fragment("r", 1).converted(Methylation::AllCpg)];
    let options = MethylOptions { mask_r1_5p: 3, ..MethylOptions::default() };
    let (_dir, prefix) = run(false, reads, options);
    let m = metrics(&prefix);
    assert_eq!(m.cpg_obs, 3);
    assert_close!(m.frac_excluded_end_mask, 0.25);
}

// ─── Territory, intervals and controls ───────────────────────────────────────

#[test]
fn intervals_restrict_the_cpg_territory() {
    // The first 64 bp hold four CpGs.
    let bed = BedBuilder::new().interval("chr1", 0, 64).to_temp_bed().unwrap();
    let options =
        MethylOptions { intervals: Some(bed.path().to_path_buf()), ..MethylOptions::default() };
    let (_dir, prefix) = run(false, tiling_reads(&Methylation::AllCpg), options);
    let m = metrics(&prefix);
    assert_eq!(m.cpg_territory, 4);
    assert_eq!(m.cpg_obs, 4);
}

#[test]
fn control_contigs_are_reported_separately() {
    let mut reads = tiling_reads(&Methylation::AllCpg);
    reads.push(
        read()
            .name("lambda")
            .at("lambda", 1)
            .len(READ_LEN)
            .matching_ref(&reference(true))
            .converted(Methylation::None),
    );
    let options = MethylOptions {
        control_contig: Some(vec!["lambda".to_string()]),
        ..MethylOptions::default()
    };
    let (_dir, prefix) = run(true, reads, options);
    let m = metrics(&prefix);
    assert_eq!(m.cpg_territory, 40);
    assert_close!(m.cpg_methylation, 1.0);

    let contigs: Vec<MethylContigMetric> =
        read_metrics_tsv(&output(&prefix, CONTIGS_SUFFIX)).unwrap();
    assert_eq!(contigs.len(), 2);
    let lambda = &contigs[1];
    assert_eq!(lambda.contig, "lambda");
    assert!(lambda.is_control);
    assert_eq!(lambda.cpg_sites, 10);
    assert_eq!(lambda.cpg_obs, 4);
    assert_close!(lambda.cpg_methylation, 0.0);
    assert_close!(lambda.cph_conversion, 1.0);
}

#[test]
fn an_unknown_control_contig_is_an_error() {
    let fasta = reference(false).to_temp_fasta().unwrap();
    let mut sam = builder(false);
    sam.add(fragment("r", 1).converted(Methylation::AllCpg));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let options = MethylOptions {
        control_contig: Some(vec!["phiX".to_string()]),
        ..MethylOptions::default()
    };
    let err = Methyl {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: dir.path().join("out") },
        reference: ReferenceOptions { reference: fasta.path().to_path_buf() },
        options,
    }
    .execute(None)
    .unwrap_err();
    assert!(err.to_string().contains("phiX"), "unexpected error: {err}");
}

// ─── M-bias, methylation levels and read conversion ──────────────────────────

#[test]
fn mbias_reports_calls_by_sequencing_cycle() {
    let reads =
        (0..3).map(|i| fragment(&format!("r{i}"), 1).converted(Methylation::AllCpg)).collect();
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let rows: Vec<MethylMbiasEntry> = read_metrics_tsv(&output(&prefix, MBIAS_SUFFIX)).unwrap();
    // The first CpG's C is at read offset 2: cycle 3.
    let first = rows.iter().find(|row| row.context == "CpG").unwrap();
    assert_eq!(first.read, "SE");
    assert_eq!(first.cycle, 3);
    assert_eq!(first.total, 3);
    assert_close!(first.frac_methylated, 1.0);
    // CpA cytosines sit at offset 7 of each motif: cycle 8 is the first CpH row.
    let first_cph = rows.iter().find(|row| row.context == "CpH").unwrap();
    assert_eq!(first_cph.cycle, 8);
    assert_close!(first_cph.frac_methylated, 0.0);
}

#[test]
fn mbias_counts_cycles_from_the_five_prime_end_of_reverse_reads() {
    // A reverse read's last stored base is its first cycle.
    let reads = vec![fragment("rev", 1).reverse().converted(Methylation::AllCpg)];
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let rows: Vec<MethylMbiasEntry> = read_metrics_tsv(&output(&prefix, MBIAS_SUFFIX)).unwrap();
    // The last CpG's G is at reference offset 51 of the read (motif 3, offset 3): stored
    // position 51 of 64, so cycle 64 - 51 = 13.
    let first = rows.iter().find(|row| row.context == "CpG").unwrap();
    assert_eq!(first.cycle, 13);
}

#[test]
fn beta_histogram_counts_sites_at_minimum_depth() {
    let reads =
        (0..10).map(|i| fragment(&format!("r{i}"), 1).converted(Methylation::AllCpg)).collect();
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let rows: Vec<MethylBetaEntry> = read_metrics_tsv(&output(&prefix, BETA_SUFFIX)).unwrap();
    assert_eq!(rows.len(), 101);
    assert_eq!(rows[100].cpgs, 4);
    let m = metrics(&prefix);
    assert_close!(m.mean_cpg_beta, 1.0);
    assert_close!(m.frac_cpgs_beta_high, 1.0);
}

#[test]
fn unconverted_reads_show_up_in_the_retained_cph_histogram() {
    // Each read spans four motifs: four CpA calls. One read is left unconverted.
    let mut reads = tiling_reads(&Methylation::AllCpg);
    reads[0] = fragment("r0", 1);
    let (_dir, prefix) = run(false, reads, MethylOptions::default());
    let rows: Vec<MethylReadConversionEntry> =
        read_metrics_tsv(&output(&prefix, READ_CONVERSION_SUFFIX)).unwrap();
    assert_eq!(rows[0].reads, 9);
    assert_eq!(rows[4].reads, 1);
}

#[test]
fn directional_reads_show_no_off_strand_conversion() {
    let (_dir, prefix) = run(false, tiling_reads(&Methylation::AllCpg), MethylOptions::default());
    assert_close!(metrics(&prefix).off_strand_cph_conversion, 0.0);
}

#[test]
fn reads_converted_on_the_other_strand_show_off_strand_conversion() {
    // A forward read whose bottom-strand CpH Gs (offsets 12 and 13 of each motif) read as A, as
    // a read from a non-directional library would.
    let mut builder = fragment("r", 1).converted(Methylation::AllCpg);
    for copy in 0..4 {
        builder = builder.sub(copy * MOTIF.len() + 12, b'A').sub(copy * MOTIF.len() + 13, b'A');
    }
    let (_dir, prefix) = run(false, vec![builder], MethylOptions::default());
    assert_close!(metrics(&prefix).off_strand_cph_conversion, 1.0);
}

#[test]
fn unsorted_input_is_rejected() {
    let fasta = reference(false).to_temp_fasta().unwrap();
    let mut sam = SamBuilder::with_contigs(&[("chr1".to_string(), CHR1_LEN)]);
    sam.add(fragment("r", 1).converted(Methylation::AllCpg));
    let bam = sam.to_temp_bam().unwrap();
    let dir = TempDir::new().unwrap();
    let err = Methyl {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: dir.path().join("out") },
        reference: ReferenceOptions { reference: fasta.path().to_path_buf() },
        options: MethylOptions::default(),
    }
    .execute(None)
    .unwrap_err();
    assert!(err.to_string().contains("coordinate-sorted"), "unexpected error: {err}");
}
