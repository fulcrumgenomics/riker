mod helpers;

use std::path::{Path, PathBuf};

use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions};
use riker_lib::commands::error::{
    Error, ErrorOptions, INDEL_SUFFIX, IndelMetric, MISMATCH_SUFFIX, MismatchMetric,
    OVERLAP_SUFFIX, OverlappingMismatchMetric,
};
use riker_lib::test_support::{BedBuilder, FastaBuilder, VcfBuilder, coord_builder, pair, read};

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// A 1000 bp reference of repeating `ACGT`. Reads derive their bases from it via
/// `read().matching_ref(&fasta)`, so a read matches the reference by default and a
/// `.sub(offset, base)` plants a single mismatch; its temp file is the `--reference`.
fn reference() -> FastaBuilder {
    FastaBuilder::new().contig("chr1", (0..1000).map(|i| b"ACGT"[i % 4]).collect::<Vec<u8>>())
}

/// Default error options: MAPQ ≥ 20, every base passes (`min_bq` 0), no masking.
/// The `reference` field is a placeholder overwritten by [`make_cmd`].
fn opts() -> ErrorOptions {
    ErrorOptions {
        reference: PathBuf::new(),
        vcf: None,
        intervals: None,
        min_mapq: 20,
        min_bq: 0,
        include_duplicates: false,
        max_isize: 1000,
        picard_compat: false,
        stratify_by: Vec::new(),
    }
}

/// Build an `error` command from a BAM, reference, output prefix, and options.
fn make_cmd(bam: &Path, ref_file: &Path, output: &Path, options: ErrorOptions) -> Error {
    Error {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: output.to_path_buf() },
        options: ErrorOptions { reference: ref_file.to_path_buf(), ..options },
    }
}

/// Run `error` with the given options and return the output prefix. The temp dir
/// is leaked so the output files persist for assertions.
fn run_opts(bam: &Path, ref_file: &Path, options: ErrorOptions) -> PathBuf {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    make_cmd(bam, ref_file, &prefix, options).execute(None).expect("error command should succeed");
    std::mem::forget(dir);
    prefix
}

/// Run `error` with default options and the given stratifiers.
fn run(bam: &Path, ref_file: &Path, stratify_by: Vec<String>) -> PathBuf {
    run_opts(bam, ref_file, ErrorOptions { stratify_by, ..opts() })
}

/// Read the mismatch metric rows from an output prefix.
fn mismatch_rows(prefix: &Path) -> Vec<MismatchMetric> {
    read_metrics_tsv(&prefix.with_file_name(format!("out{MISMATCH_SUFFIX}"))).unwrap()
}

/// Read the indel metric rows from an output prefix.
fn indel_rows(prefix: &Path) -> Vec<IndelMetric> {
    read_metrics_tsv(&prefix.with_file_name(format!("out{INDEL_SUFFIX}"))).unwrap()
}

/// Read the overlapping-mismatch metric rows from an output prefix.
fn overlap_rows(prefix: &Path) -> Vec<OverlappingMismatchMetric> {
    read_metrics_tsv(&prefix.with_file_name(format!("out{OVERLAP_SUFFIX}"))).unwrap()
}

/// The `all` stratifier row from a set of mismatch rows.
fn all_mismatch(rows: &[MismatchMetric]) -> &MismatchMetric {
    rows.iter().find(|r| r.stratifier == "all").unwrap()
}

/// The `all` stratifier row from a set of indel rows.
fn all_indel(rows: &[IndelMetric]) -> &IndelMetric {
    rows.iter().find(|r| r.stratifier == "all").unwrap()
}

// ─── Core correctness tests ─────────────────────────────────────────────────

#[test]
fn test_no_errors() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("read1").at("chr1", 100, 120).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    // Should have one row for "all" stratifier with zero errors
    assert!(!mm.is_empty());
    let all_row = all_mismatch(&mm);
    assert!(all_row.total_bases > 0);
    assert_eq!(all_row.error_bases, 0);
    assert_close!(all_row.frac_error, 0.0, 1e-9);
}

#[test]
fn test_simple_mismatches() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A read matching the reference except for 2 planted mismatches at offsets 0 and 5.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G').sub(5, b'T'));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 2);
    assert_close!(all_row.frac_error, 0.2, 1e-6);
}

#[test]
fn test_stratify_by_strand() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    // Forward read with 1 mismatch
    sam.add(read().name("read_fwd").at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G'));
    // Reverse read with 0 mismatches
    sam.add(read().name("read_rev").at("chr1", 201).len(10).matching_ref(&fasta).reverse());
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["strand".to_string()]);
    let mm = mismatch_rows(&prefix);

    let fwd = mm.iter().find(|r| r.stratifier == "strand" && r.covariate == "+").unwrap();
    let rev = mm.iter().find(|r| r.stratifier == "strand" && r.covariate == "-").unwrap();

    assert_eq!(fwd.total_bases, 10);
    assert_eq!(fwd.error_bases, 1);
    assert_eq!(rev.total_bases, 10);
    assert_eq!(rev.error_bases, 0);
}

#[test]
fn test_insertion_detection() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read with 10M3I (10 aligned + 3 inserted bases)
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("10M3I").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    assert_eq!(all_row.total_bases, 10); // only aligned bases count for total
    assert_eq!(all_row.num_insertions, 1);
    assert_eq!(all_row.num_inserted_bases, 3);
    assert_eq!(all_row.num_deletions, 0);
}

#[test]
fn test_deletion_detection() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read with 5M2D5M (5 match, 2 deleted, 5 match)
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("5M2D5M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    assert_eq!(all_row.total_bases, 10); // 5 + 5 aligned bases
    assert_eq!(all_row.num_deletions, 1);
    assert_eq!(all_row.num_deleted_bases, 2);
    assert_eq!(all_row.num_insertions, 0);
}

#[test]
fn test_min_mapq_filter() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read with MAPQ 10 (below default min of 20)
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).mapq(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    // All bases should be filtered out — empty or zero
    if let Some(all_row) = mm.iter().find(|r| r.stratifier == "all") {
        assert_eq!(all_row.total_bases, 0);
    }
    // Or the file may have no rows at all, which is also valid
}

#[test]
fn test_duplicate_exclusion() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    // Duplicate read (matches reference)
    sam.add(read().name("read_dup").at("chr1", 101).len(10).matching_ref(&fasta).duplicate());
    // Non-duplicate read with 1 mismatch
    sam.add(read().name("read_ok").at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G'));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    // Only the non-duplicate read should be counted
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 1);
}

#[test]
fn test_all_group_always_present() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    // User specifies only "bq" but "all" should still be present
    let prefix = run(bam.path(), ref_file.path(), vec!["bq".to_string()]);
    let mm = mismatch_rows(&prefix);

    // Should have both "all" and "bq" rows
    assert!(mm.iter().any(|r| r.stratifier == "all"));
    assert!(mm.iter().any(|r| r.stratifier == "bq"));
}

#[test]
fn test_composite_stratifier() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["strand,mapq".to_string()]);
    let mm = mismatch_rows(&prefix);

    // Should have rows for the composite "strand,mapq" group
    let composite_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "strand,mapq").collect();
    assert!(!composite_rows.is_empty());
    // The covariate should be comma-separated, e.g. "+,60"
    assert!(composite_rows.iter().any(|r| r.covariate.contains(',')));
}

#[test]
fn test_overlapping_reads_mismatching_ref_and_mate() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Two overlapping reads at positions 101-110 and 106-115.
    // Read 1 matches the reference; read 2 mismatches at 0-based ref 107 (read2 offset 2).
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("pair1").at("chr1", 101, 106).len(10).matching_ref(&fasta).r2(|r| r.sub(2, b'A')));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    // We should have some overlap data
    if let Some(all_row) = ov.iter().find(|r| r.stratifier == "all") {
        // There should be overlapping bases examined
        assert!(all_row.overlapping_read_bases > 0, "Expected overlapping bases");
    }
}

#[test]
fn test_three_output_files_always_created() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);

    // All three output files should exist
    assert!(prefix.with_file_name("out.error-mismatch.txt").exists());
    assert!(prefix.with_file_name("out.error-overlap.txt").exists());
    assert!(prefix.with_file_name("out.error-indel.txt").exists());
}

#[test]
fn test_min_bq_filter() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A read where bases 0..5 have quality 10 and bases 5..10 have quality 30.
    let mut quals = vec![10u8; 10];
    for q in quals.iter_mut().skip(5) {
        *q = 30;
    }

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta).quals(quals));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run_opts(bam.path(), ref_file.path(), ErrorOptions { min_bq: 20, ..opts() });
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    // Only the 5 high-quality bases should be counted
    assert_eq!(all_row.total_bases, 5);
}

#[test]
fn test_n_bases_excluded() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A 10bp read matching the reference but with N bases at offsets 4 and 9.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta).sub(4, b'N').sub(9, b'N'));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    // 2 N bases should be excluded, leaving 8 valid bases
    assert_eq!(all_row.total_bases, 8);
    assert_eq!(all_row.error_bases, 0);
}

#[test]
fn test_soft_clip_excluded() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // 10bp read with CIGAR 3S7M: 3 soft-clipped + 7 matched, aligned at 101.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("3S7M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    assert_eq!(all_row.total_bases, 7);
}

#[test]
fn test_secondary_supplementary_excluded() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    // Primary read with 1 mismatch
    sam.add(read().name("read1").at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G'));
    // Secondary read with 1 mismatch
    sam.add(
        read().name("read2").at("chr1", 101).len(10).matching_ref(&fasta).sub(1, b'T').secondary(),
    );
    // Supplementary read with 1 mismatch
    sam.add(
        read()
            .name("read3")
            .at("chr1", 101)
            .len(10)
            .matching_ref(&fasta)
            .sub(2, b'T')
            .supplementary(),
    );
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    // Only the primary read should be counted
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 1);
}

#[test]
fn test_mixed_insertion_and_deletion() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // CIGAR: 5M2I5M3D5M — 15 aligned bases, one 2bp insertion, one 3bp deletion.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("5M2I5M3D5M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    assert_eq!(all_row.total_bases, 15);
    assert_eq!(all_row.num_insertions, 1);
    assert_eq!(all_row.num_inserted_bases, 2);
    assert_eq!(all_row.num_deletions, 1);
    assert_eq!(all_row.num_deleted_bases, 3);
}

#[test]
fn test_q_score_calculation() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // 100 reads of 10bp each = 1000 total bases, with exactly 1 mismatch (in the
    // first read). Stack all reads at the same position to stay within the reference.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    for i in 0..100 {
        let mut r = read().at("chr1", 101).len(10).matching_ref(&fasta);
        if i == 0 {
            r = r.sub(0, b'G');
        }
        sam.add(r);
    }
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    assert_eq!(all_row.total_bases, 1000);
    assert_eq!(all_row.error_bases, 1);
    // q_score = -10 * log10(1/1000) = 30.0
    assert_close!(all_row.q_score, 30.0, 0.01);
}

// ─── Stratifier correctness tests ───────────────────────────────────────────

#[test]
fn test_stratify_by_cycle() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    // Forward read with a mismatch at read offset 2 (cycle 3 for forward reads)
    sam.add(read().name("read_fwd").at("chr1", 101).len(10).matching_ref(&fasta).sub(2, b'T'));
    // Reverse-strand read with a mismatch at read offset 7. For reverse reads the
    // cycle counts down from read_len, so offset 7 -> cycle = 10 - 7 = 3.
    sam.add(
        read().name("read_rev").at("chr1", 201).len(10).matching_ref(&fasta).sub(7, b'A').reverse(),
    );
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["cycle".to_string()]);
    let mm = mismatch_rows(&prefix);

    let cycle3 = mm.iter().find(|r| r.stratifier == "cycle" && r.covariate == "3").unwrap();
    // Both reads contribute a mismatch at cycle 3
    assert_eq!(cycle3.error_bases, 2);
}

#[test]
fn test_stratify_by_read_num() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read1 (forward) with 1 mismatch, Read2 (reverse) with 0 mismatches.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("pair1").at("chr1", 101, 120).len(10).matching_ref(&fasta).r1(|r| r.sub(0, b'G')));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["read_num".to_string()]);
    let mm = mismatch_rows(&prefix);

    let r1 = mm.iter().find(|r| r.stratifier == "read_num" && r.covariate == "R1").unwrap();
    let r2 = mm.iter().find(|r| r.stratifier == "read_num" && r.covariate == "R2").unwrap();
    assert_eq!(r1.error_bases, 1);
    assert_eq!(r2.error_bases, 0);
}

#[test]
fn test_stratify_by_ref_base() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Forward read aligned at 101 (1-based). Reference at 0-based 100..110 is
    // A C G T A C G T A C; a mismatch at offset 0 is at ref base 'A'.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G'));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["ref_base".to_string()]);
    let mm = mismatch_rows(&prefix);

    // The mismatch at ref_base 'A' should show up
    let ref_a = mm.iter().find(|r| r.stratifier == "ref_base" && r.covariate == "A").unwrap();
    assert!(ref_a.error_bases >= 1);
    // Other ref bases should have 0 errors (we only introduced one mismatch at an 'A')
    let ref_c = mm.iter().find(|r| r.stratifier == "ref_base" && r.covariate == "C").unwrap();
    assert_eq!(ref_c.error_bases, 0);
}

#[test]
fn test_stratify_by_hp_len() {
    // Reference with a homopolymer run: ACGT pattern, then AAAA at 0-based 100-103.
    let mut custom_seq: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();
    custom_seq.extend_from_slice(b"AAAA");
    custom_seq.extend((104..1000).map(|i| b"ACGT"[i % 4]));
    let fasta = FastaBuilder::new().contig("chr1", custom_seq);
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A forward read at 1-based 98 spanning the homopolymer region, matching the ref.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 98).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["hp_len".to_string()]);
    let mm = mismatch_rows(&prefix);

    // The hp_len stratifier should produce multiple covariate values including "0"
    // and some non-zero values for bases after the AAAA homopolymer.
    let hp_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "hp_len").collect();
    assert!(!hp_rows.is_empty());
    // Should have a row with hp_len=0 (for bases not following a homopolymer)
    assert!(hp_rows.iter().any(|r| r.covariate == "0"));
    // Should have some non-zero hp_len values due to the AAAA run
    assert!(hp_rows.iter().any(|r| {
        let val: i64 = r.covariate.parse().unwrap_or(0);
        val > 0
    }));
}

#[test]
fn test_stratify_by_indel_len() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read with 5M3I5M: 5 match + 3 inserted + 5 match.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("5M3I5M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["indel_len".to_string()]);
    let indels = indel_rows(&prefix);

    // Should have an indel_len=3 row for the 3bp insertion
    let indel_3 = indels.iter().find(|r| r.stratifier == "indel_len" && r.covariate == "3");
    assert!(indel_3.is_some(), "Expected indel_len=3 row for the 3bp insertion");
    let indel_3 = indel_3.unwrap();
    assert_eq!(indel_3.num_insertions, 1);

    // Mismatch metrics should have rows with indel_len=0 for aligned bases
    let mm = mismatch_rows(&prefix);
    let mm_0 = mm.iter().find(|r| r.stratifier == "indel_len" && r.covariate == "0");
    assert!(mm_0.is_some(), "Expected indel_len=0 row for aligned bases in mismatch metrics");
    let mm_0 = mm_0.unwrap();
    assert_eq!(mm_0.total_bases, 10); // 5 + 5 aligned bases
}

// ─── Overlap detection tests ────────────────────────────────────────────────

#[test]
fn test_overlapping_reads_all_agree() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Two overlapping reads (101-110 and 106-115) that both match the reference.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("pair1").at("chr1", 101, 106).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0, "Expected overlapping bases");
    assert_eq!(all_row.bases_mismatching_ref_and_mate, 0);
    assert_eq!(all_row.bases_matching_mate_but_not_ref, 0);
    assert_eq!(all_row.bases_in_three_way_disagreement, 0);
}

#[test]
fn test_overlapping_reads_matching_mate_but_not_ref() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Both reads carry the SAME mismatch at 0-based ref 107: read1 offset 7, read2
    // offset 2. They disagree with the reference but agree with each other.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(
        pair("pair1")
            .at("chr1", 101, 106)
            .len(10)
            .matching_ref(&fasta)
            .r1(|r| r.sub(7, b'A'))
            .r2(|r| r.sub(2, b'A')),
    );
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0);
    assert!(
        all_row.bases_matching_mate_but_not_ref > 0,
        "Expected bases_matching_mate_but_not_ref > 0, got {}",
        all_row.bases_matching_mate_but_not_ref
    );
}

#[test]
fn test_overlapping_reads_three_way() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Three-way disagreement at 0-based ref 107: ref='T', read1='A', read2='C'.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(
        pair("pair1")
            .at("chr1", 101, 106)
            .len(10)
            .matching_ref(&fasta)
            .r1(|r| r.sub(7, b'A'))
            .r2(|r| r.sub(2, b'C')),
    );
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0);
    assert!(
        all_row.bases_in_three_way_disagreement > 0,
        "Expected bases_in_three_way_disagreement > 0, got {}",
        all_row.bases_in_three_way_disagreement
    );
}

#[test]
fn test_overlapping_reads_mismatching_ref_and_mate_exact() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read 1 matches the reference; read 2 mismatches at 0-based ref 107 (offset 2).
    // In the overlap region read1 agrees with ref there, so this is exactly one
    // base_mismatching_ref_and_mate.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("pair1").at("chr1", 101, 106).len(10).matching_ref(&fasta).r2(|r| r.sub(2, b'A')));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0, "Expected overlapping bases");
    // Exactly 1 base where read2 disagrees with ref but read1 (mate) agrees
    assert_eq!(all_row.bases_mismatching_ref_and_mate, 1);
    assert_eq!(all_row.bases_matching_mate_but_not_ref, 0);
    assert_eq!(all_row.bases_in_three_way_disagreement, 0);
}

// ─── Edge case tests ────────────────────────────────────────────────────────

#[test]
fn test_multi_contig() {
    let acgt = |n: usize| (0..n).map(|i| b"ACGT"[i % 4]).collect::<Vec<u8>>();
    let fasta = FastaBuilder::new().contig("chr1", acgt(500)).contig("chr2", acgt(500));
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 500), ("chr2", 500)]);
    sam.add(read().name("read_chr1").at("chr1", 101).len(10).matching_ref(&fasta));
    sam.add(read().name("read_chr2").at("chr2", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    // Both reads should be counted: 10 + 10 = 20 total bases
    assert_eq!(all_row.total_bases, 20);
}

#[test]
fn test_deletion_at_read_start() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // CIGAR: 2D10M — deletion before any aligned base, so no anchor.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("2D10M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    assert_eq!(all_row.total_bases, 10);
    // No anchor before the deletion, so it should be silently skipped
    assert_eq!(all_row.num_deletions, 0);
}

#[test]
fn test_insertion_at_read_start() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // CIGAR: 3I10M — insertion before any aligned base, so no anchor.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("3I10M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    // No anchor before the insertion, so it should be skipped
    assert_eq!(all_row.num_insertions, 0);
}

#[test]
fn test_stratifier_parse_error() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let options = ErrorOptions { stratify_by: vec!["invalid_name".to_string()], ..opts() };
    let cmd = make_cmd(bam.path(), ref_file.path(), &prefix, options);

    let result = cmd.execute(None);
    assert!(result.is_err(), "Expected error for invalid stratifier name");
}

// ─── Additional stratifier tests ────────────────────────────────────────────

#[test]
fn test_gc_stratification() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read at 1-based 101 matching ref[100..110] = A C G T A C G T A C: 5 GC bases.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["gc".to_string()]);
    let mm = mismatch_rows(&prefix);

    // Sequence ACGTACGTAC has 5 GC bases out of 10 -> GC% = (5*100 + 5)/10 = 50.
    let gc_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "gc").collect();
    assert!(!gc_rows.is_empty());
    // All 10 bases should be in the gc=50 bucket
    let gc50 = gc_rows.iter().find(|r| r.covariate == "50");
    assert!(
        gc50.is_some(),
        "Expected gc=50 row, found: {:?}",
        gc_rows.iter().map(|r| &r.covariate).collect::<Vec<_>>()
    );
    assert_eq!(gc50.unwrap().total_bases, 10);
}

#[test]
fn test_pre_dinuc_stratification() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Forward read at 1-based 101 matching ref[100..110] = A C G T A C G T A C.
    // pre_dinuc = previous read base + current ref base; at offset 1 that is "AC".
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["pre_dinuc".to_string()]);
    let mm = mismatch_rows(&prefix);

    let dinuc_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "pre_dinuc").collect();
    assert!(!dinuc_rows.is_empty());
    // "AC" should appear (prev read base 'A' + current ref base 'C')
    assert!(
        dinuc_rows.iter().any(|r| r.covariate == "AC"),
        "Expected 'AC' dinuc, found: {:?}",
        dinuc_rows.iter().map(|r| &r.covariate).collect::<Vec<_>>()
    );
}

#[test]
fn test_context_3bp_stratification() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Forward read at 1-based 101 matching ref[100..110] = A C G T A C G T A C.
    // context_3bp at offset 1 = prev 'A' + ref 'C' + next 'G' -> "ACG".
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec!["context_3bp".to_string()]);
    let mm = mismatch_rows(&prefix);

    let ctx_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "context_3bp").collect();
    assert!(!ctx_rows.is_empty());
    // "ACG" should appear
    assert!(
        ctx_rows.iter().any(|r| r.covariate == "ACG"),
        "Expected 'ACG' context, found: {:?}",
        ctx_rows.iter().map(|r| &r.covariate).collect::<Vec<_>>()
    );
}

/// A read that looks like it should overlap its mate (paired, same contig, small insert) but
/// whose mate is not present in the BAM should still be processed for mismatch/indel errors.
#[test]
fn test_orphaned_buffered_read_still_counted() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A single read of a pair (mate at 106, small insert) whose mate is absent from
    // the BAM. classify_overlap buffers it, and it must still be counted. It carries
    // a mismatch at 0-based ref 105 (offset 5).
    let (orphan, _mate) = pair("orphan")
        .at("chr1", 101, 106)
        .len(10)
        .matching_ref(&fasta)
        .r1(|r| r.sub(5, b'G'))
        .build();
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(orphan);
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let mm = mismatch_rows(&prefix);

    let all_row = all_mismatch(&mm);
    assert_eq!(all_row.total_bases, 10, "Orphaned read's 10 bases should be counted");
    assert_eq!(all_row.error_bases, 1, "Orphaned read's mismatch should be counted");
}

// ─── Additional tests ────────────────────────────────────────────────────────

/// An insertion where the first inserted base has BQ below min_bq should cause the entire
/// insertion to be skipped, while aligned bases are still counted.
#[test]
fn test_insertion_low_bq_first_base_excluded() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // CIGAR: 30M2I44M = 76 read bases, 74 aligned. All BQ=30 except the first
    // inserted base (read offset 30) at BQ=5.
    let mut quals = vec![30u8; 76];
    quals[30] = 5;

    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("30M2I44M").matching_ref(&fasta).quals(quals));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run_opts(bam.path(), ref_file.path(), ErrorOptions { min_bq: 20, ..opts() });
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    // The insertion should be entirely skipped because the first inserted base has BQ < min_bq
    assert_eq!(
        all_row.num_insertions, 0,
        "Insertion should be skipped when first base BQ < min_bq"
    );
    assert_eq!(all_row.num_inserted_bases, 0);
    // The 74 aligned bases (all BQ=30) should still be counted.
    assert_eq!(all_row.total_bases, 74);
}

/// When both reads in an overlapping pair have a mismatch at the same position (both differ from
/// reference but agree with each other), `bases_matching_mate_but_not_ref` should be > 0 and
/// `overlapping_read_bases` should count both reads' contributions to the overlap.
#[test]
fn test_overlap_double_counted() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Both reads mismatch at 0-based ref 106: read1 offset 6, read2 offset 1.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(
        pair("pair1")
            .at("chr1", 101, 106)
            .len(10)
            .matching_ref(&fasta)
            .r1(|r| r.sub(6, b'T'))
            .r2(|r| r.sub(1, b'T')),
    );
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    // The overlap region is 106-110 (5 positions); both reads are counted, so
    // overlapping_read_bases should be 2 * 5 = 10.
    assert_eq!(
        all_row.overlapping_read_bases, 10,
        "Both reads should contribute to overlapping_read_bases"
    );
    // Both reads disagree with reference at the same position but agree with each other
    assert!(
        all_row.bases_matching_mate_but_not_ref > 0,
        "Expected bases_matching_mate_but_not_ref > 0, got {}",
        all_row.bases_matching_mate_but_not_ref
    );
}

/// Reads with insert size exceeding max_isize should be excluded from the isize stratifier
/// but still counted in the "all" group.
#[test]
fn test_max_isize_exclusion() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    let mut sam = coord_builder(&[("chr1", 1000)]);
    // Pair 1: small insert size (tlen=15, overlapping).
    sam.add(pair("small_pair").at("chr1", 101, 106).len(10).matching_ref(&fasta));
    // Pair 2: large insert size (tlen=500).
    sam.add(pair("large_pair").at("chr1", 101, 591).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let options =
        ErrorOptions { max_isize: 100, stratify_by: vec!["all".into(), "isize".into()], ..opts() };
    let prefix = run_opts(bam.path(), ref_file.path(), options);
    let mm = mismatch_rows(&prefix);

    // The "all" row should count bases from BOTH pairs (4 reads x 10 bases = 40)
    let all_row = all_mismatch(&mm);
    assert_eq!(all_row.total_bases, 40, "All 4 reads should be counted in the 'all' group");

    // The isize stratifier should have a row for tlen=15 (small pair)
    let isize_15 = mm.iter().find(|r| r.stratifier == "isize" && r.covariate == "15");
    assert!(isize_15.is_some(), "Expected isize=15 row for the small pair");

    // The isize stratifier should NOT have a row for tlen=500 (exceeds max_isize=100)
    let isize_500 = mm.iter().find(|r| r.stratifier == "isize" && r.covariate == "500");
    assert!(isize_500.is_none(), "tlen=500 should be excluded when max_isize=100");
}

/// A pair where read1 and read2 are well-separated with no overlap should produce zero
/// overlapping bases in the overlap metrics.
#[test]
fn test_non_overlapping_pair_no_overlap_metrics() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read1 at 101-110, read2 at 200-209. No overlap.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(pair("pair1").at("chr1", 101, 200).len(10).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let ov = overlap_rows(&prefix);

    // If there is an "all" overlap row, it should have zero overlapping bases
    if let Some(all_row) = ov.iter().find(|r| r.stratifier == "all") {
        assert_eq!(
            all_row.overlapping_read_bases, 0,
            "Non-overlapping pair should have 0 overlapping read bases"
        );
    }
    // Otherwise, no rows at all is also acceptable
}

/// An insertion at the very start of a read (CIGAR: 3I50M) has no preceding anchor base,
/// so the insertion should not be counted, but the aligned bases should still be counted.
#[test]
fn test_insertion_at_read_start_no_anchor() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // CIGAR: 3I50M — insertion before any aligned base, so no anchor.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("3I50M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    let prefix = run(bam.path(), ref_file.path(), vec![]);
    let indels = indel_rows(&prefix);

    let all_row = all_indel(&indels);
    // No anchor before the insertion, so it should be skipped
    assert_eq!(
        all_row.num_insertions, 0,
        "Insertion at read start should have no anchor and be skipped"
    );
    // The 50 aligned bases should still be counted
    assert_eq!(all_row.total_bases, 50);

    // Verify mismatch metrics also count the aligned bases
    let mm = mismatch_rows(&prefix);
    let mm_all = all_mismatch(&mm);
    assert_eq!(mm_all.total_bases, 50, "Aligned M bases should still be counted for mismatches");
}

/// Regression: an orphan buffered mid-read on contig A must not be processed against
/// contig B's region when the next interval crosses contigs.
///
/// Setup: chr1 is all A, chr2 is all G. A coordinate-sorted pair lives on chr1, with
/// read1 at pos 1 (1-based) and read2 at pos 30, each 50M. The intervals file covers
/// `chr1:0-15` and all of `chr2`. Read1 is returned by the chr1 query and buffered
/// (mate's expected start at pos 30 falls inside read1's [1, 50] span). Read2 lives
/// outside the chr1 interval, so it never appears in the query and the buffered
/// mate stays orphaned across the contig boundary.
///
/// Under the bug, the end-of-chr2 `flush_behind(chr2, 200)` drains read1 (its
/// `mate_ref_id=chr1 < chr2`) and runs it through the chr2 region. All 50 A-bases
/// then get compared against chr2's G reference, inflating `error_bases` to 50.
///
/// Correctly handled, the orphan is flushed against chr1's *last* region before
/// swapping to chr2; the 15 positions inside `[0, 15)` match A-vs-A and the remaining
/// 35 positions are skipped by `RegionContext::contains`, so `error_bases` is 0.
#[test]
fn test_cross_contig_orphan_not_processed_against_wrong_contig() {
    let fasta = FastaBuilder::new().contig("chr1", vec![b'A'; 200]).contig("chr2", vec![b'G'; 200]);
    let ref_file = fasta.to_temp_fasta().unwrap();

    // read1 at 1-based pos 1, mate at pos 30 — mate_pos falls inside read1's [1, 50]
    // span so the buffer keeps it pending.
    let mut sam = coord_builder(&[("chr1", 200), ("chr2", 200)]);
    sam.add(pair("read").at("chr1", 1, 30).len(50).matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    // BED: chr1 partial interval ending *before* read1's mate_pos (29); chr2 full.
    let dir = tempfile::tempdir().unwrap();
    let bed =
        BedBuilder::new().interval("chr1", 0, 15).interval("chr2", 0, 200).to_temp_bed().unwrap();

    let prefix = dir.path().join("out");
    let options = ErrorOptions { intervals: Some(bed.path().to_path_buf()), ..opts() };
    make_cmd(bam.path(), ref_file.path(), &prefix, options)
        .execute(None)
        .expect("error command should succeed");

    let mm = mismatch_rows(&prefix);
    let all_row = all_mismatch(&mm);
    assert_eq!(
        all_row.error_bases, 0,
        "orphan on chr1 must not be compared against chr2's reference bases"
    );
}

// ─── Known-sites (VCF) masking tests ────────────────────────────────────────

#[test]
fn test_vcf_masks_known_variant_mismatch() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // A read matching the reference except a mismatch at offset 0 (1-based pos 101).
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).len(10).matching_ref(&fasta).sub(0, b'G'));
    let bam = sam.to_temp_indexed_bam().unwrap();

    // Without a VCF, the mismatch counts as an error over all 10 bases.
    let plain = run(bam.path(), ref_file.path(), vec![]);
    let plain_rows = mismatch_rows(&plain);
    let all = all_mismatch(&plain_rows);
    assert_eq!(all.total_bases, 10);
    assert_eq!(all.error_bases, 1);

    // A known variant at chr1:101 masks that position: the mismatch is excluded and
    // the masked base drops out of the denominator too.
    let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
    vcf.add("chr1", 101, "A").alt("G");
    let vcf_file = vcf.to_temp_vcf().unwrap();

    let options = ErrorOptions { vcf: Some(vcf_file.path().to_path_buf()), ..opts() };
    let masked = run_opts(bam.path(), ref_file.path(), options);
    let masked_rows = mismatch_rows(&masked);
    let all = all_mismatch(&masked_rows);
    assert_eq!(all.total_bases, 9, "the masked position drops from the denominator");
    assert_eq!(all.error_bases, 0, "the mismatch at the known-variant site is not counted");
}

#[test]
fn test_vcf_masks_deletion_in_variant_span() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read 5M2D5M at 101 deletes 1-based reference positions 106 and 107.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("5M2D5M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    // Without a VCF the deletion is reported.
    let plain = run(bam.path(), ref_file.path(), vec![]);
    let plain_rows = indel_rows(&plain);
    let all = all_indel(&plain_rows);
    assert_eq!(all.num_deletions, 1);
    assert_eq!(all.num_deleted_bases, 2);

    // A known 2 bp deletion (REF "CG" at 106) masks the whole deleted span, so the
    // observed deletion is excluded.
    let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
    vcf.add("chr1", 106, "CG").alt("C");
    let vcf_file = vcf.to_temp_vcf().unwrap();

    let options = ErrorOptions { vcf: Some(vcf_file.path().to_path_buf()), ..opts() };
    let masked = run_opts(bam.path(), ref_file.path(), options);
    let masked_rows = indel_rows(&masked);
    let all = all_indel(&masked_rows);
    assert_eq!(all.num_deletions, 0, "a deletion fully within a known-variant span is excluded");
}

#[test]
fn test_vcf_masks_insertion_at_variant_anchor() {
    let fasta = reference();
    let ref_file = fasta.to_temp_fasta().unwrap();

    // Read 5M3I5M at 101: the insertion is anchored at 1-based reference pos 105.
    let mut sam = coord_builder(&[("chr1", 1000)]);
    sam.add(read().at("chr1", 101).cigar("5M3I5M").matching_ref(&fasta));
    let bam = sam.to_temp_indexed_bam().unwrap();

    // Without a VCF the insertion is reported.
    let plain = run(bam.path(), ref_file.path(), vec![]);
    let plain_rows = indel_rows(&plain);
    let all = all_indel(&plain_rows);
    assert_eq!(all.num_insertions, 1);
    assert_eq!(all.num_inserted_bases, 3);

    // A known variant at the insertion's anchor (chr1:105) excludes the insertion.
    let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
    vcf.add("chr1", 105, "A").alt("G");
    let vcf_file = vcf.to_temp_vcf().unwrap();

    let options = ErrorOptions { vcf: Some(vcf_file.path().to_path_buf()), ..opts() };
    let masked = run_opts(bam.path(), ref_file.path(), options);
    let masked_rows = indel_rows(&masked);
    let all = all_indel(&masked_rows);
    assert_eq!(all.num_insertions, 0, "an insertion anchored at a known-variant site is excluded");
}
