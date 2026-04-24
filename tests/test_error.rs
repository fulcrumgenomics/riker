mod helpers;

use std::path::PathBuf;

use helpers::{FastaBuilder, SamBuilder, coord_builder, read_metrics_tsv};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::{
    Flags, MappingQuality,
    cigar::{Op, op::Kind},
};
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OutputOptions};
use riker_lib::commands::error::{
    Error, ErrorOptions, IndelMetric, MismatchMetric, OverlappingMismatchMetric,
};

/// Build a simple FR pair where read1 is forward and read2 is reverse.
#[expect(clippy::too_many_arguments)]
fn make_pair(
    builder: &mut SamBuilder,
    name: &str,
    ref_id: usize,
    pos1: usize,
    pos2: usize,
    tlen: i32,
    seq1: &[u8],
    seq2: &[u8],
    quals: &[u8],
) {
    let read_len = seq1.len();
    let mapq = MappingQuality::new(60).unwrap();
    let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();

    // Read 1: forward
    let r1 = RecordBuf::builder()
        .set_name(name)
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(pos1).unwrap())
        .set_mapping_quality(mapq)
        .set_cigar(cigar.clone())
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(pos2).unwrap())
        .set_template_length(tlen)
        .set_sequence(Sequence::from(seq1.to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build();

    // Read 2: reverse
    let r2 = RecordBuf::builder()
        .set_name(name)
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
        )
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(pos2).unwrap())
        .set_mapping_quality(mapq)
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(pos1).unwrap())
        .set_template_length(-tlen)
        .set_sequence(Sequence::from(seq2.to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build();

    builder.add_record(r1);
    builder.add_record(r2);
}

/// Run the error command and return the output directory for assertions.
fn run_error(
    bam_path: &std::path::Path,
    fasta_path: &std::path::Path,
    stratify_by: Vec<String>,
) -> PathBuf {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");

    let cmd = Error {
        input: InputOptions { input: bam_path.to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        options: ErrorOptions {
            reference: fasta_path.to_path_buf(),
            vcf: None,
            intervals: None,
            min_mapq: 20,
            min_bq: 0, // set to 0 for tests so all bases pass
            include_duplicates: false,
            max_isize: 1000,
            picard_compat: false,
            stratify_by,
        },
    };

    cmd.execute().expect("error command should succeed");

    // Leak the tempdir so files persist for assertions
    let path = prefix.clone();
    std::mem::forget(dir);
    path
}

/// Create test reference: chr1 = 1000bp of repeating ACGTACGT...
fn test_reference() -> (tempfile::NamedTempFile, Vec<u8>) {
    let seq: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
    let fasta = FastaBuilder::new().add_contig("chr1", &seq).to_temp_fasta().unwrap();
    (fasta, seq)
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[test]
fn test_no_errors() {
    let (fasta, ref_seq) = test_reference();
    let seq = ref_seq[99..109].to_vec(); // 10bp matching reference at pos 100-109
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    make_pair(&mut builder, "read1", 0, 100, 120, 30, &seq, &ref_seq[119..129], &quals);
    let bam = builder.to_temp_indexed_bam().unwrap();

    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // Should have one row for "all" stratifier with zero errors
    assert!(!mm.is_empty());
    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.total_bases > 0);
    assert_eq!(all_row.error_bases, 0);
    assert_float_eq!(all_row.frac_error, 0.0, 1e-9);
}

#[test]
fn test_simple_mismatches() {
    let (fasta, ref_seq) = test_reference();

    // Create a read that matches reference except for 2 bases
    let mut seq = ref_seq[99..109].to_vec(); // pos 100-109 (1-based)
    seq[0] = b'T'; // mismatch at pos 100 (ref is 'A' at 0-based index 99 → 99%4=3 → 'T', wait)
    // ref at 0-based 99 is ref_seq[99] = b"ACGT"[99%4] = b"ACGT"[3] = b'T'
    // So seq[0] matching ref_seq[99] = T. Let me pick a position where ref is 'A'.
    // ref at 0-based 100 is ref_seq[100] = b"ACGT"[100%4] = b"ACGT"[0] = b'A'
    // So let's use pos 101 (1-based) through 110 (1-based)
    let mut seq = ref_seq[100..110].to_vec(); // matching reference
    seq[0] = b'G'; // mismatch at index 100, ref='A', read='G'
    seq[5] = b'T'; // mismatch at index 105, ref='C' (105%4=1), read='T'
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Only add an unpaired read to keep it simple
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap()) // 1-based
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 2);
    assert!((all_row.frac_error - 0.2).abs() < 1e-6);
}

#[test]
fn test_stratify_by_strand() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Forward read with 1 mismatch
    let mut seq_fwd = ref_seq[100..110].to_vec();
    seq_fwd[0] = b'G'; // mismatch
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let r_fwd = RecordBuf::builder()
        .set_name("read_fwd")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(seq_fwd))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Reverse read with 0 mismatches
    let seq_rev = ref_seq[200..210].to_vec();
    let r_rev = RecordBuf::builder()
        .set_name("read_rev")
        .set_flags(Flags::REVERSE_COMPLEMENTED)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(201).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq_rev))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(r_fwd);
    builder.add_record(r_rev);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["strand".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let fwd = mm.iter().find(|r| r.stratifier == "strand" && r.covariate == "+").unwrap();
    let rev = mm.iter().find(|r| r.stratifier == "strand" && r.covariate == "-").unwrap();

    assert_eq!(fwd.total_bases, 10);
    assert_eq!(fwd.error_bases, 1);
    assert_eq!(rev.total_bases, 10);
    assert_eq!(rev.error_bases, 0);
}

#[test]
fn test_insertion_detection() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 13]; // 10M + 3I = 13 read bases

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read with 10M3I (10 aligned + 3 inserted bases)
    // Sequence: 10 bases matching ref + 3 inserted bases
    let mut seq = ref_seq[100..110].to_vec();
    seq.extend_from_slice(b"AAA"); // 3bp insertion
    let cigar: Cigar =
        [Op::new(Kind::Match, 10), Op::new(Kind::Insertion, 3)].into_iter().collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 10); // only aligned bases count for total
    assert_eq!(all_row.num_insertions, 1);
    assert_eq!(all_row.num_inserted_bases, 3);
    assert_eq!(all_row.num_deletions, 0);
}

#[test]
fn test_deletion_detection() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read with 5M2D5M (5 match, 2 deleted, 5 match)
    let mut seq = ref_seq[100..105].to_vec(); // first 5 bases
    seq.extend_from_slice(&ref_seq[107..112]); // skip 2, next 5 bases
    let cigar: Cigar =
        [Op::new(Kind::Match, 5), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 5)]
            .into_iter()
            .collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 10); // 5 + 5 aligned bases
    assert_eq!(all_row.num_deletions, 1);
    assert_eq!(all_row.num_deleted_bases, 2);
    assert_eq!(all_row.num_insertions, 0);
}

#[test]
fn test_min_mapq_filter() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read with MAPQ 10 (below default min of 20)
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(10).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // All bases should be filtered out — empty or zero
    if let Some(all_row) = mm.iter().find(|r| r.stratifier == "all") {
        assert_eq!(all_row.total_bases, 0);
    }
    // Or the file may have no rows at all, which is also valid
}

#[test]
fn test_duplicate_exclusion() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Duplicate read
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read_dup")
        .set_flags(Flags::DUPLICATE)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(seq.clone()))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Non-duplicate read with 1 mismatch
    let mut seq2 = ref_seq[100..110].to_vec();
    seq2[0] = b'G';
    let record2 = RecordBuf::builder()
        .set_name("read_ok")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq2))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(record);
    builder.add_record(record2);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    // Only the non-duplicate read should be counted
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 1);
}

#[test]
fn test_all_group_always_present() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    // User specifies only "bq" but "all" should still be present
    let prefix = run_error(bam.path(), fasta.path(), vec!["bq".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // Should have both "all" and "bq" rows
    assert!(mm.iter().any(|r| r.stratifier == "all"));
    assert!(mm.iter().any(|r| r.stratifier == "bq"));
}

#[test]
fn test_composite_stratifier() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["strand,mapq".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // Should have rows for the composite "strand,mapq" group
    let composite_rows: Vec<_> = mm.iter().filter(|r| r.stratifier == "strand,mapq").collect();
    assert!(!composite_rows.is_empty());
    // The covariate should be comma-separated, e.g. "+,60"
    assert!(composite_rows.iter().any(|r| r.covariate.contains(',')));
}

#[test]
fn test_overlapping_reads_mismatching_ref_and_mate() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Two overlapping reads at positions 101-110 and 106-115
    // Read 1: all matching reference
    let seq1 = ref_seq[100..110].to_vec();
    // Read 2: matches reference except pos 107 (0-based) which is in overlap region
    let mut seq2 = ref_seq[105..115].to_vec();
    seq2[2] = b'N'; // mismatch at position 107, but N won't count
    // Actually, N gets filtered. Let's use a real mismatch.
    let mut seq2 = ref_seq[105..115].to_vec();
    // ref at 0-based 107 = b"ACGT"[107%4] = b"ACGT"[3] = b'T'
    seq2[2] = b'A'; // mismatch at 0-based pos 107, ref='T', read='A'

    make_pair(
        &mut builder,
        "pair1",
        0,
        101, // read1 pos (1-based)
        106, // read2 pos (1-based)
        15,  // tlen
        &seq1,
        &seq2,
        &quals,
    );

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

    // We should have some overlap data
    if let Some(all_row) = ov.iter().find(|r| r.stratifier == "all") {
        // There should be overlapping bases examined
        assert!(all_row.overlapping_read_bases > 0, "Expected overlapping bases");
    }
}

#[test]
fn test_three_output_files_always_created() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    // All three output files should exist
    assert!(prefix.with_file_name("out.error-mismatch.txt").exists());
    assert!(prefix.with_file_name("out.error-overlap.txt").exists());
    assert!(prefix.with_file_name("out.error-indel.txt").exists());
}

// ─── Helper for running with custom min_bq ──────────────────────────────────

/// Run the error command with a custom `min_bq` threshold.
fn run_error_with_bq(
    bam_path: &std::path::Path,
    fasta_path: &std::path::Path,
    stratify_by: Vec<String>,
    min_bq: u8,
) -> PathBuf {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");

    let cmd = Error {
        input: InputOptions { input: bam_path.to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        options: ErrorOptions {
            reference: fasta_path.to_path_buf(),
            vcf: None,
            intervals: None,
            min_mapq: 20,
            min_bq,
            include_duplicates: false,
            max_isize: 1000,
            picard_compat: false,
            stratify_by,
        },
    };

    cmd.execute().expect("error command should succeed");
    let path = prefix.clone();
    std::mem::forget(dir);
    path
}

// ─── Core correctness tests ─────────────────────────────────────────────────

#[test]
fn test_min_bq_filter() {
    let (fasta, ref_seq) = test_reference();

    // Create a read where bases 0..5 have quality 10 and bases 5..10 have quality 30
    let seq = ref_seq[100..110].to_vec(); // all matching reference
    let mut quals = vec![10u8; 10];
    for q in quals.iter_mut().skip(5) {
        *q = 30;
    }

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error_with_bq(bam.path(), fasta.path(), vec![], 20);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    // Only the 5 high-quality bases should be counted
    assert_eq!(all_row.total_bases, 5);
}

#[test]
fn test_n_bases_excluded() {
    let (fasta, ref_seq) = test_reference();

    // Build a 10bp read: first 4 match ref, then N, then 4 match ref, then N
    // Aligned at position 101 (1-based), so 0-based indices 100..110
    let mut seq = ref_seq[100..110].to_vec();
    seq[4] = b'N';
    seq[9] = b'N';
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    // 2 N bases should be excluded, leaving 8 valid bases
    assert_eq!(all_row.total_bases, 8);
    assert_eq!(all_row.error_bases, 0);
}

#[test]
fn test_soft_clip_excluded() {
    let (fasta, ref_seq) = test_reference();

    // 10bp read with CIGAR 3S7M: 3 soft-clipped + 7 matched
    // The 7 matched bases align starting at position 101 (1-based)
    let mut seq = vec![b'A'; 3]; // soft-clipped bases (don't matter)
    seq.extend_from_slice(&ref_seq[100..107]); // 7 matching ref bases
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 7)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 7);
}

#[test]
fn test_secondary_supplementary_excluded() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();

    // Primary read with 1 mismatch
    let mut primary_seq = ref_seq[100..110].to_vec();
    primary_seq[0] = b'G'; // mismatch at ref='A'
    let primary = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(primary_seq))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Secondary read with 1 mismatch
    let mut secondary_seq = ref_seq[100..110].to_vec();
    secondary_seq[1] = b'T'; // mismatch
    let secondary = RecordBuf::builder()
        .set_name("read2")
        .set_flags(Flags::SECONDARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(secondary_seq))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Supplementary read with 1 mismatch
    let mut supplementary_seq = ref_seq[100..110].to_vec();
    supplementary_seq[2] = b'T'; // mismatch
    let supplementary = RecordBuf::builder()
        .set_name("read3")
        .set_flags(Flags::SUPPLEMENTARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(supplementary_seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(primary);
    builder.add_record(secondary);
    builder.add_record(supplementary);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    // Only the primary read should be counted
    assert_eq!(all_row.total_bases, 10);
    assert_eq!(all_row.error_bases, 1);
}

#[test]
fn test_mixed_insertion_and_deletion() {
    let (fasta, ref_seq) = test_reference();

    // CIGAR: 5M2I5M3D5M
    // Read bases: 5 (match) + 2 (ins) + 5 (match) + 5 (match) = 17 read bases
    // Ref consumed: 5 + 5 + 3 + 5 = 18 ref bases
    // Aligned bases: 5 + 5 + 5 = 15
    let mut seq = ref_seq[100..105].to_vec(); // first 5M
    seq.extend_from_slice(b"AA"); // 2I
    seq.extend_from_slice(&ref_seq[105..110]); // second 5M
    // After the 3D, reference skips to 113; next 5M covers ref 113..118
    seq.extend_from_slice(&ref_seq[113..118]); // third 5M
    let quals = vec![30u8; 17];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [
        Op::new(Kind::Match, 5),
        Op::new(Kind::Insertion, 2),
        Op::new(Kind::Match, 5),
        Op::new(Kind::Deletion, 3),
        Op::new(Kind::Match, 5),
    ]
    .into_iter()
    .collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 15);
    assert_eq!(all_row.num_insertions, 1);
    assert_eq!(all_row.num_inserted_bases, 2);
    assert_eq!(all_row.num_deletions, 1);
    assert_eq!(all_row.num_deleted_bases, 3);
}

#[test]
fn test_q_score_calculation() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();

    // Create 100 reads of 10bp each = 1000 total bases, with exactly 1 mismatch.
    // Stack all reads at the same position to stay within the 1000bp reference.
    for i in 0..100 {
        let mut seq = ref_seq[100..110].to_vec();
        if i == 0 {
            seq[0] = b'G'; // single mismatch in the first read (ref='A' at 0-based 100)
        }
        let record = RecordBuf::builder()
            .set_name(format!("read{i}").as_str())
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(101).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(quals.clone()))
            .build();
        builder.add_record(record);
    }

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 1000);
    assert_eq!(all_row.error_bases, 1);
    // q_score = -10 * log10(1/1000) = 30.0
    assert_float_eq!(all_row.q_score, 30.0, 0.01);
}

// ─── Stratifier correctness tests ───────────────────────────────────────────

#[test]
fn test_stratify_by_cycle() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();

    // Forward read with a mismatch at read offset 2 (cycle 3 for forward reads)
    let mut seq_fwd = ref_seq[100..110].to_vec();
    // ref at 0-based 102 = b"ACGT"[102%4] = b"ACGT"[2] = b'G'
    seq_fwd[2] = b'T'; // mismatch at cycle 3
    let fwd = RecordBuf::builder()
        .set_name("read_fwd")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(seq_fwd))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Reverse-strand read with a mismatch at read offset 7 (cycle 3 for reverse)
    // For reverse reads, cycle starts at read_len and counts down, so
    // read_offset 7 -> cycle = 10 - 7 = 3
    let mut seq_rev = ref_seq[200..210].to_vec();
    // ref at 0-based 207 = b"ACGT"[207%4] = b"ACGT"[3] = b'T'
    seq_rev[7] = b'A'; // mismatch at cycle 3 for reverse read
    let rev = RecordBuf::builder()
        .set_name("read_rev")
        .set_flags(Flags::REVERSE_COMPLEMENTED)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(201).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq_rev))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(fwd);
    builder.add_record(rev);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["cycle".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let cycle3 = mm.iter().find(|r| r.stratifier == "cycle" && r.covariate == "3").unwrap();
    // Both reads contribute a mismatch at cycle 3
    assert_eq!(cycle3.error_bases, 2);
}

#[test]
fn test_stratify_by_read_num() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read1 (forward) with 1 mismatch
    let mut seq1 = ref_seq[100..110].to_vec();
    seq1[0] = b'G'; // mismatch
    // Read2 (reverse) with 0 mismatches
    let seq2 = ref_seq[119..129].to_vec();

    make_pair(&mut builder, "pair1", 0, 101, 120, 29, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["read_num".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let r1 = mm.iter().find(|r| r.stratifier == "read_num" && r.covariate == "R1").unwrap();
    let r2 = mm.iter().find(|r| r.stratifier == "read_num" && r.covariate == "R2").unwrap();
    assert_eq!(r1.error_bases, 1);
    assert_eq!(r2.error_bases, 0);
}

#[test]
fn test_stratify_by_ref_base() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();

    // Forward read aligned at 101 (1-based). Reference pattern at 0-based 100..110:
    // A C G T A C G T A C
    let mut seq = ref_seq[100..110].to_vec();
    // Introduce mismatch at offset 0: ref='A', read='G'
    seq[0] = b'G';
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["ref_base".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // The mismatch at ref_base 'A' should show up
    let ref_a = mm.iter().find(|r| r.stratifier == "ref_base" && r.covariate == "A").unwrap();
    assert!(ref_a.error_bases >= 1);
    // Other ref bases should have 0 errors (we only introduced one mismatch at an 'A')
    let ref_c = mm.iter().find(|r| r.stratifier == "ref_base" && r.covariate == "C").unwrap();
    assert_eq!(ref_c.error_bases, 0);
}

#[test]
fn test_stratify_by_hp_len() {
    // Build a FASTA with a homopolymer run: ACGT then AAAA then CGTACGTACGT...
    let mut custom_seq: Vec<u8> = Vec::with_capacity(1000);
    // First 100 bases: normal ACGT pattern
    custom_seq.extend((0..100).map(|i| b"ACGT"[i % 4]));
    // Bases 100-103: AAAA (4-base homopolymer)
    custom_seq.extend_from_slice(b"AAAA");
    // Fill rest with ACGT pattern
    custom_seq.extend((104..1000).map(|i| b"ACGT"[i % 4]));

    let fasta = FastaBuilder::new().add_contig("chr1", &custom_seq).to_temp_fasta().unwrap();

    // Create a forward read spanning the homopolymer region
    // Read at positions 98..108 (0-based), which is 99-based=99 -> 1-based=99
    // Actually: positions 97..107 (0-based) = 1-based 98..107
    // This spans: ref[97..107] which includes the AAAA at positions 100-103
    let seq = custom_seq[97..107].to_vec(); // all matching reference
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(98).unwrap()) // 1-based
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["hp_len".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // The base after the 4A homopolymer (read offset where we see 4 preceding A's)
    // should have hp_len >= 3 (since the preceding identical bases form the hp run)
    // The hp_len stratifier should produce multiple covariate values including "0"
    // and some non-zero values for bases after the homopolymer
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
    let (fasta, ref_seq) = test_reference();

    // Read with 5M3I5M: 5 match + 3 inserted + 5 match = 13 read bases
    let mut seq = ref_seq[100..105].to_vec();
    seq.extend_from_slice(b"AAA"); // 3bp insertion
    seq.extend_from_slice(&ref_seq[105..110]);
    let quals = vec![30u8; 13];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar =
        [Op::new(Kind::Match, 5), Op::new(Kind::Insertion, 3), Op::new(Kind::Match, 5)]
            .into_iter()
            .collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["indel_len".to_string()]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    // Should have an indel_len=3 row for the 3bp insertion
    let indel_3 = indels.iter().find(|r| r.stratifier == "indel_len" && r.covariate == "3");
    assert!(indel_3.is_some(), "Expected indel_len=3 row for the 3bp insertion");
    let indel_3 = indel_3.unwrap();
    assert_eq!(indel_3.num_insertions, 1);

    // Mismatch metrics should have rows with indel_len=0 for aligned bases
    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();
    let mm_0 = mm.iter().find(|r| r.stratifier == "indel_len" && r.covariate == "0");
    assert!(mm_0.is_some(), "Expected indel_len=0 row for aligned bases in mismatch metrics");
    let mm_0 = mm_0.unwrap();
    assert_eq!(mm_0.total_bases, 10); // 5 + 5 aligned bases
}

// ─── Overlap detection tests ────────────────────────────────────────────────

#[test]
fn test_overlapping_reads_all_agree() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Two overlapping reads where both match reference in the overlap region
    let seq1 = ref_seq[100..110].to_vec(); // pos 101-110
    let seq2 = ref_seq[105..115].to_vec(); // pos 106-115, overlap at 106-110

    make_pair(&mut builder, "pair1", 0, 101, 106, 15, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0, "Expected overlapping bases");
    assert_eq!(all_row.bases_mismatching_ref_and_mate, 0);
    assert_eq!(all_row.bases_matching_mate_but_not_ref, 0);
    assert_eq!(all_row.bases_in_three_way_disagreement, 0);
}

#[test]
fn test_overlapping_reads_matching_mate_but_not_ref() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Overlapping pair where both reads have the SAME mismatch at position 107 (0-based)
    // ref at 107 = b"ACGT"[107%4] = b"ACGT"[3] = b'T'
    let mut seq1 = ref_seq[100..110].to_vec(); // pos 101-110
    seq1[7] = b'A'; // read offset 7 -> ref pos 107, ref='T', read='A'

    let mut seq2 = ref_seq[105..115].to_vec(); // pos 106-115
    seq2[2] = b'A'; // read offset 2 -> ref pos 107, ref='T', read='A'

    make_pair(&mut builder, "pair1", 0, 101, 106, 15, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    assert!(all_row.overlapping_read_bases > 0);
    // Both reads disagree with ref but agree with each other
    assert!(
        all_row.bases_matching_mate_but_not_ref > 0,
        "Expected bases_matching_mate_but_not_ref > 0, got {}",
        all_row.bases_matching_mate_but_not_ref
    );
}

#[test]
fn test_overlapping_reads_three_way() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Three-way disagreement at position 107 (0-based):
    // ref='T' (107%4=3), read='A', mate='C' — all different
    let mut seq1 = ref_seq[100..110].to_vec(); // pos 101-110
    seq1[7] = b'A'; // ref pos 107, ref='T', read1='A'

    let mut seq2 = ref_seq[105..115].to_vec(); // pos 106-115
    seq2[2] = b'C'; // ref pos 107, ref='T', read2='C'

    make_pair(&mut builder, "pair1", 0, 101, 106, 15, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

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
    // Strengthened version of the existing test with exact assertions
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read 1: all matching reference at pos 101-110
    let seq1 = ref_seq[100..110].to_vec();
    // Read 2: matches reference except pos 107 (0-based)
    // ref at 107 = b'T', we set read to 'A'
    // In the overlap region (106-110), only pos 107 mismatches in read2
    // Read1 agrees with ref at pos 107, so this is bases_mismatching_ref_and_mate
    let mut seq2 = ref_seq[105..115].to_vec();
    seq2[2] = b'A'; // mismatch at ref pos 107 (read2 offset 2)

    make_pair(&mut builder, "pair1", 0, 101, 106, 15, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

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
    // Build a FASTA with two contigs
    let seq1: Vec<u8> = (0..500).map(|i| b"ACGT"[i % 4]).collect();
    let seq2: Vec<u8> = (0..500).map(|i| b"ACGT"[i % 4]).collect();
    let fasta = FastaBuilder::new()
        .add_contig("chr1", &seq1)
        .add_contig("chr2", &seq2)
        .to_temp_fasta()
        .unwrap();

    let quals = vec![30u8; 10];
    let mut builder = coord_builder(&[("chr1", 500), ("chr2", 500)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();

    // Read on chr1
    let read_seq1 = seq1[100..110].to_vec();
    let r1 = RecordBuf::builder()
        .set_name("read_chr1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(read_seq1))
        .set_quality_scores(QualityScores::from(quals.clone()))
        .build();

    // Read on chr2
    let read_seq2 = seq2[100..110].to_vec();
    let r2 = RecordBuf::builder()
        .set_name("read_chr2")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(read_seq2))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(r1);
    builder.add_record(r2);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    // Both reads should be counted: 10 + 10 = 20 total bases
    assert_eq!(all_row.total_bases, 20);
}

#[test]
fn test_deletion_at_read_start() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // CIGAR: 2D10M — deletion before any aligned base, so no anchor
    let seq = ref_seq[102..112].to_vec(); // 10 bases after the 2bp deletion
    let cigar: Cigar = [Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 10)].into_iter().collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 10);
    // No anchor before the deletion, so it should be silently skipped
    assert_eq!(all_row.num_deletions, 0);
}

#[test]
fn test_insertion_at_read_start() {
    let (fasta, ref_seq) = test_reference();

    // CIGAR: 3I10M — insertion before any aligned base, so no anchor
    let mut seq = vec![b'A'; 3]; // 3 inserted bases
    seq.extend_from_slice(&ref_seq[100..110]); // 10 aligned bases
    let quals = vec![30u8; 13];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar =
        [Op::new(Kind::Insertion, 3), Op::new(Kind::Match, 10)].into_iter().collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    // No anchor before the insertion, so it should be skipped
    assert_eq!(all_row.num_insertions, 0);
}

#[test]
fn test_stratifier_parse_error() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let seq = ref_seq[100..110].to_vec();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");

    let cmd = Error {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix },
        options: ErrorOptions {
            reference: fasta.path().to_path_buf(),
            vcf: None,
            intervals: None,
            min_mapq: 20,
            min_bq: 0,
            include_duplicates: false,
            max_isize: 1000,
            picard_compat: false,
            stratify_by: vec!["invalid_name".to_string()],
        },
    };

    let result = cmd.execute();
    assert!(result.is_err(), "Expected error for invalid stratifier name");
}

// ─── Additional stratifier tests ────────────────────────────────────────────

#[test]
fn test_gc_stratification() {
    let (fasta, ref_seq) = test_reference();

    // Create a read with a known GC content.
    // Use ref_seq[100..110] = A C G T A C G T A C
    // GC count = 4 (C,G,C,G), total = 10 -> GC% = 40 (rounded = (4*100+5)/10 = 40)
    let seq = ref_seq[100..110].to_vec();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["gc".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    // Sequence ACGTACGTAC has 5 GC bases (3 C's + 2 G's) out of 10
    // GC% = (5*100 + 5)/10 = 50 (rounded)
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
    let (fasta, ref_seq) = test_reference();

    // Forward read at pos 101 (1-based), ref[100..110] = A C G T A C G T A C
    // Read matches ref exactly.
    // pre_dinuc = previous read base + current ref base (in sequencing order)
    // At read offset 1: prev_base=A (offset 0), ref_base=C (pos 101 0-based) -> "AC"
    let seq = ref_seq[100..110].to_vec();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["pre_dinuc".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

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
    let (fasta, ref_seq) = test_reference();

    // Forward read at pos 101 (1-based), ref[100..110] = A C G T A C G T A C
    // context_3bp = prev_read_base + ref_base + next_read_base
    // At read offset 1: prev=A(offset 0), ref=C(pos 101), next=G(offset 2) -> "ACG"
    let seq = ref_seq[100..110].to_vec();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec!["context_3bp".to_string()]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

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
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Create a single read that classify_overlap will put into the Buffer:
    // - paired, both mapped, same contig, tlen < 2*read_len, starts before mate
    // - BUT the mate is not in the BAM
    let mut seq = ref_seq[100..110].to_vec(); // pos 101-110
    seq[5] = b'G'; // mismatch at ref pos 105, ref='C' (105%4=1 -> 'C'), read='G'

    let mapq = MappingQuality::new(60).unwrap();
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let r1 = RecordBuf::builder()
        .set_name("orphan")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(mapq)
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(106).unwrap()) // mate starts later
        .set_template_length(15) // < 2*10, so overlap expected
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();

    builder.add_record(r1);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();

    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(all_row.total_bases, 10, "Orphaned read's 10 bases should be counted");
    assert_eq!(all_row.error_bases, 1, "Orphaned read's mismatch should be counted");
}

// ─── Additional tests ────────────────────────────────────────────────────────

/// An insertion where the first inserted base has BQ below min_bq should cause the entire
/// insertion to be skipped, while aligned bases are still counted.
#[test]
fn test_insertion_low_bq_first_base_excluded() {
    let (fasta, ref_seq) = test_reference();

    // CIGAR: 30M2I44M = 76 read bases, 74 aligned bases
    // Position 101 (1-based), so ref[100..174]
    let mut seq = ref_seq[100..130].to_vec(); // 30M
    seq.extend_from_slice(b"AA"); // 2I
    seq.extend_from_slice(&ref_seq[130..174]); // 44M
    assert_eq!(seq.len(), 76);

    // All bases get BQ=30 except the first inserted base (read offset 30) gets BQ=5
    let mut quals = vec![30u8; 76];
    quals[30] = 5;

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar =
        [Op::new(Kind::Match, 30), Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 44)]
            .into_iter()
            .collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error_with_bq(bam.path(), fasta.path(), vec![], 20);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    // The insertion should be entirely skipped because the first inserted base has BQ < min_bq
    assert_eq!(
        all_row.num_insertions, 0,
        "Insertion should be skipped when first base BQ < min_bq"
    );
    assert_eq!(all_row.num_inserted_bases, 0);
    // Aligned bases should still be counted (74 total, minus any with BQ < 20, but all aligned
    // bases have BQ=30 so all 74 should pass)
    assert_eq!(all_row.total_bases, 74);
}

/// When both reads in an overlapping pair have a mismatch at the same position (both differ from
/// reference but agree with each other), `bases_matching_mate_but_not_ref` should be > 0 and
/// `overlapping_read_bases` should count both reads' contributions to the overlap.
#[test]
fn test_overlap_double_counted() {
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read1: pos 101-110 (1-based). Mismatch at pos 107 (0-based 106).
    // ref at 0-based 106 = b"ACGT"[106%4] = b"ACGT"[2] = b'G'
    let mut seq1 = ref_seq[100..110].to_vec();
    seq1[6] = b'T'; // mismatch at ref pos 107 (0-based 106), ref='G', read='T'

    // Read2: pos 106-115 (1-based). Overlap region is 106-110. Mismatch at same position.
    // ref at 0-based 106 = 'G'. Read offset for pos 107 in read2 = 107-106 = 1
    let mut seq2 = ref_seq[105..115].to_vec();
    seq2[1] = b'T'; // mismatch at ref pos 107 (0-based 106), ref='G', read='T'

    make_pair(&mut builder, "pair1", 0, 101, 106, 15, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

    let all_row = ov.iter().find(|r| r.stratifier == "all").unwrap();
    // The overlap region is 106-110 (5 positions). Both reads are counted in the overlap,
    // so overlapping_read_bases should be 2 * 5 = 10.
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
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Pair 1: small insert size (tlen=15, overlapping)
    let small_fwd = ref_seq[100..110].to_vec();
    let small_rev = ref_seq[105..115].to_vec();
    make_pair(&mut builder, "small_pair", 0, 101, 106, 15, &small_fwd, &small_rev, &quals);

    // Pair 2: large insert size (tlen=500)
    let large_fwd = ref_seq[100..110].to_vec();
    let large_rev = ref_seq[590..600].to_vec();
    make_pair(&mut builder, "large_pair", 0, 101, 591, 500, &large_fwd, &large_rev, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");

    let cmd = Error {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        options: ErrorOptions {
            reference: fasta.path().to_path_buf(),
            vcf: None,
            intervals: None,
            min_mapq: 20,
            min_bq: 0,
            include_duplicates: false,
            max_isize: 100,
            picard_compat: false,
            stratify_by: vec!["all".into(), "isize".into()],
        },
    };

    cmd.execute().expect("error command should succeed");
    let path = prefix.clone();
    std::mem::forget(dir);

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&path.with_file_name("out.error-mismatch.txt")).unwrap();

    // The "all" row should count bases from BOTH pairs (4 reads x 10 bases = 40)
    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
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
    let (fasta, ref_seq) = test_reference();
    let quals = vec![30u8; 10];

    let mut builder = coord_builder(&[("chr1", 1000)]);

    // Read1 at pos 101-110, read2 at pos 200-209. No overlap.
    let seq1 = ref_seq[100..110].to_vec();
    let seq2 = ref_seq[199..209].to_vec();
    make_pair(&mut builder, "pair1", 0, 101, 200, 109, &seq1, &seq2, &quals);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let ov: Vec<OverlappingMismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-overlap.txt")).unwrap();

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
    let (fasta, ref_seq) = test_reference();

    // CIGAR: 3I50M — insertion before any aligned base, so no anchor
    let mut seq = vec![b'A'; 3]; // 3 inserted bases
    seq.extend_from_slice(&ref_seq[100..150]); // 50 aligned bases
    let quals = vec![30u8; 53];

    let mut builder = coord_builder(&[("chr1", 1000)]);
    let cigar: Cigar =
        [Op::new(Kind::Insertion, 3), Op::new(Kind::Match, 50)].into_iter().collect();

    let record = RecordBuf::builder()
        .set_name("read1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(101).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(quals))
        .build();
    builder.add_record(record);

    let bam = builder.to_temp_indexed_bam().unwrap();
    let prefix = run_error(bam.path(), fasta.path(), vec![]);

    let indels: Vec<IndelMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-indel.txt")).unwrap();

    let all_row = indels.iter().find(|r| r.stratifier == "all").unwrap();
    // No anchor before the insertion, so it should be skipped
    assert_eq!(
        all_row.num_insertions, 0,
        "Insertion at read start should have no anchor and be skipped"
    );
    // The 50 aligned bases should still be counted
    assert_eq!(all_row.total_bases, 50);

    // Verify mismatch metrics also count the aligned bases
    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();
    let mm_all = mm.iter().find(|r| r.stratifier == "all").unwrap();
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
    let chr1 = vec![b'A'; 200];
    let chr2 = vec![b'G'; 200];
    let fasta = FastaBuilder::new()
        .add_contig("chr1", &chr1)
        .add_contig("chr2", &chr2)
        .to_temp_fasta()
        .unwrap();

    let quals = vec![30u8; 50];
    let seq_all_a = vec![b'A'; 50];
    let mut builder = coord_builder(&[("chr1", 200), ("chr2", 200)]);
    // read1 at 1-based pos 1, mate at pos 30 (0-based 29) — mate_pos falls inside
    // read1's [1, 50] span so the buffer keeps it pending.
    make_pair(&mut builder, "read", 0, 1, 30, 79, &seq_all_a, &seq_all_a, &quals);
    let bam = builder.to_temp_indexed_bam().unwrap();

    // BED: chr1 partial interval ending *before* read1's mate_pos (29); chr2 full.
    let dir = tempfile::tempdir().unwrap();
    let bed_path = dir.path().join("regions.bed");
    std::fs::write(&bed_path, "chr1\t0\t15\nchr2\t0\t200\n").unwrap();

    let prefix = dir.path().join("out");
    let cmd = Error {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        options: ErrorOptions {
            reference: fasta.path().to_path_buf(),
            vcf: None,
            intervals: Some(bed_path),
            min_mapq: 20,
            min_bq: 0,
            include_duplicates: false,
            max_isize: 1000,
            picard_compat: false,
            stratify_by: vec![],
        },
    };
    cmd.execute().expect("error command should succeed");

    let mm: Vec<MismatchMetric> =
        read_metrics_tsv(&prefix.with_file_name("out.error-mismatch.txt")).unwrap();
    let all_row = mm.iter().find(|r| r.stratifier == "all").unwrap();
    assert_eq!(
        all_row.error_bases, 0,
        "orphan on chr1 must not be compared against chr2's reference bases"
    );
}
