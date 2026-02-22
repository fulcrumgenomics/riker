mod helpers;

use anyhow::Result;
use helpers::{SamBuilder, read_metrics_tsv};
use riker_lib::collector::Collector;
use riker_lib::commands::alignment::{
    AlignmentCollector, AlignmentOptions, AlignmentSummaryMetric, METRICS_SUFFIX,
};
use riker_lib::sam::alignment_reader::AlignmentReader;
use tempfile::NamedTempFile;

// ─── Helper ───────────────────────────────────────────────────────────────────

fn run_alignment(bam: &std::path::Path) -> Result<Vec<AlignmentSummaryMetric>> {
    let prefix = NamedTempFile::with_suffix(".alignment")?;
    let prefix_path = prefix.path().to_path_buf();
    let metrics_path =
        std::path::PathBuf::from(format!("{}{METRICS_SUFFIX}", prefix_path.display()));

    let mut collector =
        AlignmentCollector::new(bam, &prefix_path, None, &AlignmentOptions::default());

    let (mut reader, header) = AlignmentReader::new(bam, None)?;
    collector.initialize(&header)?;
    for result in reader.record_bufs(&header) {
        collector.accept(&result?, &header)?;
    }
    collector.finish()?;

    read_metrics_tsv::<AlignmentSummaryMetric>(&metrics_path)
}

fn row<'a>(metrics: &'a [AlignmentSummaryMetric], category: &str) -> &'a AlignmentSummaryMetric {
    metrics
        .iter()
        .find(|m| m.category == category)
        .unwrap_or_else(|| panic!("no row for category '{category}'"))
}

// ─── Basic paired metrics ─────────────────────────────────────────────────────

#[test]
fn test_basic_paired_reads_first_and_second() -> Result<()> {
    let mut builder = SamBuilder::new();
    // 5 FR pairs, all mapped at MAPQ=60, read_len=100.
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
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let first = row(&metrics, "read1");
    let second = row(&metrics, "read2");
    let pair = row(&metrics, "pair");

    // 5 reads per category; 10 for pair.
    assert_eq!(first.total_reads, 5);
    assert_eq!(second.total_reads, 5);
    assert_eq!(pair.total_reads, 10);

    // All mapped.
    assert_eq!(first.aligned_reads, 5);
    assert_eq!(second.aligned_reads, 5);
    assert_eq!(pair.aligned_reads, 10);

    // All properly paired → aligned_reads_in_pairs should equal aligned_reads.
    assert_eq!(first.aligned_reads_in_pairs, 5);
    assert_eq!(second.aligned_reads_in_pairs, 5);

    // frac_aligned = 1.0 (all PF reads are aligned).
    assert!((first.frac_aligned - 1.0).abs() < 1e-5);

    // Strand balance: all R1 forward, all R2 reverse.
    assert!((first.strand_balance - 1.0).abs() < 1e-5);
    assert!((second.strand_balance - 0.0).abs() < 1e-5);

    // PAIR strand balance: 5 positive + 5 negative = 0.5.
    assert!((pair.strand_balance - 0.5).abs() < 1e-5);

    // Read lengths: all 100.
    assert!((first.mean_read_length - 100.0).abs() < 1e-9);
    assert_eq!(first.min_read_length, 100);
    assert_eq!(first.max_read_length, 100);

    Ok(())
}

#[test]
fn test_no_paired_data_produces_unpaired_row() -> Result<()> {
    // Empty BAM → single read1 row with all zeros.
    let builder = SamBuilder::new();
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    assert_eq!(metrics.len(), 1);
    assert_eq!(metrics[0].category, "read1");
    assert_eq!(metrics[0].total_reads, 0);
    assert_eq!(metrics[0].aligned_reads, 0);
    Ok(())
}

// ─── read1 (unpaired) category ────────────────────────────────────────────────

#[test]
fn test_unpaired_reads_category() -> Result<()> {
    let mut builder = SamBuilder::new();
    for i in 0..4 {
        builder.add_unpaired(
            &format!("r{i}"),
            0,
            1000 + i * 100,
            60,
            75,
            i % 2 == 0, // alternating strands
            false,
            false,
            None,
        );
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    // Only a read1 row when no paired reads.
    assert_eq!(metrics.len(), 1);
    let unpaired = &metrics[0];
    assert_eq!(unpaired.category, "read1");
    assert_eq!(unpaired.total_reads, 4);
    assert_eq!(unpaired.aligned_reads, 4);
    assert!((unpaired.frac_aligned - 1.0).abs() < 1e-5);
    // 2 forward, 2 reverse → strand balance = 0.5.
    assert!((unpaired.strand_balance - 0.5).abs() < 1e-5);
    Ok(())
}

// ─── QC-fail reads ────────────────────────────────────────────────────────────

#[test]
fn test_qc_fail_counted_in_total_not_stats() -> Result<()> {
    let mut builder = SamBuilder::new();
    // 3 normal pairs + 2 QC-fail pairs.
    for i in 0..3 {
        builder.add_pair(&format!("ok{i}"), 0, 100, 300, 200, 60, 100, false, false);
    }
    for i in 0..2 {
        builder.add_pair(&format!("fail{i}"), 0, 100, 300, 200, 60, 100, false, true);
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let first = row(&metrics, "read1");
    // 5 total (3 ok + 2 fail), but only 3 aligned (PF only).
    assert_eq!(first.total_reads, 5);
    assert_eq!(first.aligned_reads, 3);
    Ok(())
}

// ─── HQ threshold ─────────────────────────────────────────────────────────────

#[test]
fn test_hq_threshold_filters_low_mapq() -> Result<()> {
    let mut builder = SamBuilder::new();
    // 3 reads at MAPQ=60 (HQ) + 2 reads at MAPQ=10 (not HQ for min_mapq=20).
    for i in 0..3 {
        builder.add_unpaired(&format!("hq{i}"), 0, 1000, 60, 100, false, false, false, None);
    }
    for i in 0..2 {
        builder.add_unpaired(&format!("lq{i}"), 0, 2000, 10, 100, false, false, false, None);
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.aligned_reads, 5);
    assert_eq!(unpaired.hq_aligned_reads, 3); // only MAPQ≥20 counted as HQ
    Ok(())
}

// ─── Mismatch rate (NM tag) ────────────────────────────────────────────────────

#[test]
fn test_mismatch_rate_from_nm_tag() -> Result<()> {
    let mut builder = SamBuilder::new();
    // 2 pairs, NM=2 for each R1, NM=3 for each R2.
    // read_len=100 → 100 aligned bases per read.
    // Total NM = 2+3+2+3 = 10 over 400 aligned bases → mismatch_rate ≈ 0.025.
    // HQ (MAPQ=60): same reads → hq_mismatch_rate ≈ 0.025.
    for i in 0..2 {
        builder.add_pair_with_nm(&format!("r{i}"), 0, 100, 300, 200, 60, 100, 2, 3);
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    // 4 reads × 100 aligned bases = 400 total aligned bases.
    // Total NM = 2+3+2+3 = 10 → mismatch_rate = 10/400 = 0.025.
    assert!((pair.mismatch_rate - 0.025).abs() < 1e-6, "mismatch_rate={}", pair.mismatch_rate);
    // All reads are HQ → hq_mismatch_rate same.
    assert!(
        (pair.hq_mismatch_rate - 0.025).abs() < 1e-6,
        "hq_mismatch_rate={}",
        pair.hq_mismatch_rate
    );
    Ok(())
}

#[test]
fn test_hq_median_mismatches() -> Result<()> {
    let mut builder = SamBuilder::new();
    // Unpaired reads with NM values 1, 2, 3, 4, 5.
    for nm in [1u8, 2, 3, 4, 5] {
        builder.add_unpaired(&format!("r{nm}"), 0, 1000, 60, 100, false, false, false, Some(nm));
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // Median of {1,2,3,4,5} = 3.0.
    assert!((unpaired.hq_median_mismatches - 3.0).abs() < 0.1, "{}", unpaired.hq_median_mismatches);
    Ok(())
}

// ─── Indel rate ───────────────────────────────────────────────────────────────

#[test]
fn test_indel_rate_from_cigar() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // Build a read with 1 insertion and 1 deletion → 2 indel events, 98 M bases.
    let cigar: Cigar = [
        Op::new(Kind::Match, 49),
        Op::new(Kind::Insertion, 2),
        Op::new(Kind::Match, 49),
        Op::new(Kind::Deletion, 1),
    ]
    .into_iter()
    .collect();

    let flags = Flags::empty();
    let seq: Sequence = vec![b'A'; 100].into(); // 49 + 2 + 49 = 100 read bases
    let qual = QualityScores::from(vec![30u8; 100]);

    let record = noodles::sam::alignment::RecordBuf::builder()
        .set_name("indel_read")
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(100).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(seq)
        .set_quality_scores(qual)
        .build();

    let mut sb = SamBuilder::new();
    sb.add_record(record);
    let bam = sb.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // 2 indel events / 98 aligned bases ≈ 0.020408...
    assert!((unpaired.indel_rate - 2.0 / 98.0).abs() < 1e-6, "indel_rate={}", unpaired.indel_rate);
    Ok(())
}

// ─── Chimera detection ────────────────────────────────────────────────────────

#[test]
fn test_chimera_different_contigs() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    let mut builder = SamBuilder::with_contigs(&[
        ("chr1".to_string(), 249_250_621),
        ("chr2".to_string(), 243_199_373),
    ]);

    let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
    let seq: Sequence = vec![b'A'; 100].into();
    let qual = QualityScores::from(vec![30u8; 100]);

    // R1 maps to chr1, R2 maps to chr2 → chimeric pair.
    let r1 = noodles::sam::alignment::RecordBuf::builder()
        .set_name("chimera1")
        .set_flags(Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED | Flags::FIRST_SEGMENT)
        .set_reference_sequence_id(0) // chr1
        .set_alignment_start(Position::new(1000).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_mate_reference_sequence_id(1) // chr2
        .set_mate_alignment_start(Position::new(2000).unwrap())
        .set_template_length(0)
        .set_sequence(seq.clone())
        .set_quality_scores(qual.clone())
        .build();

    let r2 = noodles::sam::alignment::RecordBuf::builder()
        .set_name("chimera1")
        .set_flags(Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED | Flags::LAST_SEGMENT)
        .set_reference_sequence_id(1) // chr2
        .set_alignment_start(Position::new(2000).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0) // chr1
        .set_mate_alignment_start(Position::new(1000).unwrap())
        .set_template_length(0)
        .set_sequence(seq)
        .set_quality_scores(qual)
        .build();

    builder.add_record(r1);
    builder.add_record(r2);

    // Also add 3 normal FR pairs (not chimeric) for context.
    for i in 0..3 {
        builder.add_pair(
            &format!("normal{i}"),
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
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    // 1 chimeric pair / denominator (1 chimeric + 3 normal = 4 pairs × 2 reads = 8 reads)
    // chimeras = 2 (both R1 and R2 are counted), chimeras_denominator ≥ 2.
    assert!(pair.frac_chimeras > 0.0, "frac_chimeras should be > 0");
    Ok(())
}

#[test]
fn test_chimera_large_insert() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    let mut builder = SamBuilder::new();
    let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
    let seq: Sequence = vec![b'A'; 100].into();
    let qual = QualityScores::from(vec![30u8; 100]);

    // Insert size of 200,000 (> default max 100,000) → chimeric.
    let r1 = noodles::sam::alignment::RecordBuf::builder()
        .set_name("large_insert")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1000).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar.clone())
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(200_100).unwrap())
        .set_template_length(200_000i32)
        .set_sequence(seq.clone())
        .set_quality_scores(qual.clone())
        .build();

    let r2 = noodles::sam::alignment::RecordBuf::builder()
        .set_name("large_insert")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(200_100).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(1000).unwrap())
        .set_template_length(-200_000i32)
        .set_sequence(seq)
        .set_quality_scores(qual)
        .build();

    builder.add_record(r1);
    builder.add_record(r2);

    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    assert!((pair.frac_chimeras - 1.0).abs() < 1e-5, "frac_chimeras={}", pair.frac_chimeras);
    Ok(())
}

// ─── Bad cycles ───────────────────────────────────────────────────────────────

#[test]
fn test_bad_cycles_all_n_at_same_position() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // 5 reads each with N at position 0 (cycle 1 on forward strand).
    // 5/5 = 100% no-call at cycle 1 → bad_cycles = 1.
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let mut builder = SamBuilder::new();
    for i in 0..5 {
        let mut bases = vec![b'A'; 10];
        bases[0] = b'N'; // N at cycle 1
        let seq: Sequence = bases.into();
        let qual = QualityScores::from(vec![30u8; 10]);
        let record = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("r{i}").as_str())
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_sequence(seq)
            .set_quality_scores(qual)
            .build();
        builder.add_record(record);
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.bad_cycles, 1);
    Ok(())
}

#[test]
fn test_bad_cycles_below_threshold_not_counted() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // 10 reads; only 4 have N at position 0 = 40% < 80% → bad_cycles = 0.
    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let mut builder = SamBuilder::new();
    for i in 0..10 {
        let mut bases = vec![b'A'; 10];
        if i < 4 {
            bases[0] = b'N';
        }
        let seq: Sequence = bases.into();
        let qual = QualityScores::from(vec![30u8; 10]);
        let record = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("r{i}").as_str())
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_sequence(seq)
            .set_quality_scores(qual)
            .build();
        builder.add_record(record);
    }
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.bad_cycles, 0);
    Ok(())
}

// ─── Soft-clip fractions ──────────────────────────────────────────────────────

#[test]
fn test_softclip_fraction() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // 1 read: 10S 80M 10S  (100 bases total, 20 soft-clipped).
    let cigar: Cigar =
        [Op::new(Kind::SoftClip, 10), Op::new(Kind::Match, 80), Op::new(Kind::SoftClip, 10)]
            .into_iter()
            .collect();

    let seq: Sequence = vec![b'A'; 100].into();
    let qual = QualityScores::from(vec![30u8; 100]);
    let record = noodles::sam::alignment::RecordBuf::builder()
        .set_name("sc_read")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(100).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(seq)
        .set_quality_scores(qual)
        .build();

    let mut builder = SamBuilder::new();
    builder.add_record(record);
    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // 20 soft-clipped / 100 total bases = 0.2.
    assert!(
        (unpaired.frac_softclipped_reads - 0.2).abs() < 1e-5,
        "frac_softclipped_reads={}",
        unpaired.frac_softclipped_reads
    );
    // No hard clips.
    assert!((unpaired.frac_hardclipped_reads - 0.0).abs() < 1e-5);
    Ok(())
}

// ─── Reference validation ─────────────────────────────────────────────────────

#[test]
fn test_reference_validation_missing_contig_errors() -> Result<()> {
    use helpers::FastaBuilder;

    let mut builder = SamBuilder::with_contigs(&[
        ("chr1".to_string(), 249_250_621),
        ("chrZ".to_string(), 1_000_000), // not in reference
    ]);
    builder.add_pair("r1", 0, 100, 300, 200, 60, 100, false, false);
    let bam = builder.to_temp_bam()?;

    // Reference only has chr1.
    let reference = FastaBuilder::new().add_contig("chr1", &vec![b'A'; 1000]).to_temp_fasta()?;

    let prefix = NamedTempFile::with_suffix(".alignment")?;
    let prefix_path = prefix.path().to_path_buf();

    let mut collector = AlignmentCollector::new(
        bam.path(),
        &prefix_path,
        Some(reference.path().to_path_buf()),
        &AlignmentOptions::default(),
    );

    let (_, header) = AlignmentReader::new(bam.path(), None)?;
    let result = collector.initialize(&header);
    assert!(result.is_err(), "expected error for missing contig");
    assert!(result.unwrap_err().to_string().contains("chrZ"));
    Ok(())
}

// ─── Output file naming ───────────────────────────────────────────────────────

#[test]
fn test_output_file_created_with_correct_suffix() -> Result<()> {
    let mut builder = SamBuilder::new();
    builder.add_unpaired("r1", 0, 100, 60, 100, false, false, false, None);
    let bam = builder.to_temp_bam()?;

    let prefix = NamedTempFile::with_suffix(".align_test")?;
    let prefix_path = prefix.path().to_path_buf();
    let expected_metrics =
        std::path::PathBuf::from(format!("{}{METRICS_SUFFIX}", prefix_path.display()));

    let mut collector =
        AlignmentCollector::new(bam.path(), &prefix_path, None, &AlignmentOptions::default());
    let (mut reader, header) = AlignmentReader::new(bam.path(), None)?;
    collector.initialize(&header)?;
    for result in reader.record_bufs(&header) {
        collector.accept(&result?, &header)?;
    }
    collector.finish()?;

    assert!(expected_metrics.exists(), "metrics file not created at expected path");
    Ok(())
}

// ─── PAIR bad_cycles override ─────────────────────────────────────────────────

#[test]
fn test_pair_bad_cycles_is_sum_of_first_and_second() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // 5 R1 reads: all have N at cycle 1 (forward strand, position 0).
    // 5 R2 reads: all have N at cycle 2 (reverse strand read of len 10:
    //   position 0 → cycle = 10 - 0 = 10; position 9 → cycle = 1;
    //   to get cycle 2 on reverse strand we need N at position 8 of a 10-base read).
    // With 5/5 = 100% at cycle 1 (from R1) AND cycle X (from R2), bad_cycles(PAIR) = sum.

    let cigar: Cigar = [Op::new(Kind::Match, 10)].into_iter().collect();
    let mut builder = SamBuilder::new();

    for i in 0..5 {
        let pos1 = 100 + i * 200;
        let pos2 = 300 + i * 200;

        // R1 (forward): N at index 0 → cycle 1.
        let mut seq1 = vec![b'A'; 10];
        seq1[0] = b'N';

        // R2 (reverse): N at index 0 → cycle = read_len(10) - 0 = 10 on rev strand.
        // All R2 have N at cycle 10 (100% → bad).
        let mut seq2 = vec![b'A'; 10];
        seq2[0] = b'N';

        let r1 = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("r{i}").as_str())
            .set_flags(
                Flags::SEGMENTED
                    | Flags::PROPERLY_SEGMENTED
                    | Flags::MATE_REVERSE_COMPLEMENTED
                    | Flags::FIRST_SEGMENT,
            )
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos1).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(pos2).unwrap())
            .set_template_length(200)
            .set_sequence(Sequence::from(seq1))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build();

        let r2 = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("r{i}").as_str())
            .set_flags(
                Flags::SEGMENTED
                    | Flags::PROPERLY_SEGMENTED
                    | Flags::REVERSE_COMPLEMENTED
                    | Flags::LAST_SEGMENT,
            )
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos2).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(pos1).unwrap())
            .set_template_length(-200)
            .set_sequence(Sequence::from(seq2))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build();

        builder.add_record(r1);
        builder.add_record(r2);
    }

    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let first = row(&metrics, "read1");
    let second = row(&metrics, "read2");
    let pair = row(&metrics, "pair");

    // PAIR bad_cycles must equal FIRST + SECOND.
    assert_eq!(pair.bad_cycles, first.bad_cycles + second.bad_cycles);
    // Both R1 and R2 have 100% N at their respective cycle → each contributes 1 bad cycle.
    assert_eq!(first.bad_cycles, 1, "read1 bad_cycles");
    assert_eq!(second.bad_cycles, 1, "read2 bad_cycles");
    assert_eq!(pair.bad_cycles, 2, "pair bad_cycles");
    Ok(())
}

// ─── Improper pairs ───────────────────────────────────────────────────────────

#[test]
fn test_improper_pairs_counted() -> Result<()> {
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    let mut builder = SamBuilder::new();

    // 3 proper pairs (from add_pair which sets PROPERLY_SEGMENTED).
    for i in 0..3 {
        builder.add_pair(
            &format!("ok{i}"),
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

    // 2 improper pairs: SEGMENTED but NOT PROPERLY_SEGMENTED.
    let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
    let seq: Sequence = vec![b'A'; 100].into();
    let qual = QualityScores::from(vec![30u8; 100]);

    for i in 0..2 {
        let pos1 = 10_000 + i * 500;
        let pos2 = 10_200 + i * 500;
        let r1 = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("imp{i}").as_str())
            .set_flags(Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos1).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(pos2).unwrap())
            .set_template_length(200)
            .set_sequence(seq.clone())
            .set_quality_scores(qual.clone())
            .build();
        let r2 = noodles::sam::alignment::RecordBuf::builder()
            .set_name(format!("imp{i}").as_str())
            .set_flags(Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos2).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(pos1).unwrap())
            .set_template_length(-200)
            .set_sequence(seq.clone())
            .set_quality_scores(qual.clone())
            .build();
        builder.add_record(r1);
        builder.add_record(r2);
    }

    let bam = builder.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    // 4 improper reads (2 pairs × 2) / 10 total aligned = 0.4.
    assert_eq!(pair.reads_improperly_paired, 4, "reads_improperly_paired");
    assert!(
        (pair.frac_reads_improperly_paired - 0.4).abs() < 1e-5,
        "frac_improper={}",
        pair.frac_reads_improperly_paired
    );
    Ok(())
}
