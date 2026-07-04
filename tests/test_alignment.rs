mod helpers;

use anyhow::Result;
use helpers::read_metrics_tsv;
use riker_lib::assert_close;
use riker_lib::collector::Collector;
use riker_lib::commands::alignment::{
    AlignmentCollector, AlignmentOptions, AlignmentSummaryMetric, METRICS_SUFFIX,
};
use riker_lib::sam::alignment_reader::AlignmentReader;
use riker_lib::test_support::{FastaBuilder, SamBuilder, drive_collector, pair, read};
use tempfile::NamedTempFile;

// ─── Helpers ────────────────────────────────────────────────────────────────

fn run_alignment(bam: &std::path::Path) -> Result<Vec<AlignmentSummaryMetric>> {
    let prefix = NamedTempFile::with_suffix(".alignment")?;
    let prefix_path = prefix.path().to_path_buf();
    let metrics_path =
        std::path::PathBuf::from(format!("{}{METRICS_SUFFIX}", prefix_path.display()));

    let mut collector =
        AlignmentCollector::new(bam, &prefix_path, None, &AlignmentOptions::default());
    drive_collector(&mut collector, bam)?;

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
    // 5 FR pairs, all mapped at MAPQ 60, read length 100.
    let mut sam = SamBuilder::new();
    for i in 0..5 {
        sam.add(pair(sam.next_id()).at("chr1", 100 + i * 200, 300 + i * 200));
    }
    let bam = sam.to_temp_bam()?;
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
    assert_close!(first.frac_aligned, 1.0, 1e-5);

    // Strand balance: all R1 forward, all R2 reverse.
    assert_close!(first.strand_balance, 1.0, 1e-5);
    assert_close!(second.strand_balance, 0.0, 1e-5);

    // PAIR strand balance: 5 positive + 5 negative = 0.5.
    assert_close!(pair.strand_balance, 0.5, 1e-5);

    // Read lengths: all 100.
    assert_close!(first.mean_read_length, 100.0, 1e-9);
    assert_eq!(first.min_read_length, 100);
    assert_eq!(first.max_read_length, 100);

    Ok(())
}

#[test]
fn test_no_paired_data_produces_unpaired_row() -> Result<()> {
    // Empty BAM → single read1 row with all zeros.
    let sam = SamBuilder::new();
    let bam = sam.to_temp_bam()?;
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
    // 4 unpaired reads (length 75), alternating strand.
    let mut sam = SamBuilder::new();
    for i in 0..4 {
        let mut r = read().at("chr1", 1000 + i * 100).len(75);
        if i % 2 == 0 {
            r = r.reverse();
        }
        sam.add(r);
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    // Only a read1 row when no paired reads.
    assert_eq!(metrics.len(), 1);
    let unpaired = &metrics[0];
    assert_eq!(unpaired.category, "read1");
    assert_eq!(unpaired.total_reads, 4);
    assert_eq!(unpaired.aligned_reads, 4);
    assert_close!(unpaired.frac_aligned, 1.0, 1e-5);
    // 2 forward, 2 reverse → strand balance = 0.5.
    assert_close!(unpaired.strand_balance, 0.5, 1e-5);
    Ok(())
}

// ─── QC-fail reads ────────────────────────────────────────────────────────────

#[test]
fn test_qc_fail_counted_in_total_not_stats() -> Result<()> {
    // 3 normal pairs + 2 QC-fail pairs.
    let mut sam = SamBuilder::new();
    for _ in 0..3 {
        sam.add(pair(sam.next_id()).at("chr1", 100, 300));
    }
    for _ in 0..2 {
        sam.add(pair(sam.next_id()).at("chr1", 100, 300).qc_fail());
    }
    let bam = sam.to_temp_bam()?;
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
    // 3 reads at MAPQ 60 (HQ) + 2 reads at MAPQ 10 (below the min_mapq=20 threshold).
    let mut sam = SamBuilder::new();
    for _ in 0..3 {
        sam.add(read().at("chr1", 1000));
    }
    for _ in 0..2 {
        sam.add(read().at("chr1", 2000).mapq(10));
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.aligned_reads, 5);
    assert_eq!(unpaired.hq_aligned_reads, 3); // only MAPQ≥20 counted as HQ
    Ok(())
}

// ─── Mismatch rate (NM tag) ────────────────────────────────────────────────────

#[test]
fn test_mismatch_rate_from_nm_tag() -> Result<()> {
    // 2 pairs, NM=2 for each R1, NM=3 for each R2, read length 100.
    // Total NM = 2+3+2+3 = 10 over 400 aligned bases → mismatch_rate = 0.025.
    let mut sam = SamBuilder::new();
    for _ in 0..2 {
        sam.add(pair(sam.next_id()).at("chr1", 100, 300).r1(|r| r.nm(2)).r2(|r| r.nm(3)));
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    // 4 reads × 100 aligned bases = 400 total aligned bases; total NM = 10 → 10/400 = 0.025.
    assert_close!(pair.mismatch_rate, 0.025, 1e-6);
    // All reads are HQ → hq_mismatch_rate same.
    assert_close!(pair.hq_mismatch_rate, 0.025, 1e-6);
    Ok(())
}

#[test]
fn test_hq_median_mismatches() -> Result<()> {
    // Unpaired reads with NM values 1, 2, 3, 4, 5 → median 3.
    let mut sam = SamBuilder::new();
    for nm in [1u8, 2, 3, 4, 5] {
        sam.add(read().at("chr1", 1000).nm(nm));
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // Median of {1,2,3,4,5} = 3.0.
    assert_close!(unpaired.hq_median_mismatches, 3.0, 0.1);
    Ok(())
}

// ─── Indel rate ───────────────────────────────────────────────────────────────

#[test]
fn test_indel_rate_from_cigar() -> Result<()> {
    // 1 insertion (2 bp) + 1 deletion (1 bp) → 2 indel events over 98 aligned M bases.
    let mut sam = SamBuilder::new();
    sam.add(read().at("chr1", 100).cigar("49M2I49M1D"));
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // 2 indel events / 98 aligned bases ≈ 0.020408...
    assert_close!(unpaired.indel_rate, 2.0 / 98.0, 1e-6);
    Ok(())
}

// ─── Chimera detection ────────────────────────────────────────────────────────

#[test]
fn test_chimera_different_contigs() -> Result<()> {
    let mut sam = SamBuilder::with_contigs(&[
        ("chr1".to_string(), 249_250_621),
        ("chr2".to_string(), 243_199_373),
    ]);

    // R1 maps to chr1, R2 maps to chr2 → chimeric pair.
    sam.add(pair("chimera1").r1(|r| r.at("chr1", 1000)).r2(|r| r.at("chr2", 2000).reverse()));

    // Also add 3 normal FR pairs (not chimeric) for context.
    for i in 0..3 {
        sam.add(pair(sam.next_id()).at("chr1", 100 + i * 200, 300 + i * 200));
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    assert!(pair.frac_chimeras > 0.0, "frac_chimeras should be > 0");
    Ok(())
}

#[test]
fn test_chimera_large_insert() -> Result<()> {
    // Same contig, insert size ~199 kb (> default max 100 kb) → chimeric.
    let mut sam = SamBuilder::new();
    sam.add(pair("large_insert").at("chr1", 1000, 200_100));

    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    assert_close!(pair.frac_chimeras, 1.0, 1e-5);
    Ok(())
}

// ─── Bad cycles ───────────────────────────────────────────────────────────────

#[test]
fn test_bad_cycles_all_n_at_same_position() -> Result<()> {
    // 5 reads each with an N no-call at cycle 1 (forward-strand offset 0).
    // 5/5 = 100% no-call at cycle 1 → bad_cycles = 1.
    let mut sam = SamBuilder::new();
    for _ in 0..5 {
        sam.add(read().at("chr1", 100).cigar("10M").sub(0, b'N'));
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.bad_cycles, 1);
    Ok(())
}

#[test]
fn test_bad_cycles_below_threshold_not_counted() -> Result<()> {
    // 4 of 10 reads have an N at cycle 1 = 40% < 80% threshold → bad_cycles = 0.
    let mut sam = SamBuilder::new();
    for i in 0..10 {
        let mut r = read().at("chr1", 100).cigar("10M");
        if i < 4 {
            r = r.sub(0, b'N');
        }
        sam.add(r);
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    assert_eq!(unpaired.bad_cycles, 0);
    Ok(())
}

// ─── Soft-clip fractions ──────────────────────────────────────────────────────

#[test]
fn test_softclip_fraction() -> Result<()> {
    // 10S 80M 10S: 20 of 100 bases soft-clipped.
    let mut sam = SamBuilder::new();
    sam.add(read().at("chr1", 100).cigar("10S80M10S"));
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let unpaired = row(&metrics, "read1");
    // 20 soft-clipped / 100 total bases = 0.2.
    assert_close!(unpaired.frac_softclipped_reads, 0.2, 1e-5);
    // No hard clips.
    assert_close!(unpaired.frac_hardclipped_reads, 0.0, 1e-5);
    Ok(())
}

// ─── Reference validation ─────────────────────────────────────────────────────

#[test]
fn test_reference_validation_missing_contig_errors() -> Result<()> {
    // The header declares chrZ but the reference only has chr1 → validation error.
    let mut sam = SamBuilder::with_contigs(&[
        ("chr1".to_string(), 249_250_621),
        ("chrZ".to_string(), 1_000_000), // not in reference
    ]);
    sam.add(pair("r1").at("chr1", 100, 300));
    let bam = sam.to_temp_bam()?;

    // Reference only has chr1.
    let reference = NamedTempFile::with_suffix(".fa")?;
    FastaBuilder::new().contig("chr1", vec![b'A'; 1000]).write_fasta(reference.path())?;

    let prefix = NamedTempFile::with_suffix(".alignment")?;
    let prefix_path = prefix.path().to_path_buf();

    let mut collector = AlignmentCollector::new(
        bam.path(),
        &prefix_path,
        Some(reference.path().to_path_buf()),
        &AlignmentOptions::default(),
    );

    // This test needs the initialize() error itself, so it drives the collector by hand
    // rather than through drive_collector.
    let reader = AlignmentReader::open(bam.path(), None, 0)?;
    let result = collector.initialize(reader.header());
    assert!(result.is_err(), "expected error for missing contig");
    assert!(result.unwrap_err().to_string().contains("chrZ"));
    Ok(())
}

// ─── Output file naming ───────────────────────────────────────────────────────

#[test]
fn test_output_file_created_with_correct_suffix() -> Result<()> {
    let mut sam = SamBuilder::new();
    sam.add(read().at("chr1", 100));
    let bam = sam.to_temp_bam()?;

    let prefix = NamedTempFile::with_suffix(".align_test")?;
    let prefix_path = prefix.path().to_path_buf();
    let expected_metrics =
        std::path::PathBuf::from(format!("{}{METRICS_SUFFIX}", prefix_path.display()));

    let mut collector =
        AlignmentCollector::new(bam.path(), &prefix_path, None, &AlignmentOptions::default());
    drive_collector(&mut collector, bam.path())?;

    assert!(expected_metrics.exists(), "metrics file not created at expected path");
    Ok(())
}

// ─── PAIR bad_cycles override ─────────────────────────────────────────────────

#[test]
fn test_pair_bad_cycles_is_sum_of_first_and_second() -> Result<()> {
    // 5 pairs of length-10 reads. Every R1 has an N at cycle 1 (forward, offset 0);
    // every R2 has an N at offset 0, which is cycle 10 on the reverse strand. Both are
    // 100% no-call at their cycle, so read1 and read2 each get 1 bad cycle and the PAIR
    // row is their sum.
    let mut sam = SamBuilder::new();
    for i in 0..5 {
        sam.add(
            pair(sam.next_id())
                .at("chr1", 100 + i * 200, 300 + i * 200)
                .len(10)
                .r1(|r| r.sub(0, b'N'))
                .r2(|r| r.sub(0, b'N')),
        );
    }
    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let first = row(&metrics, "read1");
    let second = row(&metrics, "read2");
    let pair = row(&metrics, "pair");

    // PAIR bad_cycles must equal FIRST + SECOND.
    assert_eq!(pair.bad_cycles, first.bad_cycles + second.bad_cycles);
    assert_eq!(first.bad_cycles, 1, "read1 bad_cycles");
    assert_eq!(second.bad_cycles, 1, "read2 bad_cycles");
    assert_eq!(pair.bad_cycles, 2, "pair bad_cycles");
    Ok(())
}

// ─── Improper pairs ───────────────────────────────────────────────────────────

#[test]
fn test_improper_pairs_counted() -> Result<()> {
    let mut sam = SamBuilder::new();

    // 3 proper pairs.
    for i in 0..3 {
        sam.add(pair(sam.next_id()).at("chr1", 100 + i * 200, 300 + i * 200));
    }

    // 2 improper pairs: segmented but not properly paired.
    for i in 0..2 {
        sam.add(pair(sam.next_id()).at("chr1", 10_000 + i * 500, 10_200 + i * 500).improper());
    }

    let bam = sam.to_temp_bam()?;
    let metrics = run_alignment(bam.path())?;

    let pair = row(&metrics, "pair");
    // 4 improper reads (2 pairs × 2) / 10 total aligned = 0.4.
    assert_eq!(pair.reads_improperly_paired, 4, "reads_improperly_paired");
    assert_close!(pair.frac_reads_improperly_paired, 0.4, 1e-5);
    Ok(())
}
