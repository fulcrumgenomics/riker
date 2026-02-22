mod helpers;

use helpers::{FastaBuilder, SamBuilder, SortOrder, coord_builder, read_metrics_tsv};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::{
    Flags, MappingQuality,
    cigar::{Op, op::Kind},
};
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use riker_lib::commands::hybcap::{
    HybCap, HybCapMetric, HybCapOptions, METRICS_SUFFIX, PER_BASE_SUFFIX, PER_TARGET_SUFFIX,
    PerTargetCoverage,
};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Write a BED file to disk. Entries are (contig, start_0based, end_exclusive).
fn write_bed(path: &Path, entries: &[(&str, u64, u64)]) {
    use std::fmt::Write;
    let content = entries.iter().fold(String::new(), |mut s, (chr, start, end)| {
        writeln!(s, "{chr}\t{start}\t{end}").unwrap();
        s
    });
    std::fs::write(path, content).unwrap();
}

/// Write a BED file with names. Entries are (contig, start_0based, end_exclusive, name).
fn write_named_bed(path: &Path, entries: &[(&str, u64, u64, &str)]) {
    use std::fmt::Write;
    let content = entries.iter().fold(String::new(), |mut s, (chr, start, end, name)| {
        writeln!(s, "{chr}\t{start}\t{end}\t{name}").unwrap();
        s
    });
    std::fs::write(path, content).unwrap();
}

/// Build a hybcap command with common defaults.
#[allow(clippy::too_many_arguments)]
fn make_cmd(
    bam: &Path,
    baits: &Path,
    targets: &Path,
    prefix: &Path,
    reference: Option<&Path>,
    options: HybCapOptions,
) -> HybCap {
    let options =
        HybCapOptions { baits: baits.to_path_buf(), targets: targets.to_path_buf(), ..options };
    HybCap {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.to_path_buf() },
        reference: OptionalReferenceOptions { reference: reference.map(Path::to_path_buf) },
        options,
    }
}

/// Create default options with specific MAPQ and BQ thresholds.
fn opts(min_mapq: u8, min_bq: u8) -> HybCapOptions {
    HybCapOptions {
        min_mapping_quality: min_mapq,
        min_base_quality: min_bq,
        ..HybCapOptions::default()
    }
}

/// Run hybcap and return the single metrics row.
fn run_and_read(cmd: &HybCap) -> HybCapMetric {
    cmd.execute().unwrap();
    let metrics_path = PathBuf::from(format!("{}{METRICS_SUFFIX}", cmd.output.output.display()));
    let rows: Vec<HybCapMetric> = read_metrics_tsv(&metrics_path).unwrap();
    assert_eq!(rows.len(), 1);
    rows.into_iter().next().unwrap()
}

/// Build a paired-end read pair with custom CIGAR and quality scores.
///
/// Both reads get the same CIGAR and quality. `pos1`/`pos2` are 1-based.
#[allow(clippy::too_many_arguments)]
fn add_custom_pair(
    builder: &mut SamBuilder,
    name: &str,
    ref_id: usize,
    pos1: usize,
    pos2: usize,
    cigar_ops: &[Op],
    qual: &[u8],
    mapq: u8,
    is_duplicate: bool,
    fails_vendor_quality: bool,
) {
    let cigar: Cigar = cigar_ops.iter().copied().collect();
    let read_len: usize = cigar_ops
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Insertion
                    | Kind::SoftClip
            )
        })
        .map(|op| op.len())
        .sum();
    let seq: Sequence = vec![b'A'; read_len].into();
    let qual_scores = QualityScores::from(qual.to_vec());
    let mq = MappingQuality::new(mapq).unwrap();

    let mut flags1 = Flags::SEGMENTED
        | Flags::PROPERLY_SEGMENTED
        | Flags::MATE_REVERSE_COMPLEMENTED
        | Flags::FIRST_SEGMENT;
    let mut flags2 = Flags::SEGMENTED
        | Flags::PROPERLY_SEGMENTED
        | Flags::REVERSE_COMPLEMENTED
        | Flags::LAST_SEGMENT;
    if is_duplicate {
        flags1 |= Flags::DUPLICATE;
        flags2 |= Flags::DUPLICATE;
    }
    if fails_vendor_quality {
        flags1 |= Flags::QC_FAIL;
        flags2 |= Flags::QC_FAIL;
    }

    let pos1_p = Position::new(pos1).unwrap();
    let pos2_p = Position::new(pos2).unwrap();
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    let tlen = (pos2 as i32 - pos1 as i32) + read_len as i32;

    let r1 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags1)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(pos1_p)
        .set_mapping_quality(mq)
        .set_cigar(cigar.clone())
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(pos2_p)
        .set_template_length(tlen)
        .set_sequence(seq.clone())
        .set_quality_scores(qual_scores.clone())
        .build();

    let r2 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags2)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(pos2_p)
        .set_mapping_quality(mq)
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(pos1_p)
        .set_template_length(-tlen)
        .set_sequence(seq)
        .set_quality_scores(qual_scores)
        .build();

    builder.add_record(r1);
    builder.add_record(r2);
}

/// Add a single mapped record with custom CIGAR, quality, and flags.
#[allow(clippy::too_many_arguments)]
fn add_single_record(
    builder: &mut SamBuilder,
    name: &str,
    ref_id: usize,
    pos: usize,
    cigar_ops: &[Op],
    qual: &[u8],
    mapq: u8,
    flags: Flags,
) {
    let cigar: Cigar = cigar_ops.iter().copied().collect();
    let read_len: usize = cigar_ops
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Insertion
                    | Kind::SoftClip
            )
        })
        .map(|op| op.len())
        .sum();
    let seq: Sequence = vec![b'A'; read_len].into();
    let qual_scores = QualityScores::from(qual.to_vec());
    let mq = MappingQuality::new(mapq).unwrap();

    let record = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(pos).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar)
        .set_sequence(seq)
        .set_quality_scores(qual_scores)
        .build();

    builder.add_record(record);
}

// ─── Tests ported from Picard ─────────────────────────────────────────────────

/// Port of Picard test #1: Two 100bp reads with different base qualities.
/// Read 1 has good quality (Q=30), Read 2 has poor quality (Q=2).
/// With min_bq=10, half the bases should be excluded for base quality.
#[test]
fn test_low_base_quality_exclusion() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Target covers 200bp: chr1:0-200
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let good_qual = vec![30u8; 100];
    let bad_qual = vec![2u8; 100];
    // Two unpaired reads at pos 1, 100bp each
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &good_qual,
        20,
        Flags::empty(),
    );
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &bad_qual,
        20,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(1, 10));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 2);
    assert_eq!(m.deduped_bases_aligned, 200);
    assert_float_eq!(m.frac_exc_baseq, 0.5, 0.001);
    assert_float_eq!(m.frac_exc_overlap, 0.0, 0.001);
    assert_float_eq!(m.frac_exc_off_target, 0.0, 0.001);
    assert_float_eq!(m.frac_target_bases_1x, 0.5, 0.001);
    assert_eq!(m.max_target_coverage, 1);
    assert_eq!(m.min_target_coverage, 0);
    assert_eq!(m.total_bases, 200);
}

/// Port of Picard test #3: Two 101bp reads with different MAPQ (100 and 1).
/// With min_mapq=2, the low MAPQ read should be excluded.
#[test]
fn test_low_mapq_filtering() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 101];
    // Read 1: MAPQ=100
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 101)],
        &qual,
        100,
        Flags::empty(),
    );
    // Read 2: MAPQ=1
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 101)],
        &qual,
        1,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(2, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 2);
    assert_eq!(m.bases_aligned, 202);
    assert_float_eq!(m.frac_exc_baseq, 0.0, 0.001);
    // One read filtered by MAPQ: 101 / 202 bases
    assert_float_eq!(m.frac_exc_mapq, 101.0 / 202.0, 0.001);
    // Only 101 bases survive, covering positions 1-101 out of 200bp target
    assert_float_eq!(m.frac_target_bases_1x, 101.0 / 200.0, 0.001);
}

/// Port of Picard test #5: Two overlapping paired-end reads with clipping.
/// Both 101bp reads start at pos 1 (insert size = 360).
/// With clipping enabled, overlap bases should be excluded.
#[test]
fn test_overlapping_reads_with_clipping() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 101];
    // Paired reads: both at pos 1, 101bp, perfectly overlapping
    add_custom_pair(
        &mut bld,
        "pair1",
        0,
        1,
        1,
        &[Op::new(Kind::Match, 101)],
        &qual,
        60,
        false,
        false,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 2);
    assert_eq!(m.bases_aligned, 202);
    // The second-of-pair (left-most tie) gets fully clipped
    assert_float_eq!(m.frac_exc_overlap, 101.0 / 202.0, 0.01);
    assert_float_eq!(m.frac_target_bases_1x, 101.0 / 200.0, 0.01);
}

/// Port of Picard test #6: Overlapping reads WITHOUT clipping.
/// Same setup but CLIP_OVERLAPPING_READS=false.
#[test]
fn test_overlapping_reads_without_clipping() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 101];
    add_custom_pair(
        &mut bld,
        "pair1",
        0,
        1,
        1,
        &[Op::new(Kind::Match, 101)],
        &qual,
        60,
        false,
        false,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_float_eq!(m.frac_exc_overlap, 0.0, 0.001);
    // Both reads cover pos 1-101, so 101 bases covered at 2x → max coverage = 2
    assert_eq!(m.max_target_coverage, 2);
    // On-target bases = both reads' bases = 202 (both count)
    assert_eq!(m.on_target_bases, 202);
}

/// Port of Picard test #9: Short read covering only one of two targets.
/// Single 10bp read at pos 1-10. Two intervals: [0-10] and [90-100].
/// Only first interval covered → frac_target_bases_1x = 0.5.
#[test]
fn test_short_read_partial_target() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Two 10bp intervals
    write_bed(&bait_path, &[("chr1", 0, 10), ("chr1", 90, 100)]);
    write_bed(&target_path, &[("chr1", 0, 10), ("chr1", 90, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 10];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 10)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 1);
    assert_eq!(m.deduped_bases_aligned, 10);
    assert_eq!(m.target_territory, 20);
    assert_float_eq!(m.frac_target_bases_1x, 0.5, 0.001);
    assert_eq!(m.max_target_coverage, 1);
    assert_eq!(m.min_target_coverage, 0);
}

/// Port of Picard F80 test: Uniform coverage should give fold_80 = 1.0.
/// 100 reads each covering a 200bp target fully → uniform depth.
#[test]
fn test_fold_80_uniform_coverage() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 999, 1199)]);
    write_bed(&target_path, &[("chr1", 999, 1199)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 200];
    for i in 0..100 {
        add_single_record(
            &mut bld,
            &format!("read{i}"),
            0,
            1000,
            &[Op::new(Kind::Match, 200)],
            &qual,
            60,
            Flags::empty(),
        );
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_float_eq!(m.mean_target_coverage, 100.0, 0.01);
    assert_float_eq!(m.median_target_coverage, 100.0, 0.01);
    assert_float_eq!(m.fold_80_base_penalty, 1.0, 0.01);
    assert_float_eq!(m.frac_target_bases_100x, 1.0, 0.001);
    assert_float_eq!(m.frac_target_bases_1000x, 0.0, 0.001);
}

/// Port of Picard indel test: Deletion without INCLUDE_INDELS reduces coverage.
/// 100 reads with 90M20D10M covering a 120bp target.
/// Without INCLUDE_INDELS, deletion region has 0 depth → mean = (90*100 + 10*100) / 120.
#[test]
fn test_deletion_without_include_indels() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 120)]);
    write_bed(&target_path, &[("chr1", 0, 120)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100]; // 90M + 10M = 100 read bases
    let cigar = [Op::new(Kind::Match, 90), Op::new(Kind::Deletion, 20), Op::new(Kind::Match, 10)];
    for i in 0..100 {
        add_single_record(&mut bld, &format!("read{i}"), 0, 1, &cigar, &qual, 60, Flags::empty());
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.include_indels = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    // 90 bases at 100x + 20 bases at 0x + 10 bases at 100x = 10000 total / 120 = 83.33
    assert_float_eq!(m.mean_target_coverage, 10_000.0 / 120.0, 0.5);
    assert_eq!(m.min_target_coverage, 0);
}

/// Port of Picard indel test: Deletion WITH INCLUDE_INDELS fills in coverage.
#[test]
fn test_deletion_with_include_indels() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 120)]);
    write_bed(&target_path, &[("chr1", 0, 120)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    let cigar = [Op::new(Kind::Match, 90), Op::new(Kind::Deletion, 20), Op::new(Kind::Match, 10)];
    for i in 0..100 {
        add_single_record(&mut bld, &format!("read{i}"), 0, 1, &cigar, &qual, 60, Flags::empty());
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.include_indels = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    // With indels: 120 bases all at 100x = 100.0 mean
    assert_float_eq!(m.mean_target_coverage, 100.0, 0.5);
    assert_eq!(m.min_target_coverage, 100);
}

// ─── Gap-filling tests ─────────────────────────────────────────────────────

/// Duplicate reads should be excluded from on-target coverage.
#[test]
fn test_duplicate_exclusion() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // One non-duplicate
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    // One duplicate
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::DUPLICATE,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 2);
    assert_eq!(m.deduped_reads, 1);
    assert_float_eq!(m.frac_exc_dupe, 0.5, 0.01);
    assert_float_eq!(m.mean_target_coverage, 1.0, 0.01);
}

/// include_duplicates=true should count duplicates toward coverage.
#[test]
fn test_include_duplicates() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::DUPLICATE,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.include_duplicates = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_float_eq!(m.frac_exc_dupe, 0.0, 0.001);
    assert_float_eq!(m.mean_target_coverage, 2.0, 0.01);
}

/// QC fail reads should be excluded from all computations.
#[test]
fn test_qc_fail_excluded() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::QC_FAIL,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 1); // QC-fail read is excluded entirely
    assert_float_eq!(m.mean_target_coverage, 1.0, 0.01);
}

/// Secondary alignments should be completely skipped.
#[test]
fn test_secondary_skipped() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Primary read
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    // Secondary alignment (should be skipped entirely)
    add_single_record(
        &mut bld,
        "read1_sec",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::SECONDARY,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 1); // secondary not counted
    assert_float_eq!(m.mean_target_coverage, 1.0, 0.01);
}

/// Off-target bases should be excluded and counted as exc_off_target.
#[test]
fn test_off_target_exclusion() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Bait and target only cover positions 0-50
    write_bed(&bait_path, &[("chr1", 0, 50)]);
    write_bed(&target_path, &[("chr1", 0, 50)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Read covers 1-100, but only 50bp on target
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_float_eq!(m.frac_target_bases_1x, 1.0, 0.001);
    assert_eq!(m.on_target_bases, 50);
    // 50 off-target bases out of 100 aligned
    assert_float_eq!(m.frac_exc_off_target, 0.5, 0.01);
}

/// Bait territory and target territory should be computed correctly after merging.
#[test]
fn test_territory_computation() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Overlapping baits: [0-60) and [40-100) → merged [0-100) = 100bp
    write_bed(&bait_path, &[("chr1", 0, 60), ("chr1", 40, 100)]);
    // Non-overlapping targets: [0-30) and [70-100) = 60bp
    write_bed(&target_path, &[("chr1", 0, 30), ("chr1", 70, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.bait_territory, 100);
    assert_eq!(m.target_territory, 60);
    assert_float_eq!(m.bait_design_efficiency, 60.0 / 100.0, 0.001);
}

/// Genome size should be the sum of all contig lengths from the BAM header.
#[test]
fn test_genome_size() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld =
        SamBuilder::with_contigs(&[("chr1".to_string(), 1000), ("chr2".to_string(), 2000)])
            .sort_order(SortOrder::Coordinate);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.genome_size, 3000);
}

/// Per-target coverage output should have one row per merged target.
#[test]
fn test_per_target_coverage_output() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Non-abutting intervals so they don't merge
    write_named_bed(&bait_path, &[("chr1", 0, 50, "bait1"), ("chr1", 60, 110, "bait2")]);
    write_named_bed(&target_path, &[("chr1", 0, 50, "target1"), ("chr1", 60, 110, "target2")]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 50];
    // Only cover first target
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.per_target_coverage = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    cmd.execute().unwrap();

    let per_target_path = PathBuf::from(format!("{}{PER_TARGET_SUFFIX}", prefix.display()));
    let rows: Vec<PerTargetCoverage> = read_metrics_tsv(&per_target_path).unwrap();
    assert_eq!(rows.len(), 2);

    // First target: covered
    assert_eq!(rows[0].name, "target1");
    assert_float_eq!(rows[0].mean_coverage, 1.0, 0.01);
    assert_eq!(rows[0].min_coverage, 1);

    // Second target: not covered
    assert_eq!(rows[1].name, "target2");
    assert_float_eq!(rows[1].mean_coverage, 0.0, 0.01);
    assert_eq!(rows[1].max_coverage, 0);
    assert_float_eq!(rows[1].frac_0x, 1.0, 0.001);
}

/// Per-base coverage output should have one row per target base position.
#[test]
fn test_per_base_coverage_output() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 5)]);
    write_bed(&target_path, &[("chr1", 0, 5)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 3];
    // Read covers positions 1-3 (0-based: 0-2) → first 3 of 5 target bases
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 3)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.per_base_coverage = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    cmd.execute().unwrap();

    let per_base_path = PathBuf::from(format!("{}{PER_BASE_SUFFIX}", prefix.display()));
    let rows: Vec<riker_lib::commands::hybcap::PerBaseCoverage> =
        read_metrics_tsv(&per_base_path).unwrap();
    assert_eq!(rows.len(), 5);

    // First 3 positions covered at 1x
    for row in &rows[0..3] {
        assert_eq!(row.coverage, 1);
    }
    // Last 2 positions at 0x
    for row in &rows[3..5] {
        assert_eq!(row.coverage, 0);
    }
    // Positions should be 1-based
    assert_eq!(rows[0].pos, 1);
    assert_eq!(rows[4].pos, 5);
}

/// Bait classification: on-bait vs near-bait vs off-bait.
#[test]
fn test_bait_classification() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Bait covers positions 0-50
    write_bed(&bait_path, &[("chr1", 0, 50)]);
    write_bed(&target_path, &[("chr1", 0, 50)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Read covers 1-100 (0-based: 0-99). First 50bp on bait, next 50bp near bait
    // (within 250bp near_distance default)
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.on_bait_bases, 50);
    assert_eq!(m.near_bait_bases, 50);
    assert_eq!(m.off_bait_bases, 0);
}

/// Off-bait bases: read entirely outside expanded bait region.
#[test]
fn test_off_bait_classification() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 50)]);
    write_bed(&target_path, &[("chr1", 0, 50)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 50];
    // Read is >250bp away from bait → off-bait
    add_single_record(
        &mut bld,
        "read1",
        0,
        500,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.on_bait_bases, 0);
    assert_eq!(m.near_bait_bases, 0);
    assert_eq!(m.off_bait_bases, 50);
}

/// Zero coverage targets: one target with reads, one without.
#[test]
fn test_zero_coverage_targets() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 50), ("chr1", 500, 550)]);
    write_bed(&target_path, &[("chr1", 0, 50), ("chr1", 500, 550)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 50];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    // 1 out of 2 targets has zero coverage
    assert_float_eq!(m.frac_uncovered_targets, 0.5, 0.001);
}

/// GC/AT dropout with reference: verify non-zero values when coverage is skewed.
#[test]
fn test_gc_dropout_with_reference() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Two non-abutting targets: one AT-rich, one GC-rich (gap so they don't merge)
    write_bed(&bait_path, &[("chr1", 0, 100), ("chr1", 200, 300)]);
    write_bed(&target_path, &[("chr1", 0, 100), ("chr1", 200, 300)]);

    // Reference: AT-rich at 0-100, padding at 100-200, GC-rich at 200-300, rest N
    let at_bases = [b'A'; 100];
    let padding1 = [b'N'; 100];
    let gc_bases = [b'G'; 100];
    let padding2 = [b'N'; 9700];
    let ref_seq: Vec<u8> = at_bases
        .iter()
        .chain(padding1.iter())
        .chain(gc_bases.iter())
        .chain(padding2.iter())
        .copied()
        .collect();
    let refa = FastaBuilder::new().add_contig("chr1", &ref_seq).to_temp_fasta().unwrap();

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Only cover the AT-rich target (GC-rich has no coverage → GC dropout)
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd =
        make_cmd(bam.path(), &bait_path, &target_path, &prefix, Some(refa.path()), opts(0, 0));
    let m = run_and_read(&cmd);

    // GC dropout should be positive (GC-rich target has less coverage than expected)
    assert!(m.gc_dropout > 0.0, "gc_dropout should be > 0, got {}", m.gc_dropout);
}

/// Panel name should default to bait filename stem.
#[test]
fn test_panel_name_default() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("my_panel.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.panel_name, "my_panel");
}

/// Custom panel name should override the default.
#[test]
fn test_panel_name_custom() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.panel_name = Some("TwistExome".to_string());
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.panel_name, "TwistExome");
}

/// Fold enrichment: on-bait fraction vs genome fraction.
#[test]
fn test_fold_enrichment() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Bait covers 100bp out of 10000bp genome
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Read is entirely on bait
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    // 100% on bait, bait is 100/10000 = 1% of genome → enrichment = 1.0 / 0.01 = 100
    assert_float_eq!(m.fold_enrichment, 100.0, 1.0);
}

/// Mixed base quality: half good, half bad within same read.
#[test]
fn test_mixed_base_quality_in_read() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    // First 50 bases Q=30, last 50 bases Q=2
    let mut qual = vec![30u8; 50];
    qual.extend(vec![2u8; 50]);
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 10));
    let m = run_and_read(&cmd);

    assert_eq!(m.on_target_bases, 50);
    assert_float_eq!(m.frac_exc_baseq, 0.5, 0.01);
}

/// Empty BAM (no reads) should produce all-zero metrics.
#[test]
fn test_empty_bam() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let bld = coord_builder(&[("chr1", 10_000)]);
    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 0);
    assert_float_eq!(m.mean_target_coverage, 0.0, 0.001);
    assert_eq!(m.max_target_coverage, 0);
    assert_eq!(m.min_target_coverage, 0);
}

/// Multiple reads building up depth: verify coverage fractions.
#[test]
fn test_coverage_fractions_at_thresholds() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // 5 reads covering all 100bp → 5x coverage
    for i in 0..5 {
        add_single_record(
            &mut bld,
            &format!("read{i}"),
            0,
            1,
            &[Op::new(Kind::Match, 100)],
            &qual,
            60,
            Flags::empty(),
        );
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_float_eq!(m.mean_target_coverage, 5.0, 0.01);
    assert_float_eq!(m.frac_target_bases_1x, 1.0, 0.001);
    assert_float_eq!(m.frac_target_bases_10x, 0.0, 0.001);
}

/// Read with insertion: insertion bases should not count without INCLUDE_INDELS.
#[test]
fn test_insertion_without_include_indels() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    // 50M + 20I + 50M = 120 read bases, 100 ref-consuming
    let qual = vec![30u8; 120];
    let cigar = [Op::new(Kind::Match, 50), Op::new(Kind::Insertion, 20), Op::new(Kind::Match, 50)];
    add_single_record(&mut bld, "read1", 0, 1, &cigar, &qual, 60, Flags::empty());

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.include_indels = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    // Only M bases counted on-target (100 out of 100bp target)
    assert_eq!(m.on_target_bases, 100);
}

/// Read with insertion and INCLUDE_INDELS: inserted bases should also count.
#[test]
fn test_insertion_with_include_indels() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 120];
    let cigar = [Op::new(Kind::Match, 50), Op::new(Kind::Insertion, 20), Op::new(Kind::Match, 50)];
    add_single_record(&mut bld, "read1", 0, 1, &cigar, &qual, 60, Flags::empty());

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.include_indels = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    // M bases + I bases: 100 + 20 = 120
    assert_eq!(m.on_target_bases, 120);
}

/// Paired reads classified as selected_pairs when overlapping expanded baits.
#[test]
fn test_selected_pairs() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 50];
    // Pair on-bait
    add_custom_pair(
        &mut bld,
        "pair1",
        0,
        1,
        51,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        false,
        false,
    );
    // Pair off-bait (far from bait, > near_distance)
    add_custom_pair(
        &mut bld,
        "pair2",
        0,
        500,
        600,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        false,
        false,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    // Only the on-bait pair should be counted as selected
    assert_eq!(m.selected_pairs, 1);
    assert_eq!(m.selected_unique_pairs, 1);
}

/// Soft-clipped bases should not contribute to aligned bases or coverage.
#[test]
fn test_soft_clip_excluded() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    // 10S80M10S = 100 read bases, 80 aligned
    let qual = vec![30u8; 100];
    let cigar =
        [Op::new(Kind::SoftClip, 10), Op::new(Kind::Match, 80), Op::new(Kind::SoftClip, 10)];
    add_single_record(&mut bld, "read1", 0, 1, &cigar, &qual, 60, Flags::empty());

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.deduped_bases_aligned, 80);
    assert_eq!(m.on_target_bases, 80);
}

/// Overlapping read pair where only part of the read overlaps: verify partial clipping.
#[test]
fn test_partial_overlap_clipping() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];

    // Read1 at pos 1 (100M), Read2 at pos 51 (100M)
    // Overlap: positions 51-100 = 50bp
    let cigar = [Op::new(Kind::Match, 100)];
    let cigar_c: Cigar = cigar.iter().copied().collect();
    let seq: Sequence = vec![b'A'; 100].into();
    let mq = MappingQuality::new(60).unwrap();

    let r1 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar_c.clone())
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(51).unwrap())
        .set_template_length(150)
        .set_sequence(seq.clone())
        .set_quality_scores(QualityScores::from(qual.clone()))
        .build();

    let r2 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(51).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar_c)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(1).unwrap())
        .set_template_length(-150)
        .set_sequence(seq)
        .set_quality_scores(QualityScores::from(qual))
        .build();

    bld.add_record(r1);
    bld.add_record(r2);

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    // Read1 starts at 1, mate at 51 → Read1 is left-most, overlap = 50bp clipped from Read1
    // Total aligned = 200, overlap = 50
    assert_float_eq!(m.frac_exc_overlap, 50.0 / 200.0, 0.01);
    // Covered positions: Read1 covers 1-50 (after clip), Read2 covers 51-150
    assert_eq!(m.on_target_bases, 150);
}

/// Library size estimation requires duplicate pairs.
#[test]
fn test_hs_library_size_no_duplicates() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // 5 unique pairs, no duplicates → selected_pairs == selected_unique_pairs
    for i in 0..5 {
        add_custom_pair(
            &mut bld,
            &format!("pair{i}"),
            0,
            1 + i * 10,
            101 + i * 10,
            &[Op::new(Kind::Match, 100)],
            &qual,
            60,
            false,
            false,
        );
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    // No duplicates → library size = None (all pairs unique)
    assert!(m.hs_library_size.is_none());
}

/// Verify frac_deduped_reads fraction.
#[test]
fn test_fraction_computations() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // 3 reads: 1 non-dup, 1 dup, 1 QC-fail
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    add_single_record(
        &mut bld,
        "read2",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::DUPLICATE,
    );
    add_single_record(
        &mut bld,
        "read3",
        0,
        1,
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::QC_FAIL,
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 2); // QC-fail read excluded; 1 non-dup + 1 dup
    assert_eq!(m.deduped_reads, 1); // 1 non-dup read
    assert_float_eq!(m.frac_deduped_reads, 1.0 / 2.0, 0.01);
}

// ─── Additional integration tests ────────────────────────────────────────────

/// A read spanning two separate (non-abutting) target intervals should
/// contribute coverage to both targets.
#[test]
fn test_multiple_overlapping_targets() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // Two non-abutting targets: bases 0-10 and 20-30 on chr1
    write_bed(&bait_path, &[("chr1", 0, 10), ("chr1", 20, 30)]);
    write_bed(&target_path, &[("chr1", 0, 10), ("chr1", 20, 30)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 50];
    // One 50bp read at pos 1, covering bases 0-49 (overlaps both targets)
    add_single_record(
        &mut bld,
        "read1",
        0,
        1,
        &[Op::new(Kind::Match, 50)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 1);
    // All 20 target bases (10 + 10) should have coverage from the read
    assert_eq!(m.on_target_bases, 20);
    // Both targets fully covered at 1x
    assert_float_eq!(m.frac_target_bases_1x, 1.0, 0.001);
    assert_float_eq!(m.mean_target_coverage, 1.0, 0.001);
    assert_eq!(m.min_target_coverage, 1);
    assert_eq!(m.max_target_coverage, 1);
}

/// When every read has MAPQ below the threshold, all bases should be excluded
/// and coverage should be zero everywhere.
#[test]
fn test_all_reads_filtered_by_mapq() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // Three reads, all with MAPQ=1 (below the min_mapq=20 threshold)
    for i in 0..3 {
        add_single_record(
            &mut bld,
            &format!("read{i}"),
            0,
            1,
            &[Op::new(Kind::Match, 100)],
            &qual,
            1,
            Flags::empty(),
        );
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(20, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 3);
    assert_eq!(m.bases_aligned, 300);
    // All aligned bases excluded by MAPQ
    assert_float_eq!(m.frac_exc_mapq, 1.0, 0.001);
    // Zero coverage on all targets
    assert_float_eq!(m.mean_target_coverage, 0.0, 0.001);
    assert_eq!(m.max_target_coverage, 0);
    assert_eq!(m.min_target_coverage, 0);
    assert_float_eq!(m.frac_target_bases_1x, 0.0, 0.001);
    assert_eq!(m.on_target_bases, 0);
}

/// High coverage on a small target to verify no overflow issues with u64
/// coverage accumulators. 1000 reads on a 10bp target.
#[test]
fn test_high_coverage_u64_safety() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // 10bp target
    write_bed(&bait_path, &[("chr1", 0, 10)]);
    write_bed(&target_path, &[("chr1", 0, 10)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 10];
    // 1000 reads, each covering the full 10bp target
    for i in 0..1000 {
        add_single_record(
            &mut bld,
            &format!("read{i}"),
            0,
            1,
            &[Op::new(Kind::Match, 10)],
            &qual,
            60,
            Flags::empty(),
        );
    }

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, opts(0, 0));
    let m = run_and_read(&cmd);

    assert_eq!(m.total_reads, 1000);
    assert_eq!(m.on_target_bases, 10_000);
    assert_float_eq!(m.mean_target_coverage, 1000.0, 0.01);
    assert_eq!(m.max_target_coverage, 1000);
    assert_eq!(m.min_target_coverage, 1000);
    assert_float_eq!(m.frac_target_bases_1x, 1.0, 0.001);
    assert_float_eq!(m.frac_target_bases_100x, 1.0, 0.001);
    assert_float_eq!(m.frac_target_bases_500x, 1.0, 0.001);
    assert_float_eq!(m.frac_target_bases_1000x, 1.0, 0.001);
}

// ─── Overlap clipping with soft clips regression tests ───────────────────────

/// Regression test: overlap clipping must account for trailing soft clips.
///
/// Read1 at pos 1 CIGAR 90M10S, Read2 at pos 51 CIGAR 90M10S.
/// Read1 covers ref 0-89, Read2 covers ref 50-139. Overlap = 40bp on ref 50-89.
/// Read1 (left-most, first-of-pair → clipped on tie break irrelevant since pos1 < pos2).
/// With trailing soft clip fix: clip_start = (100 - 10) - 40 = 50 → Read1 contributes 50 bases.
/// on_target_bases = 50 + 90 = 140, frac_exc_overlap = 40 / 180.
#[test]
fn test_overlap_clipping_with_trailing_soft_clips() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100]; // 90M + 10S = 100 read bases
    let cigar_ops = [Op::new(Kind::Match, 90), Op::new(Kind::SoftClip, 10)];
    add_custom_pair(&mut bld, "pair1", 0, 1, 51, &cigar_ops, &qual, 60, false, false);

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.on_target_bases, 140);
    assert_float_eq!(m.frac_exc_overlap, 40.0 / 180.0, 0.001);
}

/// Overlap clipping with both leading and trailing soft clips.
///
/// Read1 at pos 1 CIGAR 5S85M10S, Read2 at pos 46 CIGAR 5S85M10S.
/// Read1 covers ref 0-84, Read2 covers ref 45-129. Overlap = 40bp on ref 45-84.
/// Read1 gets 40bp clipped: clip_start = (100 - 10) - 40 = 50 → 50 read bases kept.
/// Of those 50, the first 5 are leading soft clip → 45 aligned bases contributed (ref 0-44).
/// on_target_bases = 45 + 85 = 130, frac_exc_overlap = 40 / 170.
#[test]
fn test_overlap_clipping_with_leading_and_trailing_soft_clips() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100]; // 5S + 85M + 10S = 100 read bases
    let mq = MappingQuality::new(60).unwrap();
    let cigar_ops =
        [Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 85), Op::new(Kind::SoftClip, 10)];
    let cigar: Cigar = cigar_ops.iter().copied().collect();
    let seq: Sequence = vec![b'A'; 100].into();

    let r1 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar.clone())
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(46).unwrap())
        .set_template_length(130)
        .set_sequence(seq.clone())
        .set_quality_scores(QualityScores::from(qual.clone()))
        .build();

    let r2 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(46).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(1).unwrap())
        .set_template_length(-130)
        .set_sequence(seq)
        .set_quality_scores(QualityScores::from(qual))
        .build();

    bld.add_record(r1);
    bld.add_record(r2);

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.on_target_bases, 130);
    assert_float_eq!(m.frac_exc_overlap, 40.0 / 170.0, 0.001);
}

/// Edge case: overlap is large enough that the entire aligned portion of Read1 is clipped.
///
/// Read1 at pos 1 CIGAR 50M50S, Read2 at pos 1 CIGAR 50M50S.
/// Both cover ref 0-49 → overlap = 50bp.
/// Read1 (first-of-pair at same pos) not clipped; Read2 (second-of-pair) gets clipped.
/// clip_start = (100 - 50) - 50 = 0 → all aligned bases clipped from Read2.
/// on_target_bases = 50 + 0 = 50, frac_exc_overlap = 50 / 100.
#[test]
fn test_overlap_clipping_trailing_soft_clip_fully_clips_read() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 100)]);
    write_bed(&target_path, &[("chr1", 0, 100)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100]; // 50M + 50S = 100 read bases
    let cigar_ops = [Op::new(Kind::Match, 50), Op::new(Kind::SoftClip, 50)];
    add_custom_pair(&mut bld, "pair1", 0, 1, 1, &cigar_ops, &qual, 60, false, false);

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.on_target_bases, 50);
    assert_float_eq!(m.frac_exc_overlap, 50.0 / 100.0, 0.001);
}

/// Verify that per-target output uses 1-based coordinates (not 0-based from BED).
#[test]
fn test_per_target_coordinates_are_one_based() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    // BED is 0-based half-open
    write_named_bed(&bait_path, &[("chr1", 100, 200, "bait1"), ("chr1", 500, 600, "bait2")]);
    write_named_bed(&target_path, &[("chr1", 100, 200, "target1"), ("chr1", 500, 600, "target2")]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let qual = vec![30u8; 100];
    // One read covering each target
    add_single_record(
        &mut bld,
        "read1",
        0,
        101, // 1-based pos covering target1
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );
    add_single_record(
        &mut bld,
        "read2",
        0,
        501, // 1-based pos covering target2
        &[Op::new(Kind::Match, 100)],
        &qual,
        60,
        Flags::empty(),
    );

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.per_target_coverage = true;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    cmd.execute().unwrap();

    let per_target_path = PathBuf::from(format!("{}{PER_TARGET_SUFFIX}", prefix.display()));
    let rows: Vec<PerTargetCoverage> = read_metrics_tsv(&per_target_path).unwrap();
    assert_eq!(rows.len(), 2);

    // BED 100-200 → 1-based 101-200
    assert_eq!(rows[0].start, 101);
    assert_eq!(rows[0].end, 200);
    assert_eq!(rows[0].length, 100);

    // BED 500-600 → 1-based 501-600
    assert_eq!(rows[1].start, 501);
    assert_eq!(rows[1].end, 600);
    assert_eq!(rows[1].length, 100);
}

/// Trailing soft clips followed by hard clips should still be detected correctly.
///
/// Read1 at pos 1 CIGAR 80M10S5H, Read2 at pos 51 CIGAR 80M.
/// Read1 covers ref 0-79, Read2 covers ref 50-129. Overlap = 30bp on ref 50-79.
/// Read1 gets clipped: clip_start = (90 - 10) - 30 = 50 → keeps 50 aligned bases (ref 0-49).
/// on_target_bases = 50 + 80 = 130, frac_exc_overlap = 30 / 160.
#[test]
fn test_overlap_clipping_trailing_soft_clip_with_hard_clip() {
    let dir = TempDir::new().unwrap();
    let bait_path = dir.path().join("baits.bed");
    let target_path = dir.path().join("targets.bed");
    write_bed(&bait_path, &[("chr1", 0, 200)]);
    write_bed(&target_path, &[("chr1", 0, 200)]);

    let mut bld = coord_builder(&[("chr1", 10_000)]);
    let mq = MappingQuality::new(60).unwrap();

    // Read1: 80M10S5H → 90 sequence bases (H doesn't count), 80 aligned
    let cigar1: Cigar =
        [Op::new(Kind::Match, 80), Op::new(Kind::SoftClip, 10), Op::new(Kind::HardClip, 5)]
            .into_iter()
            .collect();
    let seq1: Sequence = vec![b'A'; 90].into();
    let qual1 = QualityScores::from(vec![30u8; 90]);

    // Read2: 80M → 80 sequence bases, 80 aligned
    let cigar2: Cigar = [Op::new(Kind::Match, 80)].into_iter().collect();
    let seq2: Sequence = vec![b'A'; 80].into();
    let qual2 = QualityScores::from(vec![30u8; 80]);

    let r1 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::MATE_REVERSE_COMPLEMENTED
                | Flags::FIRST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar1)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(51).unwrap())
        .set_template_length(130)
        .set_sequence(seq1)
        .set_quality_scores(qual1)
        .build();

    let r2 = RecordBuf::builder()
        .set_name("pair1")
        .set_flags(
            Flags::SEGMENTED
                | Flags::PROPERLY_SEGMENTED
                | Flags::REVERSE_COMPLEMENTED
                | Flags::LAST_SEGMENT,
        )
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(51).unwrap())
        .set_mapping_quality(mq)
        .set_cigar(cigar2)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(1).unwrap())
        .set_template_length(-130)
        .set_sequence(seq2)
        .set_quality_scores(qual2)
        .build();

    bld.add_record(r1);
    bld.add_record(r2);

    let bam = bld.to_temp_bam().unwrap();
    let prefix = dir.path().join("out");

    let mut options = opts(0, 0);
    options.dont_clip_overlapping_reads = false;
    let cmd = make_cmd(bam.path(), &bait_path, &target_path, &prefix, None, options);
    let m = run_and_read(&cmd);

    assert_eq!(m.on_target_bases, 130);
    assert_float_eq!(m.frac_exc_overlap, 30.0 / 160.0, 0.001);
}
