use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow, ensure};
use clap::{Args, ValueEnum};
use noodles::sam::Header;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::fasta::Fasta;
use crate::intervals::{Interval, Intervals};
use crate::math::{safe_div, safe_div_f};
use crate::metrics::{
    serialize_f64_2dp, serialize_f64_5dp, serialize_f64_6dp, tsv_writer, write_tsv,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sequence_dict::SequenceDictionary;
use crate::simd;

/// `(ref_start_0based, ref_end_exclusive, aligned_bases, alignment_blocks)`.
type AlignmentInfo = (u32, u32, u64, SmallVec<[(u32, u32); 8]>);

// ─── File suffixes ─────────────────────────────────────────────────────────────

/// File suffix for the summary metrics output.
pub const METRICS_SUFFIX: &str = ".hybcap-metrics.txt";

/// File suffix for the per-target coverage output.
pub const PER_TARGET_SUFFIX: &str = ".hybcap-per-target.txt";

/// File suffix for the uncompressed per-base coverage output.
pub const PER_BASE_SUFFIX: &str = ".hybcap-per-base.txt";

/// File suffix for the gzip-compressed per-base coverage output.
pub const PER_BASE_GZ_SUFFIX: &str = ".hybcap-per-base.txt.gz";

/// Output format for the per-base coverage file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum PerBaseCoverageFormat {
    /// Plain TSV (default).
    Uncompressed,
    /// Gzip-compressed TSV (suffix `.gz`).
    Compressed,
}

// ─── CLI struct ───────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the hybcap collector.
#[riker_derive::multi_options("hybcap", "Hybrid Capture Options")]
#[derive(Args, Debug, Clone)]
#[command()]
#[allow(clippy::struct_excessive_bools)]
pub struct HybCapOptions {
    /// Bait (probe) interval file (IntervalList or BED).
    #[arg(long, value_name = "FILE", value_parser = crate::commands::common::parse_existing_file)]
    pub baits: PathBuf,

    /// Target interval file (IntervalList or BED).
    #[arg(long, value_name = "FILE", value_parser = crate::commands::common::parse_existing_file)]
    pub targets: PathBuf,

    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,

    /// Optional name for the bait/target panel. Defaults to the bait filename stem.
    #[arg(long, value_name = "NAME")]
    pub panel_name: Option<String>,

    /// Distance (bp) to consider a read "near" a bait for classification.
    #[arg(long, default_value_t = HybCapOptions::DEFAULT_NEAR_DISTANCE)]
    pub near_distance: u64,

    /// Minimum mapping quality for a read's bases to count toward coverage.
    #[arg(long, default_value_t = HybCapOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapping_quality: u8,

    /// Minimum base quality for a base to count toward on-target coverage.
    #[arg(long, default_value_t = HybCapOptions::DEFAULT_MIN_BASE_QUALITY)]
    pub min_base_quality: u8,

    /// Disable clipping of overlapping bases from read pairs. By default
    /// overlapping bases are clipped to avoid double-counting coverage.
    #[arg(long, default_value_t = false)]
    pub dont_clip_overlapping_reads: bool,

    /// Count inserted and deleted bases as on-target (matching Picard INCLUDE_INDELS).
    #[arg(long, default_value_t = HybCapOptions::DEFAULT_INCLUDE_INDELS)]
    pub include_indels: bool,

    /// Output per-target coverage metrics.
    #[arg(long, default_value_t = false)]
    pub per_target_coverage: bool,

    /// Output per-base coverage values. Pass `compressed` for a gzipped
    /// `.hybcap-per-base.txt.gz` file or `uncompressed` (the default when the
    /// flag is given without a value) for a plain `.hybcap-per-base.txt`.
    #[arg(
        long,
        value_enum,
        value_name = "FORMAT",
        num_args = 0..=1,
        default_missing_value = "uncompressed",
    )]
    pub per_base_coverage: Option<PerBaseCoverageFormat>,
}

impl HybCapOptions {
    const DEFAULT_NEAR_DISTANCE: u64 = 250;
    const DEFAULT_MIN_MAPQ: u8 = 20;
    const DEFAULT_MIN_BASE_QUALITY: u8 = 20;
    const DEFAULT_INCLUDE_INDELS: bool = false;
}

impl Default for HybCapOptions {
    fn default() -> Self {
        Self {
            baits: PathBuf::new(),
            targets: PathBuf::new(),
            include_duplicates: false,
            panel_name: None,
            near_distance: Self::DEFAULT_NEAR_DISTANCE,
            min_mapping_quality: Self::DEFAULT_MIN_MAPQ,
            min_base_quality: Self::DEFAULT_MIN_BASE_QUALITY,
            dont_clip_overlapping_reads: false,
            include_indels: Self::DEFAULT_INCLUDE_INDELS,
            per_target_coverage: false,
            per_base_coverage: None,
        }
    }
}

/// Collect hybrid capture (HS) metrics from a BAM file.
///
/// Computes metrics describing the efficiency of hybrid selection (bait capture)
/// experiments, including on-target rate, coverage uniformity, and enrichment.
/// Requires bait and target interval files in IntervalList or BED format.
///
/// The input BAM must be coordinate-sorted (by reference, then position); a record
/// that arrives out of order aborts the run with an error.
///
/// Reads that fail vendor quality checks (SAM QC_FAIL flag) are excluded from
/// all computations.
///
/// Outputs are written to <prefix>.hybcap-metrics.txt, and optionally
/// <prefix>.hybcap-per-target.txt and either <prefix>.hybcap-per-base.txt
/// or <prefix>.hybcap-per-base.txt.gz (depending on
/// `--per-base-coverage`).
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker hybcap -i input.bam -o out_prefix --baits baits.bed --targets targets.bed
  riker hybcap -i input.bam -o out_prefix --baits baits.interval_list --targets targets.interval_list -r ref.fa
  riker hybcap -i input.bam -o out_prefix --baits baits.bed --targets targets.bed --per-target-coverage
  riker hybcap -i input.bam -o out_prefix --baits baits.bed --targets targets.bed --per-base-coverage compressed"
)]
pub struct HybCap {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,

    #[command(flatten)]
    pub options: HybCapOptions,
}

impl Command for HybCap {
    fn execute(&self, threads: Option<u8>) -> Result<()> {
        let plan = self.thread_plan(threads);
        let mut reader = AlignmentReader::open(
            &self.input.input,
            self.reference.reference.as_deref(),
            plan.decode_threads,
        )?;
        let sample = derive_sample(&self.input.input, reader.header());

        // Load reference if provided (needed for GC dropout)
        let fasta = match &self.reference.reference {
            Some(ref_path) => Some(Fasta::from_path(ref_path)?),
            None => None,
        };

        let mut collector = HybCapCollector::new(&self.output.output, fasta, sample, &self.options);
        let mut progress = ProgressLogger::new("hybcap", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

/// Accumulates hybrid capture metrics.
#[allow(clippy::struct_excessive_bools)]
pub struct HybCapCollector {
    // Output paths
    metrics_path: PathBuf,
    per_target_path: PathBuf,
    per_base_path: PathBuf,
    per_base_format: Option<PerBaseCoverageFormat>,

    // Input/config
    baits_path: PathBuf,
    targets_path: PathBuf,
    fasta: Option<Fasta>,
    sample: String,
    panel_name: String,

    // Options
    include_duplicates: bool,
    near_distance: u32,
    min_mapq: u8,
    min_bq: u8,
    clip_overlapping: bool,
    include_indels: bool,
    output_per_target: bool,

    // Initialized during initialize()
    dict: Option<SequenceDictionary>,
    // Per-contig sorted, non-overlapping interval lists (indexed by ref_id), swept
    // with monotonic cursors since reads arrive in coordinate order. Targets carry
    // their flat index into `target_coverages` as the third field; bait lists leave
    // it 0.
    target_ivs: Vec<Vec<(u32, u32, u32)>>,
    bait_ivs: Vec<Vec<(u32, u32, u32)>>,
    expanded_bait_ivs: Vec<Vec<(u32, u32, u32)>>,
    target_cursor: usize,
    bait_cursor: usize,
    expanded_cursor: usize,
    // Coordinate-sort guard: last (ref_id, start) seen; a change of ref_id also
    // triggers the per-contig cursor reset.
    last_ref_id: Option<usize>,
    last_start: u32,

    // Per-target coverage tracking (keyed by (ref_id, start, end))
    target_coverages: Vec<(Interval, TargetCoverage)>,
    // GC fractions per target (same order as target_coverages)
    target_gc: Vec<f64>,

    // Merged interval sets for territory calculations
    merged_bait_territory: u64,
    merged_target_territory: u64,
    raw_target_count: usize,
    genome_size: u64,

    // ── Read-level counters (PF reads only; non-PF are filtered in accept()) ──
    total_reads: u64,
    deduped_reads: u64,
    deduped_reads_aligned: u64,

    // ── Base-level counters ──
    total_bases: u64,
    bases_aligned: u64,
    deduped_bases_aligned: u64,

    // ── Bait classification (computed BEFORE dup filtering) ──
    on_bait_bases: u64,
    near_bait_bases: u64,
    off_bait_bases: u64,
    selected_pairs: u64,
    selected_unique_pairs: u64,

    // ── Target bases (computed AFTER all filtering) ──
    on_target_bases: u64,
    on_target_bases_from_pairs: u64,

    // ── Exclusion counters (raw counts, divided by bases_aligned in finish) ──
    exc_dupe_bases: u64,
    exc_mapq_bases: u64,
    exc_baseq_bases: u64,
    exc_overlap_bases: u64,
    exc_off_target_bases: u64,
}

impl HybCapCollector {
    /// Create a new collector. Intervals and overlap detectors are built in `initialize()`.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        prefix: &Path,
        fasta: Option<Fasta>,
        sample: String,
        options: &HybCapOptions,
    ) -> Self {
        let panel_name = options.panel_name.clone().unwrap_or_else(|| {
            options.baits.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string()
        });
        let per_base_suffix = match options.per_base_coverage {
            Some(PerBaseCoverageFormat::Compressed) => PER_BASE_GZ_SUFFIX,
            // Path is unused when per-base output is disabled, but matching the
            // historical default keeps the field's debug output friendly.
            Some(PerBaseCoverageFormat::Uncompressed) | None => PER_BASE_SUFFIX,
        };
        Self {
            metrics_path: super::command::output_path(prefix, METRICS_SUFFIX),
            per_target_path: super::command::output_path(prefix, PER_TARGET_SUFFIX),
            per_base_path: super::command::output_path(prefix, per_base_suffix),
            per_base_format: options.per_base_coverage,

            baits_path: options.baits.clone(),
            targets_path: options.targets.clone(),
            fasta,
            sample,
            panel_name,

            include_duplicates: options.include_duplicates,
            #[expect(
                clippy::cast_possible_truncation,
                reason = "near_distance is a small padding value that fits in u32"
            )]
            near_distance: options.near_distance as u32,
            min_mapq: options.min_mapping_quality,
            min_bq: options.min_base_quality,
            clip_overlapping: !options.dont_clip_overlapping_reads,
            include_indels: options.include_indels,
            output_per_target: options.per_target_coverage,

            dict: None,
            target_ivs: Vec::new(),
            bait_ivs: Vec::new(),
            expanded_bait_ivs: Vec::new(),
            target_cursor: 0,
            bait_cursor: 0,
            expanded_cursor: 0,
            last_ref_id: None,
            last_start: 0,

            target_coverages: Vec::new(),
            target_gc: Vec::new(),

            merged_bait_territory: 0,
            merged_target_territory: 0,
            raw_target_count: 0,
            genome_size: 0,

            total_reads: 0,
            deduped_reads: 0,
            deduped_reads_aligned: 0,

            total_bases: 0,
            bases_aligned: 0,
            deduped_bases_aligned: 0,

            on_bait_bases: 0,
            near_bait_bases: 0,
            off_bait_bases: 0,
            selected_pairs: 0,
            selected_unique_pairs: 0,

            on_target_bases: 0,
            on_target_bases_from_pairs: 0,

            exc_dupe_bases: 0,
            exc_mapq_bases: 0,
            exc_baseq_bases: 0,
            exc_overlap_bases: 0,
            exc_off_target_bases: 0,
        }
    }

    /// Return the contig name for a `ref_id`, or `?` if unknown. Used only in errors.
    fn contig_name(&self, ref_id: usize) -> &str {
        self.dict.as_ref().and_then(|d| d.get_by_index(ref_id)).map_or("?", |m| m.name())
    }

    /// Compute alignment span, aligned base count, and alignment blocks in a single
    /// CIGAR walk.  Returns `(ref_start_0based, ref_end_exclusive, aligned_bases, blocks)`
    /// or `None` if unmapped.  The blocks are `(block_start, block_end)` pairs in 0-based
    /// half-open coordinates that can be reused for bait classification without re-walking.
    fn alignment_span_and_blocks(record: &RikerRecord) -> Option<AlignmentInfo> {
        let start = record.alignment_start()?.get() - 1; // 1-based to 0-based
        let mut end = start;
        let mut aligned: u64 = 0;
        let mut blocks: SmallVec<[(u32, u32); 8]> = SmallVec::new();
        for op in record.cigar_ops() {
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    #[expect(
                        clippy::cast_possible_truncation,
                        reason = "genomic coordinates fit in u32"
                    )]
                    let block_start = end as u32;
                    end += op.len();
                    aligned += op.len() as u64;
                    #[expect(
                        clippy::cast_possible_truncation,
                        reason = "genomic coordinates fit in u32"
                    )]
                    let block_end = end as u32;
                    blocks.push((block_start, block_end));
                }
                Kind::Deletion | Kind::Skip => {
                    end += op.len();
                }
                _ => {}
            }
        }
        #[expect(clippy::cast_possible_truncation, reason = "genomic coordinates fit in u32")]
        Some((start as u32, end as u32, aligned, blocks))
    }

    /// Classify bases in a record as on-bait, near-bait, or off-bait.
    /// Called BEFORE duplicate filtering. Uses prefetched overlap results and
    /// precomputed alignment blocks (avoiding a CIGAR re-walk).
    fn classify_bait_bases(
        &mut self,
        aligned_bases: u64,
        on_expanded_bait: bool,
        bait_hits: &[(u32, u32)],
        blocks: &[(u32, u32)],
    ) {
        if !on_expanded_bait {
            self.off_bait_bases += aligned_bases;
            return;
        }

        // Classify each alignment block's bases as on-bait or near-bait using
        // prefetched bait_hits with range math — no CIGAR re-walk needed.
        let mut on_bait: u64 = 0;
        for &(block_start, block_end) in blocks {
            for &(bait_start, bait_end) in bait_hits {
                let overlap_start = block_start.max(bait_start);
                let overlap_end = block_end.min(bait_end);
                if overlap_start < overlap_end {
                    on_bait += u64::from(overlap_end - overlap_start);
                }
            }
        }

        self.on_bait_bases += on_bait;
        self.near_bait_bases += aligned_bases.saturating_sub(on_bait);
    }

    /// Determine the 0-based reference position at which overlap clipping starts for
    /// this read, or `None` if the read should not be clipped.
    ///
    /// Mirrors `count_overlapping_bases` logic: only the left-most read is clipped,
    /// on ties the second-of-pair is clipped.  Returns the mate's 0-based alignment
    /// start so the coverage walk can clip bases at or past that position.
    #[expect(clippy::cast_possible_truncation, reason = "genomic coordinates fit in u32")]
    fn compute_mate_clip_ref_pos(record: &RikerRecord) -> Option<u32> {
        let flags = record.flags();

        if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
            return None;
        }
        if record.reference_sequence_id() != record.mate_reference_sequence_id() {
            return None;
        }

        let alignment_start = record.alignment_start()?.get();
        let mate_start = record.mate_alignment_start()?.get();

        // Only clip the left-most read
        if mate_start < alignment_start {
            return None;
        }
        if mate_start == alignment_start && flags.is_first_segment() {
            return None;
        }

        // Convert mate_start from 1-based to 0-based
        Some((mate_start - 1) as u32)
    }

    /// Walk the CIGAR for a surviving (post-filter) read, accumulating on-target
    /// bases and per-target coverage depths.  Receives the read's overlapping
    /// targets as a pre-swept, start-sorted `(start, end, coverage-index)` slice.
    ///
    /// Overlap clipping is computed per M-block (not per-base) via the effective
    /// unclipped length.  Because target intervals are merged (non-overlapping),
    /// each M-block is decomposed into alternating off-target and on-target
    /// segments by range arithmetic: off-target runs use one auto-vectorized
    /// quality scan, on-target runs increment a contiguous depth slice — avoiding
    /// the per-base target-membership search.
    #[expect(clippy::too_many_lines, reason = "CIGAR walk is inherently detailed")]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "genomic coordinates and depths fit in u32/usize"
    )]
    fn walk_cigar_for_coverage(
        &mut self,
        record: &RikerRecord,
        start: u32,
        mate_clip_ref_pos: Option<u32>,
        is_mapped_pair: bool,
        cached_targets: &[(u32, u32, u32)],
    ) {
        let qual_bytes: &[u8] = record.quality_scores();
        let min_bq = self.min_bq;
        // Sentinel value: u32::MAX means "no clipping"
        let clip_ref_pos = mate_clip_ref_pos.unwrap_or(u32::MAX);

        // Scalar counters accumulate in locals and flush to `self` once at the
        // end.  A `self.field += 1` inside the per-base loops would force a
        // reload/store through `self` on every base, since the optimizer cannot
        // prove the depth-array writes do not alias the counter fields.
        let mut on_target: u64 = 0;
        let mut on_target_from_pairs: u64 = 0;
        let mut exc_baseq: u64 = 0;
        let mut exc_off_target: u64 = 0;
        let mut exc_overlap: u64 = 0;

        let mut ref_pos = start;
        let mut read_offset: usize = 0;
        let mut clipping_active = false;

        // Dedup read_count across M-blocks that touch the same target.  Targets are
        // encountered in non-decreasing order across a read (ref_pos is monotonic
        // and targets are sorted), so remembering only the last-counted local index
        // suffices — and, unlike a bitmask, imposes no cap on overlaps per read.
        let mut last_counted: usize = usize::MAX;

        for op in record.cigar_ops() {
            let op: Op = op;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let block_end = ref_pos + op.len() as u32;

                    // Compute the effective (unclipped) length once per M-block
                    // instead of testing `pos >= clip_ref_pos` on every base.
                    // Clipping always removes a suffix, so the kept part is a prefix.
                    let (effective_len, clipped_count) = if block_end <= clip_ref_pos {
                        (op.len(), 0)
                    } else if ref_pos >= clip_ref_pos {
                        (0, op.len())
                    } else {
                        let unclipped = (clip_ref_pos - ref_pos) as usize;
                        (unclipped, op.len() - unclipped)
                    };
                    if clipped_count > 0 {
                        exc_overlap += clipped_count as u64;
                        clipping_active = true;
                    }

                    // Clamp effective_len to the quality buffer; on a
                    // malformed record the CIGAR can claim more bases
                    // than there are quality bytes for, and the slices
                    // below would panic. Truncated bases are charged to
                    // the low-quality bucket (no quality info = treat as
                    // failing the threshold).
                    let avail = qual_bytes.len().saturating_sub(read_offset);
                    let scan_len = effective_len.min(avail);
                    let truncated = (effective_len - scan_len) as u64;
                    if truncated > 0 {
                        exc_baseq += truncated;
                    }

                    if scan_len > 0 {
                        // Quality slice for this block's kept (unclipped, in-bounds)
                        // prefix; sliced once for bounds-check-free access.
                        let quals = &qual_bytes[read_offset..read_offset + scan_len];
                        let r_end = ref_pos + scan_len as u32;
                        let mut cursor = ref_pos;

                        // Sweep the block against the sorted, non-overlapping
                        // targets: off-target gaps are counted with one vectorized
                        // low-quality scan; on-target runs increment a contiguous
                        // depth slice with only a per-base quality test.
                        for (local_idx, &(tgt_start, tgt_end, idx)) in
                            cached_targets.iter().enumerate()
                        {
                            if cursor >= r_end {
                                break;
                            }
                            // Off-target gap before this target.
                            let gap_end = tgt_start.min(r_end);
                            if cursor < gap_end {
                                let lo = (cursor - ref_pos) as usize;
                                let hi = (gap_end - ref_pos) as usize;
                                let low = simd::count_bases_lt_q(&quals[lo..hi], min_bq);
                                exc_baseq += low;
                                exc_off_target += (hi - lo) as u64 - low;
                                cursor = gap_end;
                            }
                            // On-target overlap [cursor, min(tgt_end, r_end)).
                            let on_end = tgt_end.min(r_end);
                            if cursor < on_end {
                                let lo = (cursor - ref_pos) as usize;
                                let seg_len = (on_end - cursor) as usize;
                                let depth_base = (cursor - tgt_start) as usize;
                                let seg_quals = &quals[lo..lo + seg_len];
                                let cov = &mut self.target_coverages[idx as usize].1;
                                let depths = &mut cov.hq_depths[depth_base..depth_base + seg_len];
                                let mut hq: u64 = 0;
                                for (d, &qual) in depths.iter_mut().zip(seg_quals) {
                                    if qual < min_bq {
                                        continue;
                                    }
                                    *d = d.saturating_add(1);
                                    hq += 1;
                                }
                                on_target += hq;
                                if is_mapped_pair {
                                    on_target_from_pairs += hq;
                                }
                                exc_baseq += seg_len as u64 - hq;
                                if hq > 0 && local_idx != last_counted {
                                    cov.read_count += 1;
                                    last_counted = local_idx;
                                }
                                cursor = on_end;
                            }
                        }

                        // Trailing off-target run after the last target (also the
                        // whole block when the read overlaps no targets).
                        if cursor < r_end {
                            let lo = (cursor - ref_pos) as usize;
                            let low = simd::count_bases_lt_q(&quals[lo..scan_len], min_bq);
                            exc_baseq += low;
                            exc_off_target += (scan_len - lo) as u64 - low;
                        }
                    }

                    ref_pos = block_end;
                    read_offset += op.len();
                }
                Kind::Insertion => {
                    if self.include_indels && !cached_targets.is_empty() {
                        // Inserted bases count toward on_target if ref_pos is on target.
                        let on_tgt = cached_targets.iter().any(|&(tgt_start, tgt_end, _)| {
                            ref_pos >= tgt_start && ref_pos < tgt_end
                        });
                        if on_tgt {
                            for i in 0..op.len() {
                                if clipping_active {
                                    exc_overlap += 1;
                                    continue;
                                }
                                let qual = qual_bytes.get(read_offset + i).copied().unwrap_or(0);
                                if qual >= min_bq {
                                    on_target += 1;
                                    if is_mapped_pair {
                                        on_target_from_pairs += 1;
                                    }
                                }
                            }
                        }
                    }
                    read_offset += op.len();
                }
                Kind::Deletion => {
                    if self.include_indels && !cached_targets.is_empty() {
                        for i in 0..op.len() {
                            let pos = ref_pos + i as u32;
                            for &(tgt_start, tgt_end, idx) in cached_targets {
                                if pos >= tgt_start && pos < tgt_end {
                                    self.target_coverages[idx as usize]
                                        .1
                                        .add_base((pos - tgt_start) as usize);
                                }
                            }
                        }
                    }
                    ref_pos += op.len() as u32;
                }
                Kind::SoftClip => {
                    read_offset += op.len();
                }
                Kind::Skip => {
                    ref_pos += op.len() as u32;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }

        self.on_target_bases += on_target;
        self.on_target_bases_from_pairs += on_target_from_pairs;
        self.exc_baseq_bases += exc_baseq;
        self.exc_off_target_bases += exc_off_target;
        self.exc_overlap_bases += exc_overlap;
    }

    /// Compute GC fractions for all targets using region queries to fetch only
    /// the sequence needed for each target.
    fn compute_all_target_gc(&mut self, dict: &SequenceDictionary) -> Vec<f64> {
        let Some(fasta) = &mut self.fasta else {
            return vec![0.0; self.target_coverages.len()];
        };

        let mut gc_fractions = Vec::with_capacity(self.target_coverages.len());

        for (iv, _) in &self.target_coverages {
            let contig_name = dict.get_by_index(iv.ref_id).map_or("", |m| m.name());
            let gc = match fasta.fetch(contig_name, u64::from(iv.start), u64::from(iv.end)) {
                Ok(bases) => Self::gc_fraction(&bases),
                Err(_) => 0.0,
            };
            gc_fractions.push(gc);
        }

        gc_fractions
    }

    /// Compute GC fraction from a sequence of bases.
    ///
    /// Non-ACGT bases (N and any IUPAC ambiguity codes) are counted in the
    /// denominator and contribute 0.5 to the GC numerator, treating them as
    /// maximally uncertain rather than ignoring them (which would inflate GC
    /// in N-rich regions) or treating them as non-GC (which would dilute it).
    fn gc_fraction(bases: &[u8]) -> f64 {
        let mut gc = 0u64;
        let mut ambiguous = 0u64;
        for &b in bases {
            match b {
                b'G' | b'C' | b'g' | b'c' => gc += 1,
                b'A' | b'T' | b'a' | b't' => {}
                _ => ambiguous += 1,
            }
        }
        safe_div_f(gc as f64 + ambiguous as f64 * 0.5, bases.len() as f64)
    }

    /// Estimate library size using the Lander-Waterman model (bisection method).
    /// Returns None if inputs are invalid (no pairs, no duplicates, etc.).
    fn estimate_library_size(read_pairs: u64, unique_pairs: u64) -> Option<u64> {
        if read_pairs == 0 || unique_pairs == 0 || unique_pairs >= read_pairs {
            return None;
        }

        let n = read_pairs as f64;
        let c = unique_pairs as f64;

        // f(x) = c/x - 1 + exp(-n/x)
        let f = |x: f64| -> f64 { c / x - 1.0 + (-n / x).exp() };

        // Initial bracket: m=1.0, M=100.0
        let mut lo = 1.0_f64;
        let mut hi = 100.0_f64;

        // Check that f(lo * c) >= 0
        if f(lo * c) < 0.0 {
            return None;
        }

        // Expand upper bound until f(hi * c) < 0
        while f(hi * c) >= 0.0 {
            hi *= 10.0;
            if hi > 1e15 {
                return None; // give up
            }
        }

        // Bisection: up to 40 iterations
        for _ in 0..40 {
            let mid = f64::midpoint(lo, hi);
            let val = f(mid * c);
            if val == 0.0 {
                break;
            }
            if val > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        #[expect(clippy::cast_sign_loss, reason = "result is always positive")]
        #[expect(clippy::cast_possible_truncation, reason = "library size fits in u64")]
        Some((c * (lo + hi) / 2.0) as u64)
    }

    /// Calculate the HS penalty for a given coverage goal.
    fn calculate_hs_penalty(
        library_size: Option<u64>,
        mean_coverage: f64,
        fold_80: f64,
        on_target_pct: f64,
        pairs: u64,
        unique_pairs: u64,
        coverage_goal: u32,
    ) -> f64 {
        let lib_size = match library_size {
            Some(ls) if ls > 0 => ls,
            _ => return 0.0,
        };

        if mean_coverage == 0.0 || on_target_pct == 0.0 || pairs == 0 || unique_pairs == 0 {
            return 0.0;
        }

        let unique_pair_goal_multiplier = (f64::from(coverage_goal) / mean_coverage) * fold_80;
        let mut pair_multiplier = unique_pair_goal_multiplier;
        let mut increment = 1.0_f64;
        let mut going_up = unique_pair_goal_multiplier >= 1.0;

        let mut final_pair_multiplier = pair_multiplier;

        for _ in 0..10_000 {
            let unique_pair_multiplier =
                Self::estimate_roi(lib_size, pair_multiplier, pairs, unique_pairs);

            let diff = (unique_pair_multiplier - unique_pair_goal_multiplier).abs();
            if diff / unique_pair_goal_multiplier <= 0.001 {
                final_pair_multiplier = pair_multiplier;
                break;
            }

            let overshot = if going_up {
                unique_pair_multiplier > unique_pair_goal_multiplier
            } else {
                unique_pair_multiplier < unique_pair_goal_multiplier
            };

            if overshot {
                increment /= 2.0;
                going_up = !going_up;
            }

            if going_up {
                pair_multiplier += increment;
            } else {
                pair_multiplier -= increment;
            }

            final_pair_multiplier = pair_multiplier;
        }

        let unique_fraction = (unique_pairs as f64 * unique_pair_goal_multiplier)
            / (pairs as f64 * final_pair_multiplier);

        if unique_fraction <= 0.0 {
            return 0.0;
        }

        (1.0 / unique_fraction) * fold_80 * (1.0 / on_target_pct)
    }

    /// Estimate ROI (return on investment) for additional sequencing.
    fn estimate_roi(library_size: u64, x: f64, pairs: u64, unique_pairs: u64) -> f64 {
        let l = library_size as f64;
        let n = pairs as f64;
        let c = unique_pairs as f64;
        l * (1.0 - (-(x * n) / l).exp()) / c
    }

    /// Compute AT and GC dropout from per-target GC fractions and coverage.
    fn compute_dropout(&self) -> (f64, f64) {
        if self.fasta.is_none() || self.target_coverages.is_empty() {
            return (0.0, 0.0);
        }

        // Bin targets by integer GC percentage (0-100)
        let mut target_bases_by_gc = [0u64; 101];
        let mut aligned_bases_by_gc = [0u64; 101];

        for (i, (iv, cov)) in self.target_coverages.iter().enumerate() {
            let gc = self.target_gc.get(i).copied().unwrap_or(0.0);
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "GC is 0.0-1.0, rounded to 0-100"
            )]
            let gc_bin = (gc * 100.0).round().min(100.0) as usize;
            target_bases_by_gc[gc_bin] += u64::from(iv.len());
            aligned_bases_by_gc[gc_bin] += cov.total_depth();
        }

        let total_target_bases: u64 = target_bases_by_gc.iter().sum();
        let total_aligned_bases: u64 = aligned_bases_by_gc.iter().sum();

        if total_target_bases == 0 || total_aligned_bases == 0 {
            return (0.0, 0.0);
        }

        let total_target_f = total_target_bases as f64;
        let total_aligned_f = total_aligned_bases as f64;

        let mut at_dropout = 0.0_f64;
        let mut gc_dropout = 0.0_f64;

        for gc_bin in 0..=100 {
            let target_pct = target_bases_by_gc[gc_bin] as f64 / total_target_f;
            let aligned_pct = aligned_bases_by_gc[gc_bin] as f64 / total_aligned_f;
            let dropout = (target_pct - aligned_pct).max(0.0) * 100.0;

            // AT dropout: bins 0-50 (inclusive)
            if gc_bin <= 50 {
                at_dropout += dropout;
            }
            // GC dropout: bins 50-100 (inclusive) — note bin 50 counted in BOTH
            if gc_bin >= 50 {
                gc_dropout += dropout;
            }
        }

        (at_dropout, gc_dropout)
    }
}

impl Collector for HybCapCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        let dict = SequenceDictionary::from(header);

        // Compute genome size from sequence dictionary
        self.genome_size = dict.iter().map(|s| s.length() as u64).sum();

        // Load bait and target intervals
        let raw_baits = Intervals::from_path(&self.baits_path, dict.clone())?;
        let raw_targets = Intervals::from_path(&self.targets_path, dict.clone())?;
        self.raw_target_count = raw_targets.raw_count();

        ensure!(
            raw_baits.count() > 0,
            "No bait intervals loaded from {}",
            self.baits_path.display()
        );
        ensure!(
            raw_targets.count() > 0,
            "No target intervals loaded from {}",
            self.targets_path.display()
        );

        // Merge overlapping intervals and compute territories
        let merged_baits = raw_baits.merged();
        let merged_targets = raw_targets.merged();
        self.merged_bait_territory = merged_baits.territory();
        self.merged_target_territory = merged_targets.territory();

        log::info!(
            "Loaded {} bait intervals ({} bp territory) and {} target intervals ({} bp territory)",
            merged_baits.count(),
            self.merged_bait_territory,
            merged_targets.count(),
            self.merged_target_territory
        );

        // Build per-contig sorted interval lists for the coordinate-order sweep.
        // Inputs are merged (sorted, non-overlapping), so each read's overlaps are a
        // contiguous run reachable by a monotonic cursor. Targets carry their flat
        // index into `target_coverages` (built in the same merged order below).
        let n_contigs = dict.len();
        self.target_ivs = bucket_by_contig(
            n_contigs,
            merged_targets.iter().enumerate().map(|(idx, iv)| {
                (iv.ref_id, iv.start, iv.end, u32::try_from(idx).expect("target index fits in u32"))
            }),
        );
        self.bait_ivs = bucket_by_contig(
            n_contigs,
            merged_baits.iter().map(|iv| (iv.ref_id, iv.start, iv.end, 0u32)),
        );

        // Expanded baits (padded by near_distance, re-merged) drive near-bait
        // classification; only run-nonemptiness is queried.
        let expanded_baits = merged_baits.padded(self.near_distance).merged();
        self.expanded_bait_ivs = bucket_by_contig(
            n_contigs,
            expanded_baits.iter().map(|iv| (iv.ref_id, iv.start, iv.end, 0u32)),
        );

        // Initialize per-target coverage arrays
        for iv in merged_targets.iter() {
            let len = iv.len() as usize;
            self.target_coverages.push((iv.clone(), TargetCoverage::new(len)));
        }

        // Compute GC fractions per target via one indexed region fetch each.
        // Targets are sparse relative to the genome, so per-region fetches read
        // far fewer bytes than loading whole contigs would.
        self.target_gc = self.compute_all_target_gc(&dict);

        self.dict = Some(dict);
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // ── Gate 0: Skip secondary alignments and non-PF (QC fail) reads ──
        if flags.is_secondary() || flags.is_qc_fail() {
            return Ok(());
        }

        let is_supplementary = flags.is_supplementary();
        let is_mapped = !flags.is_unmapped();
        let is_duplicate = flags.is_duplicate();

        // ── Read-level counters (primary only) ──
        if !is_supplementary {
            self.total_reads += 1;

            if !is_duplicate || self.include_duplicates {
                self.deduped_reads += 1;
                if is_mapped {
                    self.deduped_reads_aligned += 1;
                }
            }
        }

        // ── Compute alignment span, aligned bases, and alignment blocks in one CIGAR walk ──
        // Also used for base-level counters, bait classification, and coverage.
        let span_info = if is_mapped { Self::alignment_span_and_blocks(record) } else { None };

        // ── Base-level counters ──
        let read_len = record.sequence_len() as u64;
        if !is_supplementary {
            self.total_bases += read_len;
        }

        if let Some((_, _, aligned_bases, _)) = &span_info {
            self.bases_aligned += aligned_bases;
            if !is_duplicate || self.include_duplicates {
                self.deduped_bases_aligned += aligned_bases;
            }
        }

        // Need ref_id and alignment span for overlap detection
        let Some(ref_id) = record.reference_sequence_id() else {
            return Ok(()); // unmapped
        };

        let Some((start, end, aligned_bases, ref blocks)) = span_info else {
            return Ok(());
        };

        // ── Coordinate-sort gate + per-contig cursor reset ──
        // The cursor sweeps below require reads ordered by (ref_id, position).
        // Reject any regression (also catches a BAM whose header claims coordinate
        // sort but isn't) rather than silently emitting wrong metrics.
        if let Some(last_ref_id) = self.last_ref_id {
            if ref_id < last_ref_id || (ref_id == last_ref_id && start < self.last_start) {
                return Err(anyhow!(
                    "hybcap requires a coordinate-sorted BAM (ordered by reference, then \
                     position); encountered {}:{} after {}:{}. Sort with `samtools sort`.",
                    self.contig_name(ref_id),
                    u64::from(start) + 1,
                    self.contig_name(last_ref_id),
                    u64::from(self.last_start) + 1,
                ));
            }
            if ref_id != last_ref_id {
                self.target_cursor = 0;
                self.bait_cursor = 0;
                self.expanded_cursor = 0;
            }
        }
        self.last_ref_id = Some(ref_id);
        self.last_start = start;

        // ── Sweep the sorted per-contig interval lists (amortized O(1) per read vs.
        //    a from-scratch binary search, since reads arrive in coordinate order) ──
        let (tlo, thi) = sweep_run(&self.target_ivs[ref_id], &mut self.target_cursor, start, end);
        let cached_targets: SmallVec<[(u32, u32, u32); 4]> =
            self.target_ivs[ref_id][tlo..thi].iter().copied().collect();

        let (blo, bhi) = sweep_run(&self.bait_ivs[ref_id], &mut self.bait_cursor, start, end);
        let bait_hits: SmallVec<[(u32, u32); 4]> =
            self.bait_ivs[ref_id][blo..bhi].iter().map(|&(s, e, _)| (s, e)).collect();

        // Skip the expanded-bait sweep when the read already overlaps a bait (then
        // it's on expanded bait too); its cursor catches up lazily on the next read
        // that needs it, since prunes are start-driven.
        let on_expanded_bait = if bait_hits.is_empty() {
            let (elo, ehi) =
                sweep_run(&self.expanded_bait_ivs[ref_id], &mut self.expanded_cursor, start, end);
            ehi > elo
        } else {
            true
        };

        // ── Bait classification (BEFORE dup filtering, includes supplementary) ──
        self.classify_bait_bases(aligned_bases, on_expanded_bait, &bait_hits, blocks);

        // Count selected pairs: primary alignments only
        if !is_supplementary {
            let is_first_of_pair =
                flags.is_segmented() && flags.is_first_segment() && !flags.is_mate_unmapped();
            if is_first_of_pair && on_expanded_bait {
                self.selected_pairs += 1;
                if !is_duplicate || self.include_duplicates {
                    self.selected_unique_pairs += 1;
                }
            }
        }

        // ── Gate 2: Skip duplicates ──
        if is_duplicate && !self.include_duplicates {
            self.exc_dupe_bases += aligned_bases;
            return Ok(());
        }

        // ── Gate 3: Skip low MAPQ ──
        let mapq = record.mapping_quality().map_or(0, |mq| mq.get());
        if mapq < self.min_mapq {
            self.exc_mapq_bases += aligned_bases;
            return Ok(());
        }

        // ── Overlap clipping (computed inline during CIGAR walk) ──
        // Determine the reference position at which overlap clipping starts.
        // The left-most read of a pair is clipped at mate_alignment_start.
        let mate_clip_ref_pos =
            if self.clip_overlapping { Self::compute_mate_clip_ref_pos(record) } else { None };

        // ── Base-by-base CIGAR walk using prefetched targets ──
        let is_mapped_pair =
            flags.is_segmented() && !flags.is_unmapped() && !flags.is_mate_unmapped();
        self.walk_cigar_for_coverage(
            record,
            start,
            mate_clip_ref_pos,
            is_mapped_pair,
            &cached_targets,
        );

        Ok(())
    }

    #[expect(clippy::too_many_lines, reason = "metric computation is inherently sequential")]
    fn finish(&mut self) -> Result<()> {
        let dict = self.dict.as_ref().unwrap();

        // ── Compute coverage metrics from per-target depth arrays ──
        let mut depth_histogram: Vec<u64> = Vec::new(); // index = depth, value = count of bases
        let mut total_target_depth: u64 = 0;
        let mut min_depth = u64::MAX;
        let mut max_depth: u64 = 0;
        let mut zero_cvg_targets: usize = 0;

        for (iv, cov) in &self.target_coverages {
            let target_len = iv.len() as usize;
            let target_max = cov.max_depth();
            let target_min = cov.min_depth();

            if target_max == 0 {
                zero_cvg_targets += 1;
            }

            total_target_depth += cov.total_depth();
            min_depth = min_depth.min(target_min);
            max_depth = max_depth.max(target_max);

            // Accumulate into global depth histogram
            for &d in &cov.hq_depths {
                let d_idx = d as usize;
                if d_idx >= depth_histogram.len() {
                    depth_histogram.resize(d_idx + 1, 0);
                }
                depth_histogram[d_idx] += 1;
            }

            // If target has zero length, it doesn't contribute to coverage
            if target_len == 0 {
                zero_cvg_targets += 1;
            }
        }

        // Handle edge case where no targets have any coverage
        if min_depth == u64::MAX {
            min_depth = 0;
        }

        let mean_target_coverage = safe_div(total_target_depth, self.merged_target_territory);

        // Compute median from depth histogram
        let median_target_coverage = compute_median(&depth_histogram, self.merged_target_territory);

        // Compute fold_80_base_penalty = mean / 20th percentile
        let p20 = compute_percentile(&depth_histogram, self.merged_target_territory, 0.20);
        let fold_80_base_penalty = if p20 > 0.0 { mean_target_coverage / p20 } else { 0.0 };

        // Compute frac_target_bases at each threshold
        let bases_at_or_above = compute_bases_at_or_above(&depth_histogram);
        let frac_at_threshold = |threshold: u64| -> f64 {
            #[expect(clippy::cast_possible_truncation, reason = "coverage thresholds fit in usize")]
            let t = threshold as usize;
            if t < bases_at_or_above.len() {
                safe_div(bases_at_or_above[t], self.merged_target_territory)
            } else {
                0.0
            }
        };

        // ── GC / AT dropout ──
        let (at_dropout, gc_dropout) = self.compute_dropout();

        // ── HS library size and penalties ──
        let hs_library_size =
            Self::estimate_library_size(self.selected_pairs, self.selected_unique_pairs);

        let on_target_pct = safe_div(self.on_target_bases, self.deduped_bases_aligned);

        // Use on_target_bases_from_pairs for mean coverage in penalty calc (matching Picard)
        let mean_cov_for_penalty =
            safe_div(self.on_target_bases_from_pairs, self.merged_target_territory);

        let calc_penalty = |goal: u32| -> f64 {
            Self::calculate_hs_penalty(
                hs_library_size,
                mean_cov_for_penalty,
                fold_80_base_penalty,
                on_target_pct,
                self.selected_pairs,
                self.selected_unique_pairs,
                goal,
            )
        };

        // ── Bait classification fractions ──
        let total_bait_bases = self.on_bait_bases + self.near_bait_bases + self.off_bait_bases;
        let selected_bases = self.on_bait_bases + self.near_bait_bases;

        // ── Fold enrichment ──
        let on_bait_frac = safe_div(self.on_bait_bases, total_bait_bases);
        let bait_genome_frac =
            safe_div_f(self.merged_bait_territory as f64, self.genome_size as f64);
        let fold_enrichment = safe_div_f(on_bait_frac, bait_genome_frac);

        // ── Build metric struct ──
        let metric = HybCapMetric {
            sample: self.sample.clone(),
            panel_name: self.panel_name.clone(),
            genome_size: self.genome_size,
            bait_territory: self.merged_bait_territory,
            target_territory: self.merged_target_territory,
            bait_design_efficiency: safe_div(
                self.merged_target_territory,
                self.merged_bait_territory,
            ),
            total_reads: self.total_reads,
            deduped_reads: self.deduped_reads,
            deduped_reads_aligned: self.deduped_reads_aligned,
            total_bases: self.total_bases,
            bases_aligned: self.bases_aligned,
            deduped_bases_aligned: self.deduped_bases_aligned,
            on_bait_bases: self.on_bait_bases,
            near_bait_bases: self.near_bait_bases,
            off_bait_bases: self.off_bait_bases,
            on_target_bases: self.on_target_bases,
            on_target_bases_from_pairs: self.on_target_bases_from_pairs,
            selected_pairs: self.selected_pairs,
            selected_unique_pairs: self.selected_unique_pairs,
            frac_deduped_reads: safe_div(self.deduped_reads, self.total_reads),
            frac_deduped_reads_aligned: safe_div(self.deduped_reads_aligned, self.deduped_reads),
            frac_selected_bases: safe_div(selected_bases, total_bait_bases),
            frac_off_bait: safe_div(self.off_bait_bases, total_bait_bases),
            on_bait_vs_selected: safe_div(self.on_bait_bases, selected_bases),
            mean_bait_coverage: safe_div(self.on_bait_bases, self.merged_bait_territory),
            mean_target_coverage,
            median_target_coverage,
            max_target_coverage: max_depth,
            min_target_coverage: min_depth,
            frac_uncovered_targets: safe_div(zero_cvg_targets as u64, self.raw_target_count as u64),
            fold_enrichment,
            fold_80_base_penalty,
            frac_exc_dupe: safe_div(self.exc_dupe_bases, self.bases_aligned),
            frac_exc_mapq: safe_div(self.exc_mapq_bases, self.bases_aligned),
            frac_exc_overlap: safe_div(self.exc_overlap_bases, self.bases_aligned),
            frac_exc_baseq: safe_div(self.exc_baseq_bases, self.bases_aligned),
            frac_exc_off_target: safe_div(self.exc_off_target_bases, self.bases_aligned),
            frac_target_bases_1x: frac_at_threshold(1),
            frac_target_bases_10x: frac_at_threshold(10),
            frac_target_bases_20x: frac_at_threshold(20),
            frac_target_bases_30x: frac_at_threshold(30),
            frac_target_bases_50x: frac_at_threshold(50),
            frac_target_bases_100x: frac_at_threshold(100),
            frac_target_bases_250x: frac_at_threshold(250),
            frac_target_bases_500x: frac_at_threshold(500),
            frac_target_bases_1000x: frac_at_threshold(1000),
            at_dropout,
            gc_dropout,
            frac_usable_bases_on_bait: safe_div(self.on_bait_bases, self.total_bases),
            frac_usable_bases_on_target: safe_div(self.on_target_bases, self.total_bases),
            hs_library_size,
            hs_penalty_10x: calc_penalty(10),
            hs_penalty_20x: calc_penalty(20),
            hs_penalty_30x: calc_penalty(30),
            hs_penalty_40x: calc_penalty(40),
            hs_penalty_50x: calc_penalty(50),
            hs_penalty_100x: calc_penalty(100),
        };

        // ── Write main metrics ──
        log::info!("Writing metrics to {}", self.metrics_path.display());
        write_tsv(&self.metrics_path, &[metric])?;

        // ── Per-target coverage ──
        // Stream rows row-by-row so we don't build a Vec the size of the target
        // panel just to hand it to a writer.
        if self.output_per_target {
            log::info!("Writing per-target coverage to {}", self.per_target_path.display());
            let mut writer = tsv_writer::<PerTargetCoverage>(&self.per_target_path)?;
            for (i, (iv, cov)) in self.target_coverages.iter().enumerate() {
                let contig_name =
                    dict.get_by_index(iv.ref_id).map_or("?", |m| m.name()).to_string();
                let gc_frac = self.target_gc.get(i).copied().unwrap_or(0.0);
                let normalized = if mean_target_coverage > 0.0 {
                    cov.mean_depth() / mean_target_coverage
                } else {
                    0.0
                };
                writer.write(&PerTargetCoverage {
                    sample: self.sample.clone(),
                    chrom: contig_name,
                    start: u64::from(iv.start) + 1,
                    end: u64::from(iv.end),
                    length: u64::from(iv.len()),
                    name: iv.name().to_string(),
                    gc_frac,
                    mean_coverage: cov.mean_depth(),
                    normalized_coverage: normalized,
                    min_coverage: cov.min_depth(),
                    max_coverage: cov.max_depth(),
                    frac_0x: cov.frac_at_0x(),
                    read_count: cov.read_count,
                })?;
            }
            writer.flush()?;
        }

        // ── Per-base coverage ──
        // Stream rows row-by-row — at production target sizes this file can
        // contain tens of millions of entries. Path-driven `.gz` dispatch
        // (when `--per-base-coverage compressed`) is handled inside
        // [`tsv_writer`].
        if self.per_base_format.is_some() {
            log::info!("Writing per-base coverage to {}", self.per_base_path.display());
            let mut writer = tsv_writer::<PerBaseCoverage>(&self.per_base_path)?;
            for (iv, cov) in &self.target_coverages {
                let contig_name =
                    dict.get_by_index(iv.ref_id).map_or("?", |m| m.name()).to_string();
                let target_name = iv.name().to_string();
                for (offset, &depth) in cov.hq_depths.iter().enumerate() {
                    writer.write(&PerBaseCoverage {
                        sample: self.sample.clone(),
                        chrom: contig_name.clone(),
                        pos: u64::from(iv.start) + offset as u64 + 1, // 1-based output
                        target: target_name.clone(),
                        coverage: u64::from(depth),
                    })?;
                }
            }
            writer.flush()?;
        }

        Ok(())
    }

    fn name(&self) -> &'static str {
        "hybcap"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        RikerRecordRequirements::NONE
    }

    /// Relative per-record worker cost (see [`Collector::cost_hint`]); ~140,
    /// the heaviest collector — dual bait+target per-base pileup. Measured on an
    /// exome BAM (samply worker-compute) at roughly 2x `basic` per read (basic's
    /// hint is 76, alignment's 55), which lands it near 140 on the shared hint
    /// scale. Approximate — per-read cost varies with read length across data
    /// sets — but its ordering as the heaviest is robust.
    fn cost_hint(&self) -> u32 {
        140
    }
}

// ─── Metric structs ───────────────────────────────────────────────────────────

/// Hybrid capture (HS) sequencing metrics.
#[derive(Debug, Serialize, Deserialize, MetricDocs, Clone)]
pub struct HybCapMetric {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Name of the bait/target panel.
    pub panel_name: String,
    /// Total number of bases in the reference genome.
    pub genome_size: u64,
    /// Total number of unique bases covered by bait intervals.
    pub bait_territory: u64,
    /// Total number of unique bases covered by target intervals.
    pub target_territory: u64,
    /// Ratio of target territory to bait territory.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub bait_design_efficiency: f64,
    /// Total primary reads in the BAM file (excludes non-PF reads).
    pub total_reads: u64,
    /// Non-duplicate reads.
    pub deduped_reads: u64,
    /// Non-duplicate aligned reads.
    pub deduped_reads_aligned: u64,
    /// Total bases (read length × primary reads).
    pub total_bases: u64,
    /// Total aligned bases (from primary + supplementary).
    pub bases_aligned: u64,
    /// Non-duplicate aligned bases.
    pub deduped_bases_aligned: u64,
    /// Aligned bases falling on bait intervals (includes duplicates).
    pub on_bait_bases: u64,
    /// Aligned bases near but not on baits (includes duplicates).
    pub near_bait_bases: u64,
    /// Aligned bases not near any bait (includes duplicates).
    pub off_bait_bases: u64,
    /// High-quality non-duplicate bases on target intervals.
    pub on_target_bases: u64,
    /// High-quality non-duplicate bases on target from mapped pairs.
    pub on_target_bases_from_pairs: u64,
    /// First-of-pair reads overlapping any bait (both ends mapped).
    pub selected_pairs: u64,
    /// Same as selected_pairs but excluding duplicates.
    pub selected_unique_pairs: u64,
    /// Fraction of total reads that are non-duplicate.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_deduped_reads: f64,
    /// Fraction of non-duplicate reads that are aligned.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_deduped_reads_aligned: f64,
    /// Fraction of aligned bases that are on or near a bait.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_selected_bases: f64,
    /// Fraction of aligned bases that are off-bait.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_off_bait: f64,
    /// Fraction of selected bases that are strictly on (not near) bait.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub on_bait_vs_selected: f64,
    /// Mean coverage depth across bait bases (from aligned bases incl. duplicates).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_bait_coverage: f64,
    /// Mean coverage depth across target bases (from HQ non-dup bases).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_target_coverage: f64,
    /// Median per-base coverage depth across target bases.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub median_target_coverage: f64,
    /// Maximum per-base coverage depth across target bases.
    pub max_target_coverage: u64,
    /// Minimum per-base coverage depth across target bases.
    pub min_target_coverage: u64,
    /// Fraction of raw (pre-merge) targets with zero coverage.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_uncovered_targets: f64,
    /// Fold enrichment of on-bait bases vs. random genomic distribution.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub fold_enrichment: f64,
    /// Mean / 20th-percentile coverage depth (uniformity penalty).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub fold_80_base_penalty: f64,
    /// Fraction of aligned bases excluded because the read was a duplicate.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_exc_dupe: f64,
    /// Fraction of aligned bases excluded due to low mapping quality.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_exc_mapq: f64,
    /// Fraction of aligned bases excluded due to read-pair overlap clipping.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_exc_overlap: f64,
    /// Fraction of aligned bases excluded due to low base quality.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_exc_baseq: f64,
    /// Fraction of aligned bases excluded because they were off-target (but high quality).
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_exc_off_target: f64,
    /// Fraction of target bases with coverage >= 1x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_1x: f64,
    /// Fraction of target bases with coverage >= 10x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_10x: f64,
    /// Fraction of target bases with coverage >= 20x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_20x: f64,
    /// Fraction of target bases with coverage >= 30x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_30x: f64,
    /// Fraction of target bases with coverage >= 50x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_50x: f64,
    /// Fraction of target bases with coverage >= 100x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_100x: f64,
    /// Fraction of target bases with coverage >= 250x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_250x: f64,
    /// Fraction of target bases with coverage >= 500x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_500x: f64,
    /// Fraction of target bases with coverage >= 1000x.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_target_bases_1000x: f64,
    /// AT dropout: under-representation of AT-rich targets.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub at_dropout: f64,
    /// GC dropout: under-representation of GC-rich targets.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub gc_dropout: f64,
    /// Fraction of total bases that are on-bait (including duplicates in denominator).
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_usable_bases_on_bait: f64,
    /// Fraction of total bases that are on-target (HQ, non-dup numerator, all bases denominator).
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_usable_bases_on_target: f64,
    /// Estimated library size using Lander-Waterman model on selected pairs.
    #[serde(serialize_with = "serialize_opt_u64")]
    pub hs_library_size: Option<u64>,
    /// Fold sequencing needed to reach 80% of targets at 10x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_10x: f64,
    /// Fold sequencing needed to reach 80% of targets at 20x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_20x: f64,
    /// Fold sequencing needed to reach 80% of targets at 30x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_30x: f64,
    /// Fold sequencing needed to reach 80% of targets at 40x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_40x: f64,
    /// Fold sequencing needed to reach 80% of targets at 50x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_50x: f64,
    /// Fold sequencing needed to reach 80% of targets at 100x.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub hs_penalty_100x: f64,
}

/// Per-target coverage metrics (one row per merged target interval).
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct PerTargetCoverage {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Contig name.
    pub chrom: String,
    /// 1-based start position.
    pub start: u64,
    /// 1-based inclusive end position.
    pub end: u64,
    /// Interval length in bases.
    pub length: u64,
    /// Interval name.
    pub name: String,
    /// GC fraction of the target region (0.0 if no reference provided).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_frac: f64,
    /// Mean high-quality coverage depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_coverage: f64,
    /// Mean coverage normalized by overall mean target coverage.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub normalized_coverage: f64,
    /// Minimum per-base coverage.
    pub min_coverage: u64,
    /// Maximum per-base coverage.
    pub max_coverage: u64,
    /// Fraction of bases with zero coverage.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_0x: f64,
    /// Number of reads overlapping this target.
    pub read_count: u64,
}

/// Per-base coverage (one row per target base position).
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct PerBaseCoverage {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Contig name.
    pub chrom: String,
    /// 1-based position.
    pub pos: u64,
    /// Name of the target interval containing this position.
    pub target: String,
    /// High-quality coverage depth at this position.
    pub coverage: u64,
}

// ─── Per-target depth tracking ────────────────────────────────────────────────

/// Tracks per-base depth for a single target interval.
struct TargetCoverage {
    /// Per-base high-quality depth (indexed by offset from interval start).
    /// Uses `u32` instead of `u64` — max realistic depth is ~30,000x, well within
    /// `u32::MAX` (~4.3 billion).  For 46 Mbp of targets this saves ~185 MB.
    hq_depths: Vec<u32>,
    /// Number of reads overlapping this target (for per-target output).
    read_count: u64,
}

impl TargetCoverage {
    fn new(len: usize) -> Self {
        Self { hq_depths: vec![0; len], read_count: 0 }
    }

    fn add_base(&mut self, offset: usize) {
        if offset < self.hq_depths.len() {
            self.hq_depths[offset] = self.hq_depths[offset].saturating_add(1);
        }
    }

    fn total_depth(&self) -> u64 {
        self.hq_depths.iter().map(|&d| u64::from(d)).sum()
    }

    fn min_depth(&self) -> u64 {
        self.hq_depths.iter().copied().min().map_or(0, u64::from)
    }

    fn max_depth(&self) -> u64 {
        self.hq_depths.iter().copied().max().map_or(0, u64::from)
    }

    fn mean_depth(&self) -> f64 {
        if self.hq_depths.is_empty() {
            0.0
        } else {
            self.total_depth() as f64 / self.hq_depths.len() as f64
        }
    }

    fn frac_at_0x(&self) -> f64 {
        if self.hq_depths.is_empty() {
            0.0
        } else {
            let zeros = self.hq_depths.iter().filter(|&&d| d == 0).count();
            zeros as f64 / self.hq_depths.len() as f64
        }
    }
}

// ─── Interval sweep helpers ───────────────────────────────────────────────────

/// Bucket `(ref_id, start, end, payload)` intervals into per-contig lists indexed
/// by `ref_id`, each sorted by start.  Inputs come from merged interval sets, so
/// the per-contig lists are non-overlapping.
fn bucket_by_contig(
    n_contigs: usize,
    intervals: impl Iterator<Item = (usize, u32, u32, u32)>,
) -> Vec<Vec<(u32, u32, u32)>> {
    let mut per_contig = vec![Vec::new(); n_contigs];
    for (ref_id, start, end, payload) in intervals {
        per_contig[ref_id].push((start, end, payload));
    }
    for list in &mut per_contig {
        list.sort_unstable_by_key(|&(start, _, _)| start);
    }
    per_contig
}

/// Advance `cursor` past intervals in the sorted, non-overlapping `ivs` that end at
/// or before `start`, then return the `[lo, hi)` index range of intervals
/// overlapping `[start, end)`.  Correct only when `start` is monotonically
/// non-decreasing across calls (coordinate-sorted reads); the caller enforces that.
fn sweep_run(ivs: &[(u32, u32, u32)], cursor: &mut usize, start: u32, end: u32) -> (usize, usize) {
    while *cursor < ivs.len() && ivs[*cursor].1 <= start {
        *cursor += 1;
    }
    let lo = *cursor;
    let mut hi = lo;
    while hi < ivs.len() && ivs[hi].0 < end {
        hi += 1;
    }
    (lo, hi)
}

// ─── Histogram helper functions ───────────────────────────────────────────────

/// Compute the median from a depth histogram.
fn compute_median(histogram: &[u64], total_bases: u64) -> f64 {
    if total_bases == 0 {
        return 0.0;
    }

    let mid = total_bases / 2;
    let mut cumulative: u64 = 0;

    for (depth, &count) in histogram.iter().enumerate() {
        cumulative += count;
        if cumulative > mid {
            return depth as f64;
        }
    }

    0.0
}

/// Compute the given percentile (0.0-1.0) from a depth histogram.
fn compute_percentile(histogram: &[u64], total_bases: u64, percentile: f64) -> f64 {
    if total_bases == 0 {
        return 0.0;
    }

    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "percentile target is always in range"
    )]
    let target = (total_bases as f64 * percentile) as u64;
    let mut cumulative: u64 = 0;

    for (depth, &count) in histogram.iter().enumerate() {
        cumulative += count;
        if cumulative >= target {
            return depth as f64;
        }
    }

    0.0
}

/// Compute bases-at-or-above for each depth. Returns a vec where index i
/// gives the number of bases with depth >= i.
fn compute_bases_at_or_above(histogram: &[u64]) -> Vec<u64> {
    if histogram.is_empty() {
        return Vec::new();
    }

    let mut result = vec![0u64; histogram.len()];
    let mut cumulative: u64 = 0;

    for i in (0..histogram.len()).rev() {
        cumulative += histogram[i];
        result[i] = cumulative;
    }

    result
}

// ─── Serializer for Option<u64> ───────────────────────────────────────────────

#[allow(clippy::ref_option)]
fn serialize_opt_u64<S>(value: &Option<u64>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match value {
        Some(v) => serializer.serialize_u64(*v),
        None => serializer.serialize_str(""),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_median_simple() {
        // 5 bases at depth 0, 10 bases at depth 1
        let hist = vec![5, 10];
        assert!((compute_median(&hist, 15) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_median_empty() {
        assert!((compute_median(&[], 0)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_percentile() {
        // 10 bases at depth 0, 10 at depth 1, 10 at depth 2
        let hist = vec![10, 10, 10];
        assert!((compute_percentile(&hist, 30, 0.20)).abs() < f64::EPSILON); // 20th pctile = depth 0
        assert!((compute_percentile(&hist, 30, 0.50) - 1.0).abs() < f64::EPSILON); // 50th pctile = depth 1
    }

    #[test]
    fn test_gc_fraction_pure_acgt() {
        // 4 G/C out of 8 → 0.5.
        assert!((HybCapCollector::gc_fraction(b"ACGTACGT") - 0.5).abs() < f64::EPSILON);
        // Lowercase is treated identically.
        assert!((HybCapCollector::gc_fraction(b"acgtacgt") - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_gc_fraction_iupac_contributes_half() {
        // N and every other IUPAC ambiguity code are treated uniformly:
        // each contributes 0.5 to the numerator and 1 to the denominator.
        // 4 unambiguous ACGT (2 G/C) + 8 ambiguous → (2 + 8*0.5) / 12 = 0.5.
        let bases = b"ACGTNSWRYKMB";
        assert!((HybCapCollector::gc_fraction(bases) - 0.5).abs() < f64::EPSILON);

        // All-N pulls GC toward 0.5 regardless of case.
        assert!((HybCapCollector::gc_fraction(b"NNNNnnnn") - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_gc_fraction_empty() {
        // safe_div_f returns 0.0 for an empty slice.
        assert!(HybCapCollector::gc_fraction(b"").abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_bases_at_or_above() {
        let hist = vec![5, 10, 3, 2];
        let result = compute_bases_at_or_above(&hist);
        assert_eq!(result, vec![20, 15, 5, 2]);
    }

    #[test]
    fn test_estimate_library_size() {
        // Basic sanity: more pairs than unique → some library size
        let result = HybCapCollector::estimate_library_size(1000, 900);
        assert!(result.is_some());
        let ls = result.unwrap();
        assert!(ls > 900, "library size {ls} should be > unique pairs");

        // No duplication → None
        assert!(HybCapCollector::estimate_library_size(1000, 1000).is_none());

        // No pairs → None
        assert!(HybCapCollector::estimate_library_size(0, 0).is_none());
    }

    #[test]
    fn test_estimate_roi() {
        let lib_size = 10_000;
        let pairs = 1000;
        let unique = 900;

        // 1x sequencing should give ROI close to 1.0
        let roi = HybCapCollector::estimate_roi(lib_size, 1.0, pairs, unique);
        assert!((roi - 1.0).abs() < 0.2, "ROI at 1x = {roi}");

        // More sequencing should give higher ROI
        let roi2 = HybCapCollector::estimate_roi(lib_size, 2.0, pairs, unique);
        assert!(roi2 > roi, "ROI at 2x ({roi2}) should be > ROI at 1x ({roi})");
    }

    #[test]
    fn test_hs_penalty_basic() {
        let penalty = HybCapCollector::calculate_hs_penalty(
            Some(100_000),
            50.0,   // mean coverage
            1.5,    // fold_80
            0.8,    // on_target_pct
            10_000, // pairs
            9_000,  // unique pairs
            30,     // coverage goal
        );
        assert!(penalty > 0.0, "penalty should be positive, got {penalty}");
    }

    #[test]
    fn test_hs_penalty_no_library_size() {
        let penalty =
            HybCapCollector::calculate_hs_penalty(None, 50.0, 1.5, 0.8, 10_000, 9_000, 30);
        assert!(penalty.abs() < f64::EPSILON);
    }
}
