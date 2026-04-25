use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bitvec::prelude::*;
use clap::Args;
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use noodles::sam::alignment::record::cigar::op::Kind;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use crate::counter::Counter;
use crate::fasta::Fasta;
use crate::intervals::Intervals;
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, write_tsv};
use crate::plotting::{FG_BLUE, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_plot_pdf};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::mate_buffer::{MateBuffer, MateProbeFields, Peek};
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sequence_dict::SequenceDictionary;

// ─── File suffixes ─────────────────────────────────────────────────────────────

/// File suffix for the summary metrics output.
pub const METRICS_SUFFIX: &str = ".wgs-metrics.txt";

/// File suffix for the per-depth coverage histogram output.
pub const COVERAGE_SUFFIX: &str = ".wgs-coverage.txt";

/// File suffix for the coverage depth histogram plot.
pub const PLOT_SUFFIX: &str = ".wgs-coverage.pdf";

// ─── Options ──────────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the WGS collector.
#[riker_derive::multi_options("wgs", "WGS Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct WgsOptions {
    /// Interval file (IntervalList or BED) to restrict analysis.
    #[arg(short = 'L', long, value_name = "FILE")]
    pub intervals: Option<PathBuf>,

    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,

    /// Exclude unpaired reads from coverage.
    ///
    /// Reads without a mapped mate are excluded by default to match
    /// Picard's behavior.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_EXCLUDE_UNPAIRED)]
    pub exclude_unpaired_reads: bool,

    /// Minimum mapping quality for a base to be counted toward coverage.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Minimum base quality for a base to be counted toward coverage.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_MIN_BQ)]
    pub min_bq: u8,

    /// Maximum depth reported in the metrics and histogram.
    ///
    /// Per-position depths saturate at this value; bases arriving at a
    /// position that has already hit the cap are counted as capped exclusions.
    /// Bounded by u16 (max 65535) since the per-position depth array stores
    /// each slot as a u16.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_COVERAGE_CAP)]
    pub coverage_cap: u16,
}

impl WgsOptions {
    const DEFAULT_EXCLUDE_UNPAIRED: bool = true;
    const DEFAULT_MIN_MAPQ: u8 = 20;
    const DEFAULT_MIN_BQ: u8 = 20;
    const DEFAULT_COVERAGE_CAP: u16 = 250;
}

impl Default for WgsOptions {
    fn default() -> Self {
        Self {
            intervals: None,
            include_duplicates: false,
            exclude_unpaired_reads: Self::DEFAULT_EXCLUDE_UNPAIRED,
            min_mapq: Self::DEFAULT_MIN_MAPQ,
            min_bq: Self::DEFAULT_MIN_BQ,
            coverage_cap: Self::DEFAULT_COVERAGE_CAP,
        }
    }
}

// ─── CLI struct ───────────────────────────────────────────────────────────────

/// Collect whole-genome sequencing coverage metrics from a BAM file.
///
/// Counts aligned bases at each genomic position, accounting for
/// overlapping read pairs, and produces coverage depth statistics
/// including mean, median, standard deviation, and MAD. Also reports
/// the fraction of bases excluded by various filters. Outputs are
/// written to <prefix>.wgs-metrics.txt, <prefix>.wgs-coverage.txt,
/// and <prefix>.wgs-coverage.pdf.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker wgs -i input.bam -o out_prefix -R ref.fa
  riker wgs -i input.bam -o out_prefix -R ref.fa -L intervals.bed
  riker wgs -i input.bam -o out_prefix -R ref.fa --coverage-cap 500"
)]
pub struct Wgs {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub options: WgsOptions,
}

impl Command for Wgs {
    /// # Errors
    /// Returns an error if the BAM or reference cannot be read, or if output cannot be written.
    fn execute(&self) -> Result<()> {
        let mut reader = AlignmentReader::open(&self.input.input, Some(&self.reference.reference))?;
        let reference = Fasta::from_path(&self.reference.reference)?;

        let mut collector =
            WgsCollector::new(&self.input.input, &self.output.output, reference, &self.options)?;

        collector.initialize(reader.header())?;
        let mut progress = ProgressLogger::new("wgs", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)?;
        collector.finish()
    }
}

// ─── Collector ────────────────────────────────────────────────────────────────

/// Accumulates per-position coverage statistics using a contig-sized depth
/// array and a per-pair mate buffer for overlap detection.
pub struct WgsCollector {
    // Output paths
    metrics_path: PathBuf,
    coverage_path: PathBuf,
    plot_path: PathBuf,
    plot_title: String,

    // Input path (needed for sample name derivation in initialize)
    input_path: PathBuf,

    // Per-contig depth array (one u16 per genomic position, reused across contigs).
    depth: ContigDepth,

    // Buffer of already-seen mate records whose overlapping mate hasn't arrived yet.
    mate_buffer: MateBuffer<CachedMate>,

    // Config
    reference: Fasta,
    intervals: Option<Intervals>,
    include_duplicates: bool,
    include_unpaired: bool,
    min_mapq: u8,
    min_bq: u8,
    coverage_cap: u16,

    // BAM contig metadata
    dict: Option<SequenceDictionary>,

    // Per-contig working state
    current_ref_id: Option<usize>,
    current_non_n: BitVec,
    processed_contigs: HashSet<usize>,

    // Global accumulators
    genome_territory: u64,
    depth_histogram: Vec<u64>, // index = capped depth (0..=coverage_cap)
    bases_excl_mapq: u64,
    bases_excl_dupe: u64,
    bases_excl_unpaired: u64,
    bases_excl_baseq: u64,
    bases_excl_overlap: u64,
    bases_excl_capped: u64,

    sample: String,
}

impl WgsCollector {
    /// Create a new collector. Output paths are derived from `prefix` by appending
    /// suffixes. If `options.intervals` is set, the interval file is parsed using
    /// the FASTA index's sequence dictionary. The sample name and contig names are
    /// populated during [`Collector::initialize`] once the BAM header is available.
    ///
    /// # Errors
    /// Returns an error if the interval file cannot be read or parsed.
    pub fn new(
        input: &Path,
        prefix: &Path,
        reference: Fasta,
        options: &WgsOptions,
    ) -> Result<Self> {
        let intervals = options
            .intervals
            .as_ref()
            .map(|path| Intervals::from_path(path, reference.dict().clone()))
            .transpose()?;

        let metrics_path = super::command::output_path(prefix, METRICS_SUFFIX);
        let coverage_path = super::command::output_path(prefix, COVERAGE_SUFFIX);
        let plot_path = super::command::output_path(prefix, PLOT_SUFFIX);
        let hist_len = options.coverage_cap as usize + 1;
        Ok(Self {
            metrics_path,
            coverage_path,
            plot_path,
            plot_title: String::new(),
            input_path: input.to_path_buf(),
            depth: ContigDepth::new(options.coverage_cap),
            mate_buffer: MateBuffer::new(),
            reference,
            intervals,
            include_duplicates: options.include_duplicates,
            include_unpaired: !options.exclude_unpaired_reads,
            min_mapq: options.min_mapq,
            min_bq: options.min_bq,
            coverage_cap: options.coverage_cap,
            dict: None,
            current_ref_id: None,
            current_non_n: BitVec::EMPTY,
            processed_contigs: HashSet::new(),
            genome_territory: 0,
            depth_histogram: vec![0u64; hist_len],
            bases_excl_mapq: 0,
            bases_excl_dupe: 0,
            bases_excl_unpaired: 0,
            bases_excl_baseq: 0,
            bases_excl_overlap: 0,
            bases_excl_capped: 0,
            sample: String::new(),
        })
    }

    /// Finalize the last contig, process any un-pileup'd contigs, and write output.
    ///
    /// # Errors
    /// Returns an error if reference contigs cannot be loaded or output files cannot be written.
    #[expect(clippy::too_many_lines, reason = "sequential finalization logic")]
    fn finish_metrics(&mut self) -> Result<()> {
        // Finalize whatever contig was being processed last.
        self.finalize_contig();
        self.current_ref_id = None;

        // Process contigs that had zero reads (never appeared in pileup).
        // `dict` is set in `initialize`, which always runs before `finish`.
        let dict = self.dict.as_ref().unwrap();
        let n_contigs = dict.len();
        for ref_id in 0..n_contigs {
            if self.processed_contigs.contains(&ref_id) {
                continue;
            }

            // Skip if no intervals specified for this contig when intervals are enabled.
            if let Some(ref intervals) = self.intervals
                && !intervals.has_contig(ref_id)
            {
                continue;
            }

            let name = dict[ref_id].name();
            match self.reference.load_contig(name, false) {
                Err(e) => {
                    log::warn!("wgs: could not load contig '{name}': {e}");
                }
                Ok(seq) => {
                    let bv = build_non_n_bitvec(&seq);
                    let non_n = self.count_eligible_positions(ref_id, &bv);
                    self.genome_territory += non_n;
                    self.depth_histogram[0] += non_n;
                }
            }
        }

        let total = self.genome_territory;
        if total == 0 {
            log::warn!("wgs: genome territory is zero; writing empty metrics");
            write_tsv(&self.metrics_path, &[] as &[WgsMetrics])?;
            write_tsv(&self.coverage_path, &[] as &[WgsCoverageEntry])?;
            return Ok(());
        }

        // ── Build Counter from depth histogram for statistics helpers ────────
        let hist_counter: Counter<u64> = self
            .depth_histogram
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c > 0)
            .map(|(d, &c)| (d as u64, c))
            .collect();

        // ── Compute mean and SD from the depth histogram ────────────────────
        let (mean, sd, sum_depth) = {
            let n = total as f64;
            let mut sum = 0.0_f64;
            let mut sum_sq = 0.0_f64;
            for (depth, &count) in self.depth_histogram.iter().enumerate() {
                if count > 0 {
                    let d = depth as f64;
                    let c = count as f64;
                    sum += d * c;
                    sum_sq += d * d * c;
                }
            }
            let mean = sum / n;
            let variance = if total > 1 { (sum_sq - n * mean * mean) / (n - 1.0) } else { 0.0 };
            (mean, variance.max(0.0).sqrt(), sum)
        };

        let (median, mad) = hist_counter.median_and_mad();

        // ── Reverse cumulative sums (used for frac_at, coverage histogram) ──
        let n = self.depth_histogram.len();
        let mut at_or_above = vec![0u64; n];
        let mut running: u64 = 0;
        for d in (0..n).rev() {
            running += self.depth_histogram[d];
            at_or_above[d] = running;
        }

        let frac_at = |nx: usize| -> f64 {
            if nx >= n { 0.0 } else { at_or_above[nx] as f64 / total as f64 }
        };

        // ── Fold penalties ────────────────────────────────────────────────────
        let fold_penalty = |lower_frac: f64| -> f64 {
            let pct_depth = percentile_from_hist(&self.depth_histogram, lower_frac, total);
            if pct_depth == 0.0 { 0.0 } else { mean / pct_depth }
        };

        // ── Exclusion fractions ───────────────────────────────────────────────
        let total_excl = u128::from(
            self.bases_excl_mapq
                + self.bases_excl_dupe
                + self.bases_excl_unpaired
                + self.bases_excl_baseq
                + self.bases_excl_overlap
                + self.bases_excl_capped,
        );

        let total_raw = total_excl as f64 + sum_depth;

        let frac_excl =
            |n: u64| -> f64 { if total_raw == 0.0 { 0.0 } else { n as f64 / total_raw } };

        let frac_excl_total = if total_raw == 0.0 { 0.0 } else { total_excl as f64 / total_raw };

        let metrics = WgsMetrics {
            sample: self.sample.clone(),
            genome_territory: total,
            mean_coverage: mean,
            sd_coverage: sd,
            median_coverage: median,
            mad_coverage: mad,
            frac_excluded_mapq: frac_excl(self.bases_excl_mapq),
            frac_excluded_dupe: frac_excl(self.bases_excl_dupe),
            frac_excluded_unpaired: frac_excl(self.bases_excl_unpaired),
            frac_excluded_baseq: frac_excl(self.bases_excl_baseq),
            frac_excluded_overlap: frac_excl(self.bases_excl_overlap),
            frac_excluded_capped: frac_excl(self.bases_excl_capped),
            frac_excluded_total: frac_excl_total,
            frac_bases_at_1x: frac_at(1),
            frac_bases_at_5x: frac_at(5),
            frac_bases_at_10x: frac_at(10),
            frac_bases_at_15x: frac_at(15),
            frac_bases_at_20x: frac_at(20),
            frac_bases_at_25x: frac_at(25),
            frac_bases_at_30x: frac_at(30),
            frac_bases_at_40x: frac_at(40),
            frac_bases_at_50x: frac_at(50),
            frac_bases_at_60x: frac_at(60),
            frac_bases_at_100x: frac_at(100),
            fold_80_base_penalty: fold_penalty(0.20),
            fold_90_base_penalty: fold_penalty(0.10),
            fold_95_base_penalty: fold_penalty(0.05),
        };

        write_tsv(&self.metrics_path, &[metrics])?;

        // ── Coverage histogram ────────────────────────────────────────────────
        let coverage_rows: Vec<WgsCoverageEntry> = self
            .depth_histogram
            .iter()
            .enumerate()
            .map(|(d, &bases)| {
                let frac_bases = bases as f64 / total as f64;
                let frac_bases_at_or_above = at_or_above[d] as f64 / total as f64;
                WgsCoverageEntry {
                    depth: d as u64,
                    bases,
                    frac_bases,
                    bases_at_or_above: at_or_above[d],
                    frac_bases_at_or_above,
                }
            })
            .collect();

        write_tsv(&self.coverage_path, &coverage_rows)?;

        // ── Coverage depth plot ───────────────────────────────────────────────
        // Clip X axis at the 99.5th-percentile depth so the long tail is trimmed.
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "percentile depth ≤ coverage_cap and non-negative"
        )]
        let x_max = percentile_from_hist(&self.depth_histogram, 0.999, total) as u64;
        self.plot_histogram(x_max, mean)?;

        log::info!(
            "wgs: genome_territory={total}, mean_coverage={mean:.2}, \
             metrics={}, coverage={}, plot={}",
            self.metrics_path.display(),
            self.coverage_path.display(),
            self.plot_path.display(),
        );

        Ok(())
    }

    /// Walk the per-position depth array for the current contig and accumulate
    /// per-position statistics into the global histograms. As it goes, it
    /// zeroes each slot in the depth array (so the next contig starts clean
    /// without a separate memset) and accumulates the `bases_excl_capped`
    /// overflow counter from the per-contig buffer.
    fn finalize_contig(&mut self) {
        let Some(ref_id) = self.current_ref_id else {
            return;
        };

        // Global capped-exclusion counter absorbs the per-contig cap overflow.
        self.bases_excl_capped += self.depth.take_excl_capped();

        let non_n_bv = std::mem::take(&mut self.current_non_n);
        let intervals_bv = self.intervals.as_ref().map(|intervals| intervals.contig_bitvec(ref_id));

        let contig_len = self.depth.len();
        let capped = self.coverage_cap;
        let mut eligible: u64 = 0;
        for pos in 0..contig_len {
            let depth = self.depth.take(pos);
            // Skip positions that are N in the reference, or outside the
            // optional intervals mask.
            let non_n_ok = pos < non_n_bv.len() && non_n_bv[pos];
            let interval_ok = intervals_bv.as_ref().is_none_or(|bv| pos < bv.len() && bv[pos]);
            if !non_n_ok || !interval_ok {
                continue;
            }
            eligible += 1;
            let idx = depth.min(capped) as usize;
            self.depth_histogram[idx] += 1;
        }
        self.genome_territory += eligible;

        // All positions in the active contig have been read-and-zeroed above,
        // so the depth array is ready for the next contig without an explicit
        // memset.
    }

    /// Write a PDF area-chart of the coverage depth distribution to [`Self::plot_path`].
    ///
    /// The X axis runs from 0 to `x_max` (the 99.9th-percentile depth), so the long
    /// high-depth tail is trimmed.  The Y axis is the fraction of genome positions at
    /// each depth.  A theoretical Poisson distribution is overlaid as a dashed line
    /// using the observed `mean_coverage` as lambda.  Returns immediately when `x_max`
    /// is 0 or no territory exists.
    fn plot_histogram(&self, x_max: u64, mean_coverage: f64) -> Result<()> {
        if self.genome_territory == 0 || x_max == 0 {
            return Ok(());
        }

        let territory = self.genome_territory as f64;

        let x_max_idx = usize::try_from(x_max).expect("x_max ≤ coverage_cap (u16)");
        let xy: Vec<(f64, f64)> = self
            .depth_histogram
            .iter()
            .enumerate()
            .take(x_max_idx + 1)
            .map(|(d, &c)| (d as f64, c as f64 / territory))
            .collect();

        // Theoretical Poisson distribution scaled by the fraction of covered bases.
        let covered_bases = territory - self.depth_histogram[0] as f64;
        let covered_frac = covered_bases / territory;
        let poisson_xy: Vec<(f64, f64)> =
            (0..=x_max).map(|k| (k as f64, poisson_pmf(k, mean_coverage) * covered_frac)).collect();

        let observed = LinePlot::new()
            .with_data(xy)
            .with_color(FG_BLUE)
            .with_fill()
            .with_fill_opacity(0.3)
            .with_legend("Empirical");

        let theoretical = LinePlot::new()
            .with_data(poisson_xy)
            .with_color(FG_TEAL)
            .with_dashed()
            .with_stroke_width(1.5)
            .with_legend("Theoretical");

        let plots: Vec<Plot> = vec![observed.into(), theoretical.into()];

        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(&self.plot_title)
            .with_x_label("Coverage Depth (X)")
            .with_y_label("Fraction of Genome")
            .with_minor_ticks(5)
            .with_show_minor_grid(true)
            .with_legend_position(LegendPosition::InsideTopRight);

        write_plot_pdf(plots, layout, &self.plot_path)
    }

    /// Count non-N positions in `non_n` that fall within the configured intervals
    /// for `ref_id` (or all positions if no intervals are configured). Used
    /// only by `finish_metrics` for contigs that had zero eligible reads and
    /// therefore never triggered `finalize_contig`.
    fn count_eligible_positions(&self, ref_id: usize, non_n: &BitVec) -> u64 {
        match &self.intervals {
            None => non_n.count_ones() as u64,
            Some(intervals) => {
                let ivl = intervals.contig_bitvec(ref_id);
                if ivl.is_empty() {
                    return 0;
                }
                // AND the two bitvecs over their shared length, then count ones.
                let len = non_n.len().min(ivl.len());
                (non_n[..len].to_bitvec() & &ivl[..len]).count_ones() as u64
            }
        }
    }

    /// Begin processing a new contig: load the reference sequence to build the
    /// non-N bitvec, resize the per-position depth array, and reset the mate
    /// buffer (mates cannot straddle contigs).
    fn begin_contig(&mut self, ref_id: usize) -> Result<()> {
        // `dict` is set in `initialize` before any record arrives; `accept`
        // is the only caller that drives contig transitions.
        let name = self.dict.as_ref().unwrap().get_by_index(ref_id).map_or("", |m| m.name());
        // `contig_length` returns a u64 of genomic coordinates; platforms with
        // < 64-bit pointers are rejected by the compile_error! in `lib.rs`, so
        // the cast is lossless.
        #[allow(clippy::cast_possible_truncation)]
        let contig_len = self.reference.contig_length(name).map_or(0, |n| n as usize);
        let seq = self.reference.load_contig(name, false)?;
        self.current_non_n = build_non_n_bitvec(&seq);
        self.depth.reset_for_contig(contig_len);
        self.mate_buffer.clear();
        self.current_ref_id = Some(ref_id);
        self.processed_contigs.insert(ref_id);
        Ok(())
    }

    /// Walk the CIGAR of `record` and apply its per-base depth contribution,
    /// parameterised on a [`DepthAction`] that encodes the mate-routing case:
    ///
    /// - [`AloneAction`] — no mate overlap; always push depth.
    /// - [`BufferAction`] — buffer for a later mate; push depth and record
    ///   each position in the overlap bitmap.
    /// - [`PairAction`] — second read of a pair; skip positions the cached
    ///   mate already counted (`excl_overlap`), push depth elsewhere.
    ///
    /// The `DepthAction` methods for unused hooks are zero-cost default
    /// implementations; monomorphisation compiles each call site into
    /// straight-line per-action code with no branches on the action kind.
    #[inline]
    fn walk_depth<A: DepthAction>(&mut self, record: &RikerRecord, mut action: A) {
        let Some(pos_1based) = record.alignment_start() else {
            return;
        };
        let ref_start_u64 = pos_1based.get() as u64 - 1;
        let quals: &[u8] = record.quality_scores();
        let contig_len_u64 = self.depth.len() as u64;
        let min_bq = self.min_bq;

        for (ref_off, read_off, len) in iter_aligned_blocks(record) {
            let block_ref_start = ref_start_u64 + u64::from(ref_off);
            // CIGAR blocks that fall entirely past the contig end — a
            // malformed-BAM case the spec disallows — are silently dropped.
            if block_ref_start >= contig_len_u64 {
                continue;
            }
            // Clamp the block's effective length to the remaining contig
            // territory. Any tail past the contig end is silently dropped.
            let block_ref_end = (block_ref_start + u64::from(len)).min(contig_len_u64);
            // `usable` ≤ `len` ≤ `u32::MAX` by construction (CIGAR op length
            // is u32 in noodles), so the cast is lossless.
            #[allow(clippy::cast_possible_truncation)]
            let usable = (block_ref_end - block_ref_start) as u32;

            // Clamp both bounds to quals.len(): on a malformed BAM where the
            // CIGAR claims more read bases than the qual array provides,
            // `read_off` itself can exceed `quals.len()`, which would make the
            // naive slice panic.
            let read_start = (read_off as usize).min(quals.len());
            let read_end = (read_start + usable as usize).min(quals.len());
            let available = &quals[read_start..read_end];
            #[allow(clippy::cast_possible_truncation)]
            let base_pos = block_ref_start as u32;
            for (i, &q) in available.iter().enumerate() {
                if q < min_bq {
                    self.bases_excl_baseq += 1;
                } else {
                    #[allow(clippy::cast_possible_truncation)]
                    let ref_pos = base_pos + i as u32;
                    if action.is_mate_covered(ref_pos) {
                        self.bases_excl_overlap += 1;
                    } else {
                        self.depth.push(ref_pos);
                        action.on_depth_counted(ref_pos);
                    }
                }
            }
            // If the qual array is shorter than the CIGAR claims the block
            // spans (malformed BAM — the spec says they match), the missing
            // quals can't have passed BQ, so we fold them into the BQ-fail
            // counter to keep the exclusion totals accurate.
            self.bases_excl_baseq += (usable as usize - available.len()) as u64;
        }
    }
}

// ─── Collector trait impl ─────────────────────────────────────────────────────

impl Collector for WgsCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        self.reference.validate_bam_header(header)?;

        self.dict = Some(SequenceDictionary::from(header));

        self.sample = derive_sample(&self.input_path, header);
        self.plot_title = format!("Coverage Depth Distribution of {}", self.sample);
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // Skip reads that should never enter the pileup.
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }

        // QC-fail: skip silently (not counted in exclusions, matching Picard).
        if flags.is_qc_fail() {
            return Ok(());
        }

        // Mapping quality filter (checked first among the counted exclusions,
        // matching Picard's filter ordering so per-category fractions align).
        let mapq = record.mapping_quality().map_or(0, |m| m.get());
        if mapq < self.min_mapq {
            self.bases_excl_mapq += count_aligned_bases(record);
            return Ok(());
        }

        // Duplicate filter
        if !self.include_duplicates && flags.is_duplicate() {
            self.bases_excl_dupe += count_aligned_bases(record);
            return Ok(());
        }

        // Unpaired filter: must be segmented with a mapped mate.
        if !self.include_unpaired {
            let is_paired_mapped = flags.is_segmented() && !flags.is_mate_unmapped();
            if !is_paired_mapped {
                self.bases_excl_unpaired += count_aligned_bases(record);
                return Ok(());
            }
        }

        // ── Read passed all whole-read filters ───────────────────────────────
        let Some(ref_id) = record.reference_sequence_id() else {
            return Ok(());
        };

        // Contig transition: finalize the previous contig and start a fresh
        // depth array sized to the new contig. Any buffered mates belonged to
        // the old contig (mates can never straddle contigs) and were already
        // pushed when they arrived, so we just drop the buffer.
        if Some(ref_id) != self.current_ref_id {
            if self.current_ref_id.is_some() {
                self.finalize_contig();
            }
            self.begin_contig(ref_id)?;
        }

        // Pre-extract the fields once; the same values feed both `probe_fields`
        // (which consumes the struct) and the WouldBuffer follow-up `insert_fields`,
        // so we don't have to re-read them from `record` after the probe.
        // `name()` returns `&BStr`; deref to raw bytes — the mate buffer hashes on
        // bytes, not on string semantics.
        let ref_id = record.reference_sequence_id();
        let name: Option<&[u8]> = record.name().map(|n| &**n);
        let mate_alignment_start = record.mate_alignment_start();
        let probe_fields = MateProbeFields {
            flags,
            ref_id,
            mate_ref_id: record.mate_reference_sequence_id(),
            name,
            alignment_start: record.alignment_start(),
            alignment_end: record.alignment_end(),
            mate_alignment_start,
        };
        match self.mate_buffer.probe_fields(probe_fields) {
            Peek::Alone => self.walk_depth(record, AloneAction),
            Peek::WouldBuffer { overlap_start, overlap_len } => {
                let mut bitmap = bitvec![u64, Lsb0; 0; overlap_len as usize];
                self.walk_depth(record, BufferAction { overlap_start, bitmap: &mut bitmap });
                // Contract of `probe_fields` mirrors `probe`: on WouldBuffer
                // the ref_id, name, and mate_alignment_start are all present.
                let ref_id = ref_id.expect("probe_fields WouldBuffer guarantees ref_id is Some");
                let name = name.expect("probe_fields WouldBuffer guarantees name is Some");
                let mate_pos = mate_alignment_start.map_or(0, |p| p.get());
                self.mate_buffer.insert_fields(
                    ref_id,
                    name,
                    mate_pos,
                    CachedMate { overlap_start, bitmap },
                );
            }
            Peek::PairWith(cached) => {
                self.walk_depth(record, PairAction { cached: &cached });
            }
        }

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        // All buffered mates have already pushed their depth contributions, so
        // finalization just needs to close out the last contig and write out.
        self.finish_metrics()
    }

    fn name(&self) -> &'static str {
        "wgs"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        RikerRecordRequirements::NONE
    }
}

// ─── Metric structs ───────────────────────────────────────────────────────────

/// Whole-genome sequencing coverage summary metrics.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct WgsMetrics {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Number of non-N reference bases in the genome territory (or interval set).
    pub genome_territory: u64,
    /// Mean coverage depth across all genome territory positions.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_coverage: f64,
    /// Standard deviation of coverage depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub sd_coverage: f64,
    /// Median coverage depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub median_coverage: f64,
    /// Median absolute deviation of coverage depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mad_coverage: f64,
    /// Fraction of bases excluded due to low mapping quality.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_mapq: f64,
    /// Fraction of bases excluded because the read was marked as a duplicate.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_dupe: f64,
    /// Fraction of bases excluded because the read was unpaired or had an unmapped mate.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_unpaired: f64,
    /// Fraction of bases excluded due to low base quality.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_baseq: f64,
    /// Fraction of bases excluded due to read-pair overlap (second occurrence of a read name).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_overlap: f64,
    /// Fraction of bases excluded because the coverage cap was exceeded.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_capped: f64,
    /// Total fraction of bases excluded for any reason.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_total: f64,
    /// Fraction of genome territory positions with coverage ≥ 1×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_1x: f64,
    /// Fraction of genome territory positions with coverage ≥ 5×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_5x: f64,
    /// Fraction of genome territory positions with coverage ≥ 10×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_10x: f64,
    /// Fraction of genome territory positions with coverage ≥ 15×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_15x: f64,
    /// Fraction of genome territory positions with coverage ≥ 20×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_20x: f64,
    /// Fraction of genome territory positions with coverage ≥ 25×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_25x: f64,
    /// Fraction of genome territory positions with coverage ≥ 30×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_30x: f64,
    /// Fraction of genome territory positions with coverage ≥ 40×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_40x: f64,
    /// Fraction of genome territory positions with coverage ≥ 50×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_50x: f64,
    /// Fraction of genome territory positions with coverage ≥ 60×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_60x: f64,
    /// Fraction of genome territory positions with coverage ≥ 100×.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_100x: f64,
    /// Mean coverage divided by the 20th-percentile coverage (80% of bases are at/above that level).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub fold_80_base_penalty: f64,
    /// Mean coverage divided by the 10th-percentile coverage.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub fold_90_base_penalty: f64,
    /// Mean coverage divided by the 5th-percentile coverage.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub fold_95_base_penalty: f64,
}

/// One row of the per-depth coverage histogram.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct WgsCoverageEntry {
    /// Coverage depth (number of reads covering a position).
    pub depth: u64,
    /// Number of reference positions with exactly this depth.
    pub bases: u64,
    /// Fraction of genome territory with exactly this depth.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases: f64,
    /// Number of reference positions with depth ≥ this depth.
    pub bases_at_or_above: u64,
    /// Fraction of genome territory with depth ≥ this depth.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_at_or_above: f64,
}

// ─── ContigDepth ─────────────────────────────────────────────────────────────

/// Per-position depth array for a single contig. Allocated once and reused
/// across contigs: `reset_for_contig(n)` clears and resizes to `n` u16 slots,
/// zeroing all entries. The caller drains the array by `take(pos)` on each
/// position during `finalize_contig`, which returns the current depth and
/// zeros the slot — after the finalize pass the array is already in the
/// "fresh" state for the next contig (no separate memset needed).
struct ContigDepth {
    depth: Vec<u16>,
    coverage_cap: u16,
    /// Running count of bases that arrived at a position already at the cap.
    /// Finalize drains this into the collector-level accumulator.
    excl_capped: u64,
}

impl ContigDepth {
    fn new(coverage_cap: u16) -> Self {
        Self { depth: Vec::new(), coverage_cap, excl_capped: 0 }
    }

    /// Length of the currently-active contig's depth array.
    fn len(&self) -> usize {
        self.depth.len()
    }

    /// Resize for a new contig. `Vec::resize(n, 0)` reuses the existing backing
    /// allocation up to current capacity, truncates down or zero-fills new
    /// slots as needed; the caller is expected to have drained the previous
    /// contig via `take(pos)` so the existing slots are already zero.
    fn reset_for_contig(&mut self, contig_len: usize) {
        self.depth.resize(contig_len, 0);
    }

    /// Increment the depth at `ref_pos`, saturating at `coverage_cap`. Attempts
    /// to push past the cap bump the capped-exclusion counter instead.
    fn push(&mut self, ref_pos: u32) {
        let pos = ref_pos as usize;
        // Defensive: `walk_depth` already truncates each block to the contig
        // length, so this branch should never hit.
        if pos >= self.depth.len() {
            return;
        }
        let d = &mut self.depth[pos];
        if *d < self.coverage_cap {
            *d += 1;
        } else {
            self.excl_capped += 1;
        }
    }

    /// Read the depth at `pos` and zero the slot in one step (used by
    /// `finalize_contig` as it walks the array).
    fn take(&mut self, pos: usize) -> u16 {
        std::mem::take(&mut self.depth[pos])
    }

    /// Drain and return the running capped-exclusion counter.
    fn take_excl_capped(&mut self) -> u64 {
        std::mem::take(&mut self.excl_capped)
    }
}

// ─── Mate buffering for read-pair overlap detection ──────────────────────────

/// Per-mate state retained while the overlapping mate of a buffered read has
/// not yet arrived. The only question we need to answer at pair time is:
/// "did the first read contribute BQ-passing depth at this reference
/// position?" — a single bit per position in the overlap region. We record
/// exactly that bitmap, sized to the overlap span between the two mates'
/// alignments, rather than keeping quals + CIGAR blocks to reconstruct the
/// answer at query time.
struct CachedMate {
    /// 0-based reference position of the first bit; equals the mate's
    /// expected `alignment_start`.
    overlap_start: u32,
    /// Bit `i` (LSB-first) is set when the first read contributed
    /// BQ-passing depth at reference position `overlap_start + i`. The
    /// length of the bitmap (in bits) is the overlap region's span.
    bitmap: BitVec<u64, Lsb0>,
}

impl CachedMate {
    /// Returns `true` if the cached read contributed BQ-passing depth at
    /// `ref_pos`. Positions outside the overlap region return `false`.
    fn covered_at(&self, ref_pos: u32) -> bool {
        if ref_pos < self.overlap_start {
            return false;
        }
        let idx = (ref_pos - self.overlap_start) as usize;
        if idx >= self.bitmap.len() {
            return false;
        }
        self.bitmap[idx]
    }
}

/// Per-record strategy for the inner depth-push loop, capturing the two
/// orthogonal decisions the three mate-routing cases disagree on:
///
/// 1. [`is_mate_covered`] — whether the paired-mate already contributed
///    depth at this reference position; if so, the current record counts
///    the base as an `excl_overlap` exclusion instead of pushing depth.
/// 2. [`on_depth_counted`] — whether to record the position so that a
///    later-arriving mate can dedupe against it (populates the bitmap in
///    [`CachedMate`]).
///
/// [`AloneAction`] uses the defaults (never-covered, never-record).
/// [`BufferAction`] overrides `on_depth_counted` to populate the overlap
/// bitmap. [`PairAction`] overrides `is_mate_covered` to consult the
/// cached mate. Monomorphization keeps the unused hook in each case at
/// zero cost (a `false` return or empty body).
///
/// [`is_mate_covered`]: DepthAction::is_mate_covered
/// [`on_depth_counted`]: DepthAction::on_depth_counted
trait DepthAction {
    /// Return `true` if the other read of this pair already counted
    /// `ref_pos`. The caller should not push depth again.
    #[inline]
    fn is_mate_covered(&self, _ref_pos: u32) -> bool {
        false
    }

    /// Record that depth was pushed at `ref_pos`, so a later-arriving mate
    /// can dedupe against the current record's contribution.
    #[inline]
    fn on_depth_counted(&mut self, _ref_pos: u32) {}
}

/// No-overlap-awareness action for reads with no overlapping mate.
struct AloneAction;
impl DepthAction for AloneAction {}

/// Action for reads that will be buffered: records each depth-push
/// position into `bitmap` so the arriving mate can skip the double-count.
struct BufferAction<'a> {
    overlap_start: u32,
    bitmap: &'a mut BitVec<u64, Lsb0>,
}

impl DepthAction for BufferAction<'_> {
    #[inline]
    fn on_depth_counted(&mut self, ref_pos: u32) {
        if ref_pos >= self.overlap_start {
            let idx = (ref_pos - self.overlap_start) as usize;
            if idx < self.bitmap.len() {
                self.bitmap.set(idx, true);
            }
        }
    }
}

/// Action for the second read of a pair: consults the cached mate's
/// bitmap so overlap positions tally as `excl_overlap` instead of
/// double-pushing depth.
struct PairAction<'a> {
    cached: &'a CachedMate,
}

impl DepthAction for PairAction<'_> {
    #[inline]
    fn is_mate_covered(&self, ref_pos: u32) -> bool {
        self.cached.covered_at(ref_pos)
    }
}

/// Iterate `(ref_offset, read_offset, len)` tuples for each M/=/X CIGAR
/// operation in `record`. `ref_offset` is relative to the record's
/// `alignment_start` and `read_offset` indexes into the record's quality /
/// sequence arrays. Non-M CIGAR ops advance the ref/read cursors but do not
/// yield.
fn iter_aligned_blocks(record: &RikerRecord) -> impl Iterator<Item = (u32, u32, u32)> + '_ {
    let mut ref_off: u32 = 0;
    let mut read_off: u32 = 0;
    record.cigar_ops().filter_map(move |op| {
        // CIGAR op lengths are stored as u32 in BAM/SAM, so the cast is lossless.
        #[allow(clippy::cast_possible_truncation)]
        let len = op.len() as u32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let out = (ref_off, read_off, len);
                ref_off += len;
                read_off += len;
                Some(out)
            }
            Kind::Insertion | Kind::SoftClip => {
                read_off += len;
                None
            }
            Kind::Deletion | Kind::Skip => {
                ref_off += len;
                None
            }
            Kind::HardClip | Kind::Pad => None,
        }
    })
}

// ─── Poisson helpers ─────────────────────────────────────────────────────────

/// Compute the Poisson PMF P(X = k) for a given lambda, using log-space
/// arithmetic to avoid overflow for large k.
fn poisson_pmf(k: u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    let log_pmf = k as f64 * lambda.ln() - lambda - ln_factorial(k);
    log_pmf.exp()
}

/// Compute ln(k!) by iterative summation of ln(i) for i in 2..=k.
fn ln_factorial(k: u64) -> f64 {
    (2..=k).map(|i| (i as f64).ln()).sum()
}

// ─── CIGAR helper ─────────────────────────────────────────────────────────────

/// Count the number of reference-aligned bases (M/=/X) in a record's CIGAR.
fn count_aligned_bases(record: &RikerRecord) -> u64 {
    let mut count: u64 = 0;
    for op in record.cigar_ops() {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                count += op.len() as u64;
            }
            _ => {}
        }
    }
    count
}

// ─── Bitvec helper ────────────────────────────────────────────────────────────

/// Build a `BitVec` from a reference sequence where each bit is `true` if the
/// corresponding base is not N/n. This replaces storing the full sequence,
/// reducing per-contig memory from ~250MB to ~31MB for the largest chromosomes
/// and enabling hardware-accelerated `count_ones()` for eligible-position counts.
///
/// Builds the bitvec one word at a time (64 bases → 1 `usize`) to avoid the
/// overhead of per-bit indexing through bitvec's API.
fn build_non_n_bitvec(seq: &[u8]) -> BitVec {
    const W: usize = usize::BITS as usize;
    let len = seq.len();
    let mut words = Vec::with_capacity(len.div_ceil(W));

    for chunk in seq.chunks(W) {
        let mut word: usize = 0;
        for (j, &b) in chunk.iter().enumerate() {
            if b != b'N' && b != b'n' {
                word |= 1 << j;
            }
        }
        words.push(word);
    }

    let mut bv = BitVec::from_vec(words);
    bv.truncate(len);
    bv
}

// ─── Percentile helper ────────────────────────────────────────────────────────

/// Return the depth at which `lower_frac` of positions have depth ≤ that value
/// (i.e., `1 - lower_frac` of positions are at or above).
///
/// For `lower_frac = 0.20`, this is the 20th percentile depth (80% of positions
/// are at or above), used for `fold_80_base_penalty`.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "depth indices and totals fit safely in u64 conversions"
)]
fn percentile_from_hist(hist: &[u64], lower_frac: f64, total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let threshold = (total as f64 * lower_frac).ceil() as u64;
    let mut running: u64 = 0;
    for (d, &count) in hist.iter().enumerate() {
        running += count;
        if running >= threshold {
            return d as f64;
        }
    }
    (hist.len() - 1) as f64
}

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Percentile tests ─────────────────────────────────────────────────────

    #[test]
    fn test_percentile_from_hist_basic() {
        // hist[0]=1, hist[1]=1, hist[2]=1 → total=3
        // 20th pct: threshold=ceil(0.2*3)=1, running after d=0 is 1 ≥ 1 → depth=0
        let hist = [1u64, 1, 1];
        assert!((percentile_from_hist(&hist, 0.20, 3) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_from_hist_80pct() {
        // 10 positions at depth=10; fold_80 threshold = ceil(0.2*10) = 2
        // All at depth=10, so after depth=10 running=10 ≥ 2 → depth=10
        let mut hist = vec![0u64; 11];
        hist[10] = 10;
        assert!((percentile_from_hist(&hist, 0.20, 10) - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_empty_histogram() {
        assert!((percentile_from_hist(&[], 0.5, 0) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_single_bin() {
        // All 10 positions at depth 0
        let hist = [10u64];
        assert!((percentile_from_hist(&hist, 0.5, 10) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_frac_zero() {
        // lower_frac=0.0 → threshold=ceil(0)=0, running at d=0 is 5 ≥ 0 → depth=0
        let mut hist = vec![0u64; 11];
        hist[0] = 5;
        hist[10] = 5;
        assert!((percentile_from_hist(&hist, 0.0, 10) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_frac_one() {
        // lower_frac=1.0 → threshold=ceil(10)=10
        // hist[5]=5, hist[10]=5 → running at d=10 is 10 ≥ 10 → depth=10
        let mut hist = vec![0u64; 11];
        hist[5] = 5;
        hist[10] = 5;
        assert!((percentile_from_hist(&hist, 1.0, 10) - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_percentile_two_bins() {
        // hist[5]=5, hist[10]=5, total=10, lower_frac=0.5
        // threshold=ceil(5)=5, running at d=5 is 5 ≥ 5 → depth=5
        let mut hist = vec![0u64; 11];
        hist[5] = 5;
        hist[10] = 5;
        assert!((percentile_from_hist(&hist, 0.5, 10) - 5.0).abs() < f64::EPSILON);
    }

    // ── ContigDepth tests ────────────────────────────────────────────────────

    #[test]
    fn test_contig_depth_basic_push_and_take() {
        let mut d = ContigDepth::new(250);
        d.reset_for_contig(10);
        d.push(3);
        d.push(3);
        d.push(7);

        assert_eq!(d.take(0), 0);
        assert_eq!(d.take(3), 2);
        assert_eq!(d.take(7), 1);
        // `take` zeroes the slot:
        assert_eq!(d.take(3), 0);
    }

    #[test]
    fn test_contig_depth_saturates_at_cap() {
        let mut d = ContigDepth::new(3);
        d.reset_for_contig(5);
        for _ in 0..10 {
            d.push(2);
        }
        assert_eq!(d.take(2), 3);
        // 10 pushes, cap=3 → 7 counted as capped.
        assert_eq!(d.take_excl_capped(), 7);
    }

    #[test]
    fn test_contig_depth_resize_between_contigs() {
        let mut d = ContigDepth::new(250);
        d.reset_for_contig(10);
        d.push(5);
        d.take(5);
        // Switch to a larger contig; existing backing store grows via resize.
        d.reset_for_contig(100);
        assert_eq!(d.len(), 100);
        d.push(99);
        assert_eq!(d.take(99), 1);
    }

    #[test]
    fn test_contig_depth_ignores_out_of_range() {
        let mut d = ContigDepth::new(250);
        d.reset_for_contig(10);
        d.push(42); // out of range — silently ignored
        assert_eq!(d.take_excl_capped(), 0);
        assert_eq!(d.take(0), 0);
    }

    // ── CachedMate tests ─────────────────────────────────────────────────────

    /// Build a CachedMate by explicitly naming the ref positions where the
    /// first read contributed BQ-passing depth.
    fn cached_mate(overlap_start: u32, overlap_len: u32, covered: &[u32]) -> CachedMate {
        let mut bitmap = bitvec![u64, Lsb0; 0; overlap_len as usize];
        for &ref_pos in covered {
            assert!(ref_pos >= overlap_start);
            let idx = (ref_pos - overlap_start) as usize;
            assert!(idx < overlap_len as usize);
            bitmap.set(idx, true);
        }
        CachedMate { overlap_start, bitmap }
    }

    #[test]
    fn test_cached_mate_covered_at_contiguous_region() {
        // First read contributed depth at every position in [100, 109].
        let cached = cached_mate(100, 10, &[100, 101, 102, 103, 104, 105, 106, 107, 108, 109]);
        assert!(cached.covered_at(100));
        assert!(cached.covered_at(105));
        assert!(cached.covered_at(109));
        // Outside overlap region:
        assert!(!cached.covered_at(99));
        assert!(!cached.covered_at(110));
    }

    #[test]
    fn test_cached_mate_covered_at_with_gap() {
        // Simulates CIGAR "5M2D5M" — BQ-passing depth at [100,104] and
        // [107,111] (a 2bp gap in between from the deletion). The bitmap
        // spans [100, 111] so gap positions are represented as 0 bits.
        let cached = cached_mate(100, 12, &[100, 101, 102, 103, 104, 107, 108, 109, 110, 111]);
        assert!(cached.covered_at(100));
        assert!(cached.covered_at(104));
        // Inside the 2bp deletion — not covered:
        assert!(!cached.covered_at(105));
        assert!(!cached.covered_at(106));
        // After the deletion:
        assert!(cached.covered_at(107));
        assert!(cached.covered_at(111));
    }

    #[test]
    fn test_cached_mate_covered_at_with_bq_fail_gaps() {
        // First read aligned across [100, 109] but had a BQ-failing base at
        // 103, so the bitmap has a 0 bit there even though the position is
        // within the overlap region.
        let cached = cached_mate(100, 10, &[100, 101, 102, 104, 105, 106, 107, 108, 109]);
        assert!(cached.covered_at(102));
        assert!(!cached.covered_at(103)); // BQ-fail gap
        assert!(cached.covered_at(104));
    }

    #[test]
    fn test_cached_mate_covered_at_out_of_range() {
        let cached = cached_mate(100, 10, &[100, 109]);
        assert!(!cached.covered_at(50));
        assert!(!cached.covered_at(99));
        assert!(!cached.covered_at(110));
        assert!(!cached.covered_at(u32::MAX));
    }
}
