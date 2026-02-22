use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};

use anyhow::Result;
use bitvec::prelude::*;
use clap::Args;
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::Collector;
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use crate::counter::Counter;
use crate::fasta::Fasta;
use crate::intervals::Intervals;
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, write_tsv};
use crate::plotting::{FG_BLUE, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_plot_pdf};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::derive_sample;
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
    /// At each position, the internal buffer stores up to 2x coverage-cap
    /// read-name hashes for read-pair overlap detection before dedup. Bases
    /// exceeding the cap are counted as capped exclusions. Memory usage is
    /// approximately max-read-ref-span x 2 x coverage-cap x 8 bytes.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_COVERAGE_CAP)]
    pub coverage_cap: u32,

    /// Maximum reference span (in bases) of a single read alignment.
    ///
    /// The internal circular buffer holds this many genomic positions. It must
    /// be at least as large as the longest CIGAR reference span. The default
    /// of 1000 suits short-read Illumina data; increase for long-read
    /// technologies. Memory scales linearly with this value.
    #[arg(long, default_value_t = WgsOptions::DEFAULT_MAX_READ_REF_SPAN)]
    pub max_read_ref_span: usize,
}

impl WgsOptions {
    const DEFAULT_EXCLUDE_UNPAIRED: bool = true;
    const DEFAULT_MIN_MAPQ: u8 = 20;
    const DEFAULT_MIN_BQ: u8 = 20;
    const DEFAULT_COVERAGE_CAP: u32 = 250;
    const DEFAULT_MAX_READ_REF_SPAN: usize = 1000;
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
            max_read_ref_span: Self::DEFAULT_MAX_READ_REF_SPAN,
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
        let (mut reader, header) =
            AlignmentReader::new(&self.input.input, Some(&self.reference.reference))?;
        let reference = Fasta::from_path(&self.reference.reference)?;

        let mut collector =
            WgsCollector::new(&self.input.input, &self.output.output, reference, &self.options)?;

        collector.initialize(&header)?;
        let mut progress = ProgressLogger::new("wgs", "reads", 5_000_000);
        for result in reader.record_bufs(&header) {
            let record = result?;
            progress.record_with(&record, &header);
            collector.accept(&record, &header)?;
        }
        progress.finish();
        collector.finish()
    }
}

// ─── Collector ────────────────────────────────────────────────────────────────

/// Accumulates per-position coverage statistics using a circular depth buffer.
pub struct WgsCollector {
    // Output paths
    metrics_path: PathBuf,
    coverage_path: PathBuf,
    plot_path: PathBuf,
    plot_title: String,

    // Input path (needed for sample name derivation in initialize)
    input_path: PathBuf,

    // Circular depth buffer (replaces PileupEngine)
    buffer: DepthBuffer,

    // Config
    reference: Fasta,
    intervals: Option<Intervals>,
    include_duplicates: bool,
    include_unpaired: bool,
    min_mapq: u8,
    min_bq: u8,
    coverage_cap: u32,

    // BAM contig metadata
    dict: Option<SequenceDictionary>,

    // The ref_id of the reads currently being added to the buffer.
    buffer_ref_id: Option<usize>,

    // Per-contig working state
    current_ref_id: Option<usize>,
    current_non_n: BitVec,
    current_covered: u64,
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
    bases_excl_truncated: u64,

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
            buffer: DepthBuffer::new(options.max_read_ref_span, options.coverage_cap),
            reference,
            intervals,
            include_duplicates: options.include_duplicates,
            include_unpaired: !options.exclude_unpaired_reads,
            min_mapq: options.min_mapq,
            min_bq: options.min_bq,
            coverage_cap: options.coverage_cap,
            dict: None,
            buffer_ref_id: None,
            current_ref_id: None,
            current_non_n: BitVec::EMPTY,
            current_covered: 0,
            processed_contigs: HashSet::new(),
            genome_territory: 0,
            depth_histogram: vec![0u64; hist_len],
            bases_excl_mapq: 0,
            bases_excl_dupe: 0,
            bases_excl_unpaired: 0,
            bases_excl_baseq: 0,
            bases_excl_overlap: 0,
            bases_excl_capped: 0,
            bases_excl_truncated: 0,
            sample: String::new(),
        })
    }

    /// Flush all positions in the buffer before `genome_pos`, recording
    /// coverage statistics for each.
    fn flush_before(&mut self, genome_pos: u64) -> Result<()> {
        while self.buffer.len > 0 && self.buffer.genome_pos < genome_pos {
            self.flush_one()?;
        }
        Ok(())
    }

    /// Flush all remaining positions in the buffer.
    fn flush_all(&mut self) -> Result<()> {
        while self.buffer.len > 0 {
            self.flush_one()?;
        }
        Ok(())
    }

    /// Flush the head position of the buffer and record its coverage.
    fn flush_one(&mut self) -> Result<()> {
        let genome_pos = self.buffer.genome_pos;
        let ref_id = self.buffer_ref_id.expect("buffer_ref_id set when buffer is non-empty");
        let result = self.buffer.flush_position(self.coverage_cap);
        self.record_position(ref_id, genome_pos, &result)
    }

    /// Process a flushed position: handle contig transitions, N-base checks,
    /// interval checks, and update the histogram and accumulators.
    #[allow(clippy::cast_possible_truncation)] // compile_error! in lib.rs guarantees 64-bit platform
    fn record_position(&mut self, ref_id: usize, pos: u64, result: &FlushResult) -> Result<()> {
        // ── Contig transition ────────────────────────────────────────────────
        if Some(ref_id) != self.current_ref_id {
            self.finalize_contig();
            let name = self.dict.as_ref().unwrap().get_by_index(ref_id).map_or("", |m| m.name());
            self.current_non_n = build_non_n_bitvec(&self.reference.load_contig(name, false)?);
            self.current_ref_id = Some(ref_id);
            self.current_covered = 0;
            self.processed_contigs.insert(ref_id);
        }

        // ── Skip N bases ─────────────────────────────────────────────────────
        let pos_idx = pos as usize;
        if pos_idx >= self.current_non_n.len() || !self.current_non_n[pos_idx] {
            return Ok(());
        }

        // ── Skip positions outside intervals ─────────────────────────────────
        if !self.pos_in_intervals(ref_id, pos) {
            return Ok(());
        }

        self.current_covered += 1;

        // ── Update histogram and accumulators ────────────────────────────────
        let capped = result.depth.min(u64::from(self.coverage_cap));
        let capped_idx = capped as usize;
        self.depth_histogram[capped_idx] += 1;
        self.bases_excl_capped += result.excl_capped;
        self.bases_excl_overlap += result.excl_overlap;

        Ok(())
    }

    /// Finalize the current contig: count all eligible (non-N, in-interval) positions
    /// and add uncovered positions to `depth_histogram[0]`.
    fn finalize_contig(&mut self) {
        let Some(ref_id) = self.current_ref_id else {
            return;
        };

        let non_n_bv = std::mem::take(&mut self.current_non_n);
        let non_n = self.count_eligible_positions(ref_id, &non_n_bv);
        self.genome_territory += non_n;

        // Positions with depth=0 are the eligible positions not covered by the pileup.
        let uncovered = non_n.saturating_sub(self.current_covered);
        self.depth_histogram[0] += uncovered;

        self.current_covered = 0;
    }

    /// Count non-N positions in `non_n` that fall within the configured intervals
    /// for `ref_id` (or all positions if no intervals are configured).
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

    /// Return true if `pos` is within the configured intervals for `ref_id`,
    /// or if no intervals are configured (all positions pass).
    #[expect(clippy::cast_possible_truncation, reason = "genomic coordinates fit in u32")]
    fn pos_in_intervals(&self, ref_id: usize, pos: u64) -> bool {
        match &self.intervals {
            None => true,
            Some(intervals) => intervals.contains_pos(ref_id, pos as u32),
        }
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

        let x_max_idx = usize::try_from(x_max).expect("x_max ≤ coverage_cap (u32)");
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
                + self.bases_excl_capped
                + self.bases_excl_truncated,
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
            frac_excluded_truncated: frac_excl(self.bases_excl_truncated),
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

    #[expect(
        clippy::cast_possible_truncation,
        reason = "offset is bounded by num_positions (usize) before use"
    )]
    fn accept(&mut self, record: &RecordBuf, _header: &Header) -> Result<()> {
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
        let Some(pos_1based) = record.alignment_start() else {
            return Ok(());
        };
        let start_pos = pos_1based.get() as u64 - 1;

        // Compute name hash.
        let name_hash = hash_name(record.name().map_or(&[] as &[u8], |n| {
            let s: &[u8] = n;
            s
        }));

        // ── Contig transition: flush all buffered positions ──────────────────
        if self.buffer_ref_id.is_some() && self.buffer_ref_id != Some(ref_id) {
            self.flush_all()?;
        }
        self.buffer_ref_id = Some(ref_id);

        // ── Flush positions before this read's start (they're complete) ──────
        if self.buffer.len > 0 {
            self.flush_before(start_pos)?;
        }

        // ── Initialize buffer range if empty ─────────────────────────────────
        self.buffer.ensure_range(start_pos);

        // ── Walk CIGAR and insert name hashes for aligned bases ──────────────
        let qual: &[u8] = record.quality_scores().as_ref();
        let num_positions = self.buffer.num_positions;
        let min_bq = self.min_bq;

        let mut ref_pos = start_pos;
        let mut read_cursor: usize = 0;

        for op in CigarTrait::iter(record.cigar()).filter_map(Result::ok) {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let first_offset = (ref_pos - self.buffer.genome_pos) as usize;

                    // Determine how many bases fit within the buffer; the rest are truncated.
                    let usable = if first_offset >= num_positions {
                        self.bases_excl_truncated += len as u64;
                        ref_pos += len as u64;
                        read_cursor += len;
                        continue;
                    } else {
                        len.min(num_positions - first_offset)
                    };
                    let truncated = len - usable;
                    self.bases_excl_truncated += truncated as u64;

                    // Extend buffer len once for the whole block.
                    let last_offset = first_offset + usable - 1;
                    if last_offset >= self.buffer.len {
                        self.buffer.len = last_offset + 1;
                    }

                    // Inner loop: only BQ check + push.
                    for i in 0..usable {
                        let q = qual.get(read_cursor + i).copied().unwrap_or(0);
                        if q < min_bq {
                            self.bases_excl_baseq += 1;
                        } else {
                            self.buffer.push(first_offset + i, name_hash);
                        }
                    }
                    ref_pos += len as u64;
                    read_cursor += len;
                }
                Kind::Insertion | Kind::SoftClip => {
                    read_cursor += len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as u64;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        // Flush all remaining buffer positions.
        self.flush_all()?;
        self.finish_metrics()
    }

    fn name(&self) -> &'static str {
        "wgs"
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
    /// Fraction of bases excluded because the read's reference span exceeded --max-read-ref-len.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_truncated: f64,
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

// ─── DepthBuffer ──────────────────────────────────────────────────────────────

/// Circular buffer that accumulates name hashes per genomic position for
/// coverage depth computation.  Replaces the general-purpose `PileupEngine`
/// with a structure that stores only the data WGS needs (name hashes) and
/// performs base-quality filtering at insertion time.
///
/// Both `num_positions` and `max_depth_per_pos` are rounded up to the next
/// power of two so that the hot-path modulo and multiply can be replaced with
/// bitmask and shift operations.
struct DepthBuffer {
    /// Flat storage: position `p` uses slots `[p << depth_shift .. (p << depth_shift) + counts[p])`.
    hashes: Vec<u64>,
    /// Number of name hashes stored at each buffer position.
    counts: Vec<u32>,
    /// Number of bases that arrived at a position after it was full (>= `max_depth_per_pos`).
    excess: Vec<u64>,
    /// Number of genomic positions the buffer can hold (power-of-two, rounded up from CLI value).
    num_positions: usize,
    /// Bitmask for ring index: `num_positions - 1`.
    pos_mask: usize,
    /// Maximum name hash entries per position (power-of-two, rounded up from 2 × `coverage_cap`).
    max_depth_per_pos: usize,
    /// Shift amount equivalent to `× max_depth_per_pos`: `max_depth_per_pos.trailing_zeros()`.
    depth_shift: u32,
    /// Index of the leftmost active position in the circular buffer.
    head: usize,
    /// Genomic position corresponding to `head`.
    genome_pos: u64,
    /// Number of active positions currently in use (`head..head+len`).
    len: usize,
}

impl DepthBuffer {
    /// Create a new depth buffer.
    ///
    /// `num_positions` is the maximum reference span of a single read alignment.
    /// `coverage_cap` is the maximum reported depth; `max_depth_per_pos` is set to
    /// `2 × coverage_cap` to accommodate read-pair overlap before deduplication.
    ///
    /// Both values are rounded up to the next power of two internally so that
    /// ring-index and offset computations use bitmask/shift instead of modulo/multiply.
    fn new(num_positions: usize, coverage_cap: u32) -> Self {
        let num_positions = num_positions.next_power_of_two();
        let pos_mask = num_positions - 1;
        let max_depth_per_pos = (2 * coverage_cap as usize).next_power_of_two();
        let depth_shift = max_depth_per_pos.trailing_zeros();
        Self {
            hashes: vec![0u64; num_positions * max_depth_per_pos],
            counts: vec![0u32; num_positions],
            excess: vec![0u64; num_positions],
            num_positions,
            pos_mask,
            max_depth_per_pos,
            depth_shift,
            head: 0,
            genome_pos: 0,
            len: 0,
        }
    }

    /// Return the slice of stored hashes at buffer position `buf_pos`.
    #[cfg(test)]
    fn hashes_at(&self, buf_pos: usize) -> &[u64] {
        let ring = (self.head + buf_pos) & self.pos_mask;
        let start = ring << self.depth_shift;
        let count = self.counts[ring] as usize;
        &self.hashes[start..start + count]
    }

    /// Add a name hash at `buf_pos`. If the position is full, increments excess.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "count ≤ max_depth_per_pos = 2 × coverage_cap (u32), fits in u32"
    )]
    fn push(&mut self, buf_pos: usize, name_hash: u64) {
        let ring = (self.head + buf_pos) & self.pos_mask;
        let count = self.counts[ring] as usize;
        if count < self.max_depth_per_pos {
            let idx = (ring << self.depth_shift) + count;
            self.hashes[idx] = name_hash;
            self.counts[ring] = (count + 1) as u32;
        } else {
            self.excess[ring] += 1;
        }
    }

    /// Clear all data at `buf_pos` (count and excess).
    #[cfg(test)]
    fn clear_position(&mut self, buf_pos: usize) {
        let ring = (self.head + buf_pos) & self.pos_mask;
        self.counts[ring] = 0;
        self.excess[ring] = 0;
    }

    /// Ensure the buffer covers positions up to `genome_pos + span - 1`.
    /// If the buffer is empty, initializes `self.genome_pos`.
    fn ensure_range(&mut self, start_pos: u64) {
        if self.len == 0 {
            self.genome_pos = start_pos;
        }
    }

    /// Flush a single position at `buf_pos`, returning depth statistics.
    /// After flushing, the position is cleared, `head` advances, and `len` shrinks.
    fn flush_position(&mut self, coverage_cap: u32) -> FlushResult {
        let ring = self.head & self.pos_mask;
        let count = self.counts[ring] as usize;
        let excess = self.excess[ring];

        // Fast paths for low-count positions (common at low-to-moderate coverage).
        let (depth, excl_overlap) = match count {
            0 => (0u64, 0u64),
            1 => (1, 0),
            2 => {
                let base = ring << self.depth_shift;
                if self.hashes[base] == self.hashes[base + 1] { (1, 1) } else { (2, 0) }
            }
            _ => {
                // Sort hashes in-place for overlap detection (the slot is about to be cleared).
                let start = ring << self.depth_shift;
                self.hashes[start..start + count].sort_unstable();

                // Linear scan: count unique hashes as depth, duplicates as overlap.
                let mut d: u64 = 0;
                let mut ov: u64 = 0;
                let mut prev_hash: u64 = u64::MAX;
                for &hash in &self.hashes[start..start + count] {
                    if hash == prev_hash {
                        ov += 1;
                    } else {
                        d += 1;
                        prev_hash = hash;
                    }
                }
                (d, ov)
            }
        };

        // Apply coverage cap.
        let cap = u64::from(coverage_cap);
        let capped = depth.min(cap);
        let excl_capped = depth - capped + excess;

        // Clear and advance.
        self.counts[ring] = 0;
        self.excess[ring] = 0;
        self.head = (self.head + 1) & self.pos_mask;
        self.genome_pos += 1;
        self.len -= 1;

        FlushResult { depth, excl_overlap, excl_capped }
    }
}

/// Result of flushing a single position from the depth buffer.
struct FlushResult {
    /// Number of unique name hashes (depth after overlap removal).
    depth: u64,
    /// Number of duplicate name hashes (overlap exclusions).
    excl_overlap: u64,
    /// Bases excluded because position was at or above the coverage cap, plus
    /// excess bases that arrived after the position's hash storage was full.
    excl_capped: u64,
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

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Compute a 64-bit hash of the given byte slice using `DefaultHasher`.
fn hash_name(name: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    name.hash(&mut hasher);
    hasher.finish()
}

// ─── CIGAR helper ─────────────────────────────────────────────────────────────

/// Count the number of reference-aligned bases (M/=/X) in a record's CIGAR.
fn count_aligned_bases(record: &RecordBuf) -> u64 {
    use noodles::sam::alignment::record::Cigar;
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut count: u64 = 0;
    for op in Cigar::iter(record.cigar()).filter_map(Result::ok) {
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

    // ── DepthBuffer tests ────────────────────────────────────────────────────

    #[test]
    fn test_depth_buffer_basic_push_and_flush() {
        let mut buf = DepthBuffer::new(100, 250);
        buf.ensure_range(10);
        buf.len = 1;
        buf.push(0, 42);
        buf.push(0, 99);

        let result = buf.flush_position(250);
        assert_eq!(result.depth, 2);
        assert_eq!(result.excl_overlap, 0);
        assert_eq!(result.excl_capped, 0);
    }

    #[test]
    fn test_depth_buffer_overlap_detection() {
        let mut buf = DepthBuffer::new(100, 250);
        buf.ensure_range(0);
        buf.len = 1;
        // Same hash twice = overlap
        buf.push(0, 42);
        buf.push(0, 42);

        let result = buf.flush_position(250);
        assert_eq!(result.depth, 1);
        assert_eq!(result.excl_overlap, 1);
        assert_eq!(result.excl_capped, 0);
    }

    #[test]
    fn test_depth_buffer_coverage_cap() {
        let mut buf = DepthBuffer::new(100, 3);
        buf.ensure_range(0);
        buf.len = 1;
        for i in 0..5u64 {
            buf.push(0, i);
        }

        let result = buf.flush_position(3);
        assert_eq!(result.depth, 5);
        assert_eq!(result.excl_overlap, 0);
        assert_eq!(result.excl_capped, 2); // 5 - 3 = 2
    }

    #[test]
    fn test_depth_buffer_excess_counting() {
        // max_depth_per_pos = (2 * 3).next_power_of_two() = 8, so after 8 entries excess kicks in.
        let mut buf = DepthBuffer::new(100, 3);
        buf.ensure_range(0);
        buf.len = 1;
        // Push 10 entries: 8 stored, 2 excess
        for i in 0..10u64 {
            buf.push(0, i);
        }

        let result = buf.flush_position(3);
        // Only 8 stored hashes, all unique → depth = 8
        assert_eq!(result.depth, 8);
        assert_eq!(result.excl_overlap, 0);
        // excl_capped = (8 - 3) + 2 excess = 7
        assert_eq!(result.excl_capped, 7);
    }

    #[test]
    fn test_depth_buffer_circular_wrap() {
        // Small num_positions to test circularity.
        let mut buf = DepthBuffer::new(4, 250);
        buf.ensure_range(100);

        // Fill 4 positions, flush them all, then add more.
        for i in 0..4usize {
            if i >= buf.len {
                buf.len = i + 1;
            }
            buf.push(i, (i + 1) as u64);
        }

        // Flush all 4 positions.
        for expected_depth in 1..=4u64 {
            let result = buf.flush_position(250);
            assert_eq!(result.depth, 1, "expected depth 1 at flush {expected_depth}");
        }

        // Now head has wrapped. Add more positions.
        buf.ensure_range(104);
        buf.len = 2;
        buf.push(0, 10);
        buf.push(1, 20);
        buf.push(1, 30);

        let r1 = buf.flush_position(250);
        assert_eq!(r1.depth, 1);
        let r2 = buf.flush_position(250);
        assert_eq!(r2.depth, 2);
    }

    #[test]
    fn test_depth_buffer_hashes_at() {
        let mut buf = DepthBuffer::new(10, 250);
        buf.ensure_range(0);
        buf.len = 2;
        buf.push(0, 100);
        buf.push(0, 200);
        buf.push(1, 300);

        assert_eq!(buf.hashes_at(0), &[100, 200]);
        assert_eq!(buf.hashes_at(1), &[300]);
    }

    #[test]
    fn test_depth_buffer_clear_position() {
        let mut buf = DepthBuffer::new(10, 250);
        buf.ensure_range(0);
        buf.len = 1;
        buf.push(0, 42);
        buf.clear_position(0);
        assert_eq!(buf.hashes_at(0).len(), 0);
    }
}
