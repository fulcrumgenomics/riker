use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bitvec::vec::BitVec;
use clap::Args;
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::plot::scatter::ScatterPlot;
use kuva::render::annotations::ReferenceLine;
use kuva::render::layout::TickFormat;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use crate::fasta::Fasta;
use crate::intervals::Intervals;
use crate::math::safe_div;
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, write_tsv};
use crate::plotting::{
    FG_BLUE, FG_GRAY, FG_GREEN, FG_SKY, FG_TEAL, standard_layout, write_twin_y_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sequence_dict::SequenceDictionary;

// ─── File suffixes ─────────────────────────────────────────────────────────────

/// File suffix for the per-GC-bin detail metrics output.
pub const DETAIL_SUFFIX: &str = ".gcbias-detail.txt";

/// File suffix for the summary metrics output.
pub const SUMMARY_SUFFIX: &str = ".gcbias-summary.txt";

/// File suffix for the GC bias chart.
pub const PLOT_SUFFIX: &str = ".gcbias-chart.pdf";

// ─── Number of GC bins ────────────────────────────────────────────────────────

const NUM_GC_BINS: usize = 101;

/// Maximum number of N bases allowed in a sliding window before it is considered invalid.
const MAX_N_IN_WINDOW: u32 = 4;

// ─── Helper functions ─────────────────────────────────────────────────────────

/// Bit flag indicating a G or C base in [`BASE_FLAGS`].
const GC_FLAG: u8 = 0x01;
/// Bit flag indicating an N base in [`BASE_FLAGS`].
const N_FLAG: u8 = 0x02;

/// 256-byte lookup table mapping each ASCII byte to bit flags encoding whether
/// it is a GC base ([`GC_FLAG`]) or an N base ([`N_FLAG`]).  Case-insensitive:
/// both `G`/`g` and `C`/`c` set [`GC_FLAG`]; both `N`/`n` set [`N_FLAG`].
///
/// Using a lookup table replaces per-base branch comparisons with a single
/// indexed load, which is significantly faster for the ~2.5 billion positions
/// scanned during a whole-genome GC bias calculation.
static BASE_FLAGS: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'G' as usize] = GC_FLAG;
    table[b'g' as usize] = GC_FLAG;
    table[b'C' as usize] = GC_FLAG;
    table[b'c' as usize] = GC_FLAG;
    table[b'N' as usize] = N_FLAG;
    table[b'n' as usize] = N_FLAG;
    table
};

// ─── Options ──────────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the GC bias collector.
///
/// GC bias measures library-prep bias via read-start positions, so the
/// defaults are deliberately permissive: duplicates and supplementary
/// reads are included and the MAPQ threshold is 0.
///
/// All counting and filtering is per read (read start), not per template:
/// each read or mate is evaluated and binned independently. A pair is never
/// dropped as a unit, so any filter (MAPQ, duplicate/supplementary, or
/// interval exclusion) can drop one mate while keeping the other — they
/// occupy different start positions and are accounted for separately.
#[riker_derive::multi_options("gcbias", "GC Bias Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct GcBiasOptions {
    /// Exclude duplicate reads from GC bias calculations.
    /// By default duplicates are included because GC bias measures
    /// library-level bias and every observed fragment contributes signal.
    #[arg(long, default_value_t = false)]
    pub exclude_duplicates: bool,

    /// Sliding window size for GC content calculation.
    #[arg(long, default_value_t = GcBiasOptions::DEFAULT_WINDOW_SIZE)]
    pub window_size: u32,

    /// Minimum mapping quality for a read to be counted.
    ///
    /// Defaults to 0 (all mapped reads) because GC bias measures
    /// library-prep bias and should count every read start regardless
    /// of mapping confidence.
    #[arg(long, default_value_t = GcBiasOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Exclude supplementary alignments from GC bias calculations.
    ///
    /// By default supplementary alignments are included because they
    /// represent real molecules and contribute meaningful signal.
    #[arg(long, default_value_t = false)]
    pub exclude_supplementary: bool,

    /// Exclude reads and reference windows whose computed start position falls
    /// within these intervals (BED or IntervalList, format auto-detected).
    ///
    /// Useful for masking artifact regions — e.g. poly-G stretches or adapter
    /// constructs — that produce spurious spikes in the GC bias signal. Both
    /// the read (numerator) and the reference window (denominator) at an
    /// excluded position are dropped, keeping the normalization consistent.
    ///
    /// The "computed start position" is the read's alignment start for
    /// forward-strand reads and `alignment_end - window_size` for reverse-strand
    /// reads (matching how each read is binned). A read is excluded based on
    /// that single position, not on whether its span overlaps the interval, so
    /// pad intervals if you need to catch reads whose body — but not start —
    /// touches an artifact locus.
    ///
    /// Exclusion is per read, not per template: each mate is tested
    /// independently, so one mate of a pair may be excluded while the other is
    /// still counted (the two mates have different start positions). If you
    /// need both mates dropped together, widen the intervals to cover both or
    /// pre-filter the BAM.
    #[arg(long, value_name = "FILE", value_parser = crate::commands::common::parse_existing_file)]
    pub exclude_intervals: Option<PathBuf>,
}

impl GcBiasOptions {
    const DEFAULT_WINDOW_SIZE: u32 = 100;
    const DEFAULT_MIN_MAPQ: u8 = 0;
}

impl Default for GcBiasOptions {
    fn default() -> Self {
        Self {
            exclude_duplicates: false,
            window_size: Self::DEFAULT_WINDOW_SIZE,
            min_mapq: Self::DEFAULT_MIN_MAPQ,
            exclude_supplementary: false,
            exclude_intervals: None,
        }
    }
}

// ─── CLI struct ───────────────────────────────────────────────────────────────

/// Collect GC bias metrics from a BAM file.
///
/// Measures GC bias by sliding a window across the reference genome and
/// comparing the expected GC distribution to the observed read-start
/// distribution. Produces per-GC-bin detail metrics, a summary row, and
/// a diagnostic chart. Outputs are written to <prefix>.gcbias-detail.txt,
/// <prefix>.gcbias-summary.txt, and <prefix>.gcbias-chart.pdf.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker gcbias -i input.bam -o out_prefix -r ref.fa
  riker gcbias -i input.bam -o out_prefix -r ref.fa --window-size 150"
)]
pub struct GcBias {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub options: GcBiasOptions,
}

impl Command for GcBias {
    /// # Errors
    /// Returns an error if the BAM or reference cannot be read, or if output cannot be written.
    fn execute(&self, threads: Option<u8>) -> Result<()> {
        super::common::validate_output_prefix(&self.output.output)?;
        let plan = self.thread_plan(threads);
        let mut reader = AlignmentReader::open(
            &self.input.input,
            Some(&self.reference.reference),
            plan.decode_threads,
        )?;
        let reference = Fasta::from_path(&self.reference.reference)?;

        let mut collector =
            GcBiasCollector::new(&self.input.input, &self.output.output, reference, &self.options);

        let mut progress = ProgressLogger::new("gcbias", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

// ─── Collector ────────────────────────────────────────────────────────────────

/// Accumulates GC bias statistics from a BAM file.
pub struct GcBiasCollector {
    // Output paths
    detail_path: PathBuf,
    summary_path: PathBuf,
    plot_path: PathBuf,
    plot_title: String,

    // Input path (needed for sample name derivation)
    input_path: PathBuf,

    // Config
    reference: Fasta,
    window_size: u32,
    exclude_duplicates: bool,
    min_mapq: u8,
    exclude_supplementary: bool,
    exclude_intervals_path: Option<PathBuf>,

    // BAM contig metadata (populated in initialize)
    dict: Option<SequenceDictionary>,

    // Intervals to exclude (loaded in initialize, merged for correct masking)
    exclude_intervals: Option<Intervals>,

    // Lazy per-contig GC lookup (recomputed on contig transition)
    current_contig_id: Option<usize>,
    current_gc_at_pos: Vec<u8>, // GC% at each position, or u8::MAX for invalid
    visited_contigs: HashSet<usize>,

    // Global accumulators [101] indexed by GC%
    windows_by_gc: [u64; NUM_GC_BINS],
    reads_by_gc: [u64; NUM_GC_BINS],
    bases_by_gc: [u64; NUM_GC_BINS],
    errors_by_gc: [u64; NUM_GC_BINS],
    quality_sum_by_gc: [u64; NUM_GC_BINS],
    quality_bases_by_gc: [u64; NUM_GC_BINS],

    // Read-accounting funnel (see GcBiasSummaryMetric): every record except
    // secondary/QC-fail counts in total_reads; the mapped subset counts in
    // aligned_reads; reads actually binned are `sum(reads_by_gc)`. The filtered
    // count is derived as `aligned_reads - sum(reads_by_gc)` at finish.
    total_reads: u64,
    total_clusters: u64,
    aligned_reads: u64,
    sample: String,
}

impl GcBiasCollector {
    /// Create a new collector. Output paths are derived from `prefix` by appending suffixes.
    #[must_use]
    pub fn new(input: &Path, prefix: &Path, reference: Fasta, options: &GcBiasOptions) -> Self {
        let detail_path = super::command::output_path(prefix, DETAIL_SUFFIX);
        let summary_path = super::command::output_path(prefix, SUMMARY_SUFFIX);
        let plot_path = super::command::output_path(prefix, PLOT_SUFFIX);
        Self {
            detail_path,
            summary_path,
            plot_path,
            plot_title: String::new(),
            input_path: input.to_path_buf(),
            reference,
            window_size: options.window_size,
            exclude_duplicates: options.exclude_duplicates,
            min_mapq: options.min_mapq,
            exclude_supplementary: options.exclude_supplementary,
            exclude_intervals_path: options.exclude_intervals.clone(),
            dict: None,
            exclude_intervals: None,
            current_contig_id: None,
            current_gc_at_pos: Vec::new(),
            visited_contigs: HashSet::new(),
            windows_by_gc: [0u64; NUM_GC_BINS],
            reads_by_gc: [0u64; NUM_GC_BINS],
            bases_by_gc: [0u64; NUM_GC_BINS],
            errors_by_gc: [0u64; NUM_GC_BINS],
            quality_sum_by_gc: [0u64; NUM_GC_BINS],
            quality_bases_by_gc: [0u64; NUM_GC_BINS],
            total_reads: 0,
            total_clusters: 0,
            aligned_reads: 0,
            sample: String::new(),
        }
    }

    /// Process a single BAM record.
    fn process_record(&mut self, record: &RikerRecord) -> Result<()> {
        let flags = record.flags();

        // Secondary and QC-fail records are discarded outright: they are never
        // candidates for GC binning and are counted nowhere in the funnel.
        if flags.is_secondary() || flags.is_qc_fail() {
            return Ok(());
        }

        // total_reads counts every surviving record (incl. unmapped/dup/supp/low-mapq).
        self.total_reads += 1;

        // Cluster counting: first-of-pair or unpaired. Counted on the total_reads
        // basis so it includes unmapped reads, mirroring total_reads.
        if !flags.is_segmented() || flags.is_first_segment() {
            self.total_clusters += 1;
        }

        // Unmapped reads have no start position to bin: counted in total_reads
        // but not in aligned_reads.
        if flags.is_unmapped() {
            return Ok(());
        }

        // aligned_reads counts every mapped record. It is the denominator for
        // frac_filtered_reads, so all GC-use filters below must come after it;
        // any read that exits past this point without binning is a "filtered"
        // read, recovered as `aligned_reads - reads_used` at finish.
        self.aligned_reads += 1;

        // GC-use filters: each is a bare return (the filtered count is derived).
        if self.exclude_duplicates && flags.is_duplicate() {
            return Ok(());
        }
        if self.exclude_supplementary && flags.is_supplementary() {
            return Ok(());
        }
        let mapq = record.mapping_quality().map_or(255u8, u8::from);
        if mapq < self.min_mapq {
            return Ok(());
        }

        // Contig transition: load GC array for the new contig
        let Some(ref_id) = record.reference_sequence_id() else {
            return Ok(());
        };

        if Some(ref_id) != self.current_contig_id {
            let name = self.dict.as_ref().unwrap().get_by_index(ref_id).map_or("", |m| m.name());
            let seq = self.reference.load_contig(name, false)?;
            let mask = self.contig_exclusion_mask(ref_id);
            let (gc_at_pos, window_counts) = scan_contig_gc(&seq, self.window_size, mask.as_ref());
            self.current_gc_at_pos = gc_at_pos;
            for (bin, count) in window_counts.iter().enumerate() {
                self.windows_by_gc[bin] += count;
            }
            self.visited_contigs.insert(ref_id);
            self.current_contig_id = Some(ref_id);
        }

        // Position: forward strand → alignment_start; reverse → alignment_end - window_size.
        // The reverse branch needs alignment_end (1-based inclusive == 0-based
        // exclusive), which is only `None` for records with no CIGAR; on those
        // we have nothing better than the start, so fall through to the
        // forward formula.
        let alignment_start = record.alignment_start().map_or(0, |p| usize::from(p) - 1); // 0-based
        let pos = if flags.is_reverse_complemented()
            && let Some(end) = record.alignment_end()
        {
            usize::from(end).saturating_sub(self.window_size as usize)
        } else {
            alignment_start
        };

        // Bounds check
        if pos >= self.current_gc_at_pos.len() {
            return Ok(());
        }

        // GC lookup. A sentinel value > 100 marks an invalid window start —
        // either too many Ns or an excluded interval — so the read is dropped
        // (and counted as filtered via the funnel) here.
        let gc = self.current_gc_at_pos[pos];
        if gc > 100 {
            return Ok(());
        }
        let gc = gc as usize;

        // Accumulate
        self.reads_by_gc[gc] += 1;

        // Read length (sequence bases)
        let read_len = record.sequence_len() as u64;
        self.bases_by_gc[gc] += read_len;

        // NM tag
        let nm = record.get_integer_tag(*b"NM").unwrap_or(0);
        self.errors_by_gc[gc] += u64::from(nm);

        // Base quality accumulation (integer)
        let qual_bytes: &[u8] = record.quality_scores();
        if !qual_bytes.is_empty() {
            self.quality_sum_by_gc[gc] += qual_bytes.iter().map(|&q| u64::from(q)).sum::<u64>();
            self.quality_bases_by_gc[gc] += qual_bytes.len() as u64;
        }

        Ok(())
    }

    /// Build the per-contig exclusion mask for `ref_id`, or `None` if there are
    /// no exclusion intervals on this contig.
    ///
    /// Returning `None` (rather than an empty bitvec) for an unmasked contig
    /// lets [`scan_contig_gc`] take its no-mask fast path, avoiding a
    /// bounds-checked bit lookup at every window position.
    fn contig_exclusion_mask(&self, ref_id: usize) -> Option<BitVec> {
        let mask = self.exclude_intervals.as_ref()?.contig_bitvec(ref_id);
        if mask.is_empty() { None } else { Some(mask) }
    }

    /// Finalize metrics computation and write outputs.
    fn finish_metrics(&self) -> Result<()> {
        let reads_used: u64 = self.reads_by_gc.iter().sum();
        let total_windows: u64 = self.windows_by_gc.iter().sum();

        // Compute mean reads per window
        let mean_reads_per_window = safe_div(reads_used, total_windows);

        // Filtered reads are aligned reads that never landed in a GC bin, for any
        // reason (duplicate/supplementary when excluded, low MAPQ, excluded
        // interval, or no valid GC window). Derived rather than counted: each
        // aligned read bins at most once, so `reads_used <= aligned_reads` always
        // holds and the subtraction cannot underflow (saturating_sub is belt-and-
        // suspenders for that invariant).
        let filtered_reads = self.aligned_reads.saturating_sub(reads_used);
        let frac_filtered_reads = safe_div(filtered_reads, self.aligned_reads);

        // Build detail metrics (101 rows)
        let detail_rows: Vec<GcBiasDetailMetric> = (0..NUM_GC_BINS)
            .map(|gc| {
                let windows = self.windows_by_gc[gc];
                let reads = self.reads_by_gc[gc];
                let bases = self.bases_by_gc[gc];
                let errors = self.errors_by_gc[gc];

                let (normalized_coverage, error_bar_width) =
                    normalized_coverage_and_error(reads, windows, mean_reads_per_window);

                let reported_base_quality =
                    mean_base_quality(self.quality_sum_by_gc[gc], self.quality_bases_by_gc[gc]);

                let empirical_base_quality = phred_from_error_counts(errors, bases);

                GcBiasDetailMetric {
                    sample: self.sample.clone(),
                    gc: gc as u64,
                    windows,
                    read_starts: reads,
                    reported_base_quality,
                    empirical_base_quality,
                    normalized_coverage,
                    error_bar_width,
                }
            })
            .collect();

        let (at_dropout, gc_dropout) = compute_dropout(&self.windows_by_gc, &self.reads_by_gc);

        // Compute quintile NC (aggregate, not average)
        let quintile_nc = |start: usize, end: usize| -> f64 {
            if mean_reads_per_window == 0.0 {
                return 0.0;
            }
            let sum_reads: u64 = self.reads_by_gc[start..=end].iter().sum();
            let sum_windows: u64 = self.windows_by_gc[start..=end].iter().sum();
            if sum_windows == 0 {
                0.0
            } else {
                sum_reads as f64 / (sum_windows as f64 * mean_reads_per_window)
            }
        };

        let summary = GcBiasSummaryMetric {
            sample: self.sample.clone(),
            window_size: u64::from(self.window_size),
            total_clusters: self.total_clusters,
            total_reads: self.total_reads,
            aligned_reads: self.aligned_reads,
            filtered_reads,
            frac_filtered_reads,
            at_dropout,
            gc_dropout,
            gc_0_19_normcov: quintile_nc(0, 19),
            gc_20_39_normcov: quintile_nc(20, 39),
            gc_40_59_normcov: quintile_nc(40, 59),
            gc_60_79_normcov: quintile_nc(60, 79),
            gc_80_100_normcov: quintile_nc(80, 100),
        };

        // Write outputs
        write_tsv(&self.detail_path, &detail_rows)?;
        write_tsv(&self.summary_path, &[summary])?;

        // Generate plot
        self.plot_chart(&detail_rows)?;

        log::info!(
            "gcbias: total_reads={}, total_clusters={}, aligned_reads={}, \
             filtered_reads={filtered_reads} ({:.1}%), at_dropout={at_dropout:.3}, \
             gc_dropout={gc_dropout:.3}, detail={}, summary={}, plot={}",
            self.total_reads,
            self.total_clusters,
            self.aligned_reads,
            frac_filtered_reads * 100.0,
            self.detail_path.display(),
            self.summary_path.display(),
            self.plot_path.display(),
        );

        Ok(())
    }

    /// Write a PDF GC bias chart with dual Y-axes.
    fn plot_chart(&self, detail_rows: &[GcBiasDetailMetric]) -> Result<()> {
        let y_max = 2.0_f64;
        let max_windows = detail_rows.iter().map(|r| r.windows).max().unwrap_or(1).max(1);
        let window_scale = (y_max * 0.25) / max_windows as f64;

        // Primary series 1: Window distribution as filled step chart
        let window_xy: Vec<(f64, f64)> =
            detail_rows.iter().map(|r| (r.gc as f64, r.windows as f64 * window_scale)).collect();
        let windows_line = LinePlot::new()
            .with_data(window_xy)
            .with_color(FG_SKY)
            .with_stroke_width(1.0)
            .with_step()
            .with_fill()
            .with_fill_opacity(0.3)
            .with_legend("Genome GC");

        // Primary series 2: Normalized coverage scatter dots
        let nc_xy: Vec<(f64, f64)> = detail_rows
            .iter()
            .filter(|r| r.windows > 0)
            .map(|r| (r.gc as f64, r.normalized_coverage.min(y_max)))
            .collect();
        let nc_scatter = ScatterPlot::new()
            .with_data(nc_xy)
            .with_color(FG_BLUE)
            .with_size(4.0)
            .with_legend("Coverage");

        let primary: Vec<Plot> = vec![windows_line.into(), nc_scatter.into()];

        // Secondary series: base quality lines
        let mut secondary: Vec<Plot> = Vec::new();

        let reported_bq: Vec<(f64, f64)> = detail_rows
            .iter()
            .filter(|r| r.windows > 0 && r.reported_base_quality > 0.0)
            .map(|r| (r.gc as f64, r.reported_base_quality))
            .collect();
        if !reported_bq.is_empty() {
            secondary.push(
                LinePlot::new()
                    .with_data(reported_bq)
                    .with_color(FG_GREEN)
                    .with_stroke_width(1.0)
                    .with_legend("Reported BQ")
                    .into(),
            );
        }

        let empirical_bq: Vec<(f64, f64)> = detail_rows
            .iter()
            .filter(|r| r.windows > 0 && r.empirical_base_quality > 0.0)
            .map(|r| (r.gc as f64, r.empirical_base_quality))
            .collect();
        if !empirical_bq.is_empty() {
            secondary.push(
                LinePlot::new()
                    .with_data(empirical_bq)
                    .with_color(FG_TEAL)
                    .with_stroke_width(1.0)
                    .with_legend("Empirical BQ")
                    .into(),
            );
        }

        let layout = standard_layout(&primary)
            .with_x_axis_min(0.0)
            .with_x_axis_max(100.0)
            .with_y_axis_min(0.0)
            .with_y_axis_max(y_max)
            .with_y2_range(0.0, 40.0)
            .with_x_tick_format(TickFormat::Integer)
            .with_y_tick_format(TickFormat::Fixed(1))
            .with_y2_tick_format(TickFormat::Integer)
            .with_title(&self.plot_title)
            .with_x_label("GC%")
            .with_y_label("Normalized Coverage")
            .with_y2_label("Base Quality")
            .with_y2_label_offset(-15.0, 0.0)
            .with_reference_line(
                ReferenceLine::horizontal(1.0).with_color(FG_GRAY).with_dasharray(""),
            )
            .with_legend_position(LegendPosition::InsideBottomRight)
            .with_legend_box(false);

        write_twin_y_plot_pdf(primary, secondary, layout, &self.plot_path)
    }
}

// ─── Collector trait impl ─────────────────────────────────────────────────────

impl Collector for GcBiasCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        self.reference.validate_bam_header(header)?;

        let dict = SequenceDictionary::from(header);

        // Load exclusion intervals (if any) and merge them so the per-contig
        // bitvec mask is correct even when the input intervals overlap.
        if let Some(path) = &self.exclude_intervals_path {
            let intervals = Intervals::from_path(path, dict.clone())?.merged();
            log::info!(
                "gcbias: excluding {} interval(s) covering {} bases from GC bias",
                intervals.count(),
                intervals.territory(),
            );
            self.exclude_intervals = Some(intervals);
        }

        self.dict = Some(dict);

        self.sample = derive_sample(&self.input_path, header);
        self.plot_title = format!("GC Bias of {}", self.sample);
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        self.process_record(record)
    }

    fn finish(&mut self) -> Result<()> {
        // Scan any BAM header contigs not visited during record traversal so
        // their windows still contribute to the denominator. The exclusion mask
        // is applied here too, keeping the denominator consistent on contigs
        // that had no reads.
        let dict = self.dict.as_ref().unwrap();
        for ref_id in 0..dict.len() {
            if !self.visited_contigs.contains(&ref_id) {
                let name = dict[ref_id].name();
                let seq = self.reference.load_contig(name, false)?;
                let mask = self.contig_exclusion_mask(ref_id);
                let (_, window_counts) = scan_contig_gc(&seq, self.window_size, mask.as_ref());
                for (bin, count) in window_counts.iter().enumerate() {
                    self.windows_by_gc[bin] += count;
                }
            }
        }
        self.finish_metrics()
    }

    fn name(&self) -> &'static str {
        "gcbias"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        // Uses `sequence_len()` (no decode) + `NM` aux tag + quality scores
        // (always available). Sequence bases are never read, so we don't
        // declare `with_sequence`.
        RikerRecordRequirements::NONE.with_aux_tag(*b"NM")
    }

    /// Relative per-record worker cost (see [`Collector::cost_hint`]); ~22 from
    /// a samply worker-compute measurement on a 12x WGS BAM. gcbias is
    /// reader-bound — its worker sits idle much of the time — so its worker
    /// share is modest even though its overall runtime is not.
    fn cost_hint(&self) -> u32 {
        22
    }
}

// ─── Metric structs ───────────────────────────────────────────────────────────

/// GC bias detail metrics — one row per GC percentage bin (0-100).
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct GcBiasDetailMetric {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// GC content percentage (0-100).
    pub gc: u64,
    /// Number of reference windows at this GC percentage.
    pub windows: u64,
    /// Number of reads starting at positions with this GC percentage.
    pub read_starts: u64,
    /// Mean of actual base quality scores from reads in this bin.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub reported_base_quality: f64,
    /// Phred-scaled error rate derived from the NM tag (mismatches + indels).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub empirical_base_quality: f64,
    /// Normalized coverage: (reads/windows) / mean_reads_per_window.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub normalized_coverage: f64,
    /// Error bar width: sqrt(reads)/windows / mean_reads_per_window.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub error_bar_width: f64,
}

/// GC bias summary metrics — one row per sample.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct GcBiasSummaryMetric {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Sliding window size used for GC content calculation.
    pub window_size: u64,
    /// Total clusters: first-of-pair or unpaired reads among `total_reads`
    /// (includes unmapped reads).
    pub total_clusters: u64,
    /// Total reads considered: every read except secondary and QC-fail
    /// (includes unmapped, duplicate, supplementary, and low-MAPQ reads).
    pub total_reads: u64,
    /// Aligned reads: the mapped subset of `total_reads`.
    pub aligned_reads: u64,
    /// Filtered reads: aligned reads not used in the GC bias computation
    /// (excluded as duplicate/supplementary, low MAPQ, in an excluded interval,
    /// or lacking a valid GC window). Equals `aligned_reads` minus the total
    /// binned read starts.
    pub filtered_reads: u64,
    /// Fraction of aligned reads that were filtered out of the GC bias
    /// computation (`filtered_reads / aligned_reads`).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_filtered_reads: f64,
    /// AT dropout: deficit at GC 0-50%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub at_dropout: f64,
    /// GC dropout: deficit at GC 51-100%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_dropout: f64,
    /// Aggregate normalized coverage for GC 0-19%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_0_19_normcov: f64,
    /// Aggregate normalized coverage for GC 20-39%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_20_39_normcov: f64,
    /// Aggregate normalized coverage for GC 40-59%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_40_59_normcov: f64,
    /// Aggregate normalized coverage for GC 60-79%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_60_79_normcov: f64,
    /// Aggregate normalized coverage for GC 80-100%.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub gc_80_100_normcov: f64,
}

// ─── Extracted computation helpers ─────────────────────────────────────────────

/// Compute normalized coverage and error bar width for a single GC bin.
///
/// Returns `(normalized_coverage, error_bar_width)`. If `windows` is 0 or
/// `mean_rpw` is 0.0, both values are 0.0.
fn normalized_coverage_and_error(reads: u64, windows: u64, mean_rpw: f64) -> (f64, f64) {
    if windows == 0 || mean_rpw == 0.0 {
        return (0.0, 0.0);
    }
    let nc = (reads as f64 / windows as f64) / mean_rpw;
    let eb = ((reads as f64).sqrt() / windows as f64) / mean_rpw;
    (nc, eb)
}

/// Mean reported base quality for a GC bin: the average of the per-base quality scores of
/// the reads binned here, or `0.0` when the bin has no quality bases.
fn mean_base_quality(quality_sum: u64, quality_bases: u64) -> f64 {
    safe_div(quality_sum, quality_bases)
}

/// Empirical (Phred-scaled) base quality for a GC bin, derived from the observed error rate
/// `errors / bases`: `-10 * log10(errors / bases)`. Returns `0.0` when the bin has no bases
/// or no errors, where the Phred score would otherwise be undefined or infinite.
fn phred_from_error_counts(errors: u64, bases: u64) -> f64 {
    if bases > 0 && errors > 0 { -10.0 * (errors as f64 / bases as f64).log10() } else { 0.0 }
}

/// Compute AT and GC dropout from per-GC-bin window and read counts.
///
/// Dropout measures the deficit in read coverage relative to the reference
/// window distribution. AT dropout sums deficits for GC bins 0-50%;
/// GC dropout sums deficits for bins 51-100%.
fn compute_dropout(
    windows_by_gc: &[u64; NUM_GC_BINS],
    reads_by_gc: &[u64; NUM_GC_BINS],
) -> (f64, f64) {
    let total_windows: u64 = windows_by_gc.iter().sum();
    let total_reads: u64 = reads_by_gc.iter().sum();
    if total_windows == 0 || total_reads == 0 {
        return (0.0, 0.0);
    }
    let mut at_drop = 0.0_f64;
    let mut gc_drop = 0.0_f64;
    for gc in 0..NUM_GC_BINS {
        let relative_windows = windows_by_gc[gc] as f64 / total_windows as f64;
        let relative_reads = reads_by_gc[gc] as f64 / total_reads as f64;
        let dropout = (relative_windows - relative_reads) * 100.0;
        if dropout > 0.0 {
            if gc <= 50 {
                at_drop += dropout;
            } else {
                gc_drop += dropout;
            }
        }
    }
    (at_drop, gc_drop)
}

/// Compute GC percentage (0–100) as a `u8`.
///
/// Uses integer arithmetic: `gc_count * 100 / window_size`.  The result is
/// guaranteed to fit in a `u8` for any `gc_count ≤ window_size`.
fn gc_percentage(gc_count: u32, window_size: u32) -> u8 {
    // Safety: gc_count ≤ window_size, so the result is always in 0..=100.
    #[allow(clippy::cast_possible_truncation)]
    let pct = ((u64::from(gc_count) * 100) / u64::from(window_size)) as u8;
    pct
}

/// Scan a contig sequence with a sliding window, returning per-position GC percentages
/// (`u8::MAX` for invalid positions with too many Ns) and per-GC-bin window counts.
///
/// Uses [`BASE_FLAGS`] for branchless classification of each base as GC or N,
/// which is significantly faster than per-base `match` comparisons for
/// genome-scale contigs.  The lookup table is case-insensitive, so the input
/// sequence does not need to be uppercased first.
///
/// If `excluded` is supplied (a per-contig bitvec where set bits mark excluded
/// positions), any window whose start position is excluded is treated as
/// invalid: it is left as `u8::MAX` in `gc_at_pos` and omitted from
/// `window_counts`. Because reads are binned by their computed start position,
/// this simultaneously drops the read (numerator, via the `u8::MAX` sentinel)
/// and the reference window (denominator), keeping normalization consistent.
fn scan_contig_gc(
    seq: &[u8],
    window_size: u32,
    excluded: Option<&BitVec>,
) -> (Vec<u8>, [u64; NUM_GC_BINS]) {
    let ws = window_size as usize;
    let mut gc_at_pos = vec![u8::MAX; seq.len()];
    let mut window_counts = [0u64; NUM_GC_BINS];

    if seq.len() < ws {
        return (gc_at_pos, window_counts);
    }

    // A window start position is excluded if its bit is set in the mask.
    // `is_some_and` short-circuits to false (cheap) when no mask is supplied.
    let is_excluded = |pos: usize| excluded.is_some_and(|bv| bv.get(pos).is_some_and(|b| *b));

    // Initialize counts for the first window using the lookup table
    let mut gc_count: u32 = 0;
    let mut n_count: u32 = 0;
    for &b in &seq[..ws] {
        let flags = BASE_FLAGS[b as usize];
        gc_count += u32::from(flags & GC_FLAG);
        n_count += u32::from((flags & N_FLAG) >> 1);
    }

    // Record the first window
    if n_count <= MAX_N_IN_WINDOW && !is_excluded(0) {
        let gc_pct = gc_percentage(gc_count, window_size);
        gc_at_pos[0] = gc_pct;
        window_counts[gc_pct as usize] += 1;
    }

    // Slide the window across the contig using branchless table lookups
    for i in 1..=(seq.len() - ws) {
        let leaving_flags = BASE_FLAGS[seq[i - 1] as usize];
        let entering_flags = BASE_FLAGS[seq[i + ws - 1] as usize];

        gc_count -= u32::from(leaving_flags & GC_FLAG);
        gc_count += u32::from(entering_flags & GC_FLAG);
        n_count -= u32::from((leaving_flags & N_FLAG) >> 1);
        n_count += u32::from((entering_flags & N_FLAG) >> 1);

        if n_count <= MAX_N_IN_WINDOW && !is_excluded(i) {
            let gc_pct = gc_percentage(gc_count, window_size);
            gc_at_pos[i] = gc_pct;
            window_counts[gc_pct as usize] += 1;
        }
    }

    (gc_at_pos, window_counts)
}

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn test_base_flags_gc() {
        assert_eq!(BASE_FLAGS[b'G' as usize] & GC_FLAG, GC_FLAG);
        assert_eq!(BASE_FLAGS[b'g' as usize] & GC_FLAG, GC_FLAG);
        assert_eq!(BASE_FLAGS[b'C' as usize] & GC_FLAG, GC_FLAG);
        assert_eq!(BASE_FLAGS[b'c' as usize] & GC_FLAG, GC_FLAG);
        assert_eq!(BASE_FLAGS[b'A' as usize] & GC_FLAG, 0);
        assert_eq!(BASE_FLAGS[b'T' as usize] & GC_FLAG, 0);
        assert_eq!(BASE_FLAGS[b'N' as usize] & GC_FLAG, 0);
    }

    #[test]
    fn test_base_flags_n() {
        assert_eq!(BASE_FLAGS[b'N' as usize] & N_FLAG, N_FLAG);
        assert_eq!(BASE_FLAGS[b'n' as usize] & N_FLAG, N_FLAG);
        assert_eq!(BASE_FLAGS[b'A' as usize] & N_FLAG, 0);
        assert_eq!(BASE_FLAGS[b'G' as usize] & N_FLAG, 0);
    }

    #[test]
    fn test_scan_contig_gc_case_insensitive() {
        // Lowercase gc bases should produce the same result as uppercase
        let upper = b"GGGGAAAA";
        let lower = b"ggggaaaa";
        let mixed = b"GgGgAaAa";
        let (gc_u, counts_u) = scan_contig_gc(upper, 4, None);
        let (gc_l, counts_l) = scan_contig_gc(lower, 4, None);
        let (gc_m, counts_m) = scan_contig_gc(mixed, 4, None);
        assert_eq!(gc_u, gc_l);
        assert_eq!(gc_u, gc_m);
        assert_eq!(counts_u, counts_l);
        assert_eq!(counts_u, counts_m);
    }

    #[test]
    fn test_scan_contig_gc_simple() {
        // All-A sequence with window_size=4 → GC% = 0 at every position
        let seq = b"AAAAAAAA";
        let (gc, counts) = scan_contig_gc(seq, 4, None);
        assert_eq!(gc.len(), 8);
        assert_eq!(gc[0], 0); // first window [0..4] = AAAA → 0% GC
        assert_eq!(gc[4], 0); // last window [4..8] = AAAA → 0% GC
        assert_eq!(gc[5], u8::MAX); // no valid window starting here
        assert_eq!(counts[0], 5); // 5 valid windows all at gc=0
    }

    #[test]
    fn test_scan_contig_gc_mixed() {
        // GGGGAAAA with window_size=4
        // pos 0: GGGG → 100% GC
        // pos 1: GGGA → 75% GC
        // pos 2: GGAA → 50% GC
        // pos 3: GAAA → 25% GC
        // pos 4: AAAA → 0% GC
        let seq = b"GGGGAAAA";
        let (gc, counts) = scan_contig_gc(seq, 4, None);
        assert_eq!(gc[0], 100);
        assert_eq!(gc[1], 75);
        assert_eq!(gc[2], 50);
        assert_eq!(gc[3], 25);
        assert_eq!(gc[4], 0);
        assert_eq!(counts[0], 1);
        assert_eq!(counts[25], 1);
        assert_eq!(counts[50], 1);
        assert_eq!(counts[75], 1);
        assert_eq!(counts[100], 1);
    }

    #[test]
    fn test_scan_contig_gc_with_ns() {
        // 5 N's then 4 A's, window_size=4
        let seq = b"NNNNNAAAA";
        let (gc, _) = scan_contig_gc(seq, 4, None);
        // pos 0: NNNN → 4 Ns, valid (≤4)
        assert_eq!(gc[0], 0); // 0 GC out of 4 (Ns aren't GC)
        // pos 1: NNNA → 3 Ns, valid
        assert_eq!(gc[1], 0);
        // pos 2: NNAA → 2 Ns, valid
        assert_eq!(gc[2], 0);
    }

    #[test]
    fn test_scan_contig_gc_too_many_ns() {
        // All N's with window_size=5 → each window has 5 Ns > 4 → invalid
        let seq = b"NNNNNNNNN";
        let (gc, counts) = scan_contig_gc(seq, 5, None);
        for &v in &gc[..5] {
            assert_eq!(v, u8::MAX, "window with >4 Ns should be invalid");
        }
        // No valid windows
        assert_eq!(counts.iter().sum::<u64>(), 0);
    }

    #[test]
    fn test_scan_contig_gc_short_seq() {
        // Sequence shorter than window_size
        let seq = b"GC";
        let (gc, counts) = scan_contig_gc(seq, 4, None);
        assert_eq!(gc.len(), 2);
        assert_eq!(gc[0], u8::MAX);
        assert_eq!(gc[1], u8::MAX);
        assert_eq!(counts.iter().sum::<u64>(), 0);
    }

    // ── scan_contig_gc — additional edge cases ──────────────────────────────

    #[test]
    fn test_scan_contig_gc_window_size_1() {
        // Window size 1: each position is its own window
        let seq = b"GCAT";
        let (gc, counts) = scan_contig_gc(seq, 1, None);
        assert_eq!(gc[0], 100); // G → 100% GC
        assert_eq!(gc[1], 100); // C → 100% GC
        assert_eq!(gc[2], 0); // A → 0% GC
        assert_eq!(gc[3], 0); // T → 0% GC
        assert_eq!(counts[0], 2);
        assert_eq!(counts[100], 2);
    }

    #[test]
    fn test_scan_contig_gc_all_gc() {
        // All G/C → every window is 100% GC
        let seq = b"GCGCGCGC";
        let (gc, counts) = scan_contig_gc(seq, 4, None);
        for (pos, &val) in gc.iter().enumerate().take(5) {
            assert_eq!(val, 100, "pos {pos} should be 100% GC");
        }
        assert_eq!(counts[100], 5);
    }

    #[test]
    fn test_scan_contig_gc_boundary_n_count() {
        // Exactly MAX_N_IN_WINDOW (4) Ns → valid; 5 Ns → invalid
        // Window size=5: "NNNNA" has 4 Ns → valid; "NNNNN" has 5 → invalid
        let seq = b"NNNNANNNN";
        let (gc, _) = scan_contig_gc(seq, 5, None);
        // pos 0: NNNNA → 4 Ns, valid
        assert_ne!(gc[0], u8::MAX);
        assert_eq!(gc[0], 0); // 0 GC bases out of 5
        // pos 4: ANNNN → 4 Ns, valid
        assert_ne!(gc[4], u8::MAX);
        assert_eq!(gc[4], 0);
        // But "NNNNN" (5 Ns) would be invalid
        let seq2 = b"NNNNN";
        let (gc2, _) = scan_contig_gc(seq2, 5, None);
        assert_eq!(gc2[0], u8::MAX);
    }

    #[test]
    fn test_scan_contig_gc_with_exclusion_mask() {
        use bitvec::prelude::*;
        // 8 G's, window_size 4 → 5 windows, all gc=100.
        let seq = b"GGGGGGGG";
        // Exclude window-start positions 0 and 1.
        let mut mask: BitVec = bitvec![0; seq.len()];
        mask.set(0, true);
        mask.set(1, true);

        let (gc, counts) = scan_contig_gc(seq, 4, Some(&mask));

        // Excluded starts are invalid in gc_at_pos (so reads there are dropped)...
        assert_eq!(gc[0], u8::MAX, "excluded window start 0 marked invalid");
        assert_eq!(gc[1], u8::MAX, "excluded window start 1 marked invalid");
        assert_eq!(gc[2], 100);
        assert_eq!(gc[3], 100);
        assert_eq!(gc[4], 100);
        // ...and dropped from the window counts (denominator): 5 − 2 = 3.
        assert_eq!(counts[100], 3, "excluded windows removed from the denominator");
    }

    // ── normalized_coverage_and_error ────────────────────────────────────────

    #[test]
    fn test_nc_normal() {
        let (nc, eb) = normalized_coverage_and_error(100, 50, 1.0);
        crate::assert_close!(nc, 2.0, 1e-9);
        crate::assert_close!(eb, 0.2, 1e-9); // sqrt(100)/50 / 1.0 = 10/50 = 0.2
    }

    #[test]
    fn test_nc_zero_windows() {
        let (nc, eb) = normalized_coverage_and_error(100, 0, 1.0);
        assert_eq!(nc, 0.0);
        assert_eq!(eb, 0.0);
    }

    #[test]
    fn test_nc_zero_mean() {
        let (nc, eb) = normalized_coverage_and_error(100, 50, 0.0);
        assert_eq!(nc, 0.0);
        assert_eq!(eb, 0.0);
    }

    // ── mean_base_quality / phred_from_error_counts ──────────────────────────

    #[test]
    fn mean_base_quality_averages_the_quality_sum() {
        // 400 total quality over 10 bases → mean Q40.
        crate::assert_close!(mean_base_quality(400, 10), 40.0);
    }

    #[test]
    fn mean_base_quality_is_zero_without_quality_bases() {
        assert_eq!(mean_base_quality(999, 0), 0.0);
    }

    #[test]
    fn phred_from_error_counts_is_phred_scaled_error_rate() {
        crate::assert_close!(phred_from_error_counts(1, 100), 20.0); // 1% error → Q20
        crate::assert_close!(phred_from_error_counts(1, 1000), 30.0); // 0.1% error → Q30
    }

    #[test]
    fn phred_from_error_counts_saturated_error_rate_is_zero() {
        // Every base an error → rate 1.0 → -10*log10(1) = Q0, via the formula not the guard.
        crate::assert_close!(phred_from_error_counts(10, 10), 0.0);
    }

    #[test]
    fn phred_from_error_counts_is_zero_without_errors_or_bases() {
        assert_eq!(phred_from_error_counts(0, 100), 0.0); // no errors
        assert_eq!(phred_from_error_counts(5, 0), 0.0); // no bases
    }

    // ── compute_dropout ─────────────────────────────────────────────────────

    #[test]
    fn test_compute_dropout_uniform() {
        // Proportional reads and windows → no dropout
        let mut windows = [0u64; NUM_GC_BINS];
        let mut reads = [0u64; NUM_GC_BINS];
        windows[25] = 50;
        windows[75] = 50;
        reads[25] = 50;
        reads[75] = 50;
        let (at, gc) = compute_dropout(&windows, &reads);
        assert_eq!(at, 0.0);
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn test_compute_dropout_at_deficit() {
        // More windows at low GC than reads → AT dropout > 0
        let mut windows = [0u64; NUM_GC_BINS];
        let mut reads = [0u64; NUM_GC_BINS];
        windows[25] = 60; // AT-rich region: 60% of windows
        windows[75] = 40;
        reads[25] = 40; // underrepresented in reads
        reads[75] = 60;
        let (at, gc) = compute_dropout(&windows, &reads);
        assert!(at > 0.0, "AT dropout should be > 0, got {at}");
        assert!(gc == 0.0, "GC dropout should be 0, got {gc}");
    }

    #[test]
    fn test_compute_dropout_gc_deficit() {
        // More windows at high GC than reads → GC dropout > 0
        let mut windows = [0u64; NUM_GC_BINS];
        let mut reads = [0u64; NUM_GC_BINS];
        windows[25] = 40;
        windows[75] = 60; // GC-rich region: 60% of windows
        reads[25] = 60;
        reads[75] = 40; // underrepresented
        let (at, gc) = compute_dropout(&windows, &reads);
        assert!(at == 0.0, "AT dropout should be 0, got {at}");
        assert!(gc > 0.0, "GC dropout should be > 0, got {gc}");
    }

    #[test]
    fn test_compute_dropout_empty() {
        let windows = [0u64; NUM_GC_BINS];
        let reads = [0u64; NUM_GC_BINS];
        let (at, gc) = compute_dropout(&windows, &reads);
        assert_eq!(at, 0.0);
        assert_eq!(gc, 0.0);
    }
}
