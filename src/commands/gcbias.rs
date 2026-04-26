use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Args;
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::plot::scatter::ScatterPlot;
use kuva::render::annotations::ReferenceLine;
use kuva::render::layout::{Layout, TickFormat};
use kuva::render::plots::Plot;
use noodles::sam::Header;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use crate::fasta::Fasta;
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, write_tsv};
use crate::plotting::{
    FG_BLUE, FG_GRAY, FG_GREEN, FG_SKY, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_twin_y_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::{derive_sample, get_integer_tag};
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
/// reads are included, MAPQ threshold is 0.  This matches Picard's
/// CollectGcBiasMetrics behaviour.
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
    fn execute(&self) -> Result<()> {
        let mut reader = AlignmentReader::open(&self.input.input, Some(&self.reference.reference))?;
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

    // BAM contig metadata (populated in initialize)
    dict: Option<SequenceDictionary>,

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
            dict: None,
            current_contig_id: None,
            current_gc_at_pos: Vec::new(),
            visited_contigs: HashSet::new(),
            windows_by_gc: [0u64; NUM_GC_BINS],
            reads_by_gc: [0u64; NUM_GC_BINS],
            bases_by_gc: [0u64; NUM_GC_BINS],
            errors_by_gc: [0u64; NUM_GC_BINS],
            quality_sum_by_gc: [0u64; NUM_GC_BINS],
            quality_bases_by_gc: [0u64; NUM_GC_BINS],
            total_clusters: 0,
            aligned_reads: 0,
            sample: String::new(),
        }
    }

    /// Process a single BAM record.
    fn process_record(&mut self, record: &RikerRecord) -> Result<()> {
        let flags = record.flags();

        // Filter: skip unmapped, secondary, QC-fail; supplementary only if excluded
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_qc_fail()
            || (self.exclude_supplementary && flags.is_supplementary())
        {
            return Ok(());
        }

        // Filter: duplicates
        if self.exclude_duplicates && flags.is_duplicate() {
            return Ok(());
        }

        // Filter: MAPQ
        let mapq = record.mapping_quality().map_or(255u8, u8::from);
        if mapq < self.min_mapq {
            return Ok(());
        }

        // Cluster counting: first-of-pair or unpaired
        if !flags.is_segmented() || flags.is_first_segment() {
            self.total_clusters += 1;
        }

        self.aligned_reads += 1;

        // Contig transition: load GC array for the new contig
        let Some(ref_id) = record.reference_sequence_id() else {
            return Ok(());
        };

        if Some(ref_id) != self.current_contig_id {
            let name = self.dict.as_ref().unwrap().get_by_index(ref_id).map_or("", |m| m.name());
            let seq = self.reference.load_contig(name, false)?;
            let (gc_at_pos, window_counts) = scan_contig_gc(&seq, self.window_size);
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

        // GC lookup
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
        let nm = get_integer_tag(record, *b"NM").unwrap_or(0);
        self.errors_by_gc[gc] += u64::from(nm);

        // Base quality accumulation (integer)
        let qual_bytes: &[u8] = record.quality_scores();
        if !qual_bytes.is_empty() {
            self.quality_sum_by_gc[gc] += qual_bytes.iter().map(|&q| u64::from(q)).sum::<u64>();
            self.quality_bases_by_gc[gc] += qual_bytes.len() as u64;
        }

        Ok(())
    }

    /// Finalize metrics computation and write outputs.
    fn finish_metrics(&self) -> Result<()> {
        let total_reads: u64 = self.reads_by_gc.iter().sum();
        let total_windows: u64 = self.windows_by_gc.iter().sum();

        // Compute mean reads per window
        let mean_reads_per_window =
            if total_windows > 0 { total_reads as f64 / total_windows as f64 } else { 0.0 };

        // Build detail metrics (101 rows)
        let detail_rows: Vec<GcBiasDetailMetric> = (0..NUM_GC_BINS)
            .map(|gc| {
                let windows = self.windows_by_gc[gc];
                let reads = self.reads_by_gc[gc];
                let bases = self.bases_by_gc[gc];
                let errors = self.errors_by_gc[gc];

                let (normalized_coverage, error_bar_width) =
                    normalized_coverage_and_error(reads, windows, mean_reads_per_window);

                let reported_base_quality = {
                    let qbases = self.quality_bases_by_gc[gc];
                    if qbases > 0 { self.quality_sum_by_gc[gc] as f64 / qbases as f64 } else { 0.0 }
                };

                let empirical_base_quality = if bases > 0 && errors > 0 {
                    -10.0 * (errors as f64 / bases as f64).log10()
                } else {
                    0.0
                };

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
            aligned_reads: self.aligned_reads,
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
            "gcbias: total_clusters={}, aligned_reads={}, at_dropout={at_dropout:.3}, \
             gc_dropout={gc_dropout:.3}, detail={}, summary={}, plot={}",
            self.total_clusters,
            self.aligned_reads,
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

        let layout = Layout::auto_from_plots(&primary)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
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

        self.dict = Some(SequenceDictionary::from(header));

        self.sample = derive_sample(&self.input_path, header);
        self.plot_title = format!("GC Bias of {}", self.sample);
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        self.process_record(record)
    }

    fn finish(&mut self) -> Result<()> {
        // Scan any BAM header contigs not visited during record traversal
        let dict = self.dict.as_ref().unwrap();
        for ref_id in 0..dict.len() {
            if !self.visited_contigs.contains(&ref_id) {
                let name = dict[ref_id].name();
                let seq = self.reference.load_contig(name, false)?;
                let (_, window_counts) = scan_contig_gc(&seq, self.window_size);
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
    /// Total clusters (first-of-pair or unpaired mapped reads).
    pub total_clusters: u64,
    /// Total aligned reads counted.
    pub aligned_reads: u64,
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
fn scan_contig_gc(seq: &[u8], window_size: u32) -> (Vec<u8>, [u64; NUM_GC_BINS]) {
    let ws = window_size as usize;
    let mut gc_at_pos = vec![u8::MAX; seq.len()];
    let mut window_counts = [0u64; NUM_GC_BINS];

    if seq.len() < ws {
        return (gc_at_pos, window_counts);
    }

    // Initialize counts for the first window using the lookup table
    let mut gc_count: u32 = 0;
    let mut n_count: u32 = 0;
    for &b in &seq[..ws] {
        let flags = BASE_FLAGS[b as usize];
        gc_count += u32::from(flags & GC_FLAG);
        n_count += u32::from((flags & N_FLAG) >> 1);
    }

    // Record the first window
    if n_count <= MAX_N_IN_WINDOW {
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

        if n_count <= MAX_N_IN_WINDOW {
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
        let (gc_u, counts_u) = scan_contig_gc(upper, 4);
        let (gc_l, counts_l) = scan_contig_gc(lower, 4);
        let (gc_m, counts_m) = scan_contig_gc(mixed, 4);
        assert_eq!(gc_u, gc_l);
        assert_eq!(gc_u, gc_m);
        assert_eq!(counts_u, counts_l);
        assert_eq!(counts_u, counts_m);
    }

    #[test]
    fn test_scan_contig_gc_simple() {
        // All-A sequence with window_size=4 → GC% = 0 at every position
        let seq = b"AAAAAAAA";
        let (gc, counts) = scan_contig_gc(seq, 4);
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
        let (gc, counts) = scan_contig_gc(seq, 4);
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
        let (gc, _) = scan_contig_gc(seq, 4);
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
        let (gc, counts) = scan_contig_gc(seq, 5);
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
        let (gc, counts) = scan_contig_gc(seq, 4);
        assert_eq!(gc.len(), 2);
        assert_eq!(gc[0], u8::MAX);
        assert_eq!(gc[1], u8::MAX);
        assert_eq!(counts.iter().sum::<u64>(), 0);
    }

    #[test]
    fn test_dropout_calculation() {
        // Manually test dropout formula
        let total_windows = 100u64;
        let total_reads = 100u64;
        // Bin 25% (AT region): 40 windows, 20 reads
        // relative_windows = 0.40, relative_reads = 0.20
        // dropout = (0.40 - 0.20) * 100 = 20.0 (positive → AT dropout)
        let rel_w = 40.0 / total_windows as f64;
        let rel_r = 20.0 / total_reads as f64;
        let dropout = (rel_w - rel_r) * 100.0;
        assert!((dropout - 20.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_phred_quality_calculation() {
        // errors=1, bases=1000 → -10*log10(0.001) = 30.0
        let errors = 1.0_f64;
        let bases = 1000.0_f64;
        let phred = -10.0 * (errors / bases).log10();
        assert!((phred - 30.0).abs() < 0.001);
    }

    // ── scan_contig_gc — additional edge cases ──────────────────────────────

    #[test]
    fn test_scan_contig_gc_window_size_1() {
        // Window size 1: each position is its own window
        let seq = b"GCAT";
        let (gc, counts) = scan_contig_gc(seq, 1);
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
        let (gc, counts) = scan_contig_gc(seq, 4);
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
        let (gc, _) = scan_contig_gc(seq, 5);
        // pos 0: NNNNA → 4 Ns, valid
        assert_ne!(gc[0], u8::MAX);
        assert_eq!(gc[0], 0); // 0 GC bases out of 5
        // pos 4: ANNNN → 4 Ns, valid
        assert_ne!(gc[4], u8::MAX);
        assert_eq!(gc[4], 0);
        // But "NNNNN" (5 Ns) would be invalid
        let seq2 = b"NNNNN";
        let (gc2, _) = scan_contig_gc(seq2, 5);
        assert_eq!(gc2[0], u8::MAX);
    }

    // ── normalized_coverage_and_error ────────────────────────────────────────

    #[test]
    fn test_nc_normal() {
        let (nc, eb) = normalized_coverage_and_error(100, 50, 1.0);
        assert!((nc - 2.0).abs() < 1e-9);
        assert!((eb - 0.2).abs() < 1e-9); // sqrt(100)/50 / 1.0 = 10/50 = 0.2
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
