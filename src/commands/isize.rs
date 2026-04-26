use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use kuva::plot::LinePlot;
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::counter::Counter;
use crate::metrics::write_tsv;
use crate::plotting::{FG_BLUE, FG_GREEN, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_plot_pdf};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::pair_orientation::{PairOrientation, get_pair_orientation};
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};

/// File suffix appended to the output prefix for the metrics file.
pub const METRICS_SUFFIX: &str = ".isize-metrics.txt";

/// File suffix appended to the output prefix for the histogram file.
pub const HISTOGRAM_SUFFIX: &str = ".isize-histogram.txt";

/// File suffix appended to the output prefix for the histogram plot.
pub const PLOT_SUFFIX: &str = ".isize-histogram.pdf";

// ─── Options ─────────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the insert size collector.
#[riker_derive::multi_options("isize", "Insert Size Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct IsizeOptions {
    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,
    /// Minimum fraction of total read pairs required for an orientation to be reported.
    #[arg(long, default_value_t = IsizeOptions::DEFAULT_MIN_FRAC)]
    pub min_frac: f64,
    /// Trim to MEDIAN ± DEVIATIONS × MAD before computing mean and standard deviation.
    #[arg(long, default_value_t = IsizeOptions::DEFAULT_DEVIATIONS)]
    pub deviations: f64,
}

impl IsizeOptions {
    const DEFAULT_MIN_FRAC: f64 = 0.05;
    const DEFAULT_DEVIATIONS: f64 = 10.0;
}

impl Default for IsizeOptions {
    fn default() -> Self {
        Self {
            include_duplicates: false,
            min_frac: Self::DEFAULT_MIN_FRAC,
            deviations: Self::DEFAULT_DEVIATIONS,
        }
    }
}

/// Collect insert size distribution metrics from a paired-end BAM file.
///
/// Computes mean, median, mode, standard deviation, and MAD of insert
/// sizes for each pair orientation (FR, RF, tandem). Produces a summary
/// metrics file, a per-size histogram, and a distribution chart. Outputs
/// are written to <prefix>.isize-metrics.txt,
/// <prefix>.isize-histogram.txt, and <prefix>.isize-histogram.pdf.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker isize -i input.bam -o out_prefix
  riker isize -i input.bam -o out_prefix --include-duplicates"
)]
pub struct InsertSize {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,

    #[command(flatten)]
    pub options: IsizeOptions,
}

impl Command for InsertSize {
    /// # Errors
    /// Returns an error if the BAM file cannot be read or the output file cannot be written.
    fn execute(&self) -> Result<()> {
        let mut reader =
            AlignmentReader::open(&self.input.input, self.reference.reference.as_deref())?;
        let mut collector =
            InsertSizeCollector::new(&self.input.input, &self.output.output, &self.options);
        let mut progress = ProgressLogger::new("isize", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

// ─── Collector ───────────────────────────────────────────────────────────────

/// Collects insert size distributions from a paired-end BAM.
pub struct InsertSizeCollector {
    input: PathBuf,
    metrics_path: PathBuf,
    histogram_path: PathBuf,
    plot_path: PathBuf,
    plot_title: String,
    include_duplicates: bool,
    min_frac: f64,
    deviations: f64,
    fr: Counter<u64>,
    rf: Counter<u64>,
    tandem: Counter<u64>,
}

impl InsertSizeCollector {
    /// Create a new collector. Output paths are derived from `prefix` by appending
    /// [`METRICS_SUFFIX`], [`HISTOGRAM_SUFFIX`], and [`PLOT_SUFFIX`].  The plot title is
    /// populated during [`Collector::initialize`] once the BAM header is available.
    #[must_use]
    pub fn new(input: &std::path::Path, prefix: &std::path::Path, options: &IsizeOptions) -> Self {
        let metrics_path = super::command::output_path(prefix, METRICS_SUFFIX);
        let histogram_path = super::command::output_path(prefix, HISTOGRAM_SUFFIX);
        let plot_path = super::command::output_path(prefix, PLOT_SUFFIX);
        Self {
            input: input.to_path_buf(),
            metrics_path,
            histogram_path,
            plot_path,
            plot_title: String::new(),
            include_duplicates: options.include_duplicates,
            min_frac: options.min_frac,
            deviations: options.deviations,
            fr: Counter::new(),
            rf: Counter::new(),
            tandem: Counter::new(),
        }
    }

    /// Write a PDF area-chart of the insert size distributions to [`Self::plot_path`].
    ///
    /// Reads directly from the instance counters, applies the same minimum-fraction and
    /// outlier-trim thresholds used for metrics, and assigns FG brand colours per
    /// orientation.  A legend is shown when more than one orientation is plotted.
    /// Returns immediately (without writing a file) when no orientations qualify.
    fn plot_histogram(&self) -> Result<()> {
        let total = self.fr.total() + self.rf.total() + self.tandem.total();

        let orientations: [(&str, &str, &Counter<u64>); 3] = [
            ("FR", FG_BLUE, &self.fr),
            ("RF", FG_GREEN, &self.rf),
            ("TANDEM", FG_TEAL, &self.tandem),
        ];

        let mut plots: Vec<Plot> = Vec::new();
        for (name, color, hist) in &orientations {
            let count = hist.total();
            if count == 0 || (total > 0 && below_min_frac(count, total, self.min_frac)) {
                continue;
            }
            let (median, mad) = hist.median_and_mad();
            let trim_max = median + self.deviations * mad;
            let xy: Vec<(f64, f64)> = hist
                .sorted()
                .into_iter()
                .take_while(|&(k, _)| k as f64 <= trim_max)
                .map(|(x, y)| (x as f64, y as f64))
                .collect();
            plots.push(Plot::Line(
                LinePlot::new()
                    .with_data(xy)
                    .with_color(*color)
                    .with_fill()
                    .with_fill_opacity(0.3)
                    .with_legend(*name),
            ));
        }

        if plots.is_empty() {
            return Ok(());
        }

        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(&self.plot_title)
            .with_x_label("Insert Size (bp)")
            .with_y_label("Read Pairs")
            .with_minor_ticks(5)
            .with_show_minor_grid(true);

        write_plot_pdf(plots, layout, &self.plot_path)
    }
}

impl Collector for InsertSizeCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        let label = derive_sample(&self.input, header);
        self.plot_title = format!("Insert Size Distribution of {label}");
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // Filtering logic mirrors Picard's InsertSizeMetricsCollector.acceptRecord:
        // process second-of-pair only to avoid double-counting.
        if flags.is_unmapped()
            || flags.is_mate_unmapped()
            || !flags.is_segmented()
            || flags.is_first_segment()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
            || (!self.include_duplicates && flags.is_duplicate())
            || record.template_length() == 0
            || record.reference_sequence_id() != record.mate_reference_sequence_id()
        {
            return Ok(());
        }

        let insert_size = u64::from(record.template_length().unsigned_abs());
        let Some(orientation) = get_pair_orientation(record) else {
            return Ok(());
        };

        match orientation {
            PairOrientation::Fr => self.fr.count(insert_size),
            PairOrientation::Rf => self.rf.count(insert_size),
            PairOrientation::Tandem => self.tandem.count(insert_size),
        }

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        let total: u64 = self.fr.total() + self.rf.total() + self.tandem.total();

        if total == 0 {
            log::warn!("isize: no qualifying read pairs found; writing empty metrics file");
        }

        let orientations: [(&str, &Counter<u64>); 3] =
            [("FR", &self.fr), ("RF", &self.rf), ("TANDEM", &self.tandem)];

        let mut metrics: Vec<InsertSizeMetric> = Vec::new();

        for (name, hist) in &orientations {
            let count = hist.total();
            if count == 0 {
                continue;
            }
            if total > 0 && below_min_frac(count, total, self.min_frac) {
                continue;
            }

            let (median, mad) = hist.median_and_mad();
            let mode = hist.mode().unwrap_or(0);
            let min_is = hist.min().expect("non-empty histogram");
            let max_is = hist.max().expect("non-empty histogram");
            let trim_max = median + self.deviations * mad;
            let (mean, stddev) = hist.trimmed_mean_and_stddev(trim_max);

            metrics.push(InsertSizeMetric {
                pair_orientation: (*name).to_string(),
                read_pairs: count,
                mean_insert_size: mean,
                standard_deviation: stddev,
                median_insert_size: median,
                median_absolute_deviation: mad,
                mode_insert_size: mode,
                min_insert_size: min_is,
                max_insert_size: max_is,
            });
        }

        write_tsv(&self.metrics_path, &metrics)?;

        let all_keys: std::collections::BTreeSet<u64> =
            self.fr.keys().chain(self.rf.keys()).chain(self.tandem.keys()).copied().collect();

        let hist_entries: Vec<InsertSizeHistogramEntry> = all_keys
            .into_iter()
            .map(|k| InsertSizeHistogramEntry {
                insert_size: k,
                fr_count: self.fr.count_of(&k),
                rf_count: self.rf.count_of(&k),
                tandem_count: self.tandem.count_of(&k),
            })
            .collect();

        write_tsv(&self.histogram_path, &hist_entries)?;
        self.plot_histogram()?;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "isize"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        RikerRecordRequirements::NONE
    }
}

// ─── Metric structs ──────────────────────────────────────────────────────────

/// Insert size distribution metrics, one row per pair orientation.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct InsertSizeMetric {
    /// Pair orientation (FR, RF, or TANDEM).
    pub pair_orientation: String,
    /// Number of read pairs in this orientation.
    pub read_pairs: u64,
    /// Mean insert size after trimming outliers.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_insert_size: f64,
    /// Standard deviation of insert size after trimming outliers.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub standard_deviation: f64,
    /// Median insert size across all pairs in this orientation.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub median_insert_size: f64,
    /// Median absolute deviation of insert size.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub median_absolute_deviation: f64,
    /// Most frequent insert size.
    pub mode_insert_size: u64,
    /// Smallest insert size observed.
    pub min_insert_size: u64,
    /// Largest insert size observed.
    pub max_insert_size: u64,
}

/// Insert size histogram with counts per orientation at each size.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct InsertSizeHistogramEntry {
    /// Insert size in base pairs.
    pub insert_size: u64,
    /// Number of FR (forward-reverse) pairs at this insert size.
    pub fr_count: u64,
    /// Number of RF (reverse-forward) pairs at this insert size.
    pub rf_count: u64,
    /// Number of tandem (FF/RR) pairs at this insert size.
    pub tandem_count: u64,
}

// ─── Statistics helpers ───────────────────────────────────────────────────────

/// Returns true when `count` is below the minimum-fraction threshold.
pub(crate) fn below_min_frac(count: u64, total: u64, min_frac: f64) -> bool {
    (count as f64) < (total as f64) * min_frac
}

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    // ── below_min_frac ──────────────────────────────────────────────────────

    #[test]
    fn test_below_min_frac_below() {
        // 4% < 5% → true
        assert!(below_min_frac(4, 100, 0.05));
    }

    #[test]
    fn test_below_min_frac_at_threshold() {
        // 5% == 5% → false (not below)
        assert!(!below_min_frac(5, 100, 0.05));
    }

    #[test]
    fn test_below_min_frac_above() {
        // 10% > 5% → false
        assert!(!below_min_frac(10, 100, 0.05));
    }

    #[test]
    fn test_below_min_frac_zero_total() {
        // 0 < 0*0.05 = 0.0 → false (not strictly below)
        assert!(!below_min_frac(0, 0, 0.05));
    }
}
