use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use kuva::plot::legend::LegendPosition;
use kuva::plot::{Histogram, LinePlot};
use kuva::render::annotations::{Orientation, ShadedRegion};
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::Collector;
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::metrics::write_tsv;
use crate::plotting::{
    FG_BLUE, FG_GREEN, FG_PACIFIC, FG_RED, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::derive_sample;

// ─── Output file suffixes ────────────────────────────────────────────────────

/// File suffix for base distribution by cycle metrics.
pub const BASE_DIST_SUFFIX: &str = ".base-distribution-by-cycle.txt";
/// File suffix for mean quality by cycle metrics.
pub const MEAN_QUAL_SUFFIX: &str = ".mean-quality-by-cycle.txt";
/// File suffix for quality score distribution metrics.
pub const QUAL_DIST_SUFFIX: &str = ".quality-score-distribution.txt";
/// File suffix for base distribution by cycle plot.
pub const BASE_DIST_PLOT_SUFFIX: &str = ".base-distribution-by-cycle.pdf";
/// File suffix for mean quality by cycle plot.
pub const MEAN_QUAL_PLOT_SUFFIX: &str = ".mean-quality-by-cycle.pdf";
/// File suffix for quality score distribution plot.
pub const QUAL_DIST_PLOT_SUFFIX: &str = ".quality-score-distribution.pdf";

// ─── Command ─────────────────────────────────────────────────────────────────

/// Collect basic QC metrics: base composition by cycle, mean quality by
/// cycle, and quality score distribution.
///
/// Combines Picard's CollectBaseDistributionByCycle,
/// MeanQualityByCycle, and QualityScoreDistribution into a single BAM
/// pass. Duplicates are included by default (matching Picard). Outputs
/// six files: <prefix>.base-distribution-by-cycle.txt,
/// <prefix>.base-distribution-by-cycle.pdf,
/// <prefix>.mean-quality-by-cycle.txt,
/// <prefix>.mean-quality-by-cycle.pdf,
/// <prefix>.quality-score-distribution.txt, and
/// <prefix>.quality-score-distribution.pdf.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker basic -i input.bam -o out_prefix"
)]
pub struct Basic {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,
}

impl Command for Basic {
    /// # Errors
    /// Returns an error if the BAM file cannot be read or the output files cannot be written.
    fn execute(&self) -> Result<()> {
        let (mut reader, header) =
            AlignmentReader::new(&self.input.input, self.reference.reference.as_deref())?;
        let mut collector = BasicCollector::new(&self.input.input, &self.output.output);
        collector.initialize(&header)?;
        let mut progress = ProgressLogger::new("basic", "reads", 5_000_000);
        reader.for_each_record(&header, |record| {
            progress.record_with(record, &header);
            collector.accept(record, &header)
        })?;
        progress.finish();
        collector.finish()
    }
}

// ─── Collector ───────────────────────────────────────────────────────────────

/// Collects base distribution, mean quality, and quality score distribution
/// from a BAM file in a single pass.
pub struct BasicCollector {
    input: PathBuf,
    base_dist_path: PathBuf,
    mean_qual_path: PathBuf,
    qual_dist_path: PathBuf,
    base_dist_plot_path: PathBuf,
    mean_qual_plot_path: PathBuf,
    qual_dist_plot_path: PathBuf,
    plot_title_prefix: String,

    // Per read-end, per base (A=0,C=1,G=2,T=3,N=4), per cycle (0-indexed)
    r1_base_counts: [Vec<u64>; 5],
    r1_cycle_totals: Vec<u64>,
    r2_base_counts: [Vec<u64>; 5],
    r2_cycle_totals: Vec<u64>,

    // Per read-end, per cycle — quality sums and counts
    r1_qual_sums: Vec<f64>,
    r1_qual_counts: Vec<u64>,
    r2_qual_sums: Vec<f64>,
    r2_qual_counts: Vec<u64>,

    // Quality score distribution (index = quality score)
    qual_counts: [u64; 128],
}

impl BasicCollector {
    /// Create a new collector. Output paths are derived from `prefix`.
    #[must_use]
    pub fn new(input: &std::path::Path, prefix: &std::path::Path) -> Self {
        Self {
            input: input.to_path_buf(),
            base_dist_path: super::command::output_path(prefix, BASE_DIST_SUFFIX),
            mean_qual_path: super::command::output_path(prefix, MEAN_QUAL_SUFFIX),
            qual_dist_path: super::command::output_path(prefix, QUAL_DIST_SUFFIX),
            base_dist_plot_path: super::command::output_path(prefix, BASE_DIST_PLOT_SUFFIX),
            mean_qual_plot_path: super::command::output_path(prefix, MEAN_QUAL_PLOT_SUFFIX),
            qual_dist_plot_path: super::command::output_path(prefix, QUAL_DIST_PLOT_SUFFIX),
            plot_title_prefix: String::new(),
            r1_base_counts: [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            r1_cycle_totals: Vec::new(),
            r2_base_counts: [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            r2_cycle_totals: Vec::new(),
            r1_qual_sums: Vec::new(),
            r1_qual_counts: Vec::new(),
            r2_qual_sums: Vec::new(),
            r2_qual_counts: Vec::new(),
            qual_counts: [0u64; 128],
        }
    }

    /// Ensure all per-cycle vectors are long enough for `len` cycles.
    fn ensure_capacity(
        base_counts: &mut [Vec<u64>; 5],
        cycle_totals: &mut Vec<u64>,
        qual_sums: &mut Vec<f64>,
        qual_counts: &mut Vec<u64>,
        len: usize,
    ) {
        if len > cycle_totals.len() {
            for bc in base_counts.iter_mut() {
                bc.resize(len, 0);
            }
            cycle_totals.resize(len, 0);
            qual_sums.resize(len, 0.0);
            qual_counts.resize(len, 0);
        }
    }

    /// Build base distribution metrics from the accumulated counts.
    fn build_base_dist_metrics(&self) -> Vec<BaseDistributionByCycleMetric> {
        let mut metrics = Vec::new();

        // R1 cycles
        for i in 0..self.r1_cycle_totals.len() {
            let total = self.r1_cycle_totals[i];
            if total == 0 {
                continue;
            }
            let t = total as f64;
            metrics.push(BaseDistributionByCycleMetric {
                read_end: 1,
                cycle: i + 1,
                frac_a: self.r1_base_counts[0][i] as f64 / t,
                frac_c: self.r1_base_counts[1][i] as f64 / t,
                frac_g: self.r1_base_counts[2][i] as f64 / t,
                frac_t: self.r1_base_counts[3][i] as f64 / t,
                frac_n: self.r1_base_counts[4][i] as f64 / t,
            });
        }

        // R2 cycles
        for i in 0..self.r2_cycle_totals.len() {
            let total = self.r2_cycle_totals[i];
            if total == 0 {
                continue;
            }
            let t = total as f64;
            metrics.push(BaseDistributionByCycleMetric {
                read_end: 2,
                cycle: i + 1,
                frac_a: self.r2_base_counts[0][i] as f64 / t,
                frac_c: self.r2_base_counts[1][i] as f64 / t,
                frac_g: self.r2_base_counts[2][i] as f64 / t,
                frac_t: self.r2_base_counts[3][i] as f64 / t,
                frac_n: self.r2_base_counts[4][i] as f64 / t,
            });
        }

        metrics
    }

    /// Build mean quality by cycle metrics. R2 cycles are offset by max R1 cycle count.
    fn build_mean_qual_metrics(&self) -> Vec<MeanQualityByCycleMetric> {
        let mut metrics = Vec::new();
        let r1_max = self.r1_qual_counts.len();

        for i in 0..r1_max {
            if self.r1_qual_counts[i] == 0 {
                continue;
            }
            metrics.push(MeanQualityByCycleMetric {
                cycle: i + 1,
                mean_quality: self.r1_qual_sums[i] / self.r1_qual_counts[i] as f64,
            });
        }

        for i in 0..self.r2_qual_counts.len() {
            if self.r2_qual_counts[i] == 0 {
                continue;
            }
            metrics.push(MeanQualityByCycleMetric {
                cycle: r1_max + i + 1,
                mean_quality: self.r2_qual_sums[i] / self.r2_qual_counts[i] as f64,
            });
        }

        metrics
    }

    /// Build quality score distribution metrics (only non-zero counts).
    #[allow(clippy::cast_possible_truncation)] // quality scores are always 0..127
    fn build_qual_dist_metrics(&self) -> Vec<QualityScoreDistributionMetric> {
        let total: u64 = self.qual_counts.iter().sum();
        let total_f = total as f64;
        self.qual_counts
            .iter()
            .enumerate()
            .filter(|&(_, count)| *count > 0)
            .map(|(q, count)| QualityScoreDistributionMetric {
                quality: q as u8,
                count: *count,
                frac_bases: if total > 0 { *count as f64 / total_f } else { 0.0 },
            })
            .collect()
    }

    /// Plot base distribution by cycle.
    fn plot_base_distribution(&self, metrics: &[BaseDistributionByCycleMetric]) -> Result<()> {
        if metrics.is_empty() {
            return Ok(());
        }

        let r1_max = self.r1_cycle_totals.len();
        let max_cycle = (r1_max + self.r2_cycle_totals.len()) as f64;
        let base_names = ["A", "C", "G", "T", "N"];
        let base_colors = [FG_BLUE, FG_GREEN, FG_TEAL, FG_PACIFIC, FG_RED];

        let plots: Vec<Plot> = (0..5)
            .map(|base_idx| {
                let xy: Vec<(f64, f64)> = metrics
                    .iter()
                    .map(|m| {
                        let x = if m.read_end == 1 {
                            m.cycle as f64
                        } else {
                            (r1_max + m.cycle) as f64
                        };
                        let y = match base_idx {
                            0 => m.frac_a,
                            1 => m.frac_c,
                            2 => m.frac_g,
                            3 => m.frac_t,
                            _ => m.frac_n,
                        };
                        (x, y)
                    })
                    .collect();

                Plot::Line(
                    LinePlot::new()
                        .with_data(xy)
                        .with_color(base_colors[base_idx])
                        .with_legend(base_names[base_idx]),
                )
            })
            .collect();

        let mut layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("{} Base Distribution by Cycle", self.plot_title_prefix))
            .with_x_label("Cycle")
            .with_y_label("Fraction")
            .with_x_axis_max(max_cycle + 1.0)
            .with_legend_position(LegendPosition::OutsideRightMiddle);

        if !self.r2_cycle_totals.is_empty() {
            layout.shaded_regions.push(ShadedRegion {
                orientation: Orientation::Vertical,
                min_val: 1.0,
                max_val: r1_max as f64 + 0.5,
                color: FG_BLUE.to_owned(),
                opacity: 0.1,
            });
            layout.shaded_regions.push(ShadedRegion {
                orientation: Orientation::Vertical,
                min_val: r1_max as f64 + 0.5,
                max_val: max_cycle + 1.0,
                color: FG_TEAL.to_owned(),
                opacity: 0.1,
            });
        }

        write_plot_pdf(plots, layout, &self.base_dist_plot_path)
    }

    /// Plot mean quality by cycle.
    fn plot_mean_quality(&self, metrics: &[MeanQualityByCycleMetric]) -> Result<()> {
        if metrics.is_empty() {
            return Ok(());
        }

        let r1_max = self.r1_qual_counts.len();
        let has_r2 = !self.r2_qual_counts.is_empty();

        let r1_xy: Vec<(f64, f64)> = metrics
            .iter()
            .filter(|m| m.cycle <= r1_max)
            .map(|m| (m.cycle as f64, m.mean_quality))
            .collect();

        let mut plots = vec![Plot::Line(
            LinePlot::new().with_data(r1_xy).with_color(FG_BLUE).with_fill().with_fill_opacity(0.3),
        )];

        if has_r2 {
            let r2_xy: Vec<(f64, f64)> = metrics
                .iter()
                .filter(|m| m.cycle > r1_max)
                .map(|m| (m.cycle as f64, m.mean_quality))
                .collect();
            plots.push(Plot::Line(
                LinePlot::new()
                    .with_data(r2_xy)
                    .with_color(FG_TEAL)
                    .with_fill()
                    .with_fill_opacity(0.3),
            ));
        }

        let max_cycle = (r1_max + self.r2_qual_counts.len()) as f64;
        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("{} Mean Quality by Cycle", self.plot_title_prefix))
            .with_x_label("Cycle")
            .with_y_label("Mean Quality")
            .with_x_axis_max(max_cycle + 1.0);

        write_plot_pdf(plots, layout, &self.mean_qual_plot_path)
    }

    /// Plot quality score distribution as a histogram with a linear X axis.
    ///
    /// Each quality score is rendered as a bar using `Histogram::from_bins`. The
    /// X axis runs from 0 to max(45, highest observed quality).
    fn plot_qual_distribution(&self, metrics: &[QualityScoreDistributionMetric]) -> Result<()> {
        if metrics.is_empty() {
            return Ok(());
        }

        // Determine the upper bound: max(45, max_quality).
        let max_quality = metrics.iter().map(|m| usize::from(m.quality)).max().unwrap_or(45);
        let upper = max_quality.max(45);

        // Build a dense lookup from quality score to frac_bases.
        let mut frac_by_q = vec![0.0_f64; upper + 1];
        for m in metrics {
            let q = usize::from(m.quality);
            if q < frac_by_q.len() {
                frac_by_q[q] = m.frac_bases;
            }
        }

        // Build bin edges [0, 1, 2, ..., upper] and corresponding counts.
        let edges: Vec<f64> = (0..=upper).map(|q| q as f64).collect();
        let counts: Vec<f64> = frac_by_q[..upper].to_vec();

        let plots = vec![Plot::Histogram(Histogram::from_bins(edges, counts).with_color(FG_BLUE))];

        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("{} Quality Score Distribution", self.plot_title_prefix))
            .with_x_label("Quality Score")
            .with_y_label("Fraction of Bases")
            .with_x_axis_min(0.0)
            .with_x_axis_max(upper as f64);

        write_plot_pdf(plots, layout, &self.qual_dist_plot_path)
    }

    /// Map a base byte to an index: A=0, C=1, G=2, T=3, N/other=4.
    fn base_index(base: u8) -> usize {
        match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4, // N and anything else
        }
    }
}

impl Collector for BasicCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        let label = derive_sample(&self.input, header);
        self.plot_title_prefix = label;
        Ok(())
    }

    fn accept(&mut self, record: &RecordBuf, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // Skip secondary, supplementary, and QC-fail reads
        if flags.is_secondary() || flags.is_supplementary() || flags.is_qc_fail() {
            return Ok(());
        }

        let seq: &[u8] = record.sequence().as_ref();
        if seq.is_empty() {
            return Ok(());
        }

        let quals: &[u8] = record.quality_scores().as_ref();
        let len = seq.len();
        let is_reverse = flags.is_reverse_complemented();

        // Determine R1 vs R2: paired + last_segment = R2, else R1
        let is_r2 = flags.is_segmented() && flags.is_last_segment();

        let (base_counts, cycle_totals, qual_sums, qual_counts) = if is_r2 {
            (
                &mut self.r2_base_counts,
                &mut self.r2_cycle_totals,
                &mut self.r2_qual_sums,
                &mut self.r2_qual_counts,
            )
        } else {
            (
                &mut self.r1_base_counts,
                &mut self.r1_cycle_totals,
                &mut self.r1_qual_sums,
                &mut self.r1_qual_counts,
            )
        };

        Self::ensure_capacity(base_counts, cycle_totals, qual_sums, qual_counts, len);

        let has_quals = !quals.is_empty();

        for i in 0..len {
            // Cycle: forward reads go 0..len, reverse reads go len-1..0
            let cycle_idx = if is_reverse { len - 1 - i } else { i };

            let base = seq[i];
            let bi = Self::base_index(base);

            base_counts[bi][cycle_idx] += 1;
            cycle_totals[cycle_idx] += 1;

            if has_quals {
                let q = quals[i];
                qual_sums[cycle_idx] += f64::from(q);
                qual_counts[cycle_idx] += 1;

                // Quality distribution: exclude N bases
                if bi != 4 {
                    let qi = usize::from(q).min(127);
                    self.qual_counts[qi] += 1;
                }
            }
        }

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        let base_dist = self.build_base_dist_metrics();
        let mean_qual = self.build_mean_qual_metrics();
        let qual_dist = self.build_qual_dist_metrics();

        write_tsv(&self.base_dist_path, &base_dist)?;
        write_tsv(&self.mean_qual_path, &mean_qual)?;
        write_tsv(&self.qual_dist_path, &qual_dist)?;

        self.plot_base_distribution(&base_dist)?;
        self.plot_mean_quality(&mean_qual)?;
        self.plot_qual_distribution(&qual_dist)?;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "basic"
    }
}

// ─── Metric structs ──────────────────────────────────────────────────────────

/// Base distribution by cycle, showing the fraction of each base at each sequencing cycle.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct BaseDistributionByCycleMetric {
    /// Read end (1 or 2).
    pub read_end: u8,
    /// Sequencing cycle number (1-based).
    pub cycle: usize,
    /// Fraction of A bases at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_a: f64,
    /// Fraction of C bases at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_c: f64,
    /// Fraction of G bases at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_g: f64,
    /// Fraction of T bases at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_t: f64,
    /// Fraction of N bases at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_n: f64,
}

/// Mean base quality by sequencing cycle.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MeanQualityByCycleMetric {
    /// Sequencing cycle number (1-based).
    pub cycle: usize,
    /// Mean base quality score at this cycle.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_quality: f64,
}

/// Quality score distribution across all bases.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct QualityScoreDistributionMetric {
    /// Base quality score.
    pub quality: u8,
    /// Number of bases with this quality score.
    pub count: u64,
    /// Fraction of all bases with this quality score.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_bases: f64,
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_index_standard_bases() {
        assert_eq!(BasicCollector::base_index(b'A'), 0);
        assert_eq!(BasicCollector::base_index(b'C'), 1);
        assert_eq!(BasicCollector::base_index(b'G'), 2);
        assert_eq!(BasicCollector::base_index(b'T'), 3);
        assert_eq!(BasicCollector::base_index(b'N'), 4);
    }

    #[test]
    fn test_base_index_lowercase() {
        assert_eq!(BasicCollector::base_index(b'a'), 0);
        assert_eq!(BasicCollector::base_index(b'c'), 1);
        assert_eq!(BasicCollector::base_index(b'g'), 2);
        assert_eq!(BasicCollector::base_index(b't'), 3);
    }

    #[test]
    fn test_base_index_unknown() {
        assert_eq!(BasicCollector::base_index(b'.'), 4);
        assert_eq!(BasicCollector::base_index(b'X'), 4);
    }
}
