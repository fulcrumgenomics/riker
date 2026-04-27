use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use kuva::plot::legend::LegendPosition;
use kuva::plot::{Histogram, LinePlot};
use kuva::render::annotations::{Orientation, ShadedRegion};
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::metrics::write_tsv;
use crate::plotting::{
    FG_BLUE, FG_GREEN, FG_PACIFIC, FG_RED, FG_TEAL, PLOT_HEIGHT, PLOT_WIDTH, write_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};

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
///
/// Only bases with a corresponding quality score are counted. Well-formed
/// BAMs always have `len(SEQ) == len(QUAL)` (the BAM writer enforces this
/// at serialization time); on malformed records where the quality string
/// is shorter than the sequence, the trailing bases are skipped.
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
        let mut reader =
            AlignmentReader::open(&self.input.input, self.reference.reference.as_deref())?;
        let mut collector = BasicCollector::new(&self.input.input, &self.output.output);
        let mut progress = ProgressLogger::new("basic", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

// ─── Collector ───────────────────────────────────────────────────────────────

/// Mask used to index into `CycleStats::base_counts`. Taking the low 5 bits
/// of an ASCII letter is case-insensitive (bit 5 is the case bit) and maps
/// each of A/C/G/T/N to a distinct slot, with no collisions among the IUPAC
/// ambiguity codes (W/S/M/K/R/Y/B/D/H/V) either.
const BASE_BITS: u8 = 0x1F;
/// Number of slots in `CycleStats::base_counts`; covers all values in `0..32`.
const N_BASE_SLOTS: usize = 32;

const IDX_A: usize = (b'A' & BASE_BITS) as usize;
const IDX_C: usize = (b'C' & BASE_BITS) as usize;
const IDX_G: usize = (b'G' & BASE_BITS) as usize;
const IDX_T: usize = (b'T' & BASE_BITS) as usize;

/// Bitmask of `(base & 0x1F)` slots that correspond to canonical A/C/G/T
/// (case-insensitive). Used to gate the global quality-distribution
/// histogram so that only ACGT bases — *not* N, lowercase n, or any of the
/// IUPAC ambiguity codes (W/S/M/K/R/Y/B/D/H/V) — contribute to it. Picard
/// excludes these too; the `(ACGT_BITMASK >> bi) & 1` check restores the
/// case-insensitive symmetry between `base_counts` (which folds case via
/// `& 0x1F`) and `qual_counts` (which previously checked literal `b'N'`).
const ACGT_BITMASK: u32 = (1 << IDX_A) | (1 << IDX_C) | (1 << IDX_G) | (1 << IDX_T);

/// Per-cycle accumulators. The hot loop touches exactly one `CycleStats`
/// per base: one increment of `qual_sum` plus one increment of the slot in
/// `base_counts` selected by `(base & 0x1F)`. Both the per-cycle total and
/// the per-cycle quality count (which are equal — every counted base has
/// an associated quality observation, since `accept` only walks the
/// `min(len(SEQ), len(QUAL))` prefix) are recovered at finalize time as
/// `base_counts.iter().sum()` — cheap (32 adds per cycle, paid once).
#[derive(Clone, Default, Debug)]
struct CycleStats {
    /// Sum of quality scores observed at this cycle (for mean-quality calc).
    /// Stored as `u64` so the inner loop avoids an int→float convert and a
    /// floating-point add per base. Headroom is ample: at the worst-case
    /// `q = 255` it would take ~7.2×10^16 reads at one cycle to overflow.
    qual_sum: u64,
    /// Per-base counts indexed by `(base & 0x1F) as usize`. Use
    /// `IDX_A/C/G/T` to look up; everything else is folded into the
    /// "other / N" bucket via `total - (a + c + g + t)` at finalize time.
    base_counts: [u64; N_BASE_SLOTS],
}

impl CycleStats {
    /// Total number of bases observed at this cycle. Equivalent to
    /// `base_counts.iter().sum::<u64>()` — recomputed on demand so the hot
    /// loop doesn't pay a per-base store to maintain it.
    fn total(&self) -> u64 {
        self.base_counts.iter().sum()
    }
}

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

    // Per read-end, per cycle.
    r1_cycles: Vec<CycleStats>,
    r2_cycles: Vec<CycleStats>,

    /// Quality score distribution. Four interleaved sub-histograms
    /// (indexed by `i & 3` where `i` is the base index within a read) so
    /// that consecutive iterations of `accept`'s hot loop write to
    /// different banks. This breaks the read-after-write dependency on
    /// `qual_counts[qi]` that otherwise serializes at load-store latency
    /// when consecutive bases share the same quality score (very common
    /// in Illumina data). Summed across banks once in
    /// `build_qual_dist_metrics`.
    qual_counts: [[u64; 128]; 4],
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
            r1_cycles: Vec::new(),
            r2_cycles: Vec::new(),
            qual_counts: [[0u64; 128]; 4],
        }
    }

    /// Ensure the per-cycle vector is long enough for `len` cycles.
    fn ensure_capacity(cycles: &mut Vec<CycleStats>, len: usize) {
        if len > cycles.len() {
            cycles.resize(len, CycleStats::default());
        }
    }

    /// Per-record inner loop, monomorphized on traversal direction so the
    /// `is_reverse` choice happens once at dispatch time rather than once per
    /// base. With `REVERSE = false` the cycle index is `i` (forward); with
    /// `REVERSE = true` it is `n - 1 - i` (reverse-complemented reads in the
    /// BAM are stored in the opposite orientation from sequencing).
    ///
    /// The caller guarantees `seq.len() == quals.len()` (call site
    /// pre-slices both to `n = min(seq.len(), quals.len())`) and that
    /// `cycles` covers at least `seq.len()` entries; the function reslices
    /// `cycles[..n]` internally so the compiler has a clean bound for
    /// `cycles[cycle_idx]`.
    #[inline]
    fn process_record<const REVERSE: bool>(
        seq: &[u8],
        quals: &[u8],
        cycles: &mut [CycleStats],
        qual_counts: &mut [[u64; 128]; 4],
    ) {
        let n = seq.len();
        let cycles = &mut cycles[..n];
        for i in 0..n {
            let base = seq[i];
            let q = quals[i];
            let cycle_idx = if REVERSE { n - 1 - i } else { i };
            // The mask proves `bi < 32`, so the indexed store skips bounds
            // checks entirely.
            let bi = (base & BASE_BITS) as usize;

            let stats = &mut cycles[cycle_idx];
            stats.base_counts[bi] += 1;
            stats.qual_sum += u64::from(q);

            // Quality distribution: only ACGT (case-insensitive) bases
            // contribute. Excludes N/n and all IUPAC ambiguity codes —
            // matching Picard's QualityScoreDistribution semantics. The
            // bitmask test `(ACGT_BITMASK >> bi) & 1` is two ops vs the
            // four compares a full `bi == IDX_A || ...` check would emit.
            // Round-robin across 4 banks via `i & 3` to break the RAW
            // dependency on `qual_counts[qi]` between consecutive
            // iterations. The `q & 0x7F` mask folds the (impossible-in-
            // valid-BAM) values q ≥ 128 into bins 0..=127; a missing-
            // quality record encoded as 0xFF would land entirely in bin 127.
            if (ACGT_BITMASK >> bi) & 1 != 0 {
                let qi = (q & 0x7F) as usize;
                qual_counts[i & 3][qi] += 1;
            }
        }
    }

    /// Build base distribution metrics from the accumulated counts. Counts
    /// for A/C/G/T are read by their `(base & 0x1F)` slot; everything else
    /// (N and any non-ACGT byte) is folded into `frac_n` as a residual.
    fn build_base_dist_metrics(&self) -> Vec<BaseDistributionByCycleMetric> {
        let mut metrics = Vec::new();

        for (read_end, cycles) in [(1u8, &self.r1_cycles), (2u8, &self.r2_cycles)] {
            for (i, c) in cycles.iter().enumerate() {
                let total = c.total();
                if total == 0 {
                    continue;
                }
                let t = total as f64;
                let a = c.base_counts[IDX_A];
                let cnt_c = c.base_counts[IDX_C];
                let g = c.base_counts[IDX_G];
                let cnt_t = c.base_counts[IDX_T];
                // `total = base_counts.iter().sum()` is the sum of all 32
                // slots; the four named indices are pairwise disjoint, so
                // `a + cnt_c + g + cnt_t` is at most `total`. The
                // debug-only assertion documents this invariant.
                debug_assert!(a + cnt_c + g + cnt_t <= total);
                let n = total - (a + cnt_c + g + cnt_t);
                metrics.push(BaseDistributionByCycleMetric {
                    read_end,
                    cycle: i + 1,
                    frac_a: a as f64 / t,
                    frac_c: cnt_c as f64 / t,
                    frac_g: g as f64 / t,
                    frac_t: cnt_t as f64 / t,
                    frac_n: n as f64 / t,
                });
            }
        }

        metrics
    }

    /// Build mean quality by cycle metrics. R2 cycles are offset by max R1 cycle count.
    fn build_mean_qual_metrics(&self) -> Vec<MeanQualityByCycleMetric> {
        let mut metrics = Vec::new();
        let r1_max = self.r1_cycles.len();

        for (i, c) in self.r1_cycles.iter().enumerate() {
            let total = c.total();
            if total == 0 {
                continue;
            }
            metrics.push(MeanQualityByCycleMetric {
                cycle: i + 1,
                mean_quality: c.qual_sum as f64 / total as f64,
            });
        }

        for (i, c) in self.r2_cycles.iter().enumerate() {
            let total = c.total();
            if total == 0 {
                continue;
            }
            metrics.push(MeanQualityByCycleMetric {
                cycle: r1_max + i + 1,
                mean_quality: c.qual_sum as f64 / total as f64,
            });
        }

        metrics
    }

    /// Build quality score distribution metrics (only non-zero counts).
    /// Collapses the four interleaved sub-histograms maintained by `accept`
    /// into a single 128-bin distribution before emitting metrics.
    #[allow(clippy::cast_possible_truncation)] // quality scores are always 0..127
    fn build_qual_dist_metrics(&self) -> Vec<QualityScoreDistributionMetric> {
        let mut combined = [0u64; 128];
        for bank in &self.qual_counts {
            for (slot, &c) in bank.iter().enumerate() {
                combined[slot] += c;
            }
        }
        let total: u64 = combined.iter().sum();
        let total_f = total as f64;
        combined
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

        let r1_max = self.r1_cycles.len();
        let max_cycle = (r1_max + self.r2_cycles.len()) as f64;
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

        if !self.r2_cycles.is_empty() {
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

        let r1_max = self.r1_cycles.len();
        let has_r2 = !self.r2_cycles.is_empty();

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

        let max_cycle = (r1_max + self.r2_cycles.len()) as f64;
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
}

impl Collector for BasicCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        let label = derive_sample(&self.input, header);
        self.plot_title_prefix = label;
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // Skip secondary, supplementary, and QC-fail reads
        if flags.is_secondary() || flags.is_supplementary() || flags.is_qc_fail() {
            return Ok(());
        }

        let seq: &[u8] = record.sequence();
        let quals: &[u8] = record.quality_scores();
        // Hoist the seq/quals length reconciliation out of the per-base
        // loop. Well-formed BAMs always have `seq.len() == quals.len()`;
        // the BAM writer enforces this at serialization time. Truncated /
        // missing quality data falls out of `n` and contributes neither
        // base counts nor quality stats.
        let n = seq.len().min(quals.len());
        if n == 0 {
            return Ok(());
        }

        let is_reverse = flags.is_reverse_complemented();
        // Determine R1 vs R2: paired + last_segment = R2, else R1
        let is_r2 = flags.is_segmented() && flags.is_last_segment();

        let cycles = if is_r2 { &mut self.r2_cycles } else { &mut self.r1_cycles };
        Self::ensure_capacity(cycles, n);

        // Dispatch on is_reverse once so the inner loop is a straight-line
        // monomorphization for each direction (no per-base branch).
        let seq = &seq[..n];
        let quals = &quals[..n];
        if is_reverse {
            Self::process_record::<true>(seq, quals, cycles, &mut self.qual_counts);
        } else {
            Self::process_record::<false>(seq, quals, cycles, &mut self.qual_counts);
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

    fn field_needs(&self) -> RikerRecordRequirements {
        RikerRecordRequirements::NONE.with_sequence()
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
    fn test_base_indices_distinct() {
        let idxs = [IDX_A, IDX_C, IDX_G, IDX_T];
        for (i, a) in idxs.iter().enumerate() {
            for b in &idxs[i + 1..] {
                assert_ne!(a, b, "base index collision: {a} == {b}");
            }
            assert!(*a < N_BASE_SLOTS);
        }
    }

    #[test]
    fn test_base_mask_is_case_insensitive() {
        for (upper, lower) in [(b'A', b'a'), (b'C', b'c'), (b'G', b'g'), (b'T', b't')] {
            assert_eq!(upper & BASE_BITS, lower & BASE_BITS);
        }
    }

    /// `(base & 0x1F)` must produce a distinct slot for every BAM-encodable
    /// nucleotide character so that none of N, lowercase n, or any IUPAC
    /// ambiguity code silently collides with A/C/G/T (which would corrupt
    /// `frac_n` by inflating one of the canonical-base counts). Locks in
    /// the load-bearing claim from `BASE_BITS`'s doc comment.
    #[test]
    fn test_base_mask_no_collisions_with_iupac_or_n() {
        let canonical = [(b'A', IDX_A), (b'C', IDX_C), (b'G', IDX_G), (b'T', IDX_T)];
        // BAM 4-bit decode emits these characters: =ACMGRSVTWYHKDBN (no
        // lowercase). All non-ACGT entries must map to a slot distinct
        // from the four canonical-base slots above.
        let non_acgt = b"=MRSVWYHKDBN";
        for &b in non_acgt {
            let slot = (b & BASE_BITS) as usize;
            for &(_, canon) in &canonical {
                assert_ne!(
                    slot, canon,
                    "non-ACGT byte 0x{b:02x} ({}) collides with canonical slot {canon}",
                    b as char
                );
            }
        }
        // Lowercase n must fold to the same slot as uppercase N (case bit).
        assert_eq!(b'N' & BASE_BITS, b'n' & BASE_BITS);
    }

    /// `ACGT_BITMASK` is the bitset over `(base & 0x1F)` slots that
    /// represent canonical A/C/G/T. The hot loop tests it to gate
    /// `qual_counts` updates, so collisions or omissions corrupt the
    /// quality distribution. Lock down membership for the cases that
    /// matter.
    #[test]
    fn test_acgt_bitmask_membership() {
        for &b in b"ACGTacgt" {
            let bi = (b & BASE_BITS) as usize;
            assert_ne!(
                (ACGT_BITMASK >> bi) & 1,
                0,
                "byte 0x{b:02x} ({}) should be in ACGT_BITMASK",
                b as char
            );
        }
        for &b in b"NnWSMKRYBDHV=" {
            let bi = (b & BASE_BITS) as usize;
            assert_eq!(
                (ACGT_BITMASK >> bi) & 1,
                0,
                "non-ACGT byte 0x{b:02x} ({}) should not be in ACGT_BITMASK",
                b as char
            );
        }
    }
}
