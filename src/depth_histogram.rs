//! A histogram of per-position (or per-site) depths and the coverage statistics computed from
//! it: mean and standard deviation, median and MAD, fraction at or above a depth, percentiles,
//! fold penalties, and a PDF plot of the distribution.

use std::path::Path;

use anyhow::Result;
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::render::plots::Plot;

use crate::counter::Counter;
use crate::plotting::{FG_BLUE, FG_TEAL, standard_layout, write_plot_pdf};

/// Counts of positions by depth: `counts()[d]` is the number of positions at depth `d`.
///
/// A capped histogram has a fixed `cap + 1` bins and tallies any deeper position in the last
/// bin, so every statistic treats such positions as being at the cap. An unbounded histogram
/// grows to fit the deepest position it is given.
#[derive(Debug, Clone)]
pub struct DepthHistogram {
    counts: Vec<u64>,
    /// Largest depth tallied as itself; `None` for an unbounded histogram.
    cap: Option<u64>,
}

impl DepthHistogram {
    /// A histogram that tallies depths above `cap` at `cap`.
    #[must_use]
    pub fn capped(cap: u16) -> Self {
        Self { counts: vec![0; usize::from(cap) + 1], cap: Some(u64::from(cap)) }
    }

    /// A histogram that grows to fit any depth.
    #[must_use]
    pub fn unbounded() -> Self {
        Self { counts: Vec::new(), cap: None }
    }

    /// Tally one position at `depth`.
    #[inline]
    pub fn add(&mut self, depth: u64) {
        self.add_n(depth, 1);
    }

    /// Tally `n` positions at `depth`.
    #[inline]
    pub fn add_n(&mut self, depth: u64, n: u64) {
        let depth = self.cap.map_or(depth, |cap| depth.min(cap));
        // Depths index a histogram held in memory, so they fit in usize.
        #[allow(clippy::cast_possible_truncation)]
        let index = depth as usize;
        if index >= self.counts.len() {
            self.counts.resize(index + 1, 0);
        }
        self.counts[index] += n;
    }

    /// The count of positions at each depth, indexed by depth.
    #[must_use]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    /// The number of positions tallied.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.counts.iter().sum()
    }

    /// The sum of the tallied depths.
    #[must_use]
    pub fn sum(&self) -> f64 {
        self.counts.iter().enumerate().map(|(depth, &count)| depth as f64 * count as f64).sum()
    }

    /// The mean depth and its sample standard deviation (`n - 1` denominator), both 0 when
    /// nothing is tallied.
    #[must_use]
    pub fn mean_and_sd(&self) -> (f64, f64) {
        let total = self.total();
        if total == 0 {
            return (0.0, 0.0);
        }
        let n = total as f64;
        let mut sum = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        for (depth, &count) in self.counts.iter().enumerate() {
            if count > 0 {
                let d = depth as f64;
                let c = count as f64;
                sum += d * c;
                sum_sq += d * d * c;
            }
        }
        let mean = sum / n;
        let variance = if total > 1 { (sum_sq - n * mean * mean) / (n - 1.0) } else { 0.0 };
        (mean, variance.max(0.0).sqrt())
    }

    /// The median depth (the mean of the two middle depths for an even count) and the median
    /// absolute deviation from it.
    #[must_use]
    pub fn median_and_mad(&self) -> (f64, f64) {
        let counter: Counter<u64> = self
            .counts
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 0)
            .map(|(depth, &count)| (depth as u64, count))
            .collect();
        counter.median_and_mad()
    }

    /// For each depth `d`, the number of positions at depth `d` or above.
    #[must_use]
    pub fn at_or_above(&self) -> Vec<u64> {
        let mut at_or_above = vec![0u64; self.counts.len()];
        let mut running: u64 = 0;
        for depth in (0..self.counts.len()).rev() {
            running += self.counts[depth];
            at_or_above[depth] = running;
        }
        at_or_above
    }

    /// The fraction of positions at `depth` or above; 0 when nothing is tallied or `depth` is
    /// beyond the deepest bin.
    #[must_use]
    pub fn frac_at_or_above(&self, depth: usize) -> f64 {
        let total = self.total();
        if total == 0 || depth >= self.counts.len() {
            return 0.0;
        }
        self.counts[depth..].iter().sum::<u64>() as f64 / total as f64
    }

    /// The smallest depth at which at least `lower_frac` of positions are at or below it: for
    /// `lower_frac = 0.2`, the depth that 80% of positions reach. The target count is
    /// `ceil(total × lower_frac)`. Returns 0 when nothing is tallied.
    #[must_use]
    pub fn percentile(&self, lower_frac: f64) -> f64 {
        let total = self.total();
        if total == 0 {
            return 0.0;
        }
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "the target count is non-negative and at most the total"
        )]
        let threshold = (total as f64 * lower_frac).ceil() as u64;
        let mut running: u64 = 0;
        for (depth, &count) in self.counts.iter().enumerate() {
            running += count;
            if running >= threshold {
                return depth as f64;
            }
        }
        (self.counts.len() - 1) as f64
    }

    /// The fold penalty at `lower_frac`: `mean` divided by [`Self::percentile`], e.g. the
    /// fold-80 penalty at `lower_frac = 0.2`. 0 when that percentile is 0.
    #[must_use]
    pub fn fold_penalty(&self, mean: f64, lower_frac: f64) -> f64 {
        let percentile = self.percentile(lower_frac);
        if percentile == 0.0 { 0.0 } else { mean / percentile }
    }

    /// Write a PDF area chart of the distribution to `path`: the fraction of positions at each
    /// depth, up to the 99.9th-percentile depth so the long tail is trimmed, with a Poisson
    /// distribution at `mean` (scaled to the covered fraction) overlaid as a dashed line, and an
    /// optional subtitle under the title. Writes nothing when nothing is tallied or every
    /// position is at depth 0.
    ///
    /// # Errors
    /// Returns an error if the plot cannot be rendered or written.
    pub fn plot_pdf(
        &self,
        mean: f64,
        title: &str,
        subtitle: Option<&str>,
        x_label: &str,
        y_label: &str,
        path: &Path,
    ) -> Result<()> {
        let total = self.total();
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a percentile depth indexes the histogram and is non-negative"
        )]
        let x_max = self.percentile(0.999) as u64;
        if total == 0 || x_max == 0 {
            return Ok(());
        }
        let territory = total as f64;

        // x_max indexes the histogram, so it fits in usize.
        #[allow(clippy::cast_possible_truncation)]
        let x_max_index = x_max as usize;
        let observed_xy: Vec<(f64, f64)> = self
            .counts
            .iter()
            .enumerate()
            .take(x_max_index + 1)
            .map(|(depth, &count)| (depth as f64, count as f64 / territory))
            .collect();

        let covered_frac = (territory - self.counts[0] as f64) / territory;
        let poisson_xy: Vec<(f64, f64)> =
            (0..=x_max).map(|k| (k as f64, poisson_pmf(k, mean) * covered_frac)).collect();

        let observed = LinePlot::new()
            .with_data(observed_xy)
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

        let mut layout = standard_layout(&plots).with_title(title);
        if let Some(subtitle) = subtitle {
            layout = layout.with_subtitle(subtitle);
        }
        let layout = layout
            .with_x_label(x_label)
            .with_y_label(y_label)
            .with_minor_ticks(5)
            .with_show_minor_grid(true)
            .with_legend_position(LegendPosition::InsideTopRight);

        write_plot_pdf(plots, layout, path)
    }
}

/// The Poisson probability P(X = k) at `lambda`, computed in log space to avoid overflow for
/// large `k`.
fn poisson_pmf(k: u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    let log_pmf = k as f64 * lambda.ln() - lambda - ln_factorial(k);
    log_pmf.exp()
}

/// ln(k!), summed term by term.
fn ln_factorial(k: u64) -> f64 {
    (2..=k).map(|i| (i as f64).ln()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unbounded_from(counts: &[u64]) -> DepthHistogram {
        let mut histogram = DepthHistogram::unbounded();
        for (depth, &count) in counts.iter().enumerate() {
            histogram.add_n(depth as u64, count);
        }
        histogram
    }

    #[test]
    fn capped_histogram_tallies_deeper_positions_at_the_cap() {
        let mut histogram = DepthHistogram::capped(3);
        histogram.add(1);
        histogram.add(3);
        histogram.add(10);
        assert_eq!(histogram.counts(), &[0, 1, 0, 2]);
    }

    #[test]
    fn unbounded_histogram_grows_to_fit_the_deepest_position() {
        let mut histogram = DepthHistogram::unbounded();
        histogram.add_n(5, 2);
        assert_eq!(histogram.counts(), &[0, 0, 0, 0, 0, 2]);
        assert_eq!(histogram.total(), 2);
    }

    #[test]
    fn sum_weights_each_depth_by_its_count() {
        crate::assert_close!(unbounded_from(&[4, 2, 1]).sum(), 4.0);
    }

    #[test]
    fn mean_and_sd_use_the_sample_standard_deviation() {
        // Depths 1, 2, 3: mean 2, sample SD 1.
        let (mean, sd) = unbounded_from(&[0, 1, 1, 1]).mean_and_sd();
        crate::assert_close!(mean, 2.0);
        crate::assert_close!(sd, 1.0);
    }

    #[test]
    fn mean_and_sd_are_zero_when_empty() {
        assert_eq!(DepthHistogram::unbounded().mean_and_sd(), (0.0, 0.0));
    }

    #[test]
    fn median_of_an_odd_count_is_the_middle_depth() {
        // 5 positions at depth 0, 10 at depth 1.
        crate::assert_close!(unbounded_from(&[5, 10]).median_and_mad().0, 1.0);
    }

    #[test]
    fn median_of_an_even_count_averages_the_two_middle_depths() {
        crate::assert_close!(unbounded_from(&[0, 2, 2]).median_and_mad().0, 1.5);
    }

    #[test]
    fn median_of_an_even_count_with_equal_middle_depths_is_that_depth() {
        crate::assert_close!(unbounded_from(&[1, 3]).median_and_mad().0, 1.0);
    }

    #[test]
    fn median_of_an_empty_histogram_is_zero() {
        crate::assert_close!(DepthHistogram::unbounded().median_and_mad().0, 0.0);
    }

    #[test]
    fn mad_is_the_median_distance_from_the_median() {
        // Depths 1, 2, 2, 2, 9: median 2, deviations 1, 0, 0, 0, 7 → MAD 0.
        let mut histogram = unbounded_from(&[0, 1, 3]);
        histogram.add(9);
        assert_eq!(histogram.median_and_mad(), (2.0, 0.0));
    }

    #[test]
    fn at_or_above_is_the_reverse_cumulative_count() {
        assert_eq!(unbounded_from(&[5, 10, 3, 2]).at_or_above(), vec![20, 15, 5, 2]);
    }

    #[test]
    fn frac_at_or_above_divides_by_the_total() {
        let histogram = unbounded_from(&[5, 10, 3, 2]);
        crate::assert_close!(histogram.frac_at_or_above(2), 0.25);
        crate::assert_close!(histogram.frac_at_or_above(0), 1.0);
    }

    #[test]
    fn frac_at_or_above_beyond_the_deepest_bin_is_zero() {
        crate::assert_close!(unbounded_from(&[5, 10]).frac_at_or_above(7), 0.0);
    }

    #[test]
    fn percentile_rounds_the_target_count_up() {
        // Total 3: the 20th-percentile target is ceil(0.6) = 1, reached at depth 0.
        crate::assert_close!(unbounded_from(&[1, 1, 1]).percentile(0.20), 0.0);
    }

    #[test]
    fn percentile_with_every_position_at_one_depth_is_that_depth() {
        let mut histogram = DepthHistogram::unbounded();
        histogram.add_n(10, 10);
        crate::assert_close!(histogram.percentile(0.20), 10.0);
    }

    #[test]
    fn percentile_at_one_is_the_deepest_occupied_depth() {
        let mut histogram = DepthHistogram::unbounded();
        histogram.add_n(5, 5);
        histogram.add_n(10, 5);
        crate::assert_close!(histogram.percentile(1.0), 10.0);
        crate::assert_close!(histogram.percentile(0.5), 5.0);
    }

    #[test]
    fn percentile_of_an_empty_histogram_is_zero() {
        crate::assert_close!(DepthHistogram::unbounded().percentile(0.5), 0.0);
    }

    #[test]
    fn percentile_of_a_small_territory_finds_the_first_covered_depth() {
        // Total 4: the 20th-percentile target is ceil(0.8) = 1, so depth 0 (holding nothing)
        // does not satisfy it.
        let mut histogram = DepthHistogram::unbounded();
        histogram.add_n(7, 4);
        crate::assert_close!(histogram.percentile(0.20), 7.0);
    }

    #[test]
    fn fold_penalty_divides_the_mean_by_the_percentile() {
        let mut histogram = DepthHistogram::unbounded();
        histogram.add_n(5, 5);
        histogram.add_n(15, 5);
        crate::assert_close!(histogram.fold_penalty(10.0, 0.20), 2.0);
    }

    #[test]
    fn fold_penalty_is_zero_when_the_percentile_is_zero() {
        crate::assert_close!(unbounded_from(&[5, 5]).fold_penalty(0.5, 0.20), 0.0);
    }

    #[test]
    fn poisson_pmf_sums_to_one() {
        let total: f64 = (0..100).map(|k| poisson_pmf(k, 10.0)).sum();
        crate::assert_close!(total, 1.0, 1e-9);
    }
}
