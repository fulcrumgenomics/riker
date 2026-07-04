use std::path::Path;

use anyhow::{Result, anyhow};
use kuva::plot::{LegendPosition, LinePlot};
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use kuva::render::render::render_twin_y;
use kuva::render::render_utils::compute_tick_step;

use crate::counter::Counter;

// ─── Fulcrum Genomics brand colours (hex for kuva) ──────────────────────────

pub const FG_BLUE: &str = "#26a8e0";
pub const FG_GREEN: &str = "#38b44a";
/// Core dark blue (PMS 274) — good for a single solid series distinct from FG Blue.
pub const FG_NAVY: &str = "#160052";
/// Secondary dark green (PMS 7736). Part of the brand palette; kept for completeness
/// even though no plot currently uses it.
#[allow(dead_code)]
pub const FG_PINE: &str = "#315848";
pub const FG_PACIFIC: &str = "#1693b9";
pub const FG_SKY: &str = "#4dcce8";
pub const FG_TEAL: &str = "#2fae99";
pub const FG_EMERALD: &str = "#4dcc68";
pub const FG_FOREST: &str = "#269e2a";

/// Warm red for warning/alert data series.
pub const FG_RED: &str = "#e04040";

/// Neutral gray for reference lines and gridlines.
pub const FG_GRAY: &str = "#b4b4b4";

// ─── Standard plot dimensions (8 × 6 inches at 72 DPI) ─────────────────────

/// Standard plot width in pixels.
pub const PLOT_WIDTH: f64 = 800.0;
/// Standard plot height in pixels.
pub const PLOT_HEIGHT: f64 = 600.0;

/// A [`Layout`] auto-fitted to `plots` at the standard riker plot size
/// ([`PLOT_WIDTH`] × [`PLOT_HEIGHT`]). Chain further `.with_*` setters for the
/// per-plot title, labels, and axes.
#[must_use]
pub(crate) fn standard_layout(plots: &[Plot]) -> Layout {
    Layout::auto_from_plots(plots).with_width(PLOT_WIDTH).with_height(PLOT_HEIGHT)
}

// ─── PDF rendering helpers ──────────────────────────────────────────────────

/// Render plots to a PDF file.
///
/// Builds a single-Y-axis chart from `plots` and `layout`, renders to PDF,
/// and writes the result to `path`.
///
/// # Errors
/// Returns an error if PDF rendering fails or the file cannot be written.
pub fn write_plot_pdf(plots: Vec<Plot>, layout: Layout, path: &Path) -> Result<()> {
    let pdf_bytes = kuva::render_to_pdf(plots, layout).map_err(|e| anyhow!("{e}"))?;
    std::fs::write(path, pdf_bytes)?;
    Ok(())
}

/// Render a dual-Y-axis chart to a PDF file.
///
/// `primary` plots are drawn against the left Y axis; `secondary` plots are
/// drawn against the right Y2 axis.
///
/// # Errors
/// Returns an error if PDF rendering fails or the file cannot be written.
pub fn write_twin_y_plot_pdf(
    primary: Vec<Plot>,
    secondary: Vec<Plot>,
    layout: Layout,
    path: &Path,
) -> Result<()> {
    let scene = render_twin_y(primary, secondary, layout);
    let pdf_bytes =
        kuva::backend::pdf::PdfBackend.render_scene(&scene).map_err(|e| anyhow!("{e}"))?;
    std::fs::write(path, pdf_bytes)?;
    Ok(())
}

// ─── Insert-size histogram chart ──────────────────────────────────────────────

/// One pair-orientation series for an insert-size histogram chart.
pub struct InsertSizeSeries<'a> {
    /// Legend label (e.g. `"FR"`).
    pub name: &'a str,
    /// Fill/line colour (an FG brand hex).
    pub color: &'a str,
    /// The per-insert-size counts.
    pub histogram: &'a Counter<u64>,
}

/// Write an insert-size distribution area chart to `path`, shared by the genomic `isize`
/// command and the RNA transcript-space insert-size metrics so the two charts match.
///
/// Each series is trimmed to `median + deviations * MAD` and dropped entirely when its
/// fraction of the total falls below `min_frac`. Returns `Ok(())` without writing a file
/// when no series qualifies (mirrors the caller's empty-output behaviour).
///
/// # Errors
/// Returns an error if PDF rendering fails or the file cannot be written.
pub fn write_insert_size_histogram_pdf(
    path: &Path,
    title: &str,
    series: &[InsertSizeSeries<'_>],
    min_frac: f64,
    deviations: f64,
) -> Result<()> {
    let total: u64 = series.iter().map(|s| s.histogram.total()).sum();

    let mut plots: Vec<Plot> = Vec::new();
    let mut data_x_max: f64 = 0.0;
    for s in series {
        let count = s.histogram.total();
        if count == 0 || (total > 0 && (count as f64) < (total as f64) * min_frac) {
            continue;
        }
        let (median, mad) = s.histogram.median_and_mad();
        let trim_max = median + deviations * mad;
        let xy: Vec<(f64, f64)> = s
            .histogram
            .sorted()
            .into_iter()
            .take_while(|&(k, _)| k as f64 <= trim_max)
            .map(|(x, y)| (x as f64, y as f64))
            .collect();
        if let Some(&(x, _)) = xy.last() {
            data_x_max = data_x_max.max(x);
        }
        plots.push(Plot::Line(
            LinePlot::new()
                .with_data(xy)
                .with_color(s.color)
                .with_fill()
                .with_fill_opacity(0.3)
                .with_legend(s.name),
        ));
    }

    if plots.is_empty() {
        return Ok(());
    }

    // Snap the x-axis max up to a clean multiple of kuva's tick step so the last major
    // gridline lands on the axis edge (see the isize collector's original note).
    let x_tick_step = compute_tick_step(0.0, data_x_max, 10);
    let x_axis_max = (data_x_max / x_tick_step).ceil() * x_tick_step;

    let layout = standard_layout(&plots)
        .with_title(title)
        .with_x_label("Insert Size (bp)")
        .with_y_label("Read Pairs")
        .with_x_axis_min(0.0)
        .with_x_axis_max(x_axis_max)
        .with_x_tick_step(x_tick_step)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        // Insert-size distributions decay to the right, leaving the top-right corner clear.
        .with_legend_position(LegendPosition::InsideTopRight);

    write_plot_pdf(plots, layout, path)
}
