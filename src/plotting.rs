use std::path::Path;

use anyhow::{Result, anyhow};
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use kuva::render::render::render_twin_y;

// ─── Fulcrum Genomics brand colours (hex for kuva) ──────────────────────────

pub const FG_BLUE: &str = "#26a8e0";
pub const FG_GREEN: &str = "#38b44a";
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
