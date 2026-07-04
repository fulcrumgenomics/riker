use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use anyhow::Result;
use clap::Args;

use crate::commands::alignment::AlignmentSummaryMetric;
use crate::commands::basic::{
    BaseDistributionByCycleMetric, MeanQualityByCycleMetric, QualityScoreDistributionMetric,
};
use crate::commands::command::Command;
use crate::commands::error::{IndelMetric, MismatchMetric, OverlappingMismatchMetric};
use crate::commands::gcbias::{GcBiasDetailMetric, GcBiasSummaryMetric};
use crate::commands::hybcap::{HybCapMetric, PerBaseCoverage, PerTargetCoverage};
use crate::commands::isize::{InsertSizeHistogramEntry, InsertSizeMetric};
use crate::commands::rna::{
    RnaBiotypeMetric, RnaInsertSizeHistogramEntry, RnaInsertSizeMetric, RnaSeqMetric,
};
use crate::commands::wgs::{WgsCoverageEntry, WgsMetrics};
use crate::metrics::{render_metric_docs_markdown, render_metric_docs_text};

/// Render docs for every metric type, in a fixed order, separated by blank lines, with
/// the given per-type renderer (`render_metric_docs_text` / `render_metric_docs_markdown`).
/// The type list lives here once so the output formats can't drift — adding a metric means
/// editing one place, and every format stays in sync.
macro_rules! render_all_metric_docs {
    ($render:ident, $out:expr) => {
        render_all_metric_docs!(@seq $render, $out,
            AlignmentSummaryMetric,
            BaseDistributionByCycleMetric,
            GcBiasDetailMetric,
            GcBiasSummaryMetric,
            HybCapMetric,
            PerTargetCoverage,
            PerBaseCoverage,
            InsertSizeMetric,
            InsertSizeHistogramEntry,
            RnaSeqMetric,
            RnaBiotypeMetric,
            RnaInsertSizeMetric,
            RnaInsertSizeHistogramEntry,
            MeanQualityByCycleMetric,
            QualityScoreDistributionMetric,
            WgsMetrics,
            WgsCoverageEntry,
            MismatchMetric,
            OverlappingMismatchMetric,
            IndelMetric,
        )
    };
    // Blank line between entries, none after the last: render the first, then
    // `writeln` + render for each of the rest.
    (@seq $render:ident, $out:expr, $first:ty, $($rest:ty,)*) => {{
        $render::<$first>($out)?;
        $(
            writeln!($out)?;
            $render::<$rest>($out)?;
        )*
    }};
}

/// Print documentation for all metric types produced by riker.
///
/// Lists every metric field with its description for each output file
/// type. Supports plain text and markdown table formats. Output is
/// written to stdout, or to the file specified by -o/--output.
#[derive(Args, Debug, Clone)]
#[command(long_about)]
pub struct Docs {
    /// Output format: "text" for plain text, "markdown" for markdown tables.
    #[arg(long, default_value = "text")]
    pub format: String,

    /// Output file path. Defaults to stdout when not specified.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,
}

impl Command for Docs {
    fn execute(&self, _threads: Option<u8>) -> Result<()> {
        let stdout;
        let mut out: Box<dyn Write> = if let Some(path) = &self.output {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            stdout = io::stdout();
            Box::new(stdout.lock())
        };

        match self.format.as_str() {
            "text" => render_all_metric_docs!(render_metric_docs_text, &mut out),
            "markdown" => render_all_metric_docs!(render_metric_docs_markdown, &mut out),
            other => {
                anyhow::bail!("Unknown format '{other}'. Use 'text' or 'markdown'.");
            }
        }

        Ok(())
    }
}
