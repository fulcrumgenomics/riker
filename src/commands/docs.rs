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
use crate::commands::wgs::{WgsCoverageEntry, WgsMetrics};
use crate::metrics::{render_metric_docs_markdown, render_metric_docs_text};

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
    fn execute(&self) -> Result<()> {
        let stdout;
        let mut out: Box<dyn Write> = if let Some(path) = &self.output {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            stdout = io::stdout();
            Box::new(stdout.lock())
        };

        match self.format.as_str() {
            "text" => {
                render_metric_docs_text::<AlignmentSummaryMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<BaseDistributionByCycleMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<GcBiasDetailMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<GcBiasSummaryMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<HybCapMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<PerTargetCoverage>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<PerBaseCoverage>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<InsertSizeMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<InsertSizeHistogramEntry>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<MeanQualityByCycleMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<QualityScoreDistributionMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<WgsMetrics>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<WgsCoverageEntry>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<MismatchMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<OverlappingMismatchMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_text::<IndelMetric>(&mut out)?;
            }
            "markdown" => {
                render_metric_docs_markdown::<AlignmentSummaryMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<BaseDistributionByCycleMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<GcBiasDetailMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<GcBiasSummaryMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<HybCapMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<PerTargetCoverage>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<PerBaseCoverage>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<InsertSizeMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<InsertSizeHistogramEntry>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<MeanQualityByCycleMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<QualityScoreDistributionMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<WgsMetrics>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<WgsCoverageEntry>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<MismatchMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<OverlappingMismatchMetric>(&mut out)?;
                writeln!(out)?;
                render_metric_docs_markdown::<IndelMetric>(&mut out)?;
            }
            other => {
                anyhow::bail!("Unknown format '{other}'. Use 'text' or 'markdown'.");
            }
        }

        Ok(())
    }
}
