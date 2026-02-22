#[allow(dead_code)]
mod helpers;

use anyhow::Result;
use riker_lib::commands::command::Command;
use riker_lib::commands::docs::Docs;
use tempfile::TempDir;

/// Run `docs --format text` to a file and verify the output contains
/// expected metric type headings and field names.
#[test]
fn test_docs_text_output() -> Result<()> {
    let dir = TempDir::new()?;
    let out = dir.path().join("docs.txt");

    Docs { format: "text".to_string(), output: Some(out.clone()) }.execute()?;

    let text = std::fs::read_to_string(&out)?;

    // Metric type headings
    assert!(text.contains("Alignment summary metrics"), "missing AlignmentSummaryMetric heading");
    assert!(text.contains("Insert size distribution metrics"), "missing InsertSizeMetric heading");
    assert!(
        text.contains("Whole-genome sequencing coverage summary"),
        "missing WgsMetrics heading"
    );
    assert!(text.contains("GC bias detail metrics"), "missing GcBiasDetailMetric heading");
    assert!(text.contains("GC bias summary metrics"), "missing GcBiasSummaryMetric heading");

    // Spot-check field names from different metrics
    assert!(text.contains("frac_aligned"), "missing alignment field");
    assert!(text.contains("mean_insert_size"), "missing isize field");
    assert!(text.contains("mean_coverage"), "missing wgs field");
    assert!(text.contains("normalized_coverage"), "missing gcbias field");
    assert!(text.contains("at_dropout"), "missing gcbias summary field");

    Ok(())
}

/// Run `docs --format markdown` to a file and verify the output contains
/// markdown table structure and expected content.
#[test]
fn test_docs_markdown_output() -> Result<()> {
    let dir = TempDir::new()?;
    let out = dir.path().join("docs.md");

    Docs { format: "markdown".to_string(), output: Some(out.clone()) }.execute()?;

    let text = std::fs::read_to_string(&out)?;

    // Markdown headings
    assert!(text.contains("### Alignment summary metrics"), "missing markdown heading");
    assert!(text.contains("### Insert size distribution metrics"), "missing markdown heading");
    assert!(text.contains("### GC bias detail metrics"), "missing markdown heading");

    // Table structure
    assert!(text.contains("| Column | Description |"), "missing markdown table header");
    assert!(text.contains("| `frac_aligned`"), "missing alignment field in markdown");
    assert!(text.contains("| `median_insert_size`"), "missing isize field in markdown");
    assert!(text.contains("| `genome_territory`"), "missing wgs field in markdown");

    Ok(())
}

/// Defaults to stdout (no output path) without error.
#[test]
fn test_docs_stdout_default() -> Result<()> {
    Docs { format: "text".to_string(), output: None }.execute()
}

/// Unknown format should produce an error.
#[test]
fn test_docs_unknown_format() {
    let result = Docs { format: "csv".to_string(), output: None }.execute();
    assert!(result.is_err());
}
