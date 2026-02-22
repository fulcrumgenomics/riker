use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use csv::WriterBuilder;
use serde::Serialize;

// ─── MetricDocs trait ─────────────────────────────────────────────────────────

/// Describes a single field in a metric struct.
pub struct FieldDoc {
    pub name: &'static str,
    pub description: &'static str,
}

/// Trait implemented by `#[derive(MetricDocs)]` to expose doc comments at runtime.
pub trait MetricDocs {
    /// The struct-level doc comment.
    fn metric_description() -> &'static str;
    /// Per-field name and doc comment pairs.
    fn field_docs() -> &'static [FieldDoc];
}

/// Render metric documentation as plain text to a writer.
///
/// # Errors
/// Returns an error if writing fails.
pub fn render_metric_docs_text<T: MetricDocs>(w: &mut dyn Write) -> Result<()> {
    writeln!(w, "{}", T::metric_description())?;
    writeln!(w, "{}", "-".repeat(T::metric_description().len()))?;
    let max_name = T::field_docs().iter().map(|f| f.name.len()).max().unwrap_or(0);
    for field in T::field_docs() {
        writeln!(w, "  {:<width$}  {}", field.name, field.description, width = max_name)?;
    }
    Ok(())
}

/// Render metric documentation as a markdown table to a writer.
///
/// # Errors
/// Returns an error if writing fails.
pub fn render_metric_docs_markdown<T: MetricDocs>(w: &mut dyn Write) -> Result<()> {
    writeln!(w, "### {}", T::metric_description())?;
    writeln!(w)?;
    writeln!(w, "| Column | Description |")?;
    writeln!(w, "|--------|-------------|")?;
    for field in T::field_docs() {
        writeln!(w, "| `{}` | {} |", field.name, field.description)?;
    }
    Ok(())
}

/// Serialize an `f64` with 2 decimal places.
///
/// Use as `#[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]`.
///
/// # Errors
/// Propagates serializer errors.
pub fn serialize_f64_2dp<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("{value:.2}"))
}

/// Serialize an `f64` with 5 decimal places.
///
/// Use as `#[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]`.
///
/// # Errors
/// Propagates serializer errors.
pub fn serialize_f64_5dp<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("{value:.5}"))
}

/// Serialize an `f64` with 6 decimal places.
///
/// Use as `#[serde(serialize_with = "crate::metrics::serialize_f64_6dp")]`.
///
/// # Errors
/// Propagates serializer errors.
pub fn serialize_f64_6dp<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("{value:.6}"))
}

/// Serialize an `Option<f64>` with 5 decimal places, or empty string for `None`.
///
/// Use as `#[serde(serialize_with = "crate::metrics::serialize_opt_f64_5dp")]`.
///
/// # Errors
/// Propagates serializer errors.
pub fn serialize_opt_f64_5dp<S>(value: &Option<f64>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match value {
        Some(v) => serializer.serialize_str(&format!("{v:.5}")),
        None => serializer.serialize_str(""),
    }
}

/// Write a slice of serializable rows to a tab-separated file.
///
/// Produces lowercase headers (from serde field names), no metadata lines.
///
/// # Errors
/// Returns an error if the file cannot be created or a row cannot be serialized.
pub fn write_tsv<T: Serialize>(path: &Path, rows: &[T]) -> Result<()> {
    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("Failed to create TSV: {}", path.display()))?;
    for row in rows {
        wtr.serialize(row).context("Failed to serialize row to TSV")?;
    }
    wtr.flush().context("Failed to flush TSV writer")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use riker_derive::MetricDocs;
    use serde::{Deserialize, Serialize};
    use tempfile::NamedTempFile;

    #[test]
    fn test_serialize_f64_2dp() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "super::serialize_f64_2dp")]
            value: f64,
        }

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &[Row { value: 1.234_567_89 }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("1.23"), "got: {content}");

        write_tsv(tmp.path(), &[Row { value: 0.0 }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("0.00"), "got: {content}");
    }

    #[test]
    fn test_serialize_f64_5dp() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "super::serialize_f64_5dp")]
            value: f64,
        }

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &[Row { value: 1.234_567_89 }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("1.23457"), "got: {content}");

        write_tsv(tmp.path(), &[Row { value: 0.0 }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("0.00000"), "got: {content}");
    }

    #[test]
    fn test_write_tsv_roundtrip() {
        #[derive(Serialize, Deserialize, PartialEq, Debug)]
        struct Row {
            name: String,
            count: u64,
        }

        let rows = vec![
            Row { name: "fr".to_string(), count: 100 },
            Row { name: "rf".to_string(), count: 5 },
        ];

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &rows).unwrap();

        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(tmp.path()).unwrap();
        let result: Vec<Row> = rdr.deserialize().map(|r| r.unwrap()).collect();
        assert_eq!(result, rows);
    }

    #[test]
    fn test_write_tsv_headers_lowercase() {
        #[derive(Serialize)]
        struct Row {
            mean_insert_size: f64,
            read_pairs: u64,
        }

        let rows = vec![Row { mean_insert_size: 1.0, read_pairs: 10 }];
        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &rows).unwrap();

        let content = std::fs::read_to_string(tmp.path()).unwrap();
        let header_line = content.lines().next().unwrap();
        assert_eq!(header_line, "mean_insert_size\tread_pairs");
    }

    // ── MetricDocs tests ────────────────────────────────────────────────────

    /// A test metric with documented fields.
    #[derive(MetricDocs)]
    #[allow(dead_code)]
    struct TestMetric {
        /// The sample name.
        sample: String,
        /// Number of reads.
        count: u64,
        /// A field with
        /// multi-line doc.
        multi_line: f64,
        no_doc: u32,
    }

    #[test]
    fn test_metric_docs_description() {
        assert_eq!(TestMetric::metric_description(), "A test metric with documented fields.");
    }

    #[test]
    fn test_metric_docs_field_count() {
        assert_eq!(TestMetric::field_docs().len(), 4);
    }

    #[test]
    fn test_metric_docs_field_names_and_descs() {
        let docs = TestMetric::field_docs();
        assert_eq!(docs[0].name, "sample");
        assert_eq!(docs[0].description, "The sample name.");
        assert_eq!(docs[1].name, "count");
        assert_eq!(docs[1].description, "Number of reads.");
    }

    #[test]
    fn test_metric_docs_multi_line() {
        let docs = TestMetric::field_docs();
        assert_eq!(docs[2].name, "multi_line");
        assert_eq!(docs[2].description, "A field with multi-line doc.");
    }

    #[test]
    fn test_metric_docs_empty_doc() {
        let docs = TestMetric::field_docs();
        assert_eq!(docs[3].name, "no_doc");
        assert_eq!(docs[3].description, "");
    }

    #[test]
    fn test_render_text() {
        let mut buf = Vec::new();
        render_metric_docs_text::<TestMetric>(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("A test metric with documented fields."));
        assert!(text.contains("sample"));
        assert!(text.contains("The sample name."));
    }

    #[test]
    fn test_render_markdown() {
        let mut buf = Vec::new();
        render_metric_docs_markdown::<TestMetric>(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("### A test metric with documented fields."));
        assert!(text.contains("| `sample` | The sample name. |"));
        assert!(text.contains("| Column | Description |"));
    }

    // ── Additional serializer tests ──────────────────────────────────────────

    #[test]
    fn test_serialize_f64_6dp() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "super::serialize_f64_6dp")]
            value: f64,
        }

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &[Row { value: 1.234_567_891 }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("1.234568"), "got: {content}");
    }

    #[test]
    fn test_serialize_opt_f64_5dp_some() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "super::serialize_opt_f64_5dp")]
            value: Option<f64>,
        }

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &[Row { value: Some(1.234_567_89) }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(content.contains("1.23457"), "got: {content}");
    }

    #[test]
    fn test_serialize_opt_f64_5dp_none() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "super::serialize_opt_f64_5dp")]
            value: Option<f64>,
        }

        let tmp = NamedTempFile::new().unwrap();
        write_tsv(tmp.path(), &[Row { value: None }]).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        // Header line + data line; the empty string is quoted by the CSV writer
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1], "\"\"");
    }

    #[test]
    fn test_write_tsv_empty_rows() {
        #[derive(Serialize)]
        struct Row {
            name: String,
            count: u64,
        }

        let tmp = NamedTempFile::new().unwrap();
        let rows: Vec<Row> = Vec::new();
        write_tsv(tmp.path(), &rows).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        // Empty rows → header-only output (just the header line or empty)
        assert!(content.is_empty() || content.trim().is_empty());
    }

    #[test]
    fn test_render_metric_docs_text_no_fields() {
        /// A metric with no fields.
        #[derive(MetricDocs)]
        struct EmptyMetric {}

        let mut buf = Vec::new();
        render_metric_docs_text::<EmptyMetric>(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("A metric with no fields."));
    }
}
