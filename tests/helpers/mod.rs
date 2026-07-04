//! Shared helpers for the integration test crate.
//!
//! Test inputs (records, references, intervals, VCFs) are built with the unified
//! `riker_lib::test_support` builders; this module keeps only the output-reading
//! helper, which belongs to the test crate rather than the library.

use std::path::Path;

use anyhow::Result;
use fgoxide::io::DelimFile;
use serde::de::DeserializeOwned;

/// Read a tab-separated metrics file into a vector of typed rows.
///
/// Backed by [`DelimFile::read_tsv`], which dispatches to gzip or plain
/// reading based on the path's extension — `.gz` and `.bgz` files are
/// transparently decompressed.
///
/// # Errors
/// Returns an error if the file cannot be read or deserialized.
pub fn read_metrics_tsv<T: DeserializeOwned>(path: &Path) -> Result<Vec<T>> {
    Ok(DelimFile::default().read_tsv(path)?)
}
