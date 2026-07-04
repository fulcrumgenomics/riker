//! A tiny builder for BED interval files — baits/targets for `hybcap`, intervals
//! for `wgs`. Accumulate intervals, write them to disk.

// Builder setters return `Self`; allow the "could be `#[must_use]`" lints rather
// than annotate each.
#![allow(clippy::must_use_candidate, clippy::return_self_not_must_use)]

use std::io::{self, Write};
use std::path::Path;

/// Accumulates BED intervals and writes them to a `.bed` file.
///
/// Coordinates are **0-based, half-open** — the BED convention, and what the tools
/// parse. Note this differs from the 1-based positions [`read`](super::read) /
/// [`pair`](super::pair) take; an interval `contig 0 200` covers the reference bases
/// a read at 1-based `contig:1` with a `200M` CIGAR spans.
#[derive(Clone, Default)]
pub struct BedBuilder {
    intervals: Vec<(String, usize, usize, Option<String>)>,
}

impl BedBuilder {
    /// An empty BED.
    pub fn new() -> Self {
        Self::default()
    }

    /// A 3-column interval `contig<TAB>start<TAB>end` (0-based, half-open).
    pub fn interval(mut self, contig: impl Into<String>, start: usize, end: usize) -> Self {
        self.intervals.push((contig.into(), start, end, None));
        self
    }

    /// A 4-column interval `contig<TAB>start<TAB>end<TAB>name` (0-based, half-open).
    pub fn named(
        mut self,
        contig: impl Into<String>,
        start: usize,
        end: usize,
        name: impl Into<String>,
    ) -> Self {
        self.intervals.push((contig.into(), start, end, Some(name.into())));
        self
    }

    /// Write the BED to `path`.
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or written.
    pub fn write_bed(&self, path: &Path) -> io::Result<()> {
        let mut file = std::fs::File::create(path)?;
        for (contig, start, end, name) in &self.intervals {
            match name {
                Some(name) => writeln!(file, "{contig}\t{start}\t{end}\t{name}")?,
                None => writeln!(file, "{contig}\t{start}\t{end}")?,
            }
        }
        Ok(())
    }

    /// Write to a temporary `.bed` and return the handle, which deletes it on drop.
    /// For a test where baits and targets are identical, build once and pass the one
    /// handle's path as both.
    ///
    /// # Errors
    /// Returns an error if the temp file cannot be created or written.
    pub fn to_temp_bed(&self) -> io::Result<tempfile::NamedTempFile> {
        let tmp = tempfile::NamedTempFile::with_suffix(".bed")?;
        self.write_bed(tmp.path())?;
        Ok(tmp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_three_column_intervals() {
        let bed = BedBuilder::new().interval("chr1", 0, 200).interval("chr2", 100, 300);
        let tmp = bed.to_temp_bed().unwrap();
        assert_eq!(std::fs::read_to_string(tmp.path()).unwrap(), "chr1\t0\t200\nchr2\t100\t300\n");
    }

    #[test]
    fn writes_four_column_named_intervals() {
        let bed = BedBuilder::new().named("chr1", 0, 200, "bait1");
        let tmp = bed.to_temp_bed().unwrap();
        assert_eq!(std::fs::read_to_string(tmp.path()).unwrap(), "chr1\t0\t200\tbait1\n");
    }

    #[test]
    fn mixes_named_and_unnamed_intervals_in_order() {
        let bed = BedBuilder::new().interval("chr1", 0, 10).named("chr1", 90, 100, "t2");
        let tmp = bed.to_temp_bed().unwrap();
        assert_eq!(
            std::fs::read_to_string(tmp.path()).unwrap(),
            "chr1\t0\t10\nchr1\t90\t100\tt2\n"
        );
    }

    #[test]
    fn an_empty_bed_is_an_empty_file() {
        let tmp = BedBuilder::new().to_temp_bed().unwrap();
        assert_eq!(std::fs::read_to_string(tmp.path()).unwrap(), "");
    }
}
