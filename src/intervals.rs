use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::sequence_dict::{SequenceDictionary, SequenceMetadata};
use anyhow::{Context, Result, bail};
use bitvec::prelude::*;
use flate2::read::MultiGzDecoder;
use noodles::sam::Header;

/// The gzip magic number (first two bytes of any gzip/bgzip file).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// A 0-based half-open reference interval `[start, end)` with an optional name.
///
/// Coordinates are stored as `u32`, which supports reference contigs up to ~4.3 Gbp.
/// This is sufficient for all known reference genomes (human chr1 is ~249 Mbp).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    /// Index into the BAM header's reference sequence list.
    pub ref_id: usize,
    /// 0-based inclusive start.
    pub start: u32,
    /// 0-based exclusive end.
    pub end: u32,
    /// Interval name (always populated; synthesized from coordinates if absent in input).
    name: String,
}

impl Interval {
    /// Create a new interval with an explicit name.
    #[must_use]
    pub fn new(ref_id: usize, start: u32, end: u32, name: String) -> Self {
        Self { ref_id, start, end, name }
    }

    /// Create a new interval, synthesizing a name from the contig name and coordinates.
    #[must_use]
    pub fn with_contig_name(ref_id: usize, start: u32, end: u32, contig: &str) -> Self {
        Self { ref_id, start, end, name: format!("{contig}:{start}-{end}") }
    }

    /// Return the interval's name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the length of the interval in bases.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.end.saturating_sub(self.start)
    }

    /// Return `true` if the interval has zero length.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// A collection of intervals grouped by reference contig, supporting fast
/// point-in-interval queries and bitvec generation for per-base analysis.
pub struct Intervals {
    /// Sorted intervals per contig, indexed by `ref_id`. Empty Vec = no intervals.
    by_contig: Vec<Vec<Interval>>,
    /// Reference sequence dictionary for name/index/length lookups.
    dict: SequenceDictionary,
    /// Total number of intervals before any merging (raw count from the input file).
    raw_count: usize,
}

impl Intervals {
    /// Load intervals from a file and group them by contig using a sequence dictionary.
    ///
    /// Supports plain text, gzip, and bgzip files (detected by magic bytes).
    /// Format is auto-detected: if the first non-empty line starts with `@`, the file
    /// is treated as an IntervalList; otherwise as BED.
    ///
    /// For IntervalList files, the SAM header is parsed and validated against `bam_dict`:
    /// the interval list's `@SQ` lines must be a prefix of the BAM dictionary (same names,
    /// same lengths, same order).  Individual intervals are validated against the header.
    ///
    /// # Errors
    /// Returns an error if the file cannot be read, parsed, or fails validation.
    pub fn from_path(path: &Path, bam_dict: SequenceDictionary) -> Result<Self> {
        let raw = load_intervals(path, &bam_dict)?;
        let raw_count = raw.len();
        let n = bam_dict.len();

        let mut by_contig: Vec<Vec<Interval>> = (0..n).map(|_| Vec::new()).collect();
        for iv in raw {
            if iv.ref_id < n {
                by_contig[iv.ref_id].push(iv);
            }
        }
        for ivs in &mut by_contig {
            ivs.sort_unstable_by_key(|iv| (iv.start, iv.end));
        }

        Ok(Self { by_contig, dict: bam_dict, raw_count })
    }

    /// Return a reference to the underlying sequence dictionary.
    #[must_use]
    pub fn dict(&self) -> &SequenceDictionary {
        &self.dict
    }

    /// Return the total number of intervals loaded from the input file (before merging).
    #[must_use]
    pub fn raw_count(&self) -> usize {
        self.raw_count
    }

    /// Return the total number of intervals across all contigs.
    #[must_use]
    pub fn count(&self) -> usize {
        self.by_contig.iter().map(Vec::len).sum()
    }

    /// Return `true` if `pos` (0-based) falls within any interval for the given contig.
    #[inline]
    #[must_use]
    pub fn contains_pos(&self, ref_id: usize, pos: u32) -> bool {
        match self.by_contig.get(ref_id) {
            Some(ivs) if !ivs.is_empty() => {
                let idx = ivs.partition_point(|iv| iv.start <= pos);
                idx > 0 && ivs[idx - 1].end > pos
            }
            _ => false,
        }
    }

    /// Return `true` if the given contig has any intervals.
    #[must_use]
    pub fn has_contig(&self, ref_id: usize) -> bool {
        self.by_contig.get(ref_id).is_some_and(|ivs| !ivs.is_empty())
    }

    /// Return the sorted intervals for a contig (empty slice if none).
    #[must_use]
    pub fn get_contig(&self, ref_id: usize) -> &[Interval] {
        self.by_contig.get(ref_id).map_or(&[], Vec::as_slice)
    }

    /// Iterate over all intervals across all contigs.
    pub fn iter(&self) -> impl Iterator<Item = &Interval> {
        self.by_contig.iter().flat_map(|ivs| ivs.iter())
    }

    /// Return the total number of bases covered by all intervals (sum of lengths).
    #[must_use]
    pub fn territory(&self) -> u64 {
        self.by_contig.iter().flat_map(|ivs| ivs.iter()).map(|iv| u64::from(iv.len())).sum()
    }

    /// Return a `BitVec` with length = contig length, bits set to `true` for
    /// positions covered by one or more intervals. Returns an empty bitvec if
    /// `ref_id` is out of range or has no intervals.
    #[must_use]
    pub fn contig_bitvec(&self, ref_id: usize) -> BitVec {
        let ivs = match self.by_contig.get(ref_id) {
            Some(ivs) if !ivs.is_empty() => ivs,
            _ => return BitVec::EMPTY,
        };

        let len = self.dict.get_by_index(ref_id).map_or(0, SequenceMetadata::length);
        if len == 0 {
            return BitVec::EMPTY;
        }

        let mut bv = bitvec![0; len];
        for iv in ivs {
            let s = (iv.start as usize).min(len);
            let e = (iv.end as usize).min(len);
            if s < e {
                bv[s..e].fill(true);
            }
        }
        bv
    }

    /// Return a new `Intervals` with all overlapping/abutting intervals merged.
    ///
    /// Merged intervals concatenate names with `|`. Preserves the raw_count
    /// from the original (pre-merge) interval set.
    #[must_use]
    pub fn merged(&self) -> Self {
        let mut merged_by_contig: Vec<Vec<Interval>> = Vec::with_capacity(self.by_contig.len());

        for ivs in &self.by_contig {
            let mut merged: Vec<Interval> = Vec::new();
            for iv in ivs {
                if let Some(last) = merged.last_mut()
                    && iv.start <= last.end
                {
                    // Overlapping or abutting — extend
                    last.end = last.end.max(iv.end);
                    last.name = format!("{}|{}", last.name, iv.name);
                    continue;
                }
                merged.push(iv.clone());
            }
            merged_by_contig.push(merged);
        }

        Self { by_contig: merged_by_contig, dict: self.dict.clone(), raw_count: self.raw_count }
    }

    /// Return a new `Intervals` with each interval padded by `amount` on both sides,
    /// clamped to `[0, contig_length)`.
    #[must_use]
    pub fn padded(&self, amount: u32) -> Self {
        let mut padded_by_contig: Vec<Vec<Interval>> = Vec::with_capacity(self.by_contig.len());

        for (ref_id, ivs) in self.by_contig.iter().enumerate() {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "contig lengths fit in u32 (max ~4.3 Gbp)"
            )]
            let contig_len = self.dict.get_by_index(ref_id).map_or(u32::MAX, |m| m.length() as u32);
            let padded: Vec<Interval> = ivs
                .iter()
                .map(|iv| Interval {
                    ref_id: iv.ref_id,
                    start: iv.start.saturating_sub(amount),
                    end: iv.end.saturating_add(amount).min(contig_len),
                    name: iv.name.clone(),
                })
                .collect();
            padded_by_contig.push(padded);
        }

        Self { by_contig: padded_by_contig, dict: self.dict.clone(), raw_count: self.raw_count }
    }
}

/// Read the contents of a file as a UTF-8 string, transparently decompressing
/// gzip/bgzip files detected by magic bytes.
fn read_file_contents(path: &Path) -> Result<String> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open interval file: {}", path.display()))?;

    let mut magic = [0u8; 2];
    let bytes_read = file
        .read(&mut magic)
        .with_context(|| format!("Failed to read interval file: {}", path.display()))?;

    // Reopen the file so the reader starts from the beginning.
    let file = File::open(path)
        .with_context(|| format!("Failed to open interval file: {}", path.display()))?;

    let mut content = String::new();
    if bytes_read >= 2 && magic == GZIP_MAGIC {
        MultiGzDecoder::new(file)
            .read_to_string(&mut content)
            .with_context(|| format!("Failed to decompress interval file: {}", path.display()))?;
    } else {
        std::io::BufReader::new(file)
            .read_to_string(&mut content)
            .with_context(|| format!("Failed to read interval file: {}", path.display()))?;
    }

    Ok(content)
}

/// Load intervals from an IntervalList or BED file, resolving contig names via
/// the BAM sequence dictionary.
///
/// Format is auto-detected: if the first non-empty line starts with `@`, the file
/// is treated as an IntervalList; otherwise as BED.
///
/// * **IntervalList**: 1-based fully-closed → converted to 0-based half-open.
///   The SAM header is parsed and validated as a prefix of `bam_dict`.
/// * **BED**: already 0-based half-open.  Unknown contig names produce a warning
///   and are skipped.
///
/// The returned vec is sorted by `(ref_id, start)`.
///
/// # Errors
/// Returns an error if the file cannot be read, parsed, or fails validation.
fn load_intervals(path: &Path, bam_dict: &SequenceDictionary) -> Result<Vec<Interval>> {
    let content = read_file_contents(path)?;

    // Auto-detect format: IntervalList if the first non-empty line starts with '@'.
    let is_interval_list =
        content.lines().find(|l| !l.trim().is_empty()).is_some_and(|l| l.starts_with('@'));

    let mut intervals: Vec<Interval> = Vec::new();

    if is_interval_list {
        let (header_dict, first_data_line) = parse_and_validate_header(&content, bam_dict, path)?;
        let mut errors: Vec<String> = Vec::new();

        for (line_num, line) in content.lines().enumerate() {
            let line_num = line_num + 1; // 1-based for error messages
            let line = line.trim();
            if line.is_empty() || line.starts_with('@') {
                continue;
            }
            // Only parse data lines at or after the first data line.
            if line_num < first_data_line {
                continue;
            }

            let (contig, start, end, name) = parse_interval_list_line(line, line_num)?;
            validate_interval(&contig, start, end, line_num, &header_dict, &mut errors);
            if let Some(meta) = bam_dict.get_by_name(&contig) {
                let name = name.unwrap_or_else(|| format!("{contig}:{start}-{end}"));
                intervals.push(Interval::new(meta.index(), start, end, name));
            }
        }

        if !errors.is_empty() {
            bail!("Interval validation errors in {}:\n  {}", path.display(), errors.join("\n  "));
        }
    } else {
        for (line_num, line) in content.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let (contig, start, end, name) = parse_bed_line(line, line_num + 1)?;
            match bam_dict.get_by_name(&contig) {
                Some(meta) => {
                    let name = name.unwrap_or_else(|| format!("{contig}:{start}-{end}"));
                    intervals.push(Interval::new(meta.index(), start, end, name));
                }
                None => {
                    log::warn!("intervals: contig '{contig}' not found in BAM header — skipping");
                }
            }
        }
    }

    intervals.sort_unstable_by_key(|iv| (iv.ref_id, iv.start));
    Ok(intervals)
}

/// Parse the SAM header from an IntervalList file, validate it as a prefix of the
/// BAM dictionary, and return the header's own `SequenceDictionary` along with the
/// 1-based line number of the first data line.
fn parse_and_validate_header(
    content: &str,
    bam_dict: &SequenceDictionary,
    path: &Path,
) -> Result<(SequenceDictionary, usize)> {
    // Collect header lines and find the first data line.
    let mut header_text = String::new();
    let mut first_data_line = usize::MAX;

    for (line_num, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('@') {
            header_text.push_str(line);
            header_text.push('\n');
        } else {
            first_data_line = line_num + 1; // 1-based
            break;
        }
    }

    // If there were no data lines, set first_data_line past the end.
    if first_data_line == usize::MAX {
        first_data_line = content.lines().count() + 1;
    }

    // Parse the header.
    let header: Header = header_text.parse().with_context(|| {
        format!("Failed to parse SAM header in IntervalList file: {}", path.display())
    })?;

    let header_dict = SequenceDictionary::from(&header);

    if header_dict.is_empty() {
        bail!("IntervalList file has no @SQ lines in header: {}", path.display());
    }

    // Validate that the header dict is a prefix of the BAM dict.
    validate_dict_prefix(&header_dict, bam_dict, path)?;

    Ok((header_dict, first_data_line))
}

/// Validate that `interval_dict` is a strict prefix of `bam_dict`: same names,
/// same lengths, same order for the first N entries.
fn validate_dict_prefix(
    interval_dict: &SequenceDictionary,
    bam_dict: &SequenceDictionary,
    path: &Path,
) -> Result<()> {
    if interval_dict.len() > bam_dict.len() {
        bail!(
            "IntervalList file {} has {} contigs in its header but the BAM only has {}",
            path.display(),
            interval_dict.len(),
            bam_dict.len()
        );
    }

    let mut errors: Vec<String> = Vec::new();

    for (i, iv_meta) in interval_dict.iter().enumerate() {
        let bam_meta = &bam_dict[i];
        if iv_meta.name() != bam_meta.name() {
            errors.push(format!(
                "  contig at index {i}: IntervalList has '{}' but BAM has '{}'",
                iv_meta.name(),
                bam_meta.name()
            ));
        } else if iv_meta.length() != bam_meta.length() {
            errors.push(format!(
                "  contig '{}' at index {i}: IntervalList has length {} but BAM has length {}",
                iv_meta.name(),
                iv_meta.length(),
                bam_meta.length()
            ));
        }
    }

    if !errors.is_empty() {
        bail!(
            "IntervalList header in {} does not match BAM sequence dictionary:\n{}",
            path.display(),
            errors.join("\n")
        );
    }

    Ok(())
}

/// Validate a single interval's coordinates against the header dictionary,
/// appending any errors to the `errors` vec rather than failing immediately.
fn validate_interval(
    contig: &str,
    start: u32,
    end: u32,
    line_num: usize,
    header_dict: &SequenceDictionary,
    errors: &mut Vec<String>,
) {
    let Some(meta) = header_dict.get_by_name(contig) else {
        errors.push(format!("line {line_num}: contig '{contig}' not found in IntervalList header"));
        return;
    };

    // start and end are already 0-based half-open at this point.
    // Recover the original 1-based values for error messages.
    let start1 = u64::from(start) + 1;
    let end1 = u64::from(end);

    if start1 < 1 {
        errors.push(format!("line {line_num}: start position must be >= 1, got {start1}"));
    }

    let seq_len = meta.length() as u64;
    if end1 > seq_len {
        errors.push(format!(
            "line {line_num}: end position {end1} exceeds length of contig '{}' ({seq_len})",
            meta.name()
        ));
    }
}

/// Parse one IntervalList data line (tab-separated: contig, start1, end1, strand, name).
/// Returns `(contig, start_0based, end_exclusive, optional_name)`.
fn parse_interval_list_line(
    line: &str,
    line_num: usize,
) -> Result<(String, u32, u32, Option<String>)> {
    let mut cols = line.splitn(6, '\t');
    let contig = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("IntervalList line {line_num}: missing contig"))?
        .to_string();
    let start1: u64 = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("IntervalList line {line_num}: missing start"))?
        .parse()
        .with_context(|| format!("IntervalList line {line_num}: bad start"))?;
    let end1: u64 = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("IntervalList line {line_num}: missing end"))?
        .parse()
        .with_context(|| format!("IntervalList line {line_num}: bad end"))?;

    // Column 4 is strand (skip), column 5 is name (optional)
    let _strand = cols.next(); // skip strand
    let name = cols.next().map(|s| s.trim().to_string()).filter(|s| !s.is_empty());

    if start1 < 1 {
        bail!("IntervalList line {line_num}: start position must be >= 1, got {start1}");
    }

    // 1-based inclusive [start1, end1] → 0-based half-open [start1-1, end1)
    let start0 = u32::try_from(start1 - 1)
        .with_context(|| format!("IntervalList line {line_num}: start exceeds u32 range"))?;
    let end = u32::try_from(end1)
        .with_context(|| format!("IntervalList line {line_num}: end exceeds u32 range"))?;
    Ok((contig, start0, end, name))
}

/// Parse one BED line (tab-separated: contig, start, end, name?, ...).
/// Returns `(contig, start_0based, end_exclusive, optional_name)`.
fn parse_bed_line(line: &str, line_num: usize) -> Result<(String, u32, u32, Option<String>)> {
    let mut cols = line.splitn(5, '\t');
    let contig = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("BED line {line_num}: missing contig"))?
        .to_string();
    let start: u64 = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("BED line {line_num}: missing start"))?
        .parse()
        .with_context(|| format!("BED line {line_num}: bad start"))?;
    let end: u64 = cols
        .next()
        .ok_or_else(|| anyhow::anyhow!("BED line {line_num}: missing end"))?
        .parse()
        .with_context(|| format!("BED line {line_num}: bad end"))?;

    // Column 4 is name (optional)
    let name = cols.next().map(|s| s.trim().to_string()).filter(|s| !s.is_empty() && s != ".");

    // BED is already 0-based half-open
    let start = u32::try_from(start)
        .with_context(|| format!("BED line {line_num}: start exceeds u32 range"))?;
    let end = u32::try_from(end)
        .with_context(|| format!("BED line {line_num}: end exceeds u32 range"))?;
    Ok((contig, start, end, name))
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;
    use std::num::NonZeroUsize;

    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    /// Build a header with the given contig names and lengths.
    fn make_header(contigs: &[(&str, usize)]) -> Header {
        let mut builder = Header::builder();
        for &(name, len) in contigs {
            let ref_seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(len).expect("len > 0"));
            builder = builder.add_reference_sequence(name.as_bytes(), ref_seq);
        }
        builder.build()
    }

    /// Build a `SequenceDictionary` from contig names and lengths.
    fn make_dict(contigs: &[(&str, usize)]) -> SequenceDictionary {
        SequenceDictionary::from(&make_header(contigs))
    }

    /// Write an interval list string to a temp file and load it via `load_intervals`.
    fn load_from_string(content: &str, dict: &SequenceDictionary) -> Result<Vec<Interval>> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.interval_list");
        std::fs::write(&path, content).unwrap();
        load_intervals(&path, dict)
    }

    /// Write bytes to a temp file and load via `load_intervals`.
    fn load_from_bytes(bytes: &[u8], dict: &SequenceDictionary) -> Result<Vec<Interval>> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.interval_list");
        std::fs::write(&path, bytes).unwrap();
        load_intervals(&path, dict)
    }

    /// Helper to build an `Intervals` directly from raw data for testing.
    fn intervals_from_raw(
        n_contigs: usize,
        contig_lengths: &[usize],
        raw: &[(usize, u32, u32)],
    ) -> Intervals {
        let mut builder = Header::builder();
        for (i, &len) in contig_lengths.iter().enumerate().take(n_contigs) {
            let name = format!("seq{i}");
            let ref_seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(len).expect("len > 0"));
            builder = builder.add_reference_sequence(name.as_bytes(), ref_seq);
        }
        let header = builder.build();
        let dict = SequenceDictionary::from(&header);

        let raw_count = raw.len();
        let mut by_contig: Vec<Vec<Interval>> = (0..n_contigs).map(|_| Vec::new()).collect();
        for &(ref_id, start, end) in raw {
            if ref_id < n_contigs {
                let name = format!("seq{ref_id}:{start}-{end}");
                by_contig[ref_id].push(Interval::new(ref_id, start, end, name));
            }
        }
        for ivs in &mut by_contig {
            ivs.sort_unstable_by_key(|iv| (iv.start, iv.end));
        }
        Intervals { by_contig, dict, raw_count }
    }

    // ── Interval struct tests ───────────────────────────────────────────────

    #[test]
    fn test_interval_name() {
        let iv = Interval::new(0, 100, 200, "my_interval".to_string());
        assert_eq!(iv.name(), "my_interval");
    }

    #[test]
    fn test_interval_with_contig_name() {
        let iv = Interval::with_contig_name(0, 100, 200, "chr1");
        assert_eq!(iv.name(), "chr1:100-200");
    }

    #[test]
    fn test_interval_len() {
        let iv = Interval::new(0, 100, 200, "test".to_string());
        assert_eq!(iv.len(), 100);
        assert!(!iv.is_empty());
    }

    #[test]
    fn test_interval_empty() {
        let iv = Interval::new(0, 100, 100, "test".to_string());
        assert_eq!(iv.len(), 0);
        assert!(iv.is_empty());
    }

    // ── Intervals collection tests ──────────────────────────────────────────

    #[test]
    fn test_contains_pos_hit() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20), (0, 30, 40)]);
        assert!(ivs.contains_pos(0, 10));
        assert!(ivs.contains_pos(0, 15));
        assert!(ivs.contains_pos(0, 19));
        assert!(ivs.contains_pos(0, 30));
        assert!(ivs.contains_pos(0, 39));
    }

    #[test]
    fn test_contains_pos_miss() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20), (0, 30, 40)]);
        assert!(!ivs.contains_pos(0, 9));
        assert!(!ivs.contains_pos(0, 20)); // exclusive end
        assert!(!ivs.contains_pos(0, 25));
        assert!(!ivs.contains_pos(0, 40));
        assert!(!ivs.contains_pos(0, 100));
    }

    #[test]
    fn test_contains_pos_empty() {
        let ivs = intervals_from_raw(1, &[100], &[]);
        assert!(!ivs.contains_pos(0, 5));
    }

    #[test]
    fn test_contains_pos_bad_ref_id() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20)]);
        assert!(!ivs.contains_pos(99, 10));
    }

    #[test]
    fn test_has_contig() {
        let ivs = intervals_from_raw(3, &[100, 200, 300], &[(0, 10, 20), (2, 5, 15)]);
        assert!(ivs.has_contig(0));
        assert!(!ivs.has_contig(1));
        assert!(ivs.has_contig(2));
        assert!(!ivs.has_contig(99));
    }

    #[test]
    fn test_get_contig() {
        let ivs = intervals_from_raw(2, &[100, 200], &[(0, 10, 20), (0, 30, 40)]);
        let contig_ivs = ivs.get_contig(0);
        assert_eq!(contig_ivs.len(), 2);
        assert_eq!((contig_ivs[0].start, contig_ivs[0].end), (10, 20));
        assert_eq!((contig_ivs[1].start, contig_ivs[1].end), (30, 40));
        assert!(ivs.get_contig(1).is_empty());
        assert!(ivs.get_contig(99).is_empty());
    }

    #[test]
    fn test_iter() {
        let ivs = intervals_from_raw(3, &[100, 200, 300], &[(0, 10, 20), (2, 5, 15), (0, 30, 40)]);
        let all: Vec<_> = ivs.iter().map(|iv| (iv.ref_id, iv.start, iv.end)).collect();
        assert_eq!(all, vec![(0, 10, 20), (0, 30, 40), (2, 5, 15)]);
    }

    #[test]
    fn test_territory() {
        let ivs = intervals_from_raw(2, &[100, 200], &[(0, 10, 20), (0, 30, 40), (1, 0, 50)]);
        assert_eq!(ivs.territory(), 10 + 10 + 50);
    }

    #[test]
    fn test_raw_count() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20), (0, 30, 40)]);
        assert_eq!(ivs.raw_count(), 2);
    }

    #[test]
    fn test_merged_no_overlap() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20), (0, 30, 40)]);
        let merged = ivs.merged();
        assert_eq!(merged.count(), 2);
        assert_eq!(merged.raw_count(), 2); // raw count preserved
    }

    #[test]
    fn test_merged_with_overlap() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 25), (0, 20, 40)]);
        let merged = ivs.merged();
        let contig = merged.get_contig(0);
        assert_eq!(contig.len(), 1);
        assert_eq!(contig[0].start, 10);
        assert_eq!(contig[0].end, 40);
        assert_eq!(merged.raw_count(), 2); // raw count preserved
    }

    #[test]
    fn test_merged_abutting() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20), (0, 20, 30)]);
        let merged = ivs.merged();
        let contig = merged.get_contig(0);
        assert_eq!(contig.len(), 1);
        assert_eq!(contig[0].start, 10);
        assert_eq!(contig[0].end, 30);
    }

    #[test]
    fn test_padded() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 20, 30)]);
        let padded = ivs.padded(5);
        let contig = padded.get_contig(0);
        assert_eq!(contig.len(), 1);
        assert_eq!(contig[0].start, 15);
        assert_eq!(contig[0].end, 35);
    }

    #[test]
    fn test_padded_clamped() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 2, 98)]);
        let padded = ivs.padded(10);
        let contig = padded.get_contig(0);
        assert_eq!(contig[0].start, 0); // clamped to 0
        assert_eq!(contig[0].end, 100); // clamped to contig length
    }

    #[test]
    fn test_contig_bitvec_basic() {
        let ivs = intervals_from_raw(1, &[20], &[(0, 5, 10)]);
        let bv = ivs.contig_bitvec(0);
        assert_eq!(bv.len(), 20);
        for i in 0..20 {
            if (5..10).contains(&i) {
                assert!(bv[i], "expected set at {i}");
            } else {
                assert!(!bv[i], "expected unset at {i}");
            }
        }
    }

    #[test]
    fn test_contig_bitvec_multiple_intervals() {
        let ivs = intervals_from_raw(1, &[30], &[(0, 2, 5), (0, 10, 15)]);
        let bv = ivs.contig_bitvec(0);
        assert_eq!(bv.len(), 30);
        let set: Vec<usize> = (0..30).filter(|&i| bv[i]).collect();
        assert_eq!(set, vec![2, 3, 4, 10, 11, 12, 13, 14]);
    }

    #[test]
    fn test_contig_bitvec_empty() {
        let ivs = intervals_from_raw(1, &[100], &[]);
        let bv = ivs.contig_bitvec(0);
        assert!(bv.is_empty());
    }

    #[test]
    fn test_contig_bitvec_out_of_range() {
        let ivs = intervals_from_raw(1, &[100], &[(0, 10, 20)]);
        let bv = ivs.contig_bitvec(99);
        assert!(bv.is_empty());
    }

    #[test]
    fn test_contig_bitvec_clamped_to_length() {
        // Interval extends beyond contig length — should be clamped.
        let ivs = intervals_from_raw(1, &[10], &[(0, 5, 20)]);
        let bv = ivs.contig_bitvec(0);
        assert_eq!(bv.len(), 10);
        assert_eq!(bv.count_ones(), 5); // positions 5..10
    }

    // ── Parse function tests ────────────────────────────────────────────────

    #[test]
    fn test_parse_interval_list_line_valid() {
        let (contig, start, end, name) =
            parse_interval_list_line("chr1\t100\t200\t+\tgene1", 1).unwrap();
        assert_eq!(contig, "chr1");
        assert_eq!(start, 99); // 1-based → 0-based
        assert_eq!(end, 200);
        assert_eq!(name.as_deref(), Some("gene1"));
    }

    #[test]
    fn test_parse_interval_list_line_no_name() {
        let (contig, start, end, name) = parse_interval_list_line("chr1\t100\t200\t+", 1).unwrap();
        assert_eq!(contig, "chr1");
        assert_eq!(start, 99);
        assert_eq!(end, 200);
        assert!(name.is_none());
    }

    #[test]
    fn test_parse_interval_list_line_missing_end() {
        let result = parse_interval_list_line("chr1\t100", 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_interval_list_line_bad_number() {
        let result = parse_interval_list_line("chr1\tabc\t200\t+\tg", 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_interval_list_line_start_zero() {
        let result = parse_interval_list_line("chr1\t0\t200\t+\tname", 1);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("must be >= 1"));
    }

    #[test]
    fn test_parse_bed_line_valid() {
        let (contig, start, end, name) = parse_bed_line("chr1\t100\t200\tmy_region", 1).unwrap();
        assert_eq!(contig, "chr1");
        assert_eq!(start, 100); // BED is already 0-based
        assert_eq!(end, 200);
        assert_eq!(name.as_deref(), Some("my_region"));
    }

    #[test]
    fn test_parse_bed_line_no_name() {
        let (contig, start, end, name) = parse_bed_line("chr1\t100\t200", 1).unwrap();
        assert_eq!(contig, "chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert!(name.is_none());
    }

    #[test]
    fn test_parse_bed_line_dot_name() {
        let (_, _, _, name) = parse_bed_line("chr1\t100\t200\t.", 1).unwrap();
        assert!(name.is_none());
    }

    #[test]
    fn test_parse_bed_line_missing_end() {
        let result = parse_bed_line("chr1\t100", 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_bed_line_bad_number() {
        let result = parse_bed_line("chr1\t100\txyz", 1);
        assert!(result.is_err());
    }

    // ── Format detection tests ──────────────────────────────────────────────

    #[test]
    fn test_detect_interval_list_format() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        let ivs = result.unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0].start, 99); // 1-based → 0-based confirms IntervalList path
    }

    #[test]
    fn test_detect_bed_format() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "chr1\t100\t200\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        let ivs = result.unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0].start, 100); // stays 0-based confirms BED path
    }

    #[test]
    fn test_detect_format_skips_leading_blank_lines_for_bed() {
        // A BED file with leading blank lines should still be detected as BED.
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "\n\nchr1\t100\t200\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        assert_eq!(result.unwrap()[0].start, 100); // BED
    }

    // ── Header validation tests ─────────────────────────────────────────────

    #[test]
    fn test_header_valid_prefix() {
        // Interval list has chr1; BAM has chr1, chr2 — valid prefix.
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
    }

    #[test]
    fn test_header_exact_match() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
    }

    #[test]
    fn test_header_more_contigs_than_bam() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("2 contigs"));
    }

    #[test]
    fn test_header_name_mismatch() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chrX\tLN:2000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("chrX"), "error should mention mismatched name: {err}");
        assert!(err.contains("chr2"), "error should mention BAM name: {err}");
    }

    #[test]
    fn test_header_length_mismatch() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:999\n@SQ\tSN:chr2\tLN:2000\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("999"), "error should mention IL length: {err}");
        assert!(err.contains("1000"), "error should mention BAM length: {err}");
    }

    #[test]
    fn test_header_multiple_mismatches_reported() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:999\n@SQ\tSN:chr2\tLN:1999\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        // Both mismatches should be reported.
        assert!(err.contains("chr1"), "error should mention chr1: {err}");
        assert!(err.contains("chr2"), "error should mention chr2: {err}");
    }

    #[test]
    fn test_header_no_sq_lines() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@HD\tVN:1.6\nchr1\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("no @SQ lines"));
    }

    #[test]
    fn test_header_only_no_intervals() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    // ── Individual interval validation tests ────────────────────────────────

    #[test]
    fn test_interval_end_exceeds_contig_length() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t1001\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("exceeds length"));
    }

    #[test]
    fn test_interval_end_at_contig_length() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t1000\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
    }

    #[test]
    fn test_interval_unknown_contig_in_interval_list() {
        // Contig is in BAM dict but not in the interval list header — should error.
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr2\t100\t200\t+\tname\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found in IntervalList header"));
    }

    #[test]
    fn test_multiple_interval_errors_reported() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content =
            "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t1001\t+\tname1\nchr1\t200\t1002\t+\tname2\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("1001"), "should report first error: {err}");
        assert!(err.contains("1002"), "should report second error: {err}");
    }

    // ── Gzip/bgzip decompression tests ──────────────────────────────────────

    #[test]
    fn test_load_gzipped_interval_list() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = b"@SQ\tSN:chr1\tLN:1000\nchr1\t100\t200\t+\tname\n";

        // Compress with flate2.
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(content).unwrap();
        let compressed = encoder.finish().unwrap();

        let result = load_from_bytes(&compressed, &dict);
        assert!(result.is_ok());
        let ivs = result.unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0].start, 99);
    }

    #[test]
    fn test_load_gzipped_bed() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = b"chr1\t100\t200\tname\n";

        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(content).unwrap();
        let compressed = encoder.finish().unwrap();

        let result = load_from_bytes(&compressed, &dict);
        assert!(result.is_ok());
        let ivs = result.unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0].start, 100); // BED, no conversion
    }

    // ── BED format tests (unchanged behavior) ───────────────────────────────

    #[test]
    fn test_bed_unknown_contig_is_warning_not_error() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "chr1\t100\t200\tknown\nchrX\t100\t200\tunknown\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        // Only the known contig should be loaded.
        assert_eq!(result.unwrap().len(), 1);
    }

    #[test]
    fn test_bed_with_blank_lines() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "\nchr1\t100\t200\tname\n\nchr1\t300\t400\tname2\n\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 2);
    }

    // ── IntervalList blank line handling ─────────────────────────────────────

    #[test]
    fn test_interval_list_blank_lines_between_data() {
        let dict = make_dict(&[("chr1", 1000)]);
        let content = "@SQ\tSN:chr1\tLN:1000\nchr1\t100\t200\t+\ta\n\nchr1\t300\t400\t+\tb\n";
        let result = load_from_string(content, &dict);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 2);
    }
}
