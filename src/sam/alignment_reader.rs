use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, cram, fasta, sam};
use noodles_bgzf::io::Reader as BgzfReader;

/// A unified alignment reader that handles SAM, BAM, CRAM, and gzipped SAM files.
pub enum AlignmentReader {
    Sam(sam::io::Reader<BufReader<File>>),
    GzippedSam(Box<sam::io::Reader<BufReader<GzDecoder<File>>>>),
    Bam(bam::io::Reader<BgzfReader<File>>),
    Cram(cram::io::Reader<File>),
}

impl AlignmentReader {
    /// Open a SAM, BAM, or CRAM file and read its header.
    ///
    /// The format is detected from the file extension (`.sam`, `.sam.gz`,
    /// `.bam`, or `.cram`). A reference FASTA is required for CRAM files;
    /// it is ignored for other formats.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened, the header cannot be read,
    /// or a CRAM file is opened without a reference.
    pub fn new(path: &Path, reference: Option<&Path>) -> Result<(Self, Header)> {
        let format = detect_format(path)?;
        let context = || format!("Failed to open alignment file: {}", path.display());
        let header_context = || format!("Failed to read header from: {}", path.display());

        match (format, reference) {
            (AlignmentFormat::Sam, _) => {
                let file = File::open(path).with_context(context)?;
                let mut reader = sam::io::Reader::new(BufReader::new(file));
                let header = reader.read_header().with_context(header_context)?;
                Ok((Self::Sam(reader), header))
            }
            (AlignmentFormat::GzippedSam, _) => {
                let file = File::open(path).with_context(context)?;
                let mut reader = sam::io::Reader::new(BufReader::new(GzDecoder::new(file)));
                let header = reader.read_header().with_context(header_context)?;
                Ok((Self::GzippedSam(Box::new(reader)), header))
            }
            (AlignmentFormat::Bam, _) => {
                let file = File::open(path).with_context(context)?;
                let mut reader = bam::io::Reader::new(file);
                let header = reader.read_header().with_context(header_context)?;
                Ok((Self::Bam(reader), header))
            }
            (AlignmentFormat::Cram, Some(ref_path)) => {
                let repo = build_repository(ref_path)?;
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_path(path)
                    .with_context(context)?;
                let header = reader.read_header().with_context(header_context)?;
                Ok((Self::Cram(reader), header))
            }
            (AlignmentFormat::Cram, None) => {
                bail!("CRAM files require a reference FASTA (--reference): {}", path.display());
            }
        }
    }

    /// Returns an iterator that reads alignment records as owned `RecordBuf`s.
    pub fn record_bufs<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> Box<dyn Iterator<Item = Result<RecordBuf>> + 'a> {
        match self {
            Self::Bam(r) => Box::new(
                r.record_bufs(header).map(|result| result.context("Failed to read BAM record")),
            ),
            Self::Sam(r) => Box::new(
                r.record_bufs(header).map(|result| result.context("Failed to read SAM record")),
            ),
            Self::GzippedSam(r) => Box::new(
                r.record_bufs(header).map(|result| result.context("Failed to read SAM record")),
            ),
            Self::Cram(r) => Box::new(r.records(header).map(|result| {
                let record = result.context("Failed to read CRAM record")?;
                RecordBuf::try_from_alignment_record(header, &record)
                    .context("Failed to convert CRAM record to RecordBuf")
            })),
        }
    }
}

/// Build a `fasta::Repository` from an indexed FASTA path for CRAM decoding.
pub(crate) fn build_repository(ref_path: &Path) -> Result<fasta::Repository> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(ref_path)
        .with_context(|| format!("Failed to open reference FASTA: {}", ref_path.display()))?;

    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)))
}

/// Detected alignment file format based on extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AlignmentFormat {
    Sam,
    GzippedSam,
    Bam,
    Cram,
}

/// Detect the alignment format from a file extension.
///
/// Recognizes `.sam`, `.sam.gz`, `.bam`, and `.cram` (case-insensitive).
fn detect_format(path: &Path) -> Result<AlignmentFormat> {
    let ext = path.extension().and_then(|e| e.to_str());

    // Check for .sam.gz: extension is "gz" and the stem ends with ".sam"
    if ext.is_some_and(|e| e.eq_ignore_ascii_case("gz")) {
        let stem = path.file_stem().unwrap_or_default();
        if Path::new(stem).extension().is_some_and(|e| e.eq_ignore_ascii_case("sam")) {
            return Ok(AlignmentFormat::GzippedSam);
        }
    }

    match ext {
        Some(e) if e.eq_ignore_ascii_case("sam") => Ok(AlignmentFormat::Sam),
        Some(e) if e.eq_ignore_ascii_case("bam") => Ok(AlignmentFormat::Bam),
        Some(e) if e.eq_ignore_ascii_case("cram") => Ok(AlignmentFormat::Cram),
        _ => bail!(
            "Cannot determine alignment format from extension. \
             Expected .sam, .sam.gz, .bam, or .cram: {}",
            path.display()
        ),
    }
}
