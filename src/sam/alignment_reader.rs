//! SAM/BAM/CRAM readers — sequential and indexed flavors. Both live in
//! one module because they share the same record-shaping API
//! ([`fill_record`](AlignmentReader::fill_record),
//! [`riker_records`](AlignmentReader::riker_records)) and the same
//! `RikerRecordRequirements`-driven decode logic.
//!
//! ## Two types, not one
//!
//! Originally these were a single enum with an indexed/sequential
//! variant per format. That broke down on the parallel reader path:
//! noodles' `bam::io::IndexedReader` stores its index as
//! `Box<dyn BinningIndex>`, which is not `Send`, so a single enum that
//! could ever hold an indexed variant cannot be sent to a worker thread
//! via `std::thread::scope::spawn`. Splitting keeps
//! [`AlignmentReader`] unconditionally `Send` for the parallel pipeline
//! while still letting [`IndexedAlignmentReader`] expose region queries
//! to standalone callers (currently just `error`).
//!
//! ## Reading paths
//!
//! Both types expose:
//!
//! - [`fill_record`](AlignmentReader::fill_record) — read in place into a
//!   pre-allocated [`RikerRecord`] slot. Works for BAM and SAM (any of
//!   their indexed/streaming variants); errors on CRAM, which noodles
//!   does not let us read in place.
//! - [`riker_records`](AlignmentReader::riker_records) — owned-record
//!   iterator. Works for every format. Use for CRAM or low-volume
//!   callers.
//!
//! [`IndexedAlignmentReader`] additionally has [`query`](IndexedAlignmentReader::query)
//! for region lookups.
//!
//! Both consult a [`RikerRecordRequirements`] passed by the caller. On
//! BAM (where aux decode is lazy) the requirements gate which decoders
//! run; on SAM and CRAM noodles decodes eagerly so requirements are
//! informational.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles::core::Region;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, cram, fasta, sam};
use noodles_bgzf as bgzf;
use noodles_bgzf::io::Reader as BgzfReader;

use super::riker_record::{
    AuxTagRequirements, BamRec, FallbackRec, RikerRecord, RikerRecordRequirements,
};

// ─── AlignmentReader (sequential) ───────────────────────────────────────────

/// Sequential SAM/BAM/CRAM reader. `Send`-safe; used by every command
/// that streams the input front-to-back. Multi's parallel reader thread
/// requires this type because the alternative
/// ([`IndexedAlignmentReader`]) holds a non-`Send` boxed index.
///
/// Open via [`open`](Self::open). For region queries on indexed BAM/CRAM
/// inputs, use [`IndexedAlignmentReader::open`] instead.
pub struct AlignmentReader {
    inner: Inner,
    header: Header,
}

/// Format-specific reader state for [`AlignmentReader`]. One variant per
/// supported sequential format.
enum Inner {
    Sam(sam::io::Reader<BufReader<File>>),
    GzippedSam(Box<sam::io::Reader<BufReader<MultiGzDecoder<File>>>>),
    Bam(bam::io::Reader<BgzfReader<File>>),
    Cram(cram::io::Reader<File>),
}

impl AlignmentReader {
    /// Open `path` for streaming reads. Format is detected from the
    /// extension (`.sam`, `.sam.gz`, `.bam`, `.cram`). A reference FASTA
    /// is required for CRAM and ignored for the others.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened, the header cannot
    /// be read, or a CRAM file is opened without a reference.
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self> {
        match detect_format(path)? {
            AlignmentFormat::Sam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let mut reader = sam::io::Reader::new(BufReader::new(file));
                let header = reader.read_header().with_context(|| header_context(path))?;
                Ok(Self { inner: Inner::Sam(reader), header })
            }
            AlignmentFormat::GzippedSam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let mut reader = sam::io::Reader::new(BufReader::new(MultiGzDecoder::new(file)));
                let header = reader.read_header().with_context(|| header_context(path))?;
                Ok(Self { inner: Inner::GzippedSam(Box::new(reader)), header })
            }
            AlignmentFormat::Bam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let mut reader = bam::io::Reader::new(file);
                let header = reader.read_header().with_context(|| header_context(path))?;
                Ok(Self { inner: Inner::Bam(reader), header })
            }
            AlignmentFormat::Cram => match reference {
                Some(ref_path) => {
                    let repo = build_repository(ref_path)?;
                    let mut reader = cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(repo)
                        .build_from_path(path)
                        .with_context(|| open_context(path))?;
                    let header = reader.read_header().with_context(|| header_context(path))?;
                    Ok(Self { inner: Inner::Cram(reader), header })
                }
                None => {
                    bail!("CRAM files require a reference FASTA (--reference): {}", path.display())
                }
            },
        }
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// True iff [`fill_record`](Self::fill_record) can populate a slot
    /// without allocating. False only for CRAM, which noodles does not
    /// let us read in place.
    #[must_use]
    pub fn supports_in_place_reads(&self) -> bool {
        !matches!(self.inner, Inner::Cram(_))
    }

    /// Construct an empty [`RikerRecord`] of the variant this reader
    /// fills on its in-place path. BAM hands back [`RikerRecord::Bam`];
    /// every other format gets [`RikerRecord::Fallback`]. Used by
    /// callers that pre-allocate a pool of slots before reading.
    #[must_use]
    pub fn empty_record(&self) -> RikerRecord {
        match self.inner {
            Inner::Bam(_) => RikerRecord::Bam(BamRec::new()),
            Inner::Sam(_) | Inner::GzippedSam(_) | Inner::Cram(_) => {
                RikerRecord::Fallback(FallbackRec::empty())
            }
        }
    }

    /// Read the next record into `slot` in place, decoding only the
    /// fields named in `requirements`. Returns `Ok(true)` on success,
    /// `Ok(false)` at EOF.
    ///
    /// `slot`'s variant must match this reader's format — BAM readers
    /// fill [`RikerRecord::Bam`]; SAM readers fill
    /// [`RikerRecord::Fallback`]. Use [`empty_record`](Self::empty_record)
    /// to pre-allocate a slot of the correct variant.
    ///
    /// CRAM is not supported on this path; use the
    /// [`riker_records`](Self::riker_records) iterator instead.
    ///
    /// # Errors
    /// Returns an error if the underlying read fails, if the slot's
    /// variant doesn't match the reader's format, or if the reader is
    /// CRAM.
    pub fn fill_record(
        &mut self,
        requirements: &RikerRecordRequirements,
        slot: &mut RikerRecord,
    ) -> Result<bool> {
        match (&mut self.inner, slot) {
            (Inner::Bam(reader), RikerRecord::Bam(bam_rec)) => {
                fill_bam_slot(reader, bam_rec, requirements)
            }
            (Inner::Sam(reader), RikerRecord::Fallback(fb)) => {
                fill_sam_slot(reader, &self.header, fb)
            }
            (Inner::GzippedSam(reader), RikerRecord::Fallback(fb)) => {
                fill_sam_slot(reader, &self.header, fb)
            }
            (Inner::Cram(_), _) => {
                bail!("fill_record does not support CRAM; iterate riker_records() instead")
            }
            (Inner::Bam(_), _) => bail!("BAM reader requires a RikerRecord::Bam slot"),
            (Inner::Sam(_) | Inner::GzippedSam(_), _) => {
                bail!("SAM reader requires a RikerRecord::Fallback slot")
            }
        }
    }

    /// Returns an iterator that yields one [`RikerRecord`] per record in
    /// the file, allocating a fresh [`BamRec`] / [`FallbackRec`] per
    /// yield. Works for every format. Use
    /// [`fill_record`](Self::fill_record) for BAM/SAM hot paths instead
    /// — `fill_record` reuses one slot across the whole loop. This
    /// iterator is the right tool for CRAM (where in-place reads aren't
    /// available) or low-volume callers where per-record allocation
    /// isn't the bottleneck.
    pub fn riker_records<'a>(
        &'a mut self,
        requirements: &'a RikerRecordRequirements,
    ) -> Box<dyn Iterator<Item = Result<RikerRecord>> + 'a> {
        let header = &self.header;
        match &mut self.inner {
            Inner::Bam(reader) => Box::new(BamRikerRecordIter { reader, requirements }),
            Inner::Sam(reader) => Box::new(reader.record_bufs(header).map(|result| {
                let buf = result.context("Failed to read SAM record")?;
                Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
            })),
            Inner::GzippedSam(reader) => Box::new(reader.record_bufs(header).map(|result| {
                let buf = result.context("Failed to read SAM record")?;
                Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
            })),
            Inner::Cram(reader) => Box::new(reader.records(header).map(move |result| {
                let cram_rec = result.context("Failed to read CRAM record")?;
                let buf = RecordBuf::try_from_alignment_record(header, &cram_rec)
                    .context("Failed to convert CRAM record to RecordBuf")?;
                Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
            })),
        }
    }
}

// ─── IndexedAlignmentReader ─────────────────────────────────────────────────

/// Indexed SAM/BAM/CRAM reader supporting region queries. Not `Send` —
/// noodles stores the BAM index as `Box<dyn BinningIndex>`, which lacks
/// a `Send` bound. Use [`AlignmentReader`] for sequential / parallel
/// pipelines; this type is reserved for callers that need
/// [`query`](Self::query) (currently just `error` standalone).
///
/// Open via [`open`](Self::open). Only BAM (`.bai`/`.csi`) and CRAM
/// (`.crai`) are supported; SAM has no widely-deployed index format.
pub struct IndexedAlignmentReader {
    inner: IndexedInner,
    header: Header,
}

/// Format-specific reader state for [`IndexedAlignmentReader`].
enum IndexedInner {
    Bam(bam::io::IndexedReader<bgzf::io::Reader<File>>),
    Cram { reader: cram::io::Reader<File>, index: cram::crai::Index },
}

impl IndexedAlignmentReader {
    /// Open `path` with its index, enabling region queries via
    /// [`query`](Self::query). A reference FASTA is required for CRAM
    /// and ignored for BAM.
    ///
    /// # Errors
    /// Returns an error if the file or its index cannot be opened, the
    /// format is not BAM/CRAM, or a CRAM file is opened without a
    /// reference.
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self> {
        match detect_format(path)? {
            AlignmentFormat::Sam | AlignmentFormat::GzippedSam => {
                bail!("SAM files cannot be indexed; use a BAM or CRAM file: {}", path.display())
            }
            AlignmentFormat::Bam => {
                validate_bam_index_exists(path)?;
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(path)
                    .with_context(|| format!("Failed to open indexed BAM: {}", path.display()))?;
                let header = reader.read_header().with_context(|| header_context(path))?;
                Ok(Self { inner: IndexedInner::Bam(reader), header })
            }
            AlignmentFormat::Cram => match reference {
                Some(ref_path) => {
                    validate_cram_index_exists(path)?;
                    let repo = build_repository(ref_path)?;
                    let crai_path = append_extension(path, ".crai");
                    let index = cram::crai::fs::read(&crai_path).with_context(|| {
                        format!("Failed to read CRAM index: {}", crai_path.display())
                    })?;
                    let mut reader = cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(repo)
                        .build_from_path(path)
                        .with_context(|| {
                            format!("Failed to open indexed CRAM: {}", path.display())
                        })?;
                    let header = reader.read_header().with_context(|| header_context(path))?;
                    Ok(Self { inner: IndexedInner::Cram { reader, index }, header })
                }
                None => {
                    bail!("CRAM files require a reference FASTA (--reference): {}", path.display())
                }
            },
        }
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// True iff [`fill_record`](Self::fill_record) can populate a slot
    /// without allocating. False for CRAM.
    #[must_use]
    pub fn supports_in_place_reads(&self) -> bool {
        !matches!(self.inner, IndexedInner::Cram { .. })
    }

    /// See [`AlignmentReader::empty_record`].
    #[must_use]
    pub fn empty_record(&self) -> RikerRecord {
        match self.inner {
            IndexedInner::Bam(_) => RikerRecord::Bam(BamRec::new()),
            IndexedInner::Cram { .. } => RikerRecord::Fallback(FallbackRec::empty()),
        }
    }

    /// See [`AlignmentReader::fill_record`].
    ///
    /// # Errors
    /// Returns an error if the underlying read fails, if the slot's
    /// variant doesn't match the reader's format, or if the reader is
    /// CRAM.
    pub fn fill_record(
        &mut self,
        requirements: &RikerRecordRequirements,
        slot: &mut RikerRecord,
    ) -> Result<bool> {
        match (&mut self.inner, slot) {
            (IndexedInner::Bam(reader), RikerRecord::Bam(bam_rec)) => {
                let n = bam_rec.read_from_indexed(reader)?;
                if n == 0 {
                    return Ok(false);
                }
                apply_requirements(bam_rec, requirements)?;
                Ok(true)
            }
            (IndexedInner::Cram { .. }, _) => {
                bail!("fill_record does not support CRAM; iterate riker_records() instead")
            }
            (IndexedInner::Bam(_), _) => bail!("BAM reader requires a RikerRecord::Bam slot"),
        }
    }

    /// See [`AlignmentReader::riker_records`].
    pub fn riker_records<'a>(
        &'a mut self,
        requirements: &'a RikerRecordRequirements,
    ) -> Box<dyn Iterator<Item = Result<RikerRecord>> + 'a> {
        let header = &self.header;
        match &mut self.inner {
            IndexedInner::Bam(reader) => {
                Box::new(IndexedBamRikerRecordIter { reader, requirements })
            }
            IndexedInner::Cram { reader, .. } => {
                Box::new(reader.records(header).map(move |result| {
                    let cram_rec = result.context("Failed to read CRAM record")?;
                    let buf = RecordBuf::try_from_alignment_record(header, &cram_rec)
                        .context("Failed to convert CRAM record to RecordBuf")?;
                    Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
                }))
            }
        }
    }

    /// Stream records overlapping `region` to the callback `f`,
    /// applying `requirements` per record on the BAM fast path. CRAM
    /// records arrive as [`RikerRecord::Fallback`] with everything
    /// decoded.
    ///
    /// On BAM, a single scratch [`BamRec`] is reused across the whole
    /// region — only the inner `bam::Record` allocations come from
    /// noodles' query iterator, while the cigar / quality / aux mirror
    /// buffers are recycled. CRAM still allocates one
    /// [`FallbackRec`] per record because noodles owns the CRAM
    /// decoder's output.
    ///
    /// Callbacks return `Result<()>` so the caller can propagate
    /// errors and abort the scan early.
    ///
    /// # Errors
    /// Returns an error if the underlying query, decode, or callback
    /// fails. The callback's error short-circuits the scan.
    pub fn query_for_each<F>(
        &mut self,
        region: &Region,
        requirements: &RikerRecordRequirements,
        mut f: F,
    ) -> Result<()>
    where
        F: FnMut(&RikerRecord) -> Result<()>,
    {
        let header = &self.header;
        match &mut self.inner {
            IndexedInner::Bam(reader) => {
                let query = reader
                    .query(header, region)
                    .with_context(|| format!("Failed to query BAM for region: {region}"))?;
                // One scratch slot for the whole region — install the
                // next raw record, refresh, then hand a borrow to the
                // callback. Cigar / quality / aux Vecs in the BamRec
                // are reused across iterations.
                let mut record = RikerRecord::Bam(BamRec::new());
                for result in query.records() {
                    let raw = result.context("Failed to read BAM record during query")?;
                    let RikerRecord::Bam(ref mut bam_rec) = record else {
                        unreachable!("scratch was constructed as RikerRecord::Bam");
                    };
                    bam_rec.install(raw)?;
                    apply_requirements(bam_rec, requirements)?;
                    f(&record)?;
                }
            }
            IndexedInner::Cram { reader, index } => {
                let query = reader
                    .query(header, index, region)
                    .with_context(|| format!("Failed to query CRAM for region: {region}"))?;
                for result in query {
                    let cram_rec = result.context("Failed to read CRAM record during query")?;
                    let buf = RecordBuf::try_from_alignment_record(header, &cram_rec)
                        .context("Failed to convert CRAM record to RecordBuf")?;
                    let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));
                    f(&record)?;
                }
            }
        }
        Ok(())
    }
}

// ─── Module-level helpers ───────────────────────────────────────────────────

/// BAM-side iterator yielding [`RikerRecord::Bam`] entries from a
/// sequential reader, with the caller's requirements applied.
struct BamRikerRecordIter<'a, R: Read> {
    reader: &'a mut bam::io::Reader<R>,
    requirements: &'a RikerRecordRequirements,
}

impl<R: Read> Iterator for BamRikerRecordIter<'_, R> {
    type Item = Result<RikerRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bam_rec = BamRec::new();
        match bam_rec.read_from(self.reader) {
            Ok(0) => None,
            Ok(_) => match apply_requirements(&mut bam_rec, self.requirements) {
                Ok(()) => Some(Ok(RikerRecord::Bam(bam_rec))),
                Err(e) => Some(Err(e)),
            },
            Err(e) => Some(Err(e)),
        }
    }
}

/// BAM-side iterator yielding [`RikerRecord::Bam`] entries from an
/// indexed reader. Same shape as [`BamRikerRecordIter`] but driven by
/// `IndexedReader::read_record` (which yields a raw `bam::Record` we
/// install into a fresh [`BamRec`]).
struct IndexedBamRikerRecordIter<'a, R: Read> {
    reader: &'a mut bam::io::IndexedReader<R>,
    requirements: &'a RikerRecordRequirements,
}

impl<R: Read> Iterator for IndexedBamRikerRecordIter<'_, R> {
    type Item = Result<RikerRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bam_rec = BamRec::new();
        match bam_rec.read_from_indexed(self.reader) {
            Ok(0) => None,
            Ok(_) => match apply_requirements(&mut bam_rec, self.requirements) {
                Ok(()) => Some(Ok(RikerRecord::Bam(bam_rec))),
                Err(e) => Some(Err(e)),
            },
            Err(e) => Some(Err(e)),
        }
    }
}

/// Read one BAM record into `bam_rec` from a sequential reader and
/// apply `requirements`.
fn fill_bam_slot<R: Read>(
    reader: &mut bam::io::Reader<R>,
    bam_rec: &mut BamRec,
    requirements: &RikerRecordRequirements,
) -> Result<bool> {
    let n = bam_rec.read_from(reader)?;
    if n == 0 {
        return Ok(false);
    }
    apply_requirements(bam_rec, requirements)?;
    Ok(true)
}

/// Read one SAM record into a [`FallbackRec`] in place and refresh the
/// cached scalars. Generic over the SAM reader's underlying byte
/// source so the same body covers plain `BufReader<File>` and the
/// gzip-wrapped variant.
///
/// On `read_record_buf` error the inner `RecordBuf` may be partially
/// overwritten, leaving the slot's cached scalars (e.g. `alignment_end`)
/// describing the previous record — refresh first, then propagate the
/// error, so any later observation of the slot sees self-consistent
/// state rather than half-stale data.
fn fill_sam_slot<R: std::io::BufRead>(
    reader: &mut sam::io::Reader<R>,
    header: &Header,
    fb: &mut FallbackRec,
) -> Result<bool> {
    match reader.read_record_buf(header, fb.inner_mut()) {
        Ok(0) => Ok(false),
        Ok(_) => {
            fb.refresh_cache();
            Ok(true)
        }
        Err(e) => {
            fb.refresh_cache();
            Err(anyhow::Error::from(e).context("Failed to read SAM record"))
        }
    }
}

/// Apply caller-declared decoder requirements to a freshly-read
/// [`BamRec`]: decode sequence bytes and/or scan the requested aux tags.
fn apply_requirements(bam_rec: &mut BamRec, requirements: &RikerRecordRequirements) -> Result<()> {
    if requirements.needs_sequence() {
        bam_rec.decode_sequence();
    }
    match requirements.aux_tags() {
        AuxTagRequirements::None => {}
        AuxTagRequirements::Specific { values, presence } => {
            bam_rec.scan_aux_tags(values, presence)?;
        }
        AuxTagRequirements::All => {
            bam_rec.decode_all_aux()?;
        }
    }
    Ok(())
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

/// Returns `path` with `suffix` (e.g. `".bai"`) appended after the
/// existing extension. Used to derive index-file paths.
fn append_extension(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

/// Confirm that one of `path`.bam.bai / `path`.bai / `path`.csi exists.
fn validate_bam_index_exists(path: &Path) -> Result<()> {
    let bai_alt = path.with_extension("bam.bai");
    let bai = append_extension(path, ".bai");
    let csi = append_extension(path, ".csi");

    if bai.exists() || bai_alt.exists() || csi.exists() {
        return Ok(());
    }

    bail!(
        "BAM index not found. Expected one of:\n  {}\n  {}\n  {}\n\
         Run `samtools index {}` to create one.",
        bai.display(),
        bai_alt.display(),
        csi.display(),
        path.display(),
    );
}

/// Confirm that `path`.crai exists.
fn validate_cram_index_exists(path: &Path) -> Result<()> {
    let crai = append_extension(path, ".crai");
    if crai.exists() {
        return Ok(());
    }
    bail!(
        "CRAM index not found. Expected: {}\n\
         Run `samtools index {}` to create one.",
        crai.display(),
        path.display(),
    );
}

/// Context message for "failed to open" errors.
fn open_context(path: &Path) -> String {
    format!("Failed to open alignment file: {}", path.display())
}

/// Context message for "failed to read header" errors.
fn header_context(path: &Path) -> String {
    format!("Failed to read header from: {}", path.display())
}
