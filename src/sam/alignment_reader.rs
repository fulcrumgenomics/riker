//! SAM/BAM/CRAM readers — sequential and indexed flavors. Both live in
//! one module because they share the same `RikerRecordRequirements`-
//! driven decode logic.
//!
//! ## Two types, not one
//!
//! [`AlignmentReader`] streams front-to-back; [`IndexedAlignmentReader`]
//! answers region queries against an index. They stay separate types because
//! the two access patterns want different backends. Sequential BAM is fastest
//! through noodles (see below), but *indexed* access wants a persistent decode
//! thread pool that survives seeks — which htslib provides and a streaming
//! BGZF `MultithreadedReader` does not (it joins and respawns its threads
//! on every seek, which collapses under a many-small-region `--intervals`
//! workload). So [`IndexedAlignmentReader`] reads both BAM and CRAM through
//! htslib, and is reserved for the standalone callers that need queries
//! (currently just `error`).
//!
//! ## Reading paths
//!
//! [`AlignmentReader`] (sequential, front-to-back) exposes both:
//!
//! - [`fill_record`](AlignmentReader::fill_record) — read in place into a
//!   pre-allocated [`RikerRecord`] slot. Works for every format.
//! - [`riker_records`](AlignmentReader::riker_records) — owned-record
//!   iterator. Works for every format. Use for low-volume callers where
//!   per-record allocation isn't the bottleneck.
//!
//! [`IndexedAlignmentReader`] is region-query-only:
//!
//! - [`query_for_each`](IndexedAlignmentReader::query_for_each) — stream
//!   records overlapping a region to a callback, reusing one scratch
//!   slot for the whole region.
//!
//! Both consult a [`RikerRecordRequirements`] from the caller: given once, before
//! any record is read, to [`AlignmentReader::set_requirements`], or passed per
//! query to [`IndexedAlignmentReader::query_for_each`]. On BAM and CRAM (where
//! aux decode is lazy) the requirements gate which decoders run; on SAM noodles
//! decodes eagerly so requirements are informational. The sequential reader
//! also hands them to htslib, which then skips decoding the CRAM data series no
//! consumer reads (see [`IndexedAlignmentReader`] for why the indexed one does
//! not).
//!
//! ## BAM via noodles + bgzf; CRAM via rust-htslib
//!
//! BAM records parse through noodles. With `decode_threads == 0` it uses
//! noodles' single-threaded BGZF reader (the default, lowest RSS); with a
//! non-zero count it engages the `bgzf` crate's `MultithreadedReader` for
//! parallel BGZF decompression at no meaningful RSS cost. That reader sizes
//! its read-ahead to `worker_count.max(8)`, avoiding the one-worker stall
//! where noodles' own MT reader regressed below single-threaded.
//!
//! CRAM reads through `rust-htslib` rather than `noodles-cram`: benchmarks on
//! a 3.1 GB exome CRAM showed htslib decoding ~3.8× faster than noodles-cram
//! (1.05M vs 273k records/s). htslib also handles `.crai` indexing natively.
//! `decode_threads` sizes htslib's decode pool. SAM always stays on noodles
//! single-threaded and ignores `decode_threads`.
//!
//! htslib's textual SAM header (CRAM) is parsed into a `noodles::sam::Header`
//! at open time so every collector and helper sees the same `Header` type
//! regardless of format.

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles::core::Region;
use noodles::sam::Header;
use noodles::{bam, sam};
use rust_htslib::bam::Read as HtsRead;

use crate::path_util::append_extension;

use super::riker_record::{
    AuxTagRequirements, BamRec, FallbackRec, HtslibRec, RikerRecord, RikerRecordRequirements,
};

/// Capacity of the `BufReader` placed between the file and the BGZF decoder
/// on the sequential BAM path. The BGZF reader issues two `read_exact` calls
/// per block (18-byte header + payload up to ~64 KiB); a 256 KiB buffer covers
/// several blocks per refill and keeps the per-byte syscall overhead negligible.
const BAM_READ_BUFFER_BYTES: usize = 256 * 1024;

/// htslib's `CRAM_OPT_DECODE_MD` default for reading: `-1`, "auto", which
/// regenerates `NM`/`MD` tags the CRAM writer dropped. rust-htslib passes
/// option values as an unsigned `sam_fields`; htslib reads this one back as an
/// `int`, so all bits set is `-1`.
const CRAM_DECODE_MD_AUTO: rust_htslib::htslib::sam_fields = rust_htslib::htslib::sam_fields::MAX;

/// Error for reading an [`AlignmentReader`] before it was told what to decode.
const REQUIREMENTS_NOT_SET: &str =
    "AlignmentReader::set_requirements must be called before reading records";

// ─── AlignmentReader (sequential) ───────────────────────────────────────────

/// Sequential SAM/BAM/CRAM reader. `Send`-safe; used by every command
/// that streams the input front-to-back. Multi's parallel reader thread
/// requires this type because the alternative
/// ([`IndexedAlignmentReader`]) holds a non-`Send` boxed index.
///
/// Open via [`open`](Self::open), declare what the consumers read with
/// [`set_requirements`](Self::set_requirements), then read. For region queries
/// on indexed BAM/CRAM inputs, use [`IndexedAlignmentReader::open`] instead.
pub struct AlignmentReader {
    inner: Inner,
    header: Header,
    /// The fields every record read from here is decoded with. Reading while
    /// this is `None` is an error rather than a default, so a caller that
    /// forgets [`set_requirements`](Self::set_requirements) fails instead of
    /// silently reading records with no bases or aux tags.
    requirements: Option<RikerRecordRequirements>,
    /// Whether a record has been read, after which the requirements are fixed:
    /// htslib's CRAM decode threads may already have decoded containers ahead
    /// of the reader using them.
    started: bool,
}

/// Format-specific reader state for [`AlignmentReader`]. `Bam` is
/// single-threaded noodles; `BamMt` is noodles records over a multithreaded BGZF
/// decoder (`decode_threads > 0`). The `Htslib` arm backs CRAM; it is
/// `Box`-ed because `rust_htslib::bam::Reader` is large (a few hundred bytes
/// of fields, dwarfing the noodles BAM/SAM readers) and boxing keeps the enum
/// small and avoids per-record copies during reader-handle moves.
enum Inner {
    Sam(sam::io::Reader<BufReader<File>>),
    GzippedSam(Box<sam::io::Reader<BufReader<MultiGzDecoder<File>>>>),
    Bam(bam::io::Reader<noodles_bgzf::io::Reader<BufReader<File>>>),
    BamMt(bam::io::Reader<bgzf::MultithreadedReader<BufReader<File>>>),
    Htslib(Box<rust_htslib::bam::Reader>),
}

impl AlignmentReader {
    /// Open `path` for streaming reads. Format is detected from the
    /// extension (`.sam`, `.sam.gz`, `.bam`, `.cram`). A reference FASTA
    /// is required for CRAM and ignored for the others.
    ///
    /// `decode_threads` is the number of BGZF/CRAM decode worker threads to
    /// run *in addition to* the calling thread (derived from `--threads` by
    /// the command). `0` reads single-threaded (the default, lowest RSS). For
    /// BAM, a non-zero count engages the `bgzf` crate's `MultithreadedReader`;
    /// for CRAM it sizes htslib's decode pool. SAM ignores it.
    ///
    /// Call [`set_requirements`](Self::set_requirements) before reading any
    /// record.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened, the header cannot
    /// be read, or a CRAM file is opened without a reference.
    pub fn open(path: &Path, reference: Option<&Path>, decode_threads: usize) -> Result<Self> {
        let (inner, header) = match detect_format(path)? {
            AlignmentFormat::Sam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let mut reader = sam::io::Reader::new(BufReader::new(file));
                let header = reader.read_header().with_context(|| header_context(path))?;
                (Inner::Sam(reader), header)
            }
            AlignmentFormat::GzippedSam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let mut reader = sam::io::Reader::new(BufReader::new(MultiGzDecoder::new(file)));
                let header = reader.read_header().with_context(|| header_context(path))?;
                (Inner::GzippedSam(Box::new(reader)), header)
            }
            AlignmentFormat::Bam => {
                let file = File::open(path).with_context(|| open_context(path))?;
                let buf = BufReader::with_capacity(BAM_READ_BUFFER_BYTES, file);
                match NonZero::new(decode_threads) {
                    // Single-threaded plain BGZF reader (the default).
                    None => {
                        let mut reader = bam::io::Reader::from(noodles_bgzf::io::Reader::new(buf));
                        let header = reader.read_header().with_context(|| header_context(path))?;
                        (Inner::Bam(reader), header)
                    }
                    // Parallel BGZF decompression via the bgzf crate's reader, whose
                    // read-ahead depth is `worker_count.max(8)` (avoids noodles' one-worker stall).
                    Some(workers) => {
                        let mt = bgzf::MultithreadedReader::with_worker_count(workers, buf);
                        let mut reader = bam::io::Reader::from(mt);
                        let header = reader.read_header().with_context(|| header_context(path))?;
                        (Inner::BamMt(reader), header)
                    }
                }
            }
            AlignmentFormat::Cram => match reference {
                Some(ref_path) => {
                    let mut reader = rust_htslib::bam::Reader::from_path(path)
                        .with_context(|| open_context(path))?;
                    reader.set_reference(ref_path).with_context(|| {
                        format!(
                            "Failed to set CRAM reference {} on {}",
                            ref_path.display(),
                            path.display()
                        )
                    })?;
                    set_htslib_threads(&mut reader, decode_threads, path)?;
                    let header = parse_htslib_header(reader.header().as_bytes())
                        .with_context(|| header_context(path))?;
                    (Inner::Htslib(Box::new(reader)), header)
                }
                None => {
                    bail!("CRAM files require a reference FASTA (--reference): {}", path.display())
                }
            },
        };
        Ok(Self { inner, header, requirements: None, started: false })
    }

    /// Declare the fields the consumers will read from every record (for
    /// `multi`, the union across its collectors). Must be called before the
    /// first record is read; the requirements then hold for the rest of the
    /// file.
    ///
    /// On CRAM this also tells htslib which data series to decode, so fields
    /// left out are never decoded at all. A consumer must therefore declare
    /// everything it reads: on BAM and CRAM an undeclared sequence or aux tag
    /// reads as absent. SAM decodes every field regardless, so a SAM-based test
    /// cannot catch a missing declaration.
    ///
    /// # Errors
    /// Returns an error if a record has already been read, or if htslib
    /// rejects the CRAM options.
    pub fn set_requirements(&mut self, requirements: RikerRecordRequirements) -> Result<()> {
        if self.started {
            bail!("record requirements must be set before the first record is read");
        }
        if let Inner::Htslib(reader) = &mut self.inner {
            set_cram_required_fields(reader.as_mut(), &requirements)?;
        }
        self.requirements = Some(requirements);
        Ok(())
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Construct an empty [`RikerRecord`] of the variant this reader
    /// fills on its in-place path. BAM hands back [`RikerRecord::Bam`];
    /// CRAM hands back [`RikerRecord::Htslib`]; SAM and gzipped-SAM hand
    /// back [`RikerRecord::Fallback`]. Used by callers that pre-allocate
    /// a pool of slots before reading.
    #[must_use]
    pub fn empty_record(&self) -> RikerRecord {
        match self.inner {
            Inner::Bam(_) | Inner::BamMt(_) => RikerRecord::Bam(BamRec::new()),
            Inner::Htslib(_) => RikerRecord::Htslib(HtslibRec::new()),
            Inner::Sam(_) | Inner::GzippedSam(_) => RikerRecord::Fallback(FallbackRec::empty()),
        }
    }

    /// Read the next record into `slot` in place, decoding the fields named
    /// by [`set_requirements`](Self::set_requirements). Returns `Ok(true)`
    /// on success, `Ok(false)` at EOF. Works for every format.
    ///
    /// `slot`'s variant must match this reader's format — BAM readers
    /// fill [`RikerRecord::Bam`]; CRAM readers fill
    /// [`RikerRecord::Htslib`]; SAM readers fill
    /// [`RikerRecord::Fallback`]. Use [`empty_record`](Self::empty_record)
    /// to pre-allocate a slot of the correct variant.
    ///
    /// # Errors
    /// Returns an error if [`set_requirements`](Self::set_requirements) has
    /// not been called, if the underlying read fails, or if the slot's
    /// variant doesn't match the reader's format.
    pub fn fill_record(&mut self, slot: &mut RikerRecord) -> Result<bool> {
        let Some(requirements) = &self.requirements else { bail!(REQUIREMENTS_NOT_SET) };
        self.started = true;
        match (&mut self.inner, slot) {
            (Inner::Bam(reader), RikerRecord::Bam(bam_rec)) => {
                fill_bam_slot(reader, bam_rec, requirements)
            }
            (Inner::BamMt(reader), RikerRecord::Bam(bam_rec)) => {
                fill_bam_slot(reader, bam_rec, requirements)
            }
            (Inner::Sam(reader), RikerRecord::Fallback(fb)) => {
                fill_sam_slot(reader, &self.header, fb)
            }
            (Inner::GzippedSam(reader), RikerRecord::Fallback(fb)) => {
                fill_sam_slot(reader, &self.header, fb)
            }
            (Inner::Htslib(reader), RikerRecord::Htslib(hts_rec)) => {
                fill_htslib_slot(reader.as_mut(), hts_rec, requirements)
            }
            (Inner::Bam(_) | Inner::BamMt(_), _) => {
                bail!("BAM reader requires a RikerRecord::Bam slot")
            }
            (Inner::Htslib(_), _) => bail!("htslib reader requires a RikerRecord::Htslib slot"),
            (Inner::Sam(_) | Inner::GzippedSam(_), _) => {
                bail!("SAM reader requires a RikerRecord::Fallback slot")
            }
        }
    }

    /// Returns an iterator that yields one [`RikerRecord`] per record in
    /// the file, decoded per [`set_requirements`](Self::set_requirements)
    /// and allocating a fresh [`BamRec`] / [`FallbackRec`] /
    /// [`HtslibRec`] per yield. Works for every format. Use
    /// [`fill_record`](Self::fill_record) for hot paths instead —
    /// `fill_record` reuses one slot across the whole loop. This
    /// iterator is the right tool for low-volume callers where
    /// per-record allocation isn't the bottleneck.
    ///
    /// # Errors
    /// Returns an error if [`set_requirements`](Self::set_requirements) has
    /// not been called.
    pub fn riker_records(&mut self) -> Result<Box<dyn Iterator<Item = Result<RikerRecord>> + '_>> {
        let Some(requirements) = &self.requirements else { bail!(REQUIREMENTS_NOT_SET) };
        self.started = true;
        let header = &self.header;
        let records: Box<dyn Iterator<Item = Result<RikerRecord>> + '_> = match &mut self.inner {
            Inner::Bam(reader) => Box::new(BamRikerRecordIter { reader, requirements }),
            Inner::BamMt(reader) => Box::new(BamRikerRecordIter { reader, requirements }),
            Inner::Sam(reader) => Box::new(reader.record_bufs(header).map(|result| {
                let buf = result.context("Failed to read SAM record")?;
                Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
            })),
            Inner::GzippedSam(reader) => Box::new(reader.record_bufs(header).map(|result| {
                let buf = result.context("Failed to read SAM record")?;
                Ok(RikerRecord::Fallback(FallbackRec::from_record_buf(buf)))
            })),
            Inner::Htslib(reader) => {
                Box::new(HtslibRikerRecordIter { reader: reader.as_mut(), requirements })
            }
        };
        Ok(records)
    }
}

// ─── IndexedAlignmentReader ─────────────────────────────────────────────────

/// Indexed BAM/CRAM reader supporting region queries, backed by htslib for
/// both formats. htslib owns a persistent decode thread pool that is reused
/// across every query and freed when the reader drops — unlike a
/// per-seek-teardown streaming decoder, so it stays efficient even for the
/// many-small-region (`--intervals`) access pattern. Not `Send`; use
/// [`AlignmentReader`] for the sequential / parallel pipelines. This type is
/// reserved for callers that need
/// [`query_for_each`](Self::query_for_each) (currently just `error`
/// standalone).
///
/// Unlike [`AlignmentReader`], it does not narrow htslib's CRAM decode to the
/// requirements: its one caller, `error`, reads the sequence and `NM`, which
/// need the bases and aux tags decoded anyway.
///
/// Open via [`open`](Self::open). Only BAM (`.bai`/`.csi`) and CRAM
/// (`.crai`) are supported; SAM has no widely-deployed index format.
pub struct IndexedAlignmentReader {
    reader: Box<rust_htslib::bam::IndexedReader>,
    header: Header,
}

impl IndexedAlignmentReader {
    /// Open `path` with its index, enabling region queries via
    /// [`query_for_each`](Self::query_for_each). A reference FASTA is required
    /// for CRAM and ignored for BAM. Both formats read through htslib.
    ///
    /// `decode_threads` sizes htslib's decode thread pool (in addition to the
    /// calling thread); `0` reads single-threaded. The pool is persistent —
    /// reused across every query rather than torn down per seek — so a large
    /// `--intervals` panel doesn't pay repeated thread startup. For BAM the
    /// index is chosen by [`select_bam_index`] (most-recent of the accepted
    /// layouts) and handed to htslib explicitly.
    ///
    /// # Errors
    /// Returns an error if the file or its index cannot be opened, the
    /// format is not BAM/CRAM, or a CRAM file is opened without a
    /// reference.
    pub fn open(path: &Path, reference: Option<&Path>, decode_threads: usize) -> Result<Self> {
        let mut reader = match detect_format(path)? {
            AlignmentFormat::Sam | AlignmentFormat::GzippedSam => {
                bail!("SAM files cannot be indexed; use a BAM or CRAM file: {}", path.display())
            }
            AlignmentFormat::Bam => {
                let index_path = select_bam_index(path)?;
                rust_htslib::bam::IndexedReader::from_path_and_index(path, index_path.as_path())
                    .with_context(|| format!("Failed to open indexed BAM: {}", path.display()))?
            }
            AlignmentFormat::Cram => {
                let Some(ref_path) = reference else {
                    bail!("CRAM files require a reference FASTA (--reference): {}", path.display())
                };
                validate_cram_index_exists(path)?;
                let mut reader = rust_htslib::bam::IndexedReader::from_path(path)
                    .with_context(|| format!("Failed to open indexed CRAM: {}", path.display()))?;
                reader.set_reference(ref_path).with_context(|| {
                    format!(
                        "Failed to set CRAM reference {} on {}",
                        ref_path.display(),
                        path.display()
                    )
                })?;
                reader
            }
        };
        set_htslib_threads(&mut reader, decode_threads, path)?;
        let header = parse_htslib_header(reader.header().as_bytes())
            .with_context(|| header_context(path))?;
        Ok(Self { reader: Box::new(reader), header })
    }

    /// The parsed header.
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Stream records overlapping `region` to the callback `f`,
    /// applying `requirements` per record. A single scratch
    /// [`RikerRecord`] is reused across the whole region — the cigar /
    /// quality / aux Vecs and (for CRAM) the htslib `bam1_t` buffer are
    /// recycled.
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
        fetch_htslib_region(self.reader.as_mut(), region)?;
        let mut record = RikerRecord::Htslib(HtslibRec::new());
        // The let-else inside the loop body looks redundant (scratch is always
        // `Htslib`), but pulling the binding outside the loop would extend the
        // `&mut hts_rec` borrow across `f(&record)` at the bottom and conflict
        // with the immutable borrow there. Re-binding per iteration keeps each
        // borrow scoped tightly enough that both can coexist.
        loop {
            let RikerRecord::Htslib(ref mut hts_rec) = record else {
                unreachable!("scratch was constructed as RikerRecord::Htslib");
            };
            if !hts_rec.read_from(self.reader.as_mut())? {
                break;
            }
            apply_requirements_htslib(hts_rec, requirements)?;
            f(&record)?;
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

/// CRAM-side iterator yielding [`RikerRecord::Htslib`] entries from a
/// rust-htslib reader. Allocates one [`HtslibRec`] per yield; for hot
/// paths use [`AlignmentReader::fill_record`] which reuses a single slot.
struct HtslibRikerRecordIter<'a> {
    reader: &'a mut rust_htslib::bam::Reader,
    requirements: &'a RikerRecordRequirements,
}

impl Iterator for HtslibRikerRecordIter<'_> {
    type Item = Result<RikerRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut hts_rec = HtslibRec::new();
        match hts_rec.read_from(self.reader) {
            Ok(false) => None,
            Ok(true) => match apply_requirements_htslib(&mut hts_rec, self.requirements) {
                Ok(()) => Some(Ok(RikerRecord::Htslib(hts_rec))),
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

/// Size the htslib decode/decompression thread pool on a rust-htslib
/// reader. `threads` is the number of htslib worker threads in addition to
/// the calling thread; `0` leaves htslib single-threaded. The generic bound
/// makes it callable on either reader type.
///
/// # Errors
/// Returns an error if htslib fails to create the thread pool.
fn set_htslib_threads<R: HtsRead>(reader: &mut R, threads: usize, path: &Path) -> Result<()> {
    if threads == 0 {
        return Ok(());
    }
    reader.set_threads(threads).with_context(|| {
        format!("Failed to set {threads} htslib decode threads on {}", path.display())
    })
}

/// The index-file layouts riker accepts for the BAM at `path`, in a fixed
/// order: `<file>.bai` (Picard sibling), `<file>.bam.bai` (samtools suffix),
/// `<file>.bam.csi` (large-contig). Shared by [`select_bam_index`] so
/// selection and the "missing index" error can't drift.
fn bam_index_candidates(path: &Path) -> [PathBuf; 3] {
    [path.with_extension("bai"), append_extension(path, ".bai"), append_extension(path, ".csi")]
}

/// Choose which index to load for the BAM at `path`. When several accepted
/// layouts are present, warn and pick the **most recently modified**, so a
/// freshly re-indexed BAM is used even if a stale index of another form is
/// still lying next to it. Errors (with a `samtools index` hint) if none
/// exist — this doubles as the friendly pre-open check.
fn select_bam_index(path: &Path) -> Result<PathBuf> {
    let present: Vec<(PathBuf, std::time::SystemTime)> = bam_index_candidates(path)
        .into_iter()
        .filter_map(|candidate| {
            std::fs::metadata(&candidate)
                .and_then(|meta| meta.modified())
                .ok()
                .map(|modified| (candidate, modified))
        })
        .collect();
    pick_bam_index(path, present)
}

/// Pick the most-recently-modified index from the present `(path, mtime)`
/// candidates, warning when there is more than one. Split from the filesystem
/// enumeration in [`select_bam_index`] so the selection can be unit-tested with
/// synthetic timestamps.
fn pick_bam_index(
    path: &Path,
    mut present: Vec<(PathBuf, std::time::SystemTime)>,
) -> Result<PathBuf> {
    present.sort_by_key(|(_, modified)| std::cmp::Reverse(*modified));
    match present.as_slice() {
        [] => {
            let [sibling, suffix, csi] = bam_index_candidates(path);
            bail!(
                "BAM index not found. Expected one of:\n  {}\n  {}\n  {}\n\
                 Run `samtools index {}` to create one.",
                suffix.display(),
                sibling.display(),
                csi.display(),
                path.display(),
            )
        }
        [(chosen, _)] => Ok(chosen.clone()),
        [(chosen, _), rest @ ..] => {
            let ignored: Vec<String> = rest.iter().map(|(p, _)| p.display().to_string()).collect();
            log::warn!(
                "Multiple indexes present for {}; using the most recently modified ({}) \
                 and ignoring: {}",
                path.display(),
                chosen.display(),
                ignored.join(", "),
            );
            Ok(chosen.clone())
        }
    }
}

/// Tell htslib to decode only the CRAM data series that `requirements` need
/// (see [`cram_required_fields`]), and to regenerate dropped `NM`/`MD` tags
/// only when a consumer reads them. Must happen before the first record is
/// read: htslib's decode threads read ahead with whatever setting is in force.
///
/// # Errors
/// Returns an error if htslib rejects either option.
fn set_cram_required_fields(
    reader: &mut rust_htslib::bam::Reader,
    requirements: &RikerRecordRequirements,
) -> Result<()> {
    use rust_htslib::htslib as hts;
    reader
        .set_cram_options(
            hts::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            cram_required_fields(requirements),
        )
        .context("Failed to set the CRAM fields to decode")?;
    // htslib rebuilds the tags whenever it decodes aux tags against a loaded
    // reference, e.g. for `basic` + `rna` in `multi`: an MD string per read
    // that nothing reads.
    let decode_md = if reads_nm_or_md(requirements) { CRAM_DECODE_MD_AUTO } else { 0 };
    reader
        .set_cram_options(hts::hts_fmt_option_CRAM_OPT_DECODE_MD, decode_md)
        .context("Failed to set whether to regenerate CRAM NM/MD tags")
}

/// The htslib `CRAM_OPT_REQUIRED_FIELDS` mask for `requirements`.
///
/// The fields [`RikerRecord`] always populates (name, flags, positions, MAPQ,
/// CIGAR, mate fields, qualities) are always requested. Requesting qualities
/// also keeps the sequence length: with neither `SAM_SEQ` nor `SAM_QUAL`
/// htslib reports a zero-length read, and `sequence_len()` must stay
/// accurate without a sequence decode.
///
/// `SAM_AUX` decodes every aux tag (htslib cannot pick tags), including the
/// read group, which CRAM stores as its own data series. `SAM_SEQ` is added
/// for a declared sequence, and also when `NM` or `MD` is read: htslib
/// regenerates those while decoding the sequence against the reference, which
/// for most slices it loads only when the sequence is required.
fn cram_required_fields(requirements: &RikerRecordRequirements) -> rust_htslib::htslib::sam_fields {
    use rust_htslib::htslib as hts;
    let mut fields = hts::sam_fields_SAM_QNAME
        | hts::sam_fields_SAM_FLAG
        | hts::sam_fields_SAM_RNAME
        | hts::sam_fields_SAM_POS
        | hts::sam_fields_SAM_MAPQ
        | hts::sam_fields_SAM_CIGAR
        | hts::sam_fields_SAM_RNEXT
        | hts::sam_fields_SAM_PNEXT
        | hts::sam_fields_SAM_TLEN
        | hts::sam_fields_SAM_QUAL;
    let aux_tags = requirements.aux_tags();
    if !aux_tags.is_none() {
        fields |= hts::sam_fields_SAM_AUX;
    }
    if requirements.needs_sequence() || reads_nm_or_md(requirements) {
        fields |= hts::sam_fields_SAM_SEQ;
    }
    fields
}

/// Whether a consumer reads `NM` or `MD`. htslib's CRAM writer (samtools)
/// drops both tags by default, and the reader regenerates them from the
/// reference.
fn reads_nm_or_md(requirements: &RikerRecordRequirements) -> bool {
    let aux_tags = requirements.aux_tags();
    aux_tags.includes(*b"NM") || aux_tags.includes(*b"MD")
}

/// Read one CRAM record into `hts_rec` and apply `requirements`. Used
/// by the sequential [`AlignmentReader::fill_record`] path. The
/// generic bound makes it callable on either rust-htslib reader type;
/// the indexed path inlines the same two-line body directly inside
/// [`IndexedAlignmentReader::query_for_each`] because it has its own
/// loop shape.
fn fill_htslib_slot<R: HtsRead>(
    reader: &mut R,
    hts_rec: &mut HtslibRec,
    requirements: &RikerRecordRequirements,
) -> Result<bool> {
    if !hts_rec.read_from(reader)? {
        return Ok(false);
    }
    apply_requirements_htslib(hts_rec, requirements)?;
    Ok(true)
}

/// CRAM analogue of [`apply_requirements`]: decode sequence bytes
/// and/or scan the requested aux tags on a freshly-read [`HtslibRec`].
fn apply_requirements_htslib(
    hts_rec: &mut HtslibRec,
    requirements: &RikerRecordRequirements,
) -> Result<()> {
    if requirements.needs_sequence() {
        hts_rec.decode_sequence();
    }
    match requirements.aux_tags() {
        AuxTagRequirements::None => {}
        AuxTagRequirements::Specific { values, presence } => {
            hts_rec.scan_aux_tags(values, presence)?;
        }
        AuxTagRequirements::All => {
            hts_rec.decode_all_aux()?;
        }
    }
    Ok(())
}

/// Translate a noodles `Region` into a fetch on a rust-htslib indexed
/// reader. Uses the contig name (not a noodles-side tid lookup) so the
/// htslib reader resolves the index against its own header view.
fn fetch_htslib_region(
    reader: &mut rust_htslib::bam::IndexedReader,
    region: &Region,
) -> Result<()> {
    use noodles::core::region::Interval;
    use rust_htslib::bam::FetchDefinition;

    let name = region.name();
    let interval: Interval = region.interval();
    // htslib bounds are 0-based half-open `[start, stop)` as `i64`.
    // The unbounded-interval branches use 0 / i64::MAX (start of contig
    // / end of contig); the `try_from` failure branches `.expect()` the
    // conversion because no real genomic position approaches i64::MAX
    // (it would imply contigs >9.2 billion bases). Marking the
    // unreachable case explicitly is clearer than the prior asymmetric
    // `unwrap_or(0)` / `unwrap_or(i64::MAX)`.
    let start: i64 = interval
        .start()
        .map_or(0, |p| i64::try_from(usize::from(p) - 1).expect("genomic position fits in i64"));
    let stop: i64 = interval
        .end()
        .map_or(i64::MAX, |p| i64::try_from(usize::from(p)).expect("genomic position fits in i64"));

    reader
        .fetch(FetchDefinition::RegionString(name, start, stop))
        .with_context(|| format!("Failed to fetch CRAM region: {region}"))
}

/// Parse htslib's textual SAM header bytes into a `noodles::sam::Header`.
/// Done at open time so the rest of the pipeline can keep using the
/// noodles `Header` type — every collector, helper, and `Region` query
/// downstream is already parameterised on it.
fn parse_htslib_header(header_bytes: &[u8]) -> Result<Header> {
    // htslib's `HeaderView::as_bytes()` returns an empty slice when
    // `sam_hdr_str()` returns NULL — i.e. when the CRAM has no parsed
    // header at all. noodles' `read_header()` would silently produce
    // an empty `Header` from those zero bytes; reject it explicitly so
    // a malformed CRAM surfaces as an error instead of an empty
    // reference list that silently mis-attributes records to whatever
    // contig index 0 happens to mean elsewhere.
    if header_bytes.is_empty() {
        bail!("CRAM file has an empty SAM header — refusing to read records against it");
    }
    let mut parser = sam::io::Reader::new(Cursor::new(header_bytes));
    parser.read_header().context("Failed to parse SAM header text from CRAM")
}

/// Detected alignment file format based on extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlignmentFormat {
    Sam,
    GzippedSam,
    Bam,
    Cram,
}

/// Detect the alignment format from a file extension.
///
/// Recognizes `.sam`, `.sam.gz`, `.bam`, and `.cram` (case-insensitive).
pub(crate) fn detect_format(path: &Path) -> Result<AlignmentFormat> {
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

/// Confirm that an index for the CRAM at `path` exists in either of
/// the two layouts htslib auto-discovers:
///
/// - `<file>.cram.crai` (suffix appended) — `samtools index`'s default
/// - `<file>.crai`      (extension replaced — sibling) — Picard's default
///
/// htslib will happily open either, so rejecting one would surface a
/// confusing "CRAM index not found" error for an index htslib would have
/// used. The pre-check exists only to give a friendlier error than htslib's
/// generic open failure.
fn validate_cram_index_exists(path: &Path) -> Result<()> {
    let crai_sibling = path.with_extension("crai");
    let crai_suffix = append_extension(path, ".crai");
    if crai_suffix.exists() || crai_sibling.exists() {
        return Ok(());
    }
    bail!(
        "CRAM index not found. Expected one of:\n  {}\n  {}\n\
         Run `samtools index {}` to create one.",
        crai_suffix.display(),
        crai_sibling.display(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    /// Helper: create an empty file at `path`.
    fn touch(path: &Path) {
        fs::write(path, b"").expect("create empty file");
    }

    #[test]
    fn select_bam_index_finds_suffix_form() {
        // `samtools index` default: foo.bam → foo.bam.bai
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("sample.bam");
        touch(&bam);
        let idx = dir.path().join("sample.bam.bai");
        touch(&idx);
        assert_eq!(select_bam_index(&bam).expect("should accept suffix-form .bam.bai"), idx);
    }

    #[test]
    fn select_bam_index_finds_sibling_form() {
        // Picard default: foo.bam → foo.bai
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("sample.bam");
        touch(&bam);
        let idx = dir.path().join("sample.bai");
        touch(&idx);
        assert_eq!(select_bam_index(&bam).expect("should accept sibling-form .bai"), idx);
    }

    #[test]
    fn select_bam_index_finds_csi() {
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("sample.bam");
        touch(&bam);
        let idx = dir.path().join("sample.bam.csi");
        touch(&idx);
        assert_eq!(select_bam_index(&bam).expect("should accept .csi"), idx);
    }

    #[test]
    fn select_bam_index_missing_errors() {
        let dir = TempDir::new().unwrap();
        let bam = dir.path().join("sample.bam");
        touch(&bam);
        let err = select_bam_index(&bam).expect_err("missing index should fail");
        let msg = err.to_string();
        assert!(msg.contains("BAM index not found"), "unexpected error: {msg}");
        assert!(msg.contains("sample.bam.bai"), "should mention suffix form: {msg}");
        assert!(msg.contains("sample.bai"), "should mention sibling form: {msg}");
    }

    #[test]
    fn pick_bam_index_prefers_most_recently_modified() {
        use std::time::{Duration, UNIX_EPOCH};
        let bam = Path::new("/data/sample.bam");
        let older = (bam.with_extension("bai"), UNIX_EPOCH + Duration::from_secs(100));
        let newer = (append_extension(bam, ".bai"), UNIX_EPOCH + Duration::from_secs(200));
        // Insertion order must not matter; the newer mtime always wins.
        assert_eq!(pick_bam_index(bam, vec![older.clone(), newer.clone()]).unwrap(), newer.0);
        assert_eq!(pick_bam_index(bam, vec![newer.clone(), older.clone()]).unwrap(), newer.0);
    }

    #[test]
    fn validate_cram_index_finds_crai() {
        let dir = TempDir::new().unwrap();
        let cram = dir.path().join("sample.cram");
        touch(&cram);
        touch(&dir.path().join("sample.cram.crai"));
        validate_cram_index_exists(&cram).expect("should accept .cram.crai");
    }

    #[test]
    fn validate_cram_index_missing_errors() {
        let dir = TempDir::new().unwrap();
        let cram = dir.path().join("sample.cram");
        touch(&cram);
        let err = validate_cram_index_exists(&cram).expect_err("missing .crai should fail");
        assert!(err.to_string().contains("CRAM index not found"));
    }

    /// Whether the CRAM mask for `requirements` requests `(bases, aux tags)`.
    fn cram_mask_bases_and_aux(requirements: &RikerRecordRequirements) -> (bool, bool) {
        use rust_htslib::htslib as hts;
        let fields = cram_required_fields(requirements);
        (fields & hts::sam_fields_SAM_SEQ != 0, fields & hts::sam_fields_SAM_AUX != 0)
    }

    #[test]
    fn cram_mask_without_requirements_requests_only_the_always_populated_fields() {
        use rust_htslib::htslib as hts;
        let always_populated = hts::sam_fields_SAM_QNAME
            | hts::sam_fields_SAM_FLAG
            | hts::sam_fields_SAM_RNAME
            | hts::sam_fields_SAM_POS
            | hts::sam_fields_SAM_MAPQ
            | hts::sam_fields_SAM_CIGAR
            | hts::sam_fields_SAM_RNEXT
            | hts::sam_fields_SAM_PNEXT
            | hts::sam_fields_SAM_TLEN
            | hts::sam_fields_SAM_QUAL;
        assert_eq!(cram_required_fields(&RikerRecordRequirements::NONE), always_populated);
    }

    #[test]
    fn cram_mask_for_the_sequence_requests_the_bases_only() {
        let requirements = RikerRecordRequirements::NONE.with_sequence();
        assert_eq!(cram_mask_bases_and_aux(&requirements), (true, false));
    }

    #[test]
    fn cram_mask_for_a_stored_tag_requests_the_aux_tags_only() {
        let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"MC");
        assert_eq!(cram_mask_bases_and_aux(&requirements), (false, true));
    }

    #[test]
    fn cram_mask_for_nm_requests_the_bases_and_aux_tags() {
        let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        assert_eq!(cram_mask_bases_and_aux(&requirements), (true, true));
    }

    #[test]
    fn cram_mask_for_md_presence_requests_the_bases_and_aux_tags() {
        let requirements = RikerRecordRequirements::NONE.with_aux_tag_presence(*b"MD");
        assert_eq!(cram_mask_bases_and_aux(&requirements), (true, true));
    }

    #[test]
    fn cram_mask_for_all_aux_tags_requests_the_bases_and_aux_tags() {
        let requirements = RikerRecordRequirements::NONE.with_all_aux_tags();
        assert_eq!(cram_mask_bases_and_aux(&requirements), (true, true));
    }

    /// Read a CRAM of ten mapped reads, which htslib wrote without `NM`/`MD`,
    /// through htslib with the options `requirements` set, and count the
    /// records htslib returns with an `NM` and with an `MD` tag. Reads htslib
    /// directly because riker's accessors hide undeclared tags.
    fn nm_and_md_tags_decoded(requirements: &RikerRecordRequirements) -> Result<(usize, usize)> {
        use crate::test_support::{FastaBuilder, SamBuilder, SortOrder, read};
        let reference = FastaBuilder::hg38_mini();
        let fasta = reference.to_temp_fasta()?;
        let mut builder =
            SamBuilder::with_contigs(&reference.dict()).sort_order(SortOrder::Coordinate);
        for i in 0..10 {
            builder.add(read().name(format!("r{i}")).at("chr1", 100 + i * 20).len(50).sub(5, b'A'));
        }
        let cram = builder.to_temp_htslib_cram(fasta.path())?;

        let mut reader = rust_htslib::bam::Reader::from_path(cram.path())?;
        reader.set_reference(fasta.path())?;
        set_cram_required_fields(&mut reader, requirements)?;
        let (mut with_nm, mut with_md) = (0, 0);
        for record in reader.records() {
            let record = record?;
            with_nm += usize::from(record.aux(b"NM").is_ok());
            with_md += usize::from(record.aux(b"MD").is_ok());
        }
        Ok((with_nm, with_md))
    }

    /// Decoding both the bases and the aux tags would let htslib rebuild the
    /// dropped `NM` and `MD` tags; with neither declared it must not.
    #[test]
    fn cram_decode_skips_nm_and_md_regeneration_when_neither_is_read() -> Result<()> {
        let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"MC");
        assert_eq!(nm_and_md_tags_decoded(&requirements)?, (0, 0));
        Ok(())
    }

    #[test]
    fn cram_decode_regenerates_nm_and_md_when_nm_is_read() -> Result<()> {
        let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        assert_eq!(nm_and_md_tags_decoded(&requirements)?, (10, 10));
        Ok(())
    }
}
