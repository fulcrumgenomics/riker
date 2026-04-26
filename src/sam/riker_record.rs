//! A reusable, infallible alignment record type for riker's hot paths.
//!
//! ## Why
//!
//! Profiling `riker wgs` on a 12× BAM showed the reader thread spending ~28%
//! of its CPU (~53 s of 187 s) parsing BAM aux tags that no collector ever
//! reads, plus another ~10% (~18 s) on `Data::clear` / `Data::insert` churn
//! when we reuse a `RecordBuf` across records. The noodles `bam::Record`
//! type is already lazy — it holds the raw BAM bytes and only decodes aux
//! data on demand — but its accessors are fallible, so every call site
//! needs `?` or `.unwrap()` even for fields that are effectively never
//! malformed. That fallibility is the main reason we haven't used it.
//!
//! [`RikerRecord`] solves both problems:
//!
//! 1. It's an **enum** with one variant wrapping `bam::Record` (the BAM
//!    fast path — lazy aux data) and one wrapping `RecordBuf` (the
//!    fallback for SAM and CRAM, where noodles already decodes eagerly).
//! 2. **Scalars are validated once at construction** and cached, so the
//!    downstream accessors return `Flags`, `Option<usize>`, `Option<Position>`
//!    etc. directly — never `Result<_>`. Malformed records surface as an
//!    error at read time, not as propagating failures through the hot path.
//!
//! ## Enum vs. trait dispatch
//!
//! For two variants with small, frequently-called accessors, a `match` is
//! typically faster than a vtable indirection: each arm inlines at its
//! call site and the branch predictor pins the dominant variant in
//! homogeneous inputs (which is what riker sees — one file type per run).
//! Generics / monomorphisation would work too but would force every
//! downstream function to become generic.
//!
//! ## Reusing allocations
//!
//! [`BamRec::read_from`] reads the next BAM record directly into an
//! existing [`BamRec`], reusing the underlying `bam::Record`'s byte
//! buffer. Cached scalars are refreshed in the same call. The
//! [`RikerRecord::Fallback`] variant reuses a `RecordBuf` the same way.
//!
//! ## Expensive fields are opt-in
//!
//! The BAM fast path supports three decoder-driven fields — sequence
//! bases, aux tags (targeted), and aux tags (full). Consumers opt into
//! each via their [`RikerRecordRequirements`]; the reader
//! ([`crate::sam::alignment_reader::AlignmentReader::fill_record`] and
//! friends) consults the requirements and calls the matching fillers on
//! [`BamRec`] ([`BamRec::decode_sequence`],
//! [`BamRec::scan_aux_tags`], [`BamRec::decode_all_aux`]).

use std::collections::BTreeSet;

use anyhow::{Context, Result, anyhow};
use bstr::BStr;
use noodles::bam;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use smallvec::SmallVec;

/// BAM aux tag — two ASCII bytes, e.g. `*b"NM"` for the edit-distance tag.
pub type TagKey = [u8; 2];

// ─── RikerRecord ─────────────────────────────────────────────────────────────

/// Alignment record that can come from either a BAM fast path
/// ([`BamRec`], with lazy aux data) or a generic fallback ([`FallbackRec`])
/// used for SAM and CRAM. See the module-level docs for why this is an
/// enum and not a trait.
///
/// `BamRec` is intentionally larger than `FallbackRec` (cigar / quality
/// mirror buffers, aux store, etc) — boxing it would force a per-record
/// heap allocation and defeat the in-place fill loop the entire design
/// is built around. The size delta is a deliberate trade.
#[allow(clippy::large_enum_variant, reason = "see docstring above")]
#[derive(Clone)]
pub enum RikerRecord {
    /// BAM records: the underlying `bam::Record` holds raw bytes and
    /// defers aux-tag decoding until a filler is invoked.
    Bam(BamRec),
    /// SAM/CRAM records: a fully-decoded `RecordBuf`. No decode savings
    /// here, but the type-level split lets `Bam` stay fast in the common
    /// case.
    Fallback(FallbackRec),
}

impl RikerRecord {
    /// Build a [`RikerRecord::Fallback`] from any `sam::alignment::Record`
    /// (used for CRAM, where noodles reads records as trait objects).
    ///
    /// # Errors
    /// Returns an error if the source record can't be converted to a
    /// `RecordBuf`.
    pub fn from_alignment_record<R: sam::alignment::Record>(
        header: &sam::Header,
        record: &R,
    ) -> Result<Self> {
        let buf = RecordBuf::try_from_alignment_record(header, record)
            .context("Failed to convert alignment record to RecordBuf")?;
        Ok(Self::Fallback(FallbackRec::from_record_buf(buf)))
    }

    // ── Scalar-field accessors ──────────────────────────────────────────────

    /// Flags bitmap.
    #[inline]
    #[must_use]
    pub fn flags(&self) -> Flags {
        match self {
            Self::Bam(r) => r.flags,
            Self::Fallback(r) => r.inner.flags(),
        }
    }

    /// 0-based reference sequence index, or `None` for unmapped reads.
    #[inline]
    #[must_use]
    pub fn reference_sequence_id(&self) -> Option<usize> {
        match self {
            Self::Bam(r) => r.ref_id,
            Self::Fallback(r) => r.inner.reference_sequence_id(),
        }
    }

    /// 1-based alignment start position, or `None` for unmapped reads.
    #[inline]
    #[must_use]
    pub fn alignment_start(&self) -> Option<Position> {
        match self {
            Self::Bam(r) => r.pos,
            Self::Fallback(r) => r.inner.alignment_start(),
        }
    }

    /// 1-based inclusive alignment end, derived from `alignment_start` +
    /// the reference span of the CIGAR. `None` when unmapped or when the
    /// CIGAR has no reference-consuming ops.
    ///
    /// Cached on both variants: `RecordBuf::alignment_end()` walks the
    /// CIGAR per call, which is too expensive for the mate-buffer hot path.
    #[inline]
    #[must_use]
    pub fn alignment_end(&self) -> Option<Position> {
        match self {
            Self::Bam(r) => r.alignment_end,
            Self::Fallback(r) => r.alignment_end,
        }
    }

    /// Mapping quality, or `None` when missing (BAM stores 255).
    #[inline]
    #[must_use]
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        match self {
            Self::Bam(r) => r.mapq,
            Self::Fallback(r) => r.inner.mapping_quality(),
        }
    }

    /// Mate's 0-based reference index, or `None` when the mate is
    /// unmapped or the read is unpaired.
    #[inline]
    #[must_use]
    pub fn mate_reference_sequence_id(&self) -> Option<usize> {
        match self {
            Self::Bam(r) => r.mate_ref_id,
            Self::Fallback(r) => r.inner.mate_reference_sequence_id(),
        }
    }

    /// Mate's 1-based alignment start, or `None`.
    #[inline]
    #[must_use]
    pub fn mate_alignment_start(&self) -> Option<Position> {
        match self {
            Self::Bam(r) => r.mate_pos,
            Self::Fallback(r) => r.inner.mate_alignment_start(),
        }
    }

    /// Signed template (insert) length.
    #[inline]
    #[must_use]
    pub fn template_length(&self) -> i32 {
        match self {
            Self::Bam(r) => r.tlen,
            Self::Fallback(r) => r.inner.template_length(),
        }
    }

    // ── Byte-slice accessors ────────────────────────────────────────────────

    /// Read name, without the trailing NUL. `None` when the record has no
    /// name (BAM stores the sentinel `*`).
    ///
    /// Returned as `&BStr` for readable `Display` output — `BStr` prints
    /// as a (lossy) UTF-8 string rather than as a byte-array debug repr,
    /// which is what you want for log messages and panics.
    #[inline]
    #[must_use]
    pub fn name(&self) -> Option<&BStr> {
        match self {
            Self::Bam(r) => r.inner.name(),
            Self::Fallback(r) => r.inner.name(),
        }
    }

    /// Raw Phred base-quality bytes (0-based, *not* ASCII+33).
    #[inline]
    #[must_use]
    pub fn quality_scores(&self) -> &[u8] {
        match self {
            Self::Bam(r) => &r.quality_bytes,
            Self::Fallback(r) => r.inner.quality_scores().as_ref(),
        }
    }

    // ── CIGAR ───────────────────────────────────────────────────────────────

    /// Iterator over the CIGAR operations. For BAM we walk the raw packed
    /// bytes (no allocation); for the fallback we iterate the pre-decoded
    /// `Vec<Op>` directly.
    #[must_use]
    pub fn cigar_ops(&self) -> CigarOps<'_> {
        match self {
            Self::Bam(r) => CigarOps::Bam(r.cigar_bytes.chunks_exact(4)),
            Self::Fallback(r) => CigarOps::Fallback(r.inner.cigar().as_ref().iter()),
        }
    }

    /// Number of CIGAR ops (mainly for short-circuit checks).
    #[must_use]
    pub fn cigar_len(&self) -> usize {
        match self {
            Self::Bam(r) => r.cigar_bytes.len() / 4,
            Self::Fallback(r) => r.inner.cigar().as_ref().len(),
        }
    }

    // ── Sequence ────────────────────────────────────────────────────────────

    /// Number of bases in the sequence. Always available — doesn't
    /// require the sequence to have been decoded.
    #[must_use]
    pub fn sequence_len(&self) -> usize {
        match self {
            Self::Bam(r) => r.inner.sequence().len(),
            Self::Fallback(r) => r.inner.sequence().len(),
        }
    }

    /// Decoded ASCII sequence bases, or an empty slice if the sequence is
    /// absent. "Absent" conflates two cases — matching noodles' own
    /// convention for `RecordBuf::sequence()`:
    ///
    /// 1. The consumer's [`RikerRecordRequirements`] didn't include
    ///    `sequence`, so the reader helper never called
    ///    [`BamRec::decode_sequence`]. (BAM fast path only.)
    /// 2. The file genuinely had no sequence bases (`l_seq == 0`).
    ///
    /// Collectors that need to distinguish the two should either ensure
    /// they've declared a need for sequence (so the empty case only ever
    /// means "file had none") or check [`Self::sequence_len`] directly —
    /// it reports the file-level `l_seq` regardless of whether the
    /// decoder has run.
    ///
    /// [`RikerRecordRequirements`]: crate::sam::riker_record::RikerRecordRequirements
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        match self {
            Self::Bam(r) => {
                // Empty slice covers two cases per the doc comment above:
                // (1) the consumer didn't declare `with_sequence()` and
                // (2) the file genuinely had no bases (`l_seq == 0`).
                // The `if` is also load-bearing because `sequence_buf` is
                // a reused allocation: without it, a record where the
                // decode didn't run would expose the *previous* record's
                // bases.
                if r.sequence_populated { r.sequence_buf.as_slice() } else { &[] }
            }
            Self::Fallback(r) => r.inner.sequence().as_ref(),
        }
    }

    // ── Aux tags ────────────────────────────────────────────────────────────

    /// Look up an aux tag.
    ///
    /// Returns an *owned* [`AuxValue`] because the two paths can't share
    /// a reference cleanly — BAM stores decoded values in the record's
    /// aux store, while `Fallback` converts from noodles' pre-decoded
    /// `RecordBuf::Data` on the fly. For scalar tags (`Int`, `Float`,
    /// `Character`) the return-by-value is free; only `String`/`Hex`
    /// clone bytes, and no hot-path collector reads string aux tags
    /// today.
    ///
    /// For the BAM variant, returns `Some(...)` only when the reader has
    /// run a filler that covered this tag (either
    /// [`BamRec::scan_aux_tags`] for targeted
    /// tags, or [`BamRec::decode_all_aux`] for any tag). Returns `None`
    /// both when the tag wasn't requested *and* when it was requested
    /// but not present on the record — that ambiguity is intentional.
    ///
    /// For the `Fallback` variant, we consult the decoded `Data` on the
    /// underlying `RecordBuf` and convert on the fly.
    #[must_use]
    pub fn aux_tag(&self, tag: TagKey) -> Option<AuxValue> {
        match self {
            Self::Bam(r) => {
                if r.aux_populated {
                    r.aux_store.get(tag).cloned()
                } else {
                    None
                }
            }
            Self::Fallback(r) => {
                let tag_obj = sam::alignment::record::data::field::Tag::from(tag);
                let value = r.inner.data().get(&tag_obj)?;
                record_buf_value_to_aux(value)
            }
        }
    }

    /// Cheaper companion to [`Self::aux_tag`] when only existence
    /// matters. On the BAM fast path returns `true` if the scanner
    /// observed the tag — either because the consumer requested its
    /// value (`with_aux_tag`) or its presence
    /// (`with_aux_tag_presence`). On `Fallback` consults the decoded
    /// `Data` directly. Returns `false` for unrequested tags on BAM
    /// (matching the [`aux_tag`] behaviour).
    ///
    /// [`aux_tag`]: Self::aux_tag
    #[must_use]
    pub fn has_aux_tag(&self, tag: TagKey) -> bool {
        match self {
            Self::Bam(r) => {
                r.aux_populated && (r.aux_store.get(tag).is_some() || r.aux_present.contains(&tag))
            }
            Self::Fallback(r) => {
                let tag_obj = sam::alignment::record::data::field::Tag::from(tag);
                r.inner.data().get(&tag_obj).is_some()
            }
        }
    }
}

// ─── BamRec ─────────────────────────────────────────────────────────────────

/// BAM-backed record. Wraps [`noodles::bam::Record`] and caches the fixed
/// scalar fields so that our accessors can be infallible and plain.
///
/// The cigar and quality bytes are mirrored into reusable `Vec`s owned by
/// this struct. noodles' `bam::record::Cigar` / `QualityScores` view
/// types go through `AsRef<[u8]>` which erases the inner `'a` lifetime,
/// so we can't return a slice into the BAM buffer with a lifetime tied
/// to `&BamRec`. The per-record memcpy is bounded — a few hundred bytes
/// per record, reused allocations — and still orders of magnitude
/// cheaper than the aux-tag parse we're avoiding.
///
/// The read *name* is NOT mirrored: `bam::Record::name()` returns
/// `Option<&BStr>` whose lifetime does flow through correctly via
/// dereference, so we pass it through zero-copy (and `BStr` has a
/// human-readable `Display` impl, which `&[u8]` doesn't).
#[derive(Clone)]
pub struct BamRec {
    /// Underlying noodles record. Owns the raw BAM bytes.
    inner: bam::Record,

    // ── Cached scalar fields ────────────────────────────────────────────────
    // Set in `refresh_cache` after each read; kept in sync with `inner`.
    // The cache lets us hand callers plain values rather than the
    // `Option<io::Result<…>>` shapes noodles exposes on `bam::Record`.
    flags: Flags,
    ref_id: Option<usize>,
    pos: Option<Position>,
    mapq: Option<MappingQuality>,
    mate_ref_id: Option<usize>,
    mate_pos: Option<Position>,
    tlen: i32,
    /// 1-based inclusive alignment end, computed from `pos` + cigar
    /// reference span. `None` when `pos` is `None` or the span can't be
    /// computed.
    alignment_end: Option<Position>,

    // ── Reusable byte buffers ───────────────────────────────────────────────
    /// Raw CIGAR bytes, 4 bytes per op.
    cigar_bytes: Vec<u8>,
    /// Raw Phred base-quality bytes (0-based, not ASCII+33).
    quality_bytes: Vec<u8>,

    // ── Optional, filler-driven fields ──────────────────────────────────────
    // These are populated on demand by `decode_sequence`,
    // `scan_aux_tags`, and `decode_all_aux`. Their `*_populated` flags
    // are reset to false
    // inside `refresh_cache` so accessors can distinguish "not populated
    // this record" from a stale hit. `BamRec` itself doesn't know about
    // `RikerRecordRequirements` — a helper on `AlignmentReader` reads the
    // requirements and calls the matching fillers after each read.
    /// Decoded ASCII sequence bases. Backing Vec is kept even when
    /// `sequence_populated` is false so allocation is reused.
    sequence_buf: Vec<u8>,
    sequence_populated: bool,
    /// Store for decoded aux tag values. Cleared on each read; callers
    /// invoke `scan_aux_tags` (targeted) or `decode_all_aux` (full) to fill.
    aux_store: AuxTagStore,
    aux_populated: bool,
    /// Tags the scanner observed for presence-only requests, *not* in
    /// `aux_store`. Lookups via `has_aux_tag` consult both lists.
    /// `SmallVec` chosen because typical cardinality is 0-2 and matches
    /// `aux_store`'s reuse semantics.
    aux_present: SmallVec<[TagKey; 4]>,
}

impl BamRec {
    /// Build an empty, default-initialised record. Intended to be
    /// populated via [`Self::read_from`] in a reuse loop.
    #[must_use]
    pub(crate) fn new() -> Self {
        Self {
            inner: bam::Record::default(),
            flags: Flags::empty(),
            ref_id: None,
            pos: None,
            mapq: None,
            mate_ref_id: None,
            mate_pos: None,
            tlen: 0,
            alignment_end: None,
            cigar_bytes: Vec::new(),
            quality_bytes: Vec::new(),
            sequence_buf: Vec::new(),
            sequence_populated: false,
            aux_store: AuxTagStore::new(),
            aux_populated: false,
            aux_present: SmallVec::new(),
        }
    }

    /// Read the next record directly from a noodles BAM reader into
    /// `self`, reusing all of this struct's buffer allocations. Returns
    /// the BAM block size, or `0` at EOF.
    ///
    /// # Errors
    /// Returns an error if the underlying read fails or the record's
    /// scalars can't be validated.
    pub(crate) fn read_from<R: std::io::Read>(
        &mut self,
        reader: &mut bam::io::Reader<R>,
    ) -> Result<usize> {
        let n = reader.read_record(&mut self.inner).context("Failed to read BAM record")?;
        if n == 0 {
            return Ok(0);
        }
        self.refresh_cache()?;
        Ok(n)
    }

    /// Replace the inner `bam::Record` with the supplied one and refresh
    /// the cached scalars. Useful when the record is produced by a higher-
    /// level iterator (e.g. `bam::io::indexed_reader::Reader::query`) that
    /// already yields owned `bam::Record` values.
    ///
    /// # Errors
    /// Returns an error if the record's scalars can't be validated.
    pub(crate) fn install(&mut self, record: bam::Record) -> Result<()> {
        self.inner = record;
        self.refresh_cache()
    }

    /// Decode the sequence bases into `self.sequence_buf` using the SIMD
    /// nibble decoder.
    pub(crate) fn decode_sequence(&mut self) {
        let seq_view = self.inner.sequence();
        let packed: &[u8] = seq_view.as_ref();
        super::simd_seq::decode_packed_sequence_into(
            packed,
            seq_view.len(),
            &mut self.sequence_buf,
        );
        self.sequence_populated = true;
    }

    /// Scan the aux byte block once. Tags in `values` are decoded into
    /// `self.aux_store`; tags in `presence` are recorded in
    /// `self.aux_present` without their payloads being parsed (skipping
    /// the per-record `String`/`Hex` allocation when only existence
    /// matters).
    ///
    /// # Errors
    /// Returns an error if the aux byte block is malformed.
    pub(crate) fn scan_aux_tags(
        &mut self,
        values: &BTreeSet<TagKey>,
        presence: &BTreeSet<TagKey>,
    ) -> Result<()> {
        let aux = self.inner.data();
        let aux_bytes: &[u8] = aux.as_ref();
        scan_specific(aux_bytes, values, presence, &mut self.aux_store, &mut self.aux_present)
            .context("Failed to scan aux tags")?;
        self.aux_populated = true;
        Ok(())
    }

    /// Decode *every* aux tag in the record into the aux store.
    /// Expensive — we walk the entire aux block and materialise each
    /// field. Prefer [`Self::scan_aux_tags`] when the consumer knows
    /// the tags it wants.
    ///
    /// # Errors
    /// Returns an error if any aux field can't be decoded.
    pub(crate) fn decode_all_aux(&mut self) -> Result<()> {
        let aux = self.inner.data();
        let aux_bytes: &[u8] = aux.as_ref();
        scan_all(aux_bytes, &mut self.aux_store).context("Failed to decode all aux tags")?;
        self.aux_populated = true;
        Ok(())
    }

    /// Re-read the cached scalars from the current `inner` bytes. Called
    /// by [`Self::read_from`] after every successful underlying read so
    /// downstream accessors see up-to-date values.
    fn refresh_cache(&mut self) -> Result<()> {
        self.flags = self.inner.flags();
        self.ref_id = match self.inner.reference_sequence_id() {
            None => None,
            Some(res) => Some(res.context("invalid reference sequence id")?),
        };
        self.pos = match self.inner.alignment_start() {
            None => None,
            Some(res) => Some(res.context("invalid alignment start")?),
        };
        self.mapq = self.inner.mapping_quality();
        self.mate_ref_id = match self.inner.mate_reference_sequence_id() {
            None => None,
            Some(res) => Some(res.context("invalid mate reference sequence id")?),
        };
        self.mate_pos = match self.inner.mate_alignment_start() {
            None => None,
            Some(res) => Some(res.context("invalid mate alignment start")?),
        };
        self.tlen = self.inner.template_length();

        // Mirror the variable-length byte fields into our own owned
        // buffers so accessors can hand out slices with a lifetime tied
        // to `&BamRec`. `name` is NOT mirrored — `bam::Record::name()`
        // returns `Option<&BStr>` whose lifetime does flow through via
        // dereference.
        self.cigar_bytes.clear();
        self.cigar_bytes.extend_from_slice(self.inner.cigar().as_ref());

        self.quality_bytes.clear();
        self.quality_bytes.extend_from_slice(self.inner.quality_scores().as_ref());

        // Validate the packed CIGAR codes once at read time; downstream
        // consumers (the `CigarOps` iterator) then call `code_to_kind`
        // infallibly. This is what stops a corrupted BAM from silently
        // producing wrong span/coverage numbers.
        validate_cigar_codes(&self.cigar_bytes)?;

        self.alignment_end = compute_alignment_end(self.pos, &self.cigar_bytes)?;

        // Invalidate decoder-driven caches so accessors see "not
        // populated for this record" until a filler is called.
        self.sequence_populated = false;
        self.aux_populated = false;
        self.aux_store.clear();
        self.aux_present.clear();

        Ok(())
    }
}

// No public `Default` for `BamRec`: the constructors are pub(crate) and
// an externally-defaulted record can never be filled or refreshed.

// ─── FallbackRec ────────────────────────────────────────────────────────────

/// SAM/CRAM-backed record. Wraps a fully-decoded `RecordBuf` — the
/// decode cost is already paid by noodles for these formats, so we
/// don't pretend there's a fast path here.
///
/// `RecordBuf`'s own accessors are all infallible (`Option<usize>`,
/// `Option<Position>`, etc.), so the enum arm just delegates straight
/// through. The only field we cache is `alignment_end`, since
/// `RecordBuf::alignment_end()` walks the CIGAR on each call and the
/// mate buffer hits that accessor per record.
#[derive(Clone)]
pub struct FallbackRec {
    inner: RecordBuf,
    alignment_end: Option<Position>,
}

impl FallbackRec {
    /// Build a `FallbackRec` from a freshly-read `RecordBuf`.
    #[must_use]
    pub(crate) fn from_record_buf(inner: RecordBuf) -> Self {
        let alignment_end = inner.alignment_end();
        Self { inner, alignment_end }
    }

    /// Build an empty `FallbackRec`. Intended to be filled in-place via
    /// [`inner_mut`] + [`refresh_cache`] in a reuse loop.
    ///
    /// [`inner_mut`]: Self::inner_mut
    /// [`refresh_cache`]: Self::refresh_cache
    #[must_use]
    pub(crate) fn empty() -> Self {
        Self { inner: RecordBuf::default(), alignment_end: None }
    }

    /// Mutable access to the inner `RecordBuf` — used by readers that
    /// want to fill the record in place (e.g. the `multi` pool). Call
    /// [`refresh_cache`] after mutating to bring the cached scalars back
    /// in sync. Restricted to the crate so external callers can't bypass
    /// `refresh_cache`.
    ///
    /// [`refresh_cache`]: Self::refresh_cache
    pub(crate) fn inner_mut(&mut self) -> &mut RecordBuf {
        &mut self.inner
    }

    /// Re-read the cached scalars (currently just `alignment_end`) from
    /// the inner `RecordBuf`. Cheap — a single CIGAR walk. Call after
    /// filling the record in-place via [`inner_mut`]. Crate-internal so
    /// external callers can't observe the half-stale state between an
    /// `inner_mut` mutation and the matching refresh.
    ///
    /// [`inner_mut`]: Self::inner_mut
    pub(crate) fn refresh_cache(&mut self) {
        self.alignment_end = self.inner.alignment_end();
    }
}

// No public `Default` for `FallbackRec`: same reasoning as `BamRec`
// above — `empty()` is `pub(crate)` so external callers can't fill it.

// ─── CigarOps iterator ──────────────────────────────────────────────────────

/// Iterator over the CIGAR ops of a [`RikerRecord`]. For `Bam` we walk
/// the raw packed bytes; for `Fallback` we re-iterate the pre-decoded
/// ops from the `RecordBuf`.
pub enum CigarOps<'a> {
    /// Variant wrapping chunks of 4 bytes each for a BAM-packed CIGAR.
    Bam(std::slice::ChunksExact<'a, u8>),
    /// Variant iterating pre-decoded ops from a `RecordBuf`.
    Fallback(std::slice::Iter<'a, Op>),
}

impl Iterator for CigarOps<'_> {
    type Item = Op;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Bam(chunks) => chunks.next().map(|c| {
                // BAM-packed: bits [31..4] = length, bits [3..0] = op
                // kind.
                let word = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                let kind = code_to_kind(word & 0xF);
                #[allow(clippy::cast_possible_truncation)]
                let len = (word >> 4) as usize;
                Op::new(kind, len)
            }),
            Self::Fallback(iter) => iter.next().copied(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Bam(chunks) => chunks.size_hint(),
            Self::Fallback(iter) => iter.size_hint(),
        }
    }
}

// ─── AuxValue ───────────────────────────────────────────────────────────────

/// Decoded value of a BAM auxiliary tag.
///
/// Integer widths from the BAM wire format (`Int8`/`UInt8`/`Int16`/
/// `UInt16`/`Int32`/`UInt32`) are unified into a single `Int(i64)`
/// variant: the BAM spec doesn't mandate a width per tag and no caller
/// in riker cares which encoding was used. The `Int(i64)` range covers
/// every wire type without loss. Round-tripping a record's exact
/// integer width back to BAM is not supported by this design.
///
/// Array-typed aux values (BAM `B:<type>`) are not yet modelled — any
/// requested `B:`-typed tag will appear as "not present". Add an
/// `Array` variant when a consumer needs one.
#[derive(Debug, Clone, PartialEq)]
pub enum AuxValue {
    /// `A:` single printable character (ASCII).
    Character(u8),
    /// Any of BAM's six integer encodings, sign-extended into `i64`.
    Int(i64),
    /// `f:` single-precision float.
    Float(f32),
    /// `Z:` string value. Stored owned because our aux store outlives
    /// any single call into the record.
    String(Vec<u8>),
    /// `H:` hex-encoded byte string. Same owned-storage rationale.
    Hex(Vec<u8>),
}

impl AuxValue {
    /// Returns the integer value if this tag is `Int`.
    #[inline]
    #[must_use]
    pub fn as_int(&self) -> Option<i64> {
        if let Self::Int(n) = self { Some(*n) } else { None }
    }

    /// Returns the float if this tag is `Float`.
    #[inline]
    #[must_use]
    pub fn as_float(&self) -> Option<f32> {
        if let Self::Float(f) = self { Some(*f) } else { None }
    }

    /// Returns the character (ASCII byte) if this tag is `Character`.
    #[inline]
    #[must_use]
    pub fn as_char(&self) -> Option<u8> {
        if let Self::Character(c) = self { Some(*c) } else { None }
    }

    /// Returns a `BStr` view over a `String`-typed value. Callers can
    /// then `Display` it as text without needing a byte-to-str conversion.
    #[inline]
    #[must_use]
    pub fn as_str(&self) -> Option<&BStr> {
        if let Self::String(bytes) = self { Some(BStr::new(bytes.as_slice())) } else { None }
    }

    /// Returns a `BStr` view over an `H:`-typed value.
    #[inline]
    #[must_use]
    pub fn as_hex(&self) -> Option<&BStr> {
        if let Self::Hex(bytes) = self { Some(BStr::new(bytes.as_slice())) } else { None }
    }
}

// ─── AuxTagStore ────────────────────────────────────────────────────────────

/// Small inline-backed store of decoded aux tags. Linear-scan lookup is
/// used because typical cardinality is 1-3 and fits in the inline buffer
/// without touching the heap.
#[derive(Debug, Default, Clone)]
pub(crate) struct AuxTagStore {
    tags: SmallVec<[(TagKey, AuxValue); 4]>,
}

impl AuxTagStore {
    /// Create an empty store.
    #[must_use]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Drop all stored tags. Keeps the backing capacity for reuse.
    pub(crate) fn clear(&mut self) {
        self.tags.clear();
    }

    /// Insert or replace a tag value.
    pub(crate) fn insert(&mut self, tag: TagKey, value: AuxValue) {
        if let Some((_, existing)) = self.tags.iter_mut().find(|(k, _)| *k == tag) {
            *existing = value;
        } else {
            self.tags.push((tag, value));
        }
    }

    /// Look up a tag.
    #[must_use]
    pub(crate) fn get(&self, tag: TagKey) -> Option<&AuxValue> {
        self.tags.iter().find(|(k, _)| *k == tag).map(|(_, v)| v)
    }

    /// Number of stored tags.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.tags.len()
    }
}

// ─── RikerRecordRequirements ────────────────────────────────────────────────
//
// Declarative "what does this tool need to read from each record" types.
// Most fields (flags, positions, mapq, mate info, name, cigar, qualities)
// are populated unconditionally because they're cheap. The two expensive
// fields are the sequence (4-bit packed → ASCII decode) and aux tags
// (`Data` parse). Each consumer declares which expensive fields it needs;
// the reader unions across active consumers and runs each decoder once
// per record.

/// What a tool needs to read from each record, beyond the always-populated
/// scalars and byte slices. Passed to the reader so the per-record decode
/// work is sized to the actual consumer.
///
/// # Composition
///
/// Use [`RikerRecordRequirements::NONE`] as a starting point and add fields
/// with the builder methods, or compose independent needs with
/// [`RikerRecordRequirements::union`]:
///
/// ```ignore
/// let needs = RikerRecordRequirements::NONE
///     .with_sequence()
///     .with_aux_tag(*b"NM");
///
/// let combined = wgs_needs.union(error_needs).union(alignment_needs);
/// ```
///
/// # Defaults
///
/// [`Default`] returns [`Self::NONE`] — the cheapest possible setting. Tools
/// opt in; they never accidentally drag in expensive work.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct RikerRecordRequirements {
    /// True iff the consumer will call [`RikerRecord::sequence`] or similar
    /// accessors that require the 4-bit-packed → ASCII decode.
    pub(crate) sequence: bool,

    /// Aux-tag access pattern. See [`AuxTagRequirements`].
    pub(crate) aux_tags: AuxTagRequirements,
}

impl RikerRecordRequirements {
    /// Nothing beyond the always-populated fields. Good for tools like
    /// `wgs` and `isize` that only touch flags, positions, mapq, mate info,
    /// cigar, and quality scores.
    pub const NONE: Self = Self { sequence: false, aux_tags: AuxTagRequirements::None };

    /// Whether the consumer needs decoded sequence bases.
    #[must_use]
    pub fn needs_sequence(&self) -> bool {
        self.sequence
    }

    /// Read-only view of the aux-tag access pattern.
    #[must_use]
    pub fn aux_tags(&self) -> &AuxTagRequirements {
        &self.aux_tags
    }

    /// All expensive fields enabled. Rarely the right choice — prefer the
    /// targeted builders. Provided for cases where a tool genuinely needs
    /// every aux tag and the sequence.
    #[must_use]
    pub fn everything() -> Self {
        Self { sequence: true, aux_tags: AuxTagRequirements::All }
    }

    /// Opt in to sequence decode.
    #[must_use]
    pub fn with_sequence(mut self) -> Self {
        self.sequence = true;
        self
    }

    /// Opt in to a specific aux tag whose **decoded value** the consumer
    /// will read via [`RikerRecord::aux_tag`]. Multiple calls compose.
    #[must_use]
    pub fn with_aux_tag(mut self, tag: TagKey) -> Self {
        self.aux_tags = self.aux_tags.with(tag);
        self
    }

    /// Opt in to a specific aux tag's **presence** only — the consumer
    /// will call [`RikerRecord::has_aux_tag`] but never read the decoded
    /// value. Cheaper than [`Self::with_aux_tag`] for `Z`/`H`-typed tags
    /// (e.g. `SA`) because the scanner skips parsing the payload. If the
    /// same tag is later requested via [`Self::with_aux_tag`] the value
    /// path subsumes the presence-only request.
    #[must_use]
    pub fn with_aux_tag_presence(mut self, tag: TagKey) -> Self {
        self.aux_tags = self.aux_tags.with_presence(tag);
        self
    }

    /// Opt in to the full aux-tag `Data` parse. Callers expecting to
    /// iterate many tags (or perform lookups by unknown-at-compile-time
    /// tag) should use this; otherwise prefer [`Self::with_aux_tag`] per
    /// tag.
    #[must_use]
    pub fn with_all_aux_tags(mut self) -> Self {
        self.aux_tags = AuxTagRequirements::All;
        self
    }

    /// Union two sets of needs. Used by `multi` to aggregate across the
    /// active collectors. Deterministic and idempotent.
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        Self {
            sequence: self.sequence || other.sequence,
            aux_tags: self.aux_tags.union(other.aux_tags),
        }
    }
}

impl Default for RikerRecordRequirements {
    fn default() -> Self {
        Self::NONE
    }
}

/// How aux tags are consumed. Three levels of granularity, each enabling
/// progressively more work on the reader thread.
///
/// `Specific` is the common sweet spot: collectors that care about one or
/// two tags (e.g. `error` reads `NM`) pay a targeted scan — no HashMap,
/// no hashing, just a linear walk of the aux byte block looking for those
/// tag prefixes. `All` is strictly more expensive and should be chosen
/// only when the consumer needs iteration over all tags.
///
/// Within `Specific`, tags split into two pools:
///
/// - **`values`** — collectors that will call [`RikerRecord::aux_tag`] on
///   this tag. The scanner parses the wire bytes into an [`AuxValue`]
///   (allocating `String`/`Hex` payloads as needed) and stores it.
/// - **`presence`** — collectors that only need to ask "is this tag
///   present?" via [`RikerRecord::has_aux_tag`]. The scanner records
///   existence without parsing the payload, skipping the per-record
///   `String(Vec<u8>)` allocation that would otherwise be wasted.
#[derive(Debug, Clone, Default, Eq, PartialEq)]
pub enum AuxTagRequirements {
    /// No aux access. Reader doesn't touch the aux block.
    #[default]
    None,

    /// Reader provides targeted access to these tags. The two sub-sets
    /// are disjoint after [`union`](Self::union) — `values` wins when a
    /// tag appears in both, since decoding the value also answers
    /// presence-only queries.
    Specific {
        /// Tags whose decoded value the consumer will read.
        values: BTreeSet<TagKey>,
        /// Tags the consumer will only check for existence.
        presence: BTreeSet<TagKey>,
    },

    /// Full aux `Data` decode. Choose only when the consumer needs to
    /// iterate tags or look up by tags unknown at construction time.
    All,
}

impl AuxTagRequirements {
    /// Add a value-decoded tag. `None` becomes `Specific{values:{tag}}`;
    /// `Specific` accumulates into `values` (and removes from `presence`
    /// if it was there); `All` stays `All`.
    #[must_use]
    pub fn with(self, tag: TagKey) -> Self {
        match self {
            Self::None => {
                Self::Specific { values: BTreeSet::from([tag]), presence: BTreeSet::new() }
            }
            Self::Specific { mut values, mut presence } => {
                values.insert(tag);
                presence.remove(&tag);
                Self::Specific { values, presence }
            }
            Self::All => Self::All,
        }
    }

    /// Add a presence-only tag. `None` becomes
    /// `Specific{presence:{tag}}`; on `Specific`, no-op if the tag is
    /// already in `values` (decode subsumes presence); `All` stays `All`.
    #[must_use]
    pub fn with_presence(self, tag: TagKey) -> Self {
        match self {
            Self::None => {
                Self::Specific { values: BTreeSet::new(), presence: BTreeSet::from([tag]) }
            }
            Self::Specific { values, mut presence } => {
                if !values.contains(&tag) {
                    presence.insert(tag);
                }
                Self::Specific { values, presence }
            }
            Self::All => Self::All,
        }
    }

    /// Union. `All` absorbs anything, `None` is the identity, two
    /// `Specific` operands merge — the result keeps `values` ∪ on the
    /// value side and `(presence_a ∪ presence_b) \ (values_a ∪ values_b)`
    /// on the presence side.
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        match (self, other) {
            (Self::All, _) | (_, Self::All) => Self::All,
            (Self::None, x) | (x, Self::None) => x,
            (
                Self::Specific { values: mut va, presence: mut pa },
                Self::Specific { values: vb, presence: pb },
            ) => {
                va.extend(vb);
                pa.extend(pb);
                pa.retain(|t| !va.contains(t));
                Self::Specific { values: va, presence: pa }
            }
        }
    }

    /// Whether the reader needs to touch the aux block at all.
    #[must_use]
    pub fn is_none(&self) -> bool {
        matches!(self, Self::None)
    }
}

// ─── Module-level helpers ───────────────────────────────────────────────────
//
// These operate on primitive types or external noodles types and don't
// belong as methods on our own structs.

// ── CIGAR ──

/// Validate that every op code in a packed BAM CIGAR is in the spec's
/// 0..=8 range. Run once at read time so downstream iterators can call
/// [`code_to_kind`] infallibly. Codes 9..=15 are reserved by the BAM
/// spec and treating them as anything (including a silent fallback to
/// `SequenceMismatch`) corrupts coverage / span / mismatch metrics.
fn validate_cigar_codes(cigar_bytes: &[u8]) -> Result<()> {
    let (chunks, []) = cigar_bytes.as_chunks::<4>() else {
        return Err(anyhow!("cigar buffer length {} is not a multiple of 4", cigar_bytes.len()));
    };
    for chunk in chunks {
        let code = u32::from_le_bytes(*chunk) & 0xF;
        if code > 8 {
            return Err(anyhow!("invalid BAM CIGAR op code: {code} (valid range is 0..=8)"));
        }
    }
    Ok(())
}

/// Reference-consuming span of a packed BAM CIGAR. Assumes the bytes
/// have already been validated by [`validate_cigar_codes`] — call that
/// first.
fn cigar_reference_span(cigar_bytes: &[u8]) -> Result<u32> {
    let (chunks, []) = cigar_bytes.as_chunks::<4>() else {
        return Err(anyhow!("cigar buffer length {} is not a multiple of 4", cigar_bytes.len()));
    };
    let mut span: u32 = 0;
    for chunk in chunks {
        let word = u32::from_le_bytes(*chunk);
        let kind = code_to_kind(word & 0xF);
        let len = word >> 4;
        if consumes_reference(kind) {
            span = span.checked_add(len).ok_or_else(|| anyhow!("cigar ref span overflow"))?;
        }
    }
    Ok(span)
}

/// Compute the 1-based inclusive alignment end from `alignment_start`
/// plus the reference span of a packed BAM CIGAR.
///
/// Returns `Ok(None)` when:
/// - the record has no `alignment_start` (unmapped),
/// - the CIGAR has no reference-consuming ops (`span == 0`), or
/// - `start + span - 1` overflows `usize` or doesn't fit
///   `Position::new` (only reachable on synthetic / corrupted inputs;
///   real genomic coordinates fit in u32 with room to spare).
///
/// The mate buffer's `probe`/`accept` paths interpret `None` as
/// "alignment end unknown → don't try to overlap-buffer this read",
/// which is the safe default.
fn compute_alignment_end(pos: Option<Position>, cigar_bytes: &[u8]) -> Result<Option<Position>> {
    let Some(start) = pos else { return Ok(None) };
    let span = cigar_reference_span(cigar_bytes)?;
    if span == 0 {
        return Ok(None);
    }
    let end_1based = usize::from(start)
        .checked_add(span as usize)
        .and_then(|v| v.checked_sub(1))
        .and_then(Position::new);
    Ok(end_1based)
}

/// Map a packed BAM CIGAR op code (0..=8) to a noodles `Kind`.
///
/// Panics on out-of-range codes (9..=15). All BAM read paths run
/// [`validate_cigar_codes`] at refresh time, so by the time a record
/// reaches downstream consumers the codes are guaranteed valid; a panic
/// here means an invariant violation (someone is iterating raw bytes
/// without going through `refresh_cache`).
#[inline]
fn code_to_kind(code: u32) -> Kind {
    match code {
        0 => Kind::Match,
        1 => Kind::Insertion,
        2 => Kind::Deletion,
        3 => Kind::Skip,
        4 => Kind::SoftClip,
        5 => Kind::HardClip,
        6 => Kind::Pad,
        7 => Kind::SequenceMatch,
        8 => Kind::SequenceMismatch,
        n => {
            panic!("invalid BAM CIGAR op code: {n} (validate_cigar_codes should have rejected it)")
        }
    }
}

#[inline]
fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

// ── Aux scan + parse ──

/// Scan `aux_bytes` once. For tags in `values` the payload is parsed
/// and stored in `dst`; for tags in `presence` only the tag key is
/// recorded in `present`, skipping the per-record allocation that
/// `String`/`Hex` payloads would otherwise carry. Tags outside both
/// sets are skipped in place — their value width is computed from the
/// type byte so the scanner doesn't allocate or fully decode them.
///
/// `dst` and `present` are not cleared before insertion — callers
/// should clear them first if a fresh slot is wanted.
///
/// Returns the number of tags matched (in either set).
///
/// Array-typed tags (`B:`) are currently skipped even if requested in
/// `values`; their `dst` slot stays empty. `presence` lookups still
/// register a hit because the scanner can determine existence without
/// parsing the payload.
///
/// # Errors
/// Returns `Err` if the aux byte stream is malformed (unexpected EOF,
/// unknown type byte, etc.).
fn scan_specific(
    aux_bytes: &[u8],
    values: &BTreeSet<TagKey>,
    presence: &BTreeSet<TagKey>,
    dst: &mut AuxTagStore,
    present: &mut SmallVec<[TagKey; 4]>,
) -> Result<usize> {
    let mut pos = 0;
    let mut found = 0;
    while pos + 3 <= aux_bytes.len() {
        let tag: TagKey = [aux_bytes[pos], aux_bytes[pos + 1]];
        let type_byte = aux_bytes[pos + 2];
        pos += 3;

        let value_len = aux_value_len(aux_bytes, pos, type_byte)?;
        let end = pos + value_len;
        if end > aux_bytes.len() {
            anyhow::bail!("aux value of type '{}' extends past end of buffer", type_byte as char);
        }

        if values.contains(&tag) {
            if let Some(value) = parse_aux_value(type_byte, &aux_bytes[pos..end])? {
                // `parse_aux_value` returns `None` for types we don't yet
                // support (e.g. `B:`) — skip those silently.
                dst.insert(tag, value);
                found += 1;
            }
        } else if presence.contains(&tag) && !present.contains(&tag) {
            present.push(tag);
            found += 1;
        }
        pos = end;
    }
    if pos != aux_bytes.len() {
        anyhow::bail!("trailing garbage in aux buffer: {} stray bytes", aux_bytes.len() - pos);
    }
    Ok(found)
}

/// Variant of [`scan_specific`] that decodes *every* tag in the aux
/// block into `dst`.
///
/// `dst` is not cleared before insertion — callers should call
/// [`AuxTagStore::clear`] first if they want a fresh store.
///
/// # Errors
/// Returns `Err` if the aux byte stream is malformed.
fn scan_all(aux_bytes: &[u8], dst: &mut AuxTagStore) -> Result<usize> {
    let mut pos = 0;
    let mut found = 0;
    while pos + 3 <= aux_bytes.len() {
        let tag: TagKey = [aux_bytes[pos], aux_bytes[pos + 1]];
        let type_byte = aux_bytes[pos + 2];
        pos += 3;

        let value_len = aux_value_len(aux_bytes, pos, type_byte)?;
        let end = pos + value_len;
        if end > aux_bytes.len() {
            anyhow::bail!("aux value of type '{}' extends past end of buffer", type_byte as char);
        }

        if let Some(value) = parse_aux_value(type_byte, &aux_bytes[pos..end])? {
            dst.insert(tag, value);
            found += 1;
        }
        pos = end;
    }
    if pos != aux_bytes.len() {
        anyhow::bail!("trailing garbage in aux buffer: {} stray bytes", aux_bytes.len() - pos);
    }
    Ok(found)
}

/// Compute the byte length of an aux value given its type byte and the
/// current position within the buffer. For fixed-width types this is a
/// constant; for `Z` / `H` we scan to the NUL terminator; for `B` we
/// read the element type and count.
fn aux_value_len(aux_bytes: &[u8], pos: usize, type_byte: u8) -> Result<usize> {
    Ok(match type_byte {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => {
            let terminator = aux_bytes[pos..]
                .iter()
                .position(|&b| b == 0)
                .ok_or_else(|| anyhow!("aux string missing NUL terminator"))?;
            terminator + 1 // include the NUL
        }
        b'B' => {
            if pos + 5 > aux_bytes.len() {
                anyhow::bail!("aux B-array header truncated");
            }
            let elem_ty = aux_bytes[pos];
            let count = u32::from_le_bytes(aux_bytes[pos + 1..pos + 5].try_into().unwrap());
            let elem_size: usize = match elem_ty {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => anyhow::bail!("unsupported aux B-array element type '{}'", elem_ty as char),
            };
            // Hostile inputs can claim count = u32::MAX; checked_mul keeps
            // us from triggering a multi-GB allocation downstream when the
            // value-bytes window is materialised.
            let payload_len = (count as usize).checked_mul(elem_size).ok_or_else(|| {
                anyhow!("aux B-array size overflow: count={count}, elem_size={elem_size}")
            })?;
            payload_len.checked_add(5).ok_or_else(|| anyhow!("aux B-array length overflow"))?
        }
        _ => anyhow::bail!("unknown aux type byte '{}' ({:#x})", type_byte as char, type_byte),
    })
}

/// Parse the value bytes for a single aux field into an [`AuxValue`].
/// Returns `Ok(None)` for currently-unsupported types (just `B:`).
fn parse_aux_value(type_byte: u8, value: &[u8]) -> Result<Option<AuxValue>> {
    Ok(Some(match type_byte {
        b'A' => AuxValue::Character(value[0]),
        b'c' => AuxValue::Int(i64::from(i8::from_le_bytes([value[0]]))),
        b'C' => AuxValue::Int(i64::from(value[0])),
        b's' => AuxValue::Int(i64::from(i16::from_le_bytes(value.try_into().unwrap()))),
        b'S' => AuxValue::Int(i64::from(u16::from_le_bytes(value.try_into().unwrap()))),
        b'i' => AuxValue::Int(i64::from(i32::from_le_bytes(value.try_into().unwrap()))),
        b'I' => AuxValue::Int(i64::from(u32::from_le_bytes(value.try_into().unwrap()))),
        b'f' => AuxValue::Float(f32::from_le_bytes(value.try_into().unwrap())),
        b'Z' => {
            // Drop the trailing NUL — it's an on-wire artefact, not
            // part of the semantic string.
            let s = if value.last() == Some(&0) { &value[..value.len() - 1] } else { value };
            AuxValue::String(s.to_vec())
        }
        b'H' => {
            let s = if value.last() == Some(&0) { &value[..value.len() - 1] } else { value };
            AuxValue::Hex(s.to_vec())
        }
        b'B' => return Ok(None), // unsupported
        _ => anyhow::bail!("unknown aux type byte '{}' ({:#x})", type_byte as char, type_byte),
    }))
}

/// Convert a noodles `record_buf::data::field::Value` into our [`AuxValue`].
///
/// Integer widths collapse into a single `Int(i64)` variant per our
/// design. Returns `None` for `Array`-typed values, which we don't yet
/// support (callers see this as "tag missing" — consistent with the
/// BAM scanner).
fn record_buf_value_to_aux(v: &sam::alignment::record_buf::data::field::Value) -> Option<AuxValue> {
    use sam::alignment::record_buf::data::field::Value;
    Some(match v {
        Value::Character(c) => AuxValue::Character(*c),
        Value::Int8(n) => AuxValue::Int(i64::from(*n)),
        Value::UInt8(n) => AuxValue::Int(i64::from(*n)),
        Value::Int16(n) => AuxValue::Int(i64::from(*n)),
        Value::UInt16(n) => AuxValue::Int(i64::from(*n)),
        Value::Int32(n) => AuxValue::Int(i64::from(*n)),
        Value::UInt32(n) => AuxValue::Int(i64::from(*n)),
        Value::Float(f) => AuxValue::Float(*f),
        Value::String(s) => AuxValue::String(s.to_vec()),
        Value::Hex(s) => AuxValue::Hex(s.to_vec()),
        // Array-typed aux values aren't modelled in `AuxValue` yet;
        // signal "not available" the same way the BAM scanner does for
        // `B:*`.
        Value::Array(_) => return None,
    })
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind as OpKind;
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // ── Helpers ──

    fn make_record_buf_with_cigar(ops: Vec<Op>) -> RecordBuf {
        let cigar: Cigar = ops.into_iter().collect();
        RecordBuf::builder()
            .set_name(b"r1".to_vec())
            .set_flags(Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(101).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(201).unwrap())
            .set_template_length(200)
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 100]))
            .set_quality_scores(QualityScores::from(vec![30u8; 100]))
            .build()
    }

    /// Inverse of `code_to_kind` — BAM packed op code for a given
    /// `Kind`. Test-only; production code never re-encodes.
    fn kind_code(kind: Kind) -> u32 {
        match kind {
            Kind::Match => 0,
            Kind::Insertion => 1,
            Kind::Deletion => 2,
            Kind::Skip => 3,
            Kind::SoftClip => 4,
            Kind::HardClip => 5,
            Kind::Pad => 6,
            Kind::SequenceMatch => 7,
            Kind::SequenceMismatch => 8,
        }
    }

    /// Encode a list of (tag, "<type><value>") pairs into a BAM aux
    /// byte block. Type codes handled: `i`, `C`, `f`, `Z`, `A`.
    fn encode_aux(fields: &[(TagKey, &str)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (tag, spec) in fields {
            out.extend_from_slice(tag);
            let (type_byte, value) = spec.split_at(1);
            let tb = type_byte.as_bytes()[0];
            out.push(tb);
            match tb {
                b'i' => {
                    let n: i32 = value.parse().unwrap();
                    out.extend_from_slice(&n.to_le_bytes());
                }
                b'C' => {
                    let n: u8 = value.parse().unwrap();
                    out.push(n);
                }
                b'f' => {
                    let f: f32 = value.parse().unwrap();
                    out.extend_from_slice(&f.to_le_bytes());
                }
                b'Z' => {
                    out.extend_from_slice(value.as_bytes());
                    out.push(0);
                }
                b'A' => {
                    out.push(value.as_bytes()[0]);
                }
                _ => panic!("test encoder doesn't handle type '{}'", tb as char),
            }
        }
        out
    }

    // ── RikerRecord / FallbackRec ──

    #[test]
    fn fallback_scalar_accessors_are_infallible() {
        let ops = vec![Op::new(OpKind::Match, 100)];
        let buf = make_record_buf_with_cigar(ops);
        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));
        assert!(record.flags().is_segmented());
        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.alignment_start().unwrap().get(), 101);
        assert_eq!(record.template_length(), 200);
        assert_eq!(record.name().unwrap(), b"r1");
    }

    #[test]
    fn fallback_alignment_end_computed_from_ref_span() {
        // 50M + 5I + 45M → ref span = 50 + 45 = 95. pos=101 → end=195.
        let ops = vec![
            Op::new(OpKind::Match, 50),
            Op::new(OpKind::Insertion, 5),
            Op::new(OpKind::Match, 45),
        ];
        let buf = make_record_buf_with_cigar(ops);
        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));
        assert_eq!(record.alignment_end().unwrap().get(), 195);
    }

    #[test]
    fn fallback_sequence_is_always_populated() {
        // RecordBuf already holds ASCII bases, so Fallback's sequence
        // accessor returns the bases regardless of whether any filler
        // has run on the BAM side.
        let ops = vec![Op::new(OpKind::Match, 100)];
        let buf = make_record_buf_with_cigar(ops);
        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));

        let seq = record.sequence();
        assert_eq!(seq.len(), 100);
        assert!(seq.iter().all(|&b| b == b'A'));
        assert_eq!(record.sequence_len(), 100);
    }

    #[test]
    fn fallback_cigar_iteration_matches_input() {
        let ops = vec![
            Op::new(OpKind::Match, 50),
            Op::new(OpKind::Insertion, 5),
            Op::new(OpKind::Match, 45),
        ];
        let buf = make_record_buf_with_cigar(ops.clone());
        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));

        let round_tripped: Vec<Op> = record.cigar_ops().collect();
        assert_eq!(round_tripped, ops);
    }

    #[test]
    fn fallback_aux_tag_converts_int_and_string_on_the_fly() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        let mut buf = make_record_buf_with_cigar(vec![Op::new(OpKind::Match, 100)]);
        // `ALIGNMENT_HIT_COUNT` is the `NH` tag; `MISMATCHED_POSITIONS`
        // is `MD`. Insert via noodles' `Tag` constants; look up via the
        // raw wire bytes so the public API is the one under test.
        buf.data_mut().insert(Tag::ALIGNMENT_HIT_COUNT, Value::Int32(3));
        buf.data_mut().insert(Tag::MISMATCHED_POSITIONS, Value::String(b"50ACGT".into()));

        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));

        assert_eq!(record.aux_tag(*b"NH").and_then(|v| v.as_int()), Some(3));
        let md = record.aux_tag(*b"MD").unwrap();
        assert_eq!(md.as_str().unwrap(), "50ACGT");
    }

    #[test]
    fn fallback_aux_tag_returns_none_for_missing() {
        let buf = make_record_buf_with_cigar(vec![Op::new(OpKind::Match, 100)]);
        let record = RikerRecord::Fallback(FallbackRec::from_record_buf(buf));
        assert!(record.aux_tag(*b"NM").is_none());
    }

    // ── CIGAR helpers ──

    #[test]
    fn cigar_reference_span_handles_deletion_and_skip() {
        // 10M 5D 10N 10M → ref span = 10 + 5 + 10 + 10 = 35.
        // BAM packs each op as `(len << 4) | kind`; kind 0 = Match so
        // the M ops are just `10u32 << 4`.
        let packed = [
            10u32 << 4,       // 10M
            (5u32 << 4) | 2,  // 5D
            (10u32 << 4) | 3, // 10N
            10u32 << 4,       // 10M
        ];
        let bytes: Vec<u8> = packed.iter().flat_map(|w| w.to_le_bytes()).collect();
        assert_eq!(cigar_reference_span(&bytes).unwrap(), 35);
    }

    #[test]
    fn cigar_reference_span_rejects_non_multiple_of_four() {
        let bytes = vec![0u8; 5];
        assert!(cigar_reference_span(&bytes).is_err());
    }

    #[test]
    fn validate_cigar_codes_rejects_out_of_range_op() {
        // (10 << 4) | 9 — code 9 is reserved by the BAM spec.
        let word = (10u32 << 4) | 9;
        let bytes: Vec<u8> = word.to_le_bytes().to_vec();
        let err = validate_cigar_codes(&bytes).unwrap_err();
        assert!(
            err.to_string().contains("invalid BAM CIGAR op code"),
            "expected op-code error, got: {err}"
        );
    }

    #[test]
    fn validate_cigar_codes_accepts_all_valid_codes() {
        let bytes: Vec<u8> = (0u32..=8).flat_map(|c| ((1u32 << 4) | c).to_le_bytes()).collect();
        validate_cigar_codes(&bytes).unwrap();
    }

    #[test]
    fn kind_code_round_trips() {
        for kind in [
            OpKind::Match,
            OpKind::Insertion,
            OpKind::Deletion,
            OpKind::Skip,
            OpKind::SoftClip,
            OpKind::HardClip,
            OpKind::Pad,
            OpKind::SequenceMatch,
            OpKind::SequenceMismatch,
        ] {
            assert_eq!(code_to_kind(kind_code(kind)), kind);
        }
    }

    // ── Aux scanner ──

    /// Convenience wrapper for the scanner with no presence-only set —
    /// keeps the per-test boilerplate readable.
    fn scan_values(aux: &[u8], wanted: &BTreeSet<TagKey>, dst: &mut AuxTagStore) -> Result<usize> {
        let presence: BTreeSet<TagKey> = BTreeSet::new();
        let mut present: SmallVec<[TagKey; 4]> = SmallVec::new();
        scan_specific(aux, wanted, &presence, dst, &mut present)
    }

    #[test]
    fn scan_extracts_requested_int_tag_ignoring_others() {
        let aux = encode_aux(&[(*b"NM", "i3"), (*b"MD", "Z50"), (*b"AS", "i100")]);
        let wanted: BTreeSet<TagKey> = [*b"NM"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        let found = scan_values(&aux, &wanted, &mut dst).unwrap();
        assert_eq!(found, 1);
        assert_eq!(dst.get(*b"NM").unwrap().as_int(), Some(3));
        assert!(dst.get(*b"MD").is_none());
        assert!(dst.get(*b"AS").is_none());
    }

    #[test]
    fn scan_extracts_multiple_tags() {
        let aux = encode_aux(&[(*b"NM", "i3"), (*b"MD", "Z50ACGT"), (*b"AS", "i100")]);
        let wanted: BTreeSet<TagKey> = [*b"NM", *b"AS"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        let found = scan_values(&aux, &wanted, &mut dst).unwrap();
        assert_eq!(found, 2);
        assert_eq!(dst.get(*b"NM").unwrap().as_int(), Some(3));
        assert_eq!(dst.get(*b"AS").unwrap().as_int(), Some(100));
    }

    #[test]
    fn scan_presence_only_skips_value_decode_for_string_tag() {
        // Z-typed tag would otherwise allocate via `String(Vec<u8>)`.
        // Presence-only path should record the tag without parsing it.
        let aux = encode_aux(&[(*b"SA", "Z1,12345,+,40M,30,2;"), (*b"NM", "i7")]);
        let presence: BTreeSet<TagKey> = [*b"SA"].into_iter().collect();
        let values: BTreeSet<TagKey> = [*b"NM"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        let mut present: SmallVec<[TagKey; 4]> = SmallVec::new();
        let found = scan_specific(&aux, &values, &presence, &mut dst, &mut present).unwrap();
        assert_eq!(found, 2);
        // SA is in `present` only, not `dst`.
        assert!(present.contains(b"SA"));
        assert!(dst.get(*b"SA").is_none());
        // NM still has its value decoded.
        assert_eq!(dst.get(*b"NM").unwrap().as_int(), Some(7));
    }

    #[test]
    fn all_integer_widths_collapse_into_int() {
        // Tag-by-tag: 'c' i8 minus, 'C' u8, 's' i16, 'S' u16, 'i' i32,
        // 'I' u32.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"T1c");
        aux.push(0xFFu8); // -1 as i8
        aux.extend_from_slice(b"T2C");
        aux.push(200u8);
        aux.extend_from_slice(b"T3s");
        aux.extend_from_slice(&(-500i16).to_le_bytes());
        aux.extend_from_slice(b"T4S");
        aux.extend_from_slice(&40000u16.to_le_bytes());
        aux.extend_from_slice(b"T5i");
        aux.extend_from_slice(&(-1_000_000i32).to_le_bytes());
        aux.extend_from_slice(b"T6I");
        aux.extend_from_slice(&3_000_000_000u32.to_le_bytes());

        let wanted: BTreeSet<TagKey> =
            [*b"T1", *b"T2", *b"T3", *b"T4", *b"T5", *b"T6"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        scan_values(&aux, &wanted, &mut dst).unwrap();
        assert_eq!(dst.get(*b"T1").unwrap().as_int(), Some(-1));
        assert_eq!(dst.get(*b"T2").unwrap().as_int(), Some(200));
        assert_eq!(dst.get(*b"T3").unwrap().as_int(), Some(-500));
        assert_eq!(dst.get(*b"T4").unwrap().as_int(), Some(40000));
        assert_eq!(dst.get(*b"T5").unwrap().as_int(), Some(-1_000_000));
        assert_eq!(dst.get(*b"T6").unwrap().as_int(), Some(3_000_000_000));
    }

    #[test]
    fn string_and_char_and_float_values() {
        let aux = encode_aux(&[(*b"MD", "Z50ACGT"), (*b"XC", "AN"), (*b"XF", "f0.75")]);
        let wanted: BTreeSet<TagKey> = [*b"MD", *b"XC", *b"XF"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        scan_values(&aux, &wanted, &mut dst).unwrap();
        assert_eq!(dst.get(*b"MD").and_then(AuxValue::as_str).unwrap(), "50ACGT");
        assert_eq!(dst.get(*b"XC").and_then(AuxValue::as_char), Some(b'N'));
        let f = dst.get(*b"XF").and_then(AuxValue::as_float).unwrap();
        assert!((f - 0.75).abs() < 1e-6);
    }

    #[test]
    fn scan_skips_unsupported_array_type_silently() {
        // B:i array of 3 i32s: type_ty='i', count=3, values=1,2,3.
        // 5 header bytes + 12 data.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"BAB");
        aux.push(b'i');
        aux.extend_from_slice(&3u32.to_le_bytes());
        for v in [1i32, 2, 3] {
            aux.extend_from_slice(&v.to_le_bytes());
        }
        aux.extend_from_slice(b"NMi");
        aux.extend_from_slice(&5i32.to_le_bytes());

        let wanted: BTreeSet<TagKey> = [*b"BA", *b"NM"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        scan_values(&aux, &wanted, &mut dst).unwrap();
        // `BA` requested but unsupported type → not stored.
        assert!(dst.get(*b"BA").is_none());
        // `NM` stored fine despite the scanner walking past `BA`.
        assert_eq!(dst.get(*b"NM").unwrap().as_int(), Some(5));
    }

    #[test]
    fn aux_b_array_oversized_payload_rejected_by_scan_bounds() {
        // Hostile input: claim a B:i array with count = u32::MAX. The
        // payload length aux_value_len returns is ~16 GB on 64-bit (well
        // within usize) — what catches it is the scan_specific bounds
        // check that the value extends past the end of the buffer. The
        // critical thing is no panic, no oversize allocation, just a
        // structured error.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"BAB"); // tag + B-array type byte
        aux.push(b'i'); // element type i32 (4 bytes each)
        aux.extend_from_slice(&u32::MAX.to_le_bytes());
        // No payload. The buffer is 8 bytes total; claimed length is huge.
        let wanted: BTreeSet<TagKey> = [*b"BA"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        let err = scan_values(&aux, &wanted, &mut dst).unwrap_err();
        assert!(
            err.to_string().contains("extends past end of buffer"),
            "expected past-end error, got: {err}"
        );
    }

    #[test]
    fn malformed_aux_surface_as_error() {
        let aux = b"XXQabcd".to_vec();
        let wanted: BTreeSet<TagKey> = [*b"XX"].into_iter().collect();
        let mut dst = AuxTagStore::new();
        assert!(scan_values(&aux, &wanted, &mut dst).is_err());
    }

    #[test]
    fn aux_store_insert_replaces_on_same_tag() {
        let mut s = AuxTagStore::new();
        s.insert(*b"NM", AuxValue::Int(1));
        s.insert(*b"NM", AuxValue::Int(2));
        assert_eq!(s.len(), 1);
        assert_eq!(s.get(*b"NM").unwrap().as_int(), Some(2));
    }

    // ── RikerRecordRequirements ─────────────────────────────────────────────

    #[test]
    fn requirements_default_is_none() {
        let n = RikerRecordRequirements::default();
        assert_eq!(n, RikerRecordRequirements::NONE);
        assert!(!n.needs_sequence());
        assert!(n.aux_tags().is_none());
    }

    #[test]
    fn requirements_with_sequence_sets_flag() {
        let n = RikerRecordRequirements::NONE.with_sequence();
        assert!(n.needs_sequence());
        assert!(n.aux_tags().is_none());
    }

    #[test]
    fn requirements_with_aux_tag_transitions_none_to_specific() {
        let n = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        match n.aux_tags {
            AuxTagRequirements::Specific { values, presence } => {
                assert_eq!(values.len(), 1);
                assert!(values.contains(b"NM"));
                assert!(presence.is_empty());
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_multiple_with_aux_tag_calls_accumulate() {
        let n = RikerRecordRequirements::NONE.with_aux_tag(*b"NM").with_aux_tag(*b"MD");
        match n.aux_tags {
            AuxTagRequirements::Specific { values, presence: _ } => {
                assert_eq!(values.len(), 2);
                assert!(values.contains(b"NM"));
                assert!(values.contains(b"MD"));
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_presence_only_lands_in_presence_set() {
        let n = RikerRecordRequirements::NONE.with_aux_tag_presence(*b"SA");
        match n.aux_tags {
            AuxTagRequirements::Specific { values, presence } => {
                assert!(values.is_empty());
                assert_eq!(presence.len(), 1);
                assert!(presence.contains(b"SA"));
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_value_request_subsumes_later_presence_request() {
        // value first, presence later — should not move SA into presence.
        let n = RikerRecordRequirements::NONE.with_aux_tag(*b"SA").with_aux_tag_presence(*b"SA");
        match n.aux_tags {
            AuxTagRequirements::Specific { values, presence } => {
                assert!(values.contains(b"SA"));
                assert!(presence.is_empty());
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_value_request_evicts_earlier_presence() {
        // presence first, then upgraded to value — value wins.
        let n = RikerRecordRequirements::NONE.with_aux_tag_presence(*b"SA").with_aux_tag(*b"SA");
        match n.aux_tags {
            AuxTagRequirements::Specific { values, presence } => {
                assert!(values.contains(b"SA"));
                assert!(presence.is_empty());
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_union_value_subsumes_presence_across_collectors() {
        let a = RikerRecordRequirements::NONE.with_aux_tag(*b"SA");
        let b = RikerRecordRequirements::NONE.with_aux_tag_presence(*b"SA");
        let u = a.union(b);
        match u.aux_tags {
            AuxTagRequirements::Specific { values, presence } => {
                assert!(values.contains(b"SA"));
                assert!(presence.is_empty());
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_with_all_aux_tags_sets_all() {
        let n = RikerRecordRequirements::NONE.with_aux_tag(*b"NM").with_all_aux_tags();
        assert_eq!(n.aux_tags, AuxTagRequirements::All);
    }

    #[test]
    fn requirements_union_of_none_and_specific_is_specific() {
        let a = RikerRecordRequirements::NONE;
        let b = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        assert_eq!(a.union(b.clone()), b);
    }

    #[test]
    fn requirements_union_merges_specific_sets() {
        let a = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        let b = RikerRecordRequirements::NONE.with_aux_tag(*b"MD");
        let u = a.union(b);
        match u.aux_tags {
            AuxTagRequirements::Specific { values, presence: _ } => {
                assert!(values.contains(b"NM"));
                assert!(values.contains(b"MD"));
            }
            other => panic!("expected Specific, got {other:?}"),
        }
    }

    #[test]
    fn requirements_union_with_all_absorbs_everything() {
        let a = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        let b = RikerRecordRequirements::everything();
        assert_eq!(a.clone().union(b.clone()), b);
        assert_eq!(b.clone().union(a), b);
    }

    #[test]
    fn requirements_union_is_commutative_on_specific() {
        let a = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
        let b = RikerRecordRequirements::NONE.with_aux_tag(*b"MD");
        assert_eq!(a.clone().union(b.clone()), b.union(a));
    }

    #[test]
    fn requirements_union_is_idempotent() {
        let a = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");
        assert_eq!(a.clone().union(a.clone()), a);
    }

    #[test]
    fn requirements_sequence_flag_unions_with_or() {
        let a = RikerRecordRequirements::NONE.with_sequence();
        let b = RikerRecordRequirements::NONE;
        assert!(a.clone().union(b.clone()).needs_sequence());
        assert!(b.union(a).needs_sequence());
    }
}
