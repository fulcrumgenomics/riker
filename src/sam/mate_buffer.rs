//! Reusable mate-pair buffer for coordinate-ordered per-contig scans.
//!
//! Many per-read metrics need to resolve overlap between a read and its mate —
//! for example, to avoid double-counting depth at positions covered by both
//! reads of a pair. When records arrive in coordinate order, the mate of a
//! given read is guaranteed to come later (if at all) at a position `>=` the
//! read's own start. This module buffers records whose mate is expected later,
//! keyed by read name, and returns the cached projection when the mate arrives.
//!
//! ## Two routing APIs
//!
//! [`MateBuffer`] exposes two entry points depending on how the caller wants
//! to build its cache state:
//!
//! - [`MateBuffer::accept`] — one-step. Returns a [`MateAction`] and, for the
//!   "buffer this read" case, auto-projects the record via the [`MateCache`]
//!   trait. This fits tools like `error` that can't do meaningful work until
//!   the pair resolves — they want to cache the record, not process it yet.
//!
//! - [`MateBuffer::probe`] + [`MateBuffer::insert`] — two-step. `probe`
//!   returns a [`Peek`] describing the routing decision but does *not*
//!   construct a cache value; callers that want to compute cache state
//!   inline with per-record processing (e.g. `wgs` builds an overlap
//!   bitmap as a side-effect of its depth walk) do so between the two calls.
//!   No [`MateCache`] impl is required.
//!
//! ## Flushing
//!
//! The buffer also knows how to drain pending entries when a scan passes
//! their expected mate position (e.g. at interval boundaries):
//! [`flush_behind`] and [`clear_behind`] yield or discard entries whose mates
//! won't arrive within the remaining scan; [`flush`] and [`clear`] do the
//! same for all remaining entries.
//!
//! ## Keying contract
//!
//! The buffer is keyed on raw read name. Records without a name
//! (`record.name() == None`) cannot be paired because two unrelated unnamed
//! reads would collide on the empty key; both [`accept`] and [`probe`]
//! return `Alone` for nameless records, so unnamed reads are transparently
//! treated as singletons.
//!
//! [`accept`]: MateBuffer::accept
//! [`probe`]: MateBuffer::probe
//! [`flush_behind`]: MateBuffer::flush_behind
//! [`clear_behind`]: MateBuffer::clear_behind
//! [`flush`]: MateBuffer::flush
//! [`clear`]: MateBuffer::clear

use noodles::sam::alignment::RecordBuf;
use rustc_hash::{FxBuildHasher, FxHashMap};

/// Tool-specific projection of a [`RecordBuf`] stored in the buffer until the
/// mate arrives. Implementations should extract only the state needed for the
/// tool's pair-resolution logic.
pub trait MateCache: Sized {
    /// Build a cache entry from the record being buffered.
    fn from_record(record: &RecordBuf) -> Self;
}

/// Outcome of presenting a record to the buffer via [`MateBuffer::accept`].
pub enum MateAction<T> {
    /// No overlap possible (unpaired, unmapped mate, cross-contig, or mate
    /// falls outside this read's reference span). Caller should process the
    /// record as a single read.
    Alone,
    /// The mate is expected later within this read's span. The buffer has
    /// already stored `T::from_record(record)`. Caller decides whether to
    /// process the current record now (eager) or defer until the pair
    /// resolves (lazy).
    Buffered,
    /// The mate was previously buffered; its cache is returned. Caller should
    /// process the current record paired with the cached mate.
    PairWith(T),
}

/// Outcome of peeking at a record via [`MateBuffer::probe`] without
/// committing to a cache projection. Callers that want to build `T` inline
/// with their own per-base processing use `probe` + [`MateBuffer::insert`]
/// instead of [`MateBuffer::accept`].
pub enum Peek<T> {
    /// No overlap possible — process the record alone.
    Alone,
    /// The record will be buffered on this contig. `overlap_start` is the
    /// 0-based reference position where the arriving mate can first overlap
    /// this read (equal to the mate's expected alignment_start), and
    /// `overlap_len` is the number of reference positions in the overlap
    /// region (bounded by this read's alignment_end). Callers should size
    /// any per-position structure (e.g. a bitmap) to `overlap_len`, then
    /// follow up with [`MateBuffer::insert`] once the cache is built.
    WouldBuffer { overlap_start: u32, overlap_len: u32 },
    /// Mate was previously buffered; its cache is returned. Caller should
    /// process the current record paired with the cached mate.
    PairWith(T),
}

/// Blanket [`MateCache`] impl for [`RecordBuf`]; clones the record for
/// callers that need the full read at pair-resolution time (e.g. `error`,
/// which re-walks the cached mate's CIGAR when the pair arrives).
impl MateCache for RecordBuf {
    fn from_record(record: &RecordBuf) -> Self {
        record.clone()
    }
}

/// A mate-pair buffer keyed on read name. See module docs for the coordinate-
/// order contract and the projection model.
///
/// Callers that have a tool-specific projection type can build `T` inline via
/// [`probe`] + [`insert`] without implementing [`MateCache`]. The trait is
/// only required by [`accept`], which auto-projects via `T::from_record`.
///
/// [`probe`]: MateBuffer::probe
/// [`insert`]: MateBuffer::insert
/// [`accept`]: MateBuffer::accept
pub struct MateBuffer<T> {
    buffer: FxHashMap<Vec<u8>, Entry<T>>,
}

/// Internal buffer entry. Stores the cache value alongside the mate's expected
/// reference coordinates so `flush_behind` / `clear_behind` can filter without
/// requiring `T` to expose coordinates.
struct Entry<T> {
    value: T,
    mate_ref_id: usize,
    /// 0-based mate alignment start.
    mate_pos: u32,
}

impl<T> MateBuffer<T> {
    /// Initial buffer capacity; sized to hold a few thousand pending pairs
    /// without re-allocating on typical WGS / targeted inputs.
    const INITIAL_CAPACITY: usize = 4096;

    /// Create an empty buffer.
    #[must_use]
    pub fn new() -> Self {
        Self { buffer: FxHashMap::with_capacity_and_hasher(Self::INITIAL_CAPACITY, FxBuildHasher) }
    }

    /// Two-step variant of [`accept`]: decide how to route `record` without
    /// constructing or inserting a cache value. On `WouldBuffer` the caller
    /// is responsible for following up with [`insert`] to deposit the value
    /// (typically after building it inline with per-base processing).
    ///
    /// [`accept`]: MateBuffer::accept
    /// [`insert`]: MateBuffer::insert
    #[allow(clippy::cast_possible_truncation, reason = "1-based genomic positions fit in u32")]
    pub fn probe(&mut self, record: &RecordBuf) -> Peek<T> {
        let flags = record.flags();

        if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
            return Peek::Alone;
        }
        let Some(ref_id) = record.reference_sequence_id() else {
            return Peek::Alone;
        };
        if Some(ref_id) != record.mate_reference_sequence_id() {
            return Peek::Alone;
        }
        // The buffer is keyed on read name; unnamed records can't be paired
        // because two unrelated unnamed reads would collide on the empty key.
        let Some(name) = record.name() else {
            return Peek::Alone;
        };

        let name_bytes: &[u8] = name.as_ref();
        if let Some(entry) = self.buffer.remove(name_bytes) {
            return Peek::PairWith(entry.value);
        }

        let read_start = record.alignment_start().map_or(0, |p| p.get());
        let read_end = record.alignment_end().map_or(0, |e| e.get());
        let mate_start_1based = record.mate_alignment_start().map_or(0, |p| p.get());
        if mate_start_1based >= read_start && mate_start_1based <= read_end {
            let overlap_start = mate_start_1based.saturating_sub(1) as u32;
            // read_end is 1-based inclusive; overlap region spans mate_start..=read_end.
            let overlap_len = (read_end + 1 - mate_start_1based) as u32;
            return Peek::WouldBuffer { overlap_start, overlap_len };
        }

        Peek::Alone
    }

    /// Deposit a caller-constructed cache `value` for `record`. Intended as
    /// the follow-up to a [`probe`] call that returned [`Peek::WouldBuffer`].
    /// The record's mate coordinates are re-read here so the entry can be
    /// filtered by [`flush_behind`] / [`clear_behind`] later.
    ///
    /// Silently no-ops if `record` is missing a reference sequence id or read
    /// name, matching [`probe`]'s contract — both cases return [`Peek::Alone`]
    /// there, so a paired probe/insert sequence will never reach this path.
    ///
    /// [`probe`]: MateBuffer::probe
    /// [`flush_behind`]: MateBuffer::flush_behind
    /// [`clear_behind`]: MateBuffer::clear_behind
    #[allow(clippy::cast_possible_truncation, reason = "1-based genomic positions fit in u32")]
    pub fn insert(&mut self, record: &RecordBuf, value: T) {
        // `insert` is contracted as a follow-up to a `probe` that returned
        // `WouldBuffer`, which guarantees both fields are present. A debug
        // assertion catches callers that skip `probe` and insert cold.
        debug_assert!(record.reference_sequence_id().is_some() && record.name().is_some());
        let Some(ref_id) = record.reference_sequence_id() else { return };
        let Some(name) = record.name() else { return };
        let mate_pos_1based = record.mate_alignment_start().map_or(0, |p| p.get());
        let entry = Entry {
            value,
            mate_ref_id: ref_id,
            mate_pos: mate_pos_1based.saturating_sub(1) as u32,
        };
        self.buffer.insert(<[u8]>::to_vec(name), entry);
    }

    /// Drain and return entries whose expected mate position is strictly before
    /// `(ref_id, pos)`. Typical use: at an interval boundary, yield orphans
    /// whose mates won't appear within the remaining scan.
    pub fn flush_behind(&mut self, ref_id: usize, pos: u32) -> Vec<T> {
        self.buffer
            .extract_if(|_, entry| Self::is_behind(entry, ref_id, pos))
            .map(|(_, e)| e.value)
            .collect()
    }

    /// Drain and return all remaining buffered entries.
    pub fn flush(&mut self) -> Vec<T> {
        self.buffer.drain().map(|(_, e)| e.value).collect()
    }

    /// Drop (without returning) entries whose expected mate position is
    /// strictly before `(ref_id, pos)`. Equivalent to `flush_behind` but
    /// discards the drained values.
    pub fn clear_behind(&mut self, ref_id: usize, pos: u32) {
        self.buffer.retain(|_, entry| !Self::is_behind(entry, ref_id, pos));
    }

    /// Drop all buffered entries.
    pub fn clear(&mut self) {
        self.buffer.clear();
    }

    /// Number of buffered entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.buffer.len()
    }

    /// True if no entries are buffered.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }

    fn is_behind(entry: &Entry<T>, ref_id: usize, pos: u32) -> bool {
        entry.mate_ref_id < ref_id || (entry.mate_ref_id == ref_id && entry.mate_pos < pos)
    }
}

impl<T: MateCache> MateBuffer<T> {
    /// Decide how to process `record`:
    /// 1. Unpaired, unmapped, mate-unmapped, or cross-contig → [`MateAction::Alone`].
    /// 2. Mate already buffered (by name) → remove and return [`MateAction::PairWith`].
    /// 3. Mate's alignment start falls within this read's reference span →
    ///    project via `T::from_record` into the buffer and return
    ///    [`MateAction::Buffered`].
    /// 4. Otherwise (mate starts before or outside this read's span) →
    ///    [`MateAction::Alone`].
    #[allow(clippy::cast_possible_truncation, reason = "1-based genomic positions fit in u32")]
    pub fn accept(&mut self, record: &RecordBuf) -> MateAction<T> {
        let flags = record.flags();

        // Fast reject: can this read even have an overlapping mate?
        if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
            return MateAction::Alone;
        }
        let Some(ref_id) = record.reference_sequence_id() else {
            return MateAction::Alone;
        };
        if Some(ref_id) != record.mate_reference_sequence_id() {
            return MateAction::Alone;
        }
        // The buffer is keyed on read name; unnamed records can't be paired
        // because two unrelated unnamed reads would collide on the empty key.
        let Some(name) = record.name() else {
            return MateAction::Alone;
        };

        // Mate already in the buffer? Remove and return the pair.
        let name_bytes: &[u8] = name.as_ref();
        if let Some(entry) = self.buffer.remove(name_bytes) {
            return MateAction::PairWith(entry.value);
        }

        // Buffer iff the mate's start falls within this read's reference span.
        let read_start = record.alignment_start().map_or(0, |p| p.get());
        let read_end = record.alignment_end().map_or(0, |e| e.get());
        let mate_start_1based = record.mate_alignment_start().map_or(0, |p| p.get());
        if mate_start_1based >= read_start && mate_start_1based <= read_end {
            let entry = Entry {
                value: T::from_record(record),
                mate_ref_id: ref_id,
                // Convert to 0-based; `saturating_sub` guards the record-without-
                // alignment-start case (in which `read_end >= mate_start_1based`
                // also evaluated above is only reachable at 0, producing 0 here).
                mate_pos: mate_start_1based.saturating_sub(1) as u32,
            };
            self.buffer.insert(<[u8]>::to_vec(name), entry);
            return MateAction::Buffered;
        }

        MateAction::Alone
    }
}

impl<T> Default for MateBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::{Flags, cigar::op::Kind};
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    // ── Helpers ──────────────────────────────────────────────────────────────

    /// Build a minimal paired, mapped, mate-mapped `RecordBuf` on a single
    /// contig with a 100M CIGAR and the requested positions.
    fn paired_record(
        name: &[u8],
        ref_id: usize,
        pos_1based: usize,
        mate_pos_1based: usize,
    ) -> RecordBuf {
        let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
        RecordBuf::builder()
            .set_name(name.to_vec())
            .set_flags(Flags::SEGMENTED)
            .set_reference_sequence_id(ref_id)
            .set_mate_reference_sequence_id(ref_id)
            .set_alignment_start(Position::new(pos_1based).expect("pos"))
            .set_mate_alignment_start(Position::new(mate_pos_1based).expect("mate_pos"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 100]))
            .set_quality_scores(QualityScores::from(vec![30u8; 100]))
            .build()
    }

    /// Build an unpaired record (no SEGMENTED flag).
    fn unpaired_record(pos_1based: usize) -> RecordBuf {
        let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        RecordBuf::builder()
            .set_name(b"unpaired".to_vec())
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos_1based).expect("pos"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 50]))
            .set_quality_scores(QualityScores::from(vec![30u8; 50]))
            .build()
    }

    /// A minimal cache that just remembers the mate's read name for assertions.
    #[derive(Debug)]
    struct NameCache(Vec<u8>);
    impl MateCache for NameCache {
        fn from_record(record: &RecordBuf) -> Self {
            Self(record.name().map_or_else(Vec::new, |n| <[u8]>::to_vec(n)))
        }
    }

    // ── accept: reject paths ─────────────────────────────────────────────────

    #[test]
    fn test_accept_unpaired_is_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = unpaired_record(100);
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_accept_unmapped_is_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let mut r = paired_record(b"q", 0, 100, 120);
        *r.flags_mut() = Flags::SEGMENTED | Flags::UNMAPPED;
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_accept_mate_unmapped_is_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let mut r = paired_record(b"q", 0, 100, 120);
        *r.flags_mut() = Flags::SEGMENTED | Flags::MATE_UNMAPPED;
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_accept_cross_contig_is_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let mut r = paired_record(b"q", 0, 100, 120);
        *r.mate_reference_sequence_id_mut() = Some(1);
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    /// Build a paired record with no `set_name` call so `record.name()`
    /// returns `None`. Used to exercise the nameless-record branch.
    fn nameless_paired_record(
        ref_id: usize,
        pos_1based: usize,
        mate_pos_1based: usize,
    ) -> RecordBuf {
        let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
        RecordBuf::builder()
            .set_flags(Flags::SEGMENTED)
            .set_reference_sequence_id(ref_id)
            .set_mate_reference_sequence_id(ref_id)
            .set_alignment_start(Position::new(pos_1based).expect("pos"))
            .set_mate_alignment_start(Position::new(mate_pos_1based).expect("mate_pos"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 100]))
            .set_quality_scores(QualityScores::from(vec![30u8; 100]))
            .build()
    }

    #[test]
    fn test_accept_nameless_records_are_alone() {
        // Two nameless paired records with overlapping mate positions would
        // collide on the empty key if we let them into the buffer; make sure
        // both are routed to Alone so they're processed independently.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r1 = nameless_paired_record(0, 100, 150);
        let r2 = nameless_paired_record(0, 100, 150);
        assert!(matches!(buf.accept(&r1), MateAction::Alone));
        assert!(matches!(buf.accept(&r2), MateAction::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_probe_nameless_records_are_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = nameless_paired_record(0, 100, 150);
        assert!(matches!(buf.probe(&r), Peek::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_accept_mate_before_read_is_alone() {
        // mate_start (50) < read_start (100) — shouldn't buffer; the mate has
        // already gone past us (or we've gone past it) and won't be seen later.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 50);
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_accept_mate_past_read_end_is_alone() {
        // mate_start (500) > read_end (199) — no overlap possible, don't buffer.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 500);
        assert!(matches!(buf.accept(&r), MateAction::Alone));
        assert!(buf.is_empty());
    }

    // ── accept: buffer / pair paths ──────────────────────────────────────────

    #[test]
    fn test_accept_mate_in_span_is_buffered() {
        // mate_start (150) within read's [100, 199] → buffer.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 150);
        assert!(matches!(buf.accept(&r), MateAction::Buffered));
        assert_eq!(buf.len(), 1);
    }

    #[test]
    fn test_accept_mate_at_read_start_is_buffered() {
        // mate_start == read_start (boundary): within span (inclusive).
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 100);
        assert!(matches!(buf.accept(&r), MateAction::Buffered));
    }

    #[test]
    fn test_accept_mate_at_read_end_is_buffered() {
        // mate_start at exact read_end (1-based inclusive): within span.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 199); // 100-bp read spans 100..=199
        assert!(matches!(buf.accept(&r), MateAction::Buffered));
    }

    #[test]
    fn test_accept_pair_returns_cached_mate() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let first = paired_record(b"qname", 0, 100, 150);
        assert!(matches!(buf.accept(&first), MateAction::Buffered));

        let second = paired_record(b"qname", 0, 150, 100);
        match buf.accept(&second) {
            MateAction::PairWith(cache) => assert_eq!(cache.0, b"qname"),
            other => panic!("expected PairWith, got {:?}", other_kind(&other)),
        }
        assert!(buf.is_empty());
    }

    // ── flush / clear semantics ──────────────────────────────────────────────

    #[test]
    fn test_flush_behind_yields_orphans_on_earlier_contig() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q1", 0, 100, 150);
        let _ = buf.accept(&r);

        // We're now scanning contig 1 — anything on contig 0 is behind.
        let orphans = buf.flush_behind(1, 0);
        assert_eq!(orphans.len(), 1);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_flush_behind_yields_orphans_past_pos_on_same_contig() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        // Two buffered reads with different mate positions: 150 (0-based 149)
        // and 750 (0-based 749).
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150));
        let _ = buf.accept(&paired_record(b"q2", 0, 700, 750));

        // Scanning past 500: q1's mate (149) is behind, q2's (749) is not.
        let orphans = buf.flush_behind(0, 500);
        assert_eq!(orphans.len(), 1);
        assert_eq!(buf.len(), 1);
    }

    #[test]
    fn test_flush_behind_empty_when_none_match() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150));

        // Scanning at pos 0 of the same contig — nothing is strictly behind.
        let orphans = buf.flush_behind(0, 0);
        assert!(orphans.is_empty());
        assert_eq!(buf.len(), 1);
    }

    #[test]
    fn test_flush_behind_is_strictly_before() {
        // Entry with mate_pos == 149 (0-based) should NOT be flushed when
        // pos == 149, only when pos > 149. `flush_behind` is documented as
        // "strictly before" the given position.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150)); // mate_pos_0based = 149

        let stay = buf.flush_behind(0, 149);
        assert!(stay.is_empty());
        assert_eq!(buf.len(), 1);

        let gone = buf.flush_behind(0, 150);
        assert_eq!(gone.len(), 1);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_flush_drains_all() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150));
        let _ = buf.accept(&paired_record(b"q2", 0, 200, 250));
        let all = buf.flush();
        assert_eq!(all.len(), 2);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_clear_behind_discards_without_yield() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150));
        let _ = buf.accept(&paired_record(b"q2", 0, 700, 750));

        buf.clear_behind(0, 500);
        assert_eq!(buf.len(), 1);
    }

    #[test]
    fn test_clear_empties_buffer() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let _ = buf.accept(&paired_record(b"q1", 0, 100, 150));
        let _ = buf.accept(&paired_record(b"q2", 0, 200, 250));
        buf.clear();
        assert!(buf.is_empty());
    }

    // ── RecordBuf MateCache impl ─────────────────────────────────────────────

    #[test]
    fn test_recordbuf_cache_roundtrip() {
        let mut buf: MateBuffer<RecordBuf> = MateBuffer::new();
        let first = paired_record(b"qname", 0, 100, 150);
        assert!(matches!(buf.accept(&first), MateAction::Buffered));

        let second = paired_record(b"qname", 0, 150, 100);
        match buf.accept(&second) {
            MateAction::PairWith(cached) => {
                let name: &[u8] = cached.name().unwrap();
                assert_eq!(name, b"qname");
                assert_eq!(cached.alignment_start().unwrap().get(), 100);
            }
            _ => panic!("expected PairWith"),
        }
    }

    /// Debug-formatting helper since `MateAction<T>` doesn't require `Debug`.
    fn other_kind<T>(action: &MateAction<T>) -> &'static str {
        match action {
            MateAction::Alone => "Alone",
            MateAction::Buffered => "Buffered",
            MateAction::PairWith(_) => "PairWith",
        }
    }

    // ── probe / insert (two-step API) ────────────────────────────────────────

    #[test]
    fn test_probe_unpaired_is_alone() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = unpaired_record(100);
        assert!(matches!(buf.probe(&r), Peek::Alone));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_probe_mate_in_span_reports_overlap_coords() {
        // 100-bp read at pos 100 (1-based inclusive: 100..=199), mate at 150.
        // Overlap region: mate_start..=read_end == 150..=199 in 1-based,
        // i.e. overlap_start=149 (0-based), overlap_len=50.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let r = paired_record(b"q", 0, 100, 150);
        match buf.probe(&r) {
            Peek::WouldBuffer { overlap_start, overlap_len } => {
                assert_eq!(overlap_start, 149);
                assert_eq!(overlap_len, 50);
            }
            _ => panic!("expected WouldBuffer"),
        }
        // probe does not insert
        assert!(buf.is_empty());
    }

    #[test]
    fn test_probe_then_insert_makes_entry_visible_to_accept() {
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let first = paired_record(b"qn", 0, 100, 150);
        let Peek::WouldBuffer { .. } = buf.probe(&first) else { panic!("expected WouldBuffer") };
        buf.insert(&first, NameCache(b"qn".to_vec()));
        assert_eq!(buf.len(), 1);

        // Now the mate arriving via accept should retrieve the buffered entry.
        let second = paired_record(b"qn", 0, 150, 100);
        match buf.accept(&second) {
            MateAction::PairWith(cache) => assert_eq!(cache.0, b"qn"),
            _ => panic!("expected PairWith"),
        }
        assert!(buf.is_empty());
    }

    #[test]
    fn test_probe_consumes_buffered_mate() {
        // probe on the second read should find and remove the buffered mate.
        let mut buf: MateBuffer<NameCache> = MateBuffer::new();
        let first = paired_record(b"qn", 0, 100, 150);
        let Peek::WouldBuffer { .. } = buf.probe(&first) else { panic!("expected WouldBuffer") };
        buf.insert(&first, NameCache(b"qn".to_vec()));

        let second = paired_record(b"qn", 0, 150, 100);
        assert!(matches!(buf.probe(&second), Peek::PairWith(_)));
        assert!(buf.is_empty());
    }
}
