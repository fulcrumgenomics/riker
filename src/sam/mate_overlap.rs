//! First-read-wins deduplication of the reference positions two mates of a pair share.
//!
//! In coordinate order the leftmost mate of an overlapping pair arrives first. As a collector
//! walks that read it records, in an [`OverlapBitmap`], each position of the overlap region
//! where the read's base counted; [`MateBuffer`](super::mate_buffer::MateBuffer) holds the
//! bitmap until the mate arrives, and the mate's walk then skips those positions. A position
//! where the first read's base did not count (for example, it failed base quality) stays
//! available to the mate.
//!
//! The walk is generic over an [`OverlapAction`], one per routing case from
//! [`MateBuffer::probe`](super::mate_buffer::MateBuffer::probe):
//!
//! - [`Alone`] — no overlapping mate; nothing to record or skip.
//! - [`Buffering`] — the mate will overlap later; record counted positions.
//! - [`Pairing`] — the mate already passed; skip the positions it counted.
//!
//! Monomorphisation compiles each case to straight-line code, and [`OverlapAction::IS_ALONE`]
//! lets a walker take a vectorised path when there is no per-base bookkeeping to do.

use bitvec::prelude::*;

/// The positions of an overlap region at which the first-arriving read of a pair counted a
/// base. Covers `[overlap_start, overlap_start + len)`, which runs from the mate's alignment
/// start to the first read's alignment end.
#[derive(Debug, Clone)]
pub struct OverlapBitmap {
    /// 0-based reference position of the first bit.
    overlap_start: u32,
    /// Bit `i` is set when the first read counted a base at `overlap_start + i`.
    bitmap: BitVec<u64, Lsb0>,
}

impl OverlapBitmap {
    /// An empty bitmap over the `overlap_len` positions starting at 0-based `overlap_start`.
    #[must_use]
    pub fn new(overlap_start: u32, overlap_len: u32) -> Self {
        Self { overlap_start, bitmap: bitvec![u64, Lsb0; 0; overlap_len as usize] }
    }

    /// Whether the first read counted a base at `ref_pos`. Positions outside the overlap region
    /// return `false`.
    #[must_use]
    pub fn covered_at(&self, ref_pos: u32) -> bool {
        match ref_pos.checked_sub(self.overlap_start) {
            Some(index) => {
                let index = index as usize;
                index < self.bitmap.len() && self.bitmap[index]
            }
            None => false,
        }
    }

    /// Record that the first read counted a base at `ref_pos`. Positions outside the overlap
    /// region are ignored.
    pub fn mark(&mut self, ref_pos: u32) {
        if let Some(index) = ref_pos.checked_sub(self.overlap_start) {
            let index = index as usize;
            if index < self.bitmap.len() {
                self.bitmap.set(index, true);
            }
        }
    }

    /// [`Self::covered_at`] for the 64 positions from `ref_pos`: bit `i` is set when the first
    /// read counted a base at `ref_pos + i`.
    #[must_use]
    pub fn covered_bits(&self, ref_pos: u32) -> u64 {
        let offset = self.offset_of(ref_pos);
        if offset <= -64 || offset >= self.signed_len() {
            return 0;
        }
        // Storage bits past the region's end are never set (see `mark_bits`), so whole words
        // can be read.
        let words = self.bitmap.as_raw_slice();
        if offset < 0 {
            return words[0] << -offset;
        }
        let (index, shift) = word_and_shift(offset);
        let low = words[index] >> shift;
        let high = match words.get(index + 1) {
            Some(&word) if shift > 0 => word << (64 - shift),
            _ => 0,
        };
        low | high
    }

    /// [`Self::mark`] for each set bit `i` of `bits`, at `ref_pos + i`. Positions outside the
    /// overlap region are ignored.
    pub fn mark_bits(&mut self, ref_pos: u32, bits: u64) {
        let offset = self.offset_of(ref_pos);
        // Bit `i` lands inside the region when `0 <= offset + i < len`.
        let first = (-offset).clamp(0, 64);
        let end = (self.signed_len() - offset).clamp(0, 64);
        if first >= end {
            return;
        }
        let keep = if end == 64 { u64::MAX } else { (1u64 << end) - 1 } & !((1u64 << first) - 1);
        let bits = bits & keep;
        let words = self.bitmap.as_raw_mut_slice();
        if offset < 0 {
            words[0] |= bits >> -offset;
            return;
        }
        let (index, shift) = word_and_shift(offset);
        words[index] |= bits << shift;
        if shift > 0
            && let Some(word) = words.get_mut(index + 1)
        {
            *word |= bits >> (64 - shift);
        }
    }

    /// The region's length, signed for offset arithmetic.
    fn signed_len(&self) -> i64 {
        i64::try_from(self.bitmap.len()).unwrap_or(i64::MAX)
    }

    /// The signed offset of `ref_pos` from the region's start.
    fn offset_of(&self, ref_pos: u32) -> i64 {
        i64::from(ref_pos) - i64::from(self.overlap_start)
    }
}

/// The storage word holding bit `offset` (non-negative) and its position within the word.
fn word_and_shift(offset: i64) -> (usize, u32) {
    #[expect(clippy::cast_possible_truncation, clippy::cast_sign_loss, reason = "offset >= 0")]
    let offset = offset as usize;
    #[expect(clippy::cast_possible_truncation, reason = "a value below 64")]
    let shift = (offset % 64) as u32;
    (offset / 64, shift)
}

/// Per-record strategy for a walk over a read's aligned bases, capturing the two decisions the
/// three mate-routing cases disagree on. The defaults (never covered, nothing recorded) are the
/// [`Alone`] case.
pub trait OverlapAction {
    /// `true` only for [`Alone`]: no per-base mate bookkeeping, so a walker may process whole
    /// blocks at once. Resolved at compile time per action.
    const IS_ALONE: bool = false;

    /// Whether the other read of the pair already counted a base at `ref_pos`; if so the walker
    /// excludes this read's base there as overlap.
    #[inline]
    fn is_mate_covered(&self, _ref_pos: u32) -> bool {
        false
    }

    /// Record that this read counted a base at `ref_pos`, so a later-arriving mate skips it.
    #[inline]
    fn on_counted(&mut self, _ref_pos: u32) {}

    /// [`Self::is_mate_covered`] for the 64 positions from `ref_pos`: bit `i` for `ref_pos + i`.
    #[inline]
    fn mate_covered_bits(&self, _ref_pos: u32) -> u64 {
        0
    }

    /// [`Self::on_counted`] at `ref_pos + i` for each set bit `i` of `bits`.
    #[inline]
    fn on_counted_bits(&mut self, _ref_pos: u32, _bits: u64) {}
}

/// Action for a read with no overlapping mate.
pub struct Alone;

impl OverlapAction for Alone {
    const IS_ALONE: bool = true;
}

/// Action for a read whose overlapping mate arrives later: records counted positions.
pub struct Buffering<'a> {
    /// The bitmap to fill, inserted into the mate buffer after the walk.
    pub bitmap: &'a mut OverlapBitmap,
}

impl OverlapAction for Buffering<'_> {
    #[inline]
    fn on_counted(&mut self, ref_pos: u32) {
        self.bitmap.mark(ref_pos);
    }

    #[inline]
    fn on_counted_bits(&mut self, ref_pos: u32, bits: u64) {
        self.bitmap.mark_bits(ref_pos, bits);
    }
}

/// Action for the second read of an overlapping pair: skips positions the first read counted.
pub struct Pairing<'a> {
    /// The first read's bitmap, taken from the mate buffer.
    pub cached: &'a OverlapBitmap,
}

impl OverlapAction for Pairing<'_> {
    #[inline]
    fn is_mate_covered(&self, ref_pos: u32) -> bool {
        self.cached.covered_at(ref_pos)
    }

    #[inline]
    fn mate_covered_bits(&self, ref_pos: u32) -> u64 {
        self.cached.covered_bits(ref_pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bitmap(overlap_start: u32, overlap_len: u32, covered: &[u32]) -> OverlapBitmap {
        let mut bitmap = OverlapBitmap::new(overlap_start, overlap_len);
        for &ref_pos in covered {
            bitmap.mark(ref_pos);
        }
        bitmap
    }

    #[test]
    fn covered_at_reports_marked_positions_in_a_contiguous_region() {
        let cached = bitmap(100, 10, &[100, 101, 102, 103, 104, 105, 106, 107, 108, 109]);
        assert!(cached.covered_at(100));
        assert!(cached.covered_at(105));
        assert!(cached.covered_at(109));
        assert!(!cached.covered_at(99));
        assert!(!cached.covered_at(110));
    }

    #[test]
    fn covered_at_is_false_inside_a_deletion_gap() {
        // CIGAR 5M2D5M: counted at [100,104] and [107,111].
        let cached = bitmap(100, 12, &[100, 101, 102, 103, 104, 107, 108, 109, 110, 111]);
        assert!(cached.covered_at(104));
        assert!(!cached.covered_at(105));
        assert!(!cached.covered_at(106));
        assert!(cached.covered_at(107));
    }

    #[test]
    fn covered_at_is_false_where_the_first_read_did_not_count() {
        let cached = bitmap(100, 10, &[100, 101, 102, 104, 105, 106, 107, 108, 109]);
        assert!(cached.covered_at(102));
        assert!(!cached.covered_at(103));
        assert!(cached.covered_at(104));
    }

    #[test]
    fn covered_at_is_false_outside_the_region() {
        let cached = bitmap(100, 10, &[100, 109]);
        assert!(!cached.covered_at(50));
        assert!(!cached.covered_at(99));
        assert!(!cached.covered_at(110));
        assert!(!cached.covered_at(u32::MAX));
    }

    #[test]
    fn mark_outside_the_region_is_ignored() {
        let cached = bitmap(100, 10, &[99, 110, 105]);
        assert!(!cached.covered_at(99));
        assert!(!cached.covered_at(110));
        assert!(cached.covered_at(105));
    }

    #[test]
    fn buffering_marks_and_pairing_reads_back() {
        let mut first = OverlapBitmap::new(100, 10);
        Buffering { bitmap: &mut first }.on_counted(103);
        let pairing = Pairing { cached: &first };
        assert!(pairing.is_mate_covered(103));
        assert!(!pairing.is_mate_covered(104));
    }

    /// A bitmap over 130 positions from 1,000 with every third position marked.
    fn sparse_bitmap() -> OverlapBitmap {
        let marked: Vec<u32> = (1_000..1_130).filter(|pos| pos % 3 == 0).collect();
        bitmap(1_000, 130, &marked)
    }

    #[test]
    fn covered_bits_agree_with_covered_at_at_every_offset() {
        let cached = sparse_bitmap();
        for ref_pos in 900..1_200 {
            let expected: u64 =
                (0..64).filter(|&i| cached.covered_at(ref_pos + i)).map(|i| 1u64 << i).sum();
            assert_eq!(cached.covered_bits(ref_pos), expected, "ref_pos {ref_pos}");
        }
    }

    #[test]
    fn mark_bits_agrees_with_marking_each_position() {
        let bits = 0x8000_0000_0000_0001 | 0b1011_0110 << 20;
        for ref_pos in 900..1_200 {
            let mut by_word = OverlapBitmap::new(1_000, 130);
            by_word.mark_bits(ref_pos, bits);
            let mut by_position = OverlapBitmap::new(1_000, 130);
            for i in (0..64).filter(|i| bits >> i & 1 == 1) {
                by_position.mark(ref_pos + i);
            }
            assert_eq!(by_word.bitmap, by_position.bitmap, "ref_pos {ref_pos}");
        }
    }

    #[test]
    fn mark_bits_leaves_storage_past_the_region_clear() {
        let mut cached = OverlapBitmap::new(1_000, 70);
        cached.mark_bits(1_060, u64::MAX);
        assert_eq!(cached.covered_bits(1_060), 0b11_1111_1111);
        assert_eq!(cached.bitmap.as_raw_slice()[1] >> 6, 0);
    }

    #[test]
    fn buffering_marks_words_and_pairing_reads_them_back() {
        let mut first = OverlapBitmap::new(100, 10);
        Buffering { bitmap: &mut first }.on_counted_bits(102, 0b101);
        let pairing = Pairing { cached: &first };
        assert_eq!(pairing.mate_covered_bits(100), 0b1_0100);
    }

    #[test]
    fn alone_never_reports_coverage() {
        const { assert!(<Alone as OverlapAction>::IS_ALONE) };
        assert!(!Alone.is_mate_covered(5));
    }
}
