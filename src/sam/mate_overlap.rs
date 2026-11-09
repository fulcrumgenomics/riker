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

    #[test]
    fn alone_never_reports_coverage() {
        const { assert!(<Alone as OverlapAction>::IS_ALONE) };
        assert!(!Alone.is_mate_covered(5));
    }
}
