use rust_lapper::{Interval as LapperInterval, Lapper};

use crate::intervals::{Interval, Intervals};

/// An efficient overlap detector backed by one `Lapper` per contig.
///
/// Built once from a set of intervals, then queried many times.  Internally
/// uses a `Vec<Option<Lapper>>` indexed by `ref_id` for O(1) contig lookup.
/// Within each contig the `Lapper` provides O(log n + k) overlap queries
/// using the BITS algorithm.
pub struct Overlapper<T: Clone + Send + Sync + Eq> {
    /// One optional Lapper per contig, indexed by ref_id.
    trees: Vec<Option<Lapper<u32, T>>>,
}

impl<T: Clone + Send + Sync + Eq> Overlapper<T> {
    /// Build an overlapper from an iterator of `(ref_id, start, end, value)` tuples.
    ///
    /// Coordinates are 0-based half-open `[start, end)`, matching rust-lapper's convention.
    pub fn new(intervals: impl Iterator<Item = (usize, u32, u32, T)>) -> Self {
        // Group intervals by ref_id
        let mut by_contig: Vec<Vec<LapperInterval<u32, T>>> = Vec::new();
        for (ref_id, start, stop, val) in intervals {
            if ref_id >= by_contig.len() {
                by_contig.resize_with(ref_id + 1, Vec::new);
            }
            by_contig[ref_id].push(LapperInterval { start, stop, val });
        }

        let trees: Vec<Option<Lapper<u32, T>>> = by_contig
            .into_iter()
            .map(|ivs| if ivs.is_empty() { None } else { Some(Lapper::new(ivs)) })
            .collect();

        Self { trees }
    }

    /// Build an overlapper from an `Intervals` collection, storing a clone of each
    /// `Interval` as the associated value.
    #[must_use]
    pub fn from_intervals(intervals: &Intervals) -> Overlapper<Interval> {
        Overlapper::<Interval>::new(
            intervals.iter().map(|iv| (iv.ref_id, iv.start, iv.end, iv.clone())),
        )
    }

    /// Return all values whose intervals overlap the query range `[start, end)`.
    ///
    /// Returns an empty iterator if the contig has no intervals.
    pub fn get_overlaps(&self, ref_id: usize, start: u32, end: u32) -> impl Iterator<Item = &T> {
        self.trees
            .get(ref_id)
            .and_then(Option::as_ref)
            .into_iter()
            .flat_map(move |lapper| lapper.find(start, end).map(|iv| &iv.val))
    }

    /// Return `true` if any stored interval overlaps the query range.
    pub fn overlaps_any(&self, ref_id: usize, start: u32, end: u32) -> bool {
        self.trees
            .get(ref_id)
            .and_then(Option::as_ref)
            .is_some_and(|lapper| lapper.find(start, end).next().is_some())
    }

    /// Return `true` if any stored interval contains the given 0-based position.
    #[must_use]
    pub fn contains_position(&self, ref_id: usize, pos: u32) -> bool {
        self.overlaps_any(ref_id, pos, pos + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_overlap() {
        let overlapper = Overlapper::new(
            vec![(0_usize, 10_u32, 20_u32, "A"), (0, 30, 40, "B"), (1, 5, 15, "C")].into_iter(),
        );

        // Query overlapping first interval on contig 0
        let hits: Vec<&&str> = overlapper.get_overlaps(0, 15, 25).collect();
        assert_eq!(hits, vec![&"A"]);

        // Query overlapping both intervals on contig 0
        let hits: Vec<&&str> = overlapper.get_overlaps(0, 15, 35).collect();
        assert_eq!(hits.len(), 2);

        // Query contig 1
        let hits: Vec<&&str> = overlapper.get_overlaps(1, 0, 10).collect();
        assert_eq!(hits, vec![&"C"]);

        // Query empty contig
        let hits: Vec<&&str> = overlapper.get_overlaps(2, 0, 100).collect();
        assert!(hits.is_empty());
    }

    #[test]
    fn test_overlaps_any() {
        let overlapper = Overlapper::new(vec![(0_usize, 10_u32, 20_u32, 1)].into_iter());

        assert!(overlapper.overlaps_any(0, 15, 25));
        assert!(overlapper.overlaps_any(0, 5, 15));
        assert!(!overlapper.overlaps_any(0, 20, 30)); // half-open: [10,20) does not overlap [20,30)
        assert!(!overlapper.overlaps_any(0, 0, 10)); // [0,10) does not overlap [10,20)
        assert!(!overlapper.overlaps_any(1, 10, 20)); // no intervals on contig 1
    }

    #[test]
    fn test_contains_position() {
        let overlapper = Overlapper::new(vec![(0_usize, 10_u32, 20_u32, "X")].into_iter());

        assert!(!overlapper.contains_position(0, 9));
        assert!(overlapper.contains_position(0, 10));
        assert!(overlapper.contains_position(0, 19));
        assert!(!overlapper.contains_position(0, 20)); // exclusive end
    }

    #[test]
    fn test_empty_overlapper() {
        let overlapper: Overlapper<i32> = Overlapper::new(std::iter::empty());
        assert!(!overlapper.overlaps_any(0, 0, 100));
        assert!(overlapper.get_overlaps(0, 0, 100).next().is_none());
    }

    #[test]
    fn test_sparse_contigs() {
        // Intervals only on contig 5 — contigs 0-4 should be None
        let overlapper = Overlapper::new(vec![(5_usize, 100_u32, 200_u32, "Z")].into_iter());

        assert!(!overlapper.overlaps_any(0, 0, 1000));
        assert!(!overlapper.overlaps_any(4, 0, 1000));
        assert!(overlapper.overlaps_any(5, 150, 160));
    }
}
