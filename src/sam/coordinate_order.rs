//! Enforcement of coordinate order on a stream of mapped records.

use anyhow::{Result, bail};

use crate::sequence_dict::SequenceDictionary;

/// Rejects a mapped record that arrives before the previous one in `(reference, position)`
/// order. Collectors that tally a contig once its records have passed, or that expect a read's
/// overlapping mate to arrive after it, silently miscount on unsorted input; this catches that,
/// including a BAM whose header claims coordinate sort but isn't.
pub struct CoordinateOrderGuard {
    /// Command name used in the error message.
    tool: &'static str,
    /// Reference id and 0-based start of the previous record checked.
    last: Option<(usize, u32)>,
}

impl CoordinateOrderGuard {
    /// A guard for `tool` that has seen no records yet.
    #[must_use]
    pub fn new(tool: &'static str) -> Self {
        Self { tool, last: None }
    }

    /// Check a record at reference `ref_id` and 0-based `start` against the previous one.
    /// Returns whether the record is on a different reference from the previous one (always
    /// `true` for the first record). `dict` supplies contig names for the error message.
    ///
    /// # Errors
    /// Returns an error if the record precedes the previous one.
    pub fn check(&mut self, ref_id: usize, start: u32, dict: &SequenceDictionary) -> Result<bool> {
        let new_contig = match self.last {
            Some((last_ref_id, last_start)) => {
                if ref_id < last_ref_id || (ref_id == last_ref_id && start < last_start) {
                    let name = |id: usize| dict.get_by_index(id).map_or("?", |meta| meta.name());
                    bail!(
                        "{} requires coordinate-sorted input (ordered by reference, then \
                         position); encountered {}:{} after {}:{}. Sort with `samtools sort`.",
                        self.tool,
                        name(ref_id),
                        u64::from(start) + 1,
                        name(last_ref_id),
                        u64::from(last_start) + 1,
                    );
                }
                ref_id != last_ref_id
            }
            None => true,
        };
        self.last = Some((ref_id, start));
        Ok(new_contig)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dict() -> SequenceDictionary {
        let header =
            crate::test_support::coord_builder(&[("chr1", 1000), ("chr2", 1000)]).header().clone();
        SequenceDictionary::from(&header)
    }

    #[test]
    fn first_record_starts_a_new_contig() {
        let mut guard = CoordinateOrderGuard::new("test");
        assert!(guard.check(0, 10, &dict()).unwrap());
    }

    #[test]
    fn same_contig_at_same_or_later_position_is_accepted() {
        let dict = dict();
        let mut guard = CoordinateOrderGuard::new("test");
        guard.check(0, 10, &dict).unwrap();
        assert!(!guard.check(0, 10, &dict).unwrap());
        assert!(!guard.check(0, 11, &dict).unwrap());
    }

    #[test]
    fn moving_to_a_later_contig_reports_a_new_contig() {
        let dict = dict();
        let mut guard = CoordinateOrderGuard::new("test");
        guard.check(0, 500, &dict).unwrap();
        assert!(guard.check(1, 0, &dict).unwrap());
    }

    #[test]
    fn earlier_position_on_the_same_contig_is_rejected() {
        let dict = dict();
        let mut guard = CoordinateOrderGuard::new("test");
        guard.check(0, 10, &dict).unwrap();
        let err = guard.check(0, 9, &dict).unwrap_err().to_string();
        assert!(err.contains("test requires coordinate-sorted input"), "{err}");
        assert!(err.contains("chr1:10 after chr1:11"), "{err}");
    }

    #[test]
    fn earlier_contig_is_rejected() {
        let dict = dict();
        let mut guard = CoordinateOrderGuard::new("test");
        guard.check(1, 10, &dict).unwrap();
        let err = guard.check(0, 50, &dict).unwrap_err().to_string();
        assert!(err.contains("chr1:51 after chr2:11"), "{err}");
    }
}
