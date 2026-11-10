//! CIGAR-derived helpers shared by collectors that walk a record's aligned bases.

use noodles::sam::alignment::record::cigar::op::Kind;

use crate::sam::riker_record::RikerRecord;

/// One run of reference-aligned bases (an `M`, `=` or `X` CIGAR operation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignedBlock {
    /// Offset of the block's first base from the record's alignment start, in reference bases.
    pub ref_offset: u32,
    /// Offset of the block's first base into the record's stored sequence and qualities.
    pub read_offset: u32,
    /// Number of bases in the block.
    pub len: u32,
}

/// Iterate the record's reference-aligned blocks in CIGAR order. Insertions and soft clips
/// advance the read offset, deletions and skips advance the reference offset, and hard clips and
/// padding advance neither.
pub fn aligned_blocks(record: &RikerRecord) -> impl Iterator<Item = AlignedBlock> + '_ {
    let mut ref_offset: u32 = 0;
    let mut read_offset: u32 = 0;
    record.cigar_ops().filter_map(move |op| {
        // CIGAR op lengths are stored as u32 in BAM/SAM, so the cast is lossless.
        #[allow(clippy::cast_possible_truncation)]
        let len = op.len() as u32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let block = AlignedBlock { ref_offset, read_offset, len };
                ref_offset += len;
                read_offset += len;
                Some(block)
            }
            Kind::Insertion | Kind::SoftClip => {
                read_offset += len;
                None
            }
            Kind::Deletion | Kind::Skip => {
                ref_offset += len;
                None
            }
            Kind::HardClip | Kind::Pad => None,
        }
    })
}

/// Count the record's reference-aligned bases (`M`, `=` and `X` operations).
#[must_use]
pub fn count_aligned_bases(record: &RikerRecord) -> u64 {
    record
        .cigar_ops()
        .filter(|op| {
            matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
        })
        .map(|op| op.len() as u64)
        .sum()
}

/// The lengths of the record's leading and trailing hard clips, in the order the CIGAR lists
/// them (reference orientation), as `(leading, trailing)`.
#[must_use]
pub fn hard_clip_lengths(record: &RikerRecord) -> (u32, u32) {
    let mut ops = record.cigar_ops();
    let leading = match ops.next() {
        Some(op) if op.kind() == Kind::HardClip => op.len(),
        _ => 0,
    };
    let trailing = match record.cigar_ops().last() {
        Some(op) if op.kind() == Kind::HardClip && record.cigar_len() > 1 => op.len(),
        _ => 0,
    };
    // CIGAR op lengths are stored as u32 in BAM/SAM, so the casts are lossless.
    #[allow(clippy::cast_possible_truncation)]
    (leading as u32, trailing as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::read;

    fn blocks(cigar: &str) -> Vec<AlignedBlock> {
        let record = read().cigar(cigar).into_riker_record();
        aligned_blocks(&record).collect()
    }

    #[test]
    fn aligned_blocks_of_a_simple_match_is_one_block() {
        assert_eq!(blocks("10M"), vec![AlignedBlock { ref_offset: 0, read_offset: 0, len: 10 }]);
    }

    #[test]
    fn aligned_blocks_skip_soft_clips_on_the_read_only() {
        assert_eq!(
            blocks("3S10M2S"),
            vec![AlignedBlock { ref_offset: 0, read_offset: 3, len: 10 }]
        );
    }

    #[test]
    fn aligned_blocks_advance_reference_over_deletions_and_read_over_insertions() {
        assert_eq!(
            blocks("5M2D5M3I5M"),
            vec![
                AlignedBlock { ref_offset: 0, read_offset: 0, len: 5 },
                AlignedBlock { ref_offset: 7, read_offset: 5, len: 5 },
                AlignedBlock { ref_offset: 12, read_offset: 13, len: 5 },
            ]
        );
    }

    #[test]
    fn aligned_blocks_ignore_hard_clips() {
        assert_eq!(blocks("8H10M"), vec![AlignedBlock { ref_offset: 0, read_offset: 0, len: 10 }]);
    }

    #[test]
    fn count_aligned_bases_sums_match_operations_only() {
        let record = read().cigar("2S5M2D5M3I5M1S").into_riker_record();
        assert_eq!(count_aligned_bases(&record), 15);
    }

    #[test]
    fn hard_clip_lengths_reports_leading_and_trailing_clips() {
        let record = read().cigar("8H10M3H").into_riker_record();
        assert_eq!(hard_clip_lengths(&record), (8, 3));
    }

    #[test]
    fn hard_clip_lengths_is_zero_without_hard_clips() {
        let record = read().cigar("2S10M").into_riker_record();
        assert_eq!(hard_clip_lengths(&record), (0, 0));
    }
}
