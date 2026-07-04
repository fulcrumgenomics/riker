use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;

/// Count the number of aligned bases at the 3' end of the read that overlap
/// with the mate's alignment, suitable for soft-clipping to avoid double-counting
/// coverage from overlapping read pairs.
///
/// Implements the same algorithm as Picard/htsjdk's
/// `SAMUtils.getNumOverlappingAlignedBasesToClip()`:
///
/// - Only the left-most read of a pair is clipped (lower alignment start).
/// - On ties, the second-of-pair is clipped (first-of-pair is kept).
/// - Walks the read's own CIGAR from its alignment start, counting read bases
///   that fall at or after the mate's alignment start (PNEXT).
/// - Only uses POS, PNEXT, FLAG bits, and the current read's CIGAR.
///   Does NOT require the MC (mate CIGAR) tag.
///
/// Returns 0 if the read should not be clipped (unpaired, unmapped mate,
/// read is not the left-most, etc.).
#[must_use]
pub fn count_overlapping_bases(record: &RecordBuf) -> u64 {
    let flags = record.flags();

    // Must be paired, both ends mapped
    if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
        return 0;
    }

    // Both reads must be on the same contig
    let ref_id = record.reference_sequence_id();
    let mate_ref_id = record.mate_reference_sequence_id();
    if ref_id != mate_ref_id {
        return 0;
    }

    let alignment_start = match record.alignment_start() {
        Some(pos) => pos.get(),
        None => return 0,
    };

    let mate_alignment_start = match record.mate_alignment_start() {
        Some(pos) => pos.get(),
        None => return 0,
    };

    // Only clip the left-most read (the one whose mate starts further right).
    // If both start at the same position, clip the second-of-pair.
    if mate_alignment_start < alignment_start {
        return 0;
    }
    if mate_alignment_start == alignment_start && flags.is_first_segment() {
        return 0;
    }

    // Walk the CIGAR to count read bases at or past the mate's start position.
    // mate_alignment_start is 1-based; ref_pos tracks the current 1-based ref position.
    let mut ref_pos = alignment_start;
    let mut bases_to_clip: u64 = 0;

    for op_result in record.cigar().iter() {
        let op: Op = match op_result {
            Ok(op) => op,
            Err(_) => return 0,
        };

        let kind = op.kind();
        let len = op.len();

        // Determine if this op consumes reference bases
        let consumes_ref = matches!(
            kind,
            Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip
        );
        let ref_bases_len = if consumes_ref { len } else { 0 };

        // Check if this element reaches or passes the mate start.
        // For zero-ref-length ops (I, S, H, P), ref_bases_len is 0 so
        // the condition checks if we're already past mate start.
        let element_end = ref_pos + ref_bases_len;
        let past_mate = if ref_bases_len > 0 {
            mate_alignment_start < ref_pos + ref_bases_len
        } else {
            mate_alignment_start <= ref_pos
        };

        if past_mate {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if mate_alignment_start <= ref_pos {
                        bases_to_clip += len as u64;
                    } else {
                        bases_to_clip += (element_end - mate_alignment_start) as u64;
                    }
                }
                Kind::SoftClip | Kind::HardClip | Kind::Pad | Kind::Skip | Kind::Deletion => {
                    // S/H/P/N don't count toward clip total;
                    // D consumes reference but not query — no read bases to clip
                }
                Kind::Insertion => {
                    bases_to_clip += len as u64;
                }
            }
        }

        ref_pos += ref_bases_len;
    }

    bases_to_clip
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{pair, read};

    /// A paired, both-ends-mapped record on chr1 at `alignment_start` whose mate maps at
    /// `mate_start`, with the given CIGAR. `is_first_of_pair` selects which mate the
    /// record represents; either way the returned record itself sits at `alignment_start`.
    fn make_record(
        alignment_start: usize,
        mate_start: usize,
        cigar: &str,
        is_first_of_pair: bool,
    ) -> RecordBuf {
        // `pair().at` puts read 1 (first) at the first position and read 2 (last) at the
        // second, each carrying the other's start as its mate; pick the one that lands at
        // `alignment_start`.
        if is_first_of_pair {
            let (r1, _) =
                pair("read1").at("chr1", alignment_start, mate_start).cigar(cigar).build();
            r1
        } else {
            let (_, r2) =
                pair("read1").at("chr1", mate_start, alignment_start).cigar(cigar).build();
            r2
        }
    }

    #[test]
    fn test_simple_overlap() {
        // 10M read starting at pos 1, mate starts at pos 6 → 5 bases overlap
        let rec = make_record(1, 6, "10M", true);
        assert_eq!(count_overlapping_bases(&rec), 5);
    }

    #[test]
    fn test_no_overlap() {
        // 10M read starting at pos 1, mate starts at pos 20 → no overlap
        let rec = make_record(1, 20, "10M", true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_full_overlap() {
        // 10M read starting at pos 1, mate starts at pos 1, second of pair → 10 bases
        let rec = make_record(1, 1, "10M", false);
        assert_eq!(count_overlapping_bases(&rec), 10);
    }

    #[test]
    fn test_first_of_pair_tie_not_clipped() {
        // At same position, first of pair should NOT be clipped
        let rec = make_record(1, 1, "10M", true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_right_most_read_not_clipped() {
        // Mate starts before this read → this read is right-most → not clipped
        let rec = make_record(10, 5, "10M", true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_unpaired_not_clipped() {
        let rec = read().at("chr1", 1).cigar("10M").build();
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_overlap_with_deletion() {
        // 5M2D5M at pos 1, mate starts at pos 9
        // Ref positions: M covers 1-5, D covers 6-7, M covers 8-12
        // Mate starts at 9 → overlap is positions 9-12 = 4 M bases, plus 0 D bases
        let rec = make_record(1, 9, "5M2D5M", true);
        assert_eq!(count_overlapping_bases(&rec), 4);
    }

    #[test]
    fn test_overlap_with_insertion() {
        // 5M3I5M at pos 1, mate starts at pos 4
        // Ref positions: first M covers 1-5, I has no ref positions, second M covers 6-10
        // Mate starts at 4 → from first M: 2 bases (pos 4-5), all of I: 3 bases, all of second M: 5 bases = 10
        let rec = make_record(1, 4, "5M3I5M", true);
        assert_eq!(count_overlapping_bases(&rec), 10);
    }
}
