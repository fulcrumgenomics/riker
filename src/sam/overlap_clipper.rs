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
    use noodles::core::Position;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

    /// Build a paired record with the given alignment start, mate start, cigar, and flags.
    fn make_record(
        alignment_start: usize,
        mate_start: usize,
        cigar_ops: &[Op],
        is_first_of_pair: bool,
    ) -> RecordBuf {
        let read_len: usize = cigar_ops
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                        | Kind::Insertion
                        | Kind::SoftClip
                )
            })
            .map(|op| op.len())
            .sum();
        let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
        if is_first_of_pair {
            flags |= Flags::FIRST_SEGMENT;
        } else {
            flags |= Flags::LAST_SEGMENT;
        }

        let cigar: Cigar = cigar_ops.iter().copied().collect();
        let seq: Sequence = vec![b'A'; read_len].into();
        let qual = QualityScores::from(vec![30u8; read_len]);
        let mq = MappingQuality::new(60).unwrap();

        RecordBuf::builder()
            .set_name("read1")
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(alignment_start).unwrap())
            .set_mapping_quality(mq)
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(mate_start).unwrap())
            .set_template_length(0)
            .set_sequence(seq)
            .set_quality_scores(qual)
            .build()
    }

    #[test]
    fn test_simple_overlap() {
        // 10M read starting at pos 1, mate starts at pos 6 → 5 bases overlap
        let rec = make_record(1, 6, &[Op::new(Kind::Match, 10)], true);
        assert_eq!(count_overlapping_bases(&rec), 5);
    }

    #[test]
    fn test_no_overlap() {
        // 10M read starting at pos 1, mate starts at pos 20 → no overlap
        let rec = make_record(1, 20, &[Op::new(Kind::Match, 10)], true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_full_overlap() {
        // 10M read starting at pos 1, mate starts at pos 1, second of pair → 10 bases
        let rec = make_record(1, 1, &[Op::new(Kind::Match, 10)], false);
        assert_eq!(count_overlapping_bases(&rec), 10);
    }

    #[test]
    fn test_first_of_pair_tie_not_clipped() {
        // At same position, first of pair should NOT be clipped
        let rec = make_record(1, 1, &[Op::new(Kind::Match, 10)], true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_right_most_read_not_clipped() {
        // Mate starts before this read → this read is right-most → not clipped
        let rec = make_record(10, 5, &[Op::new(Kind::Match, 10)], true);
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_unpaired_not_clipped() {
        let rec = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(1).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
            .set_sequence(vec![b'A'; 10].into())
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build();
        assert_eq!(count_overlapping_bases(&rec), 0);
    }

    #[test]
    fn test_overlap_with_deletion() {
        // 5M2D5M at pos 1, mate starts at pos 9
        // Ref positions: M covers 1-5, D covers 6-7, M covers 8-12
        // Mate starts at 9 → overlap is positions 9-12 = 4 M bases, plus 0 D bases
        let rec = make_record(
            1,
            9,
            &[Op::new(Kind::Match, 5), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 5)],
            true,
        );
        assert_eq!(count_overlapping_bases(&rec), 4);
    }

    #[test]
    fn test_overlap_with_insertion() {
        // 5M3I5M at pos 1, mate starts at pos 4
        // Ref positions: first M covers 1-5, I has no ref positions, second M covers 6-10
        // Mate starts at 4 → from first M: 2 bases (pos 4-5), all of I: 3 bases, all of second M: 5 bases = 10
        let rec = make_record(
            1,
            4,
            &[Op::new(Kind::Match, 5), Op::new(Kind::Insertion, 3), Op::new(Kind::Match, 5)],
            true,
        );
        assert_eq!(count_overlapping_bases(&rec), 10);
    }
}
