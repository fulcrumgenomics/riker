use noodles::sam::alignment::RecordBuf;

/// Pair orientation for paired-end reads.
///
/// Matches htsjdk's `SamPairUtil.PairOrientation` classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PairOrientation {
    /// Forward-Reverse ("innie"): reads face each other.
    /// Positive-strand 5' end < negative-strand 5' end.
    Fr,

    /// Reverse-Forward ("outie"): reads face away from each other.
    /// Positive-strand 5' end >= negative-strand 5' end.
    Rf,

    /// Tandem: both reads on the same strand.
    Tandem,
}

/// Classify the pair orientation of a record using htsjdk's algorithm.
///
/// Requires that the record is paired, both reads are mapped, and both are on
/// the same reference. Returns `None` if those preconditions are not met.
#[must_use]
pub fn get_pair_orientation(record: &RecordBuf) -> Option<PairOrientation> {
    let flags = record.flags();
    if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
        return None;
    }
    if record.reference_sequence_id() != record.mate_reference_sequence_id() {
        return None;
    }

    let is_reverse = flags.is_reverse_complemented();
    let mate_reverse = flags.is_mate_reverse_complemented();

    if is_reverse == mate_reverse {
        return Some(PairOrientation::Tandem);
    }

    let alignment_start = record.alignment_start().map_or(0, usize::from);
    let mate_start = record.mate_alignment_start().map_or(0, usize::from);
    let insert_size = record.template_length();

    #[expect(
        clippy::cast_possible_wrap,
        reason = "genomic positions are non-negative and fit in i64"
    )]
    let (positive_five_prime, negative_five_prime) = if is_reverse {
        let ref_len = record.cigar().alignment_span();
        let end = alignment_start + ref_len.saturating_sub(1);
        (mate_start as i64, end as i64)
    } else {
        (alignment_start as i64, alignment_start as i64 + i64::from(insert_size))
    };

    if positive_five_prime < negative_five_prime {
        Some(PairOrientation::Fr)
    } else {
        Some(PairOrientation::Rf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind},
    };
    use noodles::sam::alignment::record_buf::Cigar;

    fn make_paired_record(
        pos: usize,
        mate_pos: usize,
        tlen: i32,
        is_reverse: bool,
        mate_reverse: bool,
        read_len: usize,
    ) -> RecordBuf {
        let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
        if is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if mate_reverse {
            flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }

        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();

        RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(pos).expect("valid position"))
            .set_mapping_quality(MappingQuality::new(60).expect("valid mapq"))
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(mate_pos).expect("valid position"))
            .set_template_length(tlen)
            .build()
    }

    #[test]
    fn test_fr_orientation() {
        // Read1 forward at 100, read2 reverse at 300, tlen=300
        let record = make_paired_record(100, 300, 300, false, true, 100);
        assert_eq!(get_pair_orientation(&record), Some(PairOrientation::Fr));
    }

    #[test]
    fn test_rf_orientation() {
        // Read1 (reverse) at 100, read2 (forward/mate) at 300.
        // negative_5' = 100 + 100 - 1 = 199; positive_5' = 300 → 300 >= 199 → RF.
        let record = make_paired_record(100, 300, 300, true, false, 100);
        assert_eq!(get_pair_orientation(&record), Some(PairOrientation::Rf));
    }

    #[test]
    fn test_tandem_both_forward() {
        let record = make_paired_record(100, 300, 300, false, false, 100);
        assert_eq!(get_pair_orientation(&record), Some(PairOrientation::Tandem));
    }

    #[test]
    fn test_tandem_both_reverse() {
        let record = make_paired_record(100, 300, 300, true, true, 100);
        assert_eq!(get_pair_orientation(&record), Some(PairOrientation::Tandem));
    }

    #[test]
    fn test_unpaired_returns_none() {
        let record = RecordBuf::builder().build();
        assert_eq!(get_pair_orientation(&record), None);
    }

    #[test]
    fn test_mate_different_contig() {
        // Reads on different contigs → None
        let flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED;
        let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
        let record = RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(1) // different contig
            .set_mate_alignment_start(Position::new(200).unwrap())
            .set_template_length(200)
            .build();
        assert_eq!(get_pair_orientation(&record), None);
    }

    #[test]
    fn test_mate_unmapped() {
        let flags = Flags::SEGMENTED | Flags::MATE_UNMAPPED;
        let record = RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .build();
        assert_eq!(get_pair_orientation(&record), None);
    }

    #[test]
    fn test_record_unmapped() {
        let flags = Flags::SEGMENTED | Flags::UNMAPPED;
        let record = RecordBuf::builder()
            .set_flags(flags)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(100).unwrap())
            .build();
        assert_eq!(get_pair_orientation(&record), None);
    }

    #[test]
    fn test_fr_rf_boundary() {
        // When positive_five_prime == negative_five_prime → Rf (not Fr)
        // Forward read at pos=100, tlen=0 → negative_five_prime = 100+0 = 100
        // positive_five_prime = 100, not < 100, so → Rf
        let flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED;
        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let record = RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(Position::new(100).unwrap())
            .set_template_length(0)
            .build();
        assert_eq!(get_pair_orientation(&record), Some(PairOrientation::Rf));
    }
}
