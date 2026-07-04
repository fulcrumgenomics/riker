//! CIGAR-string parsing for the test record builders.
//!
//! Lets a test write `"5S90M5S"` instead of hand-assembling
//! `[Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 90), Op::new(Kind::SoftClip, 5)]`
//! — the exact boilerplate the review found copy-pasted across the alignment,
//! hybcap, error, and rna tests.
//!
//! Parsing is deliberately strict and *panics* on malformed input: a bad CIGAR in
//! a test is a test bug we want surfaced at the call site, not silently coerced
//! (the old `rna` helper mapped every unknown op to `M`, hiding mistakes).

use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record_buf::Cigar;

/// A parsed CIGAR: the noodles [`Cigar`] plus the query- and reference-consumed
/// lengths.
///
/// The lengths are what the record builders need to size sequence/quality vectors
/// (`query_len`) and to derive mate spans, template length, and NM (`ref_len`),
/// so we compute them once here rather than re-walking the ops at every use.
#[derive(Debug, Clone)]
pub struct ParsedCigar {
    /// The noodles record-buf CIGAR, ready to set on a record.
    pub cigar: Cigar,
    /// The operations as `(kind, length)` pairs, for walking the CIGAR against a
    /// reference when deriving a read's bases.
    pub(crate) ops: Vec<(Kind, usize)>,
    /// Number of read bases the CIGAR consumes (`M`/`I`/`S`/`=`/`X`).
    pub query_len: usize,
    /// Number of reference bases the CIGAR spans (`M`/`D`/`N`/`=`/`X`).
    pub ref_len: usize,
}

/// Parse a CIGAR string such as `"50M2000N50M"` or `"10S80M10S"` into a
/// [`ParsedCigar`].
///
/// Accepts the nine SAM operators `MIDNSHP=X`. An empty string or `"*"` yields an
/// empty CIGAR with zero lengths — the convention for an unmapped read.
///
/// # Panics
/// Panics on malformed input: an operator with no preceding length, a trailing
/// length with no operator, a zero-length operation, or an unknown operator char.
#[must_use]
pub fn parse(spec: &str) -> ParsedCigar {
    if spec.is_empty() || spec == "*" {
        return ParsedCigar {
            cigar: Vec::<Op>::new().into_iter().collect(),
            ops: Vec::new(),
            query_len: 0,
            ref_len: 0,
        };
    }

    let mut ops: Vec<(Kind, usize)> = Vec::new();
    let mut query_len = 0usize;
    let mut ref_len = 0usize;
    let mut digits = String::new();

    for ch in spec.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
            continue;
        }

        // `digits` empty here means an operator with no length (e.g. leading "M");
        // `parse` fails on the empty string, which we turn into a clear panic.
        let n: usize = digits.parse().unwrap_or_else(|_| {
            panic!("malformed CIGAR {spec:?}: operator '{ch}' has no preceding length")
        });
        assert!(n > 0, "malformed CIGAR {spec:?}: zero-length '{ch}' operation");
        digits.clear();

        let kind = kind_of(ch, spec);
        if kind.consumes_read() {
            query_len += n;
        }
        if kind.consumes_reference() {
            ref_len += n;
        }
        ops.push((kind, n));
    }

    assert!(
        digits.is_empty(),
        "malformed CIGAR {spec:?}: trailing length {digits:?} with no operator"
    );

    let cigar = ops.iter().map(|&(kind, len)| Op::new(kind, len)).collect();
    ParsedCigar { cigar, ops, query_len, ref_len }
}

/// Map a CIGAR operator character to its noodles [`Kind`].
fn kind_of(ch: char, spec: &str) -> Kind {
    match ch {
        'M' => Kind::Match,
        'I' => Kind::Insertion,
        'D' => Kind::Deletion,
        'N' => Kind::Skip,
        'S' => Kind::SoftClip,
        'H' => Kind::HardClip,
        'P' => Kind::Pad,
        '=' => Kind::SequenceMatch,
        'X' => Kind::SequenceMismatch,
        _ => panic!("malformed CIGAR {spec:?}: unknown operator '{ch}'"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an expected CIGAR from `(kind, len)` pairs for equality assertions.
    fn cig(pairs: &[(Kind, usize)]) -> Cigar {
        pairs.iter().map(|&(kind, len)| Op::new(kind, len)).collect()
    }

    #[test]
    fn parses_a_simple_match() {
        let p = parse("100M");
        assert_eq!(p.cigar, cig(&[(Kind::Match, 100)]));
        assert_eq!(p.query_len, 100);
        assert_eq!(p.ref_len, 100);
    }

    #[test]
    fn soft_clips_consume_query_but_not_reference() {
        let p = parse("10S80M10S");
        assert_eq!(p.query_len, 100);
        assert_eq!(p.ref_len, 80);
    }

    #[test]
    fn skip_consumes_reference_but_not_query() {
        // A spliced RNA read: two 50-base exons across a 2000-base intron.
        let p = parse("50M2000N50M");
        assert_eq!(p.query_len, 100);
        assert_eq!(p.ref_len, 2100);
    }

    #[test]
    fn deletion_consumes_reference_but_not_query() {
        let p = parse("50M5D45M");
        assert_eq!(p.query_len, 95);
        assert_eq!(p.ref_len, 100);
    }

    #[test]
    fn insertion_consumes_query_but_not_reference() {
        let p = parse("50M5I45M");
        assert_eq!(p.query_len, 100);
        assert_eq!(p.ref_len, 95);
    }

    #[test]
    fn hard_clip_consumes_neither() {
        let p = parse("5H100M5H");
        assert_eq!(p.query_len, 100);
        assert_eq!(p.ref_len, 100);
    }

    #[test]
    fn parses_multi_digit_lengths() {
        assert_eq!(parse("123M").cigar, cig(&[(Kind::Match, 123)]));
    }

    #[test]
    fn parses_all_nine_operators() {
        let p = parse("2H3S4M1I1D2N1P2=3X");
        // query: S3 + M4 + I1 + =2 + X3 = 13; ref: M4 + D1 + N2 + =2 + X3 = 12.
        assert_eq!(p.query_len, 13);
        assert_eq!(p.ref_len, 12);
        assert_eq!(
            p.cigar,
            cig(&[
                (Kind::HardClip, 2),
                (Kind::SoftClip, 3),
                (Kind::Match, 4),
                (Kind::Insertion, 1),
                (Kind::Deletion, 1),
                (Kind::Skip, 2),
                (Kind::Pad, 1),
                (Kind::SequenceMatch, 2),
                (Kind::SequenceMismatch, 3),
            ])
        );
    }

    #[test]
    fn empty_string_yields_empty_cigar() {
        let p = parse("");
        assert_eq!(p.cigar, cig(&[]));
        assert_eq!((p.query_len, p.ref_len), (0, 0));
    }

    #[test]
    fn star_yields_empty_cigar() {
        let p = parse("*");
        assert_eq!(p.cigar, cig(&[]));
        assert_eq!((p.query_len, p.ref_len), (0, 0));
    }

    #[test]
    #[should_panic(expected = "no preceding length")]
    fn panics_on_operator_without_length() {
        let _ = parse("M");
    }

    #[test]
    #[should_panic(expected = "trailing length")]
    fn panics_on_length_without_operator() {
        let _ = parse("100");
    }

    #[test]
    #[should_panic(expected = "unknown operator")]
    fn panics_on_unknown_operator() {
        let _ = parse("100Z");
    }

    #[test]
    #[should_panic(expected = "zero-length")]
    fn panics_on_zero_length_op() {
        let _ = parse("0M");
    }
}
