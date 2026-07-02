//! UCSC refFlat parser.
//!
//! refFlat is tab-delimited with 11 fixed columns and **0-based half-open** coordinates:
//!
//! ```text
//! geneName  name  chrom  strand  txStart  txEnd  cdsStart  cdsEnd  exonCount  exonStarts  exonEnds
//! ```
//!
//! Coordinates are converted to 1-based inclusive exactly as fgbio `RefFlatSource` does:
//! `start + 1` for starts, `end` unchanged for ends (a half-open end is already the
//! inclusive 1-based end). A transcript with `cdsStart == cdsEnd` is non-coding. An optional
//! header line (`geneName...` or Picard's `GENE_NAME...`) is skipped if present. refFlat carries
//! no biotype column, so each transcript's biotype is inferred from CDS presence: `protein_coding`
//! when a non-zero-length CDS is present, otherwise `noncoding`.

use std::io::BufRead;

use anyhow::{Context, Result, bail};

use super::{Biotype, Exon, ParsedTranscript};

/// Parse refFlat content (streamed line by line) into transcripts.
///
/// # Errors
/// Returns an error if a line has too few columns or contains malformed numbers, or if the
/// exon start/end counts disagree with the declared exon count.
pub(super) fn parse<R: BufRead>(reader: R) -> Result<Vec<ParsedTranscript>> {
    let mut transcripts = Vec::new();
    for (i, result) in reader.lines().enumerate() {
        let raw = result.with_context(|| format!("reading refFlat line {}", i + 1))?;
        let line = raw.trim_end();
        if line.is_empty() {
            continue;
        }
        // Skip a header line if present (lowercase or Picard-style).
        if i == 0 && (line.starts_with("geneName") || line.starts_with("GENE_NAME")) {
            continue;
        }
        transcripts.push(parse_line(line, i + 1)?);
    }
    Ok(transcripts)
}

/// Parse one refFlat data line (1-based line number for error messages).
fn parse_line(line: &str, line_num: usize) -> Result<ParsedTranscript> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 11 {
        bail!("refFlat line {line_num}: expected 11 columns, found {}", cols.len());
    }

    let gene_name = cols[0].to_string();
    let tx_name = cols[1].to_string();
    let contig = cols[2].to_string();
    let negative_strand = cols[3] == "-";
    let tx_start: u32 = parse_u32(cols[4], line_num, "txStart")?;
    let tx_end: u32 = parse_u32(cols[5], line_num, "txEnd")?;
    let cds_start_raw: u32 = parse_u32(cols[6], line_num, "cdsStart")?;
    let cds_end_raw: u32 = parse_u32(cols[7], line_num, "cdsEnd")?;
    let exon_count: usize = parse_u32(cols[8], line_num, "exonCount")? as usize;
    let exon_starts = parse_csv_u32(cols[9], line_num, "exonStarts")?;
    let exon_ends = parse_csv_u32(cols[10], line_num, "exonEnds")?;

    if exon_starts.len() != exon_count || exon_ends.len() != exon_count {
        bail!(
            "refFlat line {line_num}: exonCount={exon_count} but found {} exon starts and {} exon ends",
            exon_starts.len(),
            exon_ends.len()
        );
    }

    // 0-based half-open → 1-based inclusive: starts gain 1, ends are unchanged. A non-empty exon
    // needs `e > s`; reject otherwise so the 1-based span can't invert (`start > end`).
    let mut exons = Vec::with_capacity(exon_count);
    for (&s, &e) in exon_starts.iter().zip(&exon_ends) {
        if e <= s {
            bail!("refFlat line {line_num}: exon end {e} is not after start {s}");
        }
        exons.push(Exon { start: s + 1, end: e });
    }

    // cdsStart == cdsEnd marks a non-coding transcript.
    let (cds_start, cds_end) = if cds_start_raw == cds_end_raw {
        (None, None)
    } else {
        (Some(cds_start_raw + 1), Some(cds_end_raw))
    };

    // refFlat carries no biotype column, but the presence of a non-zero-length CDS distinguishes
    // protein-coding from everything else; that coding/non-coding split is all that can be inferred.
    let biotype = if cds_start.is_some() {
        Biotype::ProteinCoding
    } else {
        Biotype::Other("noncoding".to_string())
    };

    Ok(ParsedTranscript {
        gene_name,
        tx_name,
        biotype: Some(biotype),
        contig,
        negative_strand,
        start: tx_start + 1,
        end: tx_end,
        cds_start,
        cds_end,
        exons,
    })
}

/// Parse a single unsigned integer field with a contextual error message.
fn parse_u32(s: &str, line_num: usize, field: &str) -> Result<u32> {
    s.trim()
        .parse::<u32>()
        .with_context(|| format!("refFlat line {line_num}: bad {field} value '{s}'"))
}

/// Parse a comma-separated list of unsigned integers (UCSC uses a trailing comma).
fn parse_csv_u32(s: &str, line_num: usize, field: &str) -> Result<Vec<u32>> {
    s.split(',').filter(|t| !t.is_empty()).map(|t| parse_u32(t, line_num, field)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_coding_transcript_with_two_exons() {
        // 0-based: tx [100,500), cds [150,450), exons [100,200) and [300,500).
        let line = "MYGENE\tNM_1\tchr1\t+\t100\t500\t150\t450\t2\t100,300,\t200,500,";
        let txs = parse(line.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(tx.gene_name, "MYGENE");
        assert_eq!(tx.tx_name, "NM_1");
        assert_eq!(tx.contig, "chr1");
        assert!(!tx.negative_strand);
        // 1-based inclusive conversions.
        assert_eq!(tx.start, 101);
        assert_eq!(tx.end, 500);
        assert_eq!(tx.cds_start, Some(151));
        assert_eq!(tx.cds_end, Some(450));
        assert_eq!(tx.exons, vec![Exon { start: 101, end: 200 }, Exon { start: 301, end: 500 }]);
    }

    #[test]
    fn non_coding_when_cds_start_equals_cds_end() {
        let line = "GENE\tNR_1\tchr1\t-\t100\t500\t500\t500\t1\t100,\t500,";
        let tx = &parse(line.as_bytes()).unwrap()[0];
        assert!(tx.negative_strand);
        assert_eq!(tx.cds_start, None);
        assert_eq!(tx.cds_end, None);
    }

    #[test]
    fn skips_lowercase_header() {
        let content = "geneName\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\nG\tT\tchr1\t+\t10\t90\t20\t80\t1\t10,\t90,";
        let txs = parse(content.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].gene_name, "G");
    }

    #[test]
    fn errors_on_exon_count_mismatch() {
        let line = "G\tT\tchr1\t+\t10\t90\t20\t80\t3\t10,\t90,";
        assert!(parse(line.as_bytes()).is_err());
    }

    #[test]
    fn errors_on_empty_or_inverted_exon() {
        // exonStart == exonEnd (empty, 0-based half-open) would invert to start > end.
        let line = "G\tT\tchr1\t+\t10\t90\t20\t80\t1\t50,\t50,";
        assert!(parse(line.as_bytes()).is_err());
    }
}
