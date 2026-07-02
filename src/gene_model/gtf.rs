//! GTF parser (GENCODE / Ensembl / RefSeq).
//!
//! GTF is 9 tab-delimited columns with **1-based inclusive** coordinates and attributes of
//! the form `key "value"; key2 "value2";`. We assemble transcripts from their `exon` and
//! `CDS` feature lines, grouped by `transcript_id`, and group transcripts into genes by
//! `gene_id`. Provider differences are handled by trying the known attribute keys in order
//! (`gene_name`→`gene_id` for the name; `gene_type`/`gene_biotype`/`biotype` for the biotype).
//! Transcript span is derived from the exon extent; the coding region from the CDS extent.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result, bail};

use super::{Biotype, Exon, ParsedTranscript};

/// Accumulates the feature lines belonging to one transcript before it is finalized.
struct TxAccum {
    gene_name: String,
    biotype: Option<Biotype>,
    contig: String,
    negative_strand: bool,
    exons: Vec<Exon>,
    cds: Vec<(u32, u32)>,
}

/// Parse GTF content into transcripts.
///
/// # Errors
/// Returns an error if a feature line is malformed (too few columns, unparseable coordinates).
pub(super) fn parse<R: BufRead>(reader: R) -> Result<Vec<ParsedTranscript>> {
    // Preserve first-seen transcript order for deterministic output.
    let mut order: Vec<String> = Vec::new();
    let mut by_tx: HashMap<String, TxAccum> = HashMap::new();

    for (i, result) in reader.lines().enumerate() {
        let raw = result.with_context(|| format!("reading GTF line {}", i + 1))?;
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            bail!("GTF line {}: expected 9 columns, found {}", i + 1, cols.len());
        }
        let feature = cols[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }

        let attrs = parse_attrs(cols[8]);
        let Some(transcript_id) = attr(&attrs, "transcript_id") else {
            continue; // cannot group without a transcript id
        };
        let start: u32 =
            cols[3].parse().with_context(|| format!("GTF line {}: bad start", i + 1))?;
        let end: u32 = cols[4].parse().with_context(|| format!("GTF line {}: bad end", i + 1))?;
        if end < start {
            bail!("GTF line {}: feature end {end} precedes start {start}", i + 1);
        }
        let negative_strand = cols[6] == "-";

        let accum = by_tx.entry(transcript_id.to_string()).or_insert_with(|| {
            order.push(transcript_id.to_string());
            let gene_name = attr(&attrs, "gene_name")
                .or_else(|| attr(&attrs, "gene_id"))
                .unwrap_or(transcript_id)
                .to_string();
            let biotype = attr(&attrs, "gene_type")
                .or_else(|| attr(&attrs, "gene_biotype"))
                .or_else(|| attr(&attrs, "biotype"))
                .map(Biotype::from_attr);
            TxAccum {
                gene_name,
                biotype,
                contig: cols[0].to_string(),
                negative_strand,
                exons: Vec::new(),
                cds: Vec::new(),
            }
        });

        if feature == "exon" {
            accum.exons.push(Exon { start, end });
        } else {
            accum.cds.push((start, end));
        }
    }

    let mut transcripts = Vec::with_capacity(order.len());
    for tx_id in order {
        let accum = by_tx.remove(&tx_id).expect("tx id present");
        if accum.exons.is_empty() {
            continue; // CDS-only records can't define a transcript
        }
        let start = accum.exons.iter().map(|e| e.start).min().unwrap();
        let end = accum.exons.iter().map(|e| e.end).max().unwrap();
        let (cds_start, cds_end) = if accum.cds.is_empty() {
            (None, None)
        } else {
            (
                Some(accum.cds.iter().map(|c| c.0).min().unwrap()),
                Some(accum.cds.iter().map(|c| c.1).max().unwrap()),
            )
        };
        transcripts.push(ParsedTranscript {
            gene_name: accum.gene_name,
            tx_name: tx_id,
            biotype: accum.biotype,
            contig: accum.contig,
            negative_strand: accum.negative_strand,
            start,
            end,
            cds_start,
            cds_end,
            exons: accum.exons,
        });
    }
    Ok(transcripts)
}

/// Parse a GTF attribute column into borrowed `(key, value)` pairs. Each attribute is
/// `key "value"` (GENCODE/Ensembl/RefSeq); the surrounding quotes are stripped.
fn parse_attrs(s: &str) -> Vec<(&str, &str)> {
    let mut out = Vec::new();
    for field in s.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        if let Some((key, value)) = field.split_once(char::is_whitespace) {
            out.push((key.trim(), value.trim().trim_matches('"')));
        }
    }
    out
}

/// Look up the first value for `key` in parsed attributes.
fn attr<'a>(attrs: &[(&'a str, &'a str)], key: &str) -> Option<&'a str> {
    attrs.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENCODE: &str = "\
chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"MYC\"; gene_type \"protein_coding\";
chr1\tHAVANA\tCDS\t150\t180\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"MYC\"; gene_type \"protein_coding\";
chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"MYC\"; gene_type \"protein_coding\";
";

    #[test]
    fn assembles_transcript_from_exon_and_cds_lines() {
        let txs = parse(GENCODE.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(tx.gene_name, "MYC");
        assert_eq!(tx.tx_name, "T1");
        assert_eq!(tx.biotype, Some(Biotype::ProteinCoding));
        assert_eq!(tx.start, 100);
        assert_eq!(tx.end, 400);
        assert_eq!(tx.cds_start, Some(150));
        assert_eq!(tx.cds_end, Some(180));
        assert_eq!(tx.exons.len(), 2);
    }

    #[test]
    fn ensembl_uses_gene_biotype_and_unprefixed_contig() {
        let line = "1\tensembl\texon\t10\t90\t.\t-\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\"; gene_biotype \"rRNA\";\n";
        let txs = parse(line.as_bytes()).unwrap();
        assert_eq!(txs[0].contig, "1");
        assert!(txs[0].negative_strand);
        assert_eq!(txs[0].biotype, Some(Biotype::Rrna));
        // No gene_name → falls back to gene_id.
        assert_eq!(txs[0].gene_name, "ENSG1");
    }

    #[test]
    fn non_coding_transcript_has_no_cds() {
        let line = "chr1\tx\texon\t10\t90\t.\t+\t.\ttranscript_id \"T\"; gene_id \"G\";\n";
        let txs = parse(line.as_bytes()).unwrap();
        assert_eq!(txs[0].cds_start, None);
        assert_eq!(txs[0].cds_end, None);
    }

    #[test]
    fn errors_when_feature_end_precedes_start() {
        let line = "chr1\tx\texon\t90\t10\t.\t+\t.\ttranscript_id \"T\"; gene_id \"G\";\n";
        assert!(parse(line.as_bytes()).is_err());
    }
}
