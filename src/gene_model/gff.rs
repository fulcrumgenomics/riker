//! GFF3 parser (GENCODE / Ensembl / RefSeq).
//!
//! GFF3 is 9 tab-delimited columns with **1-based inclusive** coordinates and attributes of
//! the form `key=value;key2=value2`. Features are linked by `ID`/`Parent`: a gene parents
//! transcripts (`mRNA`, `lnc_RNA`, `rRNA`, …), which parent `exon` and `CDS` records.
//!
//! Parsing is two-pass and ordering-independent: first every feature is read into maps keyed
//! by ID (with exons/CDS bucketed by their `Parent`), then any feature that parents at least
//! one exon is finalized into a transcript and linked back to its gene for the name and
//! biotype. Provider differences are absorbed by trying the known attribute keys in order, so
//! GENCODE (`gene_name`/`gene_type`), Ensembl (`Name`/`biotype`), and RefSeq
//! (`Name`/`gene_biotype`) all parse without a provider flag.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result, bail};

use super::{Biotype, Exon, ParsedTranscript};

/// A non-exon/CDS feature (gene or transcript), with the fields needed to resolve names.
struct FeatureRec {
    contig: String,
    negative_strand: bool,
    parents: Vec<String>,
    name: Option<String>,
    gene_name: Option<String>,
    gene_id: Option<String>,
    transcript_id: Option<String>,
    biotype: Option<Biotype>,
}

/// Parse GFF3 content into transcripts.
///
/// # Errors
/// Returns an error if a feature line is malformed (too few columns, unparseable coordinates).
pub(super) fn parse<R: BufRead>(reader: R) -> Result<Vec<ParsedTranscript>> {
    let mut features: HashMap<String, FeatureRec> = HashMap::new();
    let mut exons_by_parent: HashMap<String, Vec<Exon>> = HashMap::new();
    let mut cds_by_parent: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    // Maps a sequence accession (e.g. `NC_000001.11`) to its common chromosome name (e.g. `1`),
    // read from `region`/`chromosome` feature attributes. RefSeq GFFs refer to contigs by
    // accession; this lets us recover the names that reconcile against a BAM header.
    let mut seqid_names: HashMap<String, String> = HashMap::new();
    // Preserve first-seen transcript order for deterministic output.
    let mut tx_order: Vec<String> = Vec::new();

    for (i, result) in reader.lines().enumerate() {
        let raw = result.with_context(|| format!("reading GFF3 line {}", i + 1))?;
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            bail!("GFF3 line {}: expected 9 columns, found {}", i + 1, cols.len());
        }
        let kind = cols[2];
        let start: u32 =
            cols[3].parse().with_context(|| format!("GFF3 line {}: bad start", i + 1))?;
        let end: u32 = cols[4].parse().with_context(|| format!("GFF3 line {}: bad end", i + 1))?;
        if end < start {
            bail!("GFF3 line {}: feature end {end} precedes start {start}", i + 1);
        }
        let negative_strand = cols[6] == "-";
        let attrs = parse_attrs(cols[8]);

        // Record the accession→common-name mapping from chromosome region features.
        if matches!(kind, "region" | "chromosome")
            && let Some(name) = attr(&attrs, "chromosome").or_else(|| attr(&attrs, "Name"))
        {
            seqid_names.entry(cols[0].to_string()).or_insert_with(|| name.to_string());
        }

        match kind {
            "exon" => {
                for parent in parents(&attrs) {
                    exons_by_parent.entry(parent).or_default().push(Exon { start, end });
                }
            }
            "CDS" => {
                for parent in parents(&attrs) {
                    cds_by_parent.entry(parent).or_default().push((start, end));
                }
            }
            _ => {
                let Some(id) = attr(&attrs, "ID").map(str::to_string) else {
                    continue; // genes/transcripts without an ID can't be linked
                };
                tx_order.push(id.clone());
                features.insert(
                    id,
                    FeatureRec {
                        contig: cols[0].to_string(),
                        negative_strand,
                        parents: parents(&attrs),
                        name: attr(&attrs, "Name").map(str::to_string),
                        gene_name: attr(&attrs, "gene_name").map(str::to_string),
                        gene_id: attr(&attrs, "gene_id").map(str::to_string),
                        transcript_id: attr(&attrs, "transcript_id").map(str::to_string),
                        biotype: biotype_of(&attrs),
                    },
                );
            }
        }
    }

    let mut transcripts = Vec::new();
    for id in tx_order {
        let Some(exons) = exons_by_parent.get(&id) else {
            continue; // not a transcript — it parents no exons
        };
        let feature = &features[&id];
        let gene = feature.parents.first().and_then(|p| features.get(p));

        let gene_name = gene
            .and_then(gene_display_name)
            .or_else(|| gene_display_name(feature))
            .unwrap_or_else(|| id.clone());
        let biotype = gene.and_then(|g| g.biotype.clone()).or_else(|| feature.biotype.clone());

        let (cds_start, cds_end) = match cds_by_parent.get(&id) {
            None => (None, None),
            Some(cds) => (
                Some(cds.iter().map(|c| c.0).min().unwrap()),
                Some(cds.iter().map(|c| c.1).max().unwrap()),
            ),
        };

        // Derive the transcript span from the exon extent (consistent with the GTF parser and
        // robust to an mRNA feature line whose declared bounds differ from its exons).
        let start = exons.iter().map(|e| e.start).min().unwrap();
        let end = exons.iter().map(|e| e.end).max().unwrap();

        // Translate an accession (RefSeq `NC_…`) to its common contig name when known.
        let contig =
            seqid_names.get(&feature.contig).cloned().unwrap_or_else(|| feature.contig.clone());

        transcripts.push(ParsedTranscript {
            gene_name,
            tx_name: transcript_name(&id, feature),
            biotype,
            contig,
            negative_strand: feature.negative_strand,
            start,
            end,
            cds_start,
            cds_end,
            exons: exons.clone(),
        });
    }
    Ok(transcripts)
}

/// Display name for a gene-like feature, trying `gene_name`, then `Name`, then `gene_id`.
/// Returns `None` if none are present; the caller then falls back to the feature `ID`.
fn gene_display_name(feature: &FeatureRec) -> Option<String> {
    feature.gene_name.clone().or_else(|| feature.name.clone()).or_else(|| feature.gene_id.clone())
}

/// Transcript name: `transcript_id` if present, else the feature ID with a leading
/// `rna-` / `transcript:` / `transcript-` prefix stripped (matching RefSeq/Ensembl ID styles).
fn transcript_name(id: &str, feature: &FeatureRec) -> String {
    if let Some(tid) = &feature.transcript_id {
        return tid.clone();
    }
    for prefix in ["rna-", "transcript:", "transcript-"] {
        if let Some(rest) = id.strip_prefix(prefix) {
            return rest.to_string();
        }
    }
    id.to_string()
}

/// Extract a biotype from a feature's attributes, trying the provider-specific keys.
fn biotype_of(attrs: &[(&str, &str)]) -> Option<Biotype> {
    attr(attrs, "gene_type")
        .or_else(|| attr(attrs, "gene_biotype"))
        .or_else(|| attr(attrs, "biotype"))
        .map(Biotype::from_attr)
}

/// Parse a GFF3 attribute column into borrowed `(key, value)` pairs (`key=value;...`).
fn parse_attrs(s: &str) -> Vec<(&str, &str)> {
    s.split(';')
        .filter_map(|field| field.trim().split_once('='))
        .map(|(k, v)| (k.trim(), v.trim()))
        .collect()
}

/// Look up the first value for `key`.
fn attr<'a>(attrs: &[(&'a str, &'a str)], key: &str) -> Option<&'a str> {
    attrs.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

/// All parent IDs (the `Parent` attribute may be comma-separated for shared features).
fn parents(attrs: &[(&str, &str)]) -> Vec<String> {
    attr(attrs, "Parent").map(|p| p.split(',').map(str::to_string).collect()).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENCODE: &str = "\
##gff-version 3
chr1\tHAVANA\tgene\t100\t400\t.\t+\t.\tID=ENSG1;gene_id=ENSG1;gene_name=MYC;gene_type=protein_coding
chr1\tHAVANA\tmRNA\t100\t400\t.\t+\t.\tID=ENST1;Parent=ENSG1;transcript_id=ENST1
chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tID=e1;Parent=ENST1
chr1\tHAVANA\tCDS\t150\t180\t.\t+\t0\tID=c1;Parent=ENST1
chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tID=e2;Parent=ENST1
";

    #[test]
    fn assembles_gencode_transcript() {
        let txs = parse(GENCODE.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(tx.gene_name, "MYC");
        assert_eq!(tx.tx_name, "ENST1");
        assert_eq!(tx.biotype, Some(Biotype::ProteinCoding));
        assert_eq!(tx.start, 100);
        assert_eq!(tx.end, 400);
        assert_eq!(tx.cds_start, Some(150));
        assert_eq!(tx.cds_end, Some(180));
        assert_eq!(tx.exons.len(), 2);
    }

    #[test]
    fn ensembl_uses_name_and_biotype() {
        let content = "\
1\tensembl\tgene\t10\t90\t.\t-\t.\tID=gene:ENSG2;Name=RNA5S;biotype=rRNA
1\tensembl\tlnc_RNA\t10\t90\t.\t-\t.\tID=transcript:ENST2;Parent=gene:ENSG2
1\tensembl\texon\t10\t90\t.\t-\t.\tParent=transcript:ENST2
";
        let txs = parse(content.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].gene_name, "RNA5S");
        assert_eq!(txs[0].biotype, Some(Biotype::Rrna));
        assert!(txs[0].negative_strand);
        // No transcript_id → ID with `transcript:` prefix stripped.
        assert_eq!(txs[0].tx_name, "ENST2");
    }

    #[test]
    fn refseq_region_maps_accession_to_common_contig_name() {
        // RefSeq refers to contigs by accession; the region line carries the common name,
        // which we use so the transcript's contig reconciles against a chr-named BAM.
        let content = "\
NC_000001.11\tRefSeq\tregion\t1\t248956422\t.\t+\t.\tID=NC_000001.11:1..248956422;Name=1;chromosome=1
NC_000001.11\tBestRefSeq\tgene\t1\t50\t.\t+\t.\tID=gene-FOO;Name=FOO;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tmRNA\t1\t50\t.\t+\t.\tID=rna-NM_1;Parent=gene-FOO;transcript_id=NM_1
NC_000001.11\tBestRefSeq\texon\t1\t50\t.\t+\t.\tParent=rna-NM_1
";
        let txs = parse(content.as_bytes()).unwrap();
        assert_eq!(txs.len(), 1);
        // Accession NC_000001.11 → common name "1" (which later reconciles to chr1).
        assert_eq!(txs[0].contig, "1");
        assert_eq!(txs[0].gene_name, "FOO");
        assert_eq!(txs[0].tx_name, "NM_1");
    }

    #[test]
    fn refseq_strips_rna_prefix_when_no_transcript_id() {
        let content = "\
NC_1\tBestRefSeq\tgene\t1\t50\t.\t+\t.\tID=gene-FOO;Name=FOO;gene_biotype=protein_coding
NC_1\tBestRefSeq\tmRNA\t1\t50\t.\t+\t.\tID=rna-NM_999;Parent=gene-FOO
NC_1\tBestRefSeq\texon\t1\t50\t.\t+\t.\tParent=rna-NM_999
";
        let txs = parse(content.as_bytes()).unwrap();
        assert_eq!(txs[0].gene_name, "FOO");
        assert_eq!(txs[0].tx_name, "NM_999");
        assert_eq!(txs[0].contig, "NC_1");
    }

    #[test]
    fn feature_with_no_exons_is_not_a_transcript() {
        let content = "\
chr1\tx\tgene\t1\t50\t.\t+\t.\tID=g1;Name=G
chr1\tx\tmRNA\t1\t50\t.\t+\t.\tID=t1;Parent=g1
";
        assert!(parse(content.as_bytes()).unwrap().is_empty());
    }

    #[test]
    fn errors_when_feature_end_precedes_start() {
        let content = "chr1\tx\texon\t90\t10\t.\t+\t.\tID=e1;Parent=t1\n";
        assert!(parse(content.as_bytes()).is_err());
    }
}
