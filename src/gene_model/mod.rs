//! Gene-model parsing for the `rna` command.
//!
//! Loads transcript/exon/CDS structure from any of the popular gene-model
//! formats and exposes it as a queryable [`GeneModel`]:
//!
//! - **UCSC refFlat** — tab-delimited, fixed columns, 0-based half-open.
//! - **GTF** — GENCODE / Ensembl / RefSeq, attributes as `key "value";`, 1-based inclusive.
//! - **GFF3** — GENCODE / Ensembl / RefSeq, attributes as `key=value;`, `ID`/`Parent`
//!   hierarchy, 1-based inclusive.
//!
//! The format is auto-detected from the file contents (see [`detect_format`]); gzip/bgzip
//! inputs are decompressed transparently.
//!
//! ## Coordinate convention
//!
//! Internally every coordinate is **1-based inclusive**, matching both reference tools
//! (Picard `Gene`/`RefFlatReader` and fgbio `GeneAnnotations`) and the native GTF/GFF
//! convention. refFlat (0-based half-open) is converted on parse exactly as fgbio's
//! `RefFlatSource` does (`start + 1`, `end` unchanged). Loci are indexed per contig sorted by
//! start, and the collector walks them with a coordinate sweep over the (coordinate-sorted) reads.
//!
//! ## Gene → locus → transcript
//!
//! Transcripts are grouped by gene name, then clustered into [`GeneLocus`]es: transcripts of
//! the same gene that share a contig and strand and whose spans overlap form one locus
//! (ported from fgbio's `RefFlatSource`). A gene whose transcripts map to two chromosomes
//! therefore yields two independent loci — each locus is the unit of overlap and analysis,
//! so multi-locus "genes" behave as separate genes, as intended.

mod gff;
mod gtf;
mod refflat;

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result, bail};
use fgoxide::io::Io;

use crate::sequence_dict::{SequenceDictionary, SequenceMetadata};

// ─── GeneModel ───────────────────────────────────────────────────────────────

/// A parsed, queryable gene model: a flat list of [`GeneLocus`]es plus a per-contig index of
/// locus indices sorted by start, used to drive a coordinate sweep over sorted reads.
pub struct GeneModel {
    loci: Vec<GeneLocus>,
    /// For each reference id (BAM ref index), the indices of loci on that contig sorted by start
    /// ascending. Drives the coordinate sweep in the collector (no interval tree needed for a
    /// coordinate-sorted scan).
    loci_by_contig: Vec<Vec<u32>>,
}

impl GeneModel {
    /// Load a gene model from `path`, resolving contig names against `dict`.
    ///
    /// The format is auto-detected. Loci on contigs that cannot be resolved to the
    /// BAM header (after `chr`/`MT` reconciliation) are dropped with a warning.
    ///
    /// # Errors
    /// Returns an error if the file cannot be read or parsed, or if no usable loci remain.
    pub fn from_path(path: &Path, dict: &SequenceDictionary) -> Result<Self> {
        let io = Io::default();
        let open = || {
            io.new_reader(path)
                .with_context(|| format!("Failed to open gene model: {}", path.display()))
        };

        // fgoxide decompresses by extension, so a gzipped file without a `.gz`/`.bgz` extension is
        // read raw — format detection would then fail on binary content with a confusing error.
        // Catch that case up front with a targeted hint.
        let mut magic = [0u8; 2];
        let looks_gzipped = std::fs::File::open(path)
            .and_then(|mut f| std::io::Read::read_exact(&mut f, &mut magic).map(|()| magic))
            .is_ok_and(|m| m == [0x1f, 0x8b]);
        if looks_gzipped && !Io::is_gzip_path(path) {
            bail!(
                "Gene model {} appears to be gzip-compressed but lacks a .gz/.bgz extension; \
                 rename it or add the extension so it can be decompressed",
                path.display()
            );
        }

        // Stream line by line (compression handled transparently by extension) rather than slurping
        // the whole — often >1 GB uncompressed — model into a String. Detection reads a throwaway
        // reader and propagates read/decompression errors; parsing re-opens a fresh one.
        let format = detect_format(open()?)?
            .with_context(|| format!("Could not detect gene-model format: {}", path.display()))?;
        log::info!("rna: parsing gene model {} as {format}", path.display());

        let transcripts = match format {
            GeneModelFormat::RefFlat => refflat::parse(open()?),
            GeneModelFormat::Gtf => gtf::parse(open()?),
            GeneModelFormat::Gff3 => gff::parse(open()?),
        }
        .with_context(|| format!("Failed to parse gene model {} as {format}", path.display()))?;

        let loci = build_loci(transcripts, dict);
        if loci.is_empty() {
            bail!(
                "Gene model {} produced no loci on contigs present in the BAM header",
                path.display()
            );
        }

        // Index loci by contig, sorted by start, for the coordinate sweep. Sized to the full
        // dictionary so every BAM ref id maps to a (possibly empty) list.
        let mut loci_by_contig: Vec<Vec<u32>> = vec![Vec::new(); dict.len()];
        for (i, locus) in loci.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation)]
            loci_by_contig[locus.ref_id].push(i as u32);
        }
        for list in &mut loci_by_contig {
            list.sort_by_key(|&i| loci[i as usize].start);
        }

        log::info!("rna: loaded {} gene loci", loci.len());
        Ok(Self { loci, loci_by_contig })
    }

    /// All loci in the model.
    #[must_use]
    pub fn loci(&self) -> &[GeneLocus] {
        &self.loci
    }

    /// The locus at `index`.
    #[must_use]
    pub fn locus(&self, index: u32) -> &GeneLocus {
        &self.loci[index as usize]
    }

    /// Indices of all loci on `ref_id`, sorted by start ascending. Empty if the contig has no
    /// loci (or `ref_id` is out of range). Used to drive the coordinate sweep.
    #[must_use]
    pub fn loci_on_contig(&self, ref_id: usize) -> &[u32] {
        self.loci_by_contig.get(ref_id).map_or(&[], Vec::as_slice)
    }

    /// Classify the aligned bases of a read's alignment `blocks` (each `(start, len)`, 1-based)
    /// against the overlapping `loci`, returning per-class base counts.
    ///
    /// Each base is counted once at the strongest function across *all* overlapping loci
    /// (`Coding > Utr > Intronic > Intergenic`). This is computed by interval arithmetic — the
    /// intersection of each block with the merged coding, exon, and span interval sets — rather
    /// than per-base iteration. For the common 0- or 1-locus read the locus's precomputed unions
    /// are used directly; multi-locus reads merge the loci's unions first (so overlapping genes
    /// compose by union, matching the per-base max).
    #[must_use]
    pub fn classify_blocks(&self, loci: &[u32], blocks: &[(u32, u32)]) -> BaseClassCounts {
        match loci {
            [] => {
                let total: u64 = blocks.iter().map(|&(_, len)| u64::from(len)).sum();
                BaseClassCounts { intergenic: total, ..BaseClassCounts::default() }
            }
            [only] => {
                let locus = &self.loci[*only as usize];
                count_blocks(
                    blocks,
                    &locus.coding_union,
                    &locus.exon_union,
                    &[(locus.start, locus.end)],
                )
            }
            many => {
                // Each locus's coding/exon unions are already sorted and disjoint, so union them
                // with a merge rather than concatenate-and-sort. The spans (one interval per locus)
                // are few and unsorted, so a plain sort-merge is cheapest there.
                let coding = union_of_sorted_disjoint(
                    many.iter().map(|&i| self.loci[i as usize].coding_union.as_slice()),
                );
                let exon = union_of_sorted_disjoint(
                    many.iter().map(|&i| self.loci[i as usize].exon_union.as_slice()),
                );
                let span = merge_intervals(
                    many.iter()
                        .map(|&i| (self.loci[i as usize].start, self.loci[i as usize].end))
                        .collect(),
                );
                count_blocks(blocks, &coding, &exon, &span)
            }
        }
    }

    /// Ribosomal intervals derived from gene biotypes (rRNA / Mt_rRNA / rRNA_pseudogene),
    /// as 0-based half-open `(ref_id, start, end)` tuples spanning each ribosomal locus.
    ///
    /// Empty for refFlat models (whose inferred biotypes are only `protein_coding` / `noncoding`,
    /// never ribosomal). Combined with any explicit `--ribosomal-intervals` in the collector.
    #[must_use]
    pub fn ribosomal_intervals(&self) -> Vec<(usize, u32, u32)> {
        self.loci
            .iter()
            .filter(|locus| locus.biotype.as_ref().is_some_and(Biotype::is_ribosomal))
            .map(|locus| (locus.ref_id, locus.start.saturating_sub(1), locus.end))
            .collect()
    }
}

// ─── GeneLocus ───────────────────────────────────────────────────────────────

/// A set of transcripts of one gene that map to the same locus (contig, strand,
/// overlapping span). The unit of overlap detection and per-read analysis.
pub struct GeneLocus {
    /// Name of the gene this locus belongs to.
    pub gene_name: String,
    /// Gene biotype: the annotated value for GTF/GFF, or for refFlat the inferred
    /// `protein_coding` / `noncoding` split (from CDS presence). `None` only if no transcript
    /// carried or implied one.
    pub biotype: Option<Biotype>,
    /// 0-based reference index into the BAM header.
    pub ref_id: usize,
    /// 1-based inclusive start (minimum transcript start in the locus).
    pub start: u32,
    /// 1-based inclusive end (maximum transcript end in the locus).
    pub end: u32,
    /// True if the locus is on the negative strand.
    pub negative_strand: bool,
    /// Transcripts belonging to this locus.
    pub transcripts: Vec<Transcript>,
    /// Merged (disjoint, sorted) exon intervals across all transcripts, 1-based inclusive.
    /// Lets per-base classification run once per locus instead of once per transcript.
    exon_union: Vec<(u32, u32)>,
    /// Merged (disjoint, sorted) coding intervals (exon ∩ CDS) across all transcripts,
    /// 1-based inclusive.
    coding_union: Vec<(u32, u32)>,
}

impl GeneLocus {
    /// Length of the locus's concatenated exon union (number of distinct exonic bases). This is
    /// the length of the per-locus coverage array.
    #[must_use]
    pub fn exon_union_len(&self) -> usize {
        self.exon_union.iter().map(|&(s, e)| (e - s + 1) as usize).sum()
    }

    /// Add 1 to the locus's exon-union coverage for each base of the 1-based inclusive block
    /// `[block_start, block_end]` that falls in an exon of any transcript. `coverage` is indexed
    /// by concatenated exon-union coordinate (length [`Self::exon_union_len`]). Accumulating once
    /// per locus dedups bases shared by overlapping transcripts; each transcript's own coverage
    /// is recovered later by [`Self::transcript_coverage`].
    pub fn add_union_coverage(&self, block_start: u32, block_end: u32, coverage: &mut [u32]) {
        let mut offset = 0u32;
        for &(exon_start, exon_end) in &self.exon_union {
            if exon_start > block_end {
                break;
            }
            let lo = block_start.max(exon_start);
            let hi = block_end.min(exon_end);
            if lo <= hi {
                let base = (offset + (lo - exon_start)) as usize;
                for c in &mut coverage[base..=base + (hi - lo) as usize] {
                    *c += 1;
                }
            }
            offset += exon_end - exon_start + 1;
        }
    }

    /// Reconstruct one transcript's per-base coverage (transcript coordinates, introns collapsed)
    /// from this locus's exon-union coverage. Equivalent to having accumulated coverage against
    /// that transcript directly: every transcript exon lies within a single exon-union interval,
    /// so its bases map to a contiguous slice of `union_cov`.
    #[must_use]
    pub fn transcript_coverage(&self, tx: &Transcript, union_cov: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(tx.transcribed_length() as usize);
        for exon in &tx.exons {
            let base = self.union_offset_of(exon.start) as usize;
            let len = (exon.end - exon.start + 1) as usize;
            out.extend_from_slice(&union_cov[base..base + len]);
        }
        out
    }

    /// Exon-union coordinate of a 1-based genomic position that lies in some exon of the locus.
    fn union_offset_of(&self, pos: u32) -> u32 {
        let mut offset = 0u32;
        for &(exon_start, exon_end) in &self.exon_union {
            if pos >= exon_start && pos <= exon_end {
                return offset + (pos - exon_start);
            }
            offset += exon_end - exon_start + 1;
        }
        unreachable!("genomic position {pos} is not within the locus exon union")
    }
}

/// Per-class aligned-base counts produced by [`GeneModel::classify_blocks`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaseClassCounts {
    /// Bases in a coding exon of any overlapping transcript.
    pub coding: u64,
    /// Bases exonic in some transcript but coding in none.
    pub utr: u64,
    /// Bases inside some gene locus's span but exonic in none.
    pub intronic: u64,
    /// Bases outside every gene locus.
    pub intergenic: u64,
}

// ─── Transcript ──────────────────────────────────────────────────────────────

/// A single transcript: its span, optional coding region, and exons. Exons are stored in genomic
/// order (sorted by start) and pairwise disjoint (see `normalize_exons`). All coordinates are
/// 1-based inclusive.
pub struct Transcript {
    /// Transcript name / id.
    pub name: String,
    /// 1-based inclusive transcript start (first exon start).
    pub start: u32,
    /// 1-based inclusive transcript end (last exon end).
    pub end: u32,
    /// 1-based inclusive coding start, or `None` for non-coding transcripts.
    pub cds_start: Option<u32>,
    /// 1-based inclusive coding end, or `None` for non-coding transcripts.
    pub cds_end: Option<u32>,
    /// Exons in genomic order (sorted by start), pairwise disjoint.
    pub exons: Vec<Exon>,
}

impl Transcript {
    /// Length of the transcribed product (sum of exon lengths, introns removed).
    #[must_use]
    pub fn transcribed_length(&self) -> u32 {
        self.exons.iter().map(Exon::length).sum()
    }
}

// ─── Exon ────────────────────────────────────────────────────────────────────

/// An exon, 1-based inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    /// 1-based inclusive start.
    pub start: u32,
    /// 1-based inclusive end.
    pub end: u32,
}

impl Exon {
    /// Number of bases spanned (inclusive).
    #[must_use]
    pub fn length(&self) -> u32 {
        self.end - self.start + 1
    }
}

// ─── Biotype ─────────────────────────────────────────────────────────────────

/// Gene biotype, extracted from GTF/GFF attributes (`gene_type` / `gene_biotype` /
/// `biotype`). Only the variants riker reasons about are named; everything else is
/// preserved verbatim in [`Biotype::Other`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Biotype {
    /// `protein_coding`.
    ProteinCoding,
    /// `rRNA`.
    Rrna,
    /// `Mt_rRNA` (mitochondrial rRNA).
    MtRrna,
    /// `rRNA_pseudogene`.
    RrnaPseudogene,
    /// Any other biotype, stored verbatim.
    Other(String),
}

impl Biotype {
    /// Classify a biotype attribute string (case-insensitive for the ribosomal variants).
    #[must_use]
    pub fn from_attr(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "protein_coding" => Self::ProteinCoding,
            "rrna" => Self::Rrna,
            "mt_rrna" => Self::MtRrna,
            "rrna_pseudogene" => Self::RrnaPseudogene,
            _ => Self::Other(s.to_string()),
        }
    }

    /// True for ribosomal biotypes (rRNA, Mt_rRNA, rRNA_pseudogene).
    #[must_use]
    pub fn is_ribosomal(&self) -> bool {
        matches!(self, Self::Rrna | Self::MtRrna | Self::RrnaPseudogene)
    }
}

// ─── ParsedTranscript (parser → grouping handoff) ──────────────────────────────

/// A transcript as produced by a parser, before grouping into loci and resolving the
/// contig name to a BAM `ref_id`. Coordinates are already normalized to 1-based inclusive.
pub(crate) struct ParsedTranscript {
    /// Gene name (or gene id when no name attribute is present).
    pub gene_name: String,
    /// Transcript name / id (may be empty; [`Self::name`] falls back to the gene name).
    pub tx_name: String,
    /// Gene biotype when the source format carries one (GTF/GFF only).
    pub biotype: Option<Biotype>,
    /// Source contig name, before reconciliation against the BAM header.
    pub contig: String,
    /// True if the transcript is on the negative strand.
    pub negative_strand: bool,
    /// 1-based inclusive transcript start (minimum exon start).
    pub start: u32,
    /// 1-based inclusive transcript end (maximum exon end).
    pub end: u32,
    /// 1-based inclusive coding start, or `None` for a non-coding transcript.
    pub cds_start: Option<u32>,
    /// 1-based inclusive coding end, or `None` for a non-coding transcript.
    pub cds_end: Option<u32>,
    /// Exons in source order (sorted into genomic order when the locus is built).
    pub exons: Vec<Exon>,
}

impl ParsedTranscript {
    /// The transcript name, preferred when present, falling back to the gene name.
    fn name(&self) -> String {
        if self.tx_name.is_empty() { self.gene_name.clone() } else { self.tx_name.clone() }
    }
}

// ─── Format detection ──────────────────────────────────────────────────────────

/// The supported gene-model file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneModelFormat {
    /// UCSC refFlat (11 tab-delimited columns).
    RefFlat,
    /// GTF (`key "value";` attributes).
    Gtf,
    /// GFF3 (`key=value;` attributes, ID/Parent hierarchy).
    Gff3,
}

impl std::fmt::Display for GeneModelFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RefFlat => write!(f, "refFlat"),
            Self::Gtf => write!(f, "GTF"),
            Self::Gff3 => write!(f, "GFF3"),
        }
    }
}

/// Detect the gene-model format from file content by sniffing the first meaningful line.
///
/// - `##gff-version 3` directive ⇒ GFF3.
/// - A refFlat header (`geneName...` or Picard `GENE_NAME...`) ⇒ refFlat.
/// - A 9-column feature line ⇒ GTF if the attributes use `key "value"`, GFF3 if `key=value`.
/// - A line with ≥10 columns whose 4th field is `+`/`-` ⇒ refFlat.
///
/// Returns `Ok(None)` if no line is decisive.
///
/// # Errors
/// Returns an error if a line cannot be read (I/O failure, or invalid UTF-8 from feeding a
/// compressed stream that was not transparently decompressed).
pub fn detect_format<R: BufRead>(reader: R) -> Result<Option<GeneModelFormat>> {
    for line in reader.lines() {
        let raw = line.context("reading gene model for format detection")?;
        let line = raw.trim_end();
        if line.is_empty() {
            continue;
        }
        if line.starts_with("##gff-version") {
            return Ok(Some(GeneModelFormat::Gff3));
        }
        if line.starts_with('#') {
            continue;
        }
        if line.starts_with("geneName") || line.starts_with("GENE_NAME") {
            return Ok(Some(GeneModelFormat::RefFlat));
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() == 9 {
            let attrs = cols[8];
            // GTF attributes are `key "value"; ...`; GFF3 are `key=value; ...`. A
            // value can legitimately contain the *other* delimiter (e.g. a GTF
            // `gene_id "AC=1"`, or a GFF3 `Note=he said "x"`), so decide by which
            // key/value separator appears first, not merely which char is present.
            return Ok(match (attrs.find('"'), attrs.find('=')) {
                (Some(q), Some(e)) => {
                    Some(if q < e { GeneModelFormat::Gtf } else { GeneModelFormat::Gff3 })
                }
                (Some(_), None) => Some(GeneModelFormat::Gtf),
                (None, Some(_)) => Some(GeneModelFormat::Gff3),
                (None, None) => None,
            });
        }
        if cols.len() >= 10 && matches!(cols[3], "+" | "-") {
            return Ok(Some(GeneModelFormat::RefFlat));
        }
        return Ok(None);
    }
    Ok(None)
}

// ─── Locus grouping + contig resolution ────────────────────────────────────────

/// Group parsed transcripts into loci and resolve contig names to BAM `ref_id`s.
///
/// Transcripts are grouped by gene name; within each gene they are clustered (ported
/// from fgbio `RefFlatSource`): processed longest-first, each transcript joins the first
/// existing cluster that contains a same-strand transcript on the same contig whose span
/// it overlaps, otherwise it starts a new cluster. Each cluster becomes a [`GeneLocus`].
/// Loci on contigs absent from `dict` (after reconciliation) are dropped and counted.
fn build_loci(transcripts: Vec<ParsedTranscript>, dict: &SequenceDictionary) -> Vec<GeneLocus> {
    // Preserve first-seen gene order for deterministic output.
    let mut order: Vec<String> = Vec::new();
    let mut by_gene: HashMap<String, Vec<ParsedTranscript>> = HashMap::new();
    for tx in transcripts {
        by_gene
            .entry(tx.gene_name.clone())
            .or_insert_with(|| {
                order.push(tx.gene_name.clone());
                Vec::new()
            })
            .push(tx);
    }

    let mut loci: Vec<GeneLocus> = Vec::new();
    let mut dropped_loci = 0usize;
    for gene_name in order {
        let mut txs = by_gene.remove(&gene_name).unwrap();
        // Longest-span transcript first, matching fgbio's clustering order.
        txs.sort_by_key(|tx| std::cmp::Reverse(tx.end - tx.start));

        let mut clusters: Vec<Vec<ParsedTranscript>> = Vec::new();
        for tx in txs {
            let slot = clusters.iter_mut().find(|cluster| {
                cluster.iter().any(|other| {
                    other.negative_strand == tx.negative_strand
                        && other.contig == tx.contig
                        && spans_overlap(other.start, other.end, tx.start, tx.end)
                })
            });
            match slot {
                Some(cluster) => cluster.push(tx),
                None => clusters.push(vec![tx]),
            }
        }

        for cluster in clusters {
            match build_locus(&gene_name, cluster, dict) {
                Some(locus) => loci.push(locus),
                None => dropped_loci += 1,
            }
        }
    }

    if dropped_loci > 0 {
        log::warn!(
            "rna: dropped {dropped_loci} gene loci on contigs not present in the BAM header"
        );
    }
    loci
}

/// Build a single [`GeneLocus`] from a cluster of transcripts, resolving the contig.
/// Returns `None` if the contig cannot be resolved against `dict`.
fn build_locus(
    gene_name: &str,
    cluster: Vec<ParsedTranscript>,
    dict: &SequenceDictionary,
) -> Option<GeneLocus> {
    let first = cluster.first()?;
    let ref_id = resolve_contig(dict, &first.contig)?;
    let negative_strand = first.negative_strand;
    // A gene is protein-coding if any of its transcripts is; otherwise take the first annotated
    // biotype. This makes the refFlat coding/non-coding inference order-independent for genes with
    // mixed transcripts, and is a no-op for GTF/GFF (whose transcripts share one gene biotype).
    let biotype = if cluster.iter().any(|tx| matches!(tx.biotype, Some(Biotype::ProteinCoding))) {
        Some(Biotype::ProteinCoding)
    } else {
        cluster.iter().find_map(|tx| tx.biotype.clone())
    };
    let start = cluster.iter().map(|tx| tx.start).min()?;
    let end = cluster.iter().map(|tx| tx.end).max()?;

    let transcripts: Vec<Transcript> = cluster
        .into_iter()
        .map(|tx| {
            let name = tx.name();
            Transcript {
                name,
                start: tx.start,
                end: tx.end,
                cds_start: tx.cds_start,
                cds_end: tx.cds_end,
                exons: normalize_exons(tx.exons),
            }
        })
        .collect();

    // Merge exon and coding (exon ∩ CDS) intervals across all transcripts so per-base
    // classification can be done once per locus rather than once per transcript.
    let mut exon_intervals = Vec::new();
    let mut coding_intervals = Vec::new();
    for tx in &transcripts {
        for exon in &tx.exons {
            exon_intervals.push((exon.start, exon.end));
            if let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) {
                let lo = exon.start.max(cds_start);
                let hi = exon.end.min(cds_end);
                if lo <= hi {
                    coding_intervals.push((lo, hi));
                }
            }
        }
    }

    Some(GeneLocus {
        gene_name: gene_name.to_string(),
        biotype,
        ref_id,
        start,
        end,
        negative_strand,
        transcripts,
        exon_union: merge_intervals(exon_intervals),
        coding_union: merge_intervals(coding_intervals),
    })
}

/// Sort a transcript's exons by start and coalesce any that overlap into a sorted, pairwise-disjoint
/// list. A transcript's exons are the disjoint spliced segments, so on real annotations this only
/// sorts; it exists to normalize malformed or duplicate exon records (from hand-edited or
/// pipeline-merged gene models) whose overlap would otherwise break the disjointness that
/// `bases_overlapping_exons` relies on. Exons that merely touch (`end + 1 == start`) are already
/// disjoint and left as separate exons.
fn normalize_exons(mut exons: Vec<Exon>) -> Vec<Exon> {
    exons.sort_by_key(|e| (e.start, e.end));
    let mut disjoint: Vec<Exon> = Vec::with_capacity(exons.len());
    for exon in exons {
        match disjoint.last_mut() {
            Some(last) if exon.start <= last.end => last.end = last.end.max(exon.end),
            _ => disjoint.push(exon),
        }
    }
    disjoint
}

/// Classify each alignment block (`(start, len)`, 1-based) against the sorted, disjoint
/// coding / exon / span interval sets, accumulating per-class base counts. `coding ⊆ exon ⊆
/// span` is guaranteed by construction, so every difference below is non-negative.
fn count_blocks(
    blocks: &[(u32, u32)],
    coding: &[(u32, u32)],
    exon: &[(u32, u32)],
    span: &[(u32, u32)],
) -> BaseClassCounts {
    let mut counts = BaseClassCounts::default();
    for &(block_start, block_len) in blocks {
        let block_end = block_start + block_len - 1;
        let coding_bases = overlap_with_intervals(block_start, block_end, coding);
        let exonic_bases = overlap_with_intervals(block_start, block_end, exon);
        let spanned_bases = overlap_with_intervals(block_start, block_end, span);
        counts.coding += u64::from(coding_bases);
        counts.utr += u64::from(exonic_bases - coding_bases);
        counts.intronic += u64::from(spanned_bases - exonic_bases);
        counts.intergenic += u64::from(block_len - spanned_bases);
    }
    counts
}

/// Total overlap (in bases) between the 1-based inclusive block `[block_start, block_end]` and a
/// sorted, disjoint set of 1-based inclusive `intervals`. Binary-searches the first candidate
/// interval, then scans until past the block.
fn overlap_with_intervals(block_start: u32, block_end: u32, intervals: &[(u32, u32)]) -> u32 {
    let lo = intervals.partition_point(|&(_, end)| end < block_start);
    let mut total = 0;
    for &(start, end) in &intervals[lo..] {
        if start > block_end {
            break;
        }
        total += block_end.min(end) - block_start.max(start) + 1;
    }
    total
}

/// Merge a set of 1-based inclusive intervals into a sorted, disjoint set. Adjacent or
/// overlapping intervals are coalesced.
fn merge_intervals(mut intervals: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    intervals.sort_unstable();
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        push_coalesced(&mut merged, interval);
    }
    merged
}

/// Union of several already-sorted, disjoint interval sets into one sorted, disjoint set — the same
/// result as [`merge_intervals`] on their concatenation, but without the `O(n log n)` sort, since
/// each input set is already ordered. Used on the hot multi-locus classification path, where each
/// overlapping gene contributes its pre-merged exon/coding union.
fn union_of_sorted_disjoint<'a>(sets: impl Iterator<Item = &'a [(u32, u32)]>) -> Vec<(u32, u32)> {
    let mut acc: Vec<(u32, u32)> = Vec::new();
    for set in sets {
        if acc.is_empty() {
            acc.extend_from_slice(set);
        } else if !set.is_empty() {
            acc = merge_two_sorted_disjoint(&acc, set);
        }
    }
    acc
}

/// Merge two already-sorted, disjoint interval sets into one, coalescing overlapping or adjacent
/// intervals. Both inputs must be sorted by start and internally disjoint.
fn merge_two_sorted_disjoint(a: &[(u32, u32)], b: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].0 <= b[j].0 {
            push_coalesced(&mut merged, a[i]);
            i += 1;
        } else {
            push_coalesced(&mut merged, b[j]);
            j += 1;
        }
    }
    for &interval in &a[i..] {
        push_coalesced(&mut merged, interval);
    }
    for &interval in &b[j..] {
        push_coalesced(&mut merged, interval);
    }
    merged
}

/// Append `interval` to an ascending run of disjoint intervals, coalescing it into the last one
/// when they overlap or are directly adjacent (`end + 1 == start`). `merged` must already be sorted
/// and `interval.0` must be `>=` every start already pushed.
fn push_coalesced(merged: &mut Vec<(u32, u32)>, interval: (u32, u32)) {
    match merged.last_mut() {
        Some(last) if interval.0 <= last.1 + 1 => last.1 = last.1.max(interval.1),
        _ => merged.push(interval),
    }
}

/// Resolve a gene-model contig name to a BAM `ref_id`, trying common naming variants:
/// exact match, then add/strip a leading `chr`, then the `MT`↔`chrM`/`M` aliases.
fn resolve_contig(dict: &SequenceDictionary, contig: &str) -> Option<usize> {
    if let Some(meta) = dict.get_by_name(contig) {
        return Some(meta.index());
    }
    // chr-prefix reconciliation.
    let alt = if let Some(stripped) = contig.strip_prefix("chr") {
        stripped.to_string()
    } else {
        format!("chr{contig}")
    };
    if let Some(meta) = dict.get_by_name(&alt) {
        return Some(meta.index());
    }
    // Mitochondrial aliases.
    let mt_aliases: &[&str] = match contig {
        "MT" | "M" => &["chrM", "chrMT"],
        "chrM" | "chrMT" => &["MT", "M"],
        _ => &[],
    };
    mt_aliases.iter().find_map(|name| dict.get_by_name(name).map(SequenceMetadata::index))
}

/// True if two 1-based inclusive spans overlap.
fn spans_overlap(a_start: u32, a_end: u32, b_start: u32, b_end: u32) -> bool {
    a_start <= b_end && b_start <= a_end
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Format detection ──────────────────────────────────────────────────────

    #[test]
    fn detect_refflat_by_data_line() {
        let line = "GENE\tNM_1\tchr1\t+\t100\t500\t150\t450\t2\t100,300,\t200,500,\n";
        assert_eq!(detect_format(line.as_bytes()).unwrap(), Some(GeneModelFormat::RefFlat));
    }

    #[test]
    fn detect_refflat_by_picard_header() {
        let content = "GENE_NAME\tTRANSCRIPT_NAME\tCHROMOSOME\tSTRAND\n";
        assert_eq!(detect_format(content.as_bytes()).unwrap(), Some(GeneModelFormat::RefFlat));
    }

    #[test]
    fn detect_gtf_by_attribute_style() {
        let line = "chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n";
        assert_eq!(detect_format(line.as_bytes()).unwrap(), Some(GeneModelFormat::Gtf));
    }

    #[test]
    fn detect_gff3_by_directive() {
        let content = "##gff-version 3\nchr1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n";
        assert_eq!(detect_format(content.as_bytes()).unwrap(), Some(GeneModelFormat::Gff3));
    }

    #[test]
    fn detect_gff3_by_attribute_style() {
        let line = "chr1\tx\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=rna1\n";
        assert_eq!(detect_format(line.as_bytes()).unwrap(), Some(GeneModelFormat::Gff3));
    }

    #[test]
    fn detect_gtf_when_a_quoted_value_contains_equals() {
        // A GTF attribute value containing '=' must not be misdetected as GFF3:
        // the quote precedes the '=', so GTF wins.
        let line = "chr1\tHAVANA\texon\t1\t9\t.\t+\t.\tgene_id \"G=1\"; transcript_id \"T1\";\n";
        assert_eq!(detect_format(line.as_bytes()).unwrap(), Some(GeneModelFormat::Gtf));
    }

    // ── Contig resolution ─────────────────────────────────────────────────────

    fn dict_of(contigs: &[&str]) -> SequenceDictionary {
        let owned: Vec<(String, usize)> = contigs.iter().map(|&n| (n.to_string(), 1000)).collect();
        SequenceDictionary::from(crate::test_support::SamBuilder::with_contigs(&owned).header())
    }

    #[test]
    fn resolve_contig_reconciles_chr_prefix_and_mt_aliases() {
        let dict = dict_of(&["chr1", "chrM"]);
        assert_eq!(resolve_contig(&dict, "chr1"), Some(0)); // exact
        assert_eq!(resolve_contig(&dict, "1"), Some(0)); // unprefixed model → chr-prefixed BAM
        assert_eq!(resolve_contig(&dict, "MT"), Some(1)); // MT alias → chrM
        assert_eq!(resolve_contig(&dict, "M"), Some(1)); // M alias → chrM
        assert_eq!(resolve_contig(&dict, "chr99"), None); // unresolvable → dropped
    }

    #[test]
    fn resolve_contig_strips_chr_prefix_for_unprefixed_dict() {
        let dict = dict_of(&["1", "MT"]);
        assert_eq!(resolve_contig(&dict, "chr1"), Some(0)); // strip chr → 1
        assert_eq!(resolve_contig(&dict, "chrM"), Some(1)); // chrM → MT alias
    }

    // ── Transcript ────────────────────────────────────────────────────────────

    #[test]
    fn transcribed_length_sums_exons_excluding_introns() {
        // Exons [100,200] and [300,400] are 101 bases each; the intron is excluded.
        let tx = Transcript {
            name: "tx".to_string(),
            start: 100,
            end: 400,
            cds_start: Some(150),
            cds_end: Some(350),
            exons: vec![Exon { start: 100, end: 200 }, Exon { start: 300, end: 400 }],
        };
        assert_eq!(tx.transcribed_length(), 101 + 101);
    }

    #[test]
    fn union_coverage_counts_every_base_and_reconstructs_per_transcript() {
        // Two transcripts of one locus sharing exon [100,200]; tx1 also has [300,400], tx2 [500,600].
        let mk = |name: &str, exons: Vec<Exon>| Transcript {
            name: name.to_string(),
            start: exons.first().unwrap().start,
            end: exons.last().unwrap().end,
            cds_start: None,
            cds_end: None,
            exons,
        };
        let tx1 = mk("t1", vec![Exon { start: 100, end: 200 }, Exon { start: 300, end: 400 }]);
        let tx2 = mk("t2", vec![Exon { start: 100, end: 200 }, Exon { start: 500, end: 600 }]);
        let locus = GeneLocus {
            gene_name: "g".to_string(),
            biotype: None,
            ref_id: 0,
            start: 100,
            end: 600,
            negative_strand: false,
            transcripts: vec![tx1, tx2],
            exon_union: vec![(100, 200), (300, 400), (500, 600)],
            coding_union: vec![],
        };

        let mut cov = vec![0u32; locus.exon_union_len()];
        assert_eq!(cov.len(), 101 + 101 + 101);
        // A read covering all of the shared exon [100,200]: every base counted once, last included.
        locus.add_union_coverage(100, 200, &mut cov);
        assert!(cov[..101].iter().all(|&c| c == 1));

        // Reconstruction: both transcripts see the shared exon's coverage (accumulated once).
        let c1 = locus.transcript_coverage(&locus.transcripts[0], &cov);
        let c2 = locus.transcript_coverage(&locus.transcripts[1], &cov);
        assert!(c1[..101].iter().all(|&c| c == 1)); // shared exon covered
        assert!(c1[101..].iter().all(|&c| c == 0)); // tx1's second exon uncovered
        assert!(c2[..101].iter().all(|&c| c == 1)); // same shared exon
        assert!(c2[101..].iter().all(|&c| c == 0)); // tx2's second exon uncovered
    }

    // ── Block-level base classification (GeneModel::classify_blocks) ───────────

    /// Build a minimal locus with the given span and merged exon/coding unions. Transcripts are
    /// not needed by `classify_blocks` (only the unions and span are read).
    fn locus_with_unions(
        start: u32,
        end: u32,
        exon_union: Vec<(u32, u32)>,
        coding_union: Vec<(u32, u32)>,
    ) -> GeneLocus {
        GeneLocus {
            gene_name: "g".to_string(),
            biotype: None,
            ref_id: 0,
            start,
            end,
            negative_strand: false,
            transcripts: Vec::new(),
            exon_union,
            coding_union,
        }
    }

    fn model_of(loci: Vec<GeneLocus>) -> GeneModel {
        let max_ref = loci.iter().map(|l| l.ref_id).max().unwrap_or(0);
        let mut loci_by_contig = vec![Vec::new(); max_ref + 1];
        for (i, locus) in loci.iter().enumerate() {
            loci_by_contig[locus.ref_id].push(u32::try_from(i).unwrap());
        }
        for list in &mut loci_by_contig {
            list.sort_by_key(|&i: &u32| loci[i as usize].start);
        }
        GeneModel { loci, loci_by_contig }
    }

    #[test]
    fn classify_blocks_no_loci_is_all_intergenic() {
        let model = model_of(Vec::new());
        let counts = model.classify_blocks(&[], &[(100, 50), (300, 20)]);
        assert_eq!(counts, BaseClassCounts { intergenic: 70, ..BaseClassCounts::default() });
    }

    #[test]
    fn classify_blocks_single_locus_splits_by_function() {
        // Span [100,400]; exons [100,200] & [300,400]; CDS [150,350] → coding [150,200] & [300,350].
        let model = model_of(vec![locus_with_unions(
            100,
            400,
            vec![(100, 200), (300, 400)],
            vec![(150, 200), (300, 350)],
        )]);
        // One block [50,251]: 50 bases before the gene, then into exon 1 and its intron.
        let counts = model.classify_blocks(&[0], &[(50, 202)]);
        assert_eq!(counts, BaseClassCounts { coding: 51, utr: 50, intronic: 51, intergenic: 50 });
    }

    #[test]
    fn classify_blocks_unions_overlapping_genes() {
        // Two overlapping loci: a base coding in either counts as coding; exonic in either as ≥UTR.
        let a = locus_with_unions(100, 300, vec![(100, 150)], vec![]); // non-coding
        let b = locus_with_unions(200, 400, vec![(250, 350)], vec![(280, 320)]);
        let model = model_of(vec![a, b]);
        // Block [120,360]: merged span [100,400], exons {[100,150],[250,350]}, coding [280,320].
        let counts = model.classify_blocks(&[0, 1], &[(120, 241)]);
        assert_eq!(counts, BaseClassCounts { coding: 41, utr: 91, intronic: 109, intergenic: 0 });
    }

    #[test]
    fn overlap_with_intervals_sums_disjoint_overlaps() {
        let intervals = [(100, 200), (300, 400)];
        assert_eq!(overlap_with_intervals(150, 350, &intervals), 51 + 51); // [150,200]+[300,350]
        assert_eq!(overlap_with_intervals(210, 290, &intervals), 0); // intronic gap
        assert_eq!(overlap_with_intervals(50, 60, &intervals), 0); // before all
    }

    #[test]
    fn merge_two_sorted_disjoint_interleaves_and_coalesces() {
        // Interleaved, non-touching intervals keep both, in order.
        assert_eq!(
            merge_two_sorted_disjoint(&[(100, 200), (500, 600)], &[(300, 400)]),
            vec![(100, 200), (300, 400), (500, 600)]
        );
        // Overlapping across the two sets coalesce into one.
        assert_eq!(merge_two_sorted_disjoint(&[(100, 200)], &[(150, 300)]), vec![(100, 300)]);
        // Directly adjacent (end + 1 == start) coalesce, matching merge_intervals.
        assert_eq!(merge_two_sorted_disjoint(&[(100, 200)], &[(201, 300)]), vec![(100, 300)]);
        // A duplicate interval present in both sets is deduplicated.
        assert_eq!(merge_two_sorted_disjoint(&[(100, 200)], &[(100, 200)]), vec![(100, 200)]);
    }

    #[test]
    fn union_of_sorted_disjoint_matches_merge_intervals() {
        // The k-way union of pre-sorted sets equals sorting-and-merging their concatenation.
        let sets: [&[(u32, u32)]; 3] =
            [&[(100, 150), (400, 500)], &[(120, 200), (600, 700)], &[(180, 190)]];
        let via_union = union_of_sorted_disjoint(sets.iter().copied());
        let via_sort = merge_intervals(sets.iter().flat_map(|s| s.iter().copied()).collect());
        assert_eq!(via_union, via_sort);
        assert_eq!(via_union, vec![(100, 200), (400, 500), (600, 700)]);
    }

    #[test]
    fn union_of_sorted_disjoint_handles_empty_and_single_sets() {
        let empty: [&[(u32, u32)]; 0] = [];
        assert!(union_of_sorted_disjoint(empty.iter().copied()).is_empty());

        // Empty sets interspersed are skipped without affecting the result.
        let sets: [&[(u32, u32)]; 3] = [&[], &[(100, 200)], &[]];
        assert_eq!(union_of_sorted_disjoint(sets.iter().copied()), vec![(100, 200)]);
    }

    #[test]
    fn normalize_exons_sorts_and_coalesces_overlaps() {
        let ex = |s, e| Exon { start: s, end: e };

        // Already-disjoint exons given out of order are only sorted, not merged.
        assert_eq!(
            normalize_exons(vec![ex(300, 400), ex(100, 200)]),
            vec![ex(100, 200), ex(300, 400)]
        );

        // Exons that merely touch (end + 1 == start) stay separate.
        assert_eq!(
            normalize_exons(vec![ex(100, 150), ex(151, 200)]),
            vec![ex(100, 150), ex(151, 200)]
        );

        // Overlapping, nested, shared-endpoint, and duplicate exons all coalesce.
        assert_eq!(normalize_exons(vec![ex(100, 200), ex(150, 300)]), vec![ex(100, 300)]);
        assert_eq!(normalize_exons(vec![ex(100, 500), ex(200, 300)]), vec![ex(100, 500)]);
        assert_eq!(normalize_exons(vec![ex(100, 150), ex(150, 200)]), vec![ex(100, 200)]);
        assert_eq!(normalize_exons(vec![ex(100, 200), ex(100, 200)]), vec![ex(100, 200)]);
    }
}
