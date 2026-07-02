//! The `rna` command: RNA-seq QC metrics.
//!
//! Ports Picard `CollectRnaSeqMetrics` (base composition, strand specificity, transcript
//! coverage/bias) and fgbio `EstimateRnaSeqInsertSize` (insert size in transcript space)
//! into a single pass over the reads sharing one gene-model load. See `CONTRIBUTING.md`
//! and the module docs in [`crate::gene_model`] for the gene-model formats supported.
//!
//! ## Filtering (unified, Picard-biased)
//!
//! In order: vendor-QC-fail reads are dropped; the read length of every primary
//! (non-secondary) read — mapped or not — is added to `bases`; reads on an `--ignore-chrom`
//! contig are counted only as `ignored_reads`; secondary and unmapped reads are then dropped.
//! Duplicates are **included by default** (RNA-seq coordinate duplicates are usually
//! independent cDNAs from abundant transcripts, not PCR artifacts); `--exclude-duplicates`
//! drops them from every metric. `--min-mapq` (default 0, Picard parity) gates the
//! alignment-based metrics but not `bases`. The observed `duplicate_rate` is always reported.
//!
//! ## Deviations from the reference tools (all bugs fixed + documented)
//!
//! - **Coverage top-N**: Picard keeps the top *1001* most-expressed transcripts
//!   (`coverages[length - 1001]`); we keep the top 1000.
//! - **Short-transcript bias**: Picard's 5'/3' bias windows (each `end_bias_bases` wide) can
//!   overlap for transcripts shorter than `2 × end_bias_bases`. We exclude such transcripts
//!   from the CV/bias set rather than altering the formula.
//! - **Coverage undercount**: Picard's `addCoverageCounts` iterates a half-open range against
//!   an inclusive end and drops the last base of every alignment block; we count every base.
//! - **Ribosomal overlap**: Picard takes the single largest rRNA-interval intersection per
//!   fragment (`Math.max`); we take the union of (merged, disjoint) overlaps.
//! - **Fragment end**: Picard derives the rRNA fragment end from TLEN (degenerate when TLEN
//!   is 0); we use the mate end from the `MC` tag, falling back to TLEN, never read length.
//! - **Coverage-vs-position bins**: Picard's per-percent windows overlap at boundaries; we
//!   sample 101 non-overlapping normalized positions.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::{Args, ValueEnum};
use noodles::sam::Header;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};
use riker_derive::MetricDocs;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{
    InputOptions, OptionalReferenceOptions, OutputOptions, parse_existing_file,
};
use crate::counter::Counter;
use crate::gene_model::{Biotype, GeneModel, encloses_inclusive, overlap_inclusive};
use crate::intervals::Intervals;
use crate::metrics::{serialize_f64_5dp, serialize_opt_f64_5dp, write_tsv};
use crate::overlapper::Overlapper;
use crate::plotting::{
    FG_BLUE, FG_GREEN, FG_NAVY, FG_TEAL, InsertSizeSeries, PLOT_HEIGHT, PLOT_WIDTH,
    write_insert_size_histogram_pdf, write_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::pair_orientation::{PairOrientation, get_pair_orientation};
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sequence_dict::SequenceDictionary;

// ─── Output file suffixes ──────────────────────────────────────────────────────

/// Suffix for the summary metrics file.
pub const METRICS_SUFFIX: &str = ".rna-metrics.txt";
/// Suffix for the per-biotype assigned-read counts file.
pub const BIOTYPE_SUFFIX: &str = ".rna-biotype.txt";
/// Suffix for the per-biotype assigned-read bar chart.
pub const BIOTYPE_PLOT_SUFFIX: &str = ".rna-biotype.pdf";
/// Suffix for the per-gene expression (RPK) distribution histogram.
pub const GENE_EXPRESSION_PLOT_SUFFIX: &str = ".rna-gene-expression.pdf";
/// Suffix for the normalized-coverage-by-position chart.
pub const COVERAGE_PLOT_SUFFIX: &str = ".rna-coverage.pdf";
/// Suffix for the transcript-space insert-size metrics file.
pub const ISIZE_SUFFIX: &str = ".rna-insert-size.txt";
/// Suffix for the transcript-space insert-size histogram.
pub const ISIZE_HISTOGRAM_SUFFIX: &str = ".rna-insert-size-histogram.txt";
/// Suffix for the transcript-space insert-size histogram chart.
pub const ISIZE_PLOT_SUFFIX: &str = ".rna-insert-size-histogram.pdf";

/// `MC` (mate CIGAR) aux tag, required to compute transcript-space insert size.
const MC_TAG: [u8; 2] = *b"MC";

// ─── Options ─────────────────────────────────────────────────────────────────

/// Library strand-specificity. `Auto` infers the orientation from the observed
/// read-1/read-2 transcript-strand counts during the pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum StrandSpec {
    /// Infer strandedness from the data (default).
    #[default]
    Auto,
    /// Unstranded library: correct/incorrect-strand metrics are zero.
    Unstranded,
    /// First read is on the transcription strand (Picard `FIRST_READ_TRANSCRIPTION_STRAND`).
    Forward,
    /// Second read is on the transcription strand (Picard `SECOND_READ_TRANSCRIPTION_STRAND`).
    Reverse,
}

impl std::fmt::Display for StrandSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Auto => "auto",
            Self::Unstranded => "unstranded",
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        };
        f.write_str(s)
    }
}

/// Tool-specific options for the RNA-seq collector.
#[riker_derive::multi_options("rna", "RNA-seq Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct RnaOptions {
    /// Gene model (required): UCSC refFlat, GTF, or GFF3 (GENCODE/Ensembl/RefSeq), plain or
    /// gzip/bgzip-compressed. Format is auto-detected from the contents. See "GENE MODEL" above.
    #[arg(long, value_name = "FILE", value_parser = parse_existing_file)]
    pub gene_model: PathBuf,

    /// Explicit ribosomal intervals (BED or Picard IntervalList), unioned with biotype-derived
    /// rRNA genes from the gene model. See "RIBOSOMAL REGIONS" above.
    #[arg(long, value_name = "FILE", value_parser = parse_existing_file)]
    pub ribosomal_intervals: Option<PathBuf>,

    /// Library strandedness. `auto` (default) infers forward/reverse/unstranded from the data and
    /// reports it as `detected_strand`; the raw read-1/read-2 strand counts are always emitted.
    #[arg(long, value_enum, default_value_t = StrandSpec::Auto)]
    pub strand: StrandSpec,

    /// Ignore reads with mapping quality below this value in the alignment-based metrics.
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Exclude duplicate-flagged reads from all metrics. Off by default (RNA-seq duplicates
    /// are usually independent fragments from abundant transcripts).
    #[arg(long, default_value_t = false)]
    pub exclude_duplicates: bool,

    /// Minimum fraction of a fragment's bases that must overlap ribosomal intervals for the whole
    /// fragment to be counted as ribosomal (toward `ribosomal_bases` / `frac_ribosomal_bases`).
    #[arg(long, default_value_t = 0.8)]
    pub rrna_fragment_fraction: f64,

    /// Minimum transcript length (bp) for a transcript to be eligible for the coverage/bias
    /// metrics: `median_cv_coverage`, `median_5prime_bias`, `median_3prime_bias`,
    /// `median_5prime_to_3prime_bias`, and the normalized coverage-vs-position table/plot. Does not
    /// affect base composition, strand, or insert size. The effective threshold is
    /// `max(minimum_length, 2 × end_bias_bases)` so the 5'/3' windows can't overlap.
    #[arg(long, default_value_t = 500)]
    pub minimum_length: u32,

    /// Width (bp) of the 5' and 3' end windows used to compute `median_5prime_bias`,
    /// `median_3prime_bias`, and `median_5prime_to_3prime_bias`.
    #[arg(long, default_value_t = 100)]
    pub end_bias_bases: u32,

    /// Minimum fraction of read+mate bases that must overlap a transcript's exons for the pair's
    /// insert size to be trusted and counted.
    #[arg(long, default_value_t = 0.95)]
    pub isize_minimum_overlap: f64,

    /// Trim insert sizes above `median + deviations × MAD` (per orientation) before computing the
    /// mean and standard deviation and before plotting the histogram.
    #[arg(long, default_value_t = 10.0)]
    pub isize_deviations: f64,

    /// Chromosome/contig name(s) to ignore: reads mapped there count only as `ignored_reads` (and
    /// toward `bases`) and contribute to no other metric. Names must match the BAM header (e.g.
    /// `chrM`). Mirrors Picard's `IGNORE_SEQUENCE`. Repeatable.
    #[arg(long, value_name = "CONTIG", num_args = 1..)]
    pub ignore_chrom: Option<Vec<String>>,

    /// Minimum number of uniquely-assigned reads for a gene to count toward `genes_detected`.
    #[arg(long, default_value_t = 5)]
    pub genes_detected_min_reads: u32,

    /// Minimum CIGAR N (skip) length, in bp, to treat a gap as a splice junction rather than a
    /// short deletion (affects the splice-junction metrics).
    #[arg(long, default_value_t = 50)]
    pub junction_min_intron: u32,

    /// Minimum mean per-base coverage for a transcript to contribute to the TIN (Transcript
    /// Integrity Number) statistics; transcripts shorter than --minimum-length are also excluded.
    /// Filters out lowly-covered, noisy transcripts whose integrity cannot be estimated reliably.
    #[arg(long, default_value_t = 10)]
    pub tin_min_coverage: u32,
}

impl Default for RnaOptions {
    fn default() -> Self {
        Self {
            gene_model: PathBuf::new(),
            ribosomal_intervals: None,
            strand: StrandSpec::Auto,
            min_mapq: 0,
            exclude_duplicates: false,
            rrna_fragment_fraction: 0.8,
            minimum_length: 500,
            end_bias_bases: 100,
            isize_minimum_overlap: 0.95,
            isize_deviations: 10.0,
            ignore_chrom: None,
            genes_detected_min_reads: 5,
            junction_min_intron: 50,
            tin_min_coverage: 10,
        }
    }
}

// ─── Command ─────────────────────────────────────────────────────────────────

/// Collect RNA-seq QC metrics in a single pass: read accounting, read genomic origin (exonic / intronic / intergenic / ambiguous / ribosomal), base composition (coding / UTR / intronic / intergenic / ribosomal), gene detection, splice-junction annotation (known / partial-novel / novel), strand specificity, transcript 5'/3' coverage bias, transcript integrity (TIN), and transcript-space insert size. Outputs are written to <prefix>.rna-metrics.txt, <prefix>.rna-biotype.txt(.pdf), <prefix>.rna-gene-expression.pdf, <prefix>.rna-coverage.pdf, and <prefix>.rna-insert-size.txt / -histogram.txt(.pdf).
///
/// The input SAM/BAM/CRAM must be coordinate-sorted (@HD SO:coordinate); duplicates are included by default (--exclude-duplicates to drop them).
///
/// GENE MODEL: UCSC refFlat, GTF, or GFF3 (GENCODE / Ensembl / RefSeq), plain or gzip/bgzip-compressed (compression is taken from a .gz/.bgz extension). The format is auto-detected from the file contents. Contig names are reconciled to the BAM header automatically (add/strip 'chr', alias MT<->chrM, and map RefSeq accessions such as NC_000001.11 to common names).
///
/// RIBOSOMAL REGIONS: the union of genes with biotype rRNA / Mt_rRNA / rRNA_pseudogene from the gene model (GTF and GFF3 carry biotypes and contribute; refFlat does not) and an explicit --ribosomal-intervals file (BED or Picard IntervalList). With neither source the ribosomal columns are left blank. Annotation-derived rRNA is a lower bound on GRCh38, where the rDNA tandem arrays are collapsed; supply curated --ribosomal-intervals for an accurate fraction.
///
/// INSERT SIZE: the 5'-to-5' fragment length in transcript coordinates (introns collapsed), per pair orientation. Uses the MC (mate CIGAR) tag when present (any orientation, exact); without it, FR pairs are estimated from a spec-compliant TLEN field (non-FR pairs are skipped). For best results add MC with `samtools fixmate -m`. A BAM with neither MC nor valid TLEN yields no/incorrect insert sizes.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker rna -i in.bam -o out --gene-model gencode.gtf.gz
  riker rna -i in.bam -o out --gene-model refFlat.txt --ribosomal-intervals rRNA.interval_list
  riker rna -i in.bam -o out --gene-model genes.gff3 --strand reverse --ignore-chrom chrM"
)]
pub struct Rna {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,

    #[command(flatten)]
    pub options: RnaOptions,
}

impl Command for Rna {
    /// # Errors
    /// Returns an error if the input or gene model cannot be read, or output cannot be written.
    fn execute(&self, threads: Option<u8>) -> Result<()> {
        super::common::validate_output_prefix(&self.output.output)?;
        let plan = self.thread_plan(threads);
        let mut reader = AlignmentReader::open(
            &self.input.input,
            self.reference.reference.as_deref(),
            plan.decode_threads,
        )?;
        let mut collector =
            RnaCollector::new(&self.input.input, &self.output.output, &self.options);
        let mut progress = ProgressLogger::new("rna", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

// ─── Collector ───────────────────────────────────────────────────────────────

/// Collects RNA-seq metrics over a single pass. Configuration is set at construction;
/// the gene model and ribosomal intervals are loaded in [`Collector::initialize`] once the
/// header (and thus the sequence dictionary) is available.
pub struct RnaCollector {
    // Configuration.
    input: PathBuf,
    metrics_path: PathBuf,
    biotype_path: PathBuf,
    biotype_plot_path: PathBuf,
    gene_expression_plot_path: PathBuf,
    coverage_plot_path: PathBuf,
    isize_metrics_path: PathBuf,
    isize_histogram_path: PathBuf,
    isize_plot_path: PathBuf,
    gene_model_path: PathBuf,
    ribosomal_path: Option<PathBuf>,
    strand: StrandSpec,
    min_mapq: u8,
    exclude_duplicates: bool,
    rrna_fragment_fraction: f64,
    minimum_length: u32,
    end_bias_bases: u32,
    isize_minimum_overlap: f64,
    isize_deviations: f64,
    ignore_chrom: Vec<String>,
    genes_detected_min_reads: u32,
    junction_min_intron: u32,
    tin_min_coverage: u32,

    // Initialized from the header.
    sample: String,
    gene_model: Option<GeneModel>,
    /// Merged, disjoint ribosomal intervals as `(start0, end0)` per the value; `None` when no
    /// ribosomal source is available (matches Picard leaving ribosomal metrics blank).
    ribosomal: Option<Overlapper<RibosomalInterval>>,
    ignored_ref_ids: HashSet<usize>,

    // Accumulated metrics.
    bases: u64,
    aligned_bases: u64,
    ribosomal_bases: u64,
    coding_bases: u64,
    utr_bases: u64,
    intronic_bases: u64,
    intergenic_bases: u64,
    ignored_reads: u64,
    // Read-level accounting and genomic-origin classification (primary, mapped, non-supplementary
    // reads reaching classification). `exonic_reads` = `assigned_reads` (exonic to one gene) +
    // `ambiguous_reads` (exonic to >1 gene); the four categories ribosomal/exonic/intronic/
    // intergenic partition those reads under default filters.
    total_reads: u64,
    unmapped_reads: u64,
    ribosomal_reads: u64,
    exonic_reads: u64,
    intronic_reads: u64,
    intergenic_reads: u64,
    ambiguous_reads: u64,
    assigned_reads: u64,
    /// Per-locus uniquely-assigned read counts, sized to the gene model in `initialize`;
    /// aggregated by gene name at `finish` for `genes_detected`.
    assigned_counts: Vec<u32>,
    // Annotated splice junctions and their boundary sites, built from the gene model in
    // `initialize`; an observed intron is known if the whole junction matches, partial if exactly
    // one boundary matches, novel otherwise.
    known_introns: FxHashSet<JunctionKey>,
    known_intron_starts: FxHashSet<SiteKey>,
    known_intron_ends: FxHashSet<SiteKey>,
    /// Distinct observed introns; classified into known/partial/novel at `finish`.
    observed_introns: FxHashSet<JunctionKey>,
    spliced_reads: u64,
    known_splice_obs: u64,
    partial_splice_obs: u64,
    novel_splice_obs: u64,
    // Strand 2×2 table over single-locus exon-overlapping non-supplementary reads, split by
    // whether the read is read-1/unpaired and whether read and transcript strands agree.
    agree_r1: u64,
    disagree_r1: u64,
    agree_r2: u64,
    disagree_r2: u64,
    num_r1: u64,
    num_r2: u64,
    num_unexplained: u64,
    mapped_reads: u64,
    duplicate_reads: u64,
    // Per-locus coverage over the locus's concatenated exon union (introns removed, exons shared
    // by transcripts stored once), keyed by locus index. Per-transcript coverage for the bias/CV
    // metrics is reconstructed from this at `finish`.
    coverage: FxHashMap<u32, Vec<u32>>,
    // Transcript-space insert sizes per orientation.
    fr: Counter<u64>,
    rf: Counter<u64>,
    tandem: Counter<u64>,
    mc_missing_warned: bool,
    // Coordinate-sweep state for locus lookup (replaces interval-tree queries; requires
    // coordinate-sorted input). `sweep_open` is the index into the current contig's sorted locus
    // list up to which loci have been opened; `sweep_active` holds opened, not-yet-ended loci.
    sweep_ref: Option<usize>,
    sweep_open: usize,
    sweep_last_start: u32,
    sweep_active: Vec<u32>,
}

/// A ribosomal interval as 0-based half-open `(start, end)` bounds, stored as the overlap
/// value so the fragment-overlap fraction can be computed.
type RibosomalInterval = (u32, u32);

/// A splice junction as `(ref_id, intron_start_1based, intron_end_1based)`.
type JunctionKey = (u32, u32, u32);
/// A single splice-site boundary as `(ref_id, position_1based)`.
type SiteKey = (u32, u32);

/// Classification of an observed splice junction against the annotated set.
enum JunctionClass {
    /// Both splice sites match an annotated junction.
    Known,
    /// Exactly one splice site is annotated.
    Partial,
    /// Neither splice site is annotated.
    Novel,
}

impl RnaCollector {
    /// Create a new collector; output paths are derived from `prefix`.
    #[must_use]
    pub fn new(input: &Path, prefix: &Path, options: &RnaOptions) -> Self {
        use super::command::output_path;
        Self {
            input: input.to_path_buf(),
            metrics_path: output_path(prefix, METRICS_SUFFIX),
            biotype_path: output_path(prefix, BIOTYPE_SUFFIX),
            biotype_plot_path: output_path(prefix, BIOTYPE_PLOT_SUFFIX),
            gene_expression_plot_path: output_path(prefix, GENE_EXPRESSION_PLOT_SUFFIX),
            coverage_plot_path: output_path(prefix, COVERAGE_PLOT_SUFFIX),
            isize_metrics_path: output_path(prefix, ISIZE_SUFFIX),
            isize_histogram_path: output_path(prefix, ISIZE_HISTOGRAM_SUFFIX),
            isize_plot_path: output_path(prefix, ISIZE_PLOT_SUFFIX),
            gene_model_path: options.gene_model.clone(),
            ribosomal_path: options.ribosomal_intervals.clone(),
            strand: options.strand,
            min_mapq: options.min_mapq,
            exclude_duplicates: options.exclude_duplicates,
            rrna_fragment_fraction: options.rrna_fragment_fraction,
            minimum_length: options.minimum_length,
            end_bias_bases: options.end_bias_bases,
            isize_minimum_overlap: options.isize_minimum_overlap,
            isize_deviations: options.isize_deviations,
            ignore_chrom: options.ignore_chrom.clone().unwrap_or_default(),
            genes_detected_min_reads: options.genes_detected_min_reads,
            junction_min_intron: options.junction_min_intron,
            tin_min_coverage: options.tin_min_coverage,
            sample: String::new(),
            gene_model: None,
            ribosomal: None,
            ignored_ref_ids: HashSet::new(),
            bases: 0,
            aligned_bases: 0,
            ribosomal_bases: 0,
            coding_bases: 0,
            utr_bases: 0,
            intronic_bases: 0,
            intergenic_bases: 0,
            ignored_reads: 0,
            total_reads: 0,
            unmapped_reads: 0,
            ribosomal_reads: 0,
            exonic_reads: 0,
            intronic_reads: 0,
            intergenic_reads: 0,
            ambiguous_reads: 0,
            assigned_reads: 0,
            assigned_counts: Vec::new(),
            known_introns: FxHashSet::default(),
            known_intron_starts: FxHashSet::default(),
            known_intron_ends: FxHashSet::default(),
            observed_introns: FxHashSet::default(),
            spliced_reads: 0,
            known_splice_obs: 0,
            partial_splice_obs: 0,
            novel_splice_obs: 0,
            agree_r1: 0,
            disagree_r1: 0,
            agree_r2: 0,
            disagree_r2: 0,
            num_r1: 0,
            num_r2: 0,
            num_unexplained: 0,
            mapped_reads: 0,
            duplicate_reads: 0,
            coverage: FxHashMap::default(),
            fr: Counter::new(),
            rf: Counter::new(),
            tandem: Counter::new(),
            mc_missing_warned: false,
            sweep_ref: None,
            sweep_open: 0,
            sweep_last_start: 0,
            sweep_active: Vec::new(),
        }
    }

    /// Build the merged, disjoint ribosomal overlapper from the explicit interval file (if
    /// any) and biotype-derived rRNA loci. Returns `None` when no ribosomal source exists.
    fn build_ribosomal(
        &self,
        model: &GeneModel,
        dict: &SequenceDictionary,
    ) -> Result<Option<Overlapper<RibosomalInterval>>> {
        let mut intervals: Vec<(usize, u32, u32)> = model.ribosomal_intervals();
        let mut have_source = !intervals.is_empty();

        if let Some(path) = &self.ribosomal_path {
            have_source = true;
            let explicit = Intervals::from_path(path, dict.clone())?;
            intervals.extend(explicit.iter().map(|iv| (iv.ref_id, iv.start, iv.end)));
        }

        if !have_source {
            return Ok(None);
        }

        // Merge into disjoint intervals per contig so the per-fragment overlap is a true union.
        intervals.sort_unstable();
        let mut merged: Vec<(usize, u32, u32)> = Vec::new();
        for (ref_id, start, end) in intervals {
            match merged.last_mut() {
                Some(last) if last.0 == ref_id && start <= last.2 => last.2 = last.2.max(end),
                _ => merged.push((ref_id, start, end)),
            }
        }

        let overlapper =
            Overlapper::new(merged.into_iter().map(|(ref_id, s, e)| (ref_id, s, e, (s, e))));
        Ok(Some(overlapper))
    }

    /// Advance the coordinate sweep to a read spanning `[read_start, read_end]` (1-based
    /// inclusive) on `ref_id` and write the indices of overlapping gene loci into `out`.
    ///
    /// Loci are opened as their start passes the read end and dropped once their end falls behind
    /// the read start; since read starts only advance (coordinate-sorted input), this is a single
    /// forward pass over the contig's loci. The retained set can over-include loci opened by an
    /// earlier long (spliced) read, so it is filtered to those actually overlapping the read.
    ///
    /// # Errors
    /// Returns an error if a record arrives out of coordinate order on the current contig.
    fn overlapping_loci(
        &mut self,
        ref_id: usize,
        read_start: u32,
        read_end: u32,
        out: &mut Vec<u32>,
    ) -> Result<()> {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        if self.sweep_ref != Some(ref_id) {
            self.sweep_ref = Some(ref_id);
            self.sweep_open = 0;
            self.sweep_last_start = 0;
            self.sweep_active.clear();
        }
        anyhow::ensure!(
            read_start >= self.sweep_last_start,
            "rna: records are out of coordinate order (start {read_start} after \
             {} on the same contig); the input must be coordinate-sorted",
            self.sweep_last_start
        );
        self.sweep_last_start = read_start;

        let contig_loci = model.loci_on_contig(ref_id);
        while self.sweep_open < contig_loci.len()
            && model.locus(contig_loci[self.sweep_open]).start <= read_end
        {
            self.sweep_active.push(contig_loci[self.sweep_open]);
            self.sweep_open += 1;
        }
        self.sweep_active.retain(|&l| model.locus(l).end >= read_start);
        out.clear();
        out.extend(self.sweep_active.iter().copied().filter(|&l| model.locus(l).start <= read_end));
        Ok(())
    }

    /// Classify the aligned bases of a primary, mapped read against the overlapping `loci`, tally
    /// per-transcript coverage, and determine the read's genomic origin. Returns whether the read
    /// overlapped any exon (for the strand logic) and its [`ReadOrigin`].
    fn classify_bases(&mut self, record: &RikerRecord, loci: &[u32]) -> (bool, ReadOrigin) {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let blocks = rec_alignment_blocks(record);

        // Functional base composition via interval arithmetic (block ∩ exon/coding/span sets).
        let counts = model.classify_blocks(loci, &blocks);
        self.coding_bases += counts.coding;
        self.utr_bases += counts.utr;
        self.intronic_bases += counts.intronic;
        self.intergenic_bases += counts.intergenic;
        self.aligned_bases += counts.coding + counts.utr + counts.intronic + counts.intergenic;
        let overlaps_exon = counts.coding > 0 || counts.utr > 0;

        // Read-level origin. The common single-locus case avoids any per-exon walk; only reads
        // spanning multiple loci pay to resolve how many distinct genes their exons touch.
        let origin = if loci.is_empty() {
            ReadOrigin::Intergenic
        } else if !overlaps_exon {
            ReadOrigin::Intronic
        } else if loci.len() == 1 {
            ReadOrigin::Exonic(loci[0])
        } else {
            let mut exon_locus = None;
            let mut exon_gene: Option<&str> = None;
            let mut ambiguous = false;
            for &locus_idx in loci {
                let locus = model.locus(locus_idx);
                let hits_exon = locus
                    .transcripts
                    .iter()
                    .any(|tx| bases_overlapping_exons(&blocks, &tx.exons) > 0);
                if hits_exon {
                    match exon_gene {
                        None => {
                            exon_gene = Some(&locus.gene_name);
                            exon_locus = Some(locus_idx);
                        }
                        Some(name) if name != locus.gene_name => ambiguous = true,
                        Some(_) => {}
                    }
                }
            }
            match (ambiguous, exon_locus) {
                (true, _) => ReadOrigin::Ambiguous,
                (false, Some(locus_idx)) => ReadOrigin::Exonic(locus_idx),
                (false, None) => ReadOrigin::Intronic,
            }
        };

        // Coverage for the bias/CV metrics: accumulate once per locus into its exon-union array
        // (per-transcript coverage is reconstructed at `finish`).
        for &locus_idx in loci {
            let locus = model.locus(locus_idx);
            let cov = self
                .coverage
                .entry(locus_idx)
                .or_insert_with(|| vec![0u32; locus.exon_union_len()]);
            for &(block_start, block_len) in &blocks {
                locus.add_union_coverage(block_start, block_start + block_len - 1, cov);
            }
        }
        (overlaps_exon, origin)
    }

    /// Accumulate the strand 2×2 table and the R1/R2-transcript-strand template counts for a
    /// read that overlaps exactly one gene locus and at least one exon. Ports the strand and
    /// template logic of Picard `RnaSeqMetricsCollector.acceptRecord`.
    fn tally_strand(&mut self, record: &RikerRecord, locus_idx: u32) {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let locus = model.locus(locus_idx);
        let flags = record.flags();
        let negative_tx = locus.negative_strand;
        let read_one_or_unpaired = !flags.is_segmented() || flags.is_first_segment();
        let negative_read = flags.is_reverse_complemented();
        let agree = negative_read == negative_tx;

        match (read_one_or_unpaired, agree) {
            (true, true) => self.agree_r1 += 1,
            (true, false) => self.disagree_r1 += 1,
            (false, true) => self.agree_r2 += 1,
            (false, false) => self.disagree_r2 += 1,
        }

        // Templates are counted once, from read 1 / unpaired reads only.
        if !read_one_or_unpaired {
            return;
        }
        let start = position_u32(record.alignment_start());
        let end = position_u32(record.alignment_end());
        let (proper, leftmost, rightmost) = if flags.is_segmented() {
            if flags.is_mate_unmapped() {
                (false, 0, 0)
            } else {
                let mate_start = position_u32(record.mate_alignment_start());
                let mate_end = mate_alignment_end(record, mate_start);
                let proper = record.reference_sequence_id() == record.mate_reference_sequence_id()
                    && get_pair_orientation(record) == Some(PairOrientation::Fr);
                (proper, start.min(mate_start), end.max(mate_end))
            }
        } else {
            (true, start, end)
        };

        if proper && encloses_inclusive(locus.start, locus.end, leftmost, rightmost) {
            if agree {
                self.num_r1 += 1;
            } else {
                self.num_r2 += 1;
            }
        } else {
            self.num_unexplained += 1;
        }
    }

    /// Resolve the mate read's reference span (see [`tally_insert_size`]). Returns `None` when the
    /// pair should be skipped: the MC path counts once at the first segment; the MC-less fallback
    /// handles only FR pairs at the forward read with a spec-compliant TLEN.
    fn resolve_mate_span(
        &mut self,
        record: &RikerRecord,
        rec_start: u32,
        mate_start: u32,
        rec_mapped: u32,
    ) -> Option<MateSpan> {
        let flags = record.flags();
        if let Some(mc) = mate_cigar(record) {
            if !flags.is_first_segment() {
                return None; // exact path counts each pair once, at the first segment
            }
            let (blocks, end) = mate_alignment_blocks(&mc, mate_start)?;
            let mapped = blocks.iter().map(|(_, l)| *l).sum();
            Some(MateSpan { blocks, end, mapped, approximate: false })
        } else {
            if !self.mc_missing_warned {
                log::warn!(
                    "rna: some read pairs lack the MC (mate CIGAR) tag; transcript-space insert \
                     size for FR pairs is estimated from the TLEN field (non-FR pairs are \
                     skipped). Run e.g. `samtools fixmate -m` to add MC for exact results. This \
                     warning is shown once."
                );
                self.mc_missing_warned = true;
            }
            // FR pairs only, at the forward (leftmost, not reverse) read so the mate's right end is
            // recoverable from TLEN; a degenerate (zero) TLEN is skipped rather than guessed.
            if flags.is_reverse_complemented()
                || !matches!(get_pair_orientation(record), Some(PairOrientation::Fr))
                || record.template_length() == 0
            {
                return None;
            }
            let mate_end = tlen_fragment_end(rec_start, mate_start, record.template_length());
            if mate_end < mate_start {
                return None;
            }
            Some(MateSpan {
                blocks: vec![(mate_start, mate_end - mate_start + 1)],
                end: mate_end,
                mapped: rec_mapped, // proxy: assume the mate's read length matches this read's
                approximate: true,
            })
        }
    }

    /// Compute and accumulate the transcript-space insert size for a read pair, adapting fgbio
    /// `EstimateRnaSeqInsertSize`.
    ///
    /// The mate's reference span is taken from the `MC` (mate CIGAR) tag when present — the exact
    /// path, used at the first segment for any orientation. When `MC` is absent we fall back to a
    /// **TLEN-based best guess for FR pairs only**, evaluated at the forward (leftmost) read: the
    /// mate's right end is `rec_start + |TLEN| - 1`, so the pair span and both 5' ends are known
    /// even though the mate's internal splicing is not. For that approximate mate we additionally
    /// require both of its ends to fall in an exon, and its exonic footprint to clear the overlap
    /// threshold (a spliced mate still passes because the intron drops out). This needs a
    /// spec-compliant TLEN; a degenerate (zero) TLEN is skipped rather than guessed.
    ///
    /// `overlaps` are the gene loci overlapping the first read (from the coordinate sweep). The
    /// pair is used only when **exactly one** of them encloses the whole pair span and the read
    /// and mate clear the per-transcript exon-overlap threshold, with all transcripts agreeing on
    /// the size. This differs from fgbio (which requires exactly one gene to overlap the pair
    /// span): a non-enclosing bystander gene near the pair no longer disqualifies it.
    fn tally_insert_size(&mut self, record: &RikerRecord, ref_id: usize, overlaps: &[u32]) {
        let flags = record.flags();
        // Gating common to both paths: paired, mate mapped, same contig, not supplementary.
        if !flags.is_segmented()
            || flags.is_mate_unmapped()
            || flags.is_supplementary()
            || record.mate_reference_sequence_id() != Some(ref_id)
        {
            return;
        }

        let rec_start = position_u32(record.alignment_start());
        let rec_end = position_u32(record.alignment_end());
        let mate_start = position_u32(record.mate_alignment_start());
        let rec_blocks = rec_alignment_blocks(record);
        let rec_mapped: u32 = rec_blocks.iter().map(|(_, l)| *l).sum();

        let Some(mate) = self.resolve_mate_span(record, rec_start, mate_start, rec_mapped) else {
            return;
        };

        let interval_start = rec_start.min(mate_start);
        let interval_end = rec_end.max(mate.end);

        let model = self.gene_model.as_ref().expect("gene model loaded");
        // The pair must be enclosed by exactly one gene. An enclosing gene necessarily overlaps
        // the first read, so all candidates are already in `overlaps`.
        let mut enclosing = overlaps.iter().copied().filter(|&l| {
            let locus = model.locus(l);
            encloses_inclusive(locus.start, locus.end, interval_start, interval_end)
        });
        let Some(locus_idx) = enclosing.next() else {
            return;
        };
        if enclosing.next().is_some() {
            return; // ambiguous: more than one gene encloses the pair
        }
        let locus = model.locus(locus_idx);

        let mapped_bases = rec_mapped + mate.mapped;
        if mapped_bases == 0 {
            return;
        }

        // Each qualifying transcript must agree on the insert size, else the pair is skipped.
        let mate_negative = flags.is_mate_reverse_complemented();
        let read_negative = flags.is_reverse_complemented();
        let mut insert_size: Option<u32> = None;
        for tx in &locus.transcripts {
            // For the approximate mate we don't know its internal splicing, so require both ends
            // to anchor in an exon before trusting its (collapsed) exonic footprint.
            if mate.approximate
                && !(point_in_exons(mate_start, &tx.exons) && point_in_exons(mate.end, &tx.exons))
            {
                continue;
            }
            let overlap = bases_overlapping_exons(&rec_blocks, &tx.exons)
                + bases_overlapping_exons(&mate.blocks, &tx.exons);
            if f64::from(overlap) / f64::from(mapped_bases) < self.isize_minimum_overlap {
                continue;
            }
            let rec5 = if read_negative { rec_end } else { rec_start };
            let mate5 = if mate_negative { mate.end } else { mate_start };
            let lower = rec5.min(mate5);
            let upper = rec5.max(mate5);
            let tx_isize: u32 =
                tx.exons.iter().map(|e| overlap_inclusive(lower, upper, e.start, e.end)).sum();
            match insert_size {
                None => insert_size = Some(tx_isize),
                Some(prev) if prev != tx_isize => return, // transcripts disagree → skip
                Some(_) => {}
            }
        }

        let Some(isize) = insert_size else {
            return;
        };
        match get_pair_orientation(record) {
            Some(PairOrientation::Fr) => self.fr.count(u64::from(isize)),
            Some(PairOrientation::Rf) => self.rf.count(u64::from(isize)),
            Some(PairOrientation::Tandem) => self.tandem.count(u64::from(isize)),
            None => {}
        }
    }

    /// Detect the library strandedness, applying the user override or inferring from the
    /// observed R1/R2-transcript-strand counts (≥75% one way ⇒ that orientation).
    fn detect_strand(&self) -> StrandSpec {
        if self.strand != StrandSpec::Auto {
            return self.strand;
        }
        let examined = self.num_r1 + self.num_r2;
        if examined == 0 {
            return StrandSpec::Unstranded;
        }
        let frac_r1 = self.num_r1 as f64 / examined as f64;
        if frac_r1 >= 0.75 {
            StrandSpec::Forward
        } else if frac_r1 <= 0.25 {
            StrandSpec::Reverse
        } else {
            StrandSpec::Unstranded
        }
    }

    /// Build the annotated splice-junction sets from every transcript with ≥2 exons: the set of
    /// whole junctions and the sets of donor/acceptor boundary sites (used to split observed
    /// junctions into known / partial-novel / novel).
    #[allow(clippy::cast_possible_truncation)]
    fn build_known_junctions(&mut self, model: &GeneModel) {
        let mut exons: Vec<(u32, u32)> = Vec::new();
        for locus in model.loci() {
            let ref_id = locus.ref_id as u32;
            for tx in &locus.transcripts {
                if tx.exons.len() < 2 {
                    continue;
                }
                exons.clear();
                exons.extend(tx.exons.iter().map(|e| (e.start, e.end)));
                exons.sort_unstable();
                for win in exons.windows(2) {
                    let (intron_start, intron_end) = (win[0].1 + 1, win[1].0 - 1);
                    if intron_end >= intron_start {
                        self.known_introns.insert((ref_id, intron_start, intron_end));
                        self.known_intron_starts.insert((ref_id, intron_start));
                        self.known_intron_ends.insert((ref_id, intron_end));
                    }
                }
            }
        }
    }

    /// Classify an observed intron against the annotated junction sets.
    fn classify_intron(&self, key: JunctionKey) -> JunctionClass {
        let (ref_id, start, end) = key;
        if self.known_introns.contains(&key) {
            JunctionClass::Known
        } else if self.known_intron_starts.contains(&(ref_id, start))
            || self.known_intron_ends.contains(&(ref_id, end))
        {
            JunctionClass::Partial
        } else {
            JunctionClass::Novel
        }
    }

    /// Tally splice junctions for a primary, mapped, non-supplementary read: count each observed
    /// intron toward the observation-level totals and record it for the distinct-junction totals.
    #[allow(clippy::cast_possible_truncation)]
    fn tally_junctions(&mut self, record: &RikerRecord, ref_id: usize) {
        let introns = rec_introns(record, self.junction_min_intron);
        if introns.is_empty() {
            return;
        }
        self.spliced_reads += 1;
        let ref_id = ref_id as u32;
        for (start, end) in introns {
            let key = (ref_id, start, end);
            match self.classify_intron(key) {
                JunctionClass::Known => self.known_splice_obs += 1,
                JunctionClass::Partial => self.partial_splice_obs += 1,
                JunctionClass::Novel => self.novel_splice_obs += 1,
            }
            self.observed_introns.insert(key);
        }
    }

    /// Aggregate the per-locus assigned-read counts by gene biotype and write the biotype file
    /// (one row per biotype, most-assigned first). `frac_reads` is over `assigned_reads`.
    #[allow(clippy::cast_possible_truncation)]
    fn finish_biotypes(&self) -> Result<()> {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let mut by_biotype: FxHashMap<String, u64> = FxHashMap::default();
        for (locus_idx, &count) in self.assigned_counts.iter().enumerate() {
            if count > 0 {
                let label = biotype_label(model.locus(locus_idx as u32).biotype.as_ref());
                *by_biotype.entry(label).or_default() += u64::from(count);
            }
        }
        let assigned = self.assigned_reads as f64;
        let mut rows: Vec<RnaBiotypeMetric> = by_biotype
            .into_iter()
            .map(|(biotype, reads)| RnaBiotypeMetric {
                sample: self.sample.clone(),
                biotype,
                reads,
                frac_reads: if assigned > 0.0 { reads as f64 / assigned } else { 0.0 },
            })
            .collect();
        rows.sort_by(|a, b| b.reads.cmp(&a.reads).then_with(|| a.biotype.cmp(&b.biotype)));
        write_tsv(&self.biotype_path, &rows)?;
        self.plot_biotypes(&rows)
    }

    /// Bar chart of assigned reads per biotype. Read counts span several orders of magnitude
    /// (protein_coding typically dwarfs everything else in poly-A data), so the y-axis is
    /// log-scaled to keep the minor biotypes visible. Only the most abundant biotypes are shown,
    /// and each bar is annotated with its share of assigned reads.
    #[allow(clippy::cast_precision_loss)]
    fn plot_biotypes(&self, rows: &[RnaBiotypeMetric]) -> Result<()> {
        use kuva::plot::BarPlot;
        use kuva::render::annotations::TextAnnotation;
        use kuva::render::layout::{Layout, TickFormat};
        use kuva::render::plots::Plot;

        /// Cap on the number of biotypes drawn (rows are sorted most-abundant first).
        const MAX_BARS: usize = 15;

        let top: Vec<&RnaBiotypeMetric> =
            rows.iter().filter(|r| r.reads > 0).take(MAX_BARS).collect();
        if top.is_empty() {
            return Ok(());
        }
        let bars: Vec<(String, f64)> =
            top.iter().map(|r| (prettify_biotype(&r.biotype), r.reads as f64)).collect();

        let plots = vec![Plot::Bar(BarPlot::new().with_bars(bars).with_color(FG_NAVY))];
        let mut layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("Assigned Reads by Biotype of {}", self.sample))
            .with_x_label("Biotype")
            .with_y_label("Assigned reads (log scale)")
            .with_x_tick_rotate(45.0)
            .with_log_y()
            // Force every y tick to scientific notation so the labels are consistent (the default
            // mixes "1000" and "1e7" across a log axis spanning many decades).
            .with_y_tick_format(TickFormat::Sci);

        // Annotate each bar's top with its share of assigned reads. Categorical bar `i` (0-based)
        // is centered at data-x `i + 1`; lift the label just above the bar on the log y-axis.
        for (i, row) in top.iter().enumerate() {
            let pct = row.frac_reads * 100.0;
            let label =
                if pct > 0.0 && pct < 0.01 { "<0.01%".to_string() } else { format!("{pct:.2}%") };
            let x = i as f64 + 1.0;
            let y = row.reads as f64 * 1.35;
            layout = layout.with_annotation(TextAnnotation::new(label, x, y).with_font_size(9));
        }
        write_plot_pdf(plots, layout, &self.biotype_plot_path)
    }

    /// Per-gene expression-distribution histogram. For every gene with at least one uniquely-
    /// assigned read compute RPK — reads per kilobase of the gene's exon-union length — and plot
    /// the distribution, overlaying all genes against the protein-coding subset. RPK normalizes
    /// raw read counts for gene length but not for library size, so it shows each gene's relative
    /// expression and the dynamic range of the sample. Bins are equal-width in log10(RPK) so the
    /// heavy right tail spreads out, and the x-axis is log-scaled but labelled in raw RPK.
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    fn finish_gene_expression(&self) -> Result<()> {
        /// Number of equal-width (in log10 space) histogram bins over the observed RPK range.
        const N_BINS: usize = 100;

        /// Per-gene aggregate over its loci: total reads, exon-union length, any-locus-coding.
        #[derive(Default)]
        struct GeneAgg {
            reads: u64,
            length: u64,
            coding: bool,
        }

        let model = self.gene_model.as_ref().expect("gene model loaded");
        // Aggregate by gene name over all loci; a gene split across disjoint loci is scored once.
        let mut per_gene: FxHashMap<&str, GeneAgg> = FxHashMap::default();
        for (locus_idx, &count) in self.assigned_counts.iter().enumerate() {
            let locus = model.locus(locus_idx as u32);
            let agg = per_gene.entry(locus.gene_name.as_str()).or_default();
            agg.reads += u64::from(count);
            agg.length += locus.exon_union_len() as u64;
            agg.coding |= matches!(locus.biotype, Some(Biotype::ProteinCoding));
        }

        // log10(RPK) for every expressed gene (>0 reads), all genes and the protein-coding subset.
        let log_rpk = |a: &GeneAgg| (a.reads as f64 / (a.length as f64 / 1000.0)).log10();
        let expressed = |a: &&GeneAgg| a.reads > 0 && a.length > 0;
        let log_all: Vec<f64> = per_gene.values().filter(expressed).map(&log_rpk).collect();
        let log_coding: Vec<f64> =
            per_gene.values().filter(|a| a.coding).filter(expressed).map(&log_rpk).collect();
        if log_all.is_empty() {
            return Ok(());
        }

        // Genes with zero assigned reads, split by coding/non-coding (for the annotation).
        let coding_zero = per_gene.values().filter(|a| a.coding && a.reads == 0).count();
        let noncoding_zero = per_gene.values().filter(|a| !a.coding && a.reads == 0).count();

        // Bin both series on shared edges (all-genes range covers the coding subset), in log10
        // space for equal-width-on-a-log-axis bars, then convert edges back to raw RPK so the
        // histogram can be drawn against a log-scaled (but RPK-labelled) x-axis.
        let (log_edges, counts_all) = histogram_bins(&log_all, N_BINS);
        let counts_coding = histogram_counts(&log_coding, &log_edges);
        let edges: Vec<f64> = log_edges.iter().map(|e| 10f64.powf(*e)).collect();
        self.plot_gene_expression(edges, counts_all, counts_coding, coding_zero, noncoding_zero)
    }

    /// Overlaid RPK histograms (all genes vs protein-coding) on a log-scaled x-axis, with a legend
    /// inset and a zero-read-gene annotation.
    fn plot_gene_expression(
        &self,
        edges: Vec<f64>,
        counts_all: Vec<f64>,
        counts_coding: Vec<f64>,
        coding_zero: usize,
        noncoding_zero: usize,
    ) -> Result<()> {
        use kuva::plot::{Histogram, LegendEntry, LegendPosition, LegendShape};
        use kuva::render::layout::Layout;
        use kuva::render::plots::Plot;

        // 8-digit hex (#RRGGBBAA) makes the fills semi-transparent so the overlap shows through.
        let plots = vec![
            Plot::Histogram(
                Histogram::from_bins(edges.clone(), counts_all).with_color(format!("{FG_BLUE}80")),
            ),
            Plot::Histogram(
                Histogram::from_bins(edges, counts_coding).with_color(format!("{FG_GREEN}80")),
            ),
        ];
        // Histogram legends aren't auto-collected by `auto_from_plots`, so supply explicit entries
        // (solid swatches for crisp colour) — this also enables the legend.
        let legend = |label: &str, color: &str| LegendEntry {
            label: label.to_string(),
            color: color.to_string(),
            shape: LegendShape::Rect,
            dasharray: None,
        };
        let annotation = vec![
            format!("{} protein coding genes with 0 reads", thousands(coding_zero as u64)),
            format!("{} non-coding genes with 0 reads", thousands(noncoding_zero as u64)),
        ];
        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("Gene Expression Distribution of {}", self.sample))
            .with_x_label("Reads per kilobase (RPK)")
            .with_y_label("Genes")
            .with_log_x()
            .with_legend_entries(vec![
                legend("All genes", FG_BLUE),
                legend("Protein-coding", FG_GREEN),
            ])
            .with_legend_position(LegendPosition::InsideTopRight)
            .with_stats_box_at(LegendPosition::InsideTopLeft, annotation)
            // Plain text annotation — no border or background box around it.
            .with_stats_box_border(false);
        write_plot_pdf(plots, layout, &self.gene_expression_plot_path)
    }

    /// Count genes (by name) with at least `genes_detected_min_reads` uniquely-assigned reads.
    /// Per-locus counts are aggregated by gene name so a gene spanning disjoint loci counts once.
    #[allow(clippy::cast_possible_truncation)]
    fn count_genes_detected(&self) -> u64 {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let mut per_gene: FxHashMap<&str, u64> = FxHashMap::default();
        for (locus_idx, &count) in self.assigned_counts.iter().enumerate() {
            if count > 0 {
                let name = model.locus(locus_idx as u32).gene_name.as_str();
                *per_gene.entry(name).or_default() += u64::from(count);
            }
        }
        let threshold = u64::from(self.genes_detected_min_reads);
        per_gene.values().filter(|&&n| n >= threshold).count() as u64
    }

    /// Compute the Transcript Integrity Number (median / mean / stddev / count). One representative
    /// transcript per gene — the most-expressed one at least `minimum_length` long — contributes,
    /// and only if its **mean** per-base coverage exceeds `tin_min_coverage`. The mean (rather than
    /// peak) gate excludes lowly-covered, noisy transcripts while keeping deeply but unevenly
    /// covered ones, since a narrow pileup barely moves the whole-transcript mean. TIN is computed
    /// over the full per-transcript coverage reconstructed from the locus's exon-union array.
    fn finish_tin(&self) -> (f64, f64, f64, u64) {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let min_mean = f64::from(self.tin_min_coverage);
        let mut tins: Vec<f64> = Vec::new();
        for (&locus_idx, union_cov) in &self.coverage {
            let locus = model.locus(locus_idx);
            // Best-expressed eligible transcript for this gene (one per locus).
            let mut best: Option<(f64, Vec<u32>)> = None;
            for tx in &locus.transcripts {
                if tx.transcribed_length() < self.minimum_length {
                    continue;
                }
                let cov = locus.transcript_coverage(tx, union_cov);
                let mean = mean_u32(&cov);
                if best.as_ref().is_none_or(|(m, _)| mean > *m) {
                    best = Some((mean, cov));
                }
            }
            if let Some((mean, cov)) = best
                && mean > min_mean
                && let Some(tin) = transcript_tin(&cov)
            {
                tins.push(tin);
            }
        }
        if tins.is_empty() {
            return (0.0, 0.0, 0.0, 0);
        }
        let count = tins.len();
        let mean = tins.iter().sum::<f64>() / count as f64;
        let variance = tins.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / count as f64;
        let median = median(&mut tins);
        (median, mean, variance.sqrt(), count as u64)
    }

    /// Compute the coverage/bias metrics over the most highly expressed transcripts and
    /// write the normalized-coverage-by-position table and chart.
    fn finish_coverage(&self) -> Result<CoverageMetrics> {
        let model = self.gene_model.as_ref().expect("gene model loaded");
        let min_len = self.minimum_length.max(2 * self.end_bias_bases);

        // Best (highest mean coverage) eligible transcript per locus. Per-transcript coverage is
        // reconstructed on demand from the locus's exon-union coverage.
        let mut best_per_locus: FxHashMap<u32, (f64, Vec<u32>, bool)> = FxHashMap::default();
        for (&locus_idx, union_cov) in &self.coverage {
            let locus = model.locus(locus_idx);
            for tx in &locus.transcripts {
                if tx.transcribed_length() < min_len {
                    continue;
                }
                let cov = locus.transcript_coverage(tx, union_cov);
                let mean = mean_u32(&cov);
                if mean < 1.0 {
                    continue;
                }
                let better = best_per_locus.get(&locus_idx).is_none_or(|(m, _, _)| mean > *m);
                if better {
                    best_per_locus.insert(locus_idx, (mean, cov, locus.negative_strand));
                }
            }
        }

        // Keep the top 1000 by mean coverage (fixes Picard's off-by-one top-1001).
        let mut means: Vec<f64> = best_per_locus.values().map(|(m, _, _)| *m).collect();
        let threshold = if means.len() <= 1000 {
            0.0
        } else {
            means.sort_unstable_by(|a, b| b.total_cmp(a));
            means[999]
        };

        let mut cvs = Vec::new();
        let mut five = Vec::new();
        let mut three = Vec::new();
        let mut five_to_three = Vec::new();
        let mut norm_sum = [0.0f64; 101];
        let mut n_tx = 0u32;

        for (mean, cov, negative) in best_per_locus.values() {
            if *mean < threshold {
                continue;
            }
            // Orient 5'→3': reverse the genomic-order coverage for negative-strand transcripts.
            let oriented: Vec<u32> =
                if *negative { cov.iter().rev().copied().collect() } else { cov.clone() };
            let m = *mean;
            cvs.push(stddev_u32(&oriented, m) / m);

            let end_bias = self.end_bias_bases as usize;
            let five_cov = mean_u32(&oriented[..end_bias]);
            let three_cov = mean_u32(&oriented[oriented.len() - end_bias..]);
            five.push(five_cov / m);
            three.push(three_cov / m);
            five_to_three.push(if three_cov > 0.0 { five_cov / three_cov } else { 0.0 });

            let last = oriented.len() - 1;
            for (p, slot) in norm_sum.iter_mut().enumerate() {
                let pos = (p * last) / 100;
                *slot += f64::from(oriented[pos]) / m;
            }
            n_tx += 1;
        }

        let entries: Vec<RnaCoverageEntry> = (0..=100u32)
            .map(|p| RnaCoverageEntry {
                normalized_position: p,
                normalized_coverage: if n_tx > 0 {
                    norm_sum[p as usize] / f64::from(n_tx)
                } else {
                    0.0
                },
            })
            .collect();
        self.plot_coverage(&entries)?;

        Ok(CoverageMetrics {
            cv: median(&mut cvs),
            five_prime: median(&mut five),
            three_prime: median(&mut three),
            five_to_three: median(&mut five_to_three),
        })
    }

    /// Write the normalized-coverage-by-position line chart.
    fn plot_coverage(&self, entries: &[RnaCoverageEntry]) -> Result<()> {
        use kuva::plot::LinePlot;
        use kuva::render::layout::Layout;
        use kuva::render::plots::Plot;

        let xy: Vec<(f64, f64)> = entries
            .iter()
            .map(|e| (f64::from(e.normalized_position), e.normalized_coverage))
            .collect();
        let plots = vec![Plot::Line(
            LinePlot::new().with_data(xy).with_color(FG_BLUE).with_fill().with_fill_opacity(0.3),
        )];
        let layout = Layout::auto_from_plots(&plots)
            .with_width(PLOT_WIDTH)
            .with_height(PLOT_HEIGHT)
            .with_title(format!("Transcript Coverage of {}", self.sample))
            .with_x_label("Normalized Position (5' → 3')")
            .with_y_label("Normalized Coverage")
            .with_x_axis_min(0.0)
            .with_x_axis_max(100.0);
        write_plot_pdf(plots, layout, &self.coverage_plot_path)
    }

    /// Assemble and write the transcript-space insert-size metrics and histogram.
    fn finish_insert_size(&self) -> Result<()> {
        let orientations: [(&str, &Counter<u64>); 3] =
            [("FR", &self.fr), ("RF", &self.rf), ("TANDEM", &self.tandem)];

        let metrics: Vec<RnaInsertSizeMetric> = orientations
            .iter()
            .map(|(name, hist)| {
                let (median, mad) = hist.median_and_mad();
                let trim_max = median + self.isize_deviations * mad;
                let (mean, stddev) = hist.trimmed_mean_and_stddev(trim_max);
                RnaInsertSizeMetric {
                    sample: self.sample.clone(),
                    pair_orientation: (*name).to_string(),
                    read_pairs: hist.total(),
                    mean_insert_size: mean,
                    standard_deviation: stddev,
                    median_insert_size: median,
                    median_absolute_deviation: mad,
                    min_insert_size: hist.min().unwrap_or(0),
                    max_insert_size: hist.max().unwrap_or(0),
                }
            })
            .collect();
        write_tsv(&self.isize_metrics_path, &metrics)?;

        let all_keys: std::collections::BTreeSet<u64> =
            self.fr.keys().chain(self.rf.keys()).chain(self.tandem.keys()).copied().collect();
        let hist_entries: Vec<RnaInsertSizeHistogramEntry> = all_keys
            .into_iter()
            .map(|k| RnaInsertSizeHistogramEntry {
                sample: self.sample.clone(),
                insert_size: k,
                fr_count: self.fr.count_of(&k),
                rf_count: self.rf.count_of(&k),
                tandem_count: self.tandem.count_of(&k),
            })
            .collect();
        write_tsv(&self.isize_histogram_path, &hist_entries)?;

        let series = [
            InsertSizeSeries { name: "FR", color: FG_BLUE, histogram: &self.fr },
            InsertSizeSeries { name: "RF", color: FG_GREEN, histogram: &self.rf },
            InsertSizeSeries { name: "TANDEM", color: FG_TEAL, histogram: &self.tandem },
        ];
        write_insert_size_histogram_pdf(
            &self.isize_plot_path,
            &format!("RNA Insert Size of {}", self.sample),
            &series,
            0.05,
            self.isize_deviations,
        )
    }

    /// Returns true if the read's fragment overlaps ribosomal intervals by at least
    /// `rrna_fragment_fraction`, tallying the read's aligned bases as ribosomal if so.
    /// Ports Picard's fragment-overlap rRNA test with the union-overlap and `MC`-fragment-end
    /// fixes.
    fn classify_ribosomal(&mut self, record: &RikerRecord, ref_id: usize) -> bool {
        let ribosomal = self.ribosomal.as_ref().expect("ribosomal source present");
        let flags = record.flags();
        let start = position_u32(record.alignment_start());
        let end = position_u32(record.alignment_end());

        // Fragment span (1-based inclusive): full template for a proper pair, else the read.
        let (frag_start, frag_end) = if !flags.is_segmented() {
            (start, end)
        } else if flags.is_mate_unmapped() || record.mate_reference_sequence_id() != Some(ref_id) {
            return false; // chimeric / half-mapped: no fragment interval (Picard)
        } else {
            let mate_start = position_u32(record.mate_alignment_start());
            let mate_end = mate_alignment_end(record, mate_start);
            (start.min(mate_start), end.max(mate_end))
        };

        // 0-based half-open query start; `saturating_sub` guards against a degenerate 0 start.
        let frag_start0 = frag_start.saturating_sub(1);
        let frag_len = u64::from(frag_end.saturating_sub(frag_start) + 1);
        // Sum overlaps with the (disjoint, merged) ribosomal intervals: a true union.
        let mut overlap: u64 = 0;
        for &(rs, re) in ribosomal.get_overlaps(ref_id, frag_start0, frag_end) {
            let lo = frag_start0.max(rs);
            let hi = frag_end.min(re);
            if lo < hi {
                overlap += u64::from(hi - lo);
            }
        }

        if frag_len > 0 && overlap as f64 / frag_len as f64 >= self.rrna_fragment_fraction {
            let aligned: u64 =
                rec_alignment_blocks(record).iter().map(|(_, l)| u64::from(*l)).sum();
            self.aligned_bases += aligned;
            self.ribosomal_bases += aligned;
            true
        } else {
            false
        }
    }
}

impl Collector for RnaCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        // The locus sweep requires coordinate-sorted input; refuse anything else up front.
        anyhow::ensure!(
            is_coordinate_sorted(header),
            "rna requires a coordinate-sorted BAM/CRAM (@HD SO:coordinate); sort with `samtools sort`"
        );

        let dict = SequenceDictionary::from(header);
        self.sample = derive_sample(&self.input, header);

        let model = GeneModel::from_path(&self.gene_model_path, &dict)?;
        self.ribosomal = self.build_ribosomal(&model, &dict)?;
        self.assigned_counts = vec![0u32; model.loci().len()];
        self.build_known_junctions(&model);
        self.gene_model = Some(model);

        for name in &self.ignore_chrom {
            match dict.get_by_name(name) {
                Some(meta) => {
                    self.ignored_ref_ids.insert(meta.index());
                }
                None => log::warn!("rna: --ignore-chrom contig '{name}' not in BAM header"),
            }
        }
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();
        if flags.is_qc_fail() {
            return Ok(());
        }

        // `bases`/`total_reads` count every primary (non-secondary) read, mapped or not
        // (Picard PF_BASES).
        if !flags.is_secondary() {
            self.bases += record.sequence_len() as u64;
            self.total_reads += 1;
        }

        let ref_id = record.reference_sequence_id();
        // Ignore-sequence reads count only as `ignored_reads`.
        if !flags.is_unmapped()
            && !flags.is_secondary()
            && ref_id.is_some_and(|r| self.ignored_ref_ids.contains(&r))
        {
            self.ignored_reads += 1;
            return Ok(());
        }
        if flags.is_secondary() {
            return Ok(());
        }
        if flags.is_unmapped() {
            self.unmapped_reads += 1;
            return Ok(());
        }
        let ref_id = ref_id.expect("mapped read has a reference id");

        // Duplicate rate is observed over primary, mapped, non-supplementary reads.
        if !flags.is_supplementary() {
            self.mapped_reads += 1;
            if flags.is_duplicate() {
                self.duplicate_reads += 1;
            }
        }
        if self.exclude_duplicates && flags.is_duplicate() {
            return Ok(());
        }
        // MAPQ gate (default 0 ⇒ no-op). A missing/255 MAPQ is treated as passing.
        if record.mapping_quality().is_some_and(|mq| mq.get() < self.min_mapq) {
            return Ok(());
        }

        // Ribosomal fragment classification (early-return if the fragment is ribosomal).
        if self.ribosomal.is_some() && self.classify_ribosomal(record, ref_id) {
            if !flags.is_supplementary() {
                self.ribosomal_reads += 1;
            }
            return Ok(());
        }

        // Gene loci overlapping this read, from the coordinate sweep. Shared by classification,
        // strand tallying, and insert size.
        let start = position_u32(record.alignment_start());
        let end = position_u32(record.alignment_end());
        let mut loci = Vec::new();
        self.overlapping_loci(ref_id, start, end, &mut loci)?;

        // Base-level functional classification + coverage, plus the read's genomic origin.
        let (overlaps_exon, origin) = self.classify_bases(record, &loci);

        // Read-level genomic-origin accounting over primary, mapped, non-supplementary reads.
        if !flags.is_supplementary() {
            match origin {
                ReadOrigin::Exonic(locus) => {
                    self.exonic_reads += 1;
                    self.assigned_reads += 1;
                    self.assigned_counts[locus as usize] += 1;
                }
                ReadOrigin::Ambiguous => {
                    self.exonic_reads += 1;
                    self.ambiguous_reads += 1;
                }
                ReadOrigin::Intronic => self.intronic_reads += 1,
                ReadOrigin::Intergenic => self.intergenic_reads += 1,
            }
            self.tally_junctions(record, ref_id);
        }

        // Strand + template-strand metrics for single-locus exon-overlapping reads.
        if !flags.is_supplementary() && overlaps_exon && loci.len() == 1 {
            self.tally_strand(record, loci[0]);
        }

        // Transcript-space insert size (own filters; counts each pair once).
        self.tally_insert_size(record, ref_id, &loci);
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn finish(&mut self) -> Result<()> {
        let coverage = self.finish_coverage()?;
        let (median_tin, mean_tin, stddev_tin, tin_transcripts) = self.finish_tin();
        self.finish_insert_size()?;

        // `detect_strand` never returns `Auto`; unstranded (and the unreachable `Auto`) → zero.
        let detected = self.detect_strand();
        let (correct, incorrect) = match detected {
            StrandSpec::Forward => {
                (self.agree_r1 + self.disagree_r2, self.disagree_r1 + self.agree_r2)
            }
            StrandSpec::Reverse => {
                (self.disagree_r1 + self.agree_r2, self.agree_r1 + self.disagree_r2)
            }
            StrandSpec::Unstranded | StrandSpec::Auto => (0, 0),
        };

        let aligned = self.aligned_bases as f64;
        let frac = |n: u64| if aligned > 0.0 { n as f64 / aligned } else { 0.0 };
        let frac_coding = frac(self.coding_bases);
        let frac_utr = frac(self.utr_bases);
        let strand_total = correct + incorrect;
        let r_total = self.num_r1 + self.num_r2;
        let has_ribosomal = self.ribosomal.is_some();

        // Read-level fractions are over mapped (primary, non-supplementary) reads.
        let mapped = self.mapped_reads as f64;
        let frac_reads = |n: u64| if mapped > 0.0 { n as f64 / mapped } else { 0.0 };
        let genes_detected = self.count_genes_detected();

        // Distinct-junction totals: classify each distinct observed intron once.
        let (mut known_juncs, mut partial_juncs, mut novel_juncs) = (0u64, 0u64, 0u64);
        for &key in &self.observed_introns {
            match self.classify_intron(key) {
                JunctionClass::Known => known_juncs += 1,
                JunctionClass::Partial => partial_juncs += 1,
                JunctionClass::Novel => novel_juncs += 1,
            }
        }
        let total_juncs = known_juncs + partial_juncs + novel_juncs;
        let frac_known_juncs =
            if total_juncs > 0 { known_juncs as f64 / total_juncs as f64 } else { 0.0 };

        let metric = RnaSeqMetric {
            sample: self.sample.clone(),

            total_reads: self.total_reads,
            mapped_reads: self.mapped_reads,
            unmapped_reads: self.unmapped_reads,
            ignored_reads: self.ignored_reads,
            duplicate_rate: if self.mapped_reads > 0 {
                self.duplicate_reads as f64 / self.mapped_reads as f64
            } else {
                0.0
            },

            ribosomal_reads: has_ribosomal.then_some(self.ribosomal_reads),
            exonic_reads: self.exonic_reads,
            intronic_reads: self.intronic_reads,
            intergenic_reads: self.intergenic_reads,
            ambiguous_reads: self.ambiguous_reads,
            assigned_reads: self.assigned_reads,
            frac_ribosomal_reads: has_ribosomal.then(|| frac_reads(self.ribosomal_reads)),
            frac_exonic_reads: frac_reads(self.exonic_reads),
            frac_intronic_reads: frac_reads(self.intronic_reads),
            frac_intergenic_reads: frac_reads(self.intergenic_reads),
            frac_ambiguous_reads: frac_reads(self.ambiguous_reads),

            genes_detected,

            spliced_reads: self.spliced_reads,
            frac_spliced_reads: frac_reads(self.spliced_reads),
            known_splice_obs: self.known_splice_obs,
            partial_splice_obs: self.partial_splice_obs,
            novel_splice_obs: self.novel_splice_obs,
            known_juncs,
            partial_juncs,
            novel_juncs,
            frac_known_juncs,

            bases: self.bases,
            aligned_bases: self.aligned_bases,
            ribosomal_bases: has_ribosomal.then_some(self.ribosomal_bases),
            coding_bases: self.coding_bases,
            utr_bases: self.utr_bases,
            intronic_bases: self.intronic_bases,
            intergenic_bases: self.intergenic_bases,
            frac_ribosomal_bases: has_ribosomal.then(|| frac(self.ribosomal_bases)),
            frac_coding_bases: frac_coding,
            frac_utr_bases: frac_utr,
            frac_intronic_bases: frac(self.intronic_bases),
            frac_intergenic_bases: frac(self.intergenic_bases),
            frac_mrna_bases: frac_coding + frac_utr,
            // Picard `PCT_USABLE_BASES` is intentionally relative to all sequenced bases.
            frac_usable_bases: if self.bases > 0 {
                (self.coding_bases + self.utr_bases) as f64 / self.bases as f64
            } else {
                0.0
            },

            detected_strand: detected.to_string(),
            correct_strand_reads: correct,
            frac_correct_strand_reads: if strand_total > 0 {
                correct as f64 / strand_total as f64
            } else {
                0.0
            },
            r1_tx_strand_reads: self.num_r1,
            r2_tx_strand_reads: self.num_r2,
            unexplained_reads: self.num_unexplained,
            frac_r1_tx_strand_reads: if r_total > 0 {
                self.num_r1 as f64 / r_total as f64
            } else {
                0.0
            },
            frac_r2_tx_strand_reads: if r_total > 0 {
                self.num_r2 as f64 / r_total as f64
            } else {
                0.0
            },

            median_cv_coverage: coverage.cv,
            median_5prime_bias: coverage.five_prime,
            median_3prime_bias: coverage.three_prime,
            median_5prime_to_3prime_bias: coverage.five_to_three,
            median_tin,
            mean_tin,
            stddev_tin,
            tin_transcripts,
        };
        write_tsv(&self.metrics_path, &[metric])?;
        self.finish_biotypes()?;
        self.finish_gene_expression()
    }

    fn name(&self) -> &'static str {
        "rna"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        // The MC tag is needed for transcript-space insert size; nothing else is expensive.
        RikerRecordRequirements::NONE.with_aux_tag(MC_TAG)
    }

    /// Relative per-record worker cost (see [`Collector::cost_hint`]); ~150
    /// from a samply worker-compute measurement at `--threads 2` on two real
    /// ENCODE RNA-seq BAMs. rna is heavy — ~2-2.85x `basic` per record (locus
    /// sweep + per-base classification + coverage). The exact multiplier is
    /// data-dependent: RNA read characteristics shift relative tool costs, so
    /// the WGS-calibrated `basic`/`alignment` ratio does not hold on RNA data
    /// and the two anchors bracket a ~90-216 range; 150 is the center.
    fn cost_hint(&self) -> u32 {
        150
    }
}

// ─── Metric structs ──────────────────────────────────────────────────────────

/// RNA-seq QC metrics (one row): read accounting, genomic origin, gene detection, splice
/// junctions, base composition, strand specificity, and coverage uniformity/integrity.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct RnaSeqMetric {
    /// Sample name from the BAM read group SM tag, or the filename.
    pub sample: String,

    // ── Read accounting ──
    /// Every QC-passing read except secondary alignments — i.e. mapped, unmapped, and supplementary
    /// reads. The denominator for `bases`. Because supplementary alignments are counted here but not
    /// in `mapped_reads`, expect `total_reads` ≥ `mapped_reads` + `unmapped_reads`.
    pub total_reads: u64,
    /// Mapped, non-secondary, non-supplementary reads — the "one record per read" count that is the
    /// denominator for `duplicate_rate` and all `frac_*_reads` genomic-origin fractions.
    pub mapped_reads: u64,
    /// Unmapped, non-secondary reads.
    pub unmapped_reads: u64,
    /// Reads on an --ignore-chrom contig. Counted in `bases`/`total_reads` but excluded from every
    /// other metric (Picard's IGNORE_SEQUENCE behavior).
    pub ignored_reads: u64,
    /// `duplicate-flagged mapped_reads / mapped_reads`. Always reported as an observed rate;
    /// duplicates are still *included* in every other metric unless `--exclude-duplicates` is set.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub duplicate_rate: f64,

    // ── Read genomic origin (read-level) ──
    /// Reads whose fragment overlaps a ribosomal interval by at least --rrna-fragment-fraction;
    /// these are tallied here and removed before the exonic/intronic/intergenic classification.
    /// Blank if no ribosomal source (gene-model biotype or --ribosomal-intervals) is available.
    pub ribosomal_reads: Option<u64>,
    /// Reads with at least one aligned base in an exon of any gene. Equals `assigned_reads`
    /// (exonic to exactly one gene) + `ambiguous_reads` (exonic to more than one gene).
    pub exonic_reads: u64,
    /// Reads contained within a gene's span but touching no exon (i.e. wholly intronic).
    pub intronic_reads: u64,
    /// Reads overlapping no gene at all.
    pub intergenic_reads: u64,
    /// Reads with exon overlap in more than one gene. A subset of `exonic_reads`; these are *not*
    /// counted in `assigned_reads` and do not contribute to `genes_detected` or the biotype file.
    pub ambiguous_reads: u64,
    /// Reads with exon overlap in exactly one gene — the "usable for quantification" set that feeds
    /// `genes_detected` and the per-biotype counts.
    pub assigned_reads: u64,
    /// `ribosomal_reads / mapped_reads`. Blank if no ribosomal source.
    #[serde(serialize_with = "serialize_opt_f64_5dp")]
    pub frac_ribosomal_reads: Option<f64>,
    /// `exonic_reads / mapped_reads`. Under default filters,
    /// ribosomal + exonic + intronic + intergenic = mapped_reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_exonic_reads: f64,
    /// `intronic_reads / mapped_reads`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_intronic_reads: f64,
    /// `intergenic_reads / mapped_reads`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_intergenic_reads: f64,
    /// `ambiguous_reads / mapped_reads` (the ambiguous subset of `frac_exonic_reads`).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_ambiguous_reads: f64,

    // ── Gene detection ──
    /// Number of distinct genes (by name, so a gene at disjoint loci counts once) with at least
    /// --genes-detected-min-reads uniquely-assigned reads. Scales with sequencing depth, so compare
    /// only across similarly-sequenced samples.
    pub genes_detected: u64,

    // ── Splice junctions ──
    /// Reads with at least one splice — a CIGAR N (skip) gap of at least --junction-min-intron bp.
    pub spliced_reads: u64,
    /// `spliced_reads / mapped_reads`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_spliced_reads: f64,
    /// Splice-junction *observations* (counted with read multiplicity: a junction seen in N reads
    /// adds N) whose intron matches an annotated junction (both splice sites annotated).
    pub known_splice_obs: u64,
    /// Splice-junction observations (with read multiplicity) where exactly one splice site is
    /// annotated (partially novel).
    pub partial_splice_obs: u64,
    /// Splice-junction observations (with read multiplicity) where neither splice site is annotated.
    pub novel_splice_obs: u64,
    /// Distinct annotated junctions observed (each unique intron counted once, regardless of how
    /// many reads spanned it) — contrast the `*_splice_obs` counts, which include multiplicity.
    pub known_juncs: u64,
    /// Distinct observed junctions with exactly one annotated splice site.
    pub partial_juncs: u64,
    /// Distinct observed junctions with neither splice site annotated.
    pub novel_juncs: u64,
    /// `known_juncs / (known + partial + novel distinct junctions)` — the annotated fraction of the
    /// junctions seen.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_known_juncs: f64,

    // ── Base composition (base-level) ──
    /// Total bases in QC-passing, non-secondary reads, including unmapped reads and soft-clips
    /// (Picard PF_BASES). A base-level total — distinct from the read-level `total_reads`.
    pub bases: u64,
    /// Reference-aligned bases (CIGAR M/=/X) of mapped reads; excludes soft-clips and insertions
    /// (Picard PF_ALIGNED_BASES). The denominator for the `frac_*_bases` composition fractions.
    pub aligned_bases: u64,
    /// Aligned bases in fragments classified ribosomal. Blank if no ribosomal source. Note the
    /// base-level composition is independent of the read-level origin counts above: it tallies
    /// *bases*, and a single read can contribute bases to several categories.
    pub ribosomal_bases: Option<u64>,
    /// Aligned bases falling in coding (CDS) exon sequence.
    pub coding_bases: u64,
    /// Aligned bases falling in UTR (non-coding exon) sequence.
    pub utr_bases: u64,
    /// Aligned bases falling within a gene span but outside any exon (intronic).
    pub intronic_bases: u64,
    /// Aligned bases overlapping no gene.
    pub intergenic_bases: u64,
    /// `ribosomal_bases / aligned_bases`. Blank if no ribosomal source.
    #[serde(serialize_with = "serialize_opt_f64_5dp")]
    pub frac_ribosomal_bases: Option<f64>,
    /// `coding_bases / aligned_bases`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_coding_bases: f64,
    /// `utr_bases / aligned_bases`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_utr_bases: f64,
    /// `intronic_bases / aligned_bases`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_intronic_bases: f64,
    /// `intergenic_bases / aligned_bases`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_intergenic_bases: f64,
    /// `frac_coding_bases + frac_utr_bases` — the mRNA (exonic) fraction of aligned bases.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_mrna_bases: f64,
    /// `(coding_bases + utr_bases) / bases`. Note the denominator is **all sequenced bases**
    /// (`bases`), not `aligned_bases` like the other `frac_*_bases` — this is intentional, matching
    /// Picard's PCT_USABLE_BASES (the usable-mRNA fraction of the whole sequencing run).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_usable_bases: f64,

    // ── Strand specificity ──
    /// Strandedness used for the strand metrics: `unstranded`, `forward`, or `reverse`. With
    /// --strand auto (the default) this is inferred from the observed R1/R2 transcript-strand counts.
    pub detected_strand: String,
    /// Reads aligning on the strand expected under `detected_strand`. Always 0 when unstranded.
    pub correct_strand_reads: u64,
    /// `correct_strand_reads / (correct + incorrect)`. (The incorrect count is not emitted; it is
    /// `mapped exon-overlapping single-gene reads` minus `correct_strand_reads`.)
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_correct_strand_reads: f64,
    /// Templates (counted from read 1 / unpaired reads, FR-enclosed within a single gene) whose
    /// geometry implies read 1 lies on the transcription strand. Used for strand auto-detection.
    pub r1_tx_strand_reads: u64,
    /// Templates whose geometry implies read 2 lies on the transcription strand.
    pub r2_tx_strand_reads: u64,
    /// Templates that overlapped a single gene's exons but whose transcription strand could not be
    /// inferred (not properly paired / not FR-enclosed in one gene).
    pub unexplained_reads: u64,
    /// `r1_tx_strand_reads / (r1_tx_strand_reads + r2_tx_strand_reads)` — the basis for strand
    /// auto-detection (≈0.5 unstranded, ≈1 forward, ≈0 reverse).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_r1_tx_strand_reads: f64,
    /// `r2_tx_strand_reads / (r1_tx_strand_reads + r2_tx_strand_reads)`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_r2_tx_strand_reads: f64,

    // ── Coverage uniformity / 5'–3' bias / integrity ──
    /// Median across the top ~1000 expressed transcripts (one per gene, length ≥
    /// max(--minimum-length, 2×--end-bias-bases), mean coverage ≥ 1) of each transcript's coverage
    /// coefficient of variation (stddev/mean of per-base coverage). Lower = more even coverage.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub median_cv_coverage: f64,
    /// Median, over the same transcripts, of mean coverage in the 5' end window (--end-bias-bases
    /// wide) divided by whole-transcript mean coverage. ~1 is unbiased; <1 indicates 5' loss.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub median_5prime_bias: f64,
    /// As `median_5prime_bias` but for the 3' end window. >1 indicates 3' enrichment (e.g. from
    /// poly-A selection of degraded RNA).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub median_3prime_bias: f64,
    /// Median ratio of 5'-window to 3'-window coverage. ~1 is balanced; far from 1 signals a 5'/3'
    /// coverage skew.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub median_5prime_to_3prime_bias: f64,
    /// Median Transcript Integrity Number (0–100) — `100·exp(−KL(coverage‖uniform))`, an
    /// entropy-based coverage-uniformity score (100 = perfectly even). Computed over one
    /// representative transcript per gene (length ≥ --minimum-length, mean coverage >
    /// --tin-min-coverage). TIN is **scale-invariant**: it reflects relative coverage *shape*, not
    /// depth, so a deeply-sequenced transcript with a localized low-coverage stretch still scores
    /// low. See ERRATA.md (rna) for the full formula and interpretation.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub median_tin: f64,
    /// Mean TIN over the same transcripts.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub mean_tin: f64,
    /// Standard deviation of TIN over the same transcripts.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub stddev_tin: f64,
    /// Number of transcripts (genes) contributing to the TIN statistics, after the length and
    /// mean-coverage gates.
    pub tin_transcripts: u64,
}

/// Per-biotype uniquely-assigned read counts (one row per biotype).
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct RnaBiotypeMetric {
    /// Sample name.
    pub sample: String,
    /// Gene biotype (e.g. protein_coding, rRNA, lincRNA). For refFlat (no biotype column) this is
    /// the inferred `protein_coding`/`noncoding` split from CDS presence; `unknown` only when no
    /// biotype could be determined.
    pub biotype: String,
    /// Reads uniquely assigned to genes of this biotype.
    pub reads: u64,
    /// `reads / assigned_reads`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_reads: f64,
}

/// One point of the normalized-coverage-by-position series (0 = 5' end, 100 = 3' end) that backs
/// the coverage plot. Internal only — not written to a file.
#[derive(Debug)]
struct RnaCoverageEntry {
    normalized_position: u32,
    normalized_coverage: f64,
}

/// Transcript-space insert-size metrics, one row per pair orientation.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct RnaInsertSizeMetric {
    /// Sample name.
    pub sample: String,
    /// Pair orientation (FR, RF, or TANDEM).
    pub pair_orientation: String,
    /// Number of read pairs contributing to this orientation.
    pub read_pairs: u64,
    /// Mean transcript-space insert size after trimming outliers.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_insert_size: f64,
    /// Standard deviation after trimming outliers.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub standard_deviation: f64,
    /// Median transcript-space insert size.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub median_insert_size: f64,
    /// Median absolute deviation of insert size.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub median_absolute_deviation: f64,
    /// Smallest insert size observed.
    pub min_insert_size: u64,
    /// Largest insert size observed.
    pub max_insert_size: u64,
}

/// Transcript-space insert-size histogram with per-orientation counts at each size.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct RnaInsertSizeHistogramEntry {
    /// Sample name.
    pub sample: String,
    /// Insert size in transcript space (bp).
    pub insert_size: u64,
    /// FR pairs at this insert size.
    pub fr_count: u64,
    /// RF pairs at this insert size.
    pub rf_count: u64,
    /// Tandem pairs at this insert size.
    pub tandem_count: u64,
}

/// Read-level genomic origin of a primary, mapped, non-supplementary read, relative to the gene
/// model. `Exonic` carries the gene locus the read is uniquely assigned to.
enum ReadOrigin {
    /// Overlaps exon sequence of exactly one gene (uniquely assigned).
    Exonic(u32),
    /// Overlaps exon sequence of more than one gene.
    Ambiguous,
    /// Within a gene's span but touches no exon.
    Intronic,
    /// Overlaps no gene.
    Intergenic,
}

/// Coverage/bias medians produced by [`RnaCollector::finish_coverage`].
struct CoverageMetrics {
    /// Median coefficient of variation of coverage across the picked transcripts.
    cv: f64,
    /// Median normalized coverage over the 5' `end_bias_bases` of the picked transcripts.
    five_prime: f64,
    /// Median normalized coverage over the 3' `end_bias_bases` of the picked transcripts.
    three_prime: f64,
    /// Median 5'-to-3' coverage ratio across the picked transcripts.
    five_to_three: f64,
}

// ─── Module-level helpers ──────────────────────────────────────────────────────

/// True if the header declares coordinate sort order (`@HD SO:coordinate`).
fn is_coordinate_sorted(header: &Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| AsRef::<[u8]>::as_ref(so) == COORDINATE)
}

/// Bin `values` into `n_bins` equal-width bins spanning the data min..max, returning the
/// `n_bins + 1` bin edges and the per-bin counts. Returns empty vectors when `values` is empty
/// (or `n_bins` is 0). A degenerate range (all values equal) is widened so the bins are non-empty;
/// values on the top edge fall in the last bin.
#[allow(clippy::cast_precision_loss)]
fn histogram_bins(values: &[f64], n_bins: usize) -> (Vec<f64>, Vec<f64>) {
    if values.is_empty() || n_bins == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &v in values {
        min = min.min(v);
        max = max.max(v);
    }
    if (max - min).abs() < f64::EPSILON {
        min -= 0.5;
        max += 0.5;
    }
    let width = (max - min) / n_bins as f64;
    let edges: Vec<f64> = (0..=n_bins).map(|i| min + width * i as f64).collect();
    let counts = histogram_counts(values, &edges);
    (edges, counts)
}

/// Count `values` into the `edges.len() - 1` bins defined by the uniform, increasing `edges`,
/// returning the per-bin counts. Values outside `[edges[0], edges[last]]` are ignored; values on
/// the top edge fall in the last bin. Used to bin a second series onto a shared set of edges.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn histogram_counts(values: &[f64], edges: &[f64]) -> Vec<f64> {
    let n_bins = edges.len().saturating_sub(1);
    let mut counts = vec![0.0f64; n_bins];
    if n_bins == 0 {
        return counts;
    }
    let (min, max) = (edges[0], edges[n_bins]);
    let width = (max - min) / n_bins as f64;
    if width <= 0.0 {
        return counts;
    }
    for &v in values {
        if v < min || v > max {
            continue;
        }
        let idx = (((v - min) / width) as usize).min(n_bins - 1);
        counts[idx] += 1.0;
    }
    counts
}

/// Format an integer with comma thousands separators (e.g. `21075` → `"21,075"`).
fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    let first = digits.len() % 3;
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && i % 3 == first {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Prettify a biotype label for axis display: underscores become spaces and each word is
/// capitalized — **except** mixed-case tokens like `lncRNA` or `rRNA` (which start lowercase but
/// carry an uppercase letter later), which are left untouched so their conventional casing stands.
/// E.g. `protein_coding` → `Protein Coding`, `Mt_rRNA` → `Mt rRNA`, `lncRNA` → `lncRNA`.
fn prettify_biotype(label: &str) -> String {
    label
        .split('_')
        .map(|word| {
            let mut chars = word.chars();
            let Some(first) = chars.next() else { return String::new() };
            let mixed_case = first.is_lowercase() && word.chars().skip(1).any(char::is_uppercase);
            if mixed_case {
                word.to_string()
            } else {
                first.to_uppercase().collect::<String>() + chars.as_str()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Human-readable label for a gene biotype, used in the biotype file. Genes without a biotype
/// (e.g. from refFlat, which carries none) are grouped under `unknown`.
fn biotype_label(biotype: Option<&Biotype>) -> String {
    match biotype {
        Some(Biotype::ProteinCoding) => "protein_coding".to_string(),
        Some(Biotype::Rrna) => "rRNA".to_string(),
        Some(Biotype::MtRrna) => "Mt_rRNA".to_string(),
        Some(Biotype::RrnaPseudogene) => "rRNA_pseudogene".to_string(),
        Some(Biotype::Other(s)) => s.clone(),
        None => "unknown".to_string(),
    }
}

/// Convert an optional 1-based [`noodles::core::Position`] to a `u32`, or 0 if absent.
fn position_u32(pos: Option<noodles::core::Position>) -> u32 {
    #[allow(clippy::cast_possible_truncation)]
    pos.map_or(0, |p| usize::from(p) as u32)
}

/// Alignment blocks `(ref_start_1based, length)` for the M/=/X runs of a record's CIGAR.
fn rec_alignment_blocks(record: &RikerRecord) -> Vec<(u32, u32)> {
    let mut pos = position_u32(record.alignment_start());
    let mut blocks = Vec::new();
    for op in record.cigar_ops() {
        #[allow(clippy::cast_possible_truncation)]
        let len = op.len() as u32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Skip zero-length ops: downstream `block_start + block_len - 1` would underflow.
                if len > 0 {
                    blocks.push((pos, len));
                }
                pos += len;
            }
            Kind::Deletion | Kind::Skip => pos += len,
            _ => {}
        }
    }
    blocks
}

/// Intron gaps `(start_1based, end_1based)` for each CIGAR N (skip) op of length `>= min_intron`.
/// Returns a non-allocating empty vec for unspliced reads (the common case).
fn rec_introns(record: &RikerRecord, min_intron: u32) -> Vec<(u32, u32)> {
    let mut pos = position_u32(record.alignment_start());
    let mut introns = Vec::new();
    for op in record.cigar_ops() {
        #[allow(clippy::cast_possible_truncation)]
        let len = op.len() as u32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                pos += len;
            }
            Kind::Skip => {
                if len >= min_intron {
                    introns.push((pos, pos + len - 1));
                }
                pos += len;
            }
            _ => {}
        }
    }
    introns
}

/// The `MC` (mate CIGAR) tag value as bytes, if present.
fn mate_cigar(record: &RikerRecord) -> Option<Vec<u8>> {
    record.aux_tag(MC_TAG)?.as_str().map(|s| s.to_vec())
}

/// Reference end (1-based inclusive) used to bound the fragment span for the ribosomal and
/// strand-template checks (both callers take `read_end.max(this)`).
///
/// Prefers the exact mate end from the `MC` tag. When `MC` is absent it falls back to the
/// TLEN-derived fragment end — `min(read_start, mate_start) + |TLEN| - 1`, exactly Picard's
/// `CoordMath.getEnd(fragmentStart, abs(inferredInsertSize))` — never the read's own length.
/// A missing/zero `TLEN` is degenerate (as in Picard); we return `mate_start` so the mate
/// contributes only a single base rather than an underflowed span. (Transcript-space insert
/// size, by contrast, requires `MC` and never reaches this fallback.)
fn mate_alignment_end(record: &RikerRecord, mate_start: u32) -> u32 {
    if let Some((_, ref_end)) =
        mate_cigar(record).and_then(|mc| mate_alignment_blocks(&mc, mate_start))
    {
        return ref_end;
    }
    tlen_fragment_end(position_u32(record.alignment_start()), mate_start, record.template_length())
}

/// Fragment end (1-based inclusive) implied by `TLEN`: `min(read_start, mate_start) + |TLEN| - 1`,
/// matching Picard's `CoordMath.getEnd(fragmentStart, abs(inferredInsertSize))`. A zero/missing
/// `TLEN` is degenerate (as in Picard) and yields `mate_start` rather than an underflowed span.
fn tlen_fragment_end(read_start: u32, mate_start: u32, template_len: i32) -> u32 {
    let len = template_len.unsigned_abs();
    if len == 0 {
        return mate_start;
    }
    read_start.min(mate_start) + len - 1
}

/// Parse a mate CIGAR string starting at `mate_start` (1-based) into its M/=/X alignment
/// blocks `(ref_start, length)` and the 1-based inclusive reference end of the alignment.
///
/// The reference end advances through **all** reference-consuming ops (M/=/X/D/N), matching
/// fgbio's `lengthOnTarget`, so a CIGAR ending in a deletion or skip (e.g. `50M2000N`) still
/// reports the correct end — using only the last alignment block would undercount it.
/// Returns `None` if the CIGAR string is malformed.
fn mate_alignment_blocks(mc: &[u8], mate_start: u32) -> Option<(Vec<(u32, u32)>, u32)> {
    let mut pos = mate_start;
    let mut blocks = Vec::new();
    let mut len: u32 = 0;
    let mut saw_digit = false;
    for &b in mc {
        if b.is_ascii_digit() {
            len = len.checked_mul(10)?.checked_add(u32::from(b - b'0'))?;
            saw_digit = true;
        } else {
            if !saw_digit {
                return None;
            }
            match b {
                b'M' | b'=' | b'X' => {
                    // Skip zero-length ops: downstream `block_start + block_len - 1` underflows.
                    if len > 0 {
                        blocks.push((pos, len));
                    }
                    pos += len;
                }
                b'D' | b'N' => pos += len,
                b'I' | b'S' | b'H' | b'P' => {}
                _ => return None,
            }
            len = 0;
            saw_digit = false;
        }
    }
    if saw_digit {
        return None;
    }
    Some((blocks, pos.saturating_sub(1)))
}

/// Number of read bases in `blocks` that overlap the transcript's exons (1-based inclusive).
/// Ported from fgbio `numReadBasesOverlappingTranscript`.
fn bases_overlapping_exons(blocks: &[(u32, u32)], exons: &[crate::gene_model::Exon]) -> u32 {
    blocks
        .iter()
        .flat_map(|&(start, len)| {
            let block_end = start + len - 1;
            exons.iter().map(move |e| overlap_inclusive(start, block_end, e.start, e.end))
        })
        .sum()
}

/// The mate read's resolved reference span for insert-size tallying. `blocks` are its M/=/X
/// alignment blocks — the exact CIGAR blocks with `MC`, or a single whole-footprint block in the
/// TLEN fallback; `end` is its 1-based inclusive reference end; `mapped` is the denominator for the
/// exon-overlap fraction; `approximate` flags the MC-less fallback (which adds the
/// both-ends-in-an-exon requirement before a transcript qualifies).
struct MateSpan {
    blocks: Vec<(u32, u32)>,
    end: u32,
    mapped: u32,
    approximate: bool,
}

/// True if the 1-based position lies within any exon (1-based inclusive) of the slice.
fn point_in_exons(pos: u32, exons: &[crate::gene_model::Exon]) -> bool {
    exons.iter().any(|e| pos >= e.start && pos <= e.end)
}

/// Mean of a `u32` slice as `f64`; 0.0 for an empty slice.
fn mean_u32(values: &[u32]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let sum: u64 = values.iter().map(|&v| u64::from(v)).sum();
    sum as f64 / values.len() as f64
}

/// Transcript Integrity Number (0–100) from a transcript's per-base coverage profile:
/// `100 · e^H / n`, where `H = −Σ pᵢ·ln pᵢ` is the Shannon entropy of the normalized coverage
/// over all `n` exonic positions. Perfectly uniform coverage gives `e^H = n` (TIN 100); coverage
/// concentrated at one end lowers the entropy and the score. Returns `None` for zero coverage.
fn transcript_tin(cov: &[u32]) -> Option<f64> {
    let total: u64 = cov.iter().map(|&c| u64::from(c)).sum();
    if total == 0 {
        return None;
    }
    let total = total as f64;
    let mut entropy = 0.0;
    for &c in cov {
        if c > 0 {
            let p = f64::from(c) / total;
            entropy -= p * p.ln();
        }
    }
    Some(100.0 * entropy.exp() / cov.len() as f64)
}

/// Population standard deviation of a `u32` slice given its mean.
fn stddev_u32(values: &[u32], mean: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let var =
        values.iter().map(|&v| (f64::from(v) - mean).powi(2)).sum::<f64>() / values.len() as f64;
    var.max(0.0).sqrt()
}

/// Median of a slice of `f64` (sorts in place); 0.0 for an empty slice.
fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable_by(f64::total_cmp);
    let mid = values.len() / 2;
    if values.len() % 2 == 1 { values[mid] } else { f64::midpoint(values[mid - 1], values[mid]) }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn histogram_bins_partitions_all_values() {
        let values = [0.0, 1.0, 2.0, 3.0, 4.0];
        let (edges, counts) = histogram_bins(&values, 4);
        // n_bins + 1 edges spanning the data range, every value accounted for exactly once.
        assert_eq!(edges.len(), 5);
        assert_eq!(edges[0], 0.0);
        assert_eq!(edges[4], 4.0);
        assert_eq!(counts.iter().sum::<f64>(), values.len() as f64);
        // The top value falls in the last bin rather than overflowing.
        assert_eq!(counts[3], 2.0);
    }

    #[test]
    fn histogram_bins_widens_a_degenerate_range() {
        // All-equal values would give zero-width bins; the range is widened so binning is valid.
        let (edges, counts) = histogram_bins(&[2.5, 2.5, 2.5], 4);
        assert_eq!(edges.len(), 5);
        assert!(edges[4] > edges[0]);
        assert_eq!(counts.iter().sum::<f64>(), 3.0);
    }

    #[test]
    fn histogram_bins_is_empty_for_no_values() {
        let (edges, counts) = histogram_bins(&[], 4);
        assert!(edges.is_empty());
        assert!(counts.is_empty());
    }

    #[test]
    fn prettify_biotype_capitalizes_but_preserves_mixed_case_tokens() {
        assert_eq!(prettify_biotype("protein_coding"), "Protein Coding");
        assert_eq!(prettify_biotype("processed_pseudogene"), "Processed Pseudogene");
        assert_eq!(prettify_biotype("lncRNA"), "lncRNA");
        assert_eq!(prettify_biotype("Mt_rRNA"), "Mt rRNA");
        assert_eq!(prettify_biotype("misc_RNA"), "Misc RNA");
        assert_eq!(prettify_biotype("snoRNA"), "snoRNA");
        assert_eq!(prettify_biotype("TEC"), "TEC");
        assert_eq!(prettify_biotype("IG_V_gene"), "IG V Gene");
    }

    #[test]
    fn thousands_groups_digits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(42), "42");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_234), "1,234");
        assert_eq!(thousands(21_075), "21,075");
        assert_eq!(thousands(1_234_567), "1,234,567");
    }

    #[test]
    fn mate_alignment_blocks_handles_spliced_cigar() {
        // 50M2000N50M starting at 1000: blocks at 1000 (len 50) and 3050 (len 50); end 3099.
        let (blocks, end) = mate_alignment_blocks(b"50M2000N50M", 1000).unwrap();
        assert_eq!(blocks, vec![(1000, 50), (3050, 50)]);
        assert_eq!(end, 3099);
    }

    #[test]
    fn mate_alignment_end_includes_trailing_skip() {
        // A CIGAR ending in N must extend the reference end past the last aligned block
        // (matches fgbio lengthOnTarget). 50M2000N starting at 1000 → end 3049.
        let (blocks, end) = mate_alignment_blocks(b"50M2000N", 1000).unwrap();
        assert_eq!(blocks, vec![(1000, 50)]);
        assert_eq!(end, 3049);
    }

    #[test]
    fn mate_alignment_blocks_rejects_malformed() {
        assert!(mate_alignment_blocks(b"50", 1).is_none());
        assert!(mate_alignment_blocks(b"M50", 1).is_none());
        assert!(mate_alignment_blocks(b"50Z", 1).is_none());
    }

    #[test]
    fn tlen_fragment_end_uses_template_length_not_read_length() {
        // rec at 1000, mate at 1200, TLEN 350 → min(1000,1200) + 350 - 1 = 1349.
        assert_eq!(tlen_fragment_end(1000, 1200, 350), 1349);
        // Negative TLEN (this read is the rightmost mate) uses |TLEN| and the same min start.
        assert_eq!(tlen_fragment_end(1200, 1000, -350), 1349);
        // Zero/missing TLEN is degenerate: fall back to mate_start, never underflow.
        assert_eq!(tlen_fragment_end(1000, 1200, 0), 1200);
    }

    #[test]
    fn median_odd_and_even() {
        assert!((median(&mut [3.0, 1.0, 2.0]) - 2.0).abs() < 1e-9);
        assert!((median(&mut [1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-9);
        assert_eq!(median(&mut []), 0.0);
    }

    #[test]
    fn bases_overlapping_exons_sums_per_exon() {
        let exons = vec![
            crate::gene_model::Exon { start: 100, end: 200 },
            crate::gene_model::Exon { start: 300, end: 400 },
        ];
        // A block 150..=350 overlaps exon1 by 51 (150..=200) and exon2 by 51 (300..=350).
        assert_eq!(bases_overlapping_exons(&[(150, 201)], &exons), 51 + 51);
    }

    #[test]
    fn strand_detection_thresholds() {
        let mut c = RnaCollector::new(Path::new("x.bam"), Path::new("out"), &RnaOptions::default());
        c.num_r1 = 90;
        c.num_r2 = 10;
        assert_eq!(c.detect_strand(), StrandSpec::Forward);
        c.num_r1 = 10;
        c.num_r2 = 90;
        assert_eq!(c.detect_strand(), StrandSpec::Reverse);
        c.num_r1 = 50;
        c.num_r2 = 50;
        assert_eq!(c.detect_strand(), StrandSpec::Unstranded);
        c.num_r1 = 0;
        c.num_r2 = 0;
        assert_eq!(c.detect_strand(), StrandSpec::Unstranded);
    }
}
