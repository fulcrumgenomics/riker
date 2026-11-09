use std::fmt;
use std::path::{Path, PathBuf};

use anyhow::{Result, ensure};
use bitvec::prelude::*;
use clap::{Args, ValueEnum};
use kuva::plot::LinePlot;
use kuva::plot::legend::LegendPosition;
use kuva::render::plots::Plot;
use noodles::sam::Header;
use noodles::sam::alignment::record::Flags;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::{Collector, drive_collector_single_threaded};
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions, ReferenceOptions};
use crate::depth_histogram::DepthHistogram;
use crate::fasta::Fasta;
use crate::intervals::Intervals;
use crate::math::safe_div_f;
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, serialize_opt_f64_5dp, write_tsv};
use crate::plotting::{
    FG_BLUE, FG_GREEN, FG_NAVY, standard_layout, write_plot_pdf, write_twin_y_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::cigar::{aligned_blocks, count_aligned_bases, hard_clip_lengths};
use crate::sam::coordinate_order::CoordinateOrderGuard;
use crate::sam::mate_buffer::{MateBuffer, Peek};
use crate::sam::mate_overlap::{Alone, Buffering, OverlapAction, OverlapBitmap, Pairing};
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sam::{derive_sample, is_coordinate_sorted};
use crate::sequence_dict::SequenceDictionary;
use crate::simd::count_n_case_insensitive;

// ─── File suffixes ─────────────────────────────────────────────────────────────

/// File suffix for the summary metrics output.
pub const METRICS_SUFFIX: &str = ".methyl-metrics.txt";

/// File suffix for the CpG depth histogram.
pub const COVERAGE_SUFFIX: &str = ".methyl-coverage.txt";

/// File suffix for the CpG depth histogram plot.
pub const COVERAGE_PLOT_SUFFIX: &str = ".methyl-coverage.pdf";

/// File suffix for the per-contig metrics.
pub const CONTIGS_SUFFIX: &str = ".methyl-contigs.txt";

/// File suffix for the M-bias table.
pub const MBIAS_SUFFIX: &str = ".methyl-mbias.txt";

/// File suffix for the M-bias plot.
pub const MBIAS_PLOT_SUFFIX: &str = ".methyl-mbias.pdf";

/// File suffix for the CpG methylation-level histogram.
pub const BETA_SUFFIX: &str = ".methyl-beta.txt";

/// File suffix for the CpG methylation-level plot.
pub const BETA_PLOT_SUFFIX: &str = ".methyl-beta.pdf";

/// File suffix for the histogram of retained CpH per read.
pub const READ_CONVERSION_SUFFIX: &str = ".methyl-read-conversion.txt";

/// Largest sequencing cycle (0-based) tallied for M-bias; a guard against corrupt read lengths
/// growing the per-cycle tables without bound, far above any real read length.
const MAX_MBIAS_CYCLE: usize = 100_000;

/// Retained CpH counts at or above this are tallied together in the read-conversion histogram.
const MAX_RETAINED_CPH_BIN: usize = 20;

/// Reads with at least this many retained CpH count towards `frac_reads_cph_retained_ge_3`,
/// matching the rule NEB's `mark-nonconverted-reads` uses.
const RETAINED_CPH_FLAG_THRESHOLD: usize = 3;

/// Number of bins in the CpG methylation-level histogram: one per percentage point, 0–100.
const BETA_BINS: usize = 101;

/// Methylation levels below this count as low in `frac_cpgs_beta_low`.
const BETA_LOW: f64 = 0.2;

/// Methylation levels above this count as high in `frac_cpgs_beta_high`.
const BETA_HIGH: f64 = 0.8;

/// Two-sided 95% normal quantile, for the Agresti–Coull intervals on M-bias rates.
const Z_95: f64 = 1.959_963_984_540_054;

// ─── Options ──────────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the methylation collector.
#[riker_derive::multi_options("methyl", "Methyl Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct MethylOptions {
    /// Library chemistry, which decides whether a cytosine read as T was methylated.
    ///
    /// em-seq and bisulfite convert unmethylated cytosines to T; taps and
    /// illumina-5base convert methylated ones.
    #[arg(long, value_enum, default_value_t = MethylOptions::DEFAULT_CHEMISTRY)]
    pub chemistry: Chemistry,

    /// Interval file (IntervalList or BED) restricting the CpG territory.
    ///
    /// A CpG counts when its C is inside an interval. The base-level funnel,
    /// CpH conversion, M-bias and per-read metrics are not restricted.
    #[arg(short = 'L', long, value_name = "FILE", value_parser = crate::commands::common::parse_existing_file)]
    pub intervals: Option<PathBuf>,

    /// Contig holding a spike-in control (e.g. lambda, pUC19). Repeatable.
    ///
    /// Reads on a control contig are reported only in the per-contig table and
    /// as `frac_bases_excluded_control`; they contribute to no genome-level
    /// metric. Interval restriction never applies to control contigs.
    #[arg(long, value_name = "CONTIG", num_args = 1..)]
    pub control_contig: Option<Vec<String>>,

    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,

    /// Minimum mapping quality for a read to count.
    #[arg(long, default_value_t = MethylOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Minimum base quality for a base to count.
    #[arg(long, default_value_t = MethylOptions::DEFAULT_MIN_BQ)]
    pub min_bq: u8,

    /// Cycles to exclude at the 5' end of read 1 and single-end reads.
    ///
    /// Cycles count from the first base sequenced, hard-clipped bases included.
    #[arg(long, default_value_t = 0)]
    pub mask_r1_5p: u16,

    /// Cycles to exclude at the 3' end of read 1 and single-end reads.
    #[arg(long, default_value_t = 0)]
    pub mask_r1_3p: u16,

    /// Cycles to exclude at the 5' end of read 2.
    #[arg(long, default_value_t = 0)]
    pub mask_r2_5p: u16,

    /// Cycles to exclude at the 3' end of read 2.
    #[arg(long, default_value_t = 0)]
    pub mask_r2_3p: u16,

    /// Largest CpG depth used when computing depth statistics.
    ///
    /// Sites deeper than the cap are tallied at the cap in the depth histogram
    /// and statistics; observation counts are never capped.
    #[arg(long, default_value_t = MethylOptions::DEFAULT_COVERAGE_CAP)]
    pub coverage_cap: u16,

    /// Minimum depth for a CpG to enter the methylation-level histogram.
    #[arg(long, default_value_t = MethylOptions::DEFAULT_MIN_BETA_DEPTH)]
    pub min_beta_depth: u32,
}

impl MethylOptions {
    const DEFAULT_CHEMISTRY: Chemistry = Chemistry::EmSeq;
    const DEFAULT_MIN_MAPQ: u8 = 20;
    const DEFAULT_MIN_BQ: u8 = 20;
    const DEFAULT_COVERAGE_CAP: u16 = 250;
    const DEFAULT_MIN_BETA_DEPTH: u32 = 10;
}

impl Default for MethylOptions {
    fn default() -> Self {
        Self {
            chemistry: Self::DEFAULT_CHEMISTRY,
            intervals: None,
            control_contig: None,
            include_duplicates: false,
            min_mapq: Self::DEFAULT_MIN_MAPQ,
            min_bq: Self::DEFAULT_MIN_BQ,
            mask_r1_5p: 0,
            mask_r1_3p: 0,
            mask_r2_5p: 0,
            mask_r2_3p: 0,
            coverage_cap: Self::DEFAULT_COVERAGE_CAP,
            min_beta_depth: Self::DEFAULT_MIN_BETA_DEPTH,
        }
    }
}

// ─── CLI struct ───────────────────────────────────────────────────────────────

/// Collect methylation QC metrics from a bisulfite, EM-seq, TAPS or Illumina
/// 5-base BAM file.
///
/// Reports how many usable CpG observations the library yielded and where the
/// rest went, broken down by reason in two funnels: one over sequenced bases
/// (unmapped, control, QC-fail, duplicate, MAPQ, unaligned, end mask, base
/// quality, overlap) and one over CpG observations (the same read and base
/// filters, plus non-informative bases). Also reports CpG depth statistics,
/// per-context conversion and methylation, per-contig metrics (spike-in
/// controls included), M-bias, the CpG methylation-level distribution and
/// retained CpH per read.
///
/// The strand each read reports is inferred from its flags, assuming a
/// directional library: read 1 (or a single-end read) forward and read 2
/// reverse report the top strand. Secondary and supplementary alignments are
/// ignored.
///
/// The input SAM/BAM/CRAM must be coordinate-sorted (@HD SO:coordinate); a
/// record that arrives out of order aborts the run with an error.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker methyl -i input.bam -o out_prefix -r ref.fa
  riker methyl -i input.bam -o out_prefix -r ref.fa --chemistry illumina-5base
  riker methyl -i input.bam -o out_prefix -r ref.fa --control-contig lambda pUC19"
)]
pub struct Methyl {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub options: MethylOptions,
}

impl Command for Methyl {
    /// # Errors
    /// Returns an error if the BAM or reference cannot be read, or if output cannot be written.
    fn execute(&self, threads: Option<u8>) -> Result<()> {
        super::common::validate_output_prefix(&self.output.output)?;
        let plan = self.thread_plan(threads);
        let mut reader = AlignmentReader::open(
            &self.input.input,
            Some(&self.reference.reference),
            plan.decode_threads,
        )?;
        let reference = Fasta::from_path(&self.reference.reference)?;
        let mut collector =
            MethylCollector::new(&self.input.input, &self.output.output, reference, &self.options);
        let mut progress = ProgressLogger::new("methyl", "reads", 5_000_000);
        drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)
    }
}

// ─── Collector ────────────────────────────────────────────────────────────────

/// Accumulates methylation QC metrics over a coordinate-sorted BAM, one contig at a time: on
/// entering a contig it loads the reference and indexes its CpG sites, tallies each read's
/// bases and cytosine calls, and on leaving the contig folds the per-site depths into the
/// genome-wide statistics.
pub struct MethylCollector {
    // Output paths
    metrics_path: PathBuf,
    coverage_path: PathBuf,
    coverage_plot_path: PathBuf,
    contigs_path: PathBuf,
    mbias_path: PathBuf,
    mbias_plot_path: PathBuf,
    beta_path: PathBuf,
    beta_plot_path: PathBuf,
    read_conversion_path: PathBuf,

    // Input path, for the sample name derived in `initialize`.
    input_path: PathBuf,
    sample: String,

    // Config
    reference: Fasta,
    chemistry: Chemistry,
    intervals_path: Option<PathBuf>,
    // Parsed in `initialize` against the BAM header, so indexed by BAM `tid`.
    intervals: Option<Intervals>,
    control_contig_names: Vec<String>,
    include_duplicates: bool,
    min_mapq: u8,
    min_bq: u8,
    masks: EndMasks,
    min_beta_depth: u32,

    // BAM contig metadata, set in `initialize`.
    dict: Option<SequenceDictionary>,
    // Whether each BAM contig (by `tid`) is a control contig.
    is_control: Vec<bool>,
    order_guard: CoordinateOrderGuard,
    mate_buffer: MateBuffer<OverlapBitmap>,

    // The contig being processed, and the BAM contigs already processed.
    current: Option<ContigState>,
    processed: Vec<bool>,
    // Reference buffer reused across contigs.
    sequence_buffer: Vec<u8>,

    // Genome-wide accumulators (every non-control contig).
    reads: u64,
    bases_sequenced: u64,
    bases_unmapped: u64,
    bases_control: u64,
    genome: Tally,
    cpg_territory: u64,
    genome_non_n_bases: u64,
    depth_histogram: DepthHistogram,
    beta_histogram: Vec<u64>,
    beta_sum: f64,
    beta_low: u64,
    beta_high: u64,
    mbias: MbiasCounts,
    retained_cph_histogram: Vec<u64>,
    contig_rows: Vec<(usize, MethylContigMetric)>,
}

impl MethylCollector {
    /// Create a new collector. Output paths are derived from `prefix` by appending suffixes.
    /// The interval file, control contigs, sample name and contig names are resolved in
    /// [`Collector::initialize`] once the BAM header is available.
    #[must_use]
    pub fn new(input: &Path, prefix: &Path, reference: Fasta, options: &MethylOptions) -> Self {
        let path = |suffix| super::command::output_path(prefix, suffix);
        Self {
            metrics_path: path(METRICS_SUFFIX),
            coverage_path: path(COVERAGE_SUFFIX),
            coverage_plot_path: path(COVERAGE_PLOT_SUFFIX),
            contigs_path: path(CONTIGS_SUFFIX),
            mbias_path: path(MBIAS_SUFFIX),
            mbias_plot_path: path(MBIAS_PLOT_SUFFIX),
            beta_path: path(BETA_SUFFIX),
            beta_plot_path: path(BETA_PLOT_SUFFIX),
            read_conversion_path: path(READ_CONVERSION_SUFFIX),
            input_path: input.to_path_buf(),
            sample: String::new(),
            reference,
            chemistry: options.chemistry,
            intervals_path: options.intervals.clone(),
            intervals: None,
            control_contig_names: options.control_contig.clone().unwrap_or_default(),
            include_duplicates: options.include_duplicates,
            min_mapq: options.min_mapq,
            min_bq: options.min_bq,
            masks: EndMasks {
                r1_5p: usize::from(options.mask_r1_5p),
                r1_3p: usize::from(options.mask_r1_3p),
                r2_5p: usize::from(options.mask_r2_5p),
                r2_3p: usize::from(options.mask_r2_3p),
            },
            min_beta_depth: options.min_beta_depth,
            dict: None,
            is_control: Vec::new(),
            order_guard: CoordinateOrderGuard::new("methyl"),
            mate_buffer: MateBuffer::new(),
            current: None,
            processed: Vec::new(),
            sequence_buffer: Vec::new(),
            reads: 0,
            bases_sequenced: 0,
            bases_unmapped: 0,
            bases_control: 0,
            genome: Tally::default(),
            cpg_territory: 0,
            genome_non_n_bases: 0,
            depth_histogram: DepthHistogram::capped(options.coverage_cap),
            beta_histogram: vec![0; BETA_BINS],
            beta_sum: 0.0,
            beta_low: 0,
            beta_high: 0,
            mbias: MbiasCounts::default(),
            retained_cph_histogram: vec![0; MAX_RETAINED_CPH_BIN + 1],
            contig_rows: Vec::new(),
        }
    }

    /// Tally a mapped primary record on contig `ref_id` starting at 0-based `start`.
    fn accept_mapped(&mut self, record: &RikerRecord, ref_id: usize, start: usize) -> Result<()> {
        if self.current.as_ref().map(|contig| contig.ref_id) != Some(ref_id) {
            self.finalize_contig();
            self.begin_contig(ref_id)?;
        }

        let flags = record.flags();
        let sequence_len = record.sequence_len() as u64;
        let is_control = self.is_control[ref_id];
        if is_control {
            self.bases_control += sequence_len;
        }

        let top_strand = monitors_top_strand(flags);
        let exclusion = if flags.is_qc_fail() {
            Some(ReadExclusion::QcFail)
        } else if flags.is_duplicate() && !self.include_duplicates {
            Some(ReadExclusion::Duplicate)
        } else if record.mapping_quality().map_or(0, |mapq| mapq.get()) < self.min_mapq {
            Some(ReadExclusion::Mapq)
        } else {
            None
        };

        // `current` was set by `begin_contig` above.
        let contig = self.current.as_mut().unwrap();
        if let Some(exclusion) = exclusion {
            let potential = contig.sites.count_potential(record, start, top_strand);
            contig.tally.exclude_read(exclusion, sequence_len, potential);
            return Ok(());
        }
        contig.tally.bases.unaligned += sequence_len.saturating_sub(count_aligned_bases(record));

        let read = ReadInfo::new(record, start, top_strand, self.masks);
        let evidence = match self.mate_buffer.probe(record) {
            Peek::Alone => self.walk_read(record, &read, is_control, Alone),
            Peek::WouldBuffer { overlap_start, overlap_len } => {
                let mut bitmap = OverlapBitmap::new(overlap_start, overlap_len);
                let evidence =
                    self.walk_read(record, &read, is_control, Buffering { bitmap: &mut bitmap });
                self.mate_buffer.insert(record, bitmap);
                evidence
            }
            Peek::PairWith(cached) => {
                self.walk_read(record, &read, is_control, Pairing { cached: &cached })
            }
        };

        if !is_control && evidence.cph_sites > 0 {
            self.retained_cph_histogram[evidence.retained_cph.min(MAX_RETAINED_CPH_BIN)] += 1;
        }
        Ok(())
    }

    /// Walk a passing read's aligned bases, tallying each into the current contig's base and
    /// CpG funnels, conversion counts and per-site depths, and (off control contigs) into
    /// M-bias. Bases are excluded in order: end mask, base quality, then (CpG funnel only)
    /// non-informative base, then mate overlap. Returns the read's CpH evidence for the
    /// retained-CpH histogram.
    #[expect(
        clippy::too_many_lines,
        reason = "the per-base loop applies every filter, in order, in one place"
    )]
    fn walk_read<A: OverlapAction>(
        &mut self,
        record: &RikerRecord,
        read: &ReadInfo,
        is_control: bool,
        mut action: A,
    ) -> ReadEvidence {
        let min_bq = self.min_bq;
        let mbias = &mut self.mbias;
        // `current` is set before any read is walked.
        let contig = self.current.as_mut().unwrap();
        let reference = contig.sequence.as_slice();
        let sites = &mut contig.sites;
        let tally = &mut contig.tally;
        let quals = record.quality_scores();
        let bases = record.sequence();
        let mut evidence = ReadEvidence::default();
        let mut site_cursor = 0usize;

        for block in aligned_blocks(record) {
            let block_start = read.ref_start + block.ref_offset as usize;
            site_cursor = sites.seek(site_cursor, block_start.saturating_sub(1));
            for offset in 0..block.len as usize {
                let ref_pos = block_start + offset;
                let read_pos = block.read_offset as usize + offset;
                // A block past the contig end or a CIGAR longer than SEQ is malformed; the
                // remaining bases are dropped.
                if ref_pos >= reference.len() || read_pos >= bases.len() {
                    break;
                }
                let cycle = read.cycle(read_pos);
                let context = monitored_context(reference, ref_pos, read.top_strand);
                let call = context.map(|_| call_base(bases[read_pos], read.top_strand));
                let quality_ok = quals.get(read_pos).is_some_and(|&qual| qual >= min_bq);

                if quality_ok
                    && !is_control
                    && let (Some(context), Some(call)) = (context, call)
                    && call != Call::Other
                {
                    mbias.record(read.role, context == Context::CpG, cycle, call == Call::Retained);
                }

                // The site index of a CpG cytosine, when that site is in scope.
                let cpg_site = if context == Some(Context::CpG) {
                    let c_pos = if read.top_strand { ref_pos } else { ref_pos - 1 };
                    site_cursor = sites.seek(site_cursor, c_pos);
                    sites.in_scope(site_cursor).then_some(site_cursor)
                } else {
                    None
                };

                if read.is_masked(cycle) {
                    tally.bases.end_mask += 1;
                    if cpg_site.is_some() {
                        tally.cpg.potential += 1;
                        tally.cpg.end_mask += 1;
                    }
                    continue;
                }
                if !quality_ok {
                    tally.bases.baseq += 1;
                    if cpg_site.is_some() {
                        tally.cpg.potential += 1;
                        tally.cpg.baseq += 1;
                    }
                    continue;
                }

                if let (Some(context), Some(call)) = (context, call)
                    && context != Context::CpG
                    && call != Call::Other
                {
                    evidence.cph_sites += 1;
                    if call == Call::Retained {
                        evidence.retained_cph += 1;
                    }
                }

                // Positions fit in u32 for every supported reference.
                #[allow(clippy::cast_possible_truncation)]
                let pos = ref_pos as u32;
                let covered = action.is_mate_covered(pos);
                if covered {
                    tally.bases.overlap += 1;
                } else {
                    tally.bases.passing += 1;
                    action.on_counted(pos);
                }

                match (context, call) {
                    (Some(Context::CpG), Some(call)) => {
                        if let Some(site) = cpg_site {
                            tally.cpg.potential += 1;
                            if call == Call::Other {
                                tally.cpg.non_informative += 1;
                            } else if covered {
                                tally.cpg.overlap += 1;
                            } else {
                                tally.cpg.good += 1;
                                tally.count_call(Context::CpG, call);
                                sites.count_call(site, call);
                            }
                        }
                    }
                    (Some(context), Some(call)) => {
                        if call != Call::Other && !covered {
                            tally.count_call(context, call);
                        }
                    }
                    _ => {
                        if !covered
                            && let Some(call) = off_strand_cph_call(
                                reference,
                                ref_pos,
                                bases[read_pos],
                                read.top_strand,
                            )
                        {
                            tally.count_off_strand_call(call);
                        }
                    }
                }
            }
        }
        evidence
    }

    /// Start contig `ref_id`: load its reference and index its CpG sites, marking which are in
    /// scope (on a control contig, all of them; otherwise those inside `-L`, or all without it).
    ///
    /// # Errors
    /// Returns an error if the contig cannot be read from the reference.
    fn begin_contig(&mut self, ref_id: usize) -> Result<()> {
        let mut sequence = std::mem::take(&mut self.sequence_buffer);
        // `dict` is set in `initialize`, which always runs before `accept`.
        let name = self.dict.as_ref().unwrap()[ref_id].name().to_string();
        self.reference.load_contig_into(&name, true, &mut sequence)?;
        let sites = CpgSites::new(&sequence, self.scope_mask(ref_id));
        let non_n_bases = sequence.len() as u64 - count_n_case_insensitive(&sequence);
        self.mate_buffer.clear();
        self.current =
            Some(ContigState { ref_id, sequence, sites, tally: Tally::default(), non_n_bases });
        Ok(())
    }

    /// Finish the current contig, if any: fold its per-site depths into the depth and
    /// methylation-level histograms (genome contigs only), record its per-contig row, and fold
    /// its tallies into the genome totals unless it is a control.
    fn finalize_contig(&mut self) {
        let Some(contig) = self.current.take() else {
            return;
        };
        let is_control = self.is_control[contig.ref_id];
        let methylated_converts = self.chemistry.methylated_converts();
        if !is_control {
            for site in 0..contig.sites.len() {
                if !contig.sites.in_scope(site) {
                    continue;
                }
                let retained = u64::from(contig.sites.retained[site]);
                let converted = u64::from(contig.sites.converted[site]);
                let depth = retained + converted;
                self.depth_histogram.add(depth);
                if depth >= u64::from(self.min_beta_depth) && depth > 0 {
                    let methylated = if methylated_converts { converted } else { retained };
                    let beta = methylated as f64 / depth as f64;
                    self.add_beta(beta);
                }
            }
            self.genome.add(&contig.tally);
            self.cpg_territory += contig.sites.n_in_scope;
            self.genome_non_n_bases += contig.non_n_bases;
        }
        let row = self.contig_row(contig.ref_id, &contig.tally, contig.sites.n_in_scope);
        self.contig_rows.push((contig.ref_id, row));
        self.processed[contig.ref_id] = true;
        self.sequence_buffer = contig.sequence;
    }

    /// Tally one CpG's methylation level into the methylation-level histogram and summary.
    fn add_beta(&mut self, beta: f64) {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "beta is in [0, 1], so the bin is in [0, 100]"
        )]
        let bin = (beta * 100.0).round() as usize;
        self.beta_histogram[bin.min(BETA_BINS - 1)] += 1;
        self.beta_sum += beta;
        if beta < BETA_LOW {
            self.beta_low += 1;
        } else if beta > BETA_HIGH {
            self.beta_high += 1;
        }
    }

    /// Account for a contig no mapped read reached: its in-scope CpGs join the territory at depth
    /// zero, its non-N bases join the genome territory, and it gets an all-zero per-contig row.
    ///
    /// # Errors
    /// Returns an error if the contig cannot be read from the reference.
    fn account_unvisited_contig(&mut self, ref_id: usize) -> Result<()> {
        let is_control = self.is_control[ref_id];
        // A genome contig with no intervals under `-L` has no CpG territory, and its non-N bases
        // are unused because the CpG:genome ratio is not reported with `-L`; skip reading it.
        let skip_load =
            !is_control && self.intervals.as_ref().is_some_and(|ivl| !ivl.has_contig(ref_id));
        let (n_in_scope, non_n_bases) = if skip_load {
            (0, 0)
        } else {
            let mut sequence = std::mem::take(&mut self.sequence_buffer);
            // `dict` is set in `initialize`, which always runs before `finish`.
            let name = self.dict.as_ref().unwrap()[ref_id].name().to_string();
            self.reference.load_contig_into(&name, true, &mut sequence)?;
            let sites = CpgSites::new(&sequence, self.scope_mask(ref_id));
            let non_n_bases = sequence.len() as u64 - count_n_case_insensitive(&sequence);
            self.sequence_buffer = sequence;
            (sites.n_in_scope, non_n_bases)
        };
        if !is_control {
            self.depth_histogram.add_n(0, n_in_scope);
            self.cpg_territory += n_in_scope;
            self.genome_non_n_bases += non_n_bases;
        }
        let row = self.contig_row(ref_id, &Tally::default(), n_in_scope);
        self.contig_rows.push((ref_id, row));
        Ok(())
    }

    /// The per-site scope of contig `ref_id`, as a mask over its positions: `None` when every
    /// CpG is in scope (control contigs, or no `-L`).
    fn scope_mask(&self, ref_id: usize) -> Option<BitVec> {
        if self.is_control[ref_id] {
            return None;
        }
        self.intervals.as_ref().map(|intervals| intervals.contig_bitvec(ref_id))
    }

    /// The per-contig metrics row for contig `ref_id` from its tallies.
    fn contig_row(&self, ref_id: usize, tally: &Tally, cpg_sites: u64) -> MethylContigMetric {
        // `dict` is set in `initialize`, which always runs before contigs are finalized.
        let contig = self.dict.as_ref().unwrap()[ref_id].name().to_string();
        MethylContigMetric {
            sample: self.sample.clone(),
            contig,
            is_control: self.is_control[ref_id],
            cpg_sites,
            cpg_observations: tally.cpg.good,
            mean_cpg_depth: safe_div_f(tally.cpg.good as f64, cpg_sites as f64),
            cpg_methylation: tally.methylation(Context::CpG, self.chemistry),
            cph_conversion: tally.cph_conversion(),
        }
    }

    /// Finalize the last contig, account for contigs no read reached, and write every output.
    ///
    /// # Errors
    /// Returns an error if a reference contig cannot be read or an output cannot be written.
    fn finish_metrics(&mut self) -> Result<()> {
        self.finalize_contig();
        for ref_id in 0..self.processed.len() {
            if !self.processed[ref_id] {
                self.account_unvisited_contig(ref_id)?;
            }
        }
        self.contig_rows.sort_by_key(|(ref_id, _)| *ref_id);

        let metrics = self.summary_metrics();
        write_tsv(&self.metrics_path, &[metrics])?;
        self.write_coverage()?;
        let contig_rows: Vec<&MethylContigMetric> =
            self.contig_rows.iter().map(|(_, row)| row).collect();
        write_tsv(&self.contigs_path, &contig_rows)?;
        self.write_mbias()?;
        self.write_beta()?;
        self.write_read_conversion()?;

        log::info!(
            "methyl: cpg_territory={}, cpg_observations={}, metrics={}",
            self.cpg_territory,
            self.genome.cpg.good,
            self.metrics_path.display(),
        );
        Ok(())
    }

    /// Build the summary metrics row from the genome-wide accumulators.
    fn summary_metrics(&self) -> MethylMetrics {
        let bases = &self.genome.bases;
        let bases_excluded = [
            self.bases_unmapped,
            self.bases_control,
            bases.qc_fail,
            bases.dupe,
            bases.mapq,
            bases.unaligned,
            bases.end_mask,
            bases.baseq,
            bases.overlap,
        ];
        let bases_total = bases_excluded.iter().sum::<u64>() + bases.passing;
        let frac_bases = |n: u64| safe_div_f(n as f64, bases_total as f64);

        let cpg = &self.genome.cpg;
        let frac_cpg = |n: u64| safe_div_f(n as f64, cpg.potential as f64);

        let (mean, sd) = self.depth_histogram.mean_and_sd();
        let (median, mad) = self.depth_histogram.median_and_mad();
        let frac_at = |depth: usize| self.depth_histogram.frac_at_or_above(depth);

        let cpg_to_genome_depth_ratio = if self.intervals.is_none() {
            let cpg_depth = safe_div_f(cpg.good as f64, self.cpg_territory as f64);
            let genome_depth = safe_div_f(bases.passing as f64, self.genome_non_n_bases as f64);
            (genome_depth > 0.0).then(|| cpg_depth / genome_depth)
        } else {
            None
        };

        let beta_sites: u64 = self.beta_histogram.iter().sum();
        let reads_with_cph: u64 = self.retained_cph_histogram.iter().sum();
        let reads_flagged: u64 =
            self.retained_cph_histogram[RETAINED_CPH_FLAG_THRESHOLD..].iter().sum();

        MethylMetrics {
            sample: self.sample.clone(),
            chemistry: self.chemistry.to_string(),
            reads: self.reads,
            bases_sequenced: self.bases_sequenced,
            bases_passing: bases.passing,
            frac_bases_excluded_unmapped: frac_bases(self.bases_unmapped),
            frac_bases_excluded_control: frac_bases(self.bases_control),
            frac_bases_excluded_qc_fail: frac_bases(bases.qc_fail),
            frac_bases_excluded_dupe: frac_bases(bases.dupe),
            frac_bases_excluded_mapq: frac_bases(bases.mapq),
            frac_bases_excluded_unaligned: frac_bases(bases.unaligned),
            frac_bases_excluded_end_mask: frac_bases(bases.end_mask),
            frac_bases_excluded_baseq: frac_bases(bases.baseq),
            frac_bases_excluded_overlap: frac_bases(bases.overlap),
            frac_bases_excluded_total: frac_bases(bases_excluded.iter().sum()),
            cpg_territory: self.cpg_territory,
            cpg_observations_potential: cpg.potential,
            cpg_observations: cpg.good,
            cpg_observations_per_read: safe_div_f(cpg.good as f64, self.reads as f64),
            cpg_observations_per_mb: safe_div_f(
                cpg.good as f64,
                self.bases_sequenced as f64 / 1_000_000.0,
            ),
            frac_cpg_obs_excluded_qc_fail: frac_cpg(cpg.qc_fail),
            frac_cpg_obs_excluded_dupe: frac_cpg(cpg.dupe),
            frac_cpg_obs_excluded_mapq: frac_cpg(cpg.mapq),
            frac_cpg_obs_excluded_end_mask: frac_cpg(cpg.end_mask),
            frac_cpg_obs_excluded_baseq: frac_cpg(cpg.baseq),
            frac_cpg_obs_excluded_non_informative: frac_cpg(cpg.non_informative),
            frac_cpg_obs_excluded_overlap: frac_cpg(cpg.overlap),
            frac_cpg_obs_excluded_total: frac_cpg(cpg.potential - cpg.good),
            mean_cpg_depth: mean,
            sd_cpg_depth: sd,
            median_cpg_depth: median,
            mad_cpg_depth: mad,
            frac_cpgs_at_1x: frac_at(1),
            frac_cpgs_at_5x: frac_at(5),
            frac_cpgs_at_10x: frac_at(10),
            frac_cpgs_at_15x: frac_at(15),
            frac_cpgs_at_20x: frac_at(20),
            frac_cpgs_at_30x: frac_at(30),
            fold_80_cpg_penalty: self.depth_histogram.fold_penalty(mean, 0.20),
            cpg_to_genome_depth_ratio,
            cpg_methylation: self.genome.methylation(Context::CpG, self.chemistry),
            cpa_conversion: self.genome.conversion(Context::CpA),
            cpc_conversion: self.genome.conversion(Context::CpC),
            cpt_conversion: self.genome.conversion(Context::CpT),
            cph_conversion: self.genome.cph_conversion(),
            off_strand_cph_conversion: safe_div_f(
                self.genome.off_strand_converted as f64,
                (self.genome.off_strand_converted + self.genome.off_strand_retained) as f64,
            ),
            mean_cpg_beta: safe_div_f(self.beta_sum, beta_sites as f64),
            frac_cpgs_beta_low: safe_div_f(self.beta_low as f64, beta_sites as f64),
            frac_cpgs_beta_intermediate: safe_div_f(
                (beta_sites - self.beta_low - self.beta_high) as f64,
                beta_sites as f64,
            ),
            frac_cpgs_beta_high: safe_div_f(self.beta_high as f64, beta_sites as f64),
            frac_reads_cph_retained_ge_3: (!self.chemistry.methylated_converts())
                .then(|| safe_div_f(reads_flagged as f64, reads_with_cph as f64)),
        }
    }

    /// Write the CpG depth histogram and its plot.
    ///
    /// # Errors
    /// Returns an error if the outputs cannot be written.
    fn write_coverage(&self) -> Result<()> {
        let total = self.depth_histogram.total();
        let at_or_above = self.depth_histogram.at_or_above();
        let rows: Vec<MethylCoverageEntry> = self
            .depth_histogram
            .counts()
            .iter()
            .enumerate()
            .map(|(depth, &cpgs)| MethylCoverageEntry {
                sample: self.sample.clone(),
                depth: depth as u64,
                cpgs,
                frac_cpgs: safe_div_f(cpgs as f64, total as f64),
                cpgs_at_or_above: at_or_above[depth],
                frac_cpgs_at_or_above: safe_div_f(at_or_above[depth] as f64, total as f64),
            })
            .collect();
        write_tsv(&self.coverage_path, &rows)?;
        let (mean, _) = self.depth_histogram.mean_and_sd();
        self.depth_histogram.plot_pdf(
            mean,
            &format!("CpG Depth Distribution of {}", self.sample),
            "CpG Depth (X)",
            "Fraction of CpGs",
            &self.coverage_plot_path,
        )
    }

    /// Write the M-bias table and plot.
    ///
    /// # Errors
    /// Returns an error if the outputs cannot be written.
    fn write_mbias(&self) -> Result<()> {
        let methylated_converts = self.chemistry.methylated_converts();
        let mut rows: Vec<MethylMbiasEntry> = Vec::new();
        for role in ReadRole::ALL {
            for cpg in [true, false] {
                for (cycle, counts) in self.mbias.cycles(role, cpg).iter().enumerate() {
                    let total = counts.retained + counts.converted;
                    if total == 0 {
                        continue;
                    }
                    let methylated =
                        if methylated_converts { counts.converted } else { counts.retained };
                    let (ci_low, ci_high) = agresti_coull_interval(methylated, total);
                    rows.push(MethylMbiasEntry {
                        sample: self.sample.clone(),
                        read: role.label().to_string(),
                        context: if cpg { "CpG" } else { "CpH" }.to_string(),
                        cycle: cycle as u64 + 1,
                        methylated,
                        total,
                        frac_methylated: methylated as f64 / total as f64,
                        ci_low,
                        ci_high,
                    });
                }
            }
        }
        write_tsv(&self.mbias_path, &rows)?;
        self.plot_mbias(&rows)
    }

    /// Plot M-bias: CpG methylation by cycle against the left axis and CpH methylation (dashed)
    /// against the right, one colour per read. Cycles holding under 1% of a series' busiest
    /// cycle's calls (the tail past most reads' length) are left off the plot, where a handful of
    /// calls would swing the line; the table keeps them. Writes nothing when there is no data.
    ///
    /// # Errors
    /// Returns an error if the plot cannot be rendered or written.
    fn plot_mbias(&self, rows: &[MethylMbiasEntry]) -> Result<()> {
        let mut primary: Vec<Plot> = Vec::new();
        let mut secondary: Vec<Plot> = Vec::new();
        let (mut cph_min, mut cph_max) = (1.0_f64, 0.0_f64);
        for (role, color) in ReadRole::ALL.iter().zip([FG_BLUE, FG_GREEN, FG_NAVY]) {
            for cpg in [true, false] {
                let context = if cpg { "CpG" } else { "CpH" };
                let series: Vec<&MethylMbiasEntry> = rows
                    .iter()
                    .filter(|row| row.read == role.label() && row.context == context)
                    .collect();
                let busiest = series.iter().map(|row| row.total).max().unwrap_or(0);
                let xy: Vec<(f64, f64)> = series
                    .iter()
                    .filter(|row| row.total * 100 >= busiest)
                    .map(|row| (row.cycle as f64, row.frac_methylated))
                    .collect();
                if xy.is_empty() {
                    continue;
                }
                if !cpg {
                    for &(_, frac) in &xy {
                        cph_min = cph_min.min(frac);
                        cph_max = cph_max.max(frac);
                    }
                }
                let line = LinePlot::new()
                    .with_data(xy)
                    .with_color(color)
                    .with_stroke_width(1.5)
                    .with_legend(format!("{} {context}", role.label()));
                if cpg {
                    primary.push(line.into());
                } else {
                    secondary.push(line.with_dashed().into());
                }
            }
        }
        if primary.is_empty() && secondary.is_empty() {
            return Ok(());
        }
        // The right axis needs an explicit range; CpH methylation sits near 0 or near 1 with the
        // chemistry, so fit it to the data.
        let y2_pad = ((cph_max - cph_min) * 0.1).max(0.005);
        let layout = standard_layout(if primary.is_empty() { &secondary } else { &primary })
            .with_title(format!("M-bias of {}", self.sample))
            .with_x_label("Cycle")
            .with_y_label("CpG Methylation")
            .with_y2_label("CpH Methylation")
            .with_y2_range((cph_min - y2_pad).max(0.0), (cph_max + y2_pad).min(1.0))
            .with_legend_position(LegendPosition::InsideBottomRight);
        write_twin_y_plot_pdf(primary, secondary, layout, &self.mbias_plot_path)
    }

    /// Write the CpG methylation-level histogram and its plot.
    ///
    /// # Errors
    /// Returns an error if the outputs cannot be written.
    fn write_beta(&self) -> Result<()> {
        let total: u64 = self.beta_histogram.iter().sum();
        let rows: Vec<MethylBetaEntry> = self
            .beta_histogram
            .iter()
            .enumerate()
            .map(|(bin, &cpgs)| MethylBetaEntry {
                sample: self.sample.clone(),
                beta: bin as f64 / 100.0,
                cpgs,
                frac_cpgs: safe_div_f(cpgs as f64, total as f64),
            })
            .collect();
        write_tsv(&self.beta_path, &rows)?;
        if total == 0 {
            return Ok(());
        }
        let xy: Vec<(f64, f64)> = rows.iter().map(|row| (row.beta, row.frac_cpgs)).collect();
        let plots: Vec<Plot> = vec![
            LinePlot::new()
                .with_data(xy)
                .with_color(FG_BLUE)
                .with_fill()
                .with_fill_opacity(0.3)
                .into(),
        ];
        let layout = standard_layout(&plots)
            .with_title(format!(
                "CpG Methylation Levels of {} (depth ≥ {})",
                self.sample, self.min_beta_depth
            ))
            .with_x_label("Methylation Level")
            .with_y_label("Fraction of CpGs")
            .with_x_axis_min(0.0)
            .with_x_axis_max(1.0);
        write_plot_pdf(plots, layout, &self.beta_plot_path)
    }

    /// Write the histogram of retained CpH per read.
    ///
    /// # Errors
    /// Returns an error if the output cannot be written.
    fn write_read_conversion(&self) -> Result<()> {
        let total: u64 = self.retained_cph_histogram.iter().sum();
        let rows: Vec<MethylReadConversionEntry> = self
            .retained_cph_histogram
            .iter()
            .enumerate()
            .map(|(retained, &reads)| MethylReadConversionEntry {
                sample: self.sample.clone(),
                retained_cph: retained as u64,
                reads,
                frac_reads: safe_div_f(reads as f64, total as f64),
            })
            .collect();
        write_tsv(&self.read_conversion_path, &rows)
    }
}

impl Collector for MethylCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        ensure!(
            is_coordinate_sorted(header),
            "methyl requires coordinate-sorted input (@HD SO:coordinate); sort with `samtools sort`"
        );
        self.reference.validate_bam_header(header)?;
        let dict = SequenceDictionary::from(header);
        self.intervals = self
            .intervals_path
            .as_ref()
            .map(|path| Intervals::from_path(path, dict.clone()))
            .transpose()?;
        self.is_control = vec![false; dict.len()];
        for name in &self.control_contig_names {
            let Some(meta) = dict.get_by_name(name) else {
                anyhow::bail!("--control-contig {name} is not a contig in the BAM header");
            };
            self.is_control[meta.index()] = true;
        }
        self.processed = vec![false; dict.len()];
        self.dict = Some(dict);
        self.sample = derive_sample(&self.input_path, header);
        Ok(())
    }

    fn accept(&mut self, record: &RikerRecord, _header: &Header) -> Result<()> {
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }
        self.reads += 1;
        let sequence_len = record.sequence_len() as u64;
        self.bases_sequenced += sequence_len;

        let (Some(ref_id), Some(start)) =
            (record.reference_sequence_id(), record.alignment_start())
        else {
            self.bases_unmapped += sequence_len;
            return Ok(());
        };
        if flags.is_unmapped() {
            self.bases_unmapped += sequence_len;
            return Ok(());
        }
        // Positions fit in u32 for every supported reference.
        #[allow(clippy::cast_possible_truncation)]
        let start = (start.get() - 1) as u32;
        // `dict` is set in `initialize`, which always runs before `accept`.
        self.order_guard.check(ref_id, start, self.dict.as_ref().unwrap())?;
        self.accept_mapped(record, ref_id, start as usize)
    }

    fn finish(&mut self) -> Result<()> {
        self.finish_metrics()
    }

    fn name(&self) -> &'static str {
        "methyl"
    }

    fn field_needs(&self) -> RikerRecordRequirements {
        RikerRecordRequirements::NONE.with_sequence()
    }

    /// Relative per-record worker cost (see [`Collector::cost_hint`]). Not yet measured; set
    /// to `wgs`'s because both walk every aligned base with mate-overlap deduplication.
    fn cost_hint(&self) -> u32 {
        70
    }
}

// ─── Metric structs ───────────────────────────────────────────────────────────

/// Methylation QC summary metrics.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylMetrics {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Library chemistry (`--chemistry`).
    pub chemistry: String,
    /// Primary records (not secondary or supplementary), mapped or not. A BAM with unmapped
    /// reads removed understates this.
    pub reads: u64,
    /// Sum of the stored sequence lengths of primary records; soft-clipped bases count,
    /// hard-clipped bases do not.
    pub bases_sequenced: u64,
    /// Aligned bases on genome contigs that passed every filter: the base-level funnel's
    /// remainder.
    pub bases_passing: u64,
    /// Fraction of sequenced bases in unmapped reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_unmapped: f64,
    /// Fraction of sequenced bases in reads aligned to a control contig.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_control: f64,
    /// Fraction of sequenced bases in reads flagged QC-fail (e.g. unconverted templates
    /// flagged by methylsieve).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_qc_fail: f64,
    /// Fraction of sequenced bases in reads marked as duplicates.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_dupe: f64,
    /// Fraction of sequenced bases in reads below the minimum mapping quality, assessed after
    /// QC-fail and duplicates.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_mapq: f64,
    /// Fraction of sequenced bases of passing reads not aligned to the reference: soft clips
    /// and insertions.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_unaligned: f64,
    /// Fraction of sequenced bases excluded by the end masks (`--mask-*`).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_end_mask: f64,
    /// Fraction of sequenced bases below the minimum base quality.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_baseq: f64,
    /// Fraction of sequenced bases at positions the overlapping mate already counted.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_overlap: f64,
    /// Fraction of sequenced bases excluded for any reason.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_bases_excluded_total: f64,
    /// Number of reference CpG sites in the territory (both strands of a CpG are one site).
    pub cpg_territory: u64,
    /// Read bases aligned to the cytosine of a territory CpG on the strand the read reports,
    /// before any filter.
    pub cpg_observations_potential: u64,
    /// Potential CpG observations that passed every filter and read C or T (G or A on the
    /// bottom strand).
    pub cpg_observations: u64,
    /// Good CpG observations per primary read.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_observations_per_read: f64,
    /// Good CpG observations per million sequenced bases.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub cpg_observations_per_mb: f64,
    /// Fraction of potential CpG observations in QC-fail reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_qc_fail: f64,
    /// Fraction of potential CpG observations in duplicate reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_dupe: f64,
    /// Fraction of potential CpG observations in reads below the minimum mapping quality,
    /// assessed after QC-fail and duplicates. The loss among non-duplicate data is this
    /// divided by `1 - frac_cpg_obs_excluded_qc_fail - frac_cpg_obs_excluded_dupe`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_mapq: f64,
    /// Fraction of potential CpG observations excluded by the end masks.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_end_mask: f64,
    /// Fraction of potential CpG observations below the minimum base quality.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_baseq: f64,
    /// Fraction of potential CpG observations whose base was neither C nor T (G nor A on the
    /// bottom strand): a variant or sequencing error.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_non_informative: f64,
    /// Fraction of potential CpG observations at positions the overlapping mate already counted.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_overlap: f64,
    /// Fraction of potential CpG observations excluded for any reason.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpg_obs_excluded_total: f64,
    /// Mean CpG depth (good observations per territory CpG, capped at `--coverage-cap`).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_cpg_depth: f64,
    /// Standard deviation of CpG depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub sd_cpg_depth: f64,
    /// Median CpG depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub median_cpg_depth: f64,
    /// Median absolute deviation of CpG depth.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mad_cpg_depth: f64,
    /// Fraction of territory CpGs with depth ≥ 1.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_1x: f64,
    /// Fraction of territory CpGs with depth ≥ 5.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_5x: f64,
    /// Fraction of territory CpGs with depth ≥ 10.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_10x: f64,
    /// Fraction of territory CpGs with depth ≥ 15.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_15x: f64,
    /// Fraction of territory CpGs with depth ≥ 20.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_20x: f64,
    /// Fraction of territory CpGs with depth ≥ 30.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_30x: f64,
    /// Mean CpG depth divided by its 20th percentile.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub fold_80_cpg_penalty: f64,
    /// Uncapped mean CpG depth divided by the mean depth of passing bases over the non-N genome.
    /// With strands merged, 1.0 means no CpG-specific loss. Empty when `-L` is given.
    #[serde(serialize_with = "serialize_opt_f64_5dp")]
    pub cpg_to_genome_depth_ratio: Option<f64>,
    /// Fraction of good territory CpG observations that were methylated.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_methylation: f64,
    /// Fraction of CpA observations read as T. Under em-seq / bisulfite this is the conversion
    /// rate; under taps / illumina-5base it is the false-conversion rate.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpa_conversion: f64,
    /// Fraction of CpC observations read as T.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpc_conversion: f64,
    /// Fraction of CpT observations read as T.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpt_conversion: f64,
    /// Fraction of CpH observations (CpA, CpC, CpT pooled) read as T.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cph_conversion: f64,
    /// Fraction of CpH cytosines on the strand a read does *not* report that it shows as
    /// converted. Near the sequencing error rate for a directional library; a high value means
    /// the library is not directional.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub off_strand_cph_conversion: f64,
    /// Mean methylation level of territory CpGs at depth ≥ `--min-beta-depth`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub mean_cpg_beta: f64,
    /// Fraction of those CpGs with methylation level below 0.2.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_low: f64,
    /// Fraction of those CpGs with methylation level from 0.2 to 0.8.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_intermediate: f64,
    /// Fraction of those CpGs with methylation level above 0.8.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_high: f64,
    /// Fraction of passing reads with at least one CpH call that retained three or more CpH,
    /// a sign of failed conversion. Empty under taps / illumina-5base, where unmethylated CpH
    /// is expected to read as C.
    #[serde(serialize_with = "serialize_opt_f64_5dp")]
    pub frac_reads_cph_retained_ge_3: Option<f64>,
}

/// One depth bin of the CpG depth histogram.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylCoverageEntry {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// CpG depth (good observations, both strands); the last bin holds every CpG at or above
    /// the coverage cap.
    pub depth: u64,
    /// Number of territory CpGs at this depth.
    pub cpgs: u64,
    /// Fraction of territory CpGs at this depth.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs: f64,
    /// Number of territory CpGs at this depth or deeper.
    pub cpgs_at_or_above: u64,
    /// Fraction of territory CpGs at this depth or deeper.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_at_or_above: f64,
}

/// Methylation metrics for one contig.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylContigMetric {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Contig name.
    pub contig: String,
    /// Whether the contig was named with `--control-contig`.
    pub is_control: bool,
    /// Number of CpG sites on the contig in scope (inside `-L` for genome contigs).
    pub cpg_sites: u64,
    /// Good CpG observations on the contig.
    pub cpg_observations: u64,
    /// Good CpG observations per CpG site.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub mean_cpg_depth: f64,
    /// Fraction of good CpG observations that were methylated.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_methylation: f64,
    /// Fraction of CpH observations read as T.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cph_conversion: f64,
}

/// Methylation at one sequencing cycle for one read and context.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylMbiasEntry {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Read: `R1`, `R2`, or `SE` for single-end reads.
    pub read: String,
    /// Cytosine context: `CpG` or `CpH`.
    pub context: String,
    /// Sequencing cycle (1-based, hard-clipped bases included).
    pub cycle: u64,
    /// Methylated cytosine calls at this cycle.
    pub methylated: u64,
    /// Cytosine calls (C or T) at this cycle.
    pub total: u64,
    /// Fraction of calls that were methylated.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_methylated: f64,
    /// Lower bound of the 95% Agresti–Coull interval on `frac_methylated`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub ci_low: f64,
    /// Upper bound of the 95% Agresti–Coull interval on `frac_methylated`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub ci_high: f64,
}

/// One bin of the CpG methylation-level histogram.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylBetaEntry {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Methylation level, rounded to the nearest 0.01.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub beta: f64,
    /// Number of territory CpGs at depth ≥ `--min-beta-depth` with this methylation level.
    pub cpgs: u64,
    /// Fraction of those CpGs with this methylation level.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs: f64,
}

/// One bin of the histogram of retained CpH per read.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MethylReadConversionEntry {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Retained (unconverted) CpH calls in the read; the last bin holds 20 or more.
    pub retained_cph: u64,
    /// Passing reads with at least one CpH call that retained this many.
    pub reads: u64,
    /// Fraction of those reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_reads: f64,
}

// ─── Chemistry ────────────────────────────────────────────────────────────────

/// Library chemistry (`--chemistry`): which cytosines the conversion turns into T.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Chemistry {
    /// Enzymatic methyl-seq: unmethylated cytosines read as T.
    #[value(name = "em-seq")]
    EmSeq,
    /// Bisulfite sequencing: unmethylated cytosines read as T.
    #[value(name = "bisulfite")]
    Bisulfite,
    /// TAPS and TAPS+: methylated cytosines read as T.
    #[value(name = "taps")]
    Taps,
    /// Illumina 5-base: methylated cytosines read as T.
    #[value(name = "illumina-5base")]
    Illumina5Base,
}

impl Chemistry {
    /// Whether methylated cytosines are the ones that read as T.
    #[must_use]
    pub fn methylated_converts(self) -> bool {
        matches!(self, Chemistry::Taps | Chemistry::Illumina5Base)
    }
}

impl fmt::Display for Chemistry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Chemistry::EmSeq => "em-seq",
            Chemistry::Bisulfite => "bisulfite",
            Chemistry::Taps => "taps",
            Chemistry::Illumina5Base => "illumina-5base",
        })
    }
}

// ─── Per-contig state ─────────────────────────────────────────────────────────

/// The contig being processed: its reference, CpG sites and tallies.
struct ContigState {
    ref_id: usize,
    /// Upper-cased reference sequence.
    sequence: Vec<u8>,
    sites: CpgSites,
    tally: Tally,
    /// Reference bases other than N.
    non_n_bases: u64,
}

/// The CpG sites of one contig, with per-site call counts (both strands merged).
struct CpgSites {
    /// 0-based position of each site's C, ascending.
    positions: Vec<u32>,
    /// Which sites are in scope; `None` when all are.
    in_scope: Option<BitVec>,
    /// Number of sites in scope.
    n_in_scope: u64,
    /// Good observations per site that read as unconverted (C, or G on the bottom strand).
    retained: Vec<u32>,
    /// Good observations per site that read as converted (T, or A on the bottom strand).
    converted: Vec<u32>,
}

impl CpgSites {
    /// Index the CpG sites of `sequence` (upper-cased). With a `scope_mask` over the contig's
    /// positions, a site is in scope when its C is set in the mask.
    fn new(sequence: &[u8], scope_mask: Option<BitVec>) -> Self {
        // Positions fit in u32 for every supported reference.
        #[allow(clippy::cast_possible_truncation)]
        let positions: Vec<u32> = sequence
            .windows(2)
            .enumerate()
            .filter(|(_, pair)| pair == b"CG")
            .map(|(pos, _)| pos as u32)
            .collect();
        let in_scope: Option<BitVec> = scope_mask.map(|mask| {
            positions.iter().map(|&pos| mask.get(pos as usize).is_some_and(|bit| *bit)).collect()
        });
        let n_in_scope = in_scope.as_ref().map_or(positions.len(), |mask| mask.count_ones()) as u64;
        let len = positions.len();
        Self { positions, in_scope, n_in_scope, retained: vec![0; len], converted: vec![0; len] }
    }

    /// Number of CpG sites on the contig.
    fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether site `site` is in scope.
    fn in_scope(&self, site: usize) -> bool {
        self.in_scope.as_ref().is_none_or(|mask| mask[site])
    }

    /// The index of the first site at or after 0-based `pos`, searching from `from`. Callers
    /// walk a read left to right, so `from` only ever moves forward.
    fn seek(&self, from: usize, pos: usize) -> usize {
        from + self.positions[from..].partition_point(|&site_pos| (site_pos as usize) < pos)
    }

    /// Count a good call at site `site`.
    fn count_call(&mut self, site: usize, call: Call) {
        match call {
            Call::Retained => self.retained[site] = self.retained[site].saturating_add(1),
            Call::Converted => self.converted[site] = self.converted[site].saturating_add(1),
            Call::Other => {}
        }
    }

    /// The number of in-scope CpG cytosines `record`, starting at 0-based `start`, has aligned
    /// bases over on the strand it reports: its potential CpG observations.
    fn count_potential(&self, record: &RikerRecord, start: usize, top_strand: bool) -> u64 {
        let mut potential = 0u64;
        let mut cursor = 0usize;
        for block in aligned_blocks(record) {
            let block_start = start + block.ref_offset as usize;
            let block_end = block_start + block.len as usize;
            // A top-strand read reports the C of a site at `pos`; a bottom-strand read the G at
            // `pos + 1`.
            let (first, end) = if top_strand {
                (block_start, block_end)
            } else {
                (block_start.saturating_sub(1), block_end.saturating_sub(1))
            };
            cursor = self.seek(cursor, first);
            let last = self.seek(cursor, end);
            potential += match &self.in_scope {
                None => (last - cursor) as u64,
                Some(mask) => mask[cursor..last].count_ones() as u64,
            };
            cursor = last;
        }
        potential
    }
}

/// Counts for one contig (or, summed, the genome).
#[derive(Debug, Default, Clone)]
struct Tally {
    bases: BaseFunnel,
    cpg: CpgFunnel,
    /// Good calls that read unconverted, by [`Context`] (CpG counted only at in-scope sites).
    retained: [u64; 4],
    /// Good calls that read converted, by [`Context`].
    converted: [u64; 4],
    /// Off-strand CpH calls that read unconverted.
    off_strand_retained: u64,
    /// Off-strand CpH calls that read converted.
    off_strand_converted: u64,
}

impl Tally {
    /// Exclude a whole read: its sequenced bases and its potential CpG observations go to the
    /// read-level category.
    fn exclude_read(&mut self, exclusion: ReadExclusion, bases: u64, potential: u64) {
        self.cpg.potential += potential;
        match exclusion {
            ReadExclusion::QcFail => {
                self.bases.qc_fail += bases;
                self.cpg.qc_fail += potential;
            }
            ReadExclusion::Duplicate => {
                self.bases.dupe += bases;
                self.cpg.dupe += potential;
            }
            ReadExclusion::Mapq => {
                self.bases.mapq += bases;
                self.cpg.mapq += potential;
            }
        }
    }

    /// Count a good call in `context`.
    fn count_call(&mut self, context: Context, call: Call) {
        match call {
            Call::Retained => self.retained[context as usize] += 1,
            Call::Converted => self.converted[context as usize] += 1,
            Call::Other => {}
        }
    }

    /// Count an off-strand CpH call.
    fn count_off_strand_call(&mut self, call: Call) {
        match call {
            Call::Retained => self.off_strand_retained += 1,
            Call::Converted => self.off_strand_converted += 1,
            Call::Other => {}
        }
    }

    /// Fraction of calls in `context` that read converted.
    fn conversion(&self, context: Context) -> f64 {
        let converted = self.converted[context as usize];
        safe_div_f(converted as f64, (converted + self.retained[context as usize]) as f64)
    }

    /// Fraction of CpH calls (CpA, CpC, CpT pooled) that read converted.
    fn cph_conversion(&self) -> f64 {
        let cph = [Context::CpA, Context::CpC, Context::CpT];
        let converted: u64 = cph.iter().map(|&context| self.converted[context as usize]).sum();
        let retained: u64 = cph.iter().map(|&context| self.retained[context as usize]).sum();
        safe_div_f(converted as f64, (converted + retained) as f64)
    }

    /// Fraction of calls in `context` that were methylated under `chemistry`.
    fn methylation(&self, context: Context, chemistry: Chemistry) -> f64 {
        let conversion = self.conversion(context);
        let total = self.converted[context as usize] + self.retained[context as usize];
        if total == 0 {
            0.0
        } else if chemistry.methylated_converts() {
            conversion
        } else {
            1.0 - conversion
        }
    }

    /// Add `other`'s counts to these.
    fn add(&mut self, other: &Tally) {
        self.bases.add(&other.bases);
        self.cpg.add(&other.cpg);
        for context in 0..4 {
            self.retained[context] += other.retained[context];
            self.converted[context] += other.converted[context];
        }
        self.off_strand_retained += other.off_strand_retained;
        self.off_strand_converted += other.off_strand_converted;
    }
}

/// Sequenced bases of mapped primary reads by the category that excluded them, or `passing`.
#[derive(Debug, Default, Clone)]
struct BaseFunnel {
    qc_fail: u64,
    dupe: u64,
    mapq: u64,
    unaligned: u64,
    end_mask: u64,
    baseq: u64,
    overlap: u64,
    passing: u64,
}

impl BaseFunnel {
    /// Add `other`'s counts to these.
    fn add(&mut self, other: &BaseFunnel) {
        self.qc_fail += other.qc_fail;
        self.dupe += other.dupe;
        self.mapq += other.mapq;
        self.unaligned += other.unaligned;
        self.end_mask += other.end_mask;
        self.baseq += other.baseq;
        self.overlap += other.overlap;
        self.passing += other.passing;
    }
}

/// Potential CpG observations by the category that excluded them, or `good`.
#[derive(Debug, Default, Clone)]
struct CpgFunnel {
    potential: u64,
    qc_fail: u64,
    dupe: u64,
    mapq: u64,
    end_mask: u64,
    baseq: u64,
    non_informative: u64,
    overlap: u64,
    good: u64,
}

impl CpgFunnel {
    /// Add `other`'s counts to these.
    fn add(&mut self, other: &CpgFunnel) {
        self.potential += other.potential;
        self.qc_fail += other.qc_fail;
        self.dupe += other.dupe;
        self.mapq += other.mapq;
        self.end_mask += other.end_mask;
        self.baseq += other.baseq;
        self.non_informative += other.non_informative;
        self.overlap += other.overlap;
        self.good += other.good;
    }
}

/// Why a whole read was excluded, in the order the filters apply.
#[derive(Debug, Clone, Copy)]
enum ReadExclusion {
    QcFail,
    Duplicate,
    Mapq,
}

// ─── Per-read state ───────────────────────────────────────────────────────────

/// Cycles to exclude at each end of read 1 (and single-end reads) and read 2.
#[derive(Debug, Clone, Copy)]
struct EndMasks {
    r1_5p: usize,
    r1_3p: usize,
    r2_5p: usize,
    r2_3p: usize,
}

/// What a read walk needs to know about the read beyond its bases.
struct ReadInfo {
    /// 0-based alignment start.
    ref_start: usize,
    top_strand: bool,
    role: ReadRole,
    reverse: bool,
    /// Stored sequence length.
    sequence_len: usize,
    /// Hard-clipped bases before the first sequenced base, in sequencing order.
    leading_hard_clip: usize,
    /// Full length of the read as sequenced: stored bases plus hard clips.
    sequenced_len: usize,
    mask_5p: usize,
    mask_3p: usize,
}

impl ReadInfo {
    /// Describe `record`, aligned from 0-based `start`.
    fn new(record: &RikerRecord, start: usize, top_strand: bool, masks: EndMasks) -> Self {
        let flags = record.flags();
        let role = ReadRole::of(flags);
        let reverse = flags.is_reverse_complemented();
        let (clip_left, clip_right) = hard_clip_lengths(record);
        let (clip_left, clip_right) = (clip_left as usize, clip_right as usize);
        let sequence_len = record.sequence_len();
        let (mask_5p, mask_3p) = match role {
            ReadRole::R2 => (masks.r2_5p, masks.r2_3p),
            ReadRole::R1 | ReadRole::Fragment => (masks.r1_5p, masks.r1_3p),
        };
        Self {
            ref_start: start,
            top_strand,
            role,
            reverse,
            sequence_len,
            leading_hard_clip: if reverse { clip_right } else { clip_left },
            sequenced_len: sequence_len + clip_left + clip_right,
            mask_5p,
            mask_3p,
        }
    }

    /// The 0-based sequencing cycle of the base at `read_pos` in the stored sequence.
    fn cycle(&self, read_pos: usize) -> usize {
        let from_5p = if self.reverse { self.sequence_len - 1 - read_pos } else { read_pos };
        self.leading_hard_clip + from_5p
    }

    /// Whether the base at 0-based `cycle` falls in an end mask.
    fn is_masked(&self, cycle: usize) -> bool {
        cycle < self.mask_5p || self.sequenced_len - 1 - cycle < self.mask_3p
    }
}

/// A read's CpH calls (after the end mask and base quality, before overlap removal), for the
/// retained-CpH histogram.
#[derive(Debug, Default)]
struct ReadEvidence {
    cph_sites: usize,
    retained_cph: usize,
}

/// Which read of a template a record is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReadRole {
    R1,
    R2,
    /// A single-end read.
    Fragment,
}

impl ReadRole {
    const ALL: [ReadRole; 3] = [ReadRole::R1, ReadRole::R2, ReadRole::Fragment];

    /// The role of a record from its flags.
    fn of(flags: Flags) -> Self {
        if !flags.is_segmented() {
            ReadRole::Fragment
        } else if flags.is_last_segment() {
            ReadRole::R2
        } else {
            ReadRole::R1
        }
    }

    /// The label used in outputs.
    fn label(self) -> &'static str {
        match self {
            ReadRole::R1 => "R1",
            ReadRole::R2 => "R2",
            ReadRole::Fragment => "SE",
        }
    }
}

/// The context of a cytosine on the strand a read reports, from the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Context {
    CpA = 0,
    CpC = 1,
    CpG = 2,
    CpT = 3,
}

/// A read's base at a cytosine it reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Call {
    /// Read as C (G on the bottom strand): not converted.
    Retained,
    /// Read as T (A on the bottom strand): converted.
    Converted,
    /// Any other base: a variant or sequencing error.
    Other,
}

// ─── M-bias ───────────────────────────────────────────────────────────────────

/// Retained and converted calls at one cycle.
#[derive(Debug, Default, Clone, Copy)]
struct CycleCounts {
    retained: u64,
    converted: u64,
}

/// Per-cycle call counts for each read role and context class (CpG or CpH).
#[derive(Debug, Default)]
struct MbiasCounts {
    /// Indexed by `role * 2 + (CpG ? 0 : 1)`, then by 0-based cycle.
    slots: [Vec<CycleCounts>; 6],
}

impl MbiasCounts {
    /// Record a call at 0-based `cycle`.
    fn record(&mut self, role: ReadRole, cpg: bool, cycle: usize, retained: bool) {
        if cycle > MAX_MBIAS_CYCLE {
            return;
        }
        let slot = &mut self.slots[Self::slot(role, cpg)];
        if cycle >= slot.len() {
            slot.resize(cycle + 1, CycleCounts::default());
        }
        if retained {
            slot[cycle].retained += 1;
        } else {
            slot[cycle].converted += 1;
        }
    }

    /// The per-cycle counts for a role and context class, by 0-based cycle.
    fn cycles(&self, role: ReadRole, cpg: bool) -> &[CycleCounts] {
        &self.slots[Self::slot(role, cpg)]
    }

    fn slot(role: ReadRole, cpg: bool) -> usize {
        let role_index = match role {
            ReadRole::R1 => 0,
            ReadRole::R2 => 1,
            ReadRole::Fragment => 2,
        };
        role_index * 2 + usize::from(!cpg)
    }
}

// ─── Functions ────────────────────────────────────────────────────────────────

/// Whether a read reports the top strand (reference C) rather than the bottom strand
/// (reference G), assuming a directional library: read 1 or a single-end read forward, or read
/// 2 reverse.
fn monitors_top_strand(flags: Flags) -> bool {
    flags.is_reverse_complemented() == (flags.is_segmented() && flags.is_last_segment())
}

/// The context of the cytosine at 0-based `pos` on the reported strand, or `None` if there is
/// no cytosine there on that strand or its neighbour is not A, C, G or T. A top-strand cytosine
/// is a reference C followed by its next base; a bottom-strand cytosine is a reference G,
/// whose next base on the bottom strand is the complement of the preceding reference base.
fn monitored_context(reference: &[u8], pos: usize, top_strand: bool) -> Option<Context> {
    if top_strand {
        if reference[pos] != b'C' {
            return None;
        }
        match reference.get(pos + 1)? {
            b'A' => Some(Context::CpA),
            b'C' => Some(Context::CpC),
            b'G' => Some(Context::CpG),
            b'T' => Some(Context::CpT),
            _ => None,
        }
    } else {
        if reference[pos] != b'G' || pos == 0 {
            return None;
        }
        match reference[pos - 1] {
            b'T' => Some(Context::CpA),
            b'G' => Some(Context::CpC),
            b'C' => Some(Context::CpG),
            b'A' => Some(Context::CpT),
            _ => None,
        }
    }
}

/// Classify a read base at a cytosine on the reported strand.
fn call_base(base: u8, top_strand: bool) -> Call {
    match (base.to_ascii_uppercase(), top_strand) {
        (b'C', true) | (b'G', false) => Call::Retained,
        (b'T', true) | (b'A', false) => Call::Converted,
        _ => Call::Other,
    }
}

/// The call at a CpH cytosine on the strand the read does *not* report, if the reference has
/// one at 0-based `pos`: a reference G not preceded by C for a top-strand read, or a reference C
/// not followed by G for a bottom-strand read.
fn off_strand_cph_call(reference: &[u8], pos: usize, base: u8, top_strand: bool) -> Option<Call> {
    let off_strand_cph = if top_strand {
        reference[pos] == b'G' && pos > 0 && matches!(reference[pos - 1], b'A' | b'G' | b'T')
    } else {
        reference[pos] == b'C' && matches!(reference.get(pos + 1), Some(b'A' | b'C' | b'T'))
    };
    off_strand_cph.then(|| call_base(base, !top_strand)).filter(|&call| call != Call::Other)
}

/// The 95% Agresti–Coull interval on a proportion of `successes` in `total` trials, clamped to
/// `[0, 1]`.
fn agresti_coull_interval(successes: u64, total: u64) -> (f64, f64) {
    let z_squared = Z_95 * Z_95;
    let adjusted_total = total as f64 + z_squared;
    let adjusted_rate = (successes as f64 + z_squared / 2.0) / adjusted_total;
    let half_width = Z_95 * (adjusted_rate * (1.0 - adjusted_rate) / adjusted_total).sqrt();
    ((adjusted_rate - half_width).max(0.0), (adjusted_rate + half_width).min(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_one_forward_and_read_two_reverse_report_the_top_strand() {
        let r1 = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let r2 = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        assert!(monitors_top_strand(r1));
        assert!(monitors_top_strand(r2 | Flags::REVERSE_COMPLEMENTED));
        assert!(!monitors_top_strand(r1 | Flags::REVERSE_COMPLEMENTED));
        assert!(!monitors_top_strand(r2));
    }

    #[test]
    fn single_end_reads_report_the_strand_they_align_to() {
        assert!(monitors_top_strand(Flags::empty()));
        assert!(!monitors_top_strand(Flags::REVERSE_COMPLEMENTED));
    }

    #[test]
    fn top_strand_context_comes_from_the_next_reference_base() {
        let reference = b"CACCCGCTCN";
        assert_eq!(monitored_context(reference, 0, true), Some(Context::CpA));
        assert_eq!(monitored_context(reference, 2, true), Some(Context::CpC));
        assert_eq!(monitored_context(reference, 4, true), Some(Context::CpG));
        assert_eq!(monitored_context(reference, 6, true), Some(Context::CpT));
        assert_eq!(monitored_context(reference, 8, true), None);
        assert_eq!(monitored_context(reference, 1, true), None);
    }

    #[test]
    fn bottom_strand_context_comes_from_the_complement_of_the_previous_base() {
        let reference = b"TGGGCGAGNG";
        assert_eq!(monitored_context(reference, 1, false), Some(Context::CpA));
        assert_eq!(monitored_context(reference, 2, false), Some(Context::CpC));
        assert_eq!(monitored_context(reference, 5, false), Some(Context::CpG));
        assert_eq!(monitored_context(reference, 7, false), Some(Context::CpT));
        assert_eq!(monitored_context(reference, 9, false), None);
        assert_eq!(monitored_context(b"GA", 0, false), None);
    }

    #[test]
    fn context_at_the_contig_end_is_unknown() {
        assert_eq!(monitored_context(b"AC", 1, true), None);
    }

    #[test]
    fn call_base_reads_c_and_t_on_top_and_g_and_a_on_bottom() {
        assert_eq!(call_base(b'C', true), Call::Retained);
        assert_eq!(call_base(b'T', true), Call::Converted);
        assert_eq!(call_base(b'A', true), Call::Other);
        assert_eq!(call_base(b'G', false), Call::Retained);
        assert_eq!(call_base(b'A', false), Call::Converted);
        assert_eq!(call_base(b'T', false), Call::Other);
    }

    #[test]
    fn off_strand_call_uses_the_other_strands_cph_cytosines() {
        // Top-strand read at a reference G preceded by A: a bottom-strand CpT.
        assert_eq!(off_strand_cph_call(b"AG", 1, b'A', true), Some(Call::Converted));
        assert_eq!(off_strand_cph_call(b"AG", 1, b'G', true), Some(Call::Retained));
        // A G preceded by C is a CpG, which is excluded.
        assert_eq!(off_strand_cph_call(b"CG", 1, b'A', true), None);
        // Bottom-strand read at a reference C followed by A.
        assert_eq!(off_strand_cph_call(b"CA", 0, b'T', false), Some(Call::Converted));
        assert_eq!(off_strand_cph_call(b"CG", 0, b'T', false), None);
    }

    #[test]
    fn cpg_sites_index_every_cg_dinucleotide() {
        let sites = CpgSites::new(b"ACGTTCGCG", None);
        assert_eq!(sites.positions, vec![1, 5, 7]);
        assert_eq!(sites.n_in_scope, 3);
    }

    #[test]
    fn cpg_sites_in_scope_follow_the_mask_at_the_c() {
        let mut mask = bitvec![0; 9];
        mask.set(5, true);
        mask.set(8, true);
        let sites = CpgSites::new(b"ACGTTCGCG", Some(mask));
        assert_eq!(sites.n_in_scope, 1);
        assert!(!sites.in_scope(0));
        assert!(sites.in_scope(1));
        assert!(!sites.in_scope(2));
    }

    #[test]
    fn seek_finds_the_first_site_at_or_after_a_position() {
        let sites = CpgSites::new(b"ACGTTCGCG", None);
        assert_eq!(sites.seek(0, 0), 0);
        assert_eq!(sites.seek(0, 2), 1);
        assert_eq!(sites.seek(1, 7), 2);
        assert_eq!(sites.seek(0, 8), 3);
    }

    #[test]
    fn cycles_count_from_the_five_prime_end_including_hard_clips() {
        let masks = EndMasks { r1_5p: 0, r1_3p: 0, r2_5p: 0, r2_3p: 0 };
        let forward = crate::test_support::read().cigar("8H10M").into_riker_record();
        let info = ReadInfo::new(&forward, 0, true, masks);
        assert_eq!(info.cycle(0), 8);
        assert_eq!(info.cycle(9), 17);

        let reverse = crate::test_support::read().cigar("10M8H").reverse().into_riker_record();
        let info = ReadInfo::new(&reverse, 0, false, masks);
        assert_eq!(info.cycle(9), 8);
        assert_eq!(info.cycle(0), 17);
    }

    #[test]
    fn end_masks_cover_cycles_from_each_end() {
        let masks = EndMasks { r1_5p: 2, r1_3p: 1, r2_5p: 0, r2_3p: 0 };
        let record = crate::test_support::read().cigar("10M").into_riker_record();
        let info = ReadInfo::new(&record, 0, true, masks);
        assert!(info.is_masked(0));
        assert!(info.is_masked(1));
        assert!(!info.is_masked(2));
        assert!(!info.is_masked(8));
        assert!(info.is_masked(9));
    }

    #[test]
    fn agresti_coull_interval_brackets_the_rate() {
        let (low, high) = agresti_coull_interval(50, 100);
        assert!(low < 0.5 && high > 0.5);
        crate::assert_close!(0.5 - low, high - 0.5, 1e-9);
    }

    #[test]
    fn agresti_coull_interval_stays_within_zero_and_one() {
        let (low, high) = agresti_coull_interval(0, 10);
        crate::assert_close!(low, 0.0);
        assert!(high > 0.0 && high < 1.0);
    }

    #[test]
    fn methylation_follows_the_chemistry() {
        let mut tally = Tally::default();
        tally.retained[Context::CpG as usize] = 3;
        tally.converted[Context::CpG as usize] = 1;
        crate::assert_close!(tally.methylation(Context::CpG, Chemistry::EmSeq), 0.75);
        crate::assert_close!(tally.methylation(Context::CpG, Chemistry::Taps), 0.25);
    }
}
