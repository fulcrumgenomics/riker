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
use crate::metrics::{serialize_f64_2dp, serialize_f64_5dp, write_tsv};
use crate::plotting::{
    FG_BLUE, FG_GREEN, FG_NAVY, standard_layout, write_plot_pdf, write_twin_y_plot_pdf,
};
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::cigar::{aligned_blocks, hard_clip_lengths};
use crate::sam::coordinate_order::CoordinateOrderGuard;
use crate::sam::mate_buffer::{MateBuffer, Peek};
use crate::sam::mate_overlap::{Alone, Buffering, OverlapAction, OverlapBitmap, Pairing};
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use crate::sam::{derive_sample, is_coordinate_sorted};
use crate::sequence_dict::SequenceDictionary;
use crate::simd;

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

/// Number of bins in the CpG methylation-level histogram: one per percentage point, 0–100.
const BETA_BINS: usize = 101;

/// Methylation levels below this count as low in `frac_cpgs_beta_low`.
const BETA_LOW: f64 = 0.2;

/// Methylation levels above this count as high in `frac_cpgs_beta_high`.
const BETA_HIGH: f64 = 0.8;

/// Two-sided 95% normal quantile, for the Agresti–Coull intervals on M-bias rates.
const Z_95: f64 = 1.959_963_984_540_054;

/// Mask for one strand's field in a context code (see [`encode_contexts`]).
const CONTEXT_MASK: u8 = 0b111;

/// Shift of the bottom-strand field in a context code.
const BOTTOM_SHIFT: u8 = 3;

/// Mask for both strand fields of a context code.
const STRAND_FIELDS_MASK: u8 = 0b11_1111;

/// Context-code bit set at both cytosines of an in-scope CpG.
const SCOPE_BIT: u8 = 0b100_0000;

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
    /// A CpG counts when its C is inside an interval. CpH conversion, M-bias
    /// and the per-read metrics are not restricted.
    #[arg(short = 'L', long, value_name = "FILE", value_parser = crate::commands::common::parse_existing_file)]
    pub intervals: Option<PathBuf>,

    /// Contig holding a spike-in control (e.g. lambda, pUC19). Repeatable.
    ///
    /// Reads on a control contig are reported only in the per-contig table;
    /// they contribute to no genome-level metric. Interval restriction never
    /// applies to control contigs.
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
/// Reports how many usable CpG observations the library yielded and why the
/// rest were excluded (QC-fail, duplicate, MAPQ, end mask, base quality,
/// non-informative base, mate overlap), normalised per read and per Mb
/// sequenced. Also reports CpG depth statistics, per-context conversion and
/// methylation, per-contig metrics (spike-in controls included), M-bias, the
/// CpG methylation-level distribution and retained CpH per read. Only
/// cytosine positions are examined; whole-genome coverage is `wgs`'s job.
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
/// cytosine calls, and on leaving the contig folds the per-site depths into the genome-wide
/// statistics.
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
    // Reference / context-code buffer reused across contigs.
    contig_buffer: Vec<u8>,

    // Genome-wide accumulators (every non-control contig).
    reads: u64,
    bases_sequenced: u64,
    genome: Tally,
    cpg_territory: u64,
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
            contig_buffer: Vec::new(),
            reads: 0,
            bases_sequenced: 0,
            genome: Tally::default(),
            cpg_territory: 0,
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
        let is_control = self.is_control[ref_id];

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
            contig.tally.exclude_read(exclusion, potential);
            return Ok(());
        }

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

    /// Walk a passing read's aligned bases 64 positions at a time: each cytosine on the strand
    /// the read reports is tallied into the current contig's CpG funnel, conversion counts and
    /// per-site depths and (off control contigs) M-bias, and each CpH cytosine on the other strand
    /// into the directionality check. Cytosine calls are excluded in order: end mask, base
    /// quality, (CpG funnel only) non-informative base, then mate overlap. Returns the read's CpH
    /// evidence for the retained-CpH histogram.
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
        let ContigState { contexts, sites, tally, .. } = self.current.as_mut().unwrap();
        let quals = record.quality_scores();
        let bases = record.sequence();
        let strand = StrandCodes::new(read.top_strand);
        let mut evidence = ReadEvidence::default();
        let mut site_cursor = 0usize;

        for block in aligned_blocks(record) {
            let block_start = read.ref_start + block.ref_offset as usize;
            let read_offset = block.read_offset as usize;
            // A block past the contig end or a CIGAR longer than SEQ is malformed; the bases
            // beyond either are dropped.
            let block_end = (block_start + block.len as usize)
                .min(contexts.len())
                .min(block_start + bases.len().saturating_sub(read_offset));
            if block_start >= block_end {
                continue;
            }
            site_cursor = sites.seek(site_cursor, block_start.saturating_sub(1));

            for window_start in (block_start..block_end).step_by(64) {
                let len = (block_end - window_start).min(64);
                let read_pos = read_offset + (window_start - block_start);
                // A missing quality (QUAL `*`, or shorter than SEQ) fails the base-quality test.
                let window_quals =
                    quals.get(read_pos..).map_or(&[][..], |q| &q[..q.len().min(len)]);
                let window = CytosineWindow::classify(
                    &contexts[window_start..window_start + len],
                    &bases[read_pos..read_pos + len],
                    window_quals,
                    min_bq,
                    &strand,
                );
                // Positions fit in u32 for every supported reference.
                #[allow(clippy::cast_possible_truncation)]
                let window_ref_pos = window_start as u32;
                let unmasked = read.unmasked_bits(read_pos, len);
                let [cpa, cpc, cpg, cpt] = window.reported;
                let cph = cpa | cpc | cpt;
                let informative = window.retained | window.converted;

                if !is_control {
                    let calls = window.quality_ok & informative;
                    mbias.record_window(read, read_pos, cpg & calls, cph & calls, window.retained);
                }

                // Calls past the end mask and base quality, which the mate may already have made.
                let counted = (cpg | cph) & unmasked & window.quality_ok;
                let covered = action.mate_covered_bits(window_ref_pos) & counted;
                action.on_counted_bits(window_ref_pos, counted & !covered);
                let kept = counted & !covered;

                evidence.cph_sites += (cph & counted & informative).count_ones() as usize;
                evidence.retained_cph += (cph & counted & window.retained).count_ones() as usize;
                for (context, cytosines) in
                    [(Context::CpA, cpa), (Context::CpC, cpc), (Context::CpT, cpt)]
                {
                    tally.count_calls(context, cytosines & kept, &window);
                }

                let scoped = cpg & window.in_scope;
                let passed = scoped & counted;
                let good = passed & informative & !covered;
                let funnel = &mut tally.cpg;
                funnel.potential += u64::from(scoped.count_ones());
                funnel.end_mask += u64::from((scoped & !unmasked).count_ones());
                funnel.baseq += u64::from((scoped & unmasked & !window.quality_ok).count_ones());
                funnel.non_informative += u64::from((passed & !informative).count_ones());
                funnel.overlap += u64::from((passed & informative & covered).count_ones());
                funnel.good += u64::from(good.count_ones());
                tally.count_calls(Context::CpG, good, &window);
                for bit in set_bits(good) {
                    // A bottom-strand read reports a site's G, one past its C.
                    let c_pos = window_start + bit - usize::from(!read.top_strand);
                    site_cursor = sites.advance(site_cursor, c_pos);
                    let retained = window.retained >> bit & 1 == 1;
                    sites.count_call(site_cursor, retained);
                }

                let off_strand = window.off_strand_cph & window.quality_ok & unmasked;
                tally.off_strand_retained +=
                    u64::from((off_strand & window.off_strand_retained).count_ones());
                tally.off_strand_converted +=
                    u64::from((off_strand & window.off_strand_converted).count_ones());
            }
        }
        evidence
    }

    /// Start contig `ref_id`: load its reference, encode it as context codes, and index its CpG
    /// sites, marking which are in scope (on a control contig, all of them; otherwise those
    /// inside `-L`, or all without it).
    ///
    /// # Errors
    /// Returns an error if the contig cannot be read from the reference.
    fn begin_contig(&mut self, ref_id: usize) -> Result<()> {
        let mut contexts = self.load_contexts(ref_id)?;
        let sites = CpgSites::new(&contexts, self.scope_mask(ref_id));
        sites.mark_in_scope(&mut contexts);
        self.mate_buffer.clear();
        self.current = Some(ContigState { ref_id, contexts, sites, tally: Tally::default() });
        Ok(())
    }

    /// Load contig `ref_id`'s reference into the reused buffer and encode it in place as context
    /// codes.
    ///
    /// # Errors
    /// Returns an error if the contig cannot be read from the reference.
    fn load_contexts(&mut self, ref_id: usize) -> Result<Vec<u8>> {
        let mut buffer = std::mem::take(&mut self.contig_buffer);
        // `dict` is set in `initialize`, which always runs before contigs are loaded.
        let name = self.dict.as_ref().unwrap()[ref_id].name().to_string();
        self.reference.load_contig_into(&name, true, &mut buffer)?;
        encode_contexts(&mut buffer);
        Ok(buffer)
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
        }
        let row = self.contig_row(contig.ref_id, &contig.tally, contig.sites.n_in_scope);
        self.contig_rows.push((contig.ref_id, row));
        self.processed[contig.ref_id] = true;
        self.contig_buffer = contig.contexts;
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
    /// zero, and it gets an all-zero per-contig row.
    ///
    /// # Errors
    /// Returns an error if the contig cannot be read from the reference.
    fn account_unvisited_contig(&mut self, ref_id: usize) -> Result<()> {
        let is_control = self.is_control[ref_id];
        // A genome contig with no intervals under `-L` has no CpG territory; skip reading it.
        let skip_load =
            !is_control && self.intervals.as_ref().is_some_and(|ivl| !ivl.has_contig(ref_id));
        let n_in_scope = if skip_load {
            0
        } else {
            let mut sequence = std::mem::take(&mut self.contig_buffer);
            // `dict` is set in `initialize`, which always runs before `finish`.
            let name = self.dict.as_ref().unwrap()[ref_id].name().to_string();
            self.reference.load_contig_into(&name, true, &mut sequence)?;
            let n_in_scope = count_cpg_sites(&sequence, self.scope_mask(ref_id).as_ref());
            self.contig_buffer = sequence;
            n_in_scope
        };
        if !is_control {
            self.depth_histogram.add_n(0, n_in_scope);
            self.cpg_territory += n_in_scope;
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
            cpg_obs: tally.cpg.good,
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
            "methyl: cpg_territory={}, cpg_obs={}, metrics={}",
            self.cpg_territory,
            self.genome.cpg.good,
            self.metrics_path.display(),
        );
        Ok(())
    }

    /// Build the summary metrics row from the genome-wide accumulators.
    fn summary_metrics(&self) -> MethylMetrics {
        let cpg = &self.genome.cpg;
        let frac_excluded = |n: u64| safe_div_f(n as f64, cpg.potential as f64);
        let (mean, sd) = self.depth_histogram.mean_and_sd();
        let (median, mad) = self.depth_histogram.median_and_mad();
        let frac_at = |depth: usize| self.depth_histogram.frac_at_or_above(depth);
        let beta_sites: u64 = self.beta_histogram.iter().sum();

        MethylMetrics {
            sample: self.sample.clone(),
            chemistry: self.chemistry.to_string(),
            reads: self.reads,
            bases_sequenced: self.bases_sequenced,
            cpg_territory: self.cpg_territory,
            cpg_obs_potential: cpg.potential,
            cpg_obs: cpg.good,
            cpg_obs_per_read: safe_div_f(cpg.good as f64, self.reads as f64),
            cpg_obs_per_mb_sequenced: safe_div_f(
                cpg.good as f64,
                self.bases_sequenced as f64 / 1_000_000.0,
            ),
            mean_cpg_depth: mean,
            sd_cpg_depth: sd,
            median_cpg_depth: median,
            mad_cpg_depth: mad,
            frac_excluded_qc_fail: frac_excluded(cpg.qc_fail),
            frac_excluded_dupe: frac_excluded(cpg.dupe),
            frac_excluded_mapq: frac_excluded(cpg.mapq),
            frac_excluded_end_mask: frac_excluded(cpg.end_mask),
            frac_excluded_baseq: frac_excluded(cpg.baseq),
            frac_excluded_non_informative: frac_excluded(cpg.non_informative),
            frac_excluded_overlap: frac_excluded(cpg.overlap),
            frac_excluded_total: frac_excluded(cpg.potential - cpg.good),
            frac_cpgs_at_1x: frac_at(1),
            frac_cpgs_at_5x: frac_at(5),
            frac_cpgs_at_10x: frac_at(10),
            frac_cpgs_at_15x: frac_at(15),
            frac_cpgs_at_20x: frac_at(20),
            frac_cpgs_at_30x: frac_at(30),
            cpg_methylation: self.genome.methylation(Context::CpG, self.chemistry),
            cpg_conversion: self.genome.conversion(Context::CpG),
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
            "CpG Depth Distribution",
            Some(&format!("{} · {}", self.sample, self.filter_summary())),
            "CpG Depth (X)",
            "Fraction of CpGs",
            &self.coverage_plot_path,
        )
    }

    /// A one-line summary of the filters behind each CpG observation, for plot subtitles.
    fn filter_summary(&self) -> String {
        let mut parts = vec![format!("MAPQ ≥ {}", self.min_mapq), format!("BQ ≥ {}", self.min_bq)];
        if self.include_duplicates {
            parts.push("duplicates included".to_string());
        }
        let masks = self.masks;
        if masks.r1_5p + masks.r1_3p + masks.r2_5p + masks.r2_3p > 0 {
            parts.push(format!(
                "masks R1 {}/{}, R2 {}/{}",
                masks.r1_5p, masks.r1_3p, masks.r2_5p, masks.r2_3p
            ));
        }
        parts.push(format!("depth capped at {}", self.depth_histogram.counts().len() - 1));
        parts.join(", ")
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
            .with_title("M-bias")
            .with_subtitle(self.sample.clone())
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
            .with_title("CpG Methylation Levels")
            .with_subtitle(format!("{} · CpGs at depth ≥ {}", self.sample, self.min_beta_depth))
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
        self.bases_sequenced += record.sequence_len() as u64;

        let (Some(ref_id), Some(start)) =
            (record.reference_sequence_id(), record.alignment_start())
        else {
            return Ok(());
        };
        if flags.is_unmapped() {
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
    /// to `wgs`'s, whose per-read walk with mate-overlap deduplication is similar.
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
    /// Sum of the stored sequence lengths of primary records: the sequence data in the BAM.
    /// Soft-clipped bases count; hard-clipped bases do not.
    pub bases_sequenced: u64,
    /// Number of reference CpG sites in the territory (both strands of a CpG are one site).
    pub cpg_territory: u64,
    /// Read bases aligned to the cytosine of a territory CpG on the strand the read reports,
    /// before any filter.
    pub cpg_obs_potential: u64,
    /// Potential CpG observations that passed every filter and read C or T (G or A on the
    /// bottom strand).
    pub cpg_obs: u64,
    /// Good CpG observations per primary read.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_obs_per_read: f64,
    /// Good CpG observations per million bases of sequence data (`bases_sequenced`), not per
    /// megabase of genome.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub cpg_obs_per_mb_sequenced: f64,
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
    /// Fraction of potential CpG observations in QC-fail reads (e.g. unconverted templates
    /// flagged by methylsieve).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_qc_fail: f64,
    /// Fraction of potential CpG observations in duplicate reads.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_dupe: f64,
    /// Fraction of potential CpG observations in reads below the minimum mapping quality,
    /// assessed after QC-fail and duplicates. The loss among non-duplicate data is this
    /// divided by `1 - frac_excluded_qc_fail - frac_excluded_dupe`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_mapq: f64,
    /// Fraction of potential CpG observations excluded by the end masks (`--mask-*`).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_end_mask: f64,
    /// Fraction of potential CpG observations below the minimum base quality.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_baseq: f64,
    /// Fraction of potential CpG observations whose base was neither C nor T (G nor A on the
    /// bottom strand): a variant or sequencing error.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_non_informative: f64,
    /// Fraction of potential CpG observations at positions the overlapping mate already counted.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_overlap: f64,
    /// Fraction of potential CpG observations excluded for any reason.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_excluded_total: f64,
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
    /// Fraction of good territory CpG observations that were methylated.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_methylation: f64,
    /// Fraction of good territory CpG observations read as T (A on the bottom strand). Under
    /// em-seq / bisulfite this is `1 - cpg_methylation`; under taps / illumina-5base it equals
    /// `cpg_methylation`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub cpg_conversion: f64,
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
    /// converted, counting both mates of an overlapping pair. Near the sequencing error rate for
    /// a directional library; a high value means the library is not directional.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub off_strand_cph_conversion: f64,
    /// Mean methylation level (β, the methylated fraction of a CpG's good observations) of
    /// territory CpGs at depth ≥ `--min-beta-depth`.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub mean_cpg_beta: f64,
    /// Fraction of those CpGs with β < 0.2 (largely unmethylated).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_low: f64,
    /// Fraction of those CpGs with 0.2 ≤ β ≤ 0.8 (partially methylated).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_intermediate: f64,
    /// Fraction of those CpGs with β > 0.8 (largely methylated).
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_cpgs_beta_high: f64,
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
    pub cpg_obs: u64,
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
    /// One context code per reference position (see [`encode_contexts`]), with [`SCOPE_BIT`] set
    /// at both cytosines of each in-scope CpG.
    contexts: Vec<u8>,
    sites: CpgSites,
    tally: Tally,
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
    /// Index the CpG sites of a contig from its context codes: a site per top-strand CpG
    /// cytosine. With a `scope_mask` over the contig's positions, a site is in scope when its C
    /// is set in the mask.
    fn new(contexts: &[u8], scope_mask: Option<BitVec>) -> Self {
        // Positions fit in u32 for every supported reference.
        #[allow(clippy::cast_possible_truncation)]
        let positions: Vec<u32> = contexts
            .iter()
            .enumerate()
            .filter(|&(_, &code)| strand_context(code, true) == Some(Context::CpG))
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

    /// As [`Self::seek`], but scanning forward one site at a time: cheaper than a binary search
    /// when `pos` is at most a few sites past `from`, as successive CpGs within a read are.
    fn advance(&self, mut from: usize, pos: usize) -> usize {
        while from < self.positions.len() && (self.positions[from] as usize) < pos {
            from += 1;
        }
        from
    }

    /// Set [`SCOPE_BIT`] in `contexts` at both cytosines (the C, and the G that the bottom strand
    /// reads) of each in-scope site.
    fn mark_in_scope(&self, contexts: &mut [u8]) {
        for (site, &pos) in self.positions.iter().enumerate() {
            if self.in_scope(site) {
                contexts[pos as usize] |= SCOPE_BIT;
                contexts[pos as usize + 1] |= SCOPE_BIT;
            }
        }
    }

    /// Count a good call at site `site` that read unconverted (`retained`) or converted.
    fn count_call(&mut self, site: usize, retained: bool) {
        let counts = if retained { &mut self.retained } else { &mut self.converted };
        counts[site] = counts[site].saturating_add(1);
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
    /// Exclude a whole read: its potential CpG observations go to the read-level category.
    fn exclude_read(&mut self, exclusion: ReadExclusion, potential: u64) {
        self.cpg.potential += potential;
        match exclusion {
            ReadExclusion::QcFail => self.cpg.qc_fail += potential,
            ReadExclusion::Duplicate => self.cpg.dupe += potential,
            ReadExclusion::Mapq => self.cpg.mapq += potential,
        }
    }

    /// Count good calls in `context` at the positions of `window` set in `cytosines`.
    fn count_calls(&mut self, context: Context, cytosines: u64, window: &CytosineWindow) {
        self.retained[context as usize] += u64::from((cytosines & window.retained).count_ones());
        self.converted[context as usize] += u64::from((cytosines & window.converted).count_ones());
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
        self.cpg.add(&other.cpg);
        for context in 0..4 {
            self.retained[context] += other.retained[context];
            self.converted[context] += other.converted[context];
        }
        self.off_strand_retained += other.off_strand_retained;
        self.off_strand_converted += other.off_strand_converted;
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
    /// The stored-sequence positions outside both end masks: `[unmasked_start, unmasked_end)`.
    unmasked_start: usize,
    unmasked_end: usize,
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
        let leading_hard_clip = if reverse { clip_right } else { clip_left };
        let (mask_5p, mask_3p) = match role {
            ReadRole::R2 => (masks.r2_5p, masks.r2_3p),
            ReadRole::R1 | ReadRole::Fragment => (masks.r1_5p, masks.r1_3p),
        };

        // Unmasked cycles are `[mask_5p, sequenced_len - mask_3p)`; map them onto stored
        // positions, which run backwards through the cycles on a reverse-strand read.
        #[expect(clippy::cast_possible_wrap, reason = "read lengths are far below i64::MAX")]
        let (len, lead, first_cycle, end_cycle) = (
            sequence_len as i64,
            leading_hard_clip as i64,
            mask_5p as i64,
            (sequence_len + clip_left + clip_right) as i64 - mask_3p as i64,
        );
        let (start_pos, end_pos) = if reverse {
            (len + lead - end_cycle, len + lead - first_cycle)
        } else {
            (first_cycle - lead, end_cycle - lead)
        };
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "clamped to [0, sequence length]"
        )]
        let clamp = |pos: i64| pos.clamp(0, len) as usize;
        let unmasked_end = clamp(end_pos);
        Self {
            ref_start: start,
            top_strand,
            role,
            reverse,
            sequence_len,
            leading_hard_clip,
            unmasked_start: clamp(start_pos).min(unmasked_end),
            unmasked_end,
        }
    }

    /// The 0-based sequencing cycle of the base at `read_pos` in the stored sequence.
    fn cycle(&self, read_pos: usize) -> usize {
        let from_5p = if self.reverse { self.sequence_len - 1 - read_pos } else { read_pos };
        self.leading_hard_clip + from_5p
    }

    /// Of the `len` (≤ 64) stored positions from `read_pos`, those outside both end masks: bit
    /// `i` for `read_pos + i`.
    fn unmasked_bits(&self, read_pos: usize, len: usize) -> u64 {
        let start = self.unmasked_start.saturating_sub(read_pos).min(len);
        let end = self.unmasked_end.saturating_sub(read_pos).min(len);
        bit_range(start, end)
    }
}

/// The context codes and read bases that identify a read's calls, by the strand it reports.
struct StrandCodes {
    /// The context code (strand fields only) of a reported-strand cytosine in each [`Context`].
    reported: [u8; 4],
    /// The context code of an off-strand CpG cytosine.
    off_strand_cpg: u8,
    /// Upper-case read bases: unconverted then converted on the reported strand, then the same
    /// on the other strand.
    bases: [u8; 4],
}

impl StrandCodes {
    /// The codes for a read reporting the top strand (or, if `top_strand` is false, the bottom).
    fn new(top_strand: bool) -> Self {
        let (reported_shift, other_shift) =
            if top_strand { (0, BOTTOM_SHIFT) } else { (BOTTOM_SHIFT, 0) };
        let contexts = [Context::CpA, Context::CpC, Context::CpG, Context::CpT];
        Self {
            reported: contexts.map(|context| context.field() << reported_shift),
            off_strand_cpg: Context::CpG.field() << other_shift,
            bases: if top_strand { *b"CTGA" } else { *b"GACT" },
        }
    }
}

/// Up to 64 consecutive aligned positions of a read, classified against the reference: bit `i`
/// of each mask describes the window's `i`th position.
struct CytosineWindow {
    /// Cytosines on the strand the read reports, by [`Context`].
    reported: [u64; 4],
    /// CpH cytosines on the other strand.
    off_strand_cph: u64,
    /// Positions of in-scope CpGs.
    in_scope: u64,
    /// The read base is the reported strand's unconverted base (C, or G on the bottom strand).
    retained: u64,
    /// The read base is the reported strand's converted base (T, or A on the bottom strand).
    converted: u64,
    /// The read base is the other strand's unconverted base.
    off_strand_retained: u64,
    /// The read base is the other strand's converted base.
    off_strand_converted: u64,
    /// The base has a quality of at least the minimum.
    quality_ok: u64,
}

impl CytosineWindow {
    /// Classify the parallel `codes`, `bases` and `quals` of up to 64 aligned positions; a
    /// position past the end of `quals` fails the quality test.
    fn classify(
        codes: &[u8],
        bases: &[u8],
        quals: &[u8],
        min_bq: u8,
        strand: &StrandCodes,
    ) -> Self {
        let [cpa, cpc, cpg, cpt] = strand.reported;
        let [cpa, cpc, cpg, cpt, no_cytosine, cpg_on_other_strand] = simd::masked_eq_bitmasks(
            codes,
            STRAND_FIELDS_MASK,
            [cpa, cpc, cpg, cpt, 0, strand.off_strand_cpg],
        );
        let [in_scope] = simd::masked_eq_bitmasks(codes, SCOPE_BIT, [SCOPE_BIT]);
        // Clearing bit 5 upper-cases ASCII letters.
        let [retained, converted, off_strand_retained, off_strand_converted] =
            simd::masked_eq_bitmasks(bases, !0x20, strand.bases);
        // A position holds a cytosine on at most one strand, so any cytosine that is neither on
        // the reported strand nor an off-strand CpG is an off-strand CpH.
        let off_strand_cph = bit_range(0, codes.len())
            & !(no_cytosine | cpa | cpc | cpg | cpt | cpg_on_other_strand);
        Self {
            reported: [cpa, cpc, cpg, cpt],
            off_strand_cph,
            in_scope,
            retained,
            converted,
            off_strand_retained,
            off_strand_converted,
            quality_ok: simd::ge_bitmask(quals, min_bq),
        }
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

impl Context {
    /// The value a strand field of a context code holds for this context.
    fn field(self) -> u8 {
        self as u8 + 1
    }

    /// The context a strand field of a context code holds, or `None` for 0 (no cytosine).
    fn from_field(field: u8) -> Option<Context> {
        match field {
            1 => Some(Context::CpA),
            2 => Some(Context::CpC),
            3 => Some(Context::CpG),
            4 => Some(Context::CpT),
            _ => None,
        }
    }
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
    /// Record the calls `read` makes over a window of stored positions starting at `read_pos`
    /// (bit `i` for `read_pos + i`): `cpg` and `cph` mark the informative calls at CpG and CpH
    /// cytosines, and `retained` which of them read unconverted.
    fn record_window(
        &mut self,
        read: &ReadInfo,
        read_pos: usize,
        cpg: u64,
        cph: u64,
        retained: u64,
    ) {
        // Re-index the window by cycle, so bit `j` is the call at cycle `first_cycle + j`. A
        // reverse-strand read was sequenced from its last stored base, so its bits run backwards
        // and the window can start before cycle 0; the wrapping arithmetic still lands every set
        // bit on its true cycle.
        let (first_cycle, cpg, cph, retained) = if read.reverse {
            let last_cycle = read.cycle(read_pos);
            (
                last_cycle.wrapping_sub(63),
                cpg.reverse_bits(),
                cph.reverse_bits(),
                retained.reverse_bits(),
            )
        } else {
            (read.cycle(read_pos), cpg, cph, retained)
        };
        let in_range = if first_cycle > usize::MAX - 63 {
            u64::MAX
        } else {
            MAX_MBIAS_CYCLE
                .checked_sub(first_cycle)
                .map_or(0, |room| bit_range(0, (room + 1).min(64)))
        };

        for (is_cpg, calls) in [(true, cpg), (false, cph)] {
            let calls = calls & in_range;
            if calls == 0 {
                continue;
            }
            let slot = &mut self.slots[Self::slot(read.role, is_cpg)];
            let top_cycle = first_cycle.wrapping_add(63 - calls.leading_zeros() as usize);
            if top_cycle >= slot.len() {
                slot.resize(top_cycle + 1, CycleCounts::default());
            }
            for bit in set_bits(calls & retained) {
                slot[first_cycle.wrapping_add(bit)].retained += 1;
            }
            for bit in set_bits(calls & !retained) {
                slot[first_cycle.wrapping_add(bit)].converted += 1;
            }
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

/// Encode `sequence` (upper-cased reference) in place as one context code per position. The low
/// three bits hold the context of a top-strand cytosine (a reference C, by its next base), the
/// next three bits that of a bottom-strand cytosine (a reference G, by the complement of its
/// preceding base); each field is 0 when there is no cytosine on that strand or its neighbour is
/// not A, C, G or T. At most one field is set, so a code of 0 means no cytosine at all.
fn encode_contexts(sequence: &mut [u8]) {
    let mut previous = b'N';
    for pos in 0..sequence.len() {
        let current = sequence[pos];
        let next = sequence.get(pos + 1).copied().unwrap_or(b'N');
        // Branch-free: each table yields 0 unless `current` is the cytosine it applies to.
        sequence[pos] = (TOP_FIELD_BY_NEXT[usize::from(next)] * u8::from(current == b'C'))
            | (BOTTOM_FIELD_BY_PREVIOUS[usize::from(previous)] * u8::from(current == b'G'));
        previous = current;
    }
}

/// The top-strand field of a context code for a reference C, by the base that follows it.
const TOP_FIELD_BY_NEXT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'G' as usize] = 3;
    table[b'T' as usize] = 4;
    table
};

/// The bottom-strand field of a context code (already shifted) for a reference G, by the base
/// before it: the bottom strand reads its complement next.
const BOTTOM_FIELD_BY_PREVIOUS: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'T' as usize] = 1 << BOTTOM_SHIFT;
    table[b'G' as usize] = 2 << BOTTOM_SHIFT;
    table[b'C' as usize] = 3 << BOTTOM_SHIFT;
    table[b'A' as usize] = 4 << BOTTOM_SHIFT;
    table
};

/// The number of CpG sites in `sequence` (upper-cased reference), or of those whose C is set in
/// `scope_mask` when one is given.
fn count_cpg_sites(sequence: &[u8], scope_mask: Option<&BitVec>) -> u64 {
    let cpgs = sequence.windows(2).enumerate().filter(|(_, pair)| *pair == b"CG");
    match scope_mask {
        None => cpgs.count() as u64,
        Some(mask) => cpgs.filter(|&(pos, _)| mask.get(pos).is_some_and(|bit| *bit)).count() as u64,
    }
}

/// The indices of the set bits of `word`, ascending.
fn set_bits(mut word: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        if word == 0 {
            return None;
        }
        let bit = word.trailing_zeros() as usize;
        word &= word - 1;
        Some(bit)
    })
}

/// A mask with bits `[start, end)` set (`end` ≤ 64); empty when `start >= end`.
fn bit_range(start: usize, end: usize) -> u64 {
    if start >= end {
        return 0;
    }
    let below_end = if end >= 64 { u64::MAX } else { (1u64 << end) - 1 };
    below_end & !((1u64 << start) - 1)
}

/// The context of the cytosine a context code records on one strand, if any.
fn strand_context(code: u8, top_strand: bool) -> Option<Context> {
    let field =
        if top_strand { code & CONTEXT_MASK } else { (code >> BOTTOM_SHIFT) & CONTEXT_MASK };
    Context::from_field(field)
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

    /// The contexts `encode_contexts` records at each position of `sequence` on one strand.
    fn strand_contexts(sequence: &[u8], top_strand: bool) -> Vec<Option<Context>> {
        let mut codes = sequence.to_vec();
        encode_contexts(&mut codes);
        codes.iter().map(|&code| strand_context(code, top_strand)).collect()
    }

    #[test]
    fn top_strand_context_comes_from_the_next_reference_base() {
        let contexts = strand_contexts(b"CACCCGCTCN", true);
        assert_eq!(contexts[0], Some(Context::CpA));
        assert_eq!(contexts[2], Some(Context::CpC));
        assert_eq!(contexts[4], Some(Context::CpG));
        assert_eq!(contexts[6], Some(Context::CpT));
        assert_eq!(contexts[8], None);
        assert_eq!(contexts[1], None);
    }

    #[test]
    fn bottom_strand_context_comes_from_the_complement_of_the_previous_base() {
        let contexts = strand_contexts(b"TGGGCGAGNG", false);
        assert_eq!(contexts[1], Some(Context::CpA));
        assert_eq!(contexts[2], Some(Context::CpC));
        assert_eq!(contexts[5], Some(Context::CpG));
        assert_eq!(contexts[7], Some(Context::CpT));
        assert_eq!(contexts[9], None);
        assert_eq!(strand_contexts(b"GA", false)[0], None);
    }

    #[test]
    fn context_at_the_contig_end_is_unknown() {
        assert_eq!(strand_contexts(b"AC", true)[1], None);
    }

    #[test]
    fn count_cpg_sites_counts_cg_dinucleotides_in_scope() {
        assert_eq!(count_cpg_sites(b"ACGTTCGCG", None), 3);
        let mut mask = bitvec![0; 9];
        mask.set(5, true);
        assert_eq!(count_cpg_sites(b"ACGTTCGCG", Some(&mask)), 1);
    }

    #[test]
    fn a_position_records_a_cytosine_on_at_most_one_strand() {
        let mut codes = b"ACGTNacgt".to_ascii_uppercase();
        encode_contexts(&mut codes);
        for code in codes {
            assert!(strand_context(code, true).is_none() || strand_context(code, false).is_none());
        }
    }

    fn encoded(sequence: &[u8]) -> Vec<u8> {
        let mut codes = sequence.to_vec();
        encode_contexts(&mut codes);
        codes
    }

    #[test]
    fn cpg_sites_index_every_cg_dinucleotide() {
        let sites = CpgSites::new(&encoded(b"ACGTTCGCG"), None);
        assert_eq!(sites.positions, vec![1, 5, 7]);
        assert_eq!(sites.n_in_scope, 3);
    }

    #[test]
    fn cpg_sites_in_scope_follow_the_mask_at_the_c() {
        let mut mask = bitvec![0; 9];
        mask.set(5, true);
        mask.set(8, true);
        let sites = CpgSites::new(&encoded(b"ACGTTCGCG"), Some(mask));
        assert_eq!(sites.n_in_scope, 1);
        assert!(!sites.in_scope(0));
        assert!(sites.in_scope(1));
        assert!(!sites.in_scope(2));
    }

    #[test]
    fn advance_agrees_with_seek() {
        let sites = CpgSites::new(&encoded(b"ACGTTCGCG"), None);
        for from in 0..3 {
            for pos in 0..10 {
                if sites.positions.get(from).is_none_or(|&site| site as usize <= pos) {
                    assert_eq!(
                        sites.advance(from, pos),
                        sites.seek(from, pos),
                        "from {from} pos {pos}"
                    );
                }
            }
        }
    }

    #[test]
    fn seek_finds_the_first_site_at_or_after_a_position() {
        let sites = CpgSites::new(&encoded(b"ACGTTCGCG"), None);
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
    fn end_masks_exclude_positions_from_each_end() {
        let masks = EndMasks { r1_5p: 2, r1_3p: 1, r2_5p: 0, r2_3p: 0 };
        let record = crate::test_support::read().cigar("10M").into_riker_record();
        let info = ReadInfo::new(&record, 0, true, masks);
        assert_eq!(info.unmasked_bits(0, 10), 0b01_1111_1100);
        assert_eq!(info.unmasked_bits(4, 6), 0b01_1111);
    }

    #[test]
    fn unmasked_positions_are_those_whose_cycles_avoid_both_end_masks() {
        for cigar in ["10M", "3H10M", "10M2H", "3H10M2H"] {
            for reverse in [false, true] {
                for (mask_5p, mask_3p) in [(0, 0), (1, 0), (0, 2), (4, 3), (9, 9), (20, 0)] {
                    let masks = EndMasks { r1_5p: mask_5p, r1_3p: mask_3p, r2_5p: 0, r2_3p: 0 };
                    let mut builder = crate::test_support::read().cigar(cigar);
                    if reverse {
                        builder = builder.reverse();
                    }
                    let info = ReadInfo::new(&builder.into_riker_record(), 0, true, masks);
                    let sequenced_len = 10 + hard_clip_total(cigar);
                    for read_pos in 0..10 {
                        let cycle = info.cycle(read_pos);
                        let masked = cycle < mask_5p || sequenced_len - 1 - cycle < mask_3p;
                        let unmasked = info.unmasked_bits(read_pos, 1) == 1;
                        assert_eq!(
                            unmasked, !masked,
                            "{cigar} reverse={reverse} masks={mask_5p}/{mask_3p} pos={read_pos}"
                        );
                    }
                }
            }
        }
    }

    /// The total hard-clipped length in a CIGAR of `<n>H` and `M` operations.
    fn hard_clip_total(cigar: &str) -> usize {
        cigar
            .split_inclusive(['H', 'M'])
            .filter_map(|op| op.strip_suffix('H'))
            .map(|n| n.parse::<usize>().unwrap())
            .sum()
    }

    #[test]
    fn record_window_files_forward_read_calls_by_cycle() {
        let masks = EndMasks { r1_5p: 0, r1_3p: 0, r2_5p: 0, r2_3p: 0 };
        let forward = crate::test_support::read().cigar("8H10M").into_riker_record();
        let info = ReadInfo::new(&forward, 0, true, masks);
        let mut mbias = MbiasCounts::default();
        mbias.record_window(&info, 2, 0b1, 0b100, 0b1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, true)[10].retained, 1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, false)[12].converted, 1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, false).len(), 13);
    }

    #[test]
    fn record_window_files_reverse_read_calls_by_cycle() {
        let masks = EndMasks { r1_5p: 0, r1_3p: 0, r2_5p: 0, r2_3p: 0 };
        let reverse = crate::test_support::read().cigar("10M8H").reverse().into_riker_record();
        let info = ReadInfo::new(&reverse, 0, false, masks);
        let mut mbias = MbiasCounts::default();
        // Stored position 0 is cycle 17 and position 9 is cycle 8.
        mbias.record_window(&info, 0, 0b1, 1 << 9, 0b1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, true)[17].retained, 1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, false)[8].converted, 1);
        assert_eq!(mbias.cycles(ReadRole::Fragment, false).len(), 9);
    }

    /// The cycles with a retained call in the single-end CpG M-bias slot.
    fn retained_cycles(mbias: &MbiasCounts) -> Vec<usize> {
        let cycles = mbias.cycles(ReadRole::Fragment, true);
        (0..cycles.len()).filter(|&cycle| cycles[cycle].retained > 0).collect()
    }

    #[test]
    fn record_window_ignores_forward_read_cycles_past_the_limit() {
        let masks = EndMasks { r1_5p: 0, r1_3p: 0, r2_5p: 0, r2_3p: 0 };
        let clip = MAX_MBIAS_CYCLE - 5;
        let forward = crate::test_support::read().cigar(format!("{clip}H10M")).into_riker_record();
        let info = ReadInfo::new(&forward, 0, true, masks);
        let mut mbias = MbiasCounts::default();
        mbias.record_window(&info, 0, 0b11_1111_1111, 0, 0b11_1111_1111);
        assert_eq!(retained_cycles(&mbias), (clip..=MAX_MBIAS_CYCLE).collect::<Vec<_>>());
    }

    #[test]
    fn record_window_ignores_reverse_read_cycles_past_the_limit() {
        let masks = EndMasks { r1_5p: 0, r1_3p: 0, r2_5p: 0, r2_3p: 0 };
        let clip = MAX_MBIAS_CYCLE - 5;
        let reverse =
            crate::test_support::read().cigar(format!("10M{clip}H")).reverse().into_riker_record();
        let info = ReadInfo::new(&reverse, 0, true, masks);
        let mut mbias = MbiasCounts::default();
        mbias.record_window(&info, 0, 0b11_1111_1111, 0, 0b11_1111_1111);
        assert_eq!(retained_cycles(&mbias), (clip..=MAX_MBIAS_CYCLE).collect::<Vec<_>>());
    }

    #[test]
    fn set_bits_visits_each_set_bit_in_ascending_order() {
        assert_eq!(set_bits(0b1001 | 1 << 63).collect::<Vec<_>>(), vec![0, 3, 63]);
        assert_eq!(set_bits(0).count(), 0);
    }

    #[test]
    fn bit_range_sets_bits_from_start_up_to_end() {
        assert_eq!(bit_range(2, 5), 0b1_1100);
        assert_eq!(bit_range(0, 64), u64::MAX);
        assert_eq!(bit_range(63, 64), 1 << 63);
        assert_eq!(bit_range(5, 5), 0);
        assert_eq!(bit_range(6, 2), 0);
    }

    /// Classify `reference` (encoded, every CpG in scope) against `bases` at quality 30.
    fn classify(reference: &[u8], bases: &[u8], top_strand: bool) -> CytosineWindow {
        let mut codes = encoded(reference);
        CpgSites::new(&codes, None).mark_in_scope(&mut codes);
        let quals = vec![30; bases.len()];
        CytosineWindow::classify(&codes, bases, &quals, 20, &StrandCodes::new(top_strand))
    }

    #[test]
    fn classify_marks_top_strand_cytosines_by_context() {
        // C at 1 (CpG) and 6 (CpA) on the top strand; G at 2 (CpG) and 5 (CpT) on the bottom.
        let window = classify(b"ACGTAGCA", b"ACGTAGCA", true);
        assert_eq!(window.reported, [1 << 6, 0, 1 << 1, 0]);
        assert_eq!(window.off_strand_cph, 1 << 5);
    }

    #[test]
    fn classify_marks_bottom_strand_cytosines_by_context() {
        let window = classify(b"ACGTAGCA", b"ACGTAGCA", false);
        assert_eq!(window.reported, [0, 0, 1 << 2, 1 << 5]);
        assert_eq!(window.off_strand_cph, 1 << 6);
    }

    #[test]
    fn classify_reads_c_and_t_on_top_and_g_and_a_on_bottom_in_either_case() {
        let top = classify(b"ACGTAGCA", b"AcGTAGtA", true);
        assert_eq!(top.retained & top.reported[Context::CpG as usize], 1 << 1);
        assert_eq!(top.converted & top.reported[Context::CpA as usize], 1 << 6);
        let bottom = classify(b"ACGTAGCA", b"ACgTAaCA", false);
        assert_eq!(bottom.retained & bottom.reported[Context::CpG as usize], 1 << 2);
        assert_eq!(bottom.converted & bottom.reported[Context::CpT as usize], 1 << 5);
    }

    #[test]
    fn classify_reads_the_off_strand_bases_of_the_other_strand() {
        let top_read = classify(b"ACGTAGCA", b"ACGTAACA", true);
        assert_eq!(top_read.off_strand_cph & top_read.off_strand_converted, 1 << 5);
        let bottom_read = classify(b"ACGTAGCA", b"ACGTAGCA", false);
        assert_eq!(bottom_read.off_strand_cph & bottom_read.off_strand_retained, 1 << 6);
    }

    #[test]
    fn classify_fails_quality_where_quals_run_out() {
        let codes = encoded(b"ACGTAGCA");
        let window = CytosineWindow::classify(
            &codes,
            b"ACGTAGCA",
            &[30, 10, 30],
            20,
            &StrandCodes::new(true),
        );
        assert_eq!(window.quality_ok, 0b101);
    }

    #[test]
    fn mark_in_scope_flags_both_cytosines_of_each_scoped_site() {
        let mut codes = encoded(b"ACGTTCGCG");
        let mut mask = bitvec![0; 9];
        mask.set(5, true);
        CpgSites::new(&codes, Some(mask)).mark_in_scope(&mut codes);
        let scoped: Vec<usize> = (0..9).filter(|&pos| codes[pos] & SCOPE_BIT != 0).collect();
        assert_eq!(scoped, vec![5, 6]);
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
