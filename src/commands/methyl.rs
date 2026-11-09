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

    /// Walk a passing read's aligned bases, examining only cytosines: each cytosine on the strand
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
        let ContigState { contexts, sites, tally, top_cytosines, bottom_cytosines, .. } =
            self.current.as_mut().unwrap();
        let (reported, other) = if read.top_strand {
            (&*top_cytosines, &*bottom_cytosines)
        } else {
            (&*bottom_cytosines, &*top_cytosines)
        };
        let quals = record.quality_scores();
        let bases = record.sequence();
        let quality_ok = |read_pos: usize| quals.get(read_pos).is_some_and(|&qual| qual >= min_bq);
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

            // Cytosines on the strand the read reports, visited via the bitset so the bases in
            // between cost nothing.
            for ref_pos in SetBits::new(reported.as_raw_slice(), block_start, block_end) {
                let read_pos = read_offset + (ref_pos - block_start);
                let Some(context) = strand_context(contexts[ref_pos], read.top_strand) else {
                    continue;
                };
                let cycle = read.cycle(read_pos);
                let call = call_base(bases[read_pos], read.top_strand);
                let quality_ok = quality_ok(read_pos);

                if quality_ok && !is_control && call != Call::Other {
                    mbias.record(read.role, context == Context::CpG, cycle, call == Call::Retained);
                }

                // The site index of a CpG cytosine, when that site is in scope.
                let cpg_site = if context == Context::CpG {
                    let c_pos = if read.top_strand { ref_pos } else { ref_pos - 1 };
                    site_cursor = sites.advance(site_cursor, c_pos);
                    sites.in_scope(site_cursor).then_some(site_cursor)
                } else {
                    None
                };

                if read.is_masked(cycle) {
                    if cpg_site.is_some() {
                        tally.cpg.potential += 1;
                        tally.cpg.end_mask += 1;
                    }
                    continue;
                }
                if !quality_ok {
                    if cpg_site.is_some() {
                        tally.cpg.potential += 1;
                        tally.cpg.baseq += 1;
                    }
                    continue;
                }

                if context != Context::CpG && call != Call::Other {
                    evidence.cph_sites += 1;
                    if call == Call::Retained {
                        evidence.retained_cph += 1;
                    }
                }

                // Positions fit in u32 for every supported reference.
                #[allow(clippy::cast_possible_truncation)]
                let pos = ref_pos as u32;
                let covered = action.is_mate_covered(pos);
                if !covered {
                    action.on_counted(pos);
                }

                if context == Context::CpG {
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
                } else if call != Call::Other && !covered {
                    tally.count_call(context, call);
                }
            }

            // CpH cytosines on the other strand: the directionality check.
            for ref_pos in SetBits::new(other.as_raw_slice(), block_start, block_end) {
                let read_pos = read_offset + (ref_pos - block_start);
                if strand_context(contexts[ref_pos], !read.top_strand) == Some(Context::CpG)
                    || !quality_ok(read_pos)
                    || read.is_masked(read.cycle(read_pos))
                {
                    continue;
                }
                let call = call_base(bases[read_pos], !read.top_strand);
                if call != Call::Other {
                    tally.count_off_strand_call(call);
                }
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
        let contexts = self.load_contexts(ref_id)?;
        let sites = CpgSites::new(&contexts, self.scope_mask(ref_id));
        let top_cytosines = cytosine_bits(&contexts, true);
        let bottom_cytosines = cytosine_bits(&contexts, false);
        self.mate_buffer.clear();
        self.current = Some(ContigState {
            ref_id,
            contexts,
            top_cytosines,
            bottom_cytosines,
            sites,
            tally: Tally::default(),
        });
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
    /// One context code per reference position (see [`encode_contexts`]).
    contexts: Vec<u8>,
    /// Positions holding a top-strand cytosine (a reference C with a known context).
    top_cytosines: BitVec<u64, Lsb0>,
    /// Positions holding a bottom-strand cytosine (a reference G with a known context).
    bottom_cytosines: BitVec<u64, Lsb0>,
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

impl Context {
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

/// The set bits of a bitset stored as `u64` words (bit `i` of word `w` is position `64w + i`)
/// within `[start, end)`, in ascending order, scanning a word at a time.
struct SetBits<'a> {
    words: &'a [u64],
    /// Index of the word `current` came from.
    word_index: usize,
    /// Index of the word holding position `end - 1`.
    last_word: usize,
    /// Mask of the positions below `end` in the last word.
    last_mask: u64,
    /// Unvisited set bits of the current word.
    current: u64,
}

impl<'a> SetBits<'a> {
    /// The set bits of `words` within `[start, end)`; empty when `start >= end`.
    fn new(words: &'a [u64], start: usize, end: usize) -> Self {
        if start >= end {
            return Self { words, word_index: 0, last_word: 0, last_mask: 0, current: 0 };
        }
        let word_index = start / 64;
        let last_word = (end - 1) / 64;
        let last_mask = u64::MAX >> (63 - (end - 1) % 64);
        let mut current = words[word_index] & (u64::MAX << (start % 64));
        if word_index == last_word {
            current &= last_mask;
        }
        Self { words, word_index, last_word, last_mask, current }
    }
}

impl Iterator for SetBits<'_> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        while self.current == 0 {
            if self.word_index >= self.last_word {
                return None;
            }
            self.word_index += 1;
            self.current = self.words[self.word_index];
            if self.word_index == self.last_word {
                self.current &= self.last_mask;
            }
        }
        let bit = self.current.trailing_zeros() as usize;
        self.current &= self.current - 1;
        Some(self.word_index * 64 + bit)
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

/// A bit per position of `contexts`, set where it records a cytosine on the top (or bottom)
/// strand.
fn cytosine_bits(contexts: &[u8], top_strand: bool) -> BitVec<u64, Lsb0> {
    let shift = if top_strand { 0 } else { BOTTOM_SHIFT };
    let mut bits = bitvec![u64, Lsb0; 0; contexts.len()];
    for (word, chunk) in bits.as_raw_mut_slice().iter_mut().zip(contexts.chunks(64)) {
        *word = chunk.iter().enumerate().fold(0u64, |word, (bit, &code)| {
            word | (u64::from((code >> shift) & CONTEXT_MASK != 0) << bit)
        });
    }
    bits
}

/// The context of the cytosine a context code records on one strand, if any.
fn strand_context(code: u8, top_strand: bool) -> Option<Context> {
    let field =
        if top_strand { code & CONTEXT_MASK } else { (code >> BOTTOM_SHIFT) & CONTEXT_MASK };
    Context::from_field(field)
}

/// Classify a read base at a cytosine on the reported strand.
fn call_base(base: u8, top_strand: bool) -> Call {
    match (base.to_ascii_uppercase(), top_strand) {
        (b'C', true) | (b'G', false) => Call::Retained,
        (b'T', true) | (b'A', false) => Call::Converted,
        _ => Call::Other,
    }
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
    fn cytosine_bits_mark_each_strands_cytosines() {
        let codes = encoded(b"ACGTAGCA");
        let top: Vec<usize> = cytosine_bits(&codes, true).iter_ones().collect();
        let bottom: Vec<usize> = cytosine_bits(&codes, false).iter_ones().collect();
        assert_eq!(top, vec![1, 6]);
        assert_eq!(bottom, vec![2, 5]);
    }

    #[test]
    fn set_bits_visits_bits_within_the_range_across_words() {
        let words = [0b1001u64, 0, 1u64 << 63 | 1];
        let all: Vec<usize> = SetBits::new(&words, 0, 192).collect();
        assert_eq!(all, vec![0, 3, 128, 191]);
        let middle: Vec<usize> = SetBits::new(&words, 1, 191).collect();
        assert_eq!(middle, vec![3, 128]);
        assert_eq!(SetBits::new(&words, 4, 4).count(), 0);
        assert_eq!(SetBits::new(&words, 3, 4).collect::<Vec<_>>(), vec![3]);
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

    #[test]
    fn call_base_reads_c_and_t_on_top_and_g_and_a_on_bottom() {
        assert_eq!(call_base(b'C', true), Call::Retained);
        assert_eq!(call_base(b'T', true), Call::Converted);
        assert_eq!(call_base(b'A', true), Call::Other);
        assert_eq!(call_base(b'G', false), Call::Retained);
        assert_eq!(call_base(b'A', false), Call::Converted);
        assert_eq!(call_base(b'T', false), Call::Other);
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
