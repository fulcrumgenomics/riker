use std::path::PathBuf;

use anyhow::Result;
use clap::Args;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};

use crate::collector::Collector;
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::counter::Counter;
use crate::fasta::Fasta;
use crate::math::{safe_div, safe_div_f};
use crate::metrics::write_tsv;
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::pair_orientation::{PairOrientation, get_pair_orientation};
use crate::sam::record_utils::{derive_sample, get_integer_tag};
use crate::simd;

/// File suffix appended to the output prefix for alignment summary metrics.
pub const METRICS_SUFFIX: &str = ".alignment-metrics.txt";

/// Base quality threshold for HQ Q20 base counting.
const BASE_QUALITY_THRESHOLD: u8 = 20;

// ─── Options ──────────────────────────────────────────────────────────────────

/// Tool-specific tuning options for the alignment collector.
#[riker_derive::multi_options("aln", "Alignment Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct AlignmentOptions {
    /// Maximum insert size before a pair is considered chimeric.
    #[arg(long, default_value_t = AlignmentOptions::DEFAULT_MAX_INSERT_SIZE)]
    pub max_insert_size: u32,
    /// Minimum mapping quality for a read to be counted as high-quality (HQ).
    #[arg(long, default_value_t = AlignmentOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,
}

impl AlignmentOptions {
    const DEFAULT_MAX_INSERT_SIZE: u32 = 100_000;
    const DEFAULT_MIN_MAPQ: u8 = 20;
}

impl Default for AlignmentOptions {
    fn default() -> Self {
        Self { max_insert_size: Self::DEFAULT_MAX_INSERT_SIZE, min_mapq: Self::DEFAULT_MIN_MAPQ }
    }
}

// ─── CLI struct ───────────────────────────────────────────────────────────────

/// Collect alignment summary metrics from a BAM file.
///
/// Counts total, aligned, and high-quality reads broken down by read1,
/// read2, and pair categories. Computes mismatch rates, strand balance,
/// chimera rates, and adapter contamination when a reference is provided.
/// Metrics are written to <prefix>.alignment-metrics.txt.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker alignment -i input.bam -o out_prefix
  riker alignment -i input.bam -o out_prefix -R ref.fa"
)]
pub struct Alignment {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,

    #[command(flatten)]
    pub options: AlignmentOptions,
}

impl Command for Alignment {
    /// # Errors
    /// Returns an error if the BAM cannot be read or the output cannot be written.
    fn execute(&self) -> Result<()> {
        let (mut reader, header) =
            AlignmentReader::new(&self.input.input, self.reference.reference.as_deref())?;

        let mut collector = AlignmentCollector::new(
            &self.input.input,
            &self.output.output,
            self.reference.reference.clone(),
            &self.options,
        );

        collector.initialize(&header)?;
        let mut progress = ProgressLogger::new("alignment", "reads", 5_000_000);
        reader.for_each_record(&header, |record| {
            progress.record_with(record, &header);
            collector.accept(record, &header)
        })?;
        progress.finish();
        collector.finish()
    }
}

// ─── Collector ────────────────────────────────────────────────────────────────

/// Collects alignment summary metrics from a BAM file.
pub struct AlignmentCollector {
    input: PathBuf,
    metrics_path: PathBuf,
    reference_path: Option<PathBuf>,
    max_insert_size: u32,
    min_mapq: u8,
    sample: String,
    first: CategoryAccumulator,
    second: CategoryAccumulator,
}

impl AlignmentCollector {
    /// Create a new collector.  Output path is derived from `prefix` by appending
    /// [`METRICS_SUFFIX`].
    #[must_use]
    pub fn new(
        input: &std::path::Path,
        prefix: &std::path::Path,
        reference_path: Option<PathBuf>,
        options: &AlignmentOptions,
    ) -> Self {
        let metrics_path = super::command::output_path(prefix, METRICS_SUFFIX);
        Self {
            input: input.to_path_buf(),
            metrics_path,
            reference_path,
            max_insert_size: options.max_insert_size,
            min_mapq: options.min_mapq,
            sample: String::new(),
            first: CategoryAccumulator::new("read1"),
            second: CategoryAccumulator::new("read2"),
        }
    }
}

impl Collector for AlignmentCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        if let Some(ref ref_path) = self.reference_path {
            let fasta = Fasta::from_path(ref_path)?;
            fasta.validate_bam_header(header)?;
        }

        self.sample = derive_sample(&self.input, header);
        Ok(())
    }

    fn accept(&mut self, record: &RecordBuf, _header: &Header) -> Result<()> {
        let flags = record.flags();

        // Skip secondary and supplementary alignments entirely.
        if flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }

        // Route to per-category accumulator.
        if flags.is_segmented() {
            if flags.is_first_segment() {
                self.first.process_record(record, self.min_mapq, self.max_insert_size);
            } else {
                self.second.process_record(record, self.min_mapq, self.max_insert_size);
            }
        } else {
            self.first.process_record(record, self.min_mapq, self.max_insert_size);
        }

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        let mut metrics: Vec<AlignmentSummaryMetric> = Vec::new();

        if self.second.total_reads > 0 {
            let mut first_metric = self.first.compute_metric();
            let mut second_metric = self.second.compute_metric();
            let combined = self.first.merge(&self.second);
            let mut pair_metric = combined.compute_metric();

            // Override pair bad_cycles = read1 + read2 (Picard convention).
            pair_metric.bad_cycles = first_metric.bad_cycles + second_metric.bad_cycles;

            first_metric.sample.clone_from(&self.sample);
            second_metric.sample.clone_from(&self.sample);
            pair_metric.sample.clone_from(&self.sample);

            metrics.push(first_metric);
            metrics.push(second_metric);
            metrics.push(pair_metric);
        } else {
            let mut m = self.first.compute_metric();
            m.sample.clone_from(&self.sample);
            metrics.push(m);
        }

        write_tsv(&self.metrics_path, &metrics)?;
        Ok(())
    }

    fn name(&self) -> &'static str {
        "alignment"
    }
}

// ─── Metric struct ────────────────────────────────────────────────────────────

/// Alignment summary metrics per read category (read1, read2, pair).
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct AlignmentSummaryMetric {
    /// Sample name derived from the BAM read group SM tag or filename.
    pub sample: String,
    /// Read category: read1, read2, or pair.
    pub category: String,
    /// Total number of reads (including QC-failed reads).
    pub total_reads: u64,
    /// Number of PF reads that aligned to the reference.
    pub aligned_reads: u64,
    /// Fraction of PF reads that aligned to the reference.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_aligned: f64,
    /// Number of aligned reads with mapping quality >= min_mapq.
    pub hq_aligned_reads: u64,
    /// Total aligned bases from high-quality reads.
    pub hq_aligned_bases: u64,
    /// Number of aligned bases with base quality >= 20 from high-quality reads.
    pub hq_aligned_q20_bases: u64,
    /// Median number of mismatches per high-quality aligned read.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub hq_median_mismatches: f64,
    /// Rate of mismatches per aligned base across all PF aligned reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_6dp")]
    pub mismatch_rate: f64,
    /// Rate of mismatches per aligned base for high-quality reads only.
    #[serde(serialize_with = "crate::metrics::serialize_f64_6dp")]
    pub hq_mismatch_rate: f64,
    /// Rate of insertion and deletion events per aligned base.
    #[serde(serialize_with = "crate::metrics::serialize_f64_6dp")]
    pub indel_rate: f64,
    /// Mean read length across all PF reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_read_length: f64,
    /// Standard deviation of read length across all PF reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub sd_read_length: f64,
    /// Median read length across all PF reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub median_read_length: f64,
    /// Median absolute deviation of read length across all PF reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mad_read_length: f64,
    /// Minimum read length observed.
    pub min_read_length: u64,
    /// Maximum read length observed.
    pub max_read_length: u64,
    /// Mean aligned read length (M, I, =, X bases) for reads that aligned.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_aligned_read_length: f64,
    /// Number of aligned reads that are part of a mapped pair.
    pub aligned_reads_in_pairs: u64,
    /// Fraction of aligned reads that are part of a mapped pair.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_aligned_in_pairs: f64,
    /// Number of reads that are aligned and paired but not properly paired.
    pub reads_improperly_paired: u64,
    /// Fraction of aligned reads that are improperly paired.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_reads_improperly_paired: f64,
    /// Number of sequencing cycles where >= 80% of reads had an N base call.
    pub bad_cycles: u64,
    /// Fraction of aligned reads on the positive strand.
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub strand_balance: f64,
    /// Fraction of read pairs that are chimeric (different contig, large insert, or unexpected orientation).
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_chimeras: f64,
    /// Fraction of total bases that are soft-clipped (mapped reads only).
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_softclipped_reads: f64,
    /// Fraction of total bases that are hard-clipped (all PF reads).
    #[serde(serialize_with = "crate::metrics::serialize_f64_5dp")]
    pub frac_hardclipped_reads: f64,
    /// Mean number of soft-clipped bases at the 3-prime end of reads.
    #[serde(serialize_with = "crate::metrics::serialize_f64_2dp")]
    pub mean_3prime_softclipped_bases: f64,
}

// ─── Internal accumulator ─────────────────────────────────────────────────────

/// Raw counters for one read category (`read1` / `read2`).
struct CategoryAccumulator {
    category: &'static str,
    // Read counts
    total_reads: u64,
    pf_reads: u64,
    // Alignment counts (PF, primary, non-supplementary, mapped reads)
    aligned_reads: u64,
    pf_aligned_bases: u64, // M/=/X bases across all PF aligned reads
    hq_aligned_reads: u64,
    hq_aligned_bases: u64,
    hq_aligned_q20_bases: u64,
    aligned_reads_in_pairs: u64,
    reads_improperly_paired: u64,
    num_positive_strand: u64,
    // Chimera tracking
    chimeras: u64,
    chimeras_denominator: u64,
    // Clip tracking
    num_soft_clipped: u64, // S bases from mapped reads
    num_hard_clipped: u64, // H bases from all PF reads
    num_3prime_soft_clipped_bases: u64,
    num_reads_with_3prime_softclips: u64,
    // Indel events (not bases)
    indels: u64,
    // Mismatch tracking via NM tag
    total_nm: u64, // sum of NM for all PF aligned reads
    hq_nm: u64,    // sum of NM for HQ aligned reads
    hq_nm_histogram: Counter<u64>,
    // Read length tracking
    read_length_histogram: Counter<u64>,
    read_length_sum: u64,
    read_length_sum_sq: u64,
    aligned_read_length_sum: u64,
    // Bad cycle tracking: cycle number -> no-call count
    bad_cycle_nocalls: Counter<u64>,
}

impl CategoryAccumulator {
    fn new(category: &'static str) -> Self {
        Self {
            category,
            total_reads: 0,
            pf_reads: 0,
            aligned_reads: 0,
            pf_aligned_bases: 0,
            hq_aligned_reads: 0,
            hq_aligned_bases: 0,
            hq_aligned_q20_bases: 0,
            aligned_reads_in_pairs: 0,
            reads_improperly_paired: 0,
            num_positive_strand: 0,
            chimeras: 0,
            chimeras_denominator: 0,
            num_soft_clipped: 0,
            num_hard_clipped: 0,
            num_3prime_soft_clipped_bases: 0,
            num_reads_with_3prime_softclips: 0,
            indels: 0,
            total_nm: 0,
            hq_nm: 0,
            hq_nm_histogram: Counter::new(),
            read_length_histogram: Counter::new(),
            read_length_sum: 0,
            read_length_sum_sq: 0,
            aligned_read_length_sum: 0,
            bad_cycle_nocalls: Counter::new(),
        }
    }

    /// Process one PF, non-secondary, non-supplementary read.
    fn process_record(&mut self, record: &RecordBuf, min_mapq: u8, max_insert_size: u32) {
        let flags = record.flags();

        self.total_reads += 1;

        // Non-PF reads only count toward total_reads.
        if flags.is_qc_fail() {
            return;
        }

        self.pf_reads += 1;

        // Read length (sequence stored in the BAM — includes soft-clipped, excludes hard-clipped).
        let seq = record.sequence();
        let read_len = seq.len() as u64;

        self.read_length_histogram.count(read_len);
        self.read_length_sum += read_len;
        self.read_length_sum_sq += read_len * read_len;

        // Hard-clipped bases from all PF reads.
        self.num_hard_clipped += sum_cigar_op(record, Kind::HardClip);

        // Bad cycles: track per-cycle no-call (N base) counts.
        track_bad_cycles(
            &mut self.bad_cycle_nocalls,
            seq.as_ref(),
            read_len,
            flags.is_reverse_complemented(),
        );

        // --- Mapped reads only beyond this point ---
        if flags.is_unmapped() {
            return;
        }

        let mapq = record.mapping_quality().map_or(255u8, u8::from);
        let is_hq = mapq >= min_mapq;
        let negative_strand = flags.is_reverse_complemented();

        // Single CIGAR pass: computes soft-clip, aligned length, indel, and
        // per-base quality stats all at once.
        let cs = self.process_cigar(record, is_hq, negative_strand);

        self.num_soft_clipped += cs.soft_clip_bases;

        if cs.three_prime_soft_clip > 0 {
            self.num_3prime_soft_clipped_bases += cs.three_prime_soft_clip;
            self.num_reads_with_3prime_softclips += 1;
        }

        self.aligned_read_length_sum += cs.aligned_read_length;
        self.aligned_reads += 1;

        if !negative_strand {
            self.num_positive_strand += 1;
        }

        if is_hq {
            self.hq_aligned_reads += 1;
        }

        // NM-tag based mismatch tracking.
        // NM = substitutions + insertion_bases + deletion_bases; subtract indel bases
        // to obtain a pure substitution count (matching Picard's reference-comparison approach).
        let nm = u64::from(get_integer_tag(record, *b"NM").unwrap_or(0));
        let indel_bases = cs.insertion_bases + cs.deletion_bases;
        let adjusted_nm = nm.saturating_sub(indel_bases);
        self.total_nm += adjusted_nm;
        if is_hq {
            self.hq_nm += adjusted_nm;
            self.hq_nm_histogram.count(adjusted_nm);
        }

        // Paired reads: pair-in-pair and chimera logic.
        // NOTE: Picard counts PF_READS_IMPROPER_PAIRS for all mapped+paired+!proper reads,
        // including those whose mate is unmapped (no is_mate_unmapped guard).  We require
        // the mate to be mapped before counting aligned_reads_in_pairs or improperly_paired,
        // which is more conservative and avoids inflating the improper-pair count with
        // reads whose mate simply failed to align.  This produces a slightly lower value.
        if flags.is_segmented() && !flags.is_mate_unmapped() {
            self.aligned_reads_in_pairs += 1;

            if !flags.is_properly_segmented() {
                self.reads_improperly_paired += 1;
            }

            // Chimera denominator: both mapped, and either no MQ tag or both HQ.
            let mate_mq = get_integer_tag(record, *b"MQ");
            let include_chimera_check = mate_mq.is_none_or(|mq| {
                mq >= u32::from(min_mapq) && u32::from(mapq) >= u32::from(min_mapq)
            });

            if include_chimera_check {
                self.chimeras_denominator += 1;
                if is_chimeric_pair(record, max_insert_size) {
                    self.chimeras += 1;
                }
            }
        } else if is_hq {
            // Fragment or pair with unmapped mate: chimeric if SA tag present.
            self.chimeras_denominator += 1;
            if record.data().get(b"SA").is_some() {
                self.chimeras += 1;
            }
        }
    }

    /// Walk the CIGAR once and accumulate per-base quality stats into `self`,
    /// while also computing and returning summary stats that the caller needs
    /// for soft-clip, alignment-length, and indel tracking.
    fn process_cigar(
        &mut self,
        record: &RecordBuf,
        is_hq: bool,
        negative_strand: bool,
    ) -> CigarStats {
        let quals = record.quality_scores();
        let qual_bytes: &[u8] = quals.as_ref();
        let has_quals = !qual_bytes.is_empty();
        let mut read_pos: usize = 0;

        let mut soft_clip_bases = 0u64;
        let mut aligned_read_length = 0u64;
        let mut insertion_bases = 0u64;
        let mut deletion_bases = 0u64;

        // Track leading and trailing soft clips for 3' computation.
        let mut seen_non_clip = false;
        let mut leading_sc = 0u64;
        let mut trailing_sc = 0u64;

        for op in CigarTrait::iter(record.cigar()).filter_map(Result::ok) {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    seen_non_clip = true;
                    aligned_read_length += len as u64;
                    self.pf_aligned_bases += len as u64;
                    if is_hq {
                        self.hq_aligned_bases += len as u64;
                        if has_quals {
                            // Malformed BAMs can declare a CIGAR longer than the
                            // quality array. Clamp the window to `qual_bytes.len()`
                            // so out-of-range positions are silently omitted —
                            // identical count to the pre-SIMD scalar form, which
                            // used `.get().unwrap_or(0)` (Q=0 never passes a
                            // non-zero BQ threshold, so the two paths agree).
                            let end = (read_pos + len).min(qual_bytes.len());
                            let window = &qual_bytes[read_pos..end];
                            self.hq_aligned_q20_bases +=
                                simd::count_bases_ge_q(window, BASE_QUALITY_THRESHOLD);
                        }
                    }
                    read_pos += len;
                }
                Kind::Insertion => {
                    seen_non_clip = true;
                    aligned_read_length += len as u64;
                    insertion_bases += len as u64;
                    self.indels += 1;
                    read_pos += len;
                }
                Kind::Deletion => {
                    seen_non_clip = true;
                    deletion_bases += len as u64;
                    self.indels += 1;
                }
                Kind::SoftClip => {
                    soft_clip_bases += len as u64;
                    if seen_non_clip {
                        trailing_sc += len as u64;
                    } else {
                        leading_sc += len as u64;
                    }
                    read_pos += len;
                }
                Kind::HardClip | Kind::Skip | Kind::Pad => {}
            }
        }

        let three_prime_soft_clip = if negative_strand { leading_sc } else { trailing_sc };

        CigarStats {
            soft_clip_bases,
            three_prime_soft_clip,
            aligned_read_length,
            insertion_bases,
            deletion_bases,
        }
    }

    /// Merge two accumulators into a new one with `category = "pair"`.
    ///
    /// All scalar counters are summed; histograms and maps are merged by summing
    /// counts per key.  Used to synthesize the `pair` metric row from `first` and
    /// `second` without a third accumulator.
    fn merge(&self, other: &Self) -> Self {
        let mut read_length_histogram = self.read_length_histogram.clone();
        read_length_histogram += &other.read_length_histogram;
        let mut hq_nm_histogram = self.hq_nm_histogram.clone();
        hq_nm_histogram += &other.hq_nm_histogram;
        let mut bad_cycle_nocalls = self.bad_cycle_nocalls.clone();
        bad_cycle_nocalls += &other.bad_cycle_nocalls;
        Self {
            category: "pair",
            total_reads: self.total_reads + other.total_reads,
            pf_reads: self.pf_reads + other.pf_reads,
            aligned_reads: self.aligned_reads + other.aligned_reads,
            pf_aligned_bases: self.pf_aligned_bases + other.pf_aligned_bases,
            hq_aligned_reads: self.hq_aligned_reads + other.hq_aligned_reads,
            hq_aligned_bases: self.hq_aligned_bases + other.hq_aligned_bases,
            hq_aligned_q20_bases: self.hq_aligned_q20_bases + other.hq_aligned_q20_bases,
            aligned_reads_in_pairs: self.aligned_reads_in_pairs + other.aligned_reads_in_pairs,
            reads_improperly_paired: self.reads_improperly_paired + other.reads_improperly_paired,
            num_positive_strand: self.num_positive_strand + other.num_positive_strand,
            chimeras: self.chimeras + other.chimeras,
            chimeras_denominator: self.chimeras_denominator + other.chimeras_denominator,
            num_soft_clipped: self.num_soft_clipped + other.num_soft_clipped,
            num_hard_clipped: self.num_hard_clipped + other.num_hard_clipped,
            num_3prime_soft_clipped_bases: self.num_3prime_soft_clipped_bases
                + other.num_3prime_soft_clipped_bases,
            num_reads_with_3prime_softclips: self.num_reads_with_3prime_softclips
                + other.num_reads_with_3prime_softclips,
            indels: self.indels + other.indels,
            total_nm: self.total_nm + other.total_nm,
            hq_nm: self.hq_nm + other.hq_nm,
            hq_nm_histogram,
            read_length_histogram,
            read_length_sum: self.read_length_sum + other.read_length_sum,
            read_length_sum_sq: self.read_length_sum_sq + other.read_length_sum_sq,
            aligned_read_length_sum: self.aligned_read_length_sum + other.aligned_read_length_sum,
            bad_cycle_nocalls,
        }
    }

    /// Compute the final metric row from accumulated counts.
    fn compute_metric(&self) -> AlignmentSummaryMetric {
        let pf_reads = self.pf_reads;
        let pf_bases = self.read_length_sum; // denominator for clip fractions

        // Read length stats.
        let mean_rl = safe_div_f(self.read_length_sum as f64, pf_reads as f64);
        let var_rl = if pf_reads > 0 {
            (self.read_length_sum_sq as f64 / pf_reads as f64) - mean_rl * mean_rl
        } else {
            0.0
        };
        let sd_rl = var_rl.max(0.0).sqrt();
        let (median_rl, mad_rl) = self.read_length_histogram.median_and_mad();
        let min_rl = self.read_length_histogram.min().unwrap_or(0);
        let max_read_len = self.read_length_histogram.max().unwrap_or(0);

        // Aligned read length.
        // NOTE: Picard computes mean_aligned_read_length over ALL PF reads (including unmapped,
        // which contribute 0), so their denominator is pf_reads.  We use aligned_reads instead,
        // which gives the mean length of reads that actually aligned — a more intuitive metric.
        // This produces a slightly higher value than Picard (~1.5 bp for typical WGS data).
        let mean_aligned_len =
            safe_div_f(self.aligned_read_length_sum as f64, self.aligned_reads as f64);

        // Mismatch metrics.
        let mismatch_rate = safe_div_f(self.total_nm as f64, self.pf_aligned_bases as f64);
        let hq_mismatch_rate = safe_div_f(self.hq_nm as f64, self.hq_aligned_bases as f64);
        let hq_median_mismatches = self.hq_nm_histogram.median();

        let indel_rate = safe_div_f(self.indels as f64, self.pf_aligned_bases as f64);

        // Fraction metrics.
        let frac_aligned = safe_div(self.aligned_reads, pf_reads);
        let frac_aligned_in_pairs = safe_div(self.aligned_reads_in_pairs, self.aligned_reads);
        let frac_improper = safe_div(self.reads_improperly_paired, self.aligned_reads);
        let strand_balance = safe_div(self.num_positive_strand, self.aligned_reads);
        let frac_chimeras = safe_div(self.chimeras, self.chimeras_denominator);

        // Clip fractions.
        let frac_softclipped_reads = safe_div_f(self.num_soft_clipped as f64, pf_bases as f64);
        let frac_hardclipped_reads = safe_div_f(self.num_hard_clipped as f64, pf_bases as f64);
        let avg_3prime = safe_div_f(
            self.num_3prime_soft_clipped_bases as f64,
            self.num_reads_with_3prime_softclips as f64,
        );

        let bad_cycles = count_bad_cycles(&self.bad_cycle_nocalls, self.total_reads);

        AlignmentSummaryMetric {
            sample: String::new(),
            category: self.category.to_string(),
            total_reads: self.total_reads,
            aligned_reads: self.aligned_reads,
            frac_aligned,
            hq_aligned_reads: self.hq_aligned_reads,
            hq_aligned_bases: self.hq_aligned_bases,
            hq_aligned_q20_bases: self.hq_aligned_q20_bases,
            hq_median_mismatches,
            mismatch_rate,
            hq_mismatch_rate,
            indel_rate,
            mean_read_length: mean_rl,
            sd_read_length: sd_rl,
            median_read_length: median_rl,
            mad_read_length: mad_rl,
            min_read_length: min_rl,
            max_read_length: max_read_len,
            mean_aligned_read_length: mean_aligned_len,
            aligned_reads_in_pairs: self.aligned_reads_in_pairs,
            frac_aligned_in_pairs,
            reads_improperly_paired: self.reads_improperly_paired,
            frac_reads_improperly_paired: frac_improper,
            bad_cycles,
            strand_balance,
            frac_chimeras,
            frac_softclipped_reads,
            frac_hardclipped_reads,
            mean_3prime_softclipped_bases: avg_3prime,
        }
    }
}

// ─── CIGAR stats ──────────────────────────────────────────────────────────────

/// Summary statistics computed from a single CIGAR pass over a mapped record.
struct CigarStats {
    soft_clip_bases: u64,
    three_prime_soft_clip: u64,
    aligned_read_length: u64,
    insertion_bases: u64,
    deletion_bases: u64,
}

// ─── CIGAR helpers ────────────────────────────────────────────────────────────

/// Sum the lengths of all CIGAR operators of a given kind.
/// Used only for `HardClip` on the pre-mapping-filter path.
fn sum_cigar_op(record: &RecordBuf, target: Kind) -> u64 {
    CigarTrait::iter(record.cigar())
        .filter_map(Result::ok)
        .filter(|op| op.kind() == target)
        .map(|op| op.len() as u64)
        .sum()
}

// ─── Chimera detection ────────────────────────────────────────────────────────

/// Return `true` if the record represents a chimeric alignment.
///
/// A properly paired read is chimeric when any of:
/// 1. The two ends map to different reference sequences.
/// 2. The absolute insert size exceeds `max_insert_size`.
/// 3. The pair orientation is not FR (the expected orientation).
/// 4. Either end has an SA (supplementary alignment) tag.
fn is_chimeric_pair(record: &RecordBuf, max_insert_size: u32) -> bool {
    // Different contigs.
    if record.reference_sequence_id() != record.mate_reference_sequence_id() {
        return true;
    }

    // Large insert.
    let tlen = record.template_length().unsigned_abs();
    if tlen > max_insert_size {
        return true;
    }

    // Unexpected orientation — anything that isn't FR.
    match get_pair_orientation(record) {
        Some(PairOrientation::Fr) | None => {}
        Some(_) => return true,
    }

    // SA tag (split-read / chimeric alignment within one end).
    if record.data().get(b"SA").is_some() {
        return true;
    }

    false
}

// ─── Bad-cycle tracking ───────────────────────────────────────────────────────

/// Accumulate per-cycle no-call (N base) counts.
///
/// Cycle numbering follows Picard/Illumina convention:
/// - Forward strand: cycle = `position_index` + 1 (1-based from 5′).
/// - Reverse strand: cycle = `read_length` - `position_index` (1-based from 5′ of sequencing).
fn track_bad_cycles(
    bad_cycle_nocalls: &mut Counter<u64>,
    bases: &[u8],
    read_len: u64,
    negative_strand: bool,
) {
    for (i, &base) in bases.iter().enumerate() {
        if base == b'N' {
            let cycle = if negative_strand { read_len - i as u64 } else { i as u64 + 1 };
            bad_cycle_nocalls.count(cycle);
        }
    }
}

// ─── Bad-cycle counting helper ────────────────────────────────────────────────

/// Count the number of sequencing cycles where ≥ 80% of reads had a no-call (N).
fn count_bad_cycles(nocall_map: &Counter<u64>, total_reads: u64) -> u64 {
    if total_reads == 0 {
        return 0;
    }
    let total = total_reads as f64;
    nocall_map.values().filter(|&&count| count as f64 / total >= 0.8).count() as u64
}

// ─── Serializer for 6 decimal places (needed for rates) ──────────────────────

// Note: serialize_f64_6dp is declared in metrics.rs — used via the
// `crate::metrics::serialize_f64_6dp` path in #[serde(serialize_with = ...)] above.

// ─── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        Flags, MappingQuality,
        cigar::{Op, op::Kind as CigarKind},
    };
    use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};

    // ── Helper: build a minimal RecordBuf for unit tests ─────────────────────

    #[allow(clippy::too_many_arguments)]
    fn make_record(
        flags: Flags,
        pos: Option<usize>,
        mapq: u8,
        cigar: Cigar,
        seq: &[u8],
        quals: &[u8],
        ref_id: Option<usize>,
        mate_ref_id: Option<usize>,
        tlen: i32,
        data: noodles::sam::alignment::record_buf::Data,
    ) -> RecordBuf {
        let mut b = RecordBuf::builder()
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .set_template_length(tlen)
            .set_data(data);

        if let Some(rid) = ref_id {
            b = b.set_reference_sequence_id(rid);
        }
        if let Some(mrid) = mate_ref_id {
            b = b.set_mate_reference_sequence_id(mrid);
        }
        if let Some(p) = pos {
            b = b.set_alignment_start(Position::new(p).expect("pos"));
        }
        b.build()
    }

    fn simple_cigar(len: usize) -> Cigar {
        [Op::new(CigarKind::Match, len)].into_iter().collect()
    }

    // ── sum_cigar_op ──────────────────────────────────────────────────────────

    #[test]
    fn test_sum_cigar_op_soft_clip() {
        let cigar: Cigar = [
            Op::new(CigarKind::SoftClip, 5),
            Op::new(CigarKind::Match, 90),
            Op::new(CigarKind::SoftClip, 5),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        assert_eq!(sum_cigar_op(&record, Kind::SoftClip), 10);
        assert_eq!(sum_cigar_op(&record, Kind::Match), 90);
    }

    // ── process_cigar ─────────────────────────────────────────────────────────

    fn cigar_stats(record: &RecordBuf, negative_strand: bool) -> CigarStats {
        let mut acc = CategoryAccumulator::new("TEST");
        acc.process_cigar(record, false, negative_strand)
    }

    #[test]
    fn test_cigar_stats_match_only() {
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            simple_cigar(100),
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let cs = cigar_stats(&record, false);
        assert_eq!(cs.aligned_read_length, 100);
        assert_eq!(cs.soft_clip_bases, 0);
        assert_eq!(cs.three_prime_soft_clip, 0);
        assert_eq!(cs.insertion_bases, 0);
        assert_eq!(cs.deletion_bases, 0);
    }

    #[test]
    fn test_cigar_stats_with_soft_clip() {
        let cigar: Cigar = [
            Op::new(CigarKind::SoftClip, 5),
            Op::new(CigarKind::Match, 90),
            Op::new(CigarKind::SoftClip, 5),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let cs = cigar_stats(&record, false);
        assert_eq!(cs.aligned_read_length, 90);
        assert_eq!(cs.soft_clip_bases, 10);
        assert_eq!(cs.three_prime_soft_clip, 5);
    }

    #[test]
    fn test_cigar_stats_with_insertion() {
        let cigar: Cigar = [
            Op::new(CigarKind::Match, 80),
            Op::new(CigarKind::Insertion, 5),
            Op::new(CigarKind::Match, 10),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 95],
            &[30u8; 95],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let cs = cigar_stats(&record, false);
        // M + I + M = 80 + 5 + 10 = 95.
        assert_eq!(cs.aligned_read_length, 95);
        assert_eq!(cs.insertion_bases, 5);
        assert_eq!(cs.deletion_bases, 0);
    }

    #[test]
    fn test_3prime_soft_clip_forward_trailing() {
        // 5S 90M 5S — forward strand: 3′ clip is the trailing 5S.
        let cigar: Cigar = [
            Op::new(CigarKind::SoftClip, 5),
            Op::new(CigarKind::Match, 90),
            Op::new(CigarKind::SoftClip, 5),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        assert_eq!(cigar_stats(&record, false).three_prime_soft_clip, 5);
    }

    #[test]
    fn test_3prime_soft_clip_reverse_leading() {
        // 5S 90M 5S — reverse strand: 3′ clip is the leading 5S.
        let cigar: Cigar = [
            Op::new(CigarKind::SoftClip, 5),
            Op::new(CigarKind::Match, 90),
            Op::new(CigarKind::SoftClip, 5),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::REVERSE_COMPLEMENTED,
            Some(1),
            60,
            cigar,
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        assert_eq!(cigar_stats(&record, true).three_prime_soft_clip, 5);
    }

    #[test]
    fn test_3prime_soft_clip_none() {
        // 5S 90M — only a leading soft clip (5′ on forward).
        let cigar: Cigar =
            [Op::new(CigarKind::SoftClip, 5), Op::new(CigarKind::Match, 90)].into_iter().collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 95],
            &[30u8; 95],
            Some(0),
            None,
            0,
            Data::default(),
        );
        assert_eq!(cigar_stats(&record, false).three_prime_soft_clip, 0);
    }

    // ── bad-cycle tracking ────────────────────────────────────────────────────

    #[test]
    fn test_bad_cycles_forward() {
        let mut map = Counter::<u64>::new();
        // Read "ACNGT" on forward strand: N is at index 2 → cycle 3.
        track_bad_cycles(&mut map, b"ACNGT", 5, false);
        assert_eq!(map.count_of(&3), 1);
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn test_bad_cycles_reverse() {
        let mut map = Counter::<u64>::new();
        // Read "ACNGT" (len=5) on reverse strand: N at index 2 → cycle = 5 - 2 = 3.
        track_bad_cycles(&mut map, b"ACNGT", 5, true);
        assert_eq!(map.count_of(&3), 1);
    }

    // ── process_cigar accumulator stats ──────────────────────────────────────

    #[test]
    fn test_process_cigar_indels_counted() {
        // 50M 2I 48M 1D — 2 indel events.
        let cigar: Cigar = [
            Op::new(CigarKind::Match, 50),
            Op::new(CigarKind::Insertion, 2),
            Op::new(CigarKind::Match, 48),
            Op::new(CigarKind::Deletion, 1),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 100],
            &[30u8; 100],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let mut acc = CategoryAccumulator::new("TEST");
        let cs = acc.process_cigar(&record, false, false);
        assert_eq!(acc.indels, 2);
        assert_eq!(acc.pf_aligned_bases, 98); // 50M + 48M = 98 aligned bases
        assert_eq!(cs.insertion_bases, 2);
        assert_eq!(cs.deletion_bases, 1);
    }

    #[test]
    fn test_process_cigar_q20_bases() {
        // 10M, qualities: 5×Q30 + 5×Q10.
        let quals: Vec<u8> = (0..10).map(|i| if i < 5 { 30 } else { 10 }).collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            simple_cigar(10),
            &[b'A'; 10],
            &quals,
            Some(0),
            None,
            0,
            Data::default(),
        );
        let mut acc = CategoryAccumulator::new("TEST");
        acc.process_cigar(&record, true, false); // is_hq = true
        assert_eq!(acc.hq_aligned_bases, 10);
        assert_eq!(acc.hq_aligned_q20_bases, 5);
    }

    // ── process_cigar — additional edge cases ───────────────────────────────

    #[test]
    fn test_cigar_stats_with_deletion() {
        // 50M 2D 48M → aligned_read_length = 50+48 = 98 (D not in read length),
        //              deletion_bases = 2, indels = 1
        let cigar: Cigar = [
            Op::new(CigarKind::Match, 50),
            Op::new(CigarKind::Deletion, 2),
            Op::new(CigarKind::Match, 48),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 98],
            &[30u8; 98],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let mut acc = CategoryAccumulator::new("TEST");
        let cs = acc.process_cigar(&record, false, false);
        assert_eq!(cs.aligned_read_length, 98);
        assert_eq!(cs.deletion_bases, 2);
        assert_eq!(cs.insertion_bases, 0);
        assert_eq!(acc.indels, 1);
    }

    #[test]
    fn test_cigar_stats_with_hard_clip() {
        // 5H 90M 5H — hard clips don't appear in sequence/CigarStats
        let cigar: Cigar = [
            Op::new(CigarKind::HardClip, 5),
            Op::new(CigarKind::Match, 90),
            Op::new(CigarKind::HardClip, 5),
        ]
        .into_iter()
        .collect();
        let record = make_record(
            Flags::empty(),
            Some(1),
            60,
            cigar,
            &[b'A'; 90],
            &[30u8; 90],
            Some(0),
            None,
            0,
            Data::default(),
        );
        let cs = cigar_stats(&record, false);
        assert_eq!(cs.aligned_read_length, 90);
        assert_eq!(cs.soft_clip_bases, 0);
        assert_eq!(cs.insertion_bases, 0);
        assert_eq!(cs.deletion_bases, 0);
    }

    // ── is_chimeric_pair ────────────────────────────────────────────────────

    use noodles::sam::alignment::record_buf::data::field::Value as DataValue;

    /// Build a paired record suitable for chimera testing.
    fn make_chimera_record(
        ref_id: usize,
        mate_ref_id: usize,
        tlen: i32,
        is_reverse: bool,
        mate_reverse: bool,
        read_len: usize,
        sa_tag: bool,
    ) -> RecordBuf {
        let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT;
        if is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if mate_reverse {
            flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        let cigar: Cigar = [Op::new(CigarKind::Match, read_len)].into_iter().collect();
        let mut data = Data::default();
        if sa_tag {
            data.insert((*b"SA").into(), DataValue::String("chr1,100,+,50M,60,0".into()));
        }
        RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(Position::new(100).expect("pos"))
            .set_mapping_quality(MappingQuality::new(60).expect("mapq"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; read_len]))
            .set_quality_scores(QualityScores::from(vec![30u8; read_len]))
            .set_mate_reference_sequence_id(mate_ref_id)
            .set_mate_alignment_start(Position::new(200).expect("pos"))
            .set_template_length(tlen)
            .set_data(data)
            .build()
    }

    #[test]
    fn test_chimeric_different_contigs() {
        let rec = make_chimera_record(0, 1, 200, false, true, 100, false);
        assert!(is_chimeric_pair(&rec, 100_000));
    }

    #[test]
    fn test_chimeric_large_insert() {
        let rec = make_chimera_record(0, 0, 200_000, false, true, 100, false);
        assert!(is_chimeric_pair(&rec, 100_000));
    }

    #[test]
    fn test_chimeric_non_fr_orientation() {
        // Tandem: both forward → not FR → chimeric
        let rec = make_chimera_record(0, 0, 200, false, false, 100, false);
        assert!(is_chimeric_pair(&rec, 100_000));
    }

    #[test]
    fn test_chimeric_sa_tag() {
        // Normal FR pair but with SA tag → chimeric
        let rec = make_chimera_record(0, 0, 200, false, true, 100, true);
        assert!(is_chimeric_pair(&rec, 100_000));
    }

    #[test]
    fn test_not_chimeric_normal_fr() {
        let rec = make_chimera_record(0, 0, 200, false, true, 100, false);
        assert!(!is_chimeric_pair(&rec, 100_000));
    }

    #[test]
    fn test_chimeric_exact_max_insert() {
        // tlen == max_insert_size → NOT chimeric (must be strictly greater)
        let rec = make_chimera_record(0, 0, 100_000, false, true, 100, false);
        assert!(!is_chimeric_pair(&rec, 100_000));
    }

    // ── count_bad_cycles ────────────────────────────────────────────────────

    #[test]
    fn test_count_bad_cycles_at_threshold() {
        let map: Counter<u64> = [(1u64, 80u64)].into_iter().collect();
        assert_eq!(count_bad_cycles(&map, 100), 1);
    }

    #[test]
    fn test_count_bad_cycles_below_threshold() {
        let map: Counter<u64> = [(1u64, 79u64)].into_iter().collect();
        assert_eq!(count_bad_cycles(&map, 100), 0);
    }

    #[test]
    fn test_count_bad_cycles_zero_reads() {
        let map: Counter<u64> = [(1u64, 10u64)].into_iter().collect();
        assert_eq!(count_bad_cycles(&map, 0), 0);
    }

    #[test]
    fn test_count_bad_cycles_multiple() {
        let map: Counter<u64> = [(1u64, 90u64), (2, 50), (3, 85)].into_iter().collect();
        assert_eq!(count_bad_cycles(&map, 100), 2); // cycles 1 and 3
    }

    // ── compute_metric ──────────────────────────────────────────────────────

    #[test]
    fn test_compute_metric_basic() {
        let mut acc = CategoryAccumulator::new("read1");
        acc.total_reads = 100;
        acc.pf_reads = 100;
        acc.aligned_reads = 80;
        acc.pf_aligned_bases = 8000;
        acc.hq_aligned_reads = 70;
        acc.hq_aligned_bases = 7000;
        acc.hq_aligned_q20_bases = 6000;
        acc.total_nm = 80;
        acc.hq_nm = 70;
        acc.num_positive_strand = 40;
        acc.read_length_sum = 10000;
        acc.read_length_sum_sq = 1_000_000;
        acc.read_length_histogram.count_n(100, 100);
        acc.aligned_read_length_sum = 8000;

        let m = acc.compute_metric();
        assert_eq!(m.total_reads, 100);
        assert_eq!(m.aligned_reads, 80);
        assert!((m.frac_aligned - 0.8).abs() < 1e-9);
        assert!((m.mismatch_rate - 0.01).abs() < 1e-9); // 80/8000
        assert!((m.strand_balance - 0.5).abs() < 1e-9); // 40/80
    }

    #[test]
    fn test_compute_metric_zero_reads() {
        let acc = CategoryAccumulator::new("read1");
        let m = acc.compute_metric();
        assert_eq!(m.total_reads, 0);
        assert_eq!(m.frac_aligned, 0.0);
        assert_eq!(m.mismatch_rate, 0.0);
        assert_eq!(m.strand_balance, 0.0);
        assert_eq!(m.bad_cycles, 0);
        assert_eq!(m.mean_read_length, 0.0);
    }

    #[test]
    fn test_compute_metric_bad_cycles_threshold() {
        let mut acc = CategoryAccumulator::new("read1");
        acc.total_reads = 100;
        acc.pf_reads = 100;
        // Cycle 1: 80% → bad; Cycle 2: 79% → not bad
        acc.bad_cycle_nocalls.count_n(1, 80);
        acc.bad_cycle_nocalls.count_n(2, 79);
        acc.read_length_histogram.count_n(100, 100);
        acc.read_length_sum = 10000;

        let m = acc.compute_metric();
        assert_eq!(m.bad_cycles, 1);
    }

    // ── merge ───────────────────────────────────────────────────────────────

    #[test]
    fn test_merge_basic() {
        let mut a = CategoryAccumulator::new("read1");
        a.total_reads = 50;
        a.pf_reads = 45;
        a.aligned_reads = 40;
        a.pf_aligned_bases = 4000;
        a.indels = 5;
        a.total_nm = 10;
        a.read_length_sum = 5000;
        a.read_length_sum_sq = 500_000;
        a.read_length_histogram.count_n(100, 50);
        a.hq_nm_histogram.count_n(1, 10);
        a.bad_cycle_nocalls.count_n(1, 20);

        let mut b = CategoryAccumulator::new("read2");
        b.total_reads = 60;
        b.pf_reads = 55;
        b.aligned_reads = 50;
        b.pf_aligned_bases = 5000;
        b.indels = 3;
        b.total_nm = 15;
        b.read_length_sum = 6000;
        b.read_length_sum_sq = 600_000;
        b.read_length_histogram.count_n(100, 60);
        b.hq_nm_histogram.count_n(1, 15);
        b.bad_cycle_nocalls.count_n(1, 25);

        let merged = a.merge(&b);
        assert_eq!(merged.category, "pair");
        assert_eq!(merged.total_reads, 110);
        assert_eq!(merged.pf_reads, 100);
        assert_eq!(merged.aligned_reads, 90);
        assert_eq!(merged.pf_aligned_bases, 9000);
        assert_eq!(merged.indels, 8);
        assert_eq!(merged.total_nm, 25);
        assert_eq!(merged.read_length_sum, 11000);
        assert_eq!(merged.read_length_sum_sq, 1_100_000);
        assert_eq!(merged.read_length_histogram.count_of(&100), 110);
        assert_eq!(merged.hq_nm_histogram.count_of(&1), 25);
        assert_eq!(merged.bad_cycle_nocalls.count_of(&1), 45);
    }
}
