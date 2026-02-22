use std::fmt;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use bitvec::prelude::*;
use bstr::{BStr, BString};
use clap::Args;
use noodles::core::Region;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as _;
use noodles::sam::alignment::record::cigar::op::Kind;
use riker_derive::MetricDocs;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use strum::EnumCount as _;

use crate::collector::Collector;
use crate::commands::command::Command;
use crate::commands::common::{InputOptions, OutputOptions};
use crate::fasta::Fasta;
use crate::metrics::{serialize_f64_2dp, serialize_f64_6dp, write_tsv};
use crate::progress::ProgressLogger;
use crate::sam::indexed_reader::IndexedAlignmentReader;
use crate::sam::pair_orientation::{PairOrientation, get_pair_orientation};
use crate::sequence_dict::SequenceDictionary;
use crate::vcf::IndexedVcf;

// ─── Constants & type aliases ────────────────────────────────────────────────────

/// Fast HashMap using fxhash (rustc-hash) for non-cryptographic hashing.
/// Benchmarked as the fastest hasher for our small covariate key types.
type HashMap<K, V> = std::collections::HashMap<K, V, rustc_hash::FxBuildHasher>;

/// File suffix for mismatch error metrics.
pub const MISMATCH_SUFFIX: &str = ".error-mismatch.txt";

/// File suffix for overlapping-read mismatch error metrics.
pub const OVERLAP_SUFFIX: &str = ".error-overlap.txt";

/// File suffix for indel error metrics.
pub const INDEL_SUFFIX: &str = ".error-indel.txt";

/// Maximum number of stratifiers per group that are stored inline (on the stack)
/// in composite covariate keys. Groups larger than this will spill to the heap.
const INLINE_STRATIFIERS: usize = 5;

// ─── ErrorOptions ────────────────────────────────────────────────────────────────

/// Tool-specific options for the error metrics collector.
#[riker_derive::multi_options("error", "Error Metrics Options")]
#[derive(Args, Debug, Clone)]
#[command()]
pub struct ErrorOptions {
    /// Reference FASTA file (must be indexed with .fai). Required.
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: PathBuf,

    /// VCF or BCF file of known variant sites to exclude. VCF files must be
    /// bgzip-compressed and indexed (.tbi or .csi). BCF files must be indexed (.csi).
    #[arg(long, value_name = "VCF/BCF")]
    pub vcf: Option<PathBuf>,

    /// Interval list or BED file to restrict analysis to specific regions.
    /// If not specified, analyzes all contigs in the sequence dictionary.
    #[arg(long, value_name = "INTERVALS")]
    pub intervals: Option<PathBuf>,

    /// Minimum mapping quality for a read to be included.
    #[arg(long, default_value_t = ErrorOptions::DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Minimum base quality for a base to be included. For aligned bases,
    /// each base is individually checked. For insertions, the first inserted
    /// base's quality determines whether the entire insertion event is counted.
    /// Deletions are not quality-filtered (they have no associated read bases).
    #[arg(long, default_value_t = ErrorOptions::DEFAULT_MIN_BQ)]
    pub min_bq: u8,

    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,

    /// Maximum insert size for the isize stratifier. Reads with absolute insert
    /// size above this value are excluded from insert size stratification (but
    /// still counted by all other stratifiers).
    #[arg(long, default_value_t = ErrorOptions::DEFAULT_MAX_ISIZE)]
    pub max_isize: u32,

    /// TEMPORARY (validation-only): Enable Picard-compatible counting behavior.
    /// When set, (1) bases at reference-N positions are included rather than
    /// skipped, and (2) insertion bases are counted in the mismatch total_bases
    /// denominator. These match Picard's behavior to enable exact-match
    /// validation. Remove once validation is complete.
    #[arg(long, default_value_t = false, hide = true)]
    pub picard_compat: bool,

    /// Stratification groups. Each value is a comma-separated list of stratifiers.
    /// The "all" group is always included regardless. When not specified, a default
    /// set of stratifications matching Picard's defaults is used.
    ///
    /// Available stratifiers: all, bq, mapq, cycle, read_num, strand,
    /// pair_orientation, isize, gc, read_base, ref_base, hp_len,
    /// pre_dinuc, post_dinuc, context_3bp, nm, indel_len.
    #[arg(long, num_args(1..), long_help = ErrorOptions::stratify_by_help())]
    pub stratify_by: Vec<String>,
}

impl ErrorOptions {
    const DEFAULT_MIN_MAPQ: u8 = 20;
    const DEFAULT_MIN_BQ: u8 = 20;
    const DEFAULT_MAX_ISIZE: u32 = 1000;

    /// Build the long help text for --stratify-by, including available stratifiers
    /// and default groups — both generated from their source-of-truth definitions.
    fn stratify_by_help() -> String {
        use strum::IntoEnumIterator;

        let available: Vec<String> = Stratifier::iter().map(|s| s.to_string()).collect();

        let mut help = String::from(
            "Stratification groups. Each value is a comma-separated list of stratifiers. \
             The \"all\" group is always included regardless.\n\n\
             Available stratifiers: ",
        );
        help.push_str(&available.join(", "));
        help.push_str(".\n\n[default: ");
        help.push_str(&Self::DEFAULT_STRATIFY_BY.join("  "));
        help.push(']');
        help
    }

    /// Default stratification groups, matching Picard's CollectSamErrorMetrics defaults
    /// (adapted for our stratifier names and excluding stratifiers we don't support).
    const DEFAULT_STRATIFY_BY: &[&str] = &[
        "bq",
        "isize",
        "gc",
        "strand",
        "pair_orientation",
        "hp_len",
        "cycle",
        "read_num",
        "read_num,cycle",
        "read_num,hp_len",
        "read_num,gc",
        "read_num,pre_dinuc",
        "mapq",
        "nm",
        "context_3bp",
    ];
}

impl Default for ErrorOptions {
    fn default() -> Self {
        Self {
            reference: PathBuf::new(),
            vcf: None,
            intervals: None,
            min_mapq: Self::DEFAULT_MIN_MAPQ,
            min_bq: Self::DEFAULT_MIN_BQ,
            include_duplicates: false,
            max_isize: Self::DEFAULT_MAX_ISIZE,
            picard_compat: false,
            stratify_by: Vec::new(),
        }
    }
}

// ─── Error (CLI command) ─────────────────────────────────────────────────────────

/// Collect base-level error metrics from a BAM file.
///
/// Computes mismatch, overlapping-read, and indel error rates by comparing
/// aligned reads to a reference genome. Supports composable stratification
/// by base quality, cycle, read number, pair orientation, homopolymer context,
/// and more. Optionally excludes known variant sites from a VCF.
///
/// Outputs three files: <prefix>.error-mismatch.txt, <prefix>.error-overlap.txt,
/// and <prefix>.error-indel.txt.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker error -i input.bam -o out -r ref.fa
  riker error -i input.bam -o out -r ref.fa -V known.vcf.gz --stratify-by read_num,cycle bq
  riker error -i input.bam -o out -r ref.fa --intervals targets.bed --min-bq 30"
)]
pub struct Error {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub options: ErrorOptions,
}

impl Error {
    /// Build intervals from a file or auto-generate one per contig from the sequence dictionary.
    #[expect(clippy::cast_possible_truncation, reason = "contig lengths fit in u32")]
    #[expect(clippy::type_complexity, reason = "internal helper, tuple is clear in context")]
    fn build_intervals(
        options: &ErrorOptions,
        dict: &SequenceDictionary,
    ) -> Result<(Option<crate::intervals::Intervals>, Vec<(String, u32, u32)>)> {
        let parsed = options
            .intervals
            .as_ref()
            .map(|p| crate::intervals::Intervals::from_path(p, dict.clone()))
            .transpose()?;

        let intervals = if let Some(ref ivs) = parsed {
            let mut result = Vec::new();
            for ref_id in 0..dict.len() {
                if let Some(meta) = dict.get_by_index(ref_id) {
                    for interval in ivs.get_contig(ref_id) {
                        result.push((meta.name().to_string(), interval.start, interval.end));
                    }
                }
            }
            result
        } else {
            dict.iter().map(|meta| (meta.name().to_string(), 0u32, meta.length() as u32)).collect()
        };

        Ok((parsed, intervals))
    }
}

impl Command for Error {
    #[expect(clippy::cast_possible_truncation, reason = "contig lengths fit in u32")]
    fn execute(&self) -> Result<()> {
        // Upfront validation
        let ref_path = &self.options.reference;
        let fai_path = {
            let mut p = ref_path.as_os_str().to_owned();
            p.push(".fai");
            PathBuf::from(p)
        };
        if !fai_path.exists() {
            bail!(
                "FASTA index not found: expected {}\n\
                 Run `samtools faidx {}` to create it.",
                fai_path.display(),
                ref_path.display(),
            );
        }

        // Open reference
        let fasta = Fasta::from_path(ref_path)?;

        // Open indexed alignment file
        let mut alignment_reader = IndexedAlignmentReader::new(&self.input.input, Some(ref_path))?;
        let header = alignment_reader.header().clone();

        // Build sequence dictionary
        let dict = SequenceDictionary::from(&header);

        // Create collector (constructor handles group parsing, partitioning, accumulators)
        let mut collector = ErrorCollector::new(&self.output.output, fasta, &self.options)?;

        let (parsed_intervals, intervals) = Self::build_intervals(&self.options, &dict)?;

        // Pre-load variant masks from VCF if provided
        if let Some(vcf_path) = &self.options.vcf {
            let mut vcf = IndexedVcf::from_path(vcf_path)?;
            collector.variant_masks =
                crate::vcf::load_variant_masks(&mut vcf, &dict, parsed_intervals.as_ref())?;
            // vcf is dropped here — never stored on collector
        }

        log::info!("Processing {} intervals across {} contigs", intervals.len(), dict.len());

        let mut progress = ProgressLogger::new("error", "reads", 1_000_000);
        let mut total_records = 0u64;
        let mut contig_cache = ContigCache { name: String::new(), bases: Vec::new() };

        for (contig_name, start, end) in &intervals {
            // Load region context (reference bases + pre-loaded variant masks)
            let region = RegionContext::new(
                contig_name,
                *start,
                *end,
                &mut collector.reference,
                collector.variant_masks.get(contig_name.as_str()),
                &mut contig_cache,
            )?;

            // Query BAM for this region
            let region_str = if *start == 0
                && *end >= dict.get_by_name(contig_name).map_or(0, |m| m.length() as u32)
            {
                contig_name.clone()
            } else {
                format!("{contig_name}:{}-{end}", start + 1)
            };
            let query_region: Region = region_str
                .parse()
                .with_context(|| format!("Failed to parse region: {region_str}"))?;

            let records = alignment_reader.query(&query_region)?;

            for record in records {
                if !collector.passes_filters(&record) {
                    continue;
                }

                total_records += 1;
                progress.record_with(&record, &header);

                match collector.mate_buffer.accept_read(record) {
                    MateResult::Buffered => {}
                    MateResult::Single(rec) => {
                        collector.process_record(&rec, &region, None);
                    }
                    MateResult::Pair(rec, mate) => {
                        let mate_bases = build_aligned_bases(&mate, collector.min_bq);
                        let record_bases = build_aligned_bases(&rec, collector.min_bq);
                        collector.process_record(&mate, &region, Some(&record_bases));
                        collector.process_record(&rec, &region, Some(&mate_bases));
                    }
                }
            }

            // Flush orphans — process for mismatch/indel only
            let ref_id = dict
                .get_by_name(contig_name)
                .map_or(0, crate::sequence_dict::SequenceMetadata::index);
            for orphan in collector.mate_buffer.flush_behind(ref_id, *end) {
                collector.process_record(&orphan, &region, None);
            }
        }

        // Discard remaining buffered mates whose mates never arrived.
        // In the standalone path these are from the last interval and have no
        // further region context to process against.
        drop(collector.mate_buffer.drain_all());

        progress.finish();
        log::info!("Processed {total_records} records passing filters.");

        // Write output
        collector.write_output()?;

        Ok(())
    }
}

// ─── ErrorCollector ──────────────────────────────────────────────────────────────

/// Collects base-level error metrics from aligned BAM records.
pub struct ErrorCollector {
    // Output paths
    mismatch_path: PathBuf,
    overlap_path: PathBuf,
    indel_path: PathBuf,

    // Config
    reference: Fasta,
    vcf_path: Option<PathBuf>,
    variant_masks: std::collections::HashMap<String, BitVec>,
    min_mapq: u8,
    min_bq: u8,
    include_duplicates: bool,
    max_isize: u32,
    picard_compat: bool, // TEMPORARY: Picard-compat mode for validation — remove later

    // Stratifier groups (always includes "all" as first group)
    groups: Vec<StratifierGroup>,

    // Stratifiers split by level for efficient caching
    read_level_stratifiers: Vec<Stratifier>,
    base_level_stratifiers: Vec<Stratifier>,

    // Mate buffer for overlap detection
    mate_buffer: MateBuffer,

    // String interner for compact covariate keys
    interner: StringInterner,

    // Accumulators — one per stratifier group
    accumulators: Vec<GroupAccumulators>,

    // For Collector trait (sequential processing)
    dict: Option<SequenceDictionary>,
    intervals: Option<crate::intervals::Intervals>,
    interval_path: Option<PathBuf>,
    current_ref_id: Option<usize>,
    current_region: Option<RegionContext>,
    contig_cache: ContigCache,
}

impl ErrorCollector {
    /// Create a new error collector.
    ///
    /// Parses stratifier groups from `options.stratify_by`, partitions them into
    /// read-level and base-level stratifiers, and builds accumulators. The
    /// `interval_path` is stored for deferred loading (requires a sequence
    /// dictionary from the BAM header).
    ///
    /// # Errors
    /// Returns an error if stratifier group names are invalid.
    pub fn new(prefix: &Path, reference: Fasta, options: &ErrorOptions) -> Result<Self> {
        let mismatch_path = super::command::output_path(prefix, MISMATCH_SUFFIX);
        let overlap_path = super::command::output_path(prefix, OVERLAP_SUFFIX);
        let indel_path = super::command::output_path(prefix, INDEL_SUFFIX);

        // Parse stratifier groups and build accumulators eagerly
        let groups = Self::parse_groups(&options.stratify_by)?;
        log::info!(
            "Stratifier groups: {}",
            groups.iter().map(|g| g.name.as_str()).collect::<Vec<_>>().join(", ")
        );
        let (read_level_stratifiers, base_level_stratifiers) = Stratifier::partition(&groups);
        let accumulators = groups.iter().map(|g| GroupAccumulators::new(g.clone())).collect();

        Ok(Self {
            mismatch_path,
            overlap_path,
            indel_path,
            reference,
            vcf_path: options.vcf.clone(),
            variant_masks: std::collections::HashMap::new(),
            min_mapq: options.min_mapq,
            min_bq: options.min_bq,
            include_duplicates: options.include_duplicates,
            max_isize: options.max_isize,
            picard_compat: options.picard_compat,
            groups,
            read_level_stratifiers,
            base_level_stratifiers,
            mate_buffer: MateBuffer::new(),
            interner: StringInterner::new(),
            accumulators,
            dict: None,
            intervals: None,
            interval_path: options.intervals.clone(),
            current_ref_id: None,
            current_region: None,
            contig_cache: ContigCache { name: String::new(), bases: Vec::new() },
        })
    }

    /// Parse stratifier groups from the options and ensure "all" is always present.
    ///
    /// When `stratify_by` is empty, uses a default set of stratifications matching
    /// Picard's `CollectSamErrorMetrics` defaults (adapted for our stratifier names).
    fn parse_groups(stratify_by: &[String]) -> Result<Vec<StratifierGroup>> {
        let mut groups = Vec::new();

        // Always include "all" as the first group
        groups
            .push(StratifierGroup { stratifiers: vec![Stratifier::All], name: "all".to_string() });

        let specs: Vec<String> = if stratify_by.is_empty() {
            ErrorOptions::DEFAULT_STRATIFY_BY.iter().map(|s| (*s).to_string()).collect()
        } else {
            stratify_by.to_vec()
        };

        for spec in &specs {
            let group = StratifierGroup::parse(spec)?;
            // Don't duplicate "all" if the user specified it
            if group.stratifiers != vec![Stratifier::All] {
                groups.push(group);
            }
        }

        Ok(groups)
    }

    /// Check if a record passes the basic filters.
    fn passes_filters(&self, record: &RecordBuf) -> bool {
        let flags = record.flags();

        // Skip unmapped, secondary, supplementary, QC-fail
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
        {
            return false;
        }

        // Skip duplicates unless include_duplicates is set
        if !self.include_duplicates && flags.is_duplicate() {
            return false;
        }

        // Skip low MAPQ
        let mapq = record.mapping_quality().map_or(0, |m| m.get());
        if mapq < self.min_mapq {
            return false;
        }

        true
    }

    /// Convert accumulators to sorted metric rows and write TSV output.
    #[expect(
        clippy::cast_precision_loss,
        reason = "f64 precision is sufficient for error rate calculations"
    )]
    fn write_output(&self) -> Result<()> {
        let mut mismatch_rows = Vec::new();
        let mut overlap_rows = Vec::new();
        let mut indel_rows = Vec::new();

        for accum in &self.accumulators {
            let group_name = &accum.group.name;

            // Sort entries by formatted key for deterministic output
            let mut entries: Vec<_> = accum.accums.iter().collect();
            entries.sort_by(|a, b| a.0.format(&self.interner).cmp(&b.0.format(&self.interner)));

            for (key, combined) in &entries {
                let formatted_key = key.format(&self.interner);

                // Mismatch metrics
                let mm = &combined.mismatch;
                if mm.total_bases > 0 {
                    let frac = mm.error_bases as f64 / mm.total_bases as f64;
                    mismatch_rows.push(MismatchMetric {
                        stratifier: group_name.clone(),
                        covariate: formatted_key.clone(),
                        total_bases: mm.total_bases,
                        error_bases: mm.error_bases,
                        frac_error: frac,
                        q_score: q_score(mm.error_bases, mm.total_bases),
                    });
                }

                // Overlap metrics
                let ov = &combined.overlap;
                if ov.overlapping_read_bases > 0 {
                    overlap_rows.push(OverlappingMismatchMetric {
                        stratifier: group_name.clone(),
                        covariate: formatted_key.clone(),
                        overlapping_read_bases: ov.overlapping_read_bases,
                        bases_mismatching_ref_and_mate: ov.bases_mismatching_ref_and_mate,
                        bases_matching_mate_but_not_ref: ov.bases_matching_mate_but_not_ref,
                        bases_in_three_way_disagreement: ov.bases_in_three_way_disagreement,
                        q_mismatching_ref_and_mate: q_score(
                            ov.bases_mismatching_ref_and_mate,
                            ov.overlapping_read_bases,
                        ),
                        q_matching_mate_but_not_ref: q_score(
                            ov.bases_matching_mate_but_not_ref,
                            ov.overlapping_read_bases,
                        ),
                        q_three_way_disagreement: q_score(
                            ov.bases_in_three_way_disagreement,
                            ov.overlapping_read_bases,
                        ),
                    });
                }

                // Indel metrics
                let ind = &combined.indel;
                if ind.total_bases > 0 || ind.num_insertions > 0 || ind.num_deletions > 0 {
                    let total_events = ind.num_insertions + ind.num_deletions;
                    let total_for_rate = ind.total_bases.max(1);
                    indel_rows.push(IndelMetric {
                        stratifier: group_name.clone(),
                        covariate: formatted_key,
                        total_bases: ind.total_bases,
                        num_insertions: ind.num_insertions,
                        num_inserted_bases: ind.num_inserted_bases,
                        num_deletions: ind.num_deletions,
                        num_deleted_bases: ind.num_deleted_bases,
                        frac_indel_error: total_events as f64 / total_for_rate as f64,
                        q_score: q_score(total_events, ind.total_bases),
                    });
                }
            }
        }

        write_tsv(&self.mismatch_path, &mismatch_rows)?;
        write_tsv(&self.overlap_path, &overlap_rows)?;
        write_tsv(&self.indel_path, &indel_rows)?;

        log::info!("Wrote mismatch metrics to: {}", self.mismatch_path.display());
        log::info!("Wrote overlap metrics to:  {}", self.overlap_path.display());
        log::info!("Wrote indel metrics to:    {}", self.indel_path.display());

        Ok(())
    }

    /// Walk a record's CIGAR and process all aligned bases and indel events.
    ///
    /// For each aligned (M/=/X) base, calls `process_aligned_base` to accumulate
    /// mismatch, indel total, and overlap metrics. For each insertion or deletion,
    /// calls `process_indel`.
    ///
    /// `mate_bases` is an optional sorted slice of `(ref_pos, read_base)` pairs
    /// from the overlapping mate read, used for overlap error classification.
    /// When `None`, overlap metrics are not accumulated for this record.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "genomic positions and read lengths fit in u32"
    )]
    #[expect(clippy::cast_sign_loss, reason = "cycle is always positive after arithmetic")]
    #[expect(clippy::cast_possible_wrap, reason = "cycle values fit in i32")]
    #[expect(
        clippy::too_many_lines,
        reason = "CIGAR walk with inline processing is clearest as one function"
    )]
    fn process_record(
        &mut self,
        record: &RecordBuf,
        region: &RegionContext,
        mate_bases: Option<&[(u32, u8)]>,
    ) {
        let is_reverse = record.flags().is_reverse_complemented();
        let read_len = record.sequence().len();
        let Some(alignment_start) = record.alignment_start() else { return };
        let alignment_start = (alignment_start.get() - 1) as u32; // convert to 0-based

        // Build per-read caches once — avoids redundant computation per base per group
        let read_cache =
            ReadLevelCache::new(record, &self.interner, self.max_isize, self.picard_compat);
        let mut base_cache = BaseCovariateCache::from_read_level(
            &self.read_level_stratifiers,
            &read_cache,
            &self.interner,
        );

        let mut ref_pos = alignment_start;
        let mut read_pos: usize = 0;
        let mut cycle: u32 = if is_reverse { read_len as u32 } else { 1 };
        let cycle_dir: i32 = if is_reverse { -1 } else { 1 };
        let mut mate_scan_idx: usize = 0;

        // Track the last valid anchor position for indels
        let mut last_anchor_read_offset: Option<usize> = None;
        let mut last_anchor_ref_pos: Option<u32> = None;
        let mut last_anchor_cycle: Option<u32> = None;

        for op_result in record.cigar().iter() {
            let Ok(op) = op_result else { return };

            let kind = op.kind();
            let len = op.len();

            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        self.process_aligned_base(
                            record,
                            read_pos,
                            ref_pos,
                            region,
                            cycle,
                            mate_bases,
                            &mut mate_scan_idx,
                            &read_cache,
                            &mut base_cache,
                        );
                        last_anchor_read_offset = Some(read_pos);
                        last_anchor_ref_pos = Some(ref_pos);
                        last_anchor_cycle = Some(cycle);
                        ref_pos += 1;
                        read_pos += 1;
                        cycle = (cycle as i32 + cycle_dir) as u32;
                    }
                }
                Kind::Insertion => {
                    // Skip the entire insertion if the first inserted base
                    // fails the base quality threshold.
                    let quals = record.quality_scores().as_ref();
                    let first_bq = quals.get(read_pos).copied().unwrap_or(0);

                    if first_bq >= self.min_bq {
                        // TEMPORARY (picard-compat): Picard counts BQ-passing insertion
                        // bases in the mismatch TOTAL_BASES. Remove later.
                        let bq_passing = if self.picard_compat { len as u32 } else { 0 };

                        // Use the last aligned base as anchor for indel stratification
                        if let (Some(anchor_ro), Some(anchor_rp), Some(anchor_cy)) =
                            (last_anchor_read_offset, last_anchor_ref_pos, last_anchor_cycle)
                        {
                            self.process_indel(
                                record,
                                anchor_ro,
                                anchor_rp,
                                region,
                                anchor_cy,
                                true,
                                len as u32,
                                bq_passing,
                                &read_cache,
                                &mut base_cache,
                            );
                        } else if bq_passing > 0 && self.picard_compat {
                            // TEMPORARY (picard-compat): no-anchor insertion (at read
                            // start or after soft clips). Add to mismatch denominator
                            // for groups that have a key from the read-level cache.
                            for (group_idx, group) in self.groups.iter().enumerate() {
                                if let Some(key) = group.key_from_cache(&base_cache) {
                                    self.accumulators[group_idx]
                                        .accums
                                        .entry(key)
                                        .or_default()
                                        .mismatch
                                        .total_bases += u64::from(bq_passing);
                                }
                            }
                        }
                    }
                    read_pos += len;
                    cycle = (cycle as i32 + cycle_dir * len as i32) as u32;
                }
                Kind::Deletion => {
                    if let (Some(anchor_ro), Some(anchor_rp), Some(anchor_cy)) =
                        (last_anchor_read_offset, last_anchor_ref_pos, last_anchor_cycle)
                    {
                        self.process_indel(
                            record,
                            anchor_ro,
                            anchor_rp,
                            region,
                            anchor_cy,
                            false,
                            len as u32,
                            0, // no BQ-passing count for deletions
                            &read_cache,
                            &mut base_cache,
                        );
                    }
                    ref_pos += len as u32;
                }
                Kind::SoftClip => {
                    read_pos += len;
                    cycle = (cycle as i32 + cycle_dir * len as i32) as u32;
                }
                Kind::HardClip | Kind::Pad => {}
                Kind::Skip => {
                    ref_pos += len as u32;
                }
            }
        }
    }

    /// Process a single aligned base for mismatch, indel total, and overlap metrics.
    ///
    /// Skips the base if it falls outside the region, is excluded (known variant),
    /// fails the base quality threshold, or is an N base. For overlap, scans forward
    /// in the sorted `mate_bases` slice using `mate_scan_idx` as a cursor — this is
    /// O(1) amortized since both reads traverse reference positions in order.
    #[expect(clippy::too_many_arguments)]
    fn process_aligned_base(
        &mut self,
        record: &RecordBuf,
        read_offset: usize,
        ref_pos: u32,
        region: &RegionContext,
        cycle: u32,
        mate_bases: Option<&[(u32, u8)]>,
        mate_scan_idx: &mut usize,
        read_cache: &ReadLevelCache,
        base_cache: &mut BaseCovariateCache,
    ) {
        // Skip if outside region or excluded
        if !region.contains(ref_pos) || region.is_excluded(ref_pos) {
            return;
        }

        // Skip if base quality too low
        let quals = record.quality_scores().as_ref();
        if read_offset < quals.len() && quals[read_offset] < self.min_bq {
            return;
        }

        let seq = record.sequence().as_ref();
        let read_base = match seq.get(read_offset) {
            Some(&b) => b.to_ascii_uppercase(),
            None => return,
        };
        let ref_base = region.ref_base(ref_pos);

        // Skip N bases — always skip read N; skip ref N unless picard-compat mode is on.
        // TEMPORARY (picard-compat): Picard counts bases at ref-N positions — remove later.
        if !is_valid_base(read_base) || (!self.picard_compat && !is_valid_base(ref_base)) {
            return;
        }

        let is_error = read_base != ref_base;

        // Find the mate's base at this ref_pos by scanning forward in the sorted slice
        let overlap_mate_base = if let Some(mate_slice) = mate_bases {
            while *mate_scan_idx < mate_slice.len() && mate_slice[*mate_scan_idx].0 < ref_pos {
                *mate_scan_idx += 1;
            }
            if *mate_scan_idx < mate_slice.len() && mate_slice[*mate_scan_idx].0 == ref_pos {
                let b = mate_slice[*mate_scan_idx].1;
                if is_valid_base(b) { Some(b) } else { None }
            } else {
                None
            }
        } else {
            None
        };

        // Update only base-level covariate slots (read-level already populated per-record)
        base_cache.populate_base_level(
            &self.base_level_stratifiers,
            record,
            read_offset,
            ref_pos,
            region,
            0, // indel_len = 0 for aligned bases
            cycle,
            read_cache,
            &self.interner,
        );

        // Then for each group, just look up precomputed values
        for (group_idx, group) in self.groups.iter().enumerate() {
            let Some(key) = group.key_from_cache(base_cache) else {
                continue;
            };

            // Single lookup into the combined accumulator
            let combined = self.accumulators[group_idx].accums.entry(key).or_default();

            // Mismatch accumulation
            combined.mismatch.total_bases += 1;
            if is_error {
                combined.mismatch.error_bases += 1;
            }

            // Indel accumulation (indel_len=0 for aligned bases counts toward total)
            combined.indel.total_bases += 1;

            // Overlap accumulation (if mate data available for this position)
            if let Some(mate_base) = overlap_mate_base {
                combined.overlap.overlapping_read_bases += 1;

                if is_error {
                    if mate_base == ref_base {
                        combined.overlap.bases_mismatching_ref_and_mate += 1;
                    } else if mate_base == read_base {
                        combined.overlap.bases_matching_mate_but_not_ref += 1;
                    } else {
                        combined.overlap.bases_in_three_way_disagreement += 1;
                    }
                }
            }
        }
    }

    /// Process an indel event (insertion or deletion).
    /// `bq_passing_insertion_bases` is the count of inserted bases with BQ >= min_bq,
    /// used only in picard-compat mode to match Picard's per-base BQ filtering of
    /// insertion bases in the mismatch denominator.
    #[expect(clippy::too_many_arguments)]
    fn process_indel(
        &mut self,
        record: &RecordBuf,
        anchor_read_offset: usize,
        anchor_ref_pos: u32,
        region: &RegionContext,
        cycle: u32,
        is_insertion: bool,
        indel_len: u32,
        bq_passing_insertion_bases: u32,
        read_cache: &ReadLevelCache,
        base_cache: &mut BaseCovariateCache,
    ) {
        if !region.contains(anchor_ref_pos) {
            return;
        }

        // For insertions, exclude if the anchor position is a known variant site.
        // For deletions, exclude if every position in the deletion span is masked
        // by a known variant. The anchor is not considered since the deletion
        // event affects the deleted reference positions, not the anchor itself.
        if is_insertion {
            if region.is_excluded(anchor_ref_pos) {
                return;
            }
        } else {
            let del_start = anchor_ref_pos + 1;
            let del_end = anchor_ref_pos + indel_len;
            let all_excluded =
                (del_start..=del_end).all(|pos| !region.contains(pos) || region.is_excluded(pos));
            if all_excluded {
                return;
            }
        }

        // Re-populate base-level slots with the anchor position and actual indel_len
        base_cache.populate_base_level(
            &self.base_level_stratifiers,
            record,
            anchor_read_offset,
            anchor_ref_pos,
            region,
            indel_len,
            cycle,
            read_cache,
            &self.interner,
        );

        for (group_idx, group) in self.groups.iter().enumerate() {
            let Some(key) = group.key_from_cache(base_cache) else {
                continue;
            };

            let combined = self.accumulators[group_idx].accums.entry(key).or_default();
            if is_insertion {
                combined.indel.num_insertions += 1;
                combined.indel.num_inserted_bases += u64::from(indel_len);
                // TEMPORARY (picard-compat): Picard counts BQ-passing insertion bases
                // in mismatch TOTAL_BASES via BaseErrorCalculator.addBase() — remove later.
                if self.picard_compat {
                    combined.mismatch.total_bases += u64::from(bq_passing_insertion_bases);
                }
            } else {
                combined.indel.num_deletions += 1;
                combined.indel.num_deleted_bases += u64::from(indel_len);
            }
        }
    }
}

impl Collector for ErrorCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> {
        let dict = SequenceDictionary::from(header);
        if let Some(interval_path) = &self.interval_path {
            self.intervals =
                Some(crate::intervals::Intervals::from_path(interval_path, dict.clone())?);
        }
        if let Some(vcf_path) = &self.vcf_path {
            let mut vcf = IndexedVcf::from_path(vcf_path)?;
            self.variant_masks =
                crate::vcf::load_variant_masks(&mut vcf, &dict, self.intervals.as_ref())?;
            // vcf is dropped here — never stored on self
        }
        self.dict = Some(dict);
        Ok(())
    }

    #[expect(clippy::cast_possible_truncation, reason = "contig lengths fit in u32")]
    fn accept(&mut self, record: &RecordBuf, _header: &Header) -> Result<()> {
        if !self.passes_filters(record) {
            return Ok(());
        }

        let Some(ref_id) = record.reference_sequence_id() else {
            return Ok(());
        };

        // Handle contig transitions
        if self.current_ref_id != Some(ref_id) {
            // Process orphaned buffered reads from previous contig without overlap data
            if let Some(region) = self.current_region.take() {
                for orphan in self.mate_buffer.drain_all() {
                    self.process_record(&orphan, &region, None);
                }
                self.current_region = Some(region);
            } else {
                drop(self.mate_buffer.drain_all());
            }

            let dict = self.dict.as_ref().unwrap();
            let Some(meta) = dict.get_by_index(ref_id) else {
                return Ok(());
            };

            // Skip contigs with no intervals (when intervals are specified)
            if let Some(intervals) = &self.intervals
                && !intervals.has_contig(ref_id)
            {
                self.current_ref_id = Some(ref_id);
                self.current_region = None;
                return Ok(());
            }

            let contig_name = meta.name();
            let contig_len = meta.length() as u32;

            self.current_region = Some(RegionContext::new_for_contig(
                contig_name,
                contig_len,
                ref_id,
                &mut self.reference,
                self.variant_masks.get(contig_name),
                self.intervals.as_ref(),
                &mut self.contig_cache,
            )?);
            self.current_ref_id = Some(ref_id);
        }

        let Some(_region) = &self.current_region else {
            return Ok(()); // no region for this contig (no intervals)
        };

        // Temporarily take the region out of self to satisfy the borrow checker:
        // process_record borrows &mut self, but we also need &region.
        // This is zero-cost (Option take/replace is just pointer swaps).
        let region = self.current_region.take().unwrap();

        // Clone required since Collector trait gives us &RecordBuf
        match self.mate_buffer.accept_read(record.clone()) {
            MateResult::Buffered => {}
            MateResult::Single(rec) => {
                self.process_record(&rec, &region, None);
            }
            MateResult::Pair(rec, mate) => {
                let mate_bases = build_aligned_bases(&mate, self.min_bq);
                let record_bases = build_aligned_bases(&rec, self.min_bq);
                self.process_record(&mate, &region, Some(&record_bases));
                self.process_record(&rec, &region, Some(&mate_bases));
            }
        }

        self.current_region = Some(region);

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        // Process remaining buffered reads without overlap data
        if let Some(region) = self.current_region.take() {
            for orphan in self.mate_buffer.drain_all() {
                self.process_record(&orphan, &region, None);
            }
            self.current_region = Some(region);
        } else {
            drop(self.mate_buffer.drain_all());
        }
        self.write_output()
    }

    fn name(&self) -> &'static str {
        "error"
    }
}

// ─── Metric structs ─────────────────────────────────────────────────────────────

/// Mismatch error metrics — one row per stratifier group + covariate combination.
///
/// Each row reports the total bases examined and the number that mismatched
/// the reference, along with the error rate and phred-scaled Q-score.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MismatchMetric {
    /// The stratification group (e.g. "read_num,cycle").
    pub stratifier: String,
    /// The covariate value(s) for this row.
    pub covariate: String,
    /// Total bases passing quality and variant-site filters.
    pub total_bases: u64,
    /// Number of bases mismatching the reference.
    pub error_bases: u64,
    /// Fraction of bases that are errors (error_bases / total_bases).
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_error: f64,
    /// Phred-scaled error rate: -10 * log10(error_rate).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub q_score: f64,
}

/// Overlapping-read mismatch metrics — one row per stratifier + covariate.
///
/// For overlapping paired-end reads, classifies each overlapping base into
/// one of three disagreement categories to distinguish sequencing errors
/// from pre-sequencing (e.g. PCR or DNA damage) errors.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct OverlappingMismatchMetric {
    /// The stratification group.
    pub stratifier: String,
    /// The covariate value(s) for this row.
    pub covariate: String,
    /// Number of bases from reads overlapping their mate at this position.
    pub overlapping_read_bases: u64,
    /// Bases where read disagrees with both reference and mate (mate agrees with ref). Likely sequencing errors.
    pub bases_mismatching_ref_and_mate: u64,
    /// Bases where read and mate agree with each other but both disagree with reference. Likely pre-sequencing errors.
    pub bases_matching_mate_but_not_ref: u64,
    /// Bases where read, mate, and reference all differ from each other.
    pub bases_in_three_way_disagreement: u64,
    /// Q-score for bases mismatching ref and mate (likely sequencing errors).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub q_mismatching_ref_and_mate: f64,
    /// Q-score for bases matching mate but not ref (likely pre-sequencing errors).
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub q_matching_mate_but_not_ref: f64,
    /// Q-score for three-way disagreements.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub q_three_way_disagreement: f64,
}

/// Indel error metrics — one row per stratifier + covariate.
///
/// Counts insertion and deletion events from CIGAR strings. The "anchor base"
/// (base immediately before the indel in sequencing order) is used for
/// per-base stratifiers like base quality and read base.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct IndelMetric {
    /// The stratification group.
    pub stratifier: String,
    /// The covariate value(s) for this row.
    pub covariate: String,
    /// Total bases considered (alignment-consuming bases).
    pub total_bases: u64,
    /// Number of insertion events.
    pub num_insertions: u64,
    /// Total bases inserted across all insertion events.
    pub num_inserted_bases: u64,
    /// Number of deletion events.
    pub num_deletions: u64,
    /// Total bases deleted across all deletion events.
    pub num_deleted_bases: u64,
    /// Fraction of indel events relative to total bases.
    #[serde(serialize_with = "serialize_f64_6dp")]
    pub frac_indel_error: f64,
    /// Phred-scaled indel error rate.
    #[serde(serialize_with = "serialize_f64_2dp")]
    pub q_score: f64,
}

// ─── Helper structs & enums ─────────────────────────────────────────────────────

/// A single stratification axis for error metrics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, strum::EnumCount, strum::EnumIter)]
pub enum Stratifier {
    /// No stratification — all bases in one bin.
    All,
    /// Base quality at the current position.
    Bq,
    /// Mapping quality of the read.
    Mapq,
    /// 1-based machine cycle (position in read, sequencing order).
    Cycle,
    /// First or second read of pair.
    ReadNum,
    /// Alignment direction (+ or -).
    Strand,
    /// Pair orientation (f1r2, f2r1, tandem).
    PairOrientation,
    /// Absolute insert size.
    Isize,
    /// GC content of the read (0-100 percent).
    Gc,
    /// The read base at this position (in sequencing orientation).
    ReadBase,
    /// The reference base at this position (in sequencing orientation).
    RefBase,
    /// Homopolymer run length preceding this base in sequencing direction.
    HpLen,
    /// Previous read base + current reference base (dinucleotide).
    PreDinuc,
    /// Next read base + current reference base (dinucleotide).
    PostDinuc,
    /// 3-base context: prev read base + ref base + next read base.
    Context3bp,
    /// NM tag value, minus 1 if the current base is itself a mismatch.
    Nm,
    /// Indel length (0 for aligned/mismatched bases).
    IndelLen,
}

impl Stratifier {
    /// Split unique stratifiers across all groups into read-level and base-level.
    fn partition(groups: &[StratifierGroup]) -> (Vec<Self>, Vec<Self>) {
        let mut seen = [false; Self::COUNT];
        let mut read_level = Vec::new();
        let mut base_level = Vec::new();
        for group in groups {
            for &strat in &group.stratifiers {
                if !seen[strat.index()] {
                    seen[strat.index()] = true;
                    if strat.is_read_level() {
                        read_level.push(strat);
                    } else {
                        base_level.push(strat);
                    }
                }
            }
        }
        (read_level, base_level)
    }

    /// Returns true if this stratifier's value is constant across all bases in a read.
    fn is_read_level(self) -> bool {
        matches!(
            self,
            Self::All
                | Self::ReadNum
                | Self::Strand
                | Self::PairOrientation
                | Self::Isize
                | Self::Gc
                | Self::Mapq
        )
    }

    /// Compute the covariate value for this stratifier at the given base position.
    ///
    /// Returns `None` if the stratifier cannot produce a value (e.g. NM tag missing,
    /// pre_dinuc at the first base of a read). Bases returning `None` are excluded
    /// from the corresponding stratifier group.
    ///
    /// Read-level stratifiers (read_num, strand, pair_orientation, isize, gc, mapq, nm)
    /// use pre-computed values from `cache` to avoid redundant per-base work.
    #[expect(clippy::too_many_arguments)]
    fn covariate(
        self,
        record: &RecordBuf,
        read_offset: usize,
        ref_pos: u32,
        region: &RegionContext,
        indel_len: u32,
        cycle: u32,
        cache: &ReadLevelCache,
        interner: &StringInterner,
    ) -> Option<CovariateValue> {
        let is_reverse = record.flags().is_reverse_complemented();
        match self {
            Self::All => Some(CovariateValue::InternedStr(interner.get("all"))),
            Self::Bq => {
                let quals = record.quality_scores().as_ref();
                if read_offset < quals.len() {
                    Some(CovariateValue::Int(u32::from(quals[read_offset])))
                } else {
                    None
                }
            }
            Self::Mapq => cache.mapq,
            Self::Cycle => Some(CovariateValue::Int(cycle)),
            Self::ReadNum => cache.read_num,
            Self::Strand => Some(cache.strand),
            Self::PairOrientation => cache.pair_orientation,
            Self::Isize => cache.isize_val,
            Self::Gc => cache.gc,
            Self::ReadBase => {
                let base = read_base_at(record, read_offset, is_reverse)?;
                Some(CovariateValue::Char(base))
            }
            Self::RefBase => {
                let base = region.ref_base(ref_pos);
                let oriented = if is_reverse { complement(base) } else { base };
                if is_valid_base(oriented) { Some(CovariateValue::Char(oriented)) } else { None }
            }
            Self::HpLen => {
                let hp = compute_hp_length(record, read_offset, is_reverse);
                Some(CovariateValue::Int(hp))
            }
            Self::PreDinuc => {
                // Previous read base (in sequencing order) + current reference base
                let prev = read_base_in_sequencing_order(record, read_offset, is_reverse, -1)?;
                let ref_b = region.ref_base(ref_pos);
                let ref_oriented = if is_reverse { complement(ref_b) } else { ref_b };
                if !is_valid_base(ref_oriented) {
                    return None;
                }
                Some(CovariateValue::Dinuc([prev, ref_oriented]))
            }
            Self::PostDinuc => {
                // Next read base (in sequencing order) + current reference base
                let next = read_base_in_sequencing_order(record, read_offset, is_reverse, 1)?;
                let ref_b = region.ref_base(ref_pos);
                let ref_oriented = if is_reverse { complement(ref_b) } else { ref_b };
                if !is_valid_base(ref_oriented) {
                    return None;
                }
                Some(CovariateValue::Dinuc([next, ref_oriented]))
            }
            Self::Context3bp => {
                let prev = read_base_in_sequencing_order(record, read_offset, is_reverse, -1)?;
                let next = read_base_in_sequencing_order(record, read_offset, is_reverse, 1)?;
                let ref_b = region.ref_base(ref_pos);
                let ref_oriented = if is_reverse { complement(ref_b) } else { ref_b };
                if !is_valid_base(ref_oriented) {
                    return None;
                }
                Some(CovariateValue::Trinuc([prev, ref_oriented, next]))
            }
            Self::Nm => {
                let nm = cache.nm_raw?;
                let ref_b = region.ref_base(ref_pos);
                let read_b = record.sequence().as_ref().get(read_offset).copied()?;
                // Subtract 1 if this base is itself a mismatch (so we get "other mismatches")
                let adjusted =
                    if read_b.to_ascii_uppercase() == ref_b { nm } else { nm.saturating_sub(1) };
                Some(CovariateValue::Int(adjusted))
            }
            Self::IndelLen => Some(CovariateValue::Int(indel_len)),
        }
    }

    /// Returns the index of this stratifier variant (based on its discriminant).
    fn index(self) -> usize {
        self as usize
    }
}

impl FromStr for Stratifier {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "all" => Ok(Self::All),
            "bq" => Ok(Self::Bq),
            "mapq" => Ok(Self::Mapq),
            "cycle" => Ok(Self::Cycle),
            "read_num" => Ok(Self::ReadNum),
            "strand" => Ok(Self::Strand),
            "pair_orientation" => Ok(Self::PairOrientation),
            "isize" => Ok(Self::Isize),
            "gc" => Ok(Self::Gc),
            "read_base" => Ok(Self::ReadBase),
            "ref_base" => Ok(Self::RefBase),
            "hp_len" => Ok(Self::HpLen),
            "pre_dinuc" => Ok(Self::PreDinuc),
            "post_dinuc" => Ok(Self::PostDinuc),
            "context_3bp" => Ok(Self::Context3bp),
            "nm" => Ok(Self::Nm),
            "indel_len" => Ok(Self::IndelLen),
            _ => bail!(
                "Unknown stratifier: '{s}'. Available: all, bq, mapq, cycle, read_num, \
                 strand, pair_orientation, isize, gc, read_base, ref_base, hp_len, \
                 pre_dinuc, post_dinuc, context_3bp, nm, indel_len"
            ),
        }
    }
}

impl fmt::Display for Stratifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::All => "all",
            Self::Bq => "bq",
            Self::Mapq => "mapq",
            Self::Cycle => "cycle",
            Self::ReadNum => "read_num",
            Self::Strand => "strand",
            Self::PairOrientation => "pair_orientation",
            Self::Isize => "isize",
            Self::Gc => "gc",
            Self::ReadBase => "read_base",
            Self::RefBase => "ref_base",
            Self::HpLen => "hp_len",
            Self::PreDinuc => "pre_dinuc",
            Self::PostDinuc => "post_dinuc",
            Self::Context3bp => "context_3bp",
            Self::Nm => "nm",
            Self::IndelLen => "indel_len",
        };
        write!(f, "{s}")
    }
}

/// A group of stratifiers that together produce a composite covariate key.
///
/// E.g. `["read_num", "cycle"]` produces keys like `["R1", "42"]`.
#[derive(Debug, Clone)]
pub struct StratifierGroup {
    /// The individual stratifiers in this group.
    pub stratifiers: Vec<Stratifier>,
    /// Precomputed display name, e.g. "read_num,cycle".
    pub name: String,
}

impl StratifierGroup {
    /// Parse a comma-separated stratifier group string (e.g., "read_num,cycle").
    ///
    /// Returns an error if any stratifier name is unknown. The parsed group
    /// contains the list of stratifiers and a precomputed flag indicating whether
    /// the group is a single stratifier (common case, optimized key type).
    fn parse(s: &str) -> Result<Self> {
        let stratifiers: Vec<Stratifier> =
            s.split(',').map(|part| part.trim().parse()).collect::<Result<_>>()?;
        if stratifiers.is_empty() {
            bail!("Empty stratifier group");
        }
        if stratifiers.len() > INLINE_STRATIFIERS {
            log::warn!("*************************************************************");
            log::warn!(
                "Stratifier group '{}' has {} stratifiers (max inline: {}).",
                s,
                stratifiers.len(),
                INLINE_STRATIFIERS
            );
            log::warn!("Groups with more than {INLINE_STRATIFIERS} stratifiers will see");
            log::warn!("a sudden and large performance degradation due to heap");
            log::warn!("allocation in the inner loop. Consider splitting into");
            log::warn!("smaller groups.");
            log::warn!("*************************************************************");
        }
        let name = stratifiers.iter().map(ToString::to_string).collect::<Vec<_>>().join(",");
        Ok(Self { stratifiers, name })
    }

    /// Assemble a covariate key from precomputed values in the cache.
    /// Returns None if any required stratifier value is None.
    fn key_from_cache(&self, base_cache: &BaseCovariateCache) -> Option<CovariateKey> {
        if self.stratifiers.len() == 1 {
            let val = base_cache.get(self.stratifiers[0])?;
            Some(CovariateKey::Single(val))
        } else {
            let mut key = SmallVec::new();
            for &strat in &self.stratifiers {
                key.push(base_cache.get(strat)?);
            }
            Some(CovariateKey::Composite(key))
        }
    }
}

/// A covariate key for a stratifier group.
///
/// Single-stratifier groups (the common case) avoid `SmallVec` overhead by
/// storing the value directly.  Composite groups use a `SmallVec` as before.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum CovariateKey {
    /// Key for a single-stratifier group — no `SmallVec` overhead.
    Single(CovariateValue),
    /// Key for a multi-stratifier group.
    Composite(SmallVec<[CovariateValue; INLINE_STRATIFIERS]>),
}

impl CovariateKey {
    /// Format a covariate key as comma-separated covariate values for TSV output.
    fn format(&self, interner: &StringInterner) -> String {
        match self {
            Self::Single(v) => v.format(interner),
            Self::Composite(vals) => {
                let mut s = String::new();
                for (i, v) in vals.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    s.push_str(&v.format(interner));
                }
                s
            }
        }
    }
}

/// A stack-allocated, Copy-able covariate value used as HashMap keys.
///
/// Avoids per-base String allocations in the hot loop. Each stratifier
/// returns one of these variants; composite keys use `SmallVec<[CovariateValue; INLINE_STRATIFIERS]>`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CovariateValue {
    /// An unsigned integer covariate (base quality, mapping quality, cycle, insert size, etc.).
    Int(u32),
    /// A single-character covariate (read base, reference base).
    Char(u8),
    /// A two-character covariate (dinucleotide context).
    Dinuc([u8; 2]),
    /// A three-character covariate (trinucleotide context).
    Trinuc([u8; 3]),
    /// An interned string covariate (read number, strand, pair orientation, "all").
    InternedStr(u16),
}

impl CovariateValue {
    /// Format a single covariate value to a string, using the interner to resolve `InternedStr`.
    fn format(self, interner: &StringInterner) -> String {
        match self {
            Self::Int(n) => n.to_string(),
            Self::Char(c) => String::from(c as char),
            Self::Dinuc(d) => format!("{}{}", d[0] as char, d[1] as char),
            Self::Trinuc(t) => {
                format!("{}{}{}", t[0] as char, t[1] as char, t[2] as char)
            }
            Self::InternedStr(id) => interner.resolve(id).to_string(),
        }
    }
}

/// Interns strings into u16 indices for compact, fast covariate keys.
/// Pre-populated with known static strings; supports dynamic strings at runtime.
struct StringInterner {
    strings: Vec<String>,
    index: HashMap<String, u16>,
}

impl StringInterner {
    /// Create a new interner pre-populated with known static covariate strings.
    fn new() -> Self {
        let mut interner = Self { strings: Vec::new(), index: HashMap::default() };
        // Pre-populate with all known static strings
        for s in ["all", "R1", "R2", "+", "-", "f1r2", "f2r1", "tandem"] {
            interner.intern(s);
        }
        interner
    }

    /// Intern a string, returning its u16 index. If already interned, returns existing index.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "interned string count is always << u16::MAX"
    )]
    fn intern(&mut self, s: &str) -> u16 {
        if let Some(&id) = self.index.get(s) {
            return id;
        }
        let id = self.strings.len() as u16;
        self.strings.push(s.to_string());
        self.index.insert(s.to_string(), id);
        id
    }

    /// Look up an already-interned string. Panics if not found.
    fn get(&self, s: &str) -> u16 {
        *self.index.get(s).expect("string not interned")
    }

    /// Resolve an interned id back to its string.
    fn resolve(&self, id: u16) -> &str {
        &self.strings[id as usize]
    }
}

/// Precomputed reference data for a genomic interval.
pub struct RegionContext {
    /// 0-based start position (inclusive) on the reference.
    pub start: u32,
    /// 0-based end position (exclusive) on the reference.
    pub end: u32,
    /// Uppercase reference bases for this region (length = end - start).
    pub bases: Vec<u8>,
    /// Variant sites within this region (set bits = excluded positions).
    pub excluded: BitVec,
}

impl RegionContext {
    /// Create a new region context by loading reference bases and variant sites.
    ///
    /// When `contig_cache` is provided and already holds the requested contig,
    /// the cached bases are reused instead of re-reading the FASTA.
    /// `contig_variant_mask` is the pre-loaded full-contig variant mask; the
    /// region `[start, end)` is sliced from it internally.
    fn new(
        contig_name: &str,
        start: u32,
        end: u32,
        fasta: &mut Fasta,
        contig_variant_mask: Option<&BitVec>,
        contig_cache: &mut ContigCache,
    ) -> Result<Self> {
        // Reload if this is a different contig than what's cached
        if contig_cache.name != contig_name {
            contig_cache.bases = fasta.load_contig(contig_name, true)?;
            contig_cache.name = contig_name.to_string();
        }

        let start_usize = start as usize;
        let end_usize = (end as usize).min(contig_cache.bases.len());
        let bases = contig_cache.bases[start_usize..end_usize].to_vec();

        // Slice variant mask for this region if available
        let excluded = if let Some(mask) = contig_variant_mask {
            mask[start_usize..end_usize].to_bitvec()
        } else {
            bitvec![0; bases.len()]
        };

        Ok(Self { start, end, bases, excluded })
    }

    /// Create a region context covering an entire contig, for the Collector path.
    ///
    /// Loads the full contig from the reference and builds a combined exclusion
    /// `BitVec` that marks positions outside intervals (if provided) and known
    /// variant sites from the pre-loaded variant mask.
    fn new_for_contig(
        contig_name: &str,
        contig_len: u32,
        ref_id: usize,
        fasta: &mut Fasta,
        contig_variant_mask: Option<&BitVec>,
        intervals: Option<&crate::intervals::Intervals>,
        contig_cache: &mut ContigCache,
    ) -> Result<Self> {
        // Load reference bases (using contig cache)
        if contig_cache.name != contig_name {
            contig_cache.bases = fasta.load_contig(contig_name, true)?;
            contig_cache.name = contig_name.to_string();
        }
        let bases = contig_cache.bases[..contig_len as usize].to_vec();

        // Build combined exclusion BitVec
        let mut excluded = if let Some(ivs) = intervals {
            !ivs.contig_bitvec(ref_id) // positions NOT in intervals are excluded
        } else {
            bitvec![0; bases.len()] // nothing excluded
        };

        // OR in pre-loaded variant mask
        if let Some(mask) = contig_variant_mask {
            excluded |= mask;
        }

        Ok(Self { start: 0, end: contig_len, bases, excluded })
    }

    /// Get the reference base at a genomic position (0-based).
    fn ref_base(&self, ref_pos: u32) -> u8 {
        let offset = (ref_pos - self.start) as usize;
        self.bases.get(offset).copied().unwrap_or(b'N')
    }

    /// Check if a genomic position is excluded (known variant site).
    fn is_excluded(&self, ref_pos: u32) -> bool {
        let offset = (ref_pos - self.start) as usize;
        self.excluded.get(offset).is_some_and(|b| *b)
    }

    /// Check if a genomic position is within this region.
    fn contains(&self, ref_pos: u32) -> bool {
        ref_pos >= self.start && ref_pos < self.end
    }
}

/// Cache for the most recently loaded reference contig, avoiding redundant FASTA
/// I/O when multiple intervals fall on the same contig.
struct ContigCache {
    name: String,
    bases: Vec<u8>,
}

/// Cached per-read covariate values that are constant across all bases in a read.
///
/// Computing these once per record (instead of once per base per stratifier group)
/// eliminates redundant work in the hot loop — especially for `gc` which iterates
/// the entire read sequence.
struct ReadLevelCache {
    /// First-of-pair vs second-of-pair ("R1" or "R2"), or None if unpaired.
    read_num: Option<CovariateValue>,
    /// Alignment direction ("+" or "-").
    strand: CovariateValue,
    /// Pair orientation ("f1r2", "f2r1", "tandem"), or None if not determinable.
    pair_orientation: Option<CovariateValue>,
    /// Absolute insert size, capped or excluded by --max-isize, or None if excluded.
    isize_val: Option<CovariateValue>,
    /// GC content of the read as integer percentage (0-100), or None if sequence empty.
    gc: Option<CovariateValue>,
    /// Mapping quality, or None if not available.
    mapq: Option<CovariateValue>,
    /// Raw NM tag value; per-base adjustment is still done in `Stratifier::covariate`.
    nm_raw: Option<u32>,
}

impl ReadLevelCache {
    /// Build the cache from a record.  Called once per record in `process_record`.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "GC percentage is always 0-100, fits in u32"
    )]
    fn new(
        record: &RecordBuf,
        interner: &StringInterner,
        max_isize: u32,
        picard_compat: bool,
    ) -> Self {
        let flags = record.flags();
        let is_reverse = flags.is_reverse_complemented();

        let read_num = if !flags.is_segmented() {
            None
        } else if flags.is_first_segment() {
            Some(CovariateValue::InternedStr(interner.get("R1")))
        } else {
            Some(CovariateValue::InternedStr(interner.get("R2")))
        };

        let strand = if is_reverse {
            CovariateValue::InternedStr(interner.get("-"))
        } else {
            CovariateValue::InternedStr(interner.get("+"))
        };

        let pair_orientation = get_pair_orientation(record).map(|po| {
            CovariateValue::InternedStr(interner.get(match po {
                PairOrientation::Fr => {
                    if flags.is_first_segment() == is_reverse {
                        "f2r1"
                    } else {
                        "f1r2"
                    }
                }
                PairOrientation::Rf => {
                    if flags.is_first_segment() == is_reverse {
                        "f1r2"
                    } else {
                        "f2r1"
                    }
                }
                PairOrientation::Tandem => "tandem",
            }))
        });

        let raw_isize = record.template_length().unsigned_abs();
        let isize_val = if picard_compat {
            // TEMPORARY (picard-compat): cap at max_isize instead of excluding
            Some(CovariateValue::Int(raw_isize.min(max_isize)))
        } else if raw_isize > max_isize {
            None // exclude from isize stratifier
        } else {
            Some(CovariateValue::Int(raw_isize))
        };

        let gc = {
            let seq = record.sequence().as_ref();
            if seq.is_empty() {
                None
            } else {
                let gc_count = seq
                    .iter()
                    .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                    .count();
                let pct = (gc_count * 100 + seq.len() / 2) / seq.len();
                Some(CovariateValue::Int(pct as u32))
            }
        };

        let mapq = record.mapping_quality().map(|m| CovariateValue::Int(u32::from(m.get())));

        let nm_raw = get_nm_tag(record);

        Self { read_num, strand, pair_orientation, isize_val, gc, mapq, nm_raw }
    }
}

/// Pre-computed covariate values for all stratifiers at a single base position.
/// Computed once per base, reused across all stratifier groups.
struct BaseCovariateCache {
    values: [Option<CovariateValue>; Stratifier::COUNT],
}

impl BaseCovariateCache {
    /// Populate the read-level slots from the `ReadLevelCache`. Called once per record.
    /// Slots not in `read_level` are left uninitialized (they'll be overwritten per-base).
    fn from_read_level(
        read_level: &[Stratifier],
        cache: &ReadLevelCache,
        interner: &StringInterner,
    ) -> Self {
        let mut values = [None; Stratifier::COUNT];
        for &strat in read_level {
            values[strat.index()] = match strat {
                Stratifier::All => Some(CovariateValue::InternedStr(interner.get("all"))),
                Stratifier::ReadNum => cache.read_num,
                Stratifier::Strand => Some(cache.strand),
                Stratifier::PairOrientation => cache.pair_orientation,
                Stratifier::Isize => cache.isize_val,
                Stratifier::Gc => cache.gc,
                Stratifier::Mapq => cache.mapq,
                _ => None, // not a read-level stratifier
            };
        }
        Self { values }
    }

    /// Update only the base-level slots. Called once per base position.
    /// Read-level slots are preserved from `from_read_level`.
    #[expect(clippy::too_many_arguments)]
    fn populate_base_level(
        &mut self,
        base_level: &[Stratifier],
        record: &RecordBuf,
        read_offset: usize,
        ref_pos: u32,
        region: &RegionContext,
        indel_len: u32,
        cycle: u32,
        cache: &ReadLevelCache,
        interner: &StringInterner,
    ) {
        for &strat in base_level {
            self.values[strat.index()] = strat.covariate(
                record,
                read_offset,
                ref_pos,
                region,
                indel_len,
                cycle,
                cache,
                interner,
            );
        }
    }

    /// Look up the precomputed value for the given stratifier.
    fn get(&self, strat: Stratifier) -> Option<CovariateValue> {
        self.values[strat.index()]
    }
}

/// Holds a stratifier group and its per-covariate-key accumulators for mismatch,
/// indel, and overlap metrics.
struct GroupAccumulators {
    group: StratifierGroup,
    accums: HashMap<CovariateKey, CombinedAccum>,
}

impl GroupAccumulators {
    fn new(group: StratifierGroup) -> Self {
        Self { group, accums: HashMap::with_capacity_and_hasher(256, rustc_hash::FxBuildHasher) }
    }
}

/// Combined accumulator holding mismatch, indel, and overlap data for a single
/// covariate bin.  Using a single struct per key means one HashMap lookup instead
/// of three per base per group.
#[derive(Debug, Default, Clone)]
struct CombinedAccum {
    mismatch: MismatchAccum,
    indel: IndelAccum,
    overlap: OverlapAccum,
}

/// Accumulates mismatch error counts for a single covariate bin.
#[derive(Debug, Default, Clone)]
struct MismatchAccum {
    total_bases: u64,
    error_bases: u64,
}

/// Accumulates overlapping-read mismatch counts for a single covariate bin.
#[derive(Debug, Default, Clone)]
struct OverlapAccum {
    overlapping_read_bases: u64,
    bases_mismatching_ref_and_mate: u64,
    bases_matching_mate_but_not_ref: u64,
    bases_in_three_way_disagreement: u64,
}

/// Accumulates indel error counts for a single covariate bin.
#[derive(Debug, Default, Clone)]
struct IndelAccum {
    total_bases: u64,
    num_insertions: u64,
    num_inserted_bases: u64,
    num_deletions: u64,
    num_deleted_bases: u64,
}

/// Result of presenting a read to the mate buffer for overlap classification.
enum MateResult {
    /// The read was buffered — its mate may arrive later for overlap comparison.
    /// No records need processing now.
    Buffered,
    /// The read should be processed alone (no overlap data).
    Single(RecordBuf),
    /// The read and its buffered mate should be processed together with overlap data.
    Pair(RecordBuf, RecordBuf),
}

/// Buffer for storing reads that may overlap their mate pair.
///
/// Reads are keyed by read name. When a read arrives, `accept_read` determines
/// whether to buffer it (mate expected later), return it alone (no overlap), or
/// return it paired with a previously buffered mate.
struct MateBuffer {
    buffer: HashMap<BString, RecordBuf>,
}

impl MateBuffer {
    fn new() -> Self {
        Self { buffer: HashMap::with_capacity_and_hasher(4096, rustc_hash::FxBuildHasher) }
    }

    /// Accept a read and decide what to do with it:
    /// 1. If the read can't overlap its mate (unpaired, unmapped mate, different
    ///    contig) → return `Single` immediately without checking the buffer.
    /// 2. If the mate is already buffered → remove it and return `Pair(read, mate)`.
    /// 3. If the mate's start falls within this read's reference span (potential
    ///    overlap) → buffer it, return `Buffered`.
    /// 4. Otherwise → return `Single(read)` for processing without overlap.
    fn accept_read(&mut self, record: RecordBuf) -> MateResult {
        let flags = record.flags();

        // Quick check: can this read possibly overlap its mate?
        if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
            return MateResult::Single(record);
        }
        if record.reference_sequence_id() != record.mate_reference_sequence_id() {
            return MateResult::Single(record);
        }

        // Check if the mate is already buffered (zero-allocation lookup via &BStr)
        if let Some(name) = record.name() {
            let name_ref: &BStr = BStr::new(name);
            if self.buffer.contains_key(name_ref) {
                let mate = self.buffer.remove(name_ref).unwrap();
                return MateResult::Pair(record, mate);
            }
        }

        // Buffer if mate's start falls within this read's reference span
        let read_start = record.alignment_start().map_or(0, |p| p.get());
        let read_end = record.alignment_end().map_or(0, |e| e.get());
        let mate_start = record.mate_alignment_start().map_or(0, |p| p.get());

        if mate_start >= read_start && mate_start <= read_end {
            let name = BString::from(record.name().map_or_else(Vec::new, |n| <[u8]>::to_vec(n)));
            self.buffer.insert(name, record);
            return MateResult::Buffered;
        }

        MateResult::Single(record)
    }

    /// Remove and return buffered reads whose expected mate position is behind `current_pos`.
    /// These reads' mates won't arrive, so they need to be processed without overlap data.
    #[expect(clippy::cast_possible_truncation, reason = "1-based genomic positions fit in u32")]
    fn flush_behind(&mut self, current_ref_id: usize, current_pos: u32) -> Vec<RecordBuf> {
        let keys: Vec<BString> = self
            .buffer
            .iter()
            .filter(|(_, record)| {
                let mate_ref_id = record.mate_reference_sequence_id().unwrap_or(usize::MAX);
                let mate_pos = record.mate_alignment_start().map_or(0, |p| (p.get() - 1) as u32);
                mate_ref_id < current_ref_id
                    || (mate_ref_id == current_ref_id && mate_pos < current_pos)
            })
            .map(|(k, _)| k.clone())
            .collect();
        keys.into_iter().filter_map(|k| self.buffer.remove(&k)).collect()
    }

    /// Drain and return all remaining buffered reads.
    fn drain_all(&mut self) -> Vec<RecordBuf> {
        self.buffer.drain().map(|(_, v)| v).collect()
    }
}

// ─── Module-level functions ─────────────────────────────────────────────────────

/// Complement a base.
const fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Returns true if the base is a standard DNA base (A, C, G, T).
const fn is_valid_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// Get the read base at the given offset, oriented to sequencing direction.
fn read_base_at(record: &RecordBuf, read_offset: usize, is_reverse: bool) -> Option<u8> {
    let seq = record.sequence().as_ref();
    let base = *seq.get(read_offset)?;
    let oriented = if is_reverse { complement(base) } else { base };
    if is_valid_base(oriented) { Some(oriented) } else { None }
}

/// Get a read base relative to current position in sequencing order.
/// `delta` is -1 for previous, +1 for next in sequencing direction.
#[expect(clippy::cast_sign_loss, reason = "array_delta is checked non-negative before cast")]
fn read_base_in_sequencing_order(
    record: &RecordBuf,
    read_offset: usize,
    is_reverse: bool,
    delta: i32,
) -> Option<u8> {
    // In sequencing order, "previous" means:
    // - For forward reads: read_offset - 1
    // - For reverse reads: read_offset + 1 (since sequencing is right-to-left in alignment coords)
    let array_delta = if is_reverse { -delta } else { delta };
    let target = if array_delta >= 0 {
        read_offset.checked_add(array_delta as usize)?
    } else {
        read_offset.checked_sub(array_delta.unsigned_abs() as usize)?
    };
    let seq = record.sequence().as_ref();
    let base = *seq.get(target)?;
    let oriented = if is_reverse { complement(base) } else { base };
    if is_valid_base(oriented) { Some(oriented) } else { None }
}

/// Compute homopolymer run length preceding the current base in sequencing direction.
///
/// Matches Picard's `stratifyHomopolymerLength`: counts consecutive identical
/// bases in the read before this position in the direction of sequencing.
fn compute_hp_length(record: &RecordBuf, read_offset: usize, is_reverse: bool) -> u32 {
    let seq = record.sequence().as_ref();
    let len = seq.len();

    if is_reverse {
        // Sequencing is right-to-left in alignment coords; "preceding" means scanning forward
        if read_offset + 1 >= len {
            return 0;
        }
        let anchor = seq[read_offset + 1];
        let mut count = 0u32;
        let mut i = read_offset + 1;
        while i < len && seq[i] == anchor {
            count += 1;
            i += 1;
        }
        // Picard returns (i - read_offset - 1), which is count
        count
    } else {
        // Forward reads: scan backward
        if read_offset == 0 {
            return 0;
        }
        let anchor = seq[read_offset - 1];
        let mut count = 0u32;
        let mut i = read_offset;
        while i > 0 {
            i -= 1;
            if seq[i] != anchor {
                break;
            }
            count += 1;
        }
        count
    }
}

/// Extract the NM tag value from a record.
fn get_nm_tag(record: &RecordBuf) -> Option<u32> {
    use noodles::sam::alignment::record_buf::data::field::Value;
    let data = record.data();
    let tag: [u8; 2] = *b"NM";
    let value = data.get(&tag)?;
    match value {
        Value::UInt8(v) => Some(u32::from(*v)),
        Value::UInt16(v) => Some(u32::from(*v)),
        Value::UInt32(v) => Some(*v),
        Value::Int8(v) => u32::try_from(*v).ok(),
        Value::Int16(v) => u32::try_from(*v).ok(),
        Value::Int32(v) => u32::try_from(*v).ok(),
        _ => None,
    }
}

/// Compute a phred-scaled Q-score from error count and total.
#[expect(
    clippy::cast_precision_loss,
    reason = "f64 precision is sufficient for error rate calculations"
)]
fn q_score(errors: u64, total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    if errors == 0 {
        return 99.0;
    }
    let rate = errors as f64 / total as f64;
    (-10.0 * rate.log10()).min(99.0)
}

/// Build a sorted vec of (ref_pos, read_base) for aligned positions passing
/// the base quality threshold. Used for overlap comparison — both reads traverse
/// ref positions in order, so the mate's vec can be scanned with a monotonically
/// advancing index. Bases with BQ below `min_bq` are excluded so that overlap
/// is only detected at positions where both reads pass quality filters.
#[expect(
    clippy::cast_possible_truncation,
    reason = "genomic positions and CIGAR op lengths fit in u32"
)]
fn build_aligned_bases(record: &RecordBuf, min_bq: u8) -> Vec<(u32, u8)> {
    let Some(alignment_start) = record.alignment_start() else { return Vec::new() };
    let alignment_start = (alignment_start.get() - 1) as u32;
    let seq = record.sequence().as_ref();
    let quals = record.quality_scores().as_ref();
    let mut bases = Vec::with_capacity(seq.len());
    let mut ref_pos = alignment_start;
    let mut read_pos: usize = 0;

    for op_result in record.cigar().iter() {
        let Ok(op) = op_result else { return bases };
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op.len() {
                    let bq = quals.get(read_pos).copied().unwrap_or(0);
                    if bq >= min_bq
                        && let Some(&base) = seq.get(read_pos)
                    {
                        bases.push((ref_pos, base.to_ascii_uppercase()));
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += op.len();
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += op.len() as u32;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    bases
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use strum::IntoEnumIterator;

    #[test]
    fn test_stratifier_count_and_indices() {
        // EnumCount::COUNT is derived from the actual enum definition by strum,
        // so adding a variant automatically updates it. This test verifies that
        // index() maps every variant to a unique slot in 0..COUNT.
        let mut seen = [false; Stratifier::COUNT];
        for strat in Stratifier::iter() {
            let idx = strat.index();
            assert!(idx < Stratifier::COUNT, "index {idx} out of range for {strat}");
            assert!(!seen[idx], "duplicate index {idx} for {strat}");
            seen[idx] = true;
        }
        assert!(seen.iter().all(|&s| s), "not all indices covered");
    }
}
