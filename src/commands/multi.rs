//! The `multi` command: one pass over the input, many collectors, run in parallel.
//!
//! ## Threading model
//!
//! multi receives the toolkit-wide `--threads` budget and lays it out in
//! [`plan_multi`] across three roles: **parallel read threads** (bgzf inflate
//! for BAM / htslib decode pool for CRAM), a **dispatch thread** (reads,
//! extracts, and fans records out to the workers), and **compute worker
//! groups** (each owns a subset of collectors, partitioned by
//! [`partition_collectors`] into `min(workers, n_collectors)` groups). The
//! layout is tuned per input format and was chosen from benchmarking (see
//! [`plan_multi`]): BAM is compute-bound (cheap inflate) so it favors workers;
//! CRAM is decode-bound (expensive decode) so it favors read threads and stays
//! serial-main at low budgets; SAM can't decode in parallel so it's all
//! workers. A budget of 1 — and CRAM's low-budget serial-main layouts — skip
//! the pipeline entirely (see "Single-threaded path" below).
//!
//! ## Collector→worker affinity
//!
//! Each collector is assigned to exactly one worker group, balanced by the
//! collector's [`Collector::cost_hint`] (greedy longest-processing-time),
//! so the heaviest collectors tend to get their own worker. Because a
//! collector is owned by a single thread there is **no per-collector
//! mutex** — no lock contention — and because the reader ships batches to
//! each group in a fixed order, every collector sees records in file order,
//! so output is deterministic regardless of thread count.
//!
//! Reader and workers are connected through these channels:
//!
//! - `batch_pool` — an unbounded kanal channel of empty `Vec<RikerRecord>`
//!   slots the reader pulls from (and which [`RecyclableBatch::drop`]
//!   returns to). Slots are pre-allocated as the record variant the reader
//!   writes into: `RikerRecord::Bam` for BAM, `RikerRecord::Fallback` for
//!   SAM, `RikerRecord::Htslib` for CRAM.
//! - one **ordered channel per worker group** — a bounded kanal queue
//!   of `Arc<RecyclableBatch>`. The reader wraps each filled batch in an
//!   `Arc` and sends a clone to every group's channel; the owning worker
//!   blocks on `recv()` (kanal briefly spins, then parks — no busy-polling).
//!   Bounded one batch above the pool's in-flight max so the recycling pool
//!   is the practical backpressure.
//! - `return` — the kanal `return_tx`/`return_rx` captured inside each
//!   `RecyclableBatch`. When the last `Arc` reference drops (every group
//!   has finished the batch), [`RecyclableBatch::drop`] sends the inner
//!   `Vec` back to the pool.
//!
//! Every format reads in place via [`AlignmentReader::fill_record`]. The
//! reader also consults the union of every active collector's
//! [`Collector::field_needs`] once up front and passes it down. On BAM and
//! CRAM (where aux decode is lazy) this gates which decoders run per
//! record, so a collector set that never reads aux tags pays zero
//! aux-decode cost. SAM decodes eagerly inside noodles, so the union is
//! informational there.
//!
//! When the reader hits EOF it drops the group senders, each channel
//! closes, and each worker hands its collectors back to the main thread,
//! which finalizes every group ([`Collector::finish`]) only once all
//! workers and the reader have joined cleanly.
//!
//! ## Error propagation
//!
//! Channel disconnection handles the happy shutdown path. An `AtomicBool`
//! poison flag handles the sad one.
//!
//! - A worker that gets `Err` from `accept_multiple` sets the flag and
//!   returns the error, so the reader aborts on its next send and sibling
//!   workers stop at their next batch instead of processing the rest of the
//!   in-flight work.
//! - If the reader itself errors it returns `Err`; it sets the poison flag
//!   before dropping the group senders so any in-flight worker stops
//!   promptly rather than processing a truncated stream.
//!
//! Errors are surfaced by `handle.join()` on the main thread — the reader's
//! error wins if present, otherwise the first worker error. Finalization
//! runs on the main thread *after* every worker and the reader have joined
//! without error, so **a failed run leaves no partial output on disk**.
//!
//! All `poison` flag accesses use `Ordering::Relaxed`. Relaxed is sufficient
//! because the flag is only a best-effort "stop early" signal and no longer
//! guards output: no-partial-output is enforced by finalizing after the
//! worker/reader joins — `join` is the happens-before edge that orders every
//! `accept_multiple` before any `finish()` — not by the flag.
//!
//! ## Single-threaded path
//!
//! When [`plan_multi`] asks for zero compute workers — a budget of 1, or CRAM
//! at budget 2-3 where the scarce threads are better spent decoding — the run
//! collapses to [`run_single_threaded`]: no channels, no worker pool, just a
//! serial loop on the main thread that drives the same reader API (in-place
//! fills via [`AlignmentReader::fill_record`]), while the reader's decode pool
//! (if any) inflates/decodes in parallel behind it.

use std::fmt;
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Result, anyhow};
use clap::{Args, ValueEnum};
use kanal::{Receiver, Sender};
use noodles::sam::Header;

use crate::collector::Collector;
use crate::commands::alignment::{AlignmentCollector, MultiAlignmentOptions};
use crate::commands::basic::BasicCollector;
use crate::commands::command::{Command, ThreadPlan, resolve_threads};
use crate::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use crate::commands::error::{ErrorCollector, MultiErrorOptions};
use crate::commands::gcbias::{GcBiasCollector, MultiGcBiasOptions};
use crate::commands::hybcap::{HybCapCollector, MultiHybCapOptions};
use crate::commands::isize::{InsertSizeCollector, MultiIsizeOptions};
use crate::commands::wgs::{MultiWgsOptions, WgsCollector};
use crate::fasta::Fasta;
use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::{AlignmentFormat, AlignmentReader, detect_format};
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};

/// Number of records per batch sent through the work queue.
///
/// Picked at 128 from a `batch ∈ {32,48,64,128} × pool ∈ {16,32,48,64} ×
/// threads ∈ {2,4,8}` sweep on M1 Max and Graviton 3. Smaller batches
/// (e.g. 32) win marginally at `--threads 2` but cause heavy sys-time
/// blowup at higher thread counts: each batch costs a channel send/recv per
/// worker group plus the reader's fan-out clone, so small batches multiply
/// that per-record overhead. 128 keeps each `accept_multiple` call long
/// enough that the channel overhead stays cheap, while the working set
/// (`128 × 16 = 2048` records, ~1.5 MB for typical 150 bp WGS reads) still
/// fits comfortably in L2 on M1 Max and shared L3 on Graviton 3.
const BATCH_SIZE: usize = 128;

/// Number of pre-allocated batches in the recycling pool.
///
/// The pool is the sole backpressure between the reader and pool workers:
/// at most this many batches can be in flight at any instant (each batch
/// holds a real `Vec<RikerRecord>` from the pool). Picked at 16 from the
/// same sweep — keeps the working set small enough to stay cache-resident.
const NUM_BATCHES_POOLED: usize = 16;

// ─── Multi command struct ────────────────────────────────────────────────────

/// Run multiple metric collectors in a single BAM pass.
///
/// Reads the BAM file once and dispatches records to all selected
/// collectors in parallel, avoiding the overhead of multiple passes.
/// Use this when you need several metric types from the same BAM.
/// Collectors that require a reference (gcbias, wgs) need --reference.
/// The hybcap collector requires --hybcap::targets and --hybcap::baits.
///
/// Output files depend on the selected --tools and use the <prefix> from
/// -o/--output. See the help for each individual tool for the full list
/// of output files produced.
#[derive(Args, Debug, Clone)]
#[command(
    long_about,
    after_long_help = "\
Examples:
  riker multi -i input.bam -o out_prefix -r ref.fa
  riker multi -i input.bam -o out_prefix -r ref.fa --tools alignment basic isize
  riker multi -i input.bam -o out_prefix -r ref.fa --threads 4
  riker multi -i input.bam -o out_prefix --tools hybcap --hybcap::targets t.bed --hybcap::baits b.bed"
)]
pub struct Multi {
    #[command(flatten)]
    pub input: InputOptions,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub reference: OptionalReferenceOptions,

    /// Tools to run. Defaults to all except hybcap (which requires targets/baits).
    /// The wgs and gcbias tools require --reference.
    #[arg(
        long,
        num_args(1..),
        default_values_t = [CollectorKind::Alignment, CollectorKind::Basic, CollectorKind::Isize],
        help_heading = "Multi Command Options"
    )]
    pub tools: Vec<CollectorKind>,

    // Per-tool options
    #[command(flatten)]
    pub alignment_opts: MultiAlignmentOptions,
    #[command(flatten)]
    pub error_opts: MultiErrorOptions,
    #[command(flatten)]
    pub gcbias_opts: MultiGcBiasOptions,
    #[command(flatten)]
    pub hybcap_opts: MultiHybCapOptions,
    #[command(flatten)]
    pub isize_opts: MultiIsizeOptions,
    #[command(flatten)]
    pub wgs_opts: MultiWgsOptions,
}

impl Multi {
    /// Build the list of collectors based on the deduplicated kinds.
    fn build_collectors(
        &self,
        kinds: &[CollectorKind],
        header: &Header,
    ) -> Result<Vec<Box<dyn Collector>>> {
        let mut collectors: Vec<Box<dyn Collector>> = Vec::new();
        for kind in kinds {
            match kind {
                CollectorKind::Alignment => {
                    let opts = self.alignment_opts.clone().validate()?;
                    collectors.push(Box::new(AlignmentCollector::new(
                        &self.input.input,
                        &self.output.output,
                        self.reference.reference.clone(),
                        &opts,
                    )));
                }
                CollectorKind::Basic => {
                    collectors.push(Box::new(BasicCollector::new(
                        &self.input.input,
                        &self.output.output,
                    )));
                }
                CollectorKind::Error => {
                    let ref_path = self.reference.reference.as_ref().unwrap();
                    let reference = Fasta::from_path(ref_path)?;
                    let mut error_opts = self.error_opts.clone();
                    // Fall back to global --reference if --error::reference not set
                    if error_opts.error_reference.is_none() {
                        error_opts.error_reference = Some(ref_path.clone());
                    }
                    let opts = error_opts.validate()?;
                    collectors.push(Box::new(ErrorCollector::new(
                        &self.input.input,
                        &self.output.output,
                        reference,
                        &opts,
                    )?));
                }
                CollectorKind::GcBias => {
                    let ref_path = self.reference.reference.as_ref().unwrap();
                    let reference = Fasta::from_path(ref_path)?;
                    let opts = self.gcbias_opts.clone().validate()?;
                    collectors.push(Box::new(GcBiasCollector::new(
                        &self.input.input,
                        &self.output.output,
                        reference,
                        &opts,
                    )));
                }
                CollectorKind::HybCap => {
                    let opts = self.hybcap_opts.clone().validate()?;
                    let sample = derive_sample(&self.input.input, header);
                    let fasta = self
                        .reference
                        .reference
                        .as_ref()
                        .map(|p| Fasta::from_path(p))
                        .transpose()?;
                    collectors.push(Box::new(HybCapCollector::new(
                        &self.output.output,
                        fasta,
                        sample,
                        &opts,
                    )));
                }
                CollectorKind::Isize => {
                    let opts = self.isize_opts.clone().validate()?;
                    collectors.push(Box::new(InsertSizeCollector::new(
                        &self.input.input,
                        &self.output.output,
                        &opts,
                    )));
                }
                CollectorKind::Wgs => {
                    let ref_path = self.reference.reference.as_ref().unwrap();
                    let reference = Fasta::from_path(ref_path)?;
                    let opts = self.wgs_opts.clone().validate()?;
                    collectors.push(Box::new(WgsCollector::new(
                        &self.input.input,
                        &self.output.output,
                        reference,
                        &opts,
                    )?));
                }
            }
        }
        Ok(collectors)
    }
}

impl Command for Multi {
    /// # Errors
    /// Returns an error if the BAM file cannot be read or any collector fails.
    fn execute(&self, threads: Option<u8>) -> Result<()> {
        let total = resolve_threads(threads, self.default_threads());

        // Deduplicate the collector list while preserving order.
        let mut seen = Vec::new();
        for kind in &self.tools {
            if !seen.contains(kind) {
                seen.push(*kind);
            }
        }

        // Validate required inputs for selected collectors.
        for kind in &seen {
            match kind {
                CollectorKind::Error if self.reference.reference.is_none() => {
                    return Err(anyhow!("Error collector requires --reference"));
                }
                CollectorKind::GcBias if self.reference.reference.is_none() => {
                    return Err(anyhow!("GC bias collector requires --reference"));
                }
                CollectorKind::Wgs if self.reference.reference.is_none() => {
                    return Err(anyhow!("WGS collector requires --reference"));
                }
                _ => {}
            }
        }

        // Plan the thread layout from the already-deduped tool count, then
        // open the reader (also needed to build interval maps for WGS) with
        // its share of decode threads.
        let plan = plan_multi(total, seen.len().max(1), InputKind::of(&self.input.input));
        let reader = AlignmentReader::open(
            &self.input.input,
            self.reference.reference.as_deref(),
            plan.decode_threads,
        )?;

        let collectors = self.build_collectors(&seen, reader.header())?;

        // `compute_workers == 0` is serial-main mode (collectors run on the main
        // thread, decode offloaded to the reader's pool); `>= 1` spins the
        // dispatch + worker pipeline.
        if plan.compute_workers == 0 {
            run_single_threaded(reader, collectors)?;
        } else {
            run_parallel(reader, collectors, plan.compute_workers)?;
        }

        Ok(())
    }

    /// multi parallelizes by default. Absent an explicit budget it uses up to
    /// four threads (its useful parallelism saturates around there for the
    /// typical collector set), clamped to the machine's core count.
    fn default_threads(&self) -> NonZero<usize> {
        let cores = std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN);
        cores.min(NonZero::new(4).expect("4 is non-zero"))
    }
}

// ─── Collector kinds ─────────────────────────────────────────────────────────

/// Available collector kinds for the multi command.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CollectorKind {
    /// Alignment summary metrics.
    Alignment,
    /// Basic QC metrics (base distribution, mean quality, quality distribution).
    Basic,
    /// Base-level error metrics (mismatch, overlap, indel).
    #[value(name = "error")]
    Error,
    /// GC bias metrics.
    #[value(name = "gcbias")]
    GcBias,
    /// Hybrid capture (bait/target) metrics.
    #[value(name = "hybcap")]
    HybCap,
    /// Insert size distribution metrics.
    Isize,
    /// Whole-genome sequencing coverage metrics.
    Wgs,
}

impl fmt::Display for CollectorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CollectorKind::Alignment => write!(f, "alignment"),
            CollectorKind::Basic => write!(f, "basic"),
            CollectorKind::Error => write!(f, "error"),
            CollectorKind::GcBias => write!(f, "gcbias"),
            CollectorKind::HybCap => write!(f, "hybcap"),
            CollectorKind::Isize => write!(f, "isize"),
            CollectorKind::Wgs => write!(f, "wgs"),
        }
    }
}

// ─── Threading helpers ───────────────────────────────────────────────────────

/// Input format, coarsened to what drives the thread layout: BAM (cheap
/// parallel bgzf inflate), CRAM (expensive, well-parallelizing htslib decode),
/// SAM (single-threaded text parse — no decode threads to hand out).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputKind {
    Sam,
    Bam,
    Cram,
}

impl InputKind {
    /// Classify by extension, reusing the reader's [`detect_format`]. The
    /// reader does the authoritative detection at open time; this only needs
    /// to be right enough to plan threads, so an unrecognized extension falls
    /// back to `Bam` (a misclassified SAM/BAM just means the reader ignores a
    /// decode count).
    fn of(path: &Path) -> Self {
        match detect_format(path) {
            Ok(AlignmentFormat::Cram) => Self::Cram,
            Ok(AlignmentFormat::Sam | AlignmentFormat::GzippedSam) => Self::Sam,
            Ok(AlignmentFormat::Bam) | Err(_) => Self::Bam,
        }
    }
}

/// Lay out a total thread budget for `multi` across reader-decode threads and
/// compute worker groups. `compute_workers == 0` means serial-main mode
/// (collectors run on the main thread, decode offloaded to the reader's pool);
/// `>= 1` means the dispatch + worker pipeline (with one uncounted dispatch
/// thread). Workers are capped at `n_tools` (a collector can't use more groups
/// than there are collectors).
///
/// The layout is tuned from benchmarking (12x WGS BAM + its CRAM 3.1
/// normal/archive transcodes, wgs panel). The formats differ because the pole
/// differs: **BAM** inflate is cheap so multi is compute-bound — favor workers,
/// with just enough read threads to keep the dispatcher fed. **CRAM** decode is
/// expensive and parallelizes well — favor read threads, staying serial-main at
/// low budgets so a scarce thread decodes rather than idling as a worker. **SAM**
/// can't decode in parallel, so everything goes to workers.
fn plan_multi(total: NonZero<usize>, n_tools: usize, kind: InputKind) -> ThreadPlan {
    let cap = n_tools.max(1);
    let serial = ThreadPlan { decode_threads: 0, compute_workers: 0 };
    let budget = total.get();
    match (kind, budget) {
        // Budget of 1 is always a single serial pass.
        (_, 1) => serial,
        // budget = 2, 3: dispatcher is counted; small per-format special cases.
        (InputKind::Cram, 2) => ThreadPlan { decode_threads: 1, compute_workers: 0 },
        (InputKind::Cram, 3) => ThreadPlan { decode_threads: 2, compute_workers: 0 },
        (_, 2) => ThreadPlan { decode_threads: 0, compute_workers: 1 },
        (_, 3) => ThreadPlan { decode_threads: 0, compute_workers: 2.min(cap) },
        // budget >= 4: the dispatcher becomes free; split the budget into readers
        // + workers, format-tuned by which dimension gets the odd thread.
        (InputKind::Sam, _) => ThreadPlan { decode_threads: 0, compute_workers: budget.min(cap) },
        (InputKind::Bam, _) => {
            let workers = budget.div_ceil(2).min(cap);
            ThreadPlan { decode_threads: budget - workers, compute_workers: workers }
        }
        (InputKind::Cram, _) => {
            // floor(budget/2) workers (readers get the odd thread); any budget freed
            // by the worker cap goes to reads, which CRAM decode makes good use of.
            let workers = (budget / 2).min(cap);
            ThreadPlan { decode_threads: budget - workers, compute_workers: workers }
        }
    }
}

/// A batch of records shared across collector channels.
///
/// Wrapping the records in `RecyclableBatch` lets us send the inner
/// `Vec<RikerRecord>` back to the reader's pool when the last `Arc` reference
/// drops (i.e. the last collector has finished with it), so the reader can
/// reuse the pre-allocated record slots on the next read.
type Batch = Arc<RecyclableBatch>;

/// Owns a `Vec<RikerRecord>` of capacity `BATCH_SIZE` plus a count of
/// valid records. On drop, the inner `Vec` is returned to the reader's
/// pool so its allocations can be reused.
struct RecyclableBatch {
    records: Vec<RikerRecord>,
    len: usize,
    return_tx: Sender<Vec<RikerRecord>>,
}

impl RecyclableBatch {
    /// Valid records in the batch, as a slice.
    fn records(&self) -> &[RikerRecord] {
        &self.records[..self.len]
    }
}

impl Drop for RecyclableBatch {
    fn drop(&mut self) {
        // Hand the inner Vec back to the reader's pool for reuse. The
        // receiver is dropped during shutdown; ignore send errors then.
        let records = std::mem::take(&mut self.records);
        let _ = self.return_tx.send(records);
    }
}

/// Run all collectors sequentially on a single thread (no threading
/// overhead). Drives the reader the same way the parallel reader thread
/// does — in-place fills via [`AlignmentReader::fill_record`] — but
/// loops over every collector for each record instead of fanning out
/// via channels.
fn run_single_threaded(
    mut reader: AlignmentReader,
    mut collectors: Vec<Box<dyn Collector>>,
) -> Result<()> {
    let header = reader.header().clone();
    for collector in &mut collectors {
        collector.initialize(&header)?;
    }

    let requirements = combined_requirements(&collectors);
    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);

    let read_result = run_single_threaded_loop(
        &mut reader,
        &mut collectors,
        &requirements,
        &header,
        &mut progress,
    );
    progress.finish();
    read_result?;

    for collector in &mut collectors {
        collector.finish()?;
    }

    Ok(())
}

fn run_single_threaded_loop(
    reader: &mut AlignmentReader,
    collectors: &mut [Box<dyn Collector>],
    requirements: &RikerRecordRequirements,
    header: &Header,
    progress: &mut ProgressLogger,
) -> Result<()> {
    let mut record = reader.empty_record();
    while reader.fill_record(requirements, &mut record)? {
        progress.record_with(&record, header);
        for collector in collectors.iter_mut() {
            collector.accept(&record, header)?;
        }
    }
    Ok(())
}

/// Compute the union of every collector's [`RikerRecordRequirements`]. The
/// reader uses this to decide which decoder steps to run on every record;
/// each collector then sees the same fully-populated record.
fn combined_requirements(collectors: &[Box<dyn Collector>]) -> RikerRecordRequirements {
    collectors.iter().fold(RikerRecordRequirements::NONE, |acc, c| acc.union(c.field_needs()))
}

/// Run collectors in parallel with a dedicated reader thread and a pool of
/// worker groups.
///
/// Architecture:
/// - One **reader thread** pulls empty [`RecyclableBatch`] slots from a pool,
///   fills each batch via [`AlignmentReader::fill_record`], wraps it in an
///   `Arc`, and sends a clone to **every worker group's** ordered channel (one
///   bounded kanal channel per group).
/// - Each of the `workers` **groups** is owned by a single thread that blocks
///   on its channel's `recv()` and calls `accept_multiple` on the collectors
///   it owns. A collector lives in exactly one group, so there is no shared
///   collector state and no per-collector lock; and because every group
///   receives batches in the same order, each collector sees records in file
///   order (deterministic output).
/// - When the last `Arc<RecyclableBatch>` drops (every group has finished the
///   batch), its `Drop` returns the inner `Vec<RikerRecord>` to the reader's
///   pool via a kanal channel.
///
/// `workers` is the number of groups to spawn (the compute-worker budget),
/// alongside the (uncounted) reader thread; collectors are partitioned into
/// `min(workers, n_collectors)` groups by [`partition_collectors`]. Caller
/// ensures `workers >= 1` — a zero-worker (serial-main) plan is routed to
/// [`run_single_threaded`] instead.
///
/// Backpressure comes from the recycling pool: the reader blocks on
/// `pool_rx.recv()` once `NUM_BATCHES_POOLED` batches are in flight. Each group
/// channel is bounded one batch above that, so it never blocks the reader in
/// practice. No busy-polling. Finalization ([`Collector::finish`]) runs on the
/// main thread after all groups and the reader join cleanly (see the module
/// docs).
fn run_parallel(
    mut reader: AlignmentReader,
    mut collectors: Vec<Box<dyn Collector>>,
    workers: usize,
) -> Result<()> {
    debug_assert!(workers >= 1, "run_parallel requires at least 1 worker group");
    // Clone the header up front so worker threads (which borrow it) and
    // the reader thread (which owns the AlignmentReader) don't fight over
    // it. Header.clone() is one shot at startup.
    let header = reader.header().clone();
    for collector in &mut collectors {
        collector.initialize(&header)?;
    }

    // Each collector declares which expensive fields it reads; the union
    // tells the reader which decoder steps to run per record.
    let requirements = combined_requirements(&collectors);

    // Assign collectors to worker groups, balancing declared cost so the
    // heaviest collectors get their own worker. Each group is owned
    // exclusively by one worker thread — no per-collector mutex, so no
    // lock contention — and the reader ships every batch to every group's
    // ordered channel, so each collector still sees records in file order
    // (deterministic output).
    let groups = partition_collectors(collectors, workers);
    let n_groups = groups.len();

    // Pool of reusable record-batch allocations (the backpressure). Slots
    // are pre-allocated as the variant the reader writes to. Unbounded so
    // `RecyclableBatch::drop` never blocks.
    let (pool_tx, pool_rx) = kanal::unbounded::<Vec<RikerRecord>>();
    for _ in 0..NUM_BATCHES_POOLED {
        let mut vec: Vec<RikerRecord> = Vec::with_capacity(BATCH_SIZE);
        vec.resize_with(BATCH_SIZE, || reader.empty_record());
        pool_tx.send(vec).expect("pool send cannot fail: channel is unbounded and rx is alive");
    }

    // One ordered channel per worker group. Bounded so a slow group applies
    // backpressure; sized one batch above the pool's in-flight max so the
    // recycling pool stays the practical limiter.
    let mut group_senders: Vec<Sender<Batch>> = Vec::with_capacity(n_groups);
    let mut group_receivers: Vec<Receiver<Batch>> = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let (tx, rx) = kanal::bounded::<Batch>(NUM_BATCHES_POOLED + 1);
        group_senders.push(tx);
        group_receivers.push(rx);
    }

    // Poison: set when any thread errors. Signals the reader to stop early.
    let poison = AtomicBool::new(false);
    let header_ref = &header;

    std::thread::scope(|scope| -> Result<()> {
        let poison_ref: &AtomicBool = &poison;

        // One worker per group; each owns its collectors and its receiver.
        // Workers do not finalize; on clean EOF they hand their collectors
        // back for main-thread finalization (see `group_worker_loop`).
        let mut worker_handles = Vec::with_capacity(n_groups);
        for (group, rx) in groups.into_iter().zip(group_receivers) {
            worker_handles
                .push(scope.spawn(move || group_worker_loop(rx, group, header_ref, poison_ref)));
        }

        // Reader owns all group senders + the batch pool. Dropping the
        // senders on return closes the channels so workers drain and exit.
        let requirements_ref = &requirements;
        let reader_handle = scope.spawn(move || {
            let reader_result = reader_thread_loop(
                &mut reader,
                header_ref,
                &group_senders,
                pool_tx,
                pool_rx,
                requirements_ref,
                poison_ref,
            );
            // On a reader error, poison before dropping the senders so any
            // in-flight worker stops at its next batch instead of processing
            // the post-error remainder of a truncated stream. (No-partial-
            // output is guaranteed regardless by main-thread finalization,
            // which never runs when this `reader_result` is an error.)
            if reader_result.is_err() {
                poison_ref.store(true, Ordering::Relaxed);
            }
            drop(group_senders);
            reader_result
        });

        let reader_result = reader_handle.join().map_err(|_| anyhow!("reader thread panicked"))?;

        // Join every worker, collecting the collectors each hands back
        // (unfinalized) and the first error, if any.
        let mut first_error: Option<anyhow::Error> = None;
        let mut finished_groups: Vec<Vec<Box<dyn Collector>>> = Vec::with_capacity(n_groups);
        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(group)) => finished_groups.push(group),
                Ok(Err(e)) => {
                    if first_error.is_none() {
                        first_error = Some(e);
                    }
                }
                Err(_) => {
                    if first_error.is_none() {
                        first_error = Some(anyhow!("worker thread panicked"));
                    }
                }
            }
        }

        // All-or-nothing finalization: only run `finish()` — which writes the
        // output files — once the reader AND every worker have completed
        // cleanly. Checking both error sources before finalizing anything
        // guarantees a failed run leaves no partial output on disk, and the
        // worker joins above establish the happens-before edge that orders
        // every `accept_multiple` before these `finish()` calls.
        reader_result?;
        if let Some(e) = first_error {
            return Err(e);
        }
        for group in &mut finished_groups {
            for collector in group {
                collector.finish()?;
            }
        }
        Ok(())
    })?;

    Ok(())
}

/// Partition collectors into at most `max_groups` worker groups, balancing
/// each collector's declared [`Collector::cost_hint`] across groups so the
/// heaviest collectors tend to land alone. Greedy longest-processing-time:
/// assign each collector (heaviest first) to the currently-lightest group.
/// Returns between 1 and `min(max_groups, n_collectors)` non-empty groups.
fn partition_collectors(
    collectors: Vec<Box<dyn Collector>>,
    max_groups: usize,
) -> Vec<Vec<Box<dyn Collector>>> {
    let n_groups = max_groups.min(collectors.len()).max(1);
    let mut indexed: Vec<(u32, Box<dyn Collector>)> =
        collectors.into_iter().map(|c| (c.cost_hint(), c)).collect();
    // Heaviest first so greedy LPT balances well.
    indexed.sort_by_key(|c| std::cmp::Reverse(c.0));

    let mut groups: Vec<Vec<Box<dyn Collector>>> = (0..n_groups).map(|_| Vec::new()).collect();
    let mut loads: Vec<u64> = vec![0; n_groups];
    for (cost, collector) in indexed {
        let lightest =
            loads.iter().enumerate().min_by_key(|(_, load)| **load).map_or(0, |(idx, _)| idx);
        loads[lightest] += u64::from(cost);
        groups[lightest].push(collector);
    }
    // Defensive only: with >= 1 collector, greedy LPT seeds each of the
    // `n_groups <= n_collectors` groups exactly once, so none end up empty. The
    // sole way to get an empty group is an empty `collectors` (n_groups == 1),
    // which callers never pass.
    groups.retain(|g| !g.is_empty());
    groups
}

/// Reader thread: pulls records into batches and ships each batch to every
/// worker group's ordered channel. Pulls a pre-allocated
/// `Vec<RikerRecord>` from `pool_rx`, fills each slot via
/// [`AlignmentReader::fill_record`], wraps in an `Arc<RecyclableBatch>`,
/// and sends a clone to each group; on drop the batch returns its `Vec` to
/// the pool. Sending to the groups in a fixed order per batch means each
/// group sees batches in file order, so every collector's output is
/// deterministic.
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "pool_tx and pool_rx move into this (scoped) thread: the reader \
              is the only thread that receives from the pool, and we want the \
              reader's handle on pool_tx to drop when the reader exits so \
              in-flight RecyclableBatch Drops (which each carry a clone of \
              pool_tx) become the last senders and the pool channel can \
              close naturally on shutdown"
)]
fn reader_thread_loop(
    reader: &mut AlignmentReader,
    header: &Header,
    group_senders: &[Sender<Batch>],
    pool_tx: Sender<Vec<RikerRecord>>,
    pool_rx: Receiver<Vec<RikerRecord>>,
    requirements: &RikerRecordRequirements,
    poison: &AtomicBool,
) -> Result<()> {
    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);

    let result = 'outer: loop {
        if poison.load(Ordering::Relaxed) {
            break Ok(());
        }

        // Blocks when every batch is in flight — natural backpressure.
        let Ok(mut records) = pool_rx.recv() else {
            break Ok(());
        };

        let mut len = 0;
        while len < records.len() {
            match reader.fill_record(requirements, &mut records[len]) {
                Ok(false) => break,
                Ok(true) => {
                    progress.record_with(&records[len], header);
                    len += 1;
                }
                Err(e) => break 'outer Err(e),
            }
        }

        if len == 0 {
            let _ = pool_tx.send(records);
            break Ok(());
        }

        let batch: Batch = Arc::new(RecyclableBatch { records, len, return_tx: pool_tx.clone() });
        for tx in group_senders {
            if poison.load(Ordering::Relaxed) {
                break 'outer Ok(());
            }
            // A group worker gone (errored, dropped its receiver) closes the
            // channel; stop the whole run.
            if tx.send(Arc::clone(&batch)).is_err() {
                break 'outer Ok(());
            }
        }
        drop(batch);
    };

    progress.finish();
    result
}

/// Group worker: owns a set of collectors and their shared ordered channel.
/// For each batch it runs every owned collector in turn, then drops the
/// batch (recycled once every group is done with it). It does **not**
/// finalize: on clean EOF it hands its collectors back to the main thread
/// (which finalizes every group only once all workers and the reader have
/// joined cleanly — see [`run_parallel`]), so a late error in a sibling
/// group can never leave this group's output on disk. On a collector error
/// it sets the poison flag and returns the error without finalizing.
#[allow(
    clippy::needless_pass_by_value,
    reason = "each worker owns its receiver + collector group; passing by \
              reference would keep the channel open and prevent shutdown"
)]
fn group_worker_loop(
    rx: Receiver<Batch>,
    mut group: Vec<Box<dyn Collector>>,
    header: &Header,
    poison: &AtomicBool,
) -> Result<Vec<Box<dyn Collector>>> {
    while let Ok(batch) = rx.recv() {
        if poison.load(Ordering::Relaxed) {
            return Ok(group);
        }
        for collector in &mut group {
            if let Err(e) = collector.accept_multiple(batch.records(), header) {
                poison.store(true, Ordering::Relaxed);
                return Err(e);
            }
        }
    }
    // Channel closed = the reader finished. Hand the collectors back
    // unfinalized; the main thread finalizes all groups only if the whole
    // run succeeded.
    Ok(group)
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::riker_record::RikerRecord;

    /// Collector that fails on the Nth `accept` call to exercise the
    /// poison-flag shutdown path of `run_parallel`.
    struct FailingCollector {
        seen: u64,
        fail_after: u64,
    }

    impl Collector for FailingCollector {
        fn initialize(&mut self, _h: &Header) -> Result<()> {
            Ok(())
        }
        fn accept(&mut self, _r: &RikerRecord, _h: &Header) -> Result<()> {
            self.seen += 1;
            if self.seen >= self.fail_after {
                return Err(anyhow!("synthetic failure after {} records", self.seen));
            }
            Ok(())
        }
        fn finish(&mut self) -> Result<()> {
            Ok(())
        }
        fn name(&self) -> &'static str {
            "failing"
        }
        fn field_needs(&self) -> RikerRecordRequirements {
            RikerRecordRequirements::NONE
        }
    }

    /// Collector that records whether `finish` was called, used to assert that
    /// a poisoned run finalizes no one.
    struct SpyCollector {
        finished: Arc<AtomicBool>,
    }

    impl Collector for SpyCollector {
        fn initialize(&mut self, _h: &Header) -> Result<()> {
            Ok(())
        }
        fn accept(&mut self, _r: &RikerRecord, _h: &Header) -> Result<()> {
            Ok(())
        }
        fn finish(&mut self) -> Result<()> {
            self.finished.store(true, Ordering::Relaxed);
            Ok(())
        }
        fn name(&self) -> &'static str {
            "spy"
        }
        fn field_needs(&self) -> RikerRecordRequirements {
            RikerRecordRequirements::NONE
        }
    }

    /// Write a tiny temp BAM with `n` mapped records on chr1 for the
    /// `run_parallel` tests. `helpers/` is only on the integration-test target,
    /// so we build the BAM inline here.
    fn tiny_bam(n: u32) -> Result<tempfile::NamedTempFile> {
        use noodles::bam;
        use noodles::sam::Header;
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

        let header = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(10_000).unwrap()),
            )
            .build();
        let tmp = tempfile::NamedTempFile::with_suffix(".bam")?;
        let file = std::fs::File::create(tmp.path())?;
        let mut writer = bam::io::Writer::new(std::io::BufWriter::new(file));
        writer.write_header(&header)?;
        let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        for i in 0..n {
            let pos =
                noodles::core::Position::new(usize::try_from(i).unwrap() % 9_000 + 1).unwrap();
            let record = RecordBuf::builder()
                .set_name(format!("r{i}").into_bytes())
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(pos)
                .set_cigar(cigar.clone())
                .set_sequence(Sequence::from(vec![b'A'; 50]))
                .set_quality_scores(QualityScores::from(vec![30u8; 50]))
                .build();
            writer.write_alignment_record(&header, &record)?;
        }
        drop(writer);
        Ok(tmp)
    }

    /// A failing collector inside `run_parallel` should propagate its
    /// error out and not deadlock or panic. Asserts on both the error
    /// type and that `run_parallel` actually returns (rather than
    /// hanging forever waiting for the reader).
    #[test]
    fn run_parallel_propagates_collector_error() -> Result<()> {
        let tmp = tiny_bam(2_000)?;
        let reader = AlignmentReader::open(tmp.path(), None, 0)?;
        let collectors: Vec<Box<dyn Collector>> =
            vec![Box::new(FailingCollector { seen: 0, fail_after: 100 })];

        // 2 worker threads + the reader thread.
        let result = run_parallel(reader, collectors, 2);
        let err = result.expect_err("run_parallel should propagate the collector error");
        assert!(
            err.to_string().contains("synthetic failure"),
            "expected the failing collector's error, got: {err}"
        );
        Ok(())
    }

    /// On a collector error, `run_parallel` must not finalize any *surviving*
    /// worker group: finalization (which writes output) runs on the main
    /// thread only after every worker joins cleanly. With two collectors in
    /// two groups, a failure in one group must leave the other group's
    /// `finish()` uncalled — regardless of how fast the surviving group drains.
    #[test]
    fn run_parallel_does_not_finalize_survivors_on_error() -> Result<()> {
        let tmp = tiny_bam(2_000)?;
        let reader = AlignmentReader::open(tmp.path(), None, 0)?;
        let finished = Arc::new(AtomicBool::new(false));
        let collectors: Vec<Box<dyn Collector>> = vec![
            Box::new(FailingCollector { seen: 0, fail_after: 100 }),
            Box::new(SpyCollector { finished: Arc::clone(&finished) }),
        ];

        // 2 workers => partition_collectors puts each collector in its own group.
        let result = run_parallel(reader, collectors, 2);
        assert!(result.is_err(), "the failing collector must poison the run");
        assert!(
            !finished.load(Ordering::Relaxed),
            "the surviving group's collector must not be finalized on a failed run"
        );
        Ok(())
    }

    fn nz(n: usize) -> NonZero<usize> {
        NonZero::new(n).expect("test uses non-zero")
    }

    fn plan(t: usize, n_tools: usize, kind: InputKind) -> (usize, usize) {
        let p = plan_multi(nz(t), n_tools, kind);
        (p.decode_threads, p.compute_workers)
    }

    #[test]
    fn plan_bam_is_worker_leaning() {
        // t=1 serial; t=2,3 dispatcher + workers (inflate is cheap, done inline);
        // t>=4 workers get the odd thread (workers=ceil(t/2), readers=rest).
        let cases = [
            (1, (0, 0)),
            (2, (0, 1)),
            (3, (0, 2)),
            (4, (2, 2)),
            (5, (2, 3)),
            (6, (3, 3)),
            (7, (3, 4)),
        ];
        for (t, want) in cases {
            assert_eq!(plan(t, 5, InputKind::Bam), want, "bam t={t}");
        }
    }

    #[test]
    fn plan_cram_is_read_leaning() {
        // t=1 serial; t=2,3 serial-main + decode threads (a scarce thread should
        // decode, not idle as a worker); t>=4 readers get the odd thread.
        let cases = [
            (1, (0, 0)),
            (2, (1, 0)),
            (3, (2, 0)),
            (4, (2, 2)),
            (5, (3, 2)),
            (6, (3, 3)),
            (7, (4, 3)),
        ];
        for (t, want) in cases {
            assert_eq!(plan(t, 5, InputKind::Cram), want, "cram t={t}");
        }
    }

    #[test]
    fn plan_sam_gives_everything_to_workers() {
        // SAM parse is single-threaded; no decode threads, all budget to workers.
        for (t, want) in [(1, (0, 0)), (2, (0, 1)), (3, (0, 2)), (4, (0, 4)), (8, (0, 5))] {
            assert_eq!(plan(t, 5, InputKind::Sam), want, "sam t={t}");
        }
    }

    #[test]
    fn plan_caps_workers_at_tool_count_and_spends_the_rest_on_reads() {
        // 2 tools -> at most 2 workers; freed budget becomes read threads.
        assert_eq!(plan(8, 2, InputKind::Bam), (6, 2), "bam capped");
        assert_eq!(plan(8, 2, InputKind::Cram), (6, 2), "cram capped");
        // A single tool never spins more than one worker.
        assert_eq!(plan(6, 1, InputKind::Bam), (5, 1), "bam single tool");
    }

    #[test]
    fn plan_bam_and_cram_spend_the_whole_budget() {
        // BAM/CRAM absorb any worker-cap slack into read threads, so readers +
        // workers always equal the budget.
        for kind in [InputKind::Bam, InputKind::Cram] {
            for t in 4..=64 {
                let (r, w) = plan(t, 5, kind);
                assert_eq!(r + w, t, "kind={kind:?} t={t} should spend the whole budget");
            }
        }
    }

    #[test]
    fn plan_sam_is_bounded_by_tool_count() {
        // SAM can't decode in parallel, so it can't use more threads than there
        // are collector workers — extra budget is legitimately unused.
        for t in 4..=64 {
            let (r, w) = plan(t, 5, InputKind::Sam);
            assert_eq!((r, w), (0, t.min(5)), "sam t={t}");
        }
    }

    #[test]
    fn plan_input_kind_classifies_by_extension() {
        assert_eq!(InputKind::of(Path::new("x.cram")), InputKind::Cram);
        assert_eq!(InputKind::of(Path::new("x.CRAM")), InputKind::Cram);
        assert_eq!(InputKind::of(Path::new("x.sam")), InputKind::Sam);
        assert_eq!(InputKind::of(Path::new("x.sam.gz")), InputKind::Sam);
        assert_eq!(InputKind::of(Path::new("x.bam")), InputKind::Bam);
    }
}
