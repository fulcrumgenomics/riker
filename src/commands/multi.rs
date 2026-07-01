//! The `multi` command: one pass over the input, many collectors, run in parallel.
//!
//! ## Threading model
//!
//! multi receives the toolkit-wide `--threads` budget and divides it in
//! [`Multi::plan_threads`] into reader-decode workers and compute **worker
//! groups**. One thread drives the reader/fan-out; the remaining budget is
//! split by input format — BAM decode is cheap and multi is compute-bound on
//! it, so nearly all of it becomes worker groups, while CRAM decode is
//! expensive and gets its own pool. Collectors are partitioned into
//! `min(workers, n_collectors)` groups by [`partition_collectors`], each
//! owned by one worker thread. A budget of 1 skips the parallel pipeline
//! entirely (see "Single-threaded path" below).
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
//! - `batch_pool` — an unbounded mpsc channel of empty `Vec<RikerRecord>`
//!   slots the reader pulls from (and which [`RecyclableBatch::drop`]
//!   returns to). Slots are pre-allocated as the record variant the reader
//!   writes into: `RikerRecord::Bam` for BAM, `RikerRecord::Fallback` for
//!   SAM, `RikerRecord::Htslib` for CRAM.
//! - one **ordered channel per worker group** — a bounded crossbeam queue
//!   of `Arc<RecyclableBatch>`. The reader wraps each filled batch in an
//!   `Arc` and sends a clone to every group's channel; the owning worker
//!   blocks on `recv()` — no polling, no condvar. Bounded one batch above
//!   the pool's in-flight max so the recycling pool is the practical
//!   backpressure.
//! - `return` — the mpsc `return_tx`/`return_rx` captured inside each
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
//! closes, and each worker finalizes its own collectors
//! ([`Collector::finish`]) before exiting.
//!
//! ## Error propagation
//!
//! Channel disconnection handles the happy shutdown path. An `AtomicBool`
//! poison flag handles the sad one.
//!
//! - A worker that gets `Err` from `accept_multiple` sets the flag before
//!   returning (without finalizing), so the reader aborts on its next send
//!   and sibling workers stop at their next batch instead of processing the
//!   rest of the in-flight work.
//! - If the reader itself errors it returns `Err`; `run_parallel` sets the
//!   poison flag after `reader_handle.join()` so workers skip finalization
//!   and return without writing partial output.
//!
//! Errors are surfaced by `handle.join()` on the main thread — the
//! reader's error wins if present, otherwise the first worker error.
//!
//! All `poison` flag accesses use `Ordering::Relaxed`. Relaxed is
//! sufficient because the flag is a hint that *some* thread errored; the
//! actual shared state (the channels) is synchronised by its own
//! primitives and provides whatever happens-before edges the program
//! needs.
//!
//! ## Single-threaded path
//!
//! For a budget of 1 the whole thing collapses to [`run_single_threaded`]:
//! no channels, no extra threads, just a serial loop that drives the
//! same reader API — in-place fills via [`AlignmentReader::fill_record`]
//! for every format. Useful for testing and small runs where threading
//! overhead isn't worth it.

use std::fmt;
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use anyhow::{Result, anyhow, bail};
use clap::{Args, ValueEnum};
use crossbeam_channel::{Receiver, Sender};
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
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::record_utils::derive_sample;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};

/// Number of records per batch sent through the work queue.
///
/// Picked at 128 from a `batch ∈ {32,48,64,128} × pool ∈ {16,32,48,64} ×
/// threads ∈ {2,4,8}` sweep on M1 Max and Graviton 3. Smaller batches
/// (e.g. 32) win marginally at `--threads 2` but cause heavy sys-time
/// blowup at higher thread counts: each batch acquires/releases the
/// per-collector mutex, so small batches make many workers fight for the
/// lock. 128 keeps each `accept_multiple` call long enough that mutex
/// churn stays cheap, while the working set (`128 × 16 = 2048` records,
/// ~1.5 MB for typical 150 bp WGS reads) still fits comfortably in L2 on
/// M1 Max and shared L3 on Graviton 3.
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

    /// DEPRECATED: use the top-level `riker --threads` instead.
    ///
    /// Sets the total thread budget for this run (which multi divides between
    /// input decoding and compute workers). Kept for backwards compatibility;
    /// it will be removed in a future release. Setting both this and
    /// `riker --threads` is an error.
    #[arg(long, value_parser = clap::value_parser!(u8).range(1..), help_heading = "Multi Command Options")]
    pub threads: Option<u8>,

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
    /// Resolve the total thread budget from the top-level `--threads`
    /// (`global`) and the deprecated `multi --threads` (`self.threads`).
    /// Setting both is an error; the deprecated flag warns; neither falls
    /// back to [`default_threads`](Command::default_threads).
    fn resolve_total_threads(&self, global: Option<u8>) -> Result<NonZero<usize>> {
        match (global, self.threads) {
            (Some(_), Some(_)) => bail!(
                "set the thread count with the top-level `riker --threads`, not both it and \
                 `multi --threads` (the latter is deprecated)"
            ),
            (Some(_), None) => Ok(resolve_threads(global, self.default_threads())),
            (None, Some(_)) => {
                log::warn!(
                    "`multi --threads` is deprecated and will be removed in a future release; \
                     use the top-level `riker --threads` instead"
                );
                Ok(resolve_threads(self.threads, self.default_threads()))
            }
            (None, None) => Ok(self.default_threads()),
        }
    }

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
        let total = self.resolve_total_threads(threads)?;

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

        // Divide the budget between input decode and compute workers, then
        // open the reader (also needed to build interval maps for WGS).
        let plan = self.plan_threads(total, &self.input.input);
        let reader = AlignmentReader::open(
            &self.input.input,
            self.reference.reference.as_deref(),
            plan.decode_threads,
        )?;

        let collectors = self.build_collectors(&seen, reader.header())?;

        if plan.compute_workers > 1 {
            run_parallel(reader, collectors, plan.compute_workers)?;
        } else {
            run_single_threaded(reader, collectors)?;
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

    /// Split the budget between reader-decode workers and compute worker
    /// groups. multi's orchestrator thread is near-idle, so the whole budget
    /// goes to the reader + workers rather than reserving one for the main
    /// thread. BAM decode is cheap and multi is compute-bound on it, so favor
    /// compute workers; CRAM decode is expensive and scales, so give the
    /// reader a decode pool. (Starting heuristic — constants to be tuned.)
    fn plan_threads(&self, total: NonZero<usize>, input: &Path) -> ThreadPlan {
        let n = total.get();
        if n <= 1 {
            return ThreadPlan { decode_threads: 0, compute_workers: 1 };
        }
        // One thread drives the reader/fan-out; the rest are split below.
        let rest = n - 1;
        let is_cram = input.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let decode_threads = if is_cram { (rest / 2).min(4) } else { 0 };
        let compute_workers = rest.saturating_sub(decode_threads).max(1);
        ThreadPlan { decode_threads, compute_workers }
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
    return_tx: mpsc::Sender<Vec<RikerRecord>>,
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

/// Run collectors in parallel with a dedicated reader thread and a shared
/// worker pool.
///
/// Architecture:
/// - One **reader thread** pulls empty [`RecyclableBatch`] slots from a
///   pool, fills each batch via [`AlignmentReader::fill_record`], wraps
///   it in an `Arc`, and fans it out onto a shared
///   `(collector_idx, Batch)` work queue — one entry per collector.
/// - `threads` **pool threads** block on the shared work queue. On receipt of
///   a work item they lock that collector's mutex and call `accept_multiple`.
///   The mutex serializes per-collector accepts (required since `Collector`
///   is stateful) while still allowing different collectors' batches to
///   process in parallel across threads.
/// - When the last `Arc<RecyclableBatch>` drops, its `Drop` returns the inner
///   `Vec<RikerRecord>` to the reader's pool via an mpsc channel.
///
/// `threads` is the worker count, so the pool spawns `threads` workers
/// alongside the (uncounted) reader thread. Caller is responsible for
/// ensuring `threads >= 1` and that the `--threads 1` serial case is
/// routed to [`run_single_threaded`] instead — this path is used for
/// `--threads >= 2`.
///
/// Workers block on the MPMC work queue; the reader fans batches through
/// it. Backpressure comes from the recycling pool (the reader blocks on
/// `pool_rx.recv()` once `NUM_BATCHES_POOLED` batches are in flight).
/// The work queue is bounded too, but at one batch worth above the
/// pool's natural in-flight max, so it never blocks the reader in
/// practice. No busy-polling, no `Condvar`.
fn run_parallel(
    mut reader: AlignmentReader,
    mut collectors: Vec<Box<dyn Collector>>,
    threads: usize,
) -> Result<()> {
    debug_assert!(threads >= 1, "run_parallel requires at least 1 worker thread");
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
    let groups = partition_collectors(collectors, threads);
    let n_groups = groups.len();

    // Pool of reusable record-batch allocations (the backpressure). Slots
    // are pre-allocated as the variant the reader writes to. Unbounded so
    // `RecyclableBatch::drop` never blocks.
    let (pool_tx, pool_rx) = mpsc::channel::<Vec<RikerRecord>>();
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
        let (tx, rx) = crossbeam_channel::bounded::<Batch>(NUM_BATCHES_POOLED + 1);
        group_senders.push(tx);
        group_receivers.push(rx);
    }

    // Poison: set when any thread errors. Signals the reader to stop early.
    let poison = AtomicBool::new(false);
    let header_ref = &header;

    std::thread::scope(|scope| -> Result<()> {
        let poison_ref: &AtomicBool = &poison;

        // One worker per group; each owns its collectors and its receiver.
        // Workers finalize their own collectors on clean EOF (see
        // `group_worker_loop`).
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
            drop(group_senders);
            reader_result
        });

        let reader_result = reader_handle.join().map_err(|_| anyhow!("reader thread panicked"))?;
        if let Err(e) = reader_result {
            poison.store(true, Ordering::Relaxed);
            for handle in worker_handles {
                let _ = handle.join();
            }
            return Err(e);
        }

        let mut first_error: Option<anyhow::Error> = None;
        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
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

        if let Some(e) = first_error {
            return Err(e);
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
    // Drop any empty groups (only possible if n_collectors < n_groups, which
    // `min` already prevents, but keep it robust).
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
    pool_tx: mpsc::Sender<Vec<RikerRecord>>,
    pool_rx: mpsc::Receiver<Vec<RikerRecord>>,
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
/// batch (recycled once every group is done with it). On clean EOF (the
/// reader closed the channel) it finalizes each owned collector — unless the
/// run was poisoned by an error elsewhere, in which case it returns without
/// finalizing so partial output isn't written.
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
) -> Result<()> {
    while let Ok(batch) = rx.recv() {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }
        for collector in &mut group {
            if let Err(e) = collector.accept_multiple(batch.records(), header) {
                poison.store(true, Ordering::Relaxed);
                return Err(e);
            }
        }
    }
    // Channel closed = the reader finished. Finalize only if nothing errored.
    if poison.load(Ordering::Relaxed) {
        return Ok(());
    }
    for collector in &mut group {
        collector.finish()?;
    }
    Ok(())
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

    /// A failing collector inside `run_parallel` should propagate its
    /// error out and not deadlock or panic. Asserts on both the error
    /// type and that `run_parallel` actually returns (rather than
    /// hanging forever waiting for the reader).
    #[test]
    fn run_parallel_propagates_collector_error() -> Result<()> {
        use std::path::Path;
        // Build a small in-memory BAM via the helpers crate. We rebuild
        // the helper inline because helpers/ is only on the integration
        // test target. A few hundred records is plenty to ensure the
        // failing collector trips before EOF.
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
        {
            let file = std::fs::File::create(tmp.path())?;
            let mut writer = bam::io::Writer::new(std::io::BufWriter::new(file));
            writer.write_header(&header)?;
            let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
            for i in 0u32..2_000 {
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
        }

        let reader = AlignmentReader::open(Path::new(tmp.path()), None, 0)?;
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
}
