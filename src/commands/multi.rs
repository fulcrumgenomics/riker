//! The `multi` command: one pass over the input, many collectors, run in parallel.
//!
//! ## Threading model
//!
//! `--threads N` is the **total** thread count: 1 reader thread plus
//! `N - 1` pool workers. So `--threads 2` means 1 reader + 1 worker,
//! `--threads 4` means 1 reader + 3 workers. `--threads 1` skips the
//! parallel pipeline entirely (see "Single-threaded path" below).
//!
//! Reader and workers are connected through three channels:
//!
//! - `batch_pool` — an unbounded mpsc channel of empty
//!   `Vec<RikerRecord>` slots the reader pulls from (and which
//!   [`RecyclableBatch::drop`] returns to). Slots are pre-allocated as
//!   the record variant the reader writes into: `RikerRecord::Bam` for
//!   BAM, `RikerRecord::Fallback` for SAM. The pool is empty for CRAM
//!   because CRAM batches are one-way (see below).
//! - `work_queue` — a bounded crossbeam MPMC queue of
//!   `(collector_idx, Arc<RecyclableBatch>)`. The reader fans each batch
//!   onto it once per collector; pool workers block on it. The bound is
//!   sized as `(NUM_BATCHES_POOLED + 1) * n_collectors`: one batch worth
//!   above the pool's natural in-flight max, so the BAM/SAM reader never
//!   actually blocks on send — the pool is still the practical
//!   backpressure there. The bound's real job is to backstop the CRAM
//!   path, which has no pool.
//! - `return` — the mpsc `return_tx`/`return_rx` captured inside each
//!   `RecyclableBatch`. When the last `Arc` reference drops,
//!   [`RecyclableBatch::drop`] sends the inner `Vec` back to the pool
//!   (BAM/SAM) or just drops it (CRAM, `return_tx: None`).
//!
//! For BAM and SAM the reader reads records in place into a batch
//! (`AlignmentReader::fill_record` overwrites each slot). For CRAM it
//! drives `AlignmentReader::riker_records` — noodles allocates each
//! record fresh, and the bounded work queue caps how far ahead the
//! reader can run when there's no recyclable pool. Either way the
//! reader wraps the batch in an `Arc`, clones it once per collector
//! onto the work queue, and workers
//! block on `work_rx.recv()` — no polling, no condvar. On receipt they
//! lock the per-collector mutex and call `accept_multiple`, which
//! serialises accesses to any single collector while still letting
//! different collectors process in parallel across workers.
//!
//! The reader also consults the union of every active collector's
//! [`Collector::field_needs`] once up front and passes it down. On BAM
//! (where aux decode is lazy) this gates which decoders run per record,
//! so a collector set that never reads aux tags pays zero aux-decode
//! cost. SAM and CRAM decode eagerly inside noodles, so the union is
//! informational there.
//!
//! When the reader hits EOF it drops `work_tx`, the queue closes, and
//! workers exit once they have drained. `Collector::finish` is called
//! on the main thread after the scope joins — no cross-thread
//! finalisation race.
//!
//! ## Error propagation
//!
//! Channel disconnection handles the happy shutdown path. An `AtomicBool`
//! poison flag handles the sad one.
//!
//! - A pool worker that gets `Err` from `accept_multiple` sets the flag
//!   before returning, so the reader can abort within one dispatch
//!   instead of waiting for every worker `Receiver` to drop. Without this,
//!   a single-worker error with N workers still leaves N-1 workers
//!   consuming, and the reader only stops when `work_tx.send` eventually
//!   fails because every receiver has been dropped.
//! - If the reader itself errors it simply returns `Err`; `run_parallel`
//!   then sets the poison flag (after `reader_handle.join()`) so queued
//!   batches still in the work queue are skipped rather than processed.
//!
//! Errors are surfaced by `handle.join()` on the main thread — the
//! reader's error wins if present, otherwise the first pool-worker error.
//!
//! All `poison` flag accesses use `Ordering::Relaxed`. Relaxed is
//! sufficient because the flag is a hint that *some* thread errored;
//! the actual shared state (per-collector mutexes, the work queue) is
//! synchronised by its own primitives and provides whatever
//! happens-before edges the program needs.
//!
//! ## Single-threaded path
//!
//! For `--threads 1` the whole thing collapses to [`run_single_threaded`]:
//! no channels, no extra threads, just a serial loop that drives the
//! same reader API (in-place fills for BAM/SAM, iterator for CRAM).
//! Useful for testing and small runs where threading overhead isn't
//! worth it.

use std::fmt;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use anyhow::{Result, anyhow};
use clap::{Args, ValueEnum};
use crossbeam_channel::{Receiver, Sender};
use noodles::sam::Header;

use crate::collector::Collector;
use crate::commands::alignment::{AlignmentCollector, MultiAlignmentOptions};
use crate::commands::basic::BasicCollector;
use crate::commands::command::Command;
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

    /// Total number of threads to use. The reader counts as one of them, so
    /// `--threads N` means 1 reader + `N - 1` pool workers. `--threads 1`
    /// disables the parallel pipeline entirely (no reader thread, no
    /// channels) and runs serially on the main thread; `--threads 2` is
    /// 1 reader + 1 pool worker; and so on.
    #[arg(long, default_value_t = 2, help_heading = "Multi Command Options")]
    pub threads: usize,

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
    fn execute(&self) -> Result<()> {
        if self.threads == 0 {
            return Err(anyhow!("--threads must be >= 1"));
        }

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

        // Open the reader so we can build interval maps for WGS if needed.
        let reader = AlignmentReader::open(&self.input.input, self.reference.reference.as_deref())?;

        let collectors = self.build_collectors(&seen, reader.header())?;

        if self.threads > 1 {
            run_parallel(reader, collectors, self.threads)?;
        } else {
            run_single_threaded(reader, collectors)?;
        }

        Ok(())
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

/// A single work item on the MPMC work queue: which collector the batch is
/// destined for, and the shared batch itself.
type WorkItem = (usize, Batch);
type WorkTx = Sender<WorkItem>;
type WorkRx = Receiver<WorkItem>;

/// Owns a `Vec<RikerRecord>` of capacity `BATCH_SIZE` plus a count of
/// valid records. On drop, if `return_tx` is `Some`, the inner `Vec` is
/// returned to the reader's pool so its allocations can be reused.
///
/// One-way batches (CRAM, where slots can't be recycled because each
/// record is a fresh allocation anyway) carry `return_tx: None`; the
/// `Vec` is just dropped on the last `Arc` release.
struct RecyclableBatch {
    records: Vec<RikerRecord>,
    len: usize,
    return_tx: Option<mpsc::Sender<Vec<RikerRecord>>>,
}

impl RecyclableBatch {
    /// Valid records in the batch, as a slice.
    fn records(&self) -> &[RikerRecord] {
        &self.records[..self.len]
    }
}

impl Drop for RecyclableBatch {
    fn drop(&mut self) {
        if let Some(tx) = self.return_tx.take() {
            // Hand the inner Vec back to the reader's pool for reuse. The
            // receiver is dropped during shutdown; ignore send errors then.
            let records = std::mem::take(&mut self.records);
            let _ = tx.send(records);
        }
    }
}

/// Run all collectors sequentially on a single thread (no threading
/// overhead). Drives the reader the same way the parallel reader thread
/// does — in-place fills for BAM/SAM, allocate-per-record iterator for
/// CRAM — but loops over every collector for each record instead of
/// fanning out via channels.
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
    if reader.supports_in_place_reads() {
        let mut record = reader.empty_record();
        while reader.fill_record(requirements, &mut record)? {
            progress.record_with(&record, header);
            for collector in collectors.iter_mut() {
                collector.accept(&record, header)?;
            }
        }
    } else {
        for result in reader.riker_records(requirements) {
            let record = result?;
            progress.record_with(&record, header);
            for collector in collectors.iter_mut() {
                collector.accept(&record, header)?;
            }
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
/// - One **reader thread** pulls empty [`RecyclableBatch`] slots from a pool
///   (BAM/SAM) or allocates them on the fly (CRAM), fills each batch via the
///   reader API, wraps it in an `Arc`, and fans it out onto a shared
///   `(collector_idx, Batch)` work queue — one entry per collector.
/// - `threads` **pool threads** block on the shared work queue. On receipt of
///   a work item they lock that collector's mutex and call `accept_multiple`.
///   The mutex serializes per-collector accepts (required since `Collector`
///   is stateful) while still allowing different collectors' batches to
///   process in parallel across threads.
/// - When the last `Arc<RecyclableBatch>` drops, its `Drop` returns the inner
///   `Vec<RikerRecord>` to the reader's pool via an mpsc channel — except
///   for one-way CRAM batches (`return_tx: None`), which just drop.
///
/// `threads` is the **total** thread count (reader + workers), so the
/// pool spawns `threads - 1` workers. Caller is responsible for ensuring
/// `threads >= 2` (the `--threads 1` case routes to
/// [`run_single_threaded`] instead).
///
/// Workers block on the MPMC work queue; the reader fans batches through
/// it. Backpressure on BAM/SAM comes from the recycling pool (the reader
/// blocks on `pool_rx.recv()` once `NUM_BATCHES_POOLED` batches are in
/// flight). The work queue is bounded too, but at one batch worth above
/// the pool's natural in-flight max, so its bound only kicks in on the
/// CRAM path (which has no pool). No busy-polling, no `Condvar`.
fn run_parallel(
    mut reader: AlignmentReader,
    mut collectors: Vec<Box<dyn Collector>>,
    threads: usize,
) -> Result<()> {
    debug_assert!(threads >= 2, "run_parallel requires at least 2 total threads");
    let pool_workers = threads - 1;
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

    // Each collector lives behind a Mutex so a pool thread can claim exclusive
    // access before calling `accept_multiple`. The mutex only serializes accesses
    // to a single collector; different collectors process in parallel.
    let slots: Vec<Mutex<Box<dyn Collector>>> = collectors.into_iter().map(Mutex::new).collect();

    // Shared MPMC work queue: reader sends `(collector_idx, batch)` items; pool
    // threads pull and dispatch to the matching collector.
    //
    // The bound is set one batch worth above the pool's natural in-flight max:
    // BAM/SAM produce at most `NUM_BATCHES_POOLED * n_collectors` items in the
    // queue (every in-flight batch fanned to every collector), so the +1 batch
    // of slack means the BAM/SAM reader never blocks on send — the pool stays
    // the practical backpressure there. The bound's real job is the CRAM path,
    // which has no pool: without it, a slow collector would let the reader
    // pull arbitrarily many records into memory ahead of workers. `.max(1)`
    // keeps a 0-collector run from creating a 0-capacity channel.
    let work_queue_bound = (NUM_BATCHES_POOLED + 1) * slots.len().max(1);
    let (work_tx, work_rx): (WorkTx, WorkRx) = crossbeam_channel::bounded(work_queue_bound);

    // Pool of reusable record-batch allocations, used for BAM and SAM where
    // we read in place. CRAM cannot recycle slots (each record is freshly
    // allocated by noodles), so its batches are one-way (`return_tx: None`)
    // and we don't pre-fill the pool. Unbounded so `RecyclableBatch::drop`
    // never blocks.
    //
    // Pool slots are pre-allocated as the variant the reader writes to —
    // `RikerRecord::Bam(BamRec::new())` for BAM (skips the eager aux-tag
    // decode) and `RikerRecord::Fallback(FallbackRec::empty())` for SAM.
    let (pool_tx, pool_rx) = mpsc::channel::<Vec<RikerRecord>>();
    if reader.supports_in_place_reads() {
        for _ in 0..NUM_BATCHES_POOLED {
            let mut vec: Vec<RikerRecord> = Vec::with_capacity(BATCH_SIZE);
            vec.resize_with(BATCH_SIZE, || reader.empty_record());
            // Channel is unbounded (mpsc) and `pool_rx` is held in this
            // function until the reader thread takes ownership later, so
            // the send cannot fail here.
            pool_tx.send(vec).expect("pool send cannot fail: channel is unbounded and rx is alive");
        }
    }

    // Poison: set when any thread errors. Signals the reader to stop early.
    let poison = AtomicBool::new(false);
    let n_collectors = slots.len();

    // Run the reader + worker threads inside a scope so they can borrow
    // `slots` and `poison`. The scope ends (and borrows release) before we
    // reclaim `slots` to call `finish()`.
    let header_ref = &header;
    std::thread::scope(|scope| -> Result<()> {
        let slots_ref: &[Mutex<Box<dyn Collector>>] = &slots;
        let poison_ref: &AtomicBool = &poison;

        // Pool workers share `work_rx` (MPMC) and the slots. We spawn
        // `threads - 1` of them; the reader counts as the Nth thread.
        let mut pool_handles = Vec::with_capacity(pool_workers);
        for _ in 0..pool_workers {
            let work_rx = work_rx.clone();
            pool_handles.push(
                scope.spawn(move || pool_worker_loop(work_rx, slots_ref, header_ref, poison_ref)),
            );
        }
        // Drop our outer receiver handle so the queue closes once the reader's
        // sender drops (otherwise pool threads would wait forever).
        drop(work_rx);

        // Reader owns the work-queue sender and the batch pool.
        let requirements_ref = &requirements;
        let reader_handle = scope.spawn(move || {
            let reader_result = reader_thread_loop(
                &mut reader,
                header_ref,
                &work_tx,
                n_collectors,
                pool_tx,
                pool_rx,
                requirements_ref,
                poison_ref,
            );
            drop(work_tx);
            reader_result
        });

        let reader_result = reader_handle.join().map_err(|_| anyhow!("reader thread panicked"))?;
        if let Err(e) = reader_result {
            poison.store(true, Ordering::Relaxed);
            for handle in pool_handles {
                let _ = handle.join();
            }
            return Err(e);
        }

        let mut first_error: Option<anyhow::Error> = None;
        for handle in pool_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => {
                    if first_error.is_none() {
                        first_error = Some(e);
                    }
                }
                Err(_) => {
                    if first_error.is_none() {
                        first_error = Some(anyhow!("pool thread panicked"));
                    }
                }
            }
        }

        if let Some(e) = first_error {
            return Err(e);
        }
        Ok(())
    })?;

    // Scope has ended; `slots` is owned again and no thread touches it.
    // Finalize each collector exactly once.
    for slot in slots {
        let mut collector = slot.into_inner().unwrap();
        collector.finish()?;
    }
    Ok(())
}

/// Reader thread: pulls records into batches and fans each batch onto the
/// shared work queue once per collector.
///
/// Two paths inside one function. BAM/SAM use the in-place pool: the reader
/// pulls a pre-allocated `Vec<RikerRecord>` from `pool_rx` and fills each
/// slot via [`AlignmentReader::fill_record`]; on drop the batch returns its
/// `Vec` to the pool. CRAM uses [`AlignmentReader::riker_records`] —
/// noodles allocates each `RikerRecord::Fallback` afresh — and the batch
/// is one-way (`return_tx: None`), so the pool channel is unused.
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "pool_tx and pool_rx move into this (scoped) thread: the reader \
              is the only thread that receives from the pool, and we want the \
              reader's handle on pool_tx to drop when the reader exits so \
              in-flight RecyclableBatch Drops (which each carry a clone of \
              pool_tx) become the last senders and the pool channel can \
              close naturally on shutdown; the split into helper fns needs \
              these"
)]
fn reader_thread_loop(
    reader: &mut AlignmentReader,
    header: &Header,
    work_tx: &WorkTx,
    n_collectors: usize,
    pool_tx: mpsc::Sender<Vec<RikerRecord>>,
    pool_rx: mpsc::Receiver<Vec<RikerRecord>>,
    requirements: &RikerRecordRequirements,
    poison: &AtomicBool,
) -> Result<()> {
    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);

    let result = if reader.supports_in_place_reads() {
        run_in_place_reader(
            reader,
            header,
            work_tx,
            n_collectors,
            &pool_tx,
            &pool_rx,
            requirements,
            poison,
            &mut progress,
        )
    } else {
        run_iterator_reader(
            reader,
            header,
            work_tx,
            n_collectors,
            requirements,
            poison,
            &mut progress,
        )
    };

    progress.finish();
    result
}

/// In-place reader path: pull a pre-allocated batch from the pool, fill
/// it via [`AlignmentReader::fill_record`], dispatch, and let the
/// `RecyclableBatch` return the slots to the pool on drop.
#[allow(clippy::too_many_arguments, reason = "all parameters are needed; this is a private helper")]
fn run_in_place_reader(
    reader: &mut AlignmentReader,
    header: &Header,
    work_tx: &WorkTx,
    n_collectors: usize,
    pool_tx: &mpsc::Sender<Vec<RikerRecord>>,
    pool_rx: &mpsc::Receiver<Vec<RikerRecord>>,
    requirements: &RikerRecordRequirements,
    poison: &AtomicBool,
    progress: &mut ProgressLogger,
) -> Result<()> {
    loop {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }

        // Blocks when every batch is in flight — natural backpressure.
        let Ok(mut records) = pool_rx.recv() else {
            return Ok(());
        };

        let mut len = 0;
        while len < records.len() {
            if !reader.fill_record(requirements, &mut records[len])? {
                break;
            }
            progress.record_with(&records[len], header);
            len += 1;
        }

        if len == 0 {
            let _ = pool_tx.send(records);
            return Ok(());
        }

        let batch = Arc::new(RecyclableBatch { records, len, return_tx: Some(pool_tx.clone()) });
        if !dispatch_batch(work_tx, &batch, n_collectors, poison) {
            return Ok(());
        }
        drop(batch);
    }
}

/// Iterator reader path: noodles owns the record allocation. Used for
/// CRAM, where the batch is one-way (no recycling) and the pool channel
/// is unused.
fn run_iterator_reader(
    reader: &mut AlignmentReader,
    header: &Header,
    work_tx: &WorkTx,
    n_collectors: usize,
    requirements: &RikerRecordRequirements,
    poison: &AtomicBool,
    progress: &mut ProgressLogger,
) -> Result<()> {
    let mut iter = reader.riker_records(requirements);
    let mut records: Vec<RikerRecord> = Vec::with_capacity(BATCH_SIZE);
    loop {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }

        records.clear();
        for _ in 0..BATCH_SIZE {
            match iter.next() {
                Some(Ok(rec)) => {
                    progress.record_with(&rec, header);
                    records.push(rec);
                }
                Some(Err(e)) => return Err(e),
                None => break,
            }
        }

        if records.is_empty() {
            return Ok(());
        }

        let len = records.len();
        let owned = std::mem::replace(&mut records, Vec::with_capacity(BATCH_SIZE));
        let batch = Arc::new(RecyclableBatch { records: owned, len, return_tx: None });
        if !dispatch_batch(work_tx, &batch, n_collectors, poison) {
            return Ok(());
        }
        drop(batch);
    }
}

/// Fan a batch out onto the shared work queue once per collector. Returns
/// `false` when poisoned or when the queue closes (caller should stop).
fn dispatch_batch(
    work_tx: &WorkTx,
    batch: &Batch,
    n_collectors: usize,
    poison: &AtomicBool,
) -> bool {
    for idx in 0..n_collectors {
        if poison.load(Ordering::Relaxed) {
            return false;
        }
        if work_tx.send((idx, Arc::clone(batch))).is_err() {
            // All pool threads gone (they errored and dropped their receivers).
            return false;
        }
    }
    true
}

/// Pool worker: blocks on the shared work queue and dispatches each
/// `(collector_idx, batch)` to the corresponding collector under its mutex.
///
/// Finalization (`Collector::finish`) is handled in [`run_parallel`] after
/// all pool workers have exited, so this loop's only job is to accept
/// batches and return cleanly when the queue closes.
#[allow(
    clippy::needless_pass_by_value,
    reason = "each worker owns its own clone of the MPMC receiver so we can \
              drop the outer handle in run_parallel; passing by reference \
              would leave that handle alive and keep the queue open"
)]
fn pool_worker_loop(
    work_rx: WorkRx,
    slots: &[Mutex<Box<dyn Collector>>],
    header: &Header,
    poison: &AtomicBool,
) -> Result<()> {
    while let Ok((idx, batch)) = work_rx.recv() {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }
        let mut collector = slots[idx].lock().unwrap();
        if let Err(e) = collector.accept_multiple(batch.records(), header) {
            // Signal the reader and any sibling workers to stop so we don't
            // process the rest of the in-flight batches before the reader's
            // send-fails-due-to-no-receivers path triggers shutdown.
            poison.store(true, Ordering::Relaxed);
            return Err(e);
        }
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

        let reader = AlignmentReader::open(Path::new(tmp.path()), None)?;
        let collectors: Vec<Box<dyn Collector>> =
            vec![Box::new(FailingCollector { seen: 0, fail_after: 100 })];

        // 2 threads = 1 reader + 1 worker, the minimal parallel config.
        let result = run_parallel(reader, collectors, 2);
        let err = result.expect_err("run_parallel should propagate the collector error");
        assert!(
            err.to_string().contains("synthetic failure"),
            "expected the failing collector's error, got: {err}"
        );
        Ok(())
    }
}
