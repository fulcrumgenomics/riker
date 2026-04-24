use std::fmt;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use anyhow::{Result, anyhow};
use clap::{Args, ValueEnum};
use crossbeam_channel::{Receiver, Sender};
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;

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

/// Number of records per batch sent through channels.
const BATCH_SIZE: usize = 256;

/// Number of batches buffered per collector channel.
const CHANNEL_DEPTH: usize = 128;

/// Extra pool slots beyond `CHANNEL_DEPTH` to absorb the in-flight batch and
/// any transient ordering between dispatch and per-collector drain.
const POOL_EXTRA: usize = 2;

/// A batch of records shared across collector channels.
///
/// Wrapping the records in `RecyclableBatch` lets us send the inner
/// `Vec<RecordBuf>` back to the reader's pool when the last `Arc` reference
/// drops (i.e. the last collector has finished with it), so the reader can
/// reuse the pre-allocated `RecordBuf` slots on the next read and avoid the
/// per-record clone that `Reader::record_bufs` would do.
type Batch = Arc<RecyclableBatch>;

/// A single work item on the MPMC work queue: which collector the batch is
/// destined for, and the shared batch itself.
type WorkItem = (usize, Batch);
type WorkTx = Sender<WorkItem>;
type WorkRx = Receiver<WorkItem>;

/// Owns a pre-allocated `Vec<RecordBuf>` of capacity `BATCH_SIZE` plus a count
/// of valid records. On drop, the inner `Vec` is returned to the reader's
/// pool via `return_tx` so its allocations can be reused.
struct RecyclableBatch {
    records: Vec<RecordBuf>,
    len: usize,
    return_tx: mpsc::Sender<Vec<RecordBuf>>,
}

impl RecyclableBatch {
    /// Valid records in the batch, as a slice.
    fn records(&self) -> &[RecordBuf] {
        &self.records[..self.len]
    }
}

impl Drop for RecyclableBatch {
    fn drop(&mut self) {
        // Hand the inner Vec back to the reader's pool for reuse. The receiver
        // is dropped during shutdown; ignore send errors in that case.
        let records = std::mem::take(&mut self.records);
        let _ = self.return_tx.send(records);
    }
}

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
  riker multi -i input.bam -o out_prefix -R ref.fa
  riker multi -i input.bam -o out_prefix -R ref.fa --tools alignment basic isize
  riker multi -i input.bam -o out_prefix -R ref.fa --threads 4
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

    /// Number of threads for parallel collection. 1 = single-threaded (no threading overhead).
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

        // Pre-read the header so we can build interval maps for WGS if needed.
        let (reader, header) =
            AlignmentReader::new(&self.input.input, self.reference.reference.as_deref())?;

        let collectors = self.build_collectors(&seen, &header)?;

        if self.threads == 1 {
            run_single_threaded(reader, &header, collectors)?;
        } else {
            run_parallel(reader, &header, collectors, self.threads)?;
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

/// Run all collectors sequentially on a single thread (no threading overhead).
fn run_single_threaded(
    mut reader: AlignmentReader,
    header: &Header,
    mut collectors: Vec<Box<dyn Collector>>,
) -> Result<()> {
    for collector in &mut collectors {
        collector.initialize(header)?;
    }

    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);
    reader.for_each_record(header, |record| {
        progress.record_with(record, header);
        for collector in &mut collectors {
            collector.accept(record, header)?;
        }
        Ok(())
    })?;
    progress.finish();

    for collector in &mut collectors {
        collector.finish()?;
    }

    Ok(())
}

/// Run collectors in parallel with a dedicated reader thread and a shared
/// worker pool.
///
/// Architecture:
/// - One **reader thread** pulls empty [`RecyclableBatch`] slots from a pool,
///   fills them in place via `read_record_buf`, wraps each in an `Arc`, and
///   fans the `Arc` out onto a shared `(collector_idx, Arc<Batch>)` work queue
///   — one entry per collector.
/// - `threads` **pool threads** block on the shared work queue. On receipt of
///   a work item they lock that collector's mutex and call `accept_multiple`.
///   The mutex serializes per-collector accepts (required since `Collector`
///   is stateful) while still allowing different collectors' batches to
///   process in parallel across threads.
/// - When the last `Arc<RecyclableBatch>` drops, its `Drop` returns the inner
///   `Vec<RecordBuf>` to the reader's pool via an mpsc channel.
///
/// Properties vs. the previous try_lock+sleep design: no busy-polling, no
/// `Condvar`, and `--threads` directly controls how many worker threads run.
fn run_parallel(
    mut reader: AlignmentReader,
    header: &Header,
    mut collectors: Vec<Box<dyn Collector>>,
    threads: usize,
) -> Result<()> {
    // Initialize all collectors before spawning threads.
    for collector in &mut collectors {
        collector.initialize(header)?;
    }

    // Each collector lives behind a Mutex so a pool thread can claim exclusive
    // access before calling `accept_multiple`. The mutex only serializes accesses
    // to a single collector; different collectors process in parallel.
    let slots: Vec<Mutex<Box<dyn Collector>>> = collectors.into_iter().map(Mutex::new).collect();

    // Shared MPMC work queue: reader sends (collector_idx, batch) items; pool
    // threads pull and dispatch to the matching collector.
    let (work_tx, work_rx): (WorkTx, WorkRx) =
        crossbeam_channel::bounded(CHANNEL_DEPTH * slots.len().max(1));

    // Pool of reusable record-batch allocations. When a `RecyclableBatch`
    // drops (all collectors done with it), its inner `Vec<RecordBuf>` is
    // sent back here so the reader can reuse those pre-allocated slots.
    // Unbounded so `RecyclableBatch::drop` never blocks.
    let (pool_tx, pool_rx) = mpsc::channel::<Vec<RecordBuf>>();
    let pool_size = CHANNEL_DEPTH + POOL_EXTRA;
    for _ in 0..pool_size {
        let mut vec: Vec<RecordBuf> = Vec::with_capacity(BATCH_SIZE);
        vec.resize_with(BATCH_SIZE, RecordBuf::default);
        pool_tx.send(vec).expect("pool receiver exists on startup, send cannot fail");
    }

    // Poison: set when any thread errors. Signals the reader to stop early.
    let poison = AtomicBool::new(false);
    let n_collectors = slots.len();

    // Run the reader + worker threads inside a scope so they can borrow
    // `slots` and `poison`. The scope ends (and borrows release) before we
    // reclaim `slots` to call `finish()`.
    std::thread::scope(|scope| -> Result<()> {
        let slots_ref: &[Mutex<Box<dyn Collector>>] = &slots;
        let poison_ref: &AtomicBool = &poison;

        // Pool workers share `work_rx` (MPMC) and the slots.
        let mut pool_handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let work_rx = work_rx.clone();
            pool_handles.push(
                scope.spawn(move || pool_worker_loop(work_rx, slots_ref, header, poison_ref)),
            );
        }
        // Drop our outer receiver handle so the queue closes once the reader's
        // sender drops (otherwise pool threads would wait forever).
        drop(work_rx);

        // Reader owns the work-queue sender and the batch pool.
        let reader_handle = scope.spawn(move || {
            let reader_result = reader_thread_loop(
                &mut reader,
                header,
                &work_tx,
                n_collectors,
                pool_tx,
                pool_rx,
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

/// Reader thread: reads BAM records in batches into pooled `RecordBuf` slots
/// and fans each batch onto the shared work queue once per collector.
#[allow(
    clippy::needless_pass_by_value,
    reason = "pool_tx and pool_rx are moved into the reader thread: in-flight \
              RecyclableBatch Drops clone pool_tx, and pool_rx is !Sync so it \
              cannot be shared by reference"
)]
fn reader_thread_loop(
    reader: &mut AlignmentReader,
    header: &Header,
    work_tx: &WorkTx,
    n_collectors: usize,
    pool_tx: mpsc::Sender<Vec<RecordBuf>>,
    pool_rx: mpsc::Receiver<Vec<RecordBuf>>,
    poison: &AtomicBool,
) -> Result<()> {
    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);

    loop {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }

        // Blocks if every batch is currently in flight — acts as backpressure.
        let Ok(mut records) = pool_rx.recv() else {
            return Ok(());
        };

        // Read directly into the pre-allocated slots. Each RecordBuf's inner
        // Vecs (name, cigar, sequence, quality, data) are reused.
        let mut len = 0;
        while len < records.len() {
            let bytes = reader.read_record_buf(header, &mut records[len])?;
            if bytes == 0 {
                break;
            }
            progress.record_with(&records[len], header);
            len += 1;
        }

        if len == 0 {
            // EOF: recycle the unused batch and stop.
            let _ = pool_tx.send(records);
            break;
        }

        // Wrap once, fan out via Arc clones on the shared work queue. The
        // last `Arc` drop (after the last pool worker processes it) triggers
        // `RecyclableBatch::drop`, which recycles `records` to the pool.
        let batch = Arc::new(RecyclableBatch { records, len, return_tx: pool_tx.clone() });
        for idx in 0..n_collectors {
            if poison.load(Ordering::Relaxed) {
                return Ok(());
            }
            if work_tx.send((idx, Arc::clone(&batch))).is_err() {
                // All pool threads gone (they errored and dropped their
                // receivers) — stop and let the parent surface the error.
                return Ok(());
            }
        }
        drop(batch);
    }

    progress.finish();
    Ok(())
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
        collector.accept_multiple(batch.records(), header)?;
    }
    Ok(())
}
