use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError, sync_channel};
use std::sync::{Arc, Condvar, Mutex};
use std::time::Duration;

use anyhow::{Result, anyhow};
use clap::{Args, ValueEnum};
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

/// Wraps a collector with its channel receiver, header, and completion flag.
struct CollectorSlot {
    collector: Box<dyn Collector>,
    rx: Receiver<Batch>,
    header: Header,
    done: bool,
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

/// Run collectors in parallel with a dedicated reader thread and a pool of worker threads.
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

    // Create per-collector channels and slots.
    let mut senders = Vec::with_capacity(collectors.len());
    let mut slots = Vec::with_capacity(collectors.len());

    for collector in collectors {
        let (tx, rx) = sync_channel(CHANNEL_DEPTH);
        senders.push(tx);
        slots.push(Mutex::new(CollectorSlot {
            collector,
            rx,
            header: header.clone(),
            done: false,
        }));
    }

    // Pool of reusable record-batch allocations. When a `RecyclableBatch`
    // drops (last collector done), its inner `Vec<RecordBuf>` is sent back
    // here so the reader can reuse those pre-allocated slots.
    //
    // Unbounded channel: receivers may arrive at arbitrary times during
    // shutdown and we must never block in `RecyclableBatch::drop`.
    let (pool_tx, pool_rx) = mpsc::channel::<Vec<RecordBuf>>();
    let pool_size = CHANNEL_DEPTH + POOL_EXTRA;
    for _ in 0..pool_size {
        let mut vec: Vec<RecordBuf> = Vec::with_capacity(BATCH_SIZE);
        vec.resize_with(BATCH_SIZE, RecordBuf::default);
        pool_tx.send(vec).expect("pool receiver exists on startup, send cannot fail");
    }

    let poison = AtomicBool::new(false);
    let notify = (Mutex::new(()), Condvar::new());

    // Use scoped threads so everything can borrow from this stack frame.
    std::thread::scope(|scope| {
        // Spawn pool threads.
        let mut pool_handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            pool_handles.push(scope.spawn(|| pool_thread_loop(&slots, &poison, &notify)));
        }

        // Reader thread runs on the current scope too. The pool ends, the
        // senders, and the reader all move into it; `poison`/`notify` are
        // shared with the pool threads so we pass references (which the
        // `move` closure captures by copy).
        let poison_ref: &AtomicBool = &poison;
        let notify_ref: &(Mutex<()>, Condvar) = &notify;
        let reader_handle = scope.spawn(move || {
            let reader_result = reader_thread_loop(
                &mut reader,
                header,
                &senders,
                pool_tx,
                pool_rx,
                poison_ref,
                notify_ref,
            );

            // Drop all senders so pool threads see channel disconnection.
            drop(senders);

            // Wake pool threads so they notice the disconnection.
            notify_ref.1.notify_all();

            reader_result
        });

        // Collect results — reader first, then pool threads.
        let reader_result = reader_handle.join().map_err(|_| anyhow!("reader thread panicked"))?;
        if let Err(e) = reader_result {
            // Ensure pool threads stop.
            poison.store(true, Ordering::Relaxed);
            notify.1.notify_all();
            // Wait for pool threads to finish before returning the error.
            for handle in pool_handles {
                let _ = handle.join();
            }
            return Err(e);
        }

        // Collect pool thread results.
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
    })
}

/// Reader thread: reads BAM records in batches and dispatches to all collector channels.
///
/// Records are read directly into pre-allocated `RecordBuf`s pulled from
/// `pool_rx`; when the last collector has processed the resulting
/// [`RecyclableBatch`], the inner `Vec` flows back to `pool_tx` via the
/// batch's `Drop`, so we never pay per-record allocation costs after startup.
#[allow(
    clippy::needless_pass_by_value,
    reason = "pool_tx and pool_rx are moved into the reader thread's scope so \
              that in-flight RecyclableBatch Drops can clone pool_tx and pool_rx \
              is !Sync"
)]
fn reader_thread_loop(
    reader: &mut AlignmentReader,
    header: &Header,
    senders: &[std::sync::mpsc::SyncSender<Batch>],
    pool_tx: mpsc::Sender<Vec<RecordBuf>>,
    pool_rx: mpsc::Receiver<Vec<RecordBuf>>,
    poison: &AtomicBool,
    notify: &(Mutex<()>, Condvar),
) -> Result<()> {
    let mut progress = ProgressLogger::new("multi", "reads", 5_000_000);

    loop {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }

        // Pull an empty batch from the pool. Blocks if every batch is currently
        // in flight, which is the backpressure mechanism.
        let Ok(mut records) = pool_rx.recv() else {
            // Reader's pool was closed — shouldn't happen during normal run
            // since we hold pool_tx ourselves.
            return Ok(());
        };

        // Read up to BATCH_SIZE records directly into the pre-allocated slots.
        // Each RecordBuf's inner Vecs (name, cigar, sequence, quality, data)
        // are reused across records.
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
            // EOF: return the unused batch to the pool and exit.
            let _ = pool_tx.send(records);
            break;
        }

        // Wrap and dispatch. When the last collector drops its Arc, the
        // batch's Drop sends `records` back to `pool_tx` for reuse.
        let batch = Arc::new(RecyclableBatch { records, len, return_tx: pool_tx.clone() });
        for tx in senders {
            if poison.load(Ordering::Relaxed) {
                return Ok(());
            }
            if tx.send(Arc::clone(&batch)).is_err() {
                // Receiver dropped — a pool thread errored.
                return Ok(());
            }
        }
        drop(batch);
        notify.1.notify_all();
    }

    progress.finish();
    Ok(())
}

/// Pool thread work loop: services any collector that has pending batches.
fn pool_thread_loop(
    slots: &[Mutex<CollectorSlot>],
    poison: &AtomicBool,
    notify: &(Mutex<()>, Condvar),
) -> Result<()> {
    loop {
        if poison.load(Ordering::Relaxed) {
            return Ok(());
        }

        let mut did_work = false;
        let mut all_done = true;

        for slot_mutex in slots {
            if let Ok(mut slot) = slot_mutex.try_lock() {
                if slot.done {
                    continue;
                }
                all_done = false;

                match slot.rx.try_recv() {
                    Ok(batch) => {
                        let CollectorSlot { collector, header, .. } = &mut *slot;
                        collector.accept_multiple(batch.records(), header)?;
                        did_work = true;
                    }
                    Err(TryRecvError::Empty) => {}
                    Err(TryRecvError::Disconnected) => {
                        // Channel closed — drain any remaining batches then finish.
                        {
                            let CollectorSlot { collector, rx, header, .. } = &mut *slot;
                            while let Ok(batch) = rx.try_recv() {
                                collector.accept_multiple(batch.records(), header)?;
                            }
                        }
                        if !poison.load(Ordering::Relaxed) {
                            slot.collector.finish()?;
                        }
                        slot.done = true;
                    }
                }
            } else {
                // Locked by another thread — not done yet.
                all_done = false;
            }
        }

        if all_done {
            return Ok(());
        }

        if !did_work {
            let guard = notify.0.lock().unwrap();
            let _ = notify.1.wait_timeout(guard, Duration::from_millis(1)).unwrap();
        }
    }
}
