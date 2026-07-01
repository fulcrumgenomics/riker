use anyhow::Result;
use noodles::sam::Header;

use crate::progress::ProgressLogger;
use crate::sam::alignment_reader::AlignmentReader;
use crate::sam::riker_record::{RikerRecord, RikerRecordRequirements};

/// Trait implemented by each metric collector.
///
/// Each collector stores its own configuration (output paths, reference handle, thresholds)
/// as fields set at construction time. The trait methods only receive the BAM header and
/// records, enabling the `multi` command to share a single BAM pass across collectors.
pub trait Collector: Send {
    /// Called once with the BAM header before any records are processed.
    ///
    /// # Errors
    /// Returns an error if the header is invalid for this collector's configuration.
    fn initialize(&mut self, header: &Header) -> Result<()>;

    /// Called once per record in the BAM file.
    ///
    /// # Errors
    /// Returns an error if the record cannot be processed.
    fn accept(&mut self, record: &RikerRecord, header: &Header) -> Result<()>;

    /// Process a batch of records. The default loops over `accept`. The
    /// existence of this method is the optimization opportunity, not the
    /// default body: a concrete collector that overrides `accept_multiple`
    /// can amortize per-batch setup (vtable thaws, allocation reuse,
    /// inlining) over `records.len()` records, since the parallel `multi`
    /// pipeline always dispatches through this entry point. The default
    /// implementation, called through `Box<dyn Collector>`, still goes
    /// through the vtable per record.
    ///
    /// # Errors
    /// Returns an error if any record cannot be processed.
    fn accept_multiple(&mut self, records: &[RikerRecord], header: &Header) -> Result<()> {
        for record in records {
            self.accept(record, header)?;
        }
        Ok(())
    }

    /// Called after all records have been processed. Should write output files.
    ///
    /// # Errors
    /// Returns an error if metrics cannot be computed or output files cannot be written.
    fn finish(&mut self) -> Result<()>;

    /// Short name identifying this collector (used for logging).
    fn name(&self) -> &'static str;

    /// Declare which expensive-to-populate fields this collector will read
    /// through [`crate::sam::riker_record::RikerRecord`] accessors — the
    /// reader thread uses the union of these across the active collectors
    /// to drive its per-record decode work.
    ///
    /// Required (no default). Every collector must explicitly say what it
    /// needs so we don't silently forget to opt into an expensive field
    /// when adding a new collector. Return
    /// [`RikerRecordRequirements::NONE`] if no expensive fields are needed.
    fn field_needs(&self) -> RikerRecordRequirements;

    /// Relative per-record processing cost, used only by `multi` to balance
    /// collectors across worker threads (heaviest collectors get their own
    /// worker). This is a self-reported hint, not a hard number — the
    /// scheduler compares hints, so only the ratios matter. Default `1`
    /// (a light, accumulate-only collector). Collectors with expensive
    /// per-record kernels (e.g. coverage pileup) should override with a
    /// larger value.
    fn cost_hint(&self) -> u32 {
        1
    }
}

/// Drive a single reader through a single collector's full lifecycle:
/// `initialize(header)` → stream every record through `accept` → `finish()`.
/// Uses [`AlignmentReader::fill_record`] for every format (no per-record
/// allocation).
///
/// If `initialize` fails the function returns immediately — no records
/// have been read, so neither `progress.finish()` nor `Collector::finish`
/// runs. Real collectors populate fields in `initialize` that `finish`
/// then unwraps, so calling finish on a partially-constructed collector
/// would turn an honest error into a panic.
///
/// Once `initialize` succeeds, `progress.finish()` is always called
/// before return (so the "Processed N total" line appears on both the
/// success and read-loop-error paths) and `Collector::finish` is always
/// called so it can flush partial output. If both the read loop and
/// `finish` error, the read-loop error wins (it's the upstream cause).
///
/// # Errors
/// Returns an error if the underlying reader, decoder, or any of
/// `initialize` / `accept` / `finish` fail.
pub fn drive_collector_single_threaded(
    reader: &mut AlignmentReader,
    collector: &mut dyn Collector,
    progress: &mut ProgressLogger,
) -> Result<()> {
    // Clone the header once so it can be borrowed alongside a `&mut reader`
    // held by the per-record iterator. Header clone is one shot at startup.
    let header = reader.header().clone();

    collector.initialize(&header)?;

    let read_result = drive_records(reader, collector, progress, &header);
    progress.finish();
    let finish_result = collector.finish();

    // Read-loop error wins — it's the upstream cause and what the user
    // most needs to see. Finish errors only surface if the loop ran clean.
    read_result?;
    finish_result
}

fn drive_records(
    reader: &mut AlignmentReader,
    collector: &mut dyn Collector,
    progress: &mut ProgressLogger,
    header: &Header,
) -> Result<()> {
    let requirements = collector.field_needs();
    let mut record = reader.empty_record();
    while reader.fill_record(&requirements, &mut record)? {
        progress.record_with(&record, header);
        collector.accept(&record, header)?;
    }
    Ok(())
}
