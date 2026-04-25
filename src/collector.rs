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
}

/// Drive a single reader through a single collector's full lifecycle:
/// `initialize(header)` → stream every record through `accept` → `finish()`.
/// Uses [`AlignmentReader::fill_record`] for BAM/SAM (no per-record
/// allocation) and falls back to [`AlignmentReader::riker_records`] for
/// CRAM (which noodles does not let us read in place).
///
/// `progress.finish()` is called unconditionally before returning, so the
/// "Processed N total" line appears on both the success and error paths.
///
/// `Collector::finish` always runs — even when initialize or the read
/// loop errored — so collectors can release resources and flush partial
/// output. If both the read loop and `finish` error, the read-loop
/// error wins (it's the upstream cause).
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

    if let Err(e) = collector.initialize(&header) {
        // Initialize failed before any records were processed; still call
        // finish in case it has cleanup to do, then propagate the
        // initialize error.
        let _ = collector.finish();
        return Err(e);
    }

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
    if reader.supports_in_place_reads() {
        let mut record = reader.empty_record();
        while reader.fill_record(&requirements, &mut record)? {
            progress.record_with(&record, header);
            collector.accept(&record, header)?;
        }
    } else {
        for result in reader.riker_records(&requirements) {
            let record = result?;
            progress.record_with(&record, header);
            collector.accept(&record, header)?;
        }
    }
    Ok(())
}
