use anyhow::Result;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;

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
    fn accept(&mut self, record: &RecordBuf, header: &Header) -> Result<()>;

    /// Process a batch of records. The default iterates calling `accept()` for each record.
    /// Since the default calls `self.accept()` on the concrete type (not through the vtable),
    /// this reduces virtual dispatch from once-per-record to once-per-batch.
    ///
    /// # Errors
    /// Returns an error if any record cannot be processed.
    fn accept_multiple(&mut self, records: &[RecordBuf], header: &Header) -> Result<()> {
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
}
