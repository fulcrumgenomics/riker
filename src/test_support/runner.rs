//! Drive a [`Collector`] over a BAM, so test files stop repeating the
//! open → initialize → accept-loop → finish boilerplate (the loop the review
//! found copy-pasted in the alignment, basic, and bam_reader tests).

use std::path::Path;

use anyhow::Result;

use crate::collector::Collector;
use crate::sam::alignment_reader::AlignmentReader;

/// Feed every record in `bam` through `collector`: open the reader, initialize
/// with the header, accept each record, then finish. The caller constructs the
/// collector (with its output paths) and reads its outputs; this owns only the
/// read loop.
///
/// # Errors
/// Propagates reader-open, record-decode, and collector errors.
pub fn drive_collector(collector: &mut impl Collector, bam: &Path) -> Result<()> {
    let mut reader = AlignmentReader::open(bam, None, 0)?;
    let header = reader.header().clone();
    collector.initialize(&header)?;
    let requirements = collector.field_needs();
    for record in reader.riker_records(&requirements) {
        collector.accept(&record?, &header)?;
    }
    collector.finish()?;
    Ok(())
}
