use std::time::Instant;

use log::Level;
use noodles::sam::Header;

use crate::sam::riker_record::RikerRecord;

/// Width of the formatted count field, enough for `"1,000,000,000"` (13 chars).
const COUNT_WIDTH: usize = 13;

/// Width of the formatted position field (`chrom:pos`), enough for a 6-char
/// contig name, a colon, and a 9-digit comma-formatted position (11 chars) = 18.
const POS_WIDTH: usize = 18;

/// Format a count with commas as thousands separators (e.g. `1000000` → `"1,000,000"`).
fn format_count(n: u64) -> String {
    let s = n.to_string();
    s.as_bytes()
        .rchunks(3)
        .rev()
        .map(|c| std::str::from_utf8(c).unwrap())
        .collect::<Vec<_>>()
        .join(",")
}

/// Format a `chrom:pos` string with commas in the position.
fn format_position(chrom: &str, pos: u64) -> String {
    format!("{chrom}:{}", format_count(pos))
}

/// Format a duration as `MMm SSs` with zero-padded two-digit minutes and seconds.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "elapsed seconds are non-negative and fit in u64 for any practical duration"
)]
fn format_elapsed(secs: f64) -> String {
    let total_secs = secs as u64;
    let m = total_secs / 60;
    let s = total_secs % 60;
    format!("{m:02}m {s:02}s")
}

/// Logs progress every N records.
pub struct ProgressLogger {
    name: &'static str,
    unit: &'static str,
    every: u64,
    count: u64,
    /// Next count at which to emit a progress message.  Avoids per-call division.
    next_milestone: u64,
    last_milestone: Instant,
    start: Instant,
}

impl ProgressLogger {
    #[must_use]
    pub fn new(name: &'static str, unit: &'static str, every: u64) -> Self {
        let now = Instant::now();
        Self { name, unit, every, count: 0, next_milestone: every, last_milestone: now, start: now }
    }

    /// Record one item and log if the interval has been reached.
    pub fn record(&mut self) {
        self.count += 1;
        if self.count >= self.next_milestone {
            self.next_milestone += self.every;
            self.emit(None);
        }
    }

    /// Record one item at a known genomic position and log if the interval has been reached.
    /// `pos` is displayed as-is (callers should pass 1-based coordinates for consistency).
    pub fn record_with_position(&mut self, chrom: &str, pos: u64) {
        self.count += 1;
        if self.count >= self.next_milestone {
            self.next_milestone += self.every;
            self.emit(Some(&format_position(chrom, pos)));
        }
    }

    /// Record one item and, if the log interval is reached, append the current genomic
    /// position derived from `record` and `header` to the log message.  Unmapped records
    /// are shown as `<unmapped>`.
    ///
    /// The header lookup and string allocation for the position are deferred until
    /// a milestone is actually hit, avoiding ~800M allocations on a 30× WGS run.
    pub fn record_with(&mut self, record: &RikerRecord, header: &Header) {
        self.count += 1;
        if self.count >= self.next_milestone {
            self.next_milestone += self.every;
            let (chrom, pos) = extract_chrom_pos(record, header);
            self.emit(Some(&format_position(&chrom, pos)));
        }
    }

    /// Advance the counter by `n` and emit a progress message if a milestone boundary is
    /// crossed, showing `chrom:pos` as the current location.  If `n` spans multiple
    /// milestone boundaries only one message is emitted (at the current position).
    pub fn record_n_with_position(&mut self, n: u64, chrom: &str, pos: u64) {
        if n == 0 {
            return;
        }
        self.count += n;
        if self.count >= self.next_milestone {
            // Advance past all crossed milestones (only one log message emitted).
            while self.next_milestone <= self.count {
                self.next_milestone += self.every;
            }
            self.emit(Some(&format_position(chrom, pos)));
        }
    }

    /// Log final totals.
    pub fn finish(&self) {
        let total = format_elapsed(self.start.elapsed().as_secs_f64());
        log::log!(
            target: self.name, Level::Info,
            "Processed {:>COUNT_WIDTH$} {} total in {total}.",
            format_count(self.count), self.unit,
        );
    }

    /// Emit a log message with an optional `@ chrom:pos` suffix and timing notes.
    fn emit(&mut self, pos: Option<&str>) {
        let milestone_secs = self.last_milestone.elapsed().as_secs_f64();
        let total_elapsed = format_elapsed(self.start.elapsed().as_secs_f64());
        let last_took = format!("last {} took {:.1}s", format_count(self.every), milestone_secs);

        match pos {
            None => log::log!(
                target: self.name, Level::Info,
                "Processed {:>COUNT_WIDTH$} {} - elapsed time {total_elapsed} - {last_took}.",
                format_count(self.count), self.unit,
            ),
            Some(p) => log::log!(
                target: self.name, Level::Info,
                "Processed {:>COUNT_WIDTH$} {} @ {:<POS_WIDTH$} - elapsed time {total_elapsed} - {last_took}.",
                format_count(self.count), self.unit, p,
            ),
        }
        self.last_milestone = Instant::now();
    }
}

/// Extract a `(chrom, 1-based-pos)` pair from a record for use in log messages.
fn extract_chrom_pos(record: &RikerRecord, header: &Header) -> (String, u64) {
    let chrom = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .map_or_else(|| "<unmapped>".to_string(), |(name, _)| name.to_string());
    let pos = record.alignment_start().map_or(0, usize::from) as u64;
    (chrom, pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_count_small() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(1), "1");
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_thousands() {
        assert_eq!(format_count(1_000), "1,000");
        assert_eq!(format_count(10_000), "10,000");
        assert_eq!(format_count(100_000), "100,000");
    }

    #[test]
    fn test_format_count_millions() {
        assert_eq!(format_count(1_000_000), "1,000,000");
        assert_eq!(format_count(10_000_000), "10,000,000");
        assert_eq!(format_count(1_234_567_890), "1,234,567,890");
    }

    // ── format_position ──────────────────────────────────────────────────────

    #[test]
    fn test_format_position_small() {
        assert_eq!(format_position("chr1", 100), "chr1:100");
    }

    #[test]
    fn test_format_position_large() {
        assert_eq!(format_position("chr1", 123_456_789), "chr1:123,456,789");
    }

    #[test]
    fn test_format_position_zero() {
        assert_eq!(format_position("chrX", 0), "chrX:0");
    }

    // ── format_elapsed ───────────────────────────────────────────────────────

    #[test]
    fn test_format_elapsed_seconds_only() {
        assert_eq!(format_elapsed(5.3), "00m 05s");
    }

    #[test]
    fn test_format_elapsed_minutes_and_seconds() {
        assert_eq!(format_elapsed(125.7), "02m 05s");
    }

    #[test]
    fn test_format_elapsed_exact_minute() {
        assert_eq!(format_elapsed(60.0), "01m 00s");
    }

    // ── ProgressLogger ───────────────────────────────────────────────────────

    #[test]
    fn test_record_no_panic() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        for _ in 0..25 {
            pl.record();
        }
        pl.finish();
    }

    #[test]
    fn test_record_with_position_no_panic() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        for i in 0..25 {
            pl.record_with_position("chr1", i);
        }
        pl.finish();
    }

    #[test]
    fn test_record_n_zero_is_noop() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        pl.record_n_with_position(0, "chr1", 1);
        pl.record_n_with_position(0, "chr1", 2);
        pl.record_n_with_position(0, "chr1", 3);
        pl.finish();
    }

    #[test]
    fn test_record_n_crosses_milestones() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        // Jump past multiple milestones in one call
        pl.record_n_with_position(35, "chr1", 100);
        pl.finish();
    }

    // ── extract_chrom_pos ────────────────────────────────────────────────────

    #[test]
    fn test_extract_chrom_pos_mapped() {
        use noodles::core::Position;
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;

        let header = Header::builder()
            .add_reference_sequence(
                b"chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build();

        let rec = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(100).unwrap())
            .build();
        let riker = RikerRecord::from_alignment_record(&header, &rec).unwrap();

        let (chrom, pos) = extract_chrom_pos(&riker, &header);
        assert_eq!(chrom, "chr1");
        assert_eq!(pos, 100);
    }

    #[test]
    fn test_extract_chrom_pos_unmapped() {
        use noodles::sam::alignment::RecordBuf;

        let header = Header::default();
        let rec = RecordBuf::default();
        let riker = RikerRecord::from_alignment_record(&header, &rec).unwrap();
        let (chrom, pos) = extract_chrom_pos(&riker, &header);
        assert_eq!(chrom, "<unmapped>");
        assert_eq!(pos, 0);
    }
}
