use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use noodles::core::Region;
use noodles::fasta;
use noodles::sam::Header;

use crate::sequence_dict::SequenceDictionary;

/// A loaded FASTA reference with random-access sequence retrieval.
///
/// Wraps a noodles `IndexedReader` for efficient seek-based access.  Contig
/// metadata (names, lengths) is pre-computed from the `.fai` index at
/// construction time.
///
/// Whole-contig loads go through [`load_contig_into`](Self::load_contig_into),
/// a bulk `read_exact` of the contig's byte span (from the `.fai` geometry)
/// followed by in-place newline stripping into the caller's buffer. This is both
/// faster than noodles' line-by-line reader and avoids its per-read zero-fill; the
/// caller owns the buffer, so `Fasta` holds no per-load scratch.
pub struct Fasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    dict: SequenceDictionary,
    /// Per-contig `.fai` geometry in index order, parallel to `dict`, copied out
    /// of the index so a load can borrow the reader mutably without holding the
    /// index borrow. Looked up by the dict's name→index map (reused rather than
    /// keeping a second name→geometry map).
    geom: Vec<FaiGeom>,
}

impl Fasta {
    /// Open a FASTA file and its `.fai` index.
    ///
    /// # Panics
    /// Panics if a `.fai` contig length or line geometry exceeds `usize`
    /// (should never happen: `lib.rs` rejects sub-64-bit targets at compile
    /// time, so every `u64` from the index fits `usize`).
    ///
    /// # Errors
    /// Returns an error if the FASTA or its index cannot be read, or if the
    /// `.fai` line geometry is corrupt (`line_width < line_bases`, which would
    /// mis-stride the bulk read).
    pub fn from_path(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed FASTA: {}", path.display()))?;

        let dict = SequenceDictionary::from(reader.index());

        // Build geometry in .fai record order so it stays index-aligned with
        // `dict` (both come from the same records); loads resolve name→index via
        // the dict, then index this Vec. Names live only in `dict`.
        let mut geom = Vec::with_capacity(dict.len());
        for rec in reader.index().as_ref() {
            let length = usize::try_from(rec.length()).expect("contig length fits in usize");
            let line_bases =
                usize::try_from(rec.line_bases()).expect("fai line_bases fits in usize");
            let line_width =
                usize::try_from(rec.line_width()).expect("fai line_width fits in usize");
            // A well-formed .fai has line_width >= line_bases > 0 for any non-empty
            // contig (the difference is the line terminator). A malformed index that
            // violates this would make the bulk read mis-stride, so reject it.
            if length > 0 && (line_bases == 0 || line_width < line_bases) {
                let name = String::from_utf8_lossy(rec.name().as_ref());
                bail!(
                    "Corrupt .fai geometry for contig '{name}': line_bases={line_bases}, \
                     line_width={line_width} (need line_width >= line_bases > 0)"
                );
            }
            geom.push(FaiGeom { length, offset: rec.offset(), line_bases, line_width });
        }

        Ok(Self { reader, dict, geom })
    }

    /// Return a reference to the underlying sequence dictionary.
    #[must_use]
    pub fn dict(&self) -> &SequenceDictionary {
        &self.dict
    }

    /// Verify that every reference sequence in `header` is present in the index.
    ///
    /// # Errors
    /// Returns an error if any BAM contig is absent from the reference index.
    pub fn validate_bam_header(&self, header: &Header) -> Result<()> {
        for (name, _) in header.reference_sequences() {
            let name_str = std::str::from_utf8(name.as_ref())
                .with_context(|| "BAM contig name is not valid UTF-8")?;
            if self.dict.get_by_name(name_str).is_none() {
                return Err(anyhow!(
                    "BAM contig '{name_str}' not found in reference FASTA index. \
                     All BAM contigs must be present in the reference."
                ));
            }
        }
        Ok(())
    }

    /// Load the full sequence of `contig_name` into `out`.
    ///
    /// `out` is filled with the contig's bases (newlines stripped) and truncated
    /// to the contig length; when `uppercase` is true the bases are upper-cased in
    /// place. `out` doubles as the read buffer — the raw contig span (bases plus
    /// line terminators) is read straight into it and the terminators are then
    /// squeezed out in place — so `Fasta` holds no scratch of its own and the
    /// caller owns the buffer. Reusing one `out` across many contigs avoids
    /// per-contig allocation, and because `read_exact` overwrites the whole raw
    /// span, `resize` only ever zero-fills genuine growth past the largest contig
    /// seen so far (no per-load `memset`).
    ///
    /// # Errors
    /// Returns an error if the contig is unknown or the FASTA cannot be read.
    pub fn load_contig_into(
        &mut self,
        contig_name: &str,
        uppercase: bool,
        out: &mut Vec<u8>,
    ) -> Result<()> {
        // Resolve name→index via the shared dict (its FxHashMap map), then index
        // the parallel geometry Vec. Copy the fields out (all `Copy`) so the
        // borrow ends before the mutable seek/read below.
        let idx = self
            .dict
            .get_by_name(contig_name)
            .ok_or_else(|| anyhow!("Contig '{contig_name}' not found in reference"))?
            .index();
        let geom = &self.geom[idx];
        let (length, offset, line_bases, line_width) =
            (geom.length, geom.offset, geom.line_bases, geom.line_width);

        if length == 0 {
            out.clear();
            return Ok(());
        }

        // Read the raw span (bases + terminators) directly into `out`. `resize`
        // reuses the existing allocation and zero-fills only growth beyond the
        // current length; `read_exact` then overwrites all of `[..total_raw]`, so a
        // reused buffer sees no memset after the largest contig.
        let total_raw = total_sequence_bytes(length, line_bases, line_width);
        out.resize(total_raw, 0);
        let reader = self.reader.get_mut();
        reader
            .seek(SeekFrom::Start(offset))
            .with_context(|| format!("Failed to seek to contig '{contig_name}' in FASTA"))?;
        reader
            .read_exact(&mut out[..total_raw])
            .with_context(|| format!("Failed to read contig '{contig_name}' from FASTA"))?;

        // Squeeze out the line terminators in place: shift each line's `line_bases`
        // bases left over the preceding terminators. The write cursor `dst` never
        // overtakes the read cursor `src` (`dst = k*line_bases <= k*line_width =
        // src`), so `copy_within` is a safe forward (memmove) copy and no source is
        // overwritten before it is read.
        // `from_path` guarantees `line_bases > 0` for any `length > 0` contig, so
        // `n` is always positive and the loop terminates (n == 0 would spin).
        debug_assert!(line_bases > 0, "line_bases must be > 0 for a non-empty contig");
        let terminator = line_width - line_bases;
        let (mut dst, mut src, mut remaining) = (0usize, 0usize, length);
        while remaining > 0 {
            let n = remaining.min(line_bases);
            out.copy_within(src..src + n, dst);
            dst += n;
            src += n + terminator;
            remaining -= n;
        }
        out.truncate(length);

        if uppercase {
            out.make_ascii_uppercase();
        }
        Ok(())
    }

    /// Load and return the full sequence of `contig_name`.
    ///
    /// Convenience wrapper over [`load_contig_into`](Self::load_contig_into) for
    /// one-shot callers; prefer `load_contig_into` with a reused buffer in loops.
    ///
    /// # Errors
    /// Returns an error if the contig is unknown, or if the FASTA cannot be read.
    pub fn load_contig(&mut self, contig_name: &str, uppercase: bool) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.load_contig_into(contig_name, uppercase, &mut out)?;
        Ok(out)
    }

    /// Fetch the sequence for a specific region of a contig.
    ///
    /// Coordinates are 0-based half-open `[start, end)`, matching the convention
    /// used throughout riker.  Internally converted to the 1-based inclusive
    /// positions that noodles expects.
    ///
    /// # Errors
    /// Returns an error if the region cannot be read from the FASTA.
    pub fn fetch(&mut self, contig_name: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        use noodles::core::Position;

        #[expect(
            clippy::cast_possible_truncation,
            reason = "genomic coordinates fit in usize on all supported platforms"
        )]
        let pos_start = Position::new((start + 1) as usize)
            .ok_or_else(|| anyhow!("Invalid start position {start} for region query"))?;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "genomic coordinates fit in usize on all supported platforms"
        )]
        let pos_end = Position::new(end as usize)
            .ok_or_else(|| anyhow!("Invalid end position {end} for region query"))?;

        let region = Region::new(contig_name, pos_start..=pos_end);
        let record = self
            .reader
            .query(&region)
            .with_context(|| format!("Failed to fetch {contig_name}:{start}-{end} from FASTA"))?;

        Ok(record.sequence().as_ref().to_vec())
    }

    /// Return the length of the named contig, or `None` if unknown.
    #[must_use]
    pub fn contig_length(&self, name: &str) -> Option<u64> {
        self.dict.get_by_name(name).map(|m| m.length() as u64)
    }

    /// Return an ordered slice of all contig names in the index.
    #[must_use]
    pub fn contig_names(&self) -> Vec<&str> {
        self.dict.names()
    }
}

/// One contig's `.fai` line geometry plus its byte offset and length.
struct FaiGeom {
    length: usize,
    offset: u64,
    line_bases: usize,
    line_width: usize,
}

/// Number of raw bytes (bases plus line terminators) a contig of `seq_len`
/// bases occupies on disk, given its `.fai` `line_bases` and `line_width`.
fn total_sequence_bytes(seq_len: usize, line_bases: usize, line_width: usize) -> usize {
    if seq_len == 0 || line_bases == 0 {
        return 0;
    }
    let complete = seq_len / line_bases;
    let remainder = seq_len % line_bases;
    if remainder > 0 {
        complete * line_width + remainder
    } else {
        // Exactly fills `complete` lines; the last has no trailing terminator.
        (complete - 1) * line_width + line_bases
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    /// Write a tiny indexed FASTA (`.fa` + `.fai`) to a temp dir and return the
    /// `.fa` path. `records` is `(name, sequence)`; `line_bases` sets wrapping.
    fn write_indexed_fasta(
        dir: &Path,
        records: &[(&str, &[u8])],
        line_bases: usize,
    ) -> std::path::PathBuf {
        let fasta_path = dir.join("ref.fa");
        let index_path = dir.join("ref.fa.fai");
        let mut fasta_file = File::create(&fasta_path).unwrap();
        let mut index_file = File::create(&index_path).unwrap();
        let mut offset: u64 = 0;
        let line_width = line_bases + 1; // single '\n' terminator
        for (name, seq) in records {
            let header = format!(">{name}\n");
            fasta_file.write_all(header.as_bytes()).unwrap();
            offset += header.len() as u64;
            writeln!(index_file, "{name}\t{}\t{offset}\t{line_bases}\t{line_width}", seq.len())
                .unwrap();
            let mut written = 0;
            for chunk in seq.chunks(line_bases) {
                fasta_file.write_all(chunk).unwrap();
                fasta_file.write_all(b"\n").unwrap();
                written += chunk.len() + 1;
            }
            offset += written as u64;
        }
        fasta_path
    }

    #[test]
    fn total_sequence_bytes_partial_and_full_last_line() {
        // 22 bp at 10 bp/line (width 11 incl. '\n'): two full lines + "AC" (no term).
        assert_eq!(total_sequence_bytes(22, 10, 11), 11 + 11 + 2);
        // Exactly two full lines: last line has no terminator.
        assert_eq!(total_sequence_bytes(20, 10, 11), 11 + 10);
        assert_eq!(total_sequence_bytes(0, 10, 11), 0);
    }

    #[test]
    fn load_contig_into_strips_newlines_across_wrapping() {
        let dir = tempfile::tempdir().unwrap();
        let seq1 = b"ACGTACGTNNNNacgtACGT"; // 20 bp, mixed case + Ns
        let seq2 = b"TTTTAAAACCCCGGGGA"; // 17 bp, partial last line at width 8
        let fa = write_indexed_fasta(dir.path(), &[("chr1", seq1), ("chr2", seq2)], 8);
        let mut fasta = Fasta::from_path(&fa).unwrap();

        let mut buf = Vec::new();
        fasta.load_contig_into("chr1", false, &mut buf).unwrap();
        assert_eq!(buf, seq1, "chr1 bases, newlines stripped, case preserved");

        // Reusing the same buffer for a differently-sized contig must be exact.
        fasta.load_contig_into("chr2", false, &mut buf).unwrap();
        assert_eq!(buf, seq2, "chr2 reuses the buffer, partial last line");

        // Uppercase option.
        fasta.load_contig_into("chr1", true, &mut buf).unwrap();
        assert_eq!(buf, b"ACGTACGTNNNNACGTACGT");
    }

    #[test]
    fn load_contig_into_agrees_with_fetch() {
        // Cross-check the bulk read against noodles' region query for the whole contig.
        let dir = tempfile::tempdir().unwrap();
        let seq = b"ACGTNNNNACGTACGTNNNNTTTTGGGGCCCCAAAA"; // 36 bp
        let fa = write_indexed_fasta(dir.path(), &[("c", seq)], 10);
        let mut fasta = Fasta::from_path(&fa).unwrap();

        let mut bulk = Vec::new();
        fasta.load_contig_into("c", false, &mut bulk).unwrap();
        let via_fetch = fasta.fetch("c", 0, seq.len() as u64).unwrap();
        assert_eq!(bulk, via_fetch, "bulk read must match noodles' region query");
    }

    #[test]
    fn load_contig_unknown_contig_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fa = write_indexed_fasta(dir.path(), &[("c", b"ACGT")], 4);
        let mut fasta = Fasta::from_path(&fa).unwrap();
        assert!(fasta.load_contig("nope", false).is_err());
    }

    #[test]
    fn load_contig_into_exact_multiple_of_line_bases_round_trips() {
        // 20 bp at 10 bp/line = exactly two full lines (no partial last line),
        // exercising the exact-multiple branch of the raw-byte-count geometry
        // through a real load_contig_into round-trip (not just the pure helper).
        let dir = tempfile::tempdir().unwrap();
        let seq = b"ACGTACGTACNNNNGGGGCC"; // 20 bp = 2 * 10
        let fa = write_indexed_fasta(dir.path(), &[("c", seq)], 10);
        let mut fasta = Fasta::from_path(&fa).unwrap();
        let mut buf = Vec::new();
        fasta.load_contig_into("c", false, &mut buf).unwrap();
        assert_eq!(buf, seq);
    }

    #[test]
    fn from_path_rejects_corrupt_fai_geometry() {
        // A .fai claiming line_width (4) < line_bases (10) would mis-stride the
        // bulk read, so from_path must reject it. Hand-write the .fai (the helper
        // only emits well-formed geometry).
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");
        std::fs::write(&fa, b">c\nACGTACGTAC\n").unwrap();
        // name \t length \t offset \t line_bases \t line_width
        std::fs::write(&fai, b"c\t10\t3\t10\t4\n").unwrap();
        // `Fasta` has no `Debug`, so bind the error via let-else rather than unwrap_err.
        let Err(err) = Fasta::from_path(&fa) else {
            panic!("expected a corrupt-geometry error");
        };
        assert!(
            err.to_string().contains("Corrupt .fai geometry"),
            "expected a corrupt-geometry error, got: {err}"
        );
    }
}
