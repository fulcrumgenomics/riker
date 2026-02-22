use std::fs::File;
use std::io::{Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use noodles::core::Region;
use noodles::fasta;
use noodles::sam::Header;

use crate::sequence_dict::SequenceDictionary;

/// A loaded FASTA reference with random-access sequence retrieval.
///
/// Wraps a noodles `IndexedReader` for efficient seek-based access.  Contig
/// metadata (names, lengths) is pre-computed from the `.fai` index at
/// construction time.
pub struct Fasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    dict: SequenceDictionary,
}

impl Fasta {
    /// Open a FASTA file and its `.fai` index.
    ///
    /// # Errors
    /// Returns an error if the FASTA or its index cannot be read.
    pub fn from_path(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed FASTA: {}", path.display()))?;

        let dict = SequenceDictionary::from(reader.index());

        Ok(Self { reader, dict })
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

    /// Load and return the full sequence of `contig_name`.
    ///
    /// When `uppercase` is true the returned bytes are converted to uppercase
    /// ASCII in place.
    ///
    /// # Errors
    /// Returns an error if the contig is unknown, or if the FASTA cannot be read.
    pub fn load_contig(&mut self, contig_name: &str, uppercase: bool) -> Result<Vec<u8>> {
        let meta = self
            .dict
            .get_by_name(contig_name)
            .ok_or_else(|| anyhow!("Contig '{contig_name}' not found in reference"))?;
        let expected_len = meta.length();

        // Seek directly to the sequence start using the FAI index, then read
        // into our buffer.  This avoids the intermediate Record allocation that
        // `query()` performs — halving peak memory for large contigs.
        let region = Region::new(contig_name, ..);
        let pos =
            self.reader.index().query(&region).with_context(|| {
                format!("Failed to look up contig '{contig_name}' in FASTA index")
            })?;
        self.reader
            .get_mut()
            .seek(SeekFrom::Start(pos))
            .with_context(|| format!("Failed to seek to contig '{contig_name}' in FASTA"))?;

        let mut seq = Vec::with_capacity(expected_len);
        self.reader
            .read_sequence(&mut seq)
            .with_context(|| format!("Failed to read contig '{contig_name}' from FASTA"))?;

        if uppercase {
            for b in &mut seq {
                *b = b.to_ascii_uppercase();
            }
        }
        Ok(seq)
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
