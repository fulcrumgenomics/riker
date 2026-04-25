use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::{
    Flags, MappingQuality,
    cigar::{Op, op::Kind},
};
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
use serde::de::DeserializeOwned;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::NamedTempFile;

// ─── SamBuilder ─────────────────────────────────────────────────────────────

/// The order in which records are written to the BAM file.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SortOrder {
    /// Records are written in insertion order (default).
    #[default]
    Unsorted,
    /// Records are sorted by (`ref_id`, `alignment_start`) before writing.
    Coordinate,
    /// Records are sorted by read name before writing.
    QueryName,
}

/// Builds in-memory SAM records and writes them to a temporary BAM file.
/// All test data is constructed programmatically — no checked-in files.
pub struct SamBuilder {
    header: Header,
    records: Vec<RecordBuf>,
    sort_order: SortOrder,
}

impl SamBuilder {
    /// Create a builder with a single `chr1` contig of 249 Mbp.
    pub fn new() -> Self {
        Self::with_contigs(&[("chr1".to_string(), 249_250_621)])
    }

    /// Create a builder with the given list of `(name, length)` contigs.
    pub fn with_contigs(contigs: &[(String, usize)]) -> Self {
        let mut builder = Header::builder();
        for (name, len) in contigs {
            let ref_seq = Map::<ReferenceSequence>::new(
                NonZeroUsize::new(*len).expect("contig length must be > 0"),
            );
            builder = builder.add_reference_sequence(name.as_bytes(), ref_seq);
        }
        Self { header: builder.build(), records: Vec::new(), sort_order: SortOrder::Unsorted }
    }

    /// Set the sort order for the BAM file written by `to_temp_bam` / `write_to_file`.
    #[allow(dead_code)]
    pub fn sort_order(mut self, order: SortOrder) -> Self {
        self.sort_order = order;
        self
    }

    /// Add a properly paired FR read pair.
    ///
    /// `pos1` and `pos2` are 1-based alignment starts. `tlen` is the (positive)
    /// template length assigned to read1; read2 gets `-tlen`.
    #[allow(dead_code, clippy::too_many_arguments)]
    pub fn add_pair(
        &mut self,
        name: &str,
        ref_id: usize,
        pos1: usize,
        pos2: usize,
        tlen: i32,
        mapq: u8,
        read_len: usize,
        is_duplicate: bool,
        fails_vendor_quality: bool,
    ) -> &mut Self {
        let mut flags1 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::MATE_REVERSE_COMPLEMENTED
            | Flags::FIRST_SEGMENT;
        let mut flags2 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::REVERSE_COMPLEMENTED
            | Flags::LAST_SEGMENT;
        if is_duplicate {
            flags1 |= Flags::DUPLICATE;
            flags2 |= Flags::DUPLICATE;
        }
        if fails_vendor_quality {
            flags1 |= Flags::QC_FAIL;
            flags2 |= Flags::QC_FAIL;
        }

        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        let seq: Sequence = vec![b'A'; read_len].into();
        let qual = noodles::sam::alignment::record_buf::QualityScores::from(vec![30u8; read_len]);
        let pos1_p = Position::new(pos1).expect("valid pos1");
        let pos2_p = Position::new(pos2).expect("valid pos2");
        let mq_val = MappingQuality::new(mapq).expect("valid mapq");

        let r1 = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags1)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos1_p)
            .set_mapping_quality(mq_val)
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(ref_id)
            .set_mate_alignment_start(pos2_p)
            .set_template_length(tlen)
            .set_sequence(seq.clone())
            .set_quality_scores(qual.clone())
            .build();

        let r2 = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags2)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos2_p)
            .set_mapping_quality(mq_val)
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(ref_id)
            .set_mate_alignment_start(pos1_p)
            .set_template_length(-tlen)
            .set_sequence(seq)
            .set_quality_scores(qual)
            .build();

        self.records.push(r1);
        self.records.push(r2);
        self
    }

    /// Add a single, arbitrary record (for edge-case tests).
    #[allow(dead_code)]
    pub fn add_record(&mut self, record: RecordBuf) -> &mut Self {
        self.records.push(record);
        self
    }

    /// Add a mapped unpaired (fragment) read.
    ///
    /// `pos` is 1-based.  `nm` sets the NM auxiliary tag when `Some`.
    #[allow(dead_code, clippy::too_many_arguments)]
    pub fn add_unpaired(
        &mut self,
        name: &str,
        ref_id: usize,
        pos: usize,
        mapq: u8,
        read_len: usize,
        is_reverse: bool,
        is_duplicate: bool,
        fails_vendor_quality: bool,
        nm: Option<u8>,
    ) -> &mut Self {
        let mut flags = Flags::empty();
        if is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if is_duplicate {
            flags |= Flags::DUPLICATE;
        }
        if fails_vendor_quality {
            flags |= Flags::QC_FAIL;
        }

        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        let seq: Sequence = vec![b'A'; read_len].into();
        let qual = QualityScores::from(vec![30u8; read_len]);
        let pos_p = Position::new(pos).expect("valid pos");
        let mq_val = MappingQuality::new(mapq).expect("valid mapq");

        let mut data = noodles::sam::alignment::record_buf::Data::default();
        if let Some(n) = nm {
            data.insert((*b"NM").into(), DataValue::UInt8(n));
        }

        let record = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos_p)
            .set_mapping_quality(mq_val)
            .set_cigar(cigar)
            .set_sequence(seq)
            .set_quality_scores(qual)
            .set_data(data)
            .build();

        self.records.push(record);
        self
    }

    /// Add a paired read pair with an NM tag and configurable insert size / MAPQ.
    ///
    /// Returns `&mut Self` for chaining.
    #[allow(dead_code, clippy::too_many_arguments)]
    pub fn add_pair_with_nm(
        &mut self,
        name: &str,
        ref_id: usize,
        pos1: usize,
        pos2: usize,
        tlen: i32,
        mapq: u8,
        read_len: usize,
        nm1: u8,
        nm2: u8,
    ) -> &mut Self {
        let flags1 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::MATE_REVERSE_COMPLEMENTED
            | Flags::FIRST_SEGMENT;
        let flags2 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::REVERSE_COMPLEMENTED
            | Flags::LAST_SEGMENT;

        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        let seq: Sequence = vec![b'A'; read_len].into();
        let qual = QualityScores::from(vec![30u8; read_len]);
        let pos1_p = Position::new(pos1).expect("valid pos1");
        let pos2_p = Position::new(pos2).expect("valid pos2");
        let mq_val = MappingQuality::new(mapq).expect("valid mapq");

        let mut data1 = noodles::sam::alignment::record_buf::Data::default();
        data1.insert((*b"NM").into(), DataValue::UInt8(nm1));
        let mut data2 = noodles::sam::alignment::record_buf::Data::default();
        data2.insert((*b"NM").into(), DataValue::UInt8(nm2));

        let r1 = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags1)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos1_p)
            .set_mapping_quality(mq_val)
            .set_cigar(cigar.clone())
            .set_mate_reference_sequence_id(ref_id)
            .set_mate_alignment_start(pos2_p)
            .set_template_length(tlen)
            .set_sequence(seq.clone())
            .set_quality_scores(qual.clone())
            .set_data(data1)
            .build();

        let r2 = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags2)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos2_p)
            .set_mapping_quality(mq_val)
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(ref_id)
            .set_mate_alignment_start(pos1_p)
            .set_template_length(-tlen)
            .set_sequence(seq)
            .set_quality_scores(qual)
            .set_data(data2)
            .build();

        self.records.push(r1);
        self.records.push(r2);
        self
    }

    /// Add an unmapped read (for adapter-detection and UNPAIRED category tests).
    #[allow(dead_code)]
    pub fn add_unmapped(
        &mut self,
        name: &str,
        read_len: usize,
        fails_vendor_quality: bool,
        bases: Option<&[u8]>,
    ) -> &mut Self {
        let mut flags = Flags::UNMAPPED;
        if fails_vendor_quality {
            flags |= Flags::QC_FAIL;
        }

        let seq_bytes = bases.map_or_else(|| vec![b'A'; read_len], <[u8]>::to_vec);
        let actual_len = seq_bytes.len();
        let seq: Sequence = seq_bytes.into();
        let qual = QualityScores::from(vec![30u8; actual_len]);

        let record = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_sequence(seq)
            .set_quality_scores(qual)
            .build();

        self.records.push(record);
        self
    }

    /// Write all records to a BAM file at `path`, respecting the configured sort order.
    ///
    /// # Errors
    /// Returns an error if the BAM file cannot be created or a record cannot be written.
    pub fn write_to_file(&self, path: &Path) -> Result<()> {
        let mut records: Vec<&RecordBuf> = self.records.iter().collect();
        match self.sort_order {
            SortOrder::Unsorted => {}
            SortOrder::Coordinate => {
                records.sort_by_key(|r| {
                    let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
                    let pos = r.alignment_start().map_or(u64::MAX, |p| p.get() as u64);
                    (rid, pos)
                });
            }
            SortOrder::QueryName => {
                records.sort_by_key(|r| {
                    r.name().map(|n| -> Vec<u8> {
                        let bytes: &[u8] = n;
                        bytes.to_vec()
                    })
                });
            }
        }

        let file = std::fs::File::create(path)?;
        let mut writer = bam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&self.header)?;
        for record in records {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        Ok(())
    }

    /// Write all records to a CRAM file at `path`, using `fasta_path` as
    /// the reference. Reference sequences in the BAM header must match
    /// the sequences in the FASTA index for CRAM encoding to round-trip.
    /// Sort order is honoured the same way [`Self::write_to_file`] does.
    ///
    /// # Errors
    /// Returns an error if the CRAM file or its references cannot be
    /// resolved, or if record writing fails.
    pub fn write_cram_to_file(&self, path: &Path, fasta_path: &Path) -> Result<()> {
        use noodles::cram;
        use noodles::fasta;

        let mut records: Vec<&RecordBuf> = self.records.iter().collect();
        match self.sort_order {
            SortOrder::Unsorted => {}
            SortOrder::Coordinate => {
                records.sort_by_key(|r| {
                    let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
                    let pos = r.alignment_start().map_or(u64::MAX, |p| p.get() as u64);
                    (rid, pos)
                });
            }
            SortOrder::QueryName => {
                records.sort_by_key(|r| {
                    r.name().map(|n| -> Vec<u8> {
                        let bytes: &[u8] = n;
                        bytes.to_vec()
                    })
                });
            }
        }

        let fasta_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let repo =
            fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta_reader));

        let file = std::fs::File::create(path)?;
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_writer(file);
        writer.write_header(&self.header)?;
        for record in records {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        writer.try_finish(&self.header)?;
        Ok(())
    }

    /// Write to a temporary CRAM file and return the handle (deleted on
    /// drop). `fasta_path` is the indexed FASTA the CRAM encodes
    /// against.
    ///
    /// # Errors
    /// Returns an error if the temp file cannot be created or written.
    #[allow(dead_code)]
    pub fn to_temp_cram(&self, fasta_path: &Path) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".cram")?;
        self.write_cram_to_file(tmp.path(), fasta_path)?;
        Ok(tmp)
    }

    /// Write to a temporary BAM file and return the handle (deleted on drop).
    ///
    /// # Errors
    /// Returns an error if the temp file cannot be created or written.
    #[allow(dead_code)]
    pub fn to_temp_bam(&self) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".bam")?;
        self.write_to_file(tmp.path())?;
        Ok(tmp)
    }

    /// Write all records to a plain (uncompressed) SAM file at `path`.
    /// Sort order is honoured the same way [`Self::write_to_file`] does.
    ///
    /// # Errors
    /// Returns an error if the SAM file cannot be created or a record
    /// cannot be written.
    #[allow(dead_code)]
    pub fn write_sam_to_file(&self, path: &Path) -> Result<()> {
        use noodles::sam;

        let mut records: Vec<&RecordBuf> = self.records.iter().collect();
        match self.sort_order {
            SortOrder::Unsorted => {}
            SortOrder::Coordinate => {
                records.sort_by_key(|r| {
                    let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
                    let pos = r.alignment_start().map_or(u64::MAX, |p| p.get() as u64);
                    (rid, pos)
                });
            }
            SortOrder::QueryName => {
                records.sort_by_key(|r| {
                    r.name().map(|n| -> Vec<u8> {
                        let bytes: &[u8] = n;
                        bytes.to_vec()
                    })
                });
            }
        }

        let file = std::fs::File::create(path)?;
        let mut writer = sam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&self.header)?;
        for record in records {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        Ok(())
    }

    /// Write to a temporary `.sam` file and return the handle. Used by
    /// integration tests that exercise the SAM in-place fill path.
    ///
    /// # Errors
    /// Returns an error if the temp file cannot be created or written.
    #[allow(dead_code)]
    pub fn to_temp_sam(&self) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".sam")?;
        self.write_sam_to_file(tmp.path())?;
        Ok(tmp)
    }

    /// Write to a temporary BAM file with a BAI index and return the handle.
    ///
    /// The BAI index is written to `<path>.bai`. This is needed for commands
    /// that use indexed BAM reading (e.g. `error`).
    ///
    /// # Errors
    /// Returns an error if the BAM or index cannot be written.
    #[allow(dead_code)]
    pub fn to_temp_indexed_bam(&self) -> Result<NamedTempFile> {
        use noodles::bam::{self as bam_crate, bai};
        use noodles::csi::binning_index::Indexer;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::sam::alignment::Record as _;

        let tmp = NamedTempFile::with_suffix(".bam")?;
        self.write_to_file(tmp.path())?;

        // Build the BAI index by reading through the BAM
        let mut reader = bam_crate::io::Reader::new(std::fs::File::open(tmp.path())?);
        let header = reader.read_header()?;

        let mut indexer = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();
        let mut record = bam_crate::Record::default();

        while reader.read_record(&mut record)? != 0 {
            let chunk_end = reader.get_ref().virtual_position();

            let alignment_context = match (
                record.reference_sequence_id().transpose()?,
                record.alignment_start().transpose()?,
                record.alignment_end().transpose()?,
            ) {
                (Some(id), Some(start), Some(end)) => {
                    let is_mapped = !record.flags().is_unmapped();
                    Some((id, start, end, is_mapped))
                }
                _ => None,
            };

            let chunk = Chunk::new(chunk_start, chunk_end);
            indexer.add_record(alignment_context, chunk)?;
            chunk_start = chunk_end;
        }

        let index = indexer.build(header.reference_sequences().len());

        // Write the BAI index
        let bai_path = {
            let mut p = tmp.path().as_os_str().to_owned();
            p.push(".bai");
            std::path::PathBuf::from(p)
        };
        bai::fs::write(bai_path, &index)?;

        Ok(tmp)
    }

    #[allow(dead_code)]
    pub fn header(&self) -> &Header {
        &self.header
    }
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ─── FastaBuilder ────────────────────────────────────────────────────────

/// Builds an in-memory FASTA reference and writes it to a temporary file with index.
#[allow(dead_code)]
pub struct FastaBuilder {
    contigs: Vec<(String, Vec<u8>)>,
}

#[allow(dead_code)]
impl FastaBuilder {
    pub fn new() -> Self {
        Self { contigs: Vec::new() }
    }

    /// Add a contig with the given name and sequence (ASCII bases).
    pub fn add_contig(mut self, name: &str, sequence: &[u8]) -> Self {
        self.contigs.push((name.to_string(), sequence.to_vec()));
        self
    }

    /// Write FASTA + .fai index to a temporary file and return the handle.
    ///
    /// The .fai index is written to `<path>.fai` alongside the FASTA.
    ///
    /// # Errors
    /// Returns an error if the temp file cannot be created or written.
    pub fn to_temp_fasta(&self) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".fa")?;
        self.write_fasta(tmp.path())?;
        Ok(tmp)
    }

    fn write_fasta(&self, path: &Path) -> Result<()> {
        use std::io::Write;
        let mut fai_path = path.to_path_buf();
        fai_path.set_extension("fa.fai");

        let mut fasta = std::fs::File::create(path)?;
        let mut fai = std::fs::File::create(&fai_path)?;

        let line_bases: u64 = 60;
        let mut offset: u64 = 0;

        for (name, seq) in &self.contigs {
            let header_line = format!(">{name}\n");
            fasta.write_all(header_line.as_bytes())?;
            let header_len = header_line.len() as u64;

            let seq_start = offset + header_len;
            let seq_len = seq.len() as u64;

            // Write sequence in lines of 60 bases
            #[allow(clippy::cast_possible_truncation)]
            for chunk in seq.chunks(line_bases as usize) {
                fasta.write_all(chunk)?;
                fasta.write_all(b"\n")?;
            }

            // FAI: name, length, offset, bases_per_line, bytes_per_line
            let bytes_per_line = line_bases + 1; // +1 for newline
            writeln!(fai, "{name}\t{seq_len}\t{seq_start}\t{line_bases}\t{bytes_per_line}")?;

            // Update offset: header + full lines + partial last line
            let full_lines = seq_len / line_bases;
            let remainder = seq_len % line_bases;
            let seq_bytes =
                full_lines * bytes_per_line + if remainder > 0 { remainder + 1 } else { 0 };
            offset = seq_start + seq_bytes;
        }

        Ok(())
    }
}

#[allow(dead_code)]
impl Default for FastaBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Metric reading ──────────────────────────────────────────────────────────

/// Read a tab-separated metrics file into a vector of typed rows.
///
/// # Errors
/// Returns an error if the file cannot be read or deserialized.
pub fn read_metrics_tsv<T: DeserializeOwned>(path: &Path) -> Result<Vec<T>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let rows = rdr.deserialize().collect::<Result<Vec<T>, _>>()?;
    Ok(rows)
}

// ─── Shared test builders ────────────────────────────────────────────────────

/// Create a [`SamBuilder`] with the given contigs, always sorted by coordinate.
#[allow(dead_code)]
pub fn coord_builder(contigs: &[(&str, usize)]) -> SamBuilder {
    let contigs: Vec<(String, usize)> = contigs.iter().map(|(n, l)| (n.to_string(), *l)).collect();
    SamBuilder::with_contigs(&contigs).sort_order(SortOrder::Coordinate)
}

// ─── Float assertions ────────────────────────────────────────────────────────

/// Assert that two `f64` values are within `epsilon` of each other.
#[macro_export]
macro_rules! assert_float_eq {
    ($a:expr, $b:expr, $eps:expr) => {{
        let a: f64 = $a;
        let b: f64 = $b;
        let eps: f64 = $eps;
        assert!(
            (a - b).abs() <= eps,
            "assertion failed: |{} - {}| = {} > {}",
            a,
            b,
            (a - b).abs(),
            eps
        );
    }};
}
