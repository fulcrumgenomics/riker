//! The SAM record builders and the [`SamBuilder`] that gathers them into a file.
//!
//! Ports the ergonomics of fgbio's `SamBuilder` to Rust — name only the fields that
//! differ from a sensible default, via chained setters, instead of 8–10 positional
//! args or the raw `RecordBuf::builder()` blocks the review found duplicated across
//! the alignment, hybcap, error, wgs, and gcbias tests:
//!
//! - [`read`] builds one record (a fragment); [`pair`] a mate pair, deriving mate
//!   fields, `MC`, and template length. Pair geometry comes from [`PairBuilder::at`]
//!   (two read positions) or [`PairBuilder::fr`] (an FR pair by fragment span, which
//!   models adapter read-through when the fragment is shorter than the reads).
//! - [`SamBuilder`] gathers records — `sam.add(read()/pair())` — and writes a temp
//!   SAM/BAM/CRAM. It hands out template names via [`SamBuilder::next_id`] and holds
//!   the read-through [`Adapters`].
//!
//! ```ignore
//! let mut sam = SamBuilder::new();
//! sam.add(pair(sam.next_id()).fr("chr1", 100..400))
//!    .add(read().at("chr1", 900).cigar("10S80M10S").reverse());
//! let bam = sam.to_temp_bam()?;
//! ```

// Builder entry points and setters return the builder or `Self`, so allow the two
// "could be `#[must_use]`" lints rather than annotate every method. `similar_names`
// is allowed because read-pair code is naturally full of `r1`/`r2`-style names.
#![allow(clippy::must_use_candidate, clippy::return_self_not_must_use, clippy::similar_names)]

use std::cell::Cell;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::ops::Range;
use std::path::Path;

use anyhow::Result;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;
use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::sam::header::record::value::{
    Map, map::Header as HeaderMap, map::ReadGroup, map::ReferenceSequence,
};
use tempfile::NamedTempFile;

use crate::sam::riker_record::RikerRecord;

use super::cigar;
use super::fasta_builder::{self, FastaBuilder};

/// Read length used when neither an explicit CIGAR nor explicit bases pin it.
const DEFAULT_READ_LEN: usize = 100;
/// Default mapping quality — a confidently-mapped read, matching fgbio's default.
const DEFAULT_MAPQ: u8 = 60;
/// Default per-base quality (Phred), matching the previous `SamBuilder` helpers.
const DEFAULT_QUAL: u8 = 30;

/// Illumina TruSeq 3' adapters — the sequence a read runs into when its fragment is
/// shorter than the read. Keyed to read *number*, not genomic strand: the first read
/// of a pair (R1) reads into `TRUSEQ_R1` and the second (R2) into `TRUSEQ_R2`,
/// whichever strand each maps to. (Strand only decides how the read-through is laid
/// into the stored `SEQ`: a forward read's at the 3' end as-is, a reverse read's at
/// the 5' end reverse-complemented.) The default read-through adapters (see
/// [`Adapters`]); byte-for-byte the TruSeq entry in chelae's `adapter_db`.
const TRUSEQ_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const TRUSEQ_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// How a read's `NM` (edit-distance) tag is resolved at build time.
enum NmPolicy {
    /// Derive NM from the read's bases versus the reference (mapped reads only).
    Auto,
    /// Emit no NM tag.
    None,
    /// Emit an explicit NM value.
    Explicit(u8),
}

// ─── SamBuilder ────────────────────────────────────────────────────────────────

/// The order records are written in (and stamped into the `@HD SO` tag).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SortOrder {
    /// Insertion order (default).
    #[default]
    Unsorted,
    /// Sorted by `(reference_sequence_id, alignment_start)`.
    Coordinate,
}

/// Accumulates `read()`/`pair()` records and writes them to a temporary SAM/BAM/CRAM
/// — the top-level entry point for building a test alignment file. Pick contigs, then
/// `add` reads and pairs, then `to_temp_bam()` (or a SAM/CRAM/indexed variant).
///
/// ```ignore
/// let mut sam = SamBuilder::new();
/// sam.add(pair(sam.next_id()).at("chr1", 100, 300))
///    .add(read().at("chr1", 500).reverse());
/// let bam = sam.to_temp_bam()?;
/// ```
pub struct SamBuilder {
    header: Header,
    records: Vec<RecordBuf>,
    sort_order: SortOrder,
    // Running counter + optional prefix behind `next_id`; a `Cell` so `next_id` takes
    // `&self` and can be called inline in `add(pair(sam.next_id())..)` without fighting
    // `add`'s `&mut self` borrow.
    counter: Cell<usize>,
    id_prefix: Option<String>,
    /// The read-through adapters a `pair().fr()` short fragment sequences into
    /// (Illumina TruSeq by default; override via [`Self::adapter_r1`] / [`Self::adapter_r2`]).
    adapters: Adapters,
}

impl SamBuilder {
    /// A builder with a single `chr1` contig of 249 Mbp — long enough that no test
    /// outgrows it. Note this length is only the header's; a mapped `read()` derives its
    /// *bases* from the frozen default reference (see [`fasta_builder::hg38_mini`]), whose
    /// `chr1` is ~10 kb, so a read placed past that pads with `A` (harmless where base
    /// content is unread, or pin it with `.bases()`/`.matching_ref()`).
    pub fn new() -> Self {
        Self::with_contigs(&[("chr1".to_string(), 249_250_621)])
    }

    /// A builder with the given `(name, length)` contigs. Their order sets each contig's
    /// `reference_sequence_id`; records added via [`Self::add`] resolve their contig names
    /// against this header, so any order works.
    ///
    /// # Panics
    /// Panics if any contig length is 0.
    pub fn with_contigs(contigs: &[(String, usize)]) -> Self {
        let mut header = Header::builder();
        for (name, length) in contigs {
            let reference = Map::<ReferenceSequence>::new(
                NonZeroUsize::new(*length).expect("contig length must be > 0"),
            );
            header = header.add_reference_sequence(name.as_bytes(), reference);
        }
        Self {
            header: header.build(),
            records: Vec::new(),
            sort_order: SortOrder::Unsorted,
            counter: Cell::new(0),
            id_prefix: None,
            adapters: Adapters::default(),
        }
    }

    /// Set the write order, stamping the matching `@HD SO` tag so tools that require a
    /// declared sort order (e.g. `hybcap`, `rna`) accept the file.
    pub fn sort_order(mut self, order: SortOrder) -> Self {
        self.sort_order = order;
        let value: &[u8] = match order {
            SortOrder::Unsorted => b"unsorted",
            SortOrder::Coordinate => b"coordinate",
        };
        self.stamp_sort_order(value);
        self
    }

    /// Stamp `@HD SO:coordinate` WITHOUT reordering records — a deliberately mislabeled
    /// BAM that exercises a tool's runtime order guard.
    pub fn declare_coordinate_sorted(mut self) -> Self {
        self.stamp_sort_order(b"coordinate");
        self
    }

    /// Add an `@RG` read-group line with the given ID and `SM` (sample) tag, so tests
    /// that read a sample name back from the header (e.g. `derive_sample`) can build one.
    pub fn read_group(mut self, id: impl Into<String>, sample: impl Into<String>) -> Self {
        use noodles::sam::header::record::value::map::read_group::tag;
        let sample: String = sample.into();
        let mut read_group = Map::<ReadGroup>::default();
        read_group.other_fields_mut().insert(tag::SAMPLE, sample.into());
        let id: String = id.into();
        self.header.read_groups_mut().insert(id.into(), read_group);
        self
    }

    /// Add an `@RG` read-group line with only an ID and no `SM` tag — for exercising the
    /// no-sample fallback path.
    pub fn read_group_without_sample(mut self, id: impl Into<String>) -> Self {
        let id: String = id.into();
        self.header.read_groups_mut().insert(id.into(), Map::<ReadGroup>::default());
        self
    }

    /// Set the prefix [`Self::next_id`] prepends to its counter (default: none, so ids
    /// are bare numbers — `"1"`, `"2"`, …; `id_prefix("q")` gives `"q1"`, `"q2"`, …).
    pub fn id_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.id_prefix = Some(prefix.into());
        self
    }

    /// Override the read-through adapter for read 1 (default: Illumina TruSeq). A
    /// `pair().fr()` fragment shorter than its reads sequences its 3' overhang into this
    /// adapter, so a test can pin an unambiguous, self-chosen read-through. See [`Adapters`].
    pub fn adapter_r1(mut self, adapter: impl Into<Vec<u8>>) -> Self {
        self.adapters.r1 = adapter.into();
        self
    }

    /// Override the read-through adapter for read 2 (default: Illumina TruSeq); the read-2
    /// counterpart of [`Self::adapter_r1`].
    pub fn adapter_r2(mut self, adapter: impl Into<Vec<u8>>) -> Self {
        self.adapters.r2 = adapter.into();
        self
    }

    /// The next unique read/template name: the running counter, prefixed per
    /// [`Self::id_prefix`]. Lets a test use `sam.next_id()` for a made-for-you name
    /// instead of hand-rolling `format!("read{i}")`.
    pub fn next_id(&self) -> String {
        self.counter.set(self.counter.get() + 1);
        let n = self.counter.get();
        match &self.id_prefix {
            Some(prefix) if !prefix.is_empty() => format!("{prefix}{n}"),
            _ => n.to_string(),
        }
    }

    /// Add a `read()`/`pair()` builder or a raw [`RecordBuf`]. A `pair().fr()` read-through
    /// fills from this builder's [`Adapters`] (TruSeq by default; override via
    /// [`Self::adapter_r1`] / [`Self::adapter_r2`]).
    ///
    /// Contig names are resolved to `reference_sequence_id`s against *this builder's own
    /// header*, so a record always lands on the contig its name names, whatever order the
    /// header declares its contigs in. (A raw [`RecordBuf`] is added as-is; it carries its
    /// own ids.)
    ///
    /// # Panics
    /// Panics if a `read()`/`pair()` names a contig this builder's header does not declare.
    pub fn add(&mut self, records: impl IntoRecords) -> &mut Self {
        let built = {
            let resolve = |name: &str| contig_id_in_header(&self.header, name);
            records.into_records(&self.adapters, &resolve)
        };
        self.records.extend(built);
        self
    }

    /// The header, for tests that also drive a collector directly.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Write to a temp `.bam` and return the handle (deletes it on drop).
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or a record cannot be written.
    pub fn to_temp_bam(&self) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".bam")?;
        self.write_bam(tmp.path())?;
        Ok(tmp)
    }

    /// Write to a temp `.sam` and return the handle.
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or a record cannot be written.
    pub fn to_temp_sam(&self) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".sam")?;
        self.write_sam(tmp.path())?;
        Ok(tmp)
    }

    /// Write to a temp `.cram` encoded against `fasta_path` and return the handle. The
    /// header's reference sequences must match the FASTA for CRAM to round-trip.
    ///
    /// # Errors
    /// Returns an error if the references cannot be resolved or writing fails.
    pub fn to_temp_cram(&self, fasta_path: &Path) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".cram")?;
        self.write_cram(tmp.path(), fasta_path)?;
        Ok(tmp)
    }

    /// Write to a temp `.bam` with a sibling `.bam.bai` index (for indexed reads, e.g. `error`).
    ///
    /// # Errors
    /// Returns an error if the BAM or its index cannot be written.
    pub fn to_temp_indexed_bam(&self) -> Result<NamedTempFile> {
        use noodles::bam::bai;
        let tmp = NamedTempFile::with_suffix(".bam")?;
        self.write_bam(tmp.path())?;
        let index: bai::Index = build_bam_index(tmp.path())?;
        bai::fs::write(with_suffix(tmp.path(), ".bai"), &index)?;
        Ok(tmp)
    }

    /// Write to a temp `.bam` with a sibling `.bam.csi` index (the large-contig format).
    ///
    /// # Errors
    /// Returns an error if the BAM or its index cannot be written.
    pub fn to_temp_csi_indexed_bam(&self) -> Result<NamedTempFile> {
        use noodles::csi;
        let tmp = NamedTempFile::with_suffix(".bam")?;
        self.write_bam(tmp.path())?;
        let index: csi::Index = build_bam_index(tmp.path())?;
        csi::fs::write(with_suffix(tmp.path(), ".csi"), &index)?;
        Ok(tmp)
    }

    /// Write all records to a BAM at `path`, honoring the sort order.
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or a record cannot be written.
    pub fn write_bam(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = bam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&self.header)?;
        for record in self.sorted_records() {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        Ok(())
    }

    /// Write all records to a plain SAM at `path`, honoring the sort order.
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or a record cannot be written.
    pub fn write_sam(&self, path: &Path) -> Result<()> {
        use noodles::sam;
        let file = std::fs::File::create(path)?;
        let mut writer = sam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&self.header)?;
        for record in self.sorted_records() {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        Ok(())
    }

    /// Write all records to a CRAM at `path` encoded against `fasta_path`.
    ///
    /// # Errors
    /// Returns an error if the references cannot be resolved or a record cannot be written.
    pub fn write_cram(&self, path: &Path, fasta_path: &Path) -> Result<()> {
        use noodles::cram;
        use noodles::fasta;

        let reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let repository =
            fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader));
        let file = std::fs::File::create(path)?;
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_writer(file);
        writer.write_header(&self.header)?;
        for record in self.sorted_records() {
            AlignmentWrite::write_alignment_record(&mut writer, &self.header, record)?;
        }
        writer.try_finish(&self.header)?;
        Ok(())
    }

    /// Set the `@HD SO` value in place, creating the `@HD` map if absent. Mutating in
    /// place keeps `@RG`/`@PG`/`@CO` (which live outside the `@HD` map) intact.
    fn stamp_sort_order(&mut self, value: &[u8]) {
        self.header
            .header_mut()
            .get_or_insert_with(|| Map::<HeaderMap>::builder().build().expect("valid @HD header"))
            .other_fields_mut()
            .insert(SORT_ORDER, value.to_vec().into());
    }

    /// The records in write order per [`Self::sort_order`].
    fn sorted_records(&self) -> Vec<&RecordBuf> {
        let mut records: Vec<&RecordBuf> = self.records.iter().collect();
        match self.sort_order {
            SortOrder::Unsorted => {}
            SortOrder::Coordinate => records.sort_by_key(|record| {
                let id = record.reference_sequence_id().unwrap_or(usize::MAX);
                let start = record.alignment_start().map_or(u64::MAX, |p| p.get() as u64);
                (id, start)
            }),
        }
        records
    }
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// A coordinate-sorted [`SamBuilder`] over `contigs` — the shorthand for tools that
/// require a sorted input (`hybcap`, `rna`, `wgs`). Takes `&str` names for brevity.
pub fn coord_builder(contigs: &[(&str, usize)]) -> SamBuilder {
    let owned: Vec<(String, usize)> =
        contigs.iter().map(|&(name, length)| (name.to_string(), length)).collect();
    SamBuilder::with_contigs(&owned).sort_order(SortOrder::Coordinate)
}

/// The read-through adapters a short read sequences into, keyed to read *number*: `r1`
/// for the first read of a pair, `r2` for the second, whichever strand each maps to.
/// Default to Illumina TruSeq; override per-SamBuilder via
/// [`SamBuilder::adapter_r1`] / [`SamBuilder::adapter_r2`].
#[derive(Clone)]
pub struct Adapters {
    /// The adapter read 1 reads into.
    pub r1: Vec<u8>,
    /// The adapter read 2 reads into.
    pub r2: Vec<u8>,
}

impl Default for Adapters {
    fn default() -> Self {
        Self { r1: TRUSEQ_R1.to_vec(), r2: TRUSEQ_R2.to_vec() }
    }
}

/// What [`SamBuilder::add`] accepts: a fluent `read()`/`pair()` builder (finished on the
/// way in, so a pair's [`PairBuilder::fr`] read-through picks up the SamBuilder's
/// adapters) or an already-built [`RecordBuf`].
pub trait IntoRecords {
    /// Consume into the records to add, finishing any pair with `adapters` and resolving
    /// each mate's contig name to a `reference_sequence_id` via `resolve` (the owning
    /// [`SamBuilder`]'s header).
    fn into_records(self, adapters: &Adapters, resolve: &dyn Fn(&str) -> usize) -> Vec<RecordBuf>;
}

impl IntoRecords for RecordBuf {
    fn into_records(
        self,
        _adapters: &Adapters,
        _resolve: &dyn Fn(&str) -> usize,
    ) -> Vec<RecordBuf> {
        vec![self]
    }
}

impl IntoRecords for ReadBuilder {
    fn into_records(self, _adapters: &Adapters, resolve: &dyn Fn(&str) -> usize) -> Vec<RecordBuf> {
        vec![self.build_with(resolve)]
    }
}

impl IntoRecords for PairBuilder {
    fn into_records(self, adapters: &Adapters, resolve: &dyn Fn(&str) -> usize) -> Vec<RecordBuf> {
        let (r1, r2) = self.build_with_adapters(adapters, resolve);
        vec![r1, r2]
    }
}

/// Build a binning index for `bam_path` by reading through it. Generic over the
/// reference-sequence index type so it serves both BAI (`LinearIndex`) and CSI
/// (`BinnedIndex`); the caller's type annotation picks which.
fn build_bam_index<I>(bam_path: &Path) -> Result<noodles::csi::binning_index::Index<I>>
where
    I: noodles::csi::binning_index::index::reference_sequence::Index + Default,
{
    use noodles::csi::binning_index::Indexer;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::sam::alignment::Record as _;

    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path)?);
    let header = reader.read_header()?;
    let mut indexer: Indexer<I> = Indexer::default();
    let mut chunk_start = reader.get_ref().virtual_position();
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        let chunk_end = reader.get_ref().virtual_position();
        let context = match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(id), Some(start), Some(end)) => {
                Some((id, start, end, !record.flags().is_unmapped()))
            }
            _ => None,
        };
        indexer.add_record(context, Chunk::new(chunk_start, chunk_end))?;
        chunk_start = chunk_end;
    }
    Ok(indexer.build(header.reference_sequences().len()))
}

/// `path` with `suffix` appended after its full name (e.g. `x.bam` + `.bai` → `x.bam.bai`).
fn with_suffix(path: &Path, suffix: &str) -> std::path::PathBuf {
    let mut os = path.as_os_str().to_owned();
    os.push(suffix);
    std::path::PathBuf::from(os)
}

/// Start a fluent builder for one SAM record. See the [module docs](self).
#[must_use]
pub fn read() -> ReadBuilder {
    ReadBuilder::new()
}

/// Fluent builder for a single [`RecordBuf`]. Construct via [`read`].
///
/// Every field has a sensible default (a mapped, forward, primary read on `chr1`
/// at position 1, 100 bp of `A` at quality 30, MAPQ 60), so a test names only
/// what it is actually testing.
pub struct ReadBuilder {
    name: String,
    contig: String,
    start: usize,
    mapped: bool,
    mapq: u8,
    read_len: usize,
    cigar: Option<String>,
    bases: Option<Vec<u8>>,
    quals: Option<Vec<u8>>,
    uniform_qual: u8,
    flags: Flags,
    tags: Vec<([u8; 2], DataValue)>,
    nm: NmPolicy,
    // Reference to derive bases from (None = the frozen default) and per-offset base
    // substitutions applied after derivation.
    matching_ref: Option<FastaBuilder>,
    subs: Vec<(usize, u8)>,
    // Mate fields, populated only by `PairBuilder`; a fragment read leaves them unset.
    mate_ref_id: Option<usize>,
    mate_start: Option<usize>,
    template_length: i32,
    // Bases that fill the read's soft-clip op instead of the default `A`. Set by
    // `pair().fr(..)` to lay adapter read-through into a short fragment's overhang.
    readthrough: Option<Vec<u8>>,
}

impl ReadBuilder {
    fn new() -> Self {
        Self {
            name: "read".to_string(),
            contig: "chr1".to_string(),
            start: 1,
            mapped: true,
            mapq: DEFAULT_MAPQ,
            read_len: DEFAULT_READ_LEN,
            cigar: None,
            bases: None,
            quals: None,
            uniform_qual: DEFAULT_QUAL,
            flags: Flags::empty(),
            tags: Vec::new(),
            nm: NmPolicy::Auto,
            matching_ref: None,
            subs: Vec::new(),
            mate_ref_id: None,
            mate_start: None,
            template_length: 0,
            readthrough: None,
        }
    }

    /// Set the read name (default `"read"`).
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Place the read as mapped at `contig` (e.g. `"chr1"`) and `start` (1-based).
    ///
    /// The contig name is resolved to a `reference_sequence_id` against the default
    /// reference (see [`fasta_builder`](super::fasta_builder)); naming a contig
    /// outside that set panics at build time.
    pub fn at(mut self, contig: impl Into<String>, start: usize) -> Self {
        self.contig = contig.into();
        self.start = start;
        self.mapped = true;
        self
    }

    /// Set the CIGAR from a string such as `"10S80M10S"`. Sequence and qualities
    /// are sized to the CIGAR's query length unless set explicitly.
    pub fn cigar(mut self, cigar: impl Into<String>) -> Self {
        self.cigar = Some(cigar.into());
        self
    }

    /// Set the read length used when neither a CIGAR nor bases pin it (default 100).
    pub fn len(mut self, len: usize) -> Self {
        self.read_len = len;
        self
    }

    /// Set the mapping quality (default 60).
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    /// Set the read bases explicitly (default: all `A` of the query length).
    pub fn bases(mut self, bases: impl Into<Vec<u8>>) -> Self {
        self.bases = Some(bases.into());
        self
    }

    /// Set per-base qualities explicitly (default: uniform, see [`Self::qual`]).
    pub fn quals(mut self, quals: impl Into<Vec<u8>>) -> Self {
        self.quals = Some(quals.into());
        self
    }

    /// Set a single uniform per-base quality (default 30).
    pub fn qual(mut self, qual: u8) -> Self {
        self.uniform_qual = qual;
        self
    }

    /// Derive the read's bases from `reference` instead of the frozen default.
    pub fn matching_ref(mut self, reference: &FastaBuilder) -> Self {
        self.matching_ref = Some(reference.clone());
        self
    }

    /// Substitute `base` at 0-based read offset `offset` after base derivation —
    /// the way to plant a mismatch against the reference.
    pub fn sub(mut self, offset: usize, base: u8) -> Self {
        self.subs.push((offset, base));
        self
    }

    /// Mark the read as reverse-strand.
    pub fn reverse(mut self) -> Self {
        self.flags |= Flags::REVERSE_COMPLEMENTED;
        self
    }

    /// Mark the read as a secondary alignment.
    pub fn secondary(mut self) -> Self {
        self.flags |= Flags::SECONDARY;
        self
    }

    /// Mark the read as a supplementary (chimeric) alignment.
    pub fn supplementary(mut self) -> Self {
        self.flags |= Flags::SUPPLEMENTARY;
        self
    }

    /// Mark the read as a PCR/optical duplicate.
    pub fn duplicate(mut self) -> Self {
        self.flags |= Flags::DUPLICATE;
        self
    }

    /// Mark the read as failing vendor quality (clears the PF bit).
    pub fn qc_fail(mut self) -> Self {
        self.flags |= Flags::QC_FAIL;
        self
    }

    /// Make the read unmapped: it gets no reference, position, MAPQ, or CIGAR.
    pub fn unmapped(mut self) -> Self {
        self.mapped = false;
        self
    }

    /// Set an explicit `NM` (edit-distance) tag, overriding the derived value.
    pub fn nm(mut self, nm: u8) -> Self {
        self.nm = NmPolicy::Explicit(nm);
        self
    }

    /// Suppress the `NM` tag otherwise derived on a mapped read.
    pub fn no_nm(mut self) -> Self {
        self.nm = NmPolicy::None;
        self
    }

    /// Attach an `MC` (mate CIGAR) tag.
    pub fn mc(mut self, cigar: impl Into<String>) -> Self {
        let cigar: String = cigar.into();
        self.tags.push((*b"MC", DataValue::String(cigar.into_bytes().into())));
        self
    }

    /// Attach an arbitrary two-character aux tag.
    pub fn tag(mut self, tag: [u8; 2], value: DataValue) -> Self {
        self.tags.push((tag, value));
        self
    }

    /// Build the [`RecordBuf`], resolving the contig name against the frozen default
    /// table (see [`fasta_builder::DEFAULT_CONTIGS`]). A record added through a
    /// [`SamBuilder`] instead resolves its `reference_sequence_id` against that
    /// builder's own header (see [`SamBuilder::add`]).
    ///
    /// # Panics
    /// Panics if explicit bases or qualities disagree with the CIGAR's query
    /// length, if a mapped read is given a start of 0, or if MAPQ is 255.
    #[must_use]
    pub fn build(self) -> RecordBuf {
        self.build_with(&super::fasta_builder::default_contig_id)
    }

    /// Build the record and adapt it into a [`RikerRecord`] against an empty header, for
    /// tests that drive the internal per-record kernels (CIGAR walks, accumulators) directly
    /// rather than through a file. The empty header mirrors how these unit tests exercise a
    /// record in isolation, independent of any reference-sequence dictionary.
    ///
    /// # Panics
    /// Panics for the same reasons as [`Self::build`], or if the built record cannot be
    /// adapted into a [`RikerRecord`].
    #[must_use]
    pub fn into_riker_record(self) -> RikerRecord {
        RikerRecord::from_alignment_record(&Header::default(), &self.build()).unwrap()
    }

    /// Build the record, resolving the contig name to a `reference_sequence_id` via
    /// `resolve` (the frozen default table for a standalone [`Self::build`], or the
    /// owning [`SamBuilder`]'s header when added through it).
    fn build_with(self, resolve: &dyn Fn(&str) -> usize) -> RecordBuf {
        assert!(
            !self.mapped || self.start >= 1,
            "a mapped read's start must be >= 1, got {}",
            self.start
        );
        let parsed = if let Some(spec) = &self.cigar {
            cigar::parse(spec)
        } else {
            let len = if let Some(bases) = &self.bases { bases.len() } else { self.read_len };
            cigar::parse(&format!("{len}M"))
        };
        let query_len = parsed.query_len;

        let reference_span = if self.mapped {
            let default = fasta_builder::default_reference();
            let reference = self.matching_ref.as_ref().unwrap_or(default);
            Some(reference.bases(&self.contig, self.start, parsed.ref_len))
        } else {
            None
        };

        let mut bases = if let Some(bases) = self.bases {
            bases
        } else if let Some(reference_span) = &reference_span {
            derive_bases(reference_span, &parsed.ops, query_len, self.readthrough.as_deref())
        } else {
            vec![b'A'; query_len]
        };
        for &(offset, base) in &self.subs {
            assert!(offset < bases.len(), "sub offset {offset} past read length {}", bases.len());
            bases[offset] = base;
        }
        let quals = self.quals.unwrap_or_else(|| vec![self.uniform_qual; query_len]);
        assert_eq!(bases.len(), query_len, "bases length disagrees with CIGAR query length");
        assert_eq!(quals.len(), query_len, "qualities length disagrees with CIGAR query length");

        let nm = match self.nm {
            NmPolicy::None => None,
            NmPolicy::Explicit(value) => self.mapped.then_some(value),
            NmPolicy::Auto => {
                reference_span.as_ref().map(|reference| compute_nm(reference, &parsed.ops, &bases))
            }
        };

        let mut flags = self.flags;
        let mut builder = RecordBuf::builder().set_name(self.name.as_str());

        if self.mapped {
            builder = builder
                .set_reference_sequence_id(resolve(&self.contig))
                .set_alignment_start(Position::new(self.start).expect("start must be >= 1"))
                .set_mapping_quality(MappingQuality::new(self.mapq).expect("MAPQ must be < 255"))
                .set_cigar(parsed.cigar);
        } else {
            flags |= Flags::UNMAPPED;
        }

        if let Some(mate_ref) = self.mate_ref_id {
            builder = builder.set_mate_reference_sequence_id(mate_ref);
        }
        if let Some(mate_start) = self.mate_start {
            let start = Position::new(mate_start).expect("mate start must be >= 1");
            builder = builder.set_mate_alignment_start(start);
        }
        if self.template_length != 0 {
            builder = builder.set_template_length(self.template_length);
        }

        builder = builder
            .set_flags(flags)
            .set_sequence(Sequence::from(bases))
            .set_quality_scores(QualityScores::from(quals));

        if !self.tags.is_empty() || nm.is_some() {
            let mut data = Data::default();
            for (tag, value) in self.tags {
                data.insert(tag.into(), value);
            }
            if let Some(nm) = nm {
                data.insert((*b"NM").into(), DataValue::UInt8(nm));
            }
            builder = builder.set_data(data);
        }

        builder.build()
    }

    /// The effective CIGAR string this read builds with — the explicit one, or the
    /// derived `{query_len}M`. Used to fill the mate's `MC` tag.
    fn effective_cigar_string(&self) -> String {
        if let Some(cigar) = &self.cigar {
            cigar.clone()
        } else {
            let len = if let Some(bases) = &self.bases { bases.len() } else { self.read_len };
            format!("{len}M")
        }
    }

    /// Reference bases this read spans, from its CIGAR — used for template length.
    fn ref_span(&self) -> usize {
        if let Some(cigar) = &self.cigar {
            cigar::parse(cigar).ref_len
        } else if let Some(bases) = &self.bases {
            bases.len()
        } else {
            self.read_len
        }
    }
}

/// Start a fluent builder for a read pair, mirroring fgbio's `SamBuilder.addPair`.
///
/// The mate reference/position, mate-strand and proper-pair flags, `MC` tag, and
/// signed template length are all *derived* from the two configured reads at build
/// time, so a test states only the reads themselves.
///
/// ```ignore
/// // FR pair, both 100 bp; mate fields, MC, and TLEN derived:
/// let (r1, r2) = pair("p1").at("chr1", 100, 300).build();
/// // Per-mate CIGARs — a spliced R1, a simple R2:
/// let (r1, r2) = pair("p1")
///     .r1(|r| r.at("chr1", 1500).cigar("50M2000N50M"))
///     .r2(|r| r.at("chr1", 3600).cigar("50M"))
///     .build();
/// ```
pub fn pair(name: impl Into<String>) -> PairBuilder {
    PairBuilder::new(name.into())
}

/// Fluent builder for a read pair. Construct via [`pair`].
pub struct PairBuilder {
    r1: ReadBuilder,
    r2: ReadBuilder,
    proper: bool,
    mate_cigar: bool,
    // `(contig, span)` from [`Self::fr`], applied at build time (when the final read
    // length is known) to position both mates and derive their CIGARs, soft-clipping
    // adapter read-through when the insert is shorter than the reads.
    fragment: Option<(String, Range<usize>)>,
}

impl PairBuilder {
    fn new(name: String) -> Self {
        Self {
            r1: read().name(name.clone()),
            r2: read().name(name),
            proper: true,
            mate_cigar: true,
            fragment: None,
        }
    }

    /// Place an FR pair on `contig`: read 1 forward at `start1`, read 2 reverse at
    /// `start2`. For other orientations or per-mate CIGARs, use [`Self::r1`] / [`Self::r2`].
    pub fn at(mut self, contig: impl Into<String>, start1: usize, start2: usize) -> Self {
        let contig = contig.into();
        self.r1 = self.r1.at(contig.clone(), start1);
        self.r2 = self.r2.at(contig, start2).reverse();
        self
    }

    /// Place an FR pair on `contig` spanning the 1-based, half-open reference range
    /// `span`: read 1 forward with its 5' end at `span.start`, read 2 reverse with its
    /// 5' end at `span.end - 1`, facing inward, with `|TLEN| == span.len()`.
    ///
    /// `.fr(c, s..e)` is the fragment-first counterpart to [`Self::at`]. When the
    /// fragment is at least a read long, both mates map fully (`{len}M`) and
    /// `.fr(c, s..e)` equals `.at(c, s, e - len)`. When the fragment is *shorter* than
    /// the reads, each mate sequences past its far end into adapter, so the read-through
    /// is soft-clipped: read 1 becomes `{n}M{overhang}S` and read 2 `{overhang}S{n}M`
    /// (`n = span.len()`), both mapping the whole fragment, with the clipped bases
    /// filled from the read-through adapters (see [`Adapters`]). This is the only way
    /// to express an insert shorter than the read length, since two fully-mapped reads
    /// always span at least one read length.
    ///
    /// The read length in effect at build time is used, so call [`Self::len`] before
    /// `fr` when overriding the default. Do not combine with [`Self::at`] or per-mate
    /// CIGAR customization.
    ///
    /// # Panics
    /// Panics if `span` is empty or inverted (`span.start >= span.end`).
    pub fn fr(mut self, contig: impl Into<String>, span: Range<usize>) -> Self {
        assert!(
            span.end > span.start,
            "fr() fragment span must be non-empty and forward (start < end), got {span:?}"
        );
        self.fragment = Some((contig.into(), span));
        self
    }

    /// Set both mates' read length (default 100).
    pub fn len(mut self, len: usize) -> Self {
        self.r1 = self.r1.len(len);
        self.r2 = self.r2.len(len);
        self
    }

    /// Set both mates' mapping quality (default 60).
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.r1 = self.r1.mapq(mapq);
        self.r2 = self.r2.mapq(mapq);
        self
    }

    /// Set both mates' CIGAR (a common shape for a symmetric pair); for different
    /// per-mate CIGARs use [`Self::r1`] / [`Self::r2`].
    pub fn cigar(mut self, cigar: impl Into<String>) -> Self {
        let cigar = cigar.into();
        self.r1 = self.r1.cigar(cigar.clone());
        self.r2 = self.r2.cigar(cigar);
        self
    }

    /// Set both mates' per-base qualities.
    pub fn quals(mut self, quals: impl Into<Vec<u8>>) -> Self {
        let quals = quals.into();
        self.r1 = self.r1.quals(quals.clone());
        self.r2 = self.r2.quals(quals);
        self
    }

    /// Derive both mates' bases from `reference` (see [`ReadBuilder::matching_ref`]).
    /// Per-mate substitutions can then be layered on via [`Self::r1`] / [`Self::r2`].
    pub fn matching_ref(mut self, reference: &FastaBuilder) -> Self {
        self.r1 = self.r1.matching_ref(reference);
        self.r2 = self.r2.matching_ref(reference);
        self
    }

    /// Customize read 1 with the full [`ReadBuilder`] API.
    pub fn r1(mut self, configure: impl FnOnce(ReadBuilder) -> ReadBuilder) -> Self {
        self.r1 = configure(self.r1);
        self
    }

    /// Customize read 2 with the full [`ReadBuilder`] API.
    pub fn r2(mut self, configure: impl FnOnce(ReadBuilder) -> ReadBuilder) -> Self {
        self.r2 = configure(self.r2);
        self
    }

    /// Mark the pair as not properly paired (clears the proper-pair flag).
    pub fn improper(mut self) -> Self {
        self.proper = false;
        self
    }

    /// Suppress the automatic `MC` (mate CIGAR) tags.
    pub fn no_mc(mut self) -> Self {
        self.mate_cigar = false;
        self
    }

    /// Suppress the derived `NM` tag on both mates.
    pub fn no_nm(mut self) -> Self {
        self.r1 = self.r1.no_nm();
        self.r2 = self.r2.no_nm();
        self
    }

    /// Mark both mates as failing vendor quality — a whole-pair property.
    pub fn qc_fail(mut self) -> Self {
        self.r1 = self.r1.qc_fail();
        self.r2 = self.r2.qc_fail();
        self
    }

    /// Build the pair as `(read1, read2)`, deriving every mate-linking field. Any
    /// [`Self::fr`] read-through fills from the default TruSeq [`Adapters`]; go through
    /// [`SamBuilder::add`] to use a SamBuilder's configured adapters instead.
    pub fn build(self) -> (RecordBuf, RecordBuf) {
        self.build_with_adapters(&Adapters::default(), &super::fasta_builder::default_contig_id)
    }

    /// Build the pair, filling any [`Self::fr`] read-through overhang with `adapters` and
    /// resolving each mate's contig name to a reference id via `resolve` (the frozen
    /// default table standalone, or the owning [`SamBuilder`]'s header when added to one).
    fn build_with_adapters(
        self,
        adapters: &Adapters,
        resolve: &dyn Fn(&str) -> usize,
    ) -> (RecordBuf, RecordBuf) {
        let PairBuilder { mut r1, mut r2, proper, mate_cigar, fragment } = self;

        // A fragment spec positions both mates and derives their CIGARs now that the
        // final read length is known. When the insert is shorter than the reads the
        // overhang is soft-clipped and filled with adapter read-through: read 1 reads
        // into the R1 adapter, read 2 into the R2 adapter (keyed to read number), each
        // oriented for its mate's strand (reverse-complemented on the reverse mate).
        if let Some((contig, span)) = fragment {
            let start = span.start;
            let insert = span.end.saturating_sub(span.start);
            let len = r1.read_len;
            let mapped = insert.min(len);
            let overhang = len - mapped;
            let (cigar1, cigar2) = if overhang == 0 {
                (format!("{len}M"), format!("{len}M"))
            } else {
                (format!("{mapped}M{overhang}S"), format!("{overhang}S{mapped}M"))
            };
            r1 = r1.at(contig.clone(), start).cigar(cigar1);
            r2 = r2.at(contig, start + (insert - mapped)).cigar(cigar2).reverse();
            if overhang > 0 {
                let r1_reverse = r1.flags.contains(Flags::REVERSE_COMPLEMENTED);
                let r2_reverse = r2.flags.contains(Flags::REVERSE_COMPLEMENTED);
                r1.readthrough = Some(readthrough_fill(&adapters.r1, overhang, r1_reverse));
                r2.readthrough = Some(readthrough_fill(&adapters.r2, overhang, r2_reverse));
            }
        }

        r1.flags |= Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        r2.flags |= Flags::SEGMENTED | Flags::LAST_SEGMENT;

        // Mate-strand / mate-unmapped flags mirror the *other* read's state.
        if r2.flags.contains(Flags::REVERSE_COMPLEMENTED) {
            r1.flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        if r1.flags.contains(Flags::REVERSE_COMPLEMENTED) {
            r2.flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        if !r2.mapped {
            r1.flags |= Flags::MATE_UNMAPPED;
        }
        if !r1.mapped {
            r2.flags |= Flags::MATE_UNMAPPED;
        }

        // Proper only when both mates map to the same contig, unless overridden. This is a
        // name comparison, independent of the reference-id mapping.
        let same_contig = r1.mapped && r2.mapped && r1.contig == r2.contig;
        if proper && same_contig {
            r1.flags |= Flags::PROPERLY_SEGMENTED;
            r2.flags |= Flags::PROPERLY_SEGMENTED;
        }

        if r2.mapped {
            r1.mate_ref_id = Some(resolve(&r2.contig));
            r1.mate_start = Some(r2.start);
        }
        if r1.mapped {
            r2.mate_ref_id = Some(resolve(&r1.contig));
            r2.mate_start = Some(r1.start);
        }

        let (tlen1, tlen2) = template_lengths(&r1, &r2);
        r1.template_length = tlen1;
        r2.template_length = tlen2;

        // Each mate's MC is the other's CIGAR, unless suppressed or already set
        // (e.g. a deliberately malformed MC a test wants to keep).
        if mate_cigar {
            if r2.mapped && !has_tag(&r1.tags, *b"MC") {
                let mate = r2.effective_cigar_string();
                r1.tags.push((*b"MC", DataValue::String(mate.into_bytes().into())));
            }
            if r1.mapped && !has_tag(&r2.tags, *b"MC") {
                let mate = r1.effective_cigar_string();
                r2.tags.push((*b"MC", DataValue::String(mate.into_bytes().into())));
            }
        }

        (r1.build_with(resolve), r2.build_with(resolve))
    }
}

/// Resolve a contig `name` to its `reference_sequence_id` in `header` — its position among
/// the header's `@SQ` lines. This is what [`SamBuilder::add`] uses so records resolve
/// against the builder's own header rather than the frozen default table.
///
/// # Panics
/// Panics if `header` declares no contig named `name`, so a misnamed contig fails loudly
/// when added instead of silently writing a record onto the wrong (or an out-of-range)
/// reference.
fn contig_id_in_header(header: &Header, name: &str) -> usize {
    header
        .reference_sequences()
        .keys()
        .position(|key| {
            let key: &[u8] = key;
            key == name.as_bytes()
        })
        .unwrap_or_else(|| panic!("contig {name:?} is not declared in the SamBuilder header"))
}

/// Whether `tags` already carries the two-byte `tag`.
fn has_tag(tags: &[([u8; 2], DataValue)], tag: [u8; 2]) -> bool {
    tags.iter().any(|(key, _)| *key == tag)
}

/// Signed template lengths `(tlen1, tlen2)` per the SAM spec: the leftmost mapped
/// segment is positive, the rightmost negative, with magnitude equal to the number
/// of reference bases the template spans. Zero across contigs or when either mate
/// is unmapped.
fn template_lengths(r1: &ReadBuilder, r2: &ReadBuilder) -> (i32, i32) {
    if !r1.mapped || !r2.mapped || r1.contig != r2.contig {
        return (0, 0);
    }
    let end1 = r1.start + r1.ref_span() - 1;
    let end2 = r2.start + r2.ref_span() - 1;
    let leftmost = r1.start.min(r2.start);
    let rightmost = end1.max(end2);
    let span = i32::try_from(rightmost - leftmost + 1).expect("template length fits in i32");
    if r1.start <= r2.start { (span, -span) } else { (-span, span) }
}

/// Derive a read's bases by copying `reference` along `ops`: reference-consuming
/// ops copy the upper-cased reference base (reads are upper-case even over a
/// soft-masked reference), an explicit mismatch (`X`) op emits a differing base, and
/// inserted bases fill with `A`. Soft-clip bases are laid in from `readthrough`
/// (adapter read-through, in order) when supplied, else `A`. `reference` starts at
/// the read's alignment position.
fn derive_bases(
    reference: &[u8],
    ops: &[(Kind, usize)],
    query_len: usize,
    readthrough: Option<&[u8]>,
) -> Vec<u8> {
    let mut bases = Vec::with_capacity(query_len);
    let mut reference_offset = 0usize;
    let mut readthrough_offset = 0usize;
    for &(kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch => {
                for _ in 0..len {
                    let base = reference.get(reference_offset).copied().unwrap_or(b'A');
                    bases.push(base.to_ascii_uppercase());
                    reference_offset += 1;
                }
            }
            Kind::SequenceMismatch => {
                for _ in 0..len {
                    let base = reference.get(reference_offset).copied().unwrap_or(b'A');
                    bases.push(mismatch_base(base));
                    reference_offset += 1;
                }
            }
            Kind::Insertion => bases.resize(bases.len() + len, b'A'),
            Kind::SoftClip => {
                for _ in 0..len {
                    let base = readthrough
                        .and_then(|fill| fill.get(readthrough_offset).copied())
                        .unwrap_or(b'A');
                    bases.push(base);
                    readthrough_offset += 1;
                }
            }
            Kind::Deletion | Kind::Skip => reference_offset += len,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    bases
}

/// A base that differs from `reference` (upper-cased), for an `X` mismatch op.
fn mismatch_base(reference: u8) -> u8 {
    if reference.eq_ignore_ascii_case(&b'A') { b'C' } else { b'A' }
}

/// The soft-clip bases for a read that sequenced past a short insert: the first
/// `overhang` bases of `adapter`, then poly-A (the sequencer keeps calling bases
/// after the adapter runs out). For a `reverse`-strand read the result is
/// reverse-complemented, because its stored `SEQ` is the reverse complement of what
/// came off the machine — so the trailing poly-A lands as leading poly-T and the
/// adapter as its reverse complement.
fn readthrough_fill(adapter: &[u8], overhang: usize, reverse: bool) -> Vec<u8> {
    let tail: Vec<u8> = (0..overhang).map(|i| adapter.get(i).copied().unwrap_or(b'A')).collect();
    if reverse { reverse_complement(&tail) } else { tail }
}

/// The reverse complement of `seq` (A/C/G/T swapped; other bytes pass through).
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement(base)).collect()
}

/// The Watson–Crick complement of an upper-cased base; non-ACGT bytes pass through.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// The `NM` edit distance for `bases` aligned to `reference` by `ops`: base
/// mismatches (case-insensitive, so a soft-masked reference still matches) plus
/// inserted and deleted bases. `N` (reference skip) does not count.
fn compute_nm(reference: &[u8], ops: &[(Kind, usize)], bases: &[u8]) -> u8 {
    let mut nm = 0usize;
    let mut reference_offset = 0usize;
    let mut query_offset = 0usize;
    for &(kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let reference_base = reference.get(reference_offset).copied().unwrap_or(b'A');
                    let read_base = bases.get(query_offset).copied().unwrap_or(b'A');
                    if !read_base.eq_ignore_ascii_case(&reference_base) {
                        nm += 1;
                    }
                    reference_offset += 1;
                    query_offset += 1;
                }
            }
            Kind::Insertion => {
                nm += len;
                query_offset += len;
            }
            Kind::Deletion => {
                nm += len;
                reference_offset += len;
            }
            Kind::Skip => reference_offset += len,
            Kind::SoftClip => query_offset += len,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    u8::try_from(nm).unwrap_or(u8::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Uppercase a byte vector — reads are upper-case even over soft-masked reference.
    fn upper(bases: &[u8]) -> Vec<u8> {
        bases.iter().map(u8::to_ascii_uppercase).collect()
    }

    /// A `Data` block carrying just an `NM` tag, for the derived-NM assertions.
    fn nm_data(nm: u8) -> Data {
        let mut data = Data::default();
        data.insert((*b"NM").into(), DataValue::UInt8(nm));
        data
    }

    #[test]
    fn defaults_to_a_mapped_forward_primary_read() {
        let rec = read().at("chr1", 100).build();
        assert_eq!(rec.flags(), Flags::empty());
        assert_eq!(rec.reference_sequence_id(), Some(0));
        assert_eq!(rec.alignment_start(), Position::new(100));
        assert_eq!(rec.mapping_quality(), MappingQuality::new(60));
        assert_eq!(*rec.cigar(), cigar::parse("100M").cigar);
        assert_eq!(rec.sequence().len(), 100); // content now derives from the reference
        assert_eq!(*rec.quality_scores(), QualityScores::from(vec![30u8; 100]));
    }

    #[test]
    fn cigar_sizes_sequence_and_quals_to_query_length() {
        let rec = read().at("chr1", 1).cigar("10S80M10S").build();
        assert_eq!(*rec.cigar(), cigar::parse("10S80M10S").cigar);
        assert_eq!(rec.sequence().len(), 100);
        assert_eq!(rec.quality_scores().len(), 100);
    }

    #[test]
    fn explicit_bases_and_quals_are_used_verbatim() {
        let rec =
            read().at("chr1", 1).cigar("4M").bases(*b"ACGT").quals(vec![10u8, 20, 30, 40]).build();
        assert_eq!(*rec.sequence(), Sequence::from(b"ACGT".to_vec()));
        assert_eq!(*rec.quality_scores(), QualityScores::from(vec![10u8, 20, 30, 40]));
    }

    #[test]
    fn a_mapped_read_derives_matching_bases_from_the_reference() {
        let rec = read().at("chr1", 100).cigar("50M").build();
        let want = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 50));
        assert_eq!(*rec.sequence(), Sequence::from(want));
    }

    #[test]
    fn a_read_over_a_soft_masked_region_is_uppercased() {
        // chr1 [2000, 2500) (0-based) is soft-masked, so 1-based 2100 lands inside it.
        let masked = FastaBuilder::hg38_mini().bases("chr1", 2100, 50);
        assert!(masked.iter().any(u8::is_ascii_lowercase), "fixture should be masked here");
        let rec = read().at("chr1", 2100).cigar("50M").build();
        assert_eq!(*rec.sequence(), Sequence::from(upper(&masked)));
    }

    #[test]
    fn sub_plants_a_mismatch_against_the_reference() {
        let rec = read().at("chr1", 100).cigar("10M").sub(3, b'N').build();
        let mut want = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 10));
        want[3] = b'N';
        assert_eq!(*rec.sequence(), Sequence::from(want));
    }

    #[test]
    fn matching_ref_overrides_the_default_reference() {
        let custom = FastaBuilder::new().contig("chr1", vec![b'G'; 100]);
        let rec = read().at("chr1", 1).cigar("10M").matching_ref(&custom).build();
        assert_eq!(*rec.sequence(), Sequence::from(vec![b'G'; 10]));
    }

    #[test]
    fn a_matching_read_gets_nm_zero() {
        assert_eq!(*read().at("chr1", 100).cigar("50M").build().data(), nm_data(0));
    }

    #[test]
    fn a_substitution_bumps_nm() {
        let rec = read().at("chr1", 100).cigar("10M").sub(3, b'N').build();
        assert_eq!(*rec.data(), nm_data(1));
    }

    #[test]
    fn insertions_and_deletions_count_toward_nm() {
        // 5 inserted + 5 deleted bases, matching M elsewhere.
        let rec = read().at("chr1", 100).cigar("20M5I20M5D20M").build();
        assert_eq!(*rec.data(), nm_data(10));
    }

    #[test]
    fn a_soft_masked_reference_does_not_inflate_nm() {
        // The read is upper-case over the lowercase [2000, 2500) region; NM must be 0.
        let rec = read().at("chr1", 2100).cigar("50M").build();
        assert_eq!(*rec.data(), nm_data(0));
    }

    #[test]
    fn no_nm_suppresses_the_nm_tag() {
        let rec = read().at("chr1", 100).cigar("50M").no_nm().build();
        assert_eq!(*rec.data(), Data::default());
    }

    #[test]
    fn explicit_nm_overrides_the_derived_value() {
        let rec = read().at("chr1", 100).cigar("50M").nm(7).build();
        assert_eq!(*rec.data(), nm_data(7));
    }

    #[test]
    fn name_defaults_and_is_settable() {
        assert_eq!(read().at("chr1", 1).build().name().unwrap().to_vec(), b"read".to_vec());
        assert_eq!(
            read().name("q1").at("chr1", 1).build().name().unwrap().to_vec(),
            b"q1".to_vec()
        );
    }

    #[test]
    fn mapq_is_settable() {
        assert_eq!(
            read().at("chr1", 1).mapq(20).build().mapping_quality(),
            MappingQuality::new(20)
        );
    }

    #[test]
    fn reverse_sets_the_reverse_flag() {
        assert!(read().at("chr1", 1).reverse().build().flags().is_reverse_complemented());
    }

    #[test]
    fn secondary_and_supplementary_set_their_flags() {
        assert!(read().at("chr1", 1).secondary().build().flags().is_secondary());
        assert!(read().at("chr1", 1).supplementary().build().flags().is_supplementary());
    }

    #[test]
    fn duplicate_and_qc_fail_set_their_flags() {
        assert!(read().at("chr1", 1).duplicate().build().flags().is_duplicate());
        assert!(read().at("chr1", 1).qc_fail().build().flags().is_qc_fail());
    }

    #[test]
    fn unmapped_clears_reference_position_and_cigar() {
        let rec = read().unmapped().build();
        assert!(rec.flags().is_unmapped());
        assert_eq!(rec.reference_sequence_id(), None);
        assert_eq!(rec.alignment_start(), None);
        assert_eq!(*rec.cigar(), cigar::parse("").cigar);
        // An unmapped read still carries sequence and qualities; with no reference
        // to derive from, the bases default to all-A.
        assert_eq!(*rec.sequence(), Sequence::from(vec![b'A'; 100]));
    }

    #[test]
    fn nm_tag_is_attached() {
        let rec = read().at("chr1", 1).cigar("100M").nm(3).build();
        let mut want = Data::default();
        want.insert((*b"NM").into(), DataValue::UInt8(3));
        assert_eq!(*rec.data(), want);
    }

    #[test]
    fn mc_tag_is_attached() {
        let rec = read().at("chr1", 1).mc("50M").no_nm().build();
        let mut want = Data::default();
        want.insert((*b"MC").into(), DataValue::String(b"50M".to_vec().into()));
        assert_eq!(*rec.data(), want);
    }

    // ─── pairs ─────────────────────────────────────────────────────────────

    #[test]
    fn fr_pair_derives_segment_and_mate_strand_flags() {
        let (r1, r2) = pair("p1").at("chr1", 100, 300).build();
        assert!(r1.flags().is_segmented() && r1.flags().is_first_segment());
        assert!(r1.flags().is_properly_segmented());
        assert!(!r1.flags().is_reverse_complemented());
        assert!(r1.flags().is_mate_reverse_complemented());
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_reverse_complemented());
        assert!(!r2.flags().is_mate_reverse_complemented());
    }

    #[test]
    fn fr_pair_derives_mate_positions_and_template_length() {
        let (r1, r2) = pair("p1").at("chr1", 100, 300).build();
        assert_eq!(r1.mate_reference_sequence_id(), Some(0));
        assert_eq!(r1.mate_alignment_start(), Position::new(300));
        assert_eq!(r2.mate_alignment_start(), Position::new(100));
        // 100 bp mates at 100 and 300 span bases 100..=399 = 300; R1 is leftmost.
        assert_eq!(r1.template_length(), 300);
        assert_eq!(r2.template_length(), -300);
    }

    #[test]
    fn matching_ref_derives_both_mates_from_the_reference() {
        let custom = FastaBuilder::new().contig("chr1", vec![b'G'; 400]);
        let (r1, r2) =
            pair("p1").at("chr1", 100, 300).matching_ref(&custom).r1(|r| r.sub(2, b'A')).build();
        let mut want_r1 = vec![b'G'; 100];
        want_r1[2] = b'A';
        assert_eq!(*r1.sequence(), Sequence::from(want_r1));
        assert_eq!(*r2.sequence(), Sequence::from(vec![b'G'; 100]));
    }

    #[test]
    fn mate_cigar_tags_are_derived_from_the_other_read() {
        let (r1, r2) = pair("p1")
            .no_nm()
            .r1(|r| r.at("chr1", 1500).cigar("50M2000N50M"))
            .r2(|r| r.at("chr1", 3600).cigar("50M"))
            .build();
        let mut want_r1 = Data::default();
        want_r1.insert((*b"MC").into(), DataValue::String(b"50M".to_vec().into()));
        assert_eq!(*r1.data(), want_r1);
        let mut want_r2 = Data::default();
        want_r2.insert((*b"MC").into(), DataValue::String(b"50M2000N50M".to_vec().into()));
        assert_eq!(*r2.data(), want_r2);
    }

    #[test]
    fn no_mc_suppresses_the_mate_cigar_tag() {
        let (r1, _) = pair("p1").at("chr1", 100, 300).no_mc().no_nm().build();
        assert_eq!(*r1.data(), Data::default());
    }

    #[test]
    fn an_explicit_mate_cigar_is_not_overwritten() {
        let (r1, _) = pair("p1")
            .no_nm()
            .r1(|r| r.at("chr1", 100).mc("garbage"))
            .r2(|r| r.at("chr1", 300).reverse())
            .build();
        let mut want = Data::default();
        want.insert((*b"MC").into(), DataValue::String(b"garbage".to_vec().into()));
        assert_eq!(*r1.data(), want);
    }

    #[test]
    fn improper_clears_the_proper_pair_flag() {
        let (r1, _) = pair("p1").at("chr1", 100, 300).improper().build();
        assert!(r1.flags().is_segmented());
        assert!(!r1.flags().is_properly_segmented());
    }

    #[test]
    fn qc_fail_marks_both_mates() {
        let (r1, r2) = pair("p1").at("chr1", 100, 300).qc_fail().build();
        assert!(r1.flags().is_qc_fail());
        assert!(r2.flags().is_qc_fail());
    }

    #[test]
    fn cigar_and_quals_apply_to_both_mates() {
        let (r1, r2) =
            pair("p1").at("chr1", 100, 300).cigar("10S90M").quals(vec![25u8; 100]).build();
        assert_eq!(*r1.cigar(), cigar::parse("10S90M").cigar);
        assert_eq!(*r2.cigar(), cigar::parse("10S90M").cigar);
        assert_eq!(*r1.quality_scores(), QualityScores::from(vec![25u8; 100]));
        assert_eq!(*r2.quality_scores(), QualityScores::from(vec![25u8; 100]));
    }

    #[test]
    fn cross_contig_pair_has_zero_tlen_and_is_not_proper() {
        let (r1, r2) =
            pair("p1").r1(|r| r.at("chr1", 100)).r2(|r| r.at("chr2", 100).reverse()).build();
        assert_eq!(r1.mate_reference_sequence_id(), Some(1));
        assert_eq!(r2.mate_reference_sequence_id(), Some(0));
        assert_eq!(r1.template_length(), 0);
        assert_eq!(r2.template_length(), 0);
        assert!(!r1.flags().is_properly_segmented());
    }

    #[test]
    fn an_unmapped_mate_sets_the_mate_unmapped_flag() {
        let (r1, r2) = pair("p1").r1(|r| r.at("chr1", 100)).r2(ReadBuilder::unmapped).build();
        assert!(r1.flags().is_mate_unmapped());
        assert!(r2.flags().is_unmapped());
        assert!(!r1.flags().is_properly_segmented());
    }

    #[test]
    fn fr_at_least_read_length_maps_both_mates_fully() {
        // A 300 bp fragment with 100 bp reads is just the plain FR pair
        // `.at("chr1", 100, 300)`: both mates map fully, R2 offset to span 300.
        let (r1, r2) = pair("p1").fr("chr1", 100..400).build();
        assert_eq!(*r1.cigar(), cigar::parse("100M").cigar);
        assert_eq!(*r2.cigar(), cigar::parse("100M").cigar);
        assert_eq!(r1.alignment_start(), Position::new(100));
        assert_eq!(r2.alignment_start(), Position::new(300));
        assert_eq!(r1.template_length(), 300);
        assert_eq!(r2.template_length(), -300);
    }

    #[test]
    fn fr_equal_to_read_length_fully_overlaps_without_clipping() {
        let (r1, r2) = pair("p1").fr("chr1", 100..200).build();
        assert_eq!(*r1.cigar(), cigar::parse("100M").cigar);
        assert_eq!(*r2.cigar(), cigar::parse("100M").cigar);
        assert_eq!(r1.alignment_start(), Position::new(100));
        assert_eq!(r2.alignment_start(), Position::new(100));
        assert_eq!(r1.template_length(), 100);
        assert_eq!(r2.template_length(), -100);
    }

    #[test]
    fn fr_shorter_than_reads_soft_clips_the_readthrough() {
        // A 90 bp fragment with 100 bp reads: each mate sequences 10 bp past the far
        // end into adapter, which is soft-clipped. Both map the whole [100, 189]
        // fragment and |TLEN| is the 90 bp insert, not the 100 bp read.
        let (r1, r2) = pair("p1").fr("chr1", 100..190).build();
        assert_eq!(*r1.cigar(), cigar::parse("90M10S").cigar);
        assert_eq!(*r2.cigar(), cigar::parse("10S90M").cigar);
        assert_eq!(r1.alignment_start(), Position::new(100));
        assert_eq!(r2.alignment_start(), Position::new(100));
        assert_eq!(r1.template_length(), 90);
        assert_eq!(r2.template_length(), -90);
    }

    #[test]
    fn fr_uses_the_read_length_in_effect_at_build_time() {
        // `.len()` after `.fr()` still takes effect (CIGARs derive at build time): a
        // 60 bp fragment with 75 bp reads clips 15 bp of read-through.
        let (r1, r2) = pair("p1").fr("chr1", 100..160).len(75).build();
        assert_eq!(*r1.cigar(), cigar::parse("60M15S").cigar);
        assert_eq!(*r2.cigar(), cigar::parse("15S60M").cigar);
        assert_eq!(r1.template_length(), 60);
        assert_eq!(r2.template_length(), -60);
    }

    #[test]
    fn fr_readthrough_fills_the_overhang_with_adapter() {
        // 90 bp fragment, 100 bp reads → a 10 bp read-through on each mate.
        let (r1, r2) = pair("p1").fr("chr1", 100..190).build();
        let genome = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 90));

        // R1 (forward): the genomic bases, then the first 10 bp of the R1 adapter, as-is.
        let mut want_r1 = genome.clone();
        want_r1.extend_from_slice(&TRUSEQ_R1[..10]);
        assert_eq!(*r1.sequence(), Sequence::from(want_r1));

        // R2 (reverse): the R2 adapter read-through reverse-complemented, then genomic.
        let mut want_r2 = reverse_complement(&TRUSEQ_R2[..10]);
        want_r2.extend_from_slice(&genome);
        assert_eq!(*r2.sequence(), Sequence::from(want_r2));
    }

    #[test]
    fn fr_readthrough_past_the_adapter_is_poly_a() {
        // 5 bp fragment, 50 bp reads → a 45 bp overhang: the ~33 bp adapter, then poly-A.
        let (r1, r2) = pair("p1").fr("chr1", 100..105).len(50).build();
        let genome = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 5));

        // R1 (forward): 5 genomic bases, the full R1 adapter, then poly-A to 50 bp.
        let mut want_r1 = genome.clone();
        want_r1.extend_from_slice(TRUSEQ_R1);
        want_r1.resize(50, b'A');
        assert_eq!(*r1.sequence(), Sequence::from(want_r1));

        // R2 (reverse): the as-sequenced tail (R2 adapter + poly-A) reverse-complemented
        // — so leading poly-T, then revcomp(adapter) — then the 5 genomic bases.
        let mut r2_tail = TRUSEQ_R2.to_vec();
        r2_tail.resize(45, b'A');
        let mut want_r2 = reverse_complement(&r2_tail);
        want_r2.extend_from_slice(&genome);
        assert_eq!(*r2.sequence(), Sequence::from(want_r2));
    }

    #[test]
    fn readthrough_adapter_is_keyed_to_read_number_orientation_to_strand() {
        // Distinct 10 bp adapters make the wiring unambiguous: R1 gets adapter r1, R2
        // gets adapter r2 (by read number), each oriented for its strand.
        let adapters = Adapters { r1: b"TTTTTTTTTT".to_vec(), r2: b"GGGGGGGGGG".to_vec() };
        let records = pair("p1")
            .fr("chr1", 100..190)
            .into_records(&adapters, &crate::test_support::default_contig_id);
        let genome = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 90));

        // R1 (forward, first read) → adapter r1, as-is.
        let mut want_r1 = genome.clone();
        want_r1.extend_from_slice(b"TTTTTTTTTT");
        assert_eq!(*records[0].sequence(), Sequence::from(want_r1));

        // R2 (reverse, second read) → adapter r2, reverse-complemented (GGGG.. → CCCC..).
        let mut want_r2 = reverse_complement(b"GGGGGGGGGG");
        want_r2.extend_from_slice(&genome);
        assert_eq!(*records[1].sequence(), Sequence::from(want_r2));
    }

    #[test]
    fn sam_builder_adapter_overrides_drive_read_through() {
        // The read-through adapters are configurable per-SamBuilder: a `pair().fr()` added
        // through the builder fills its overhang from the builder's adapters, not TruSeq.
        let mut sam =
            SamBuilder::new().adapter_r1(b"TTTTTTTTTT".to_vec()).adapter_r2(b"GGGGGGGGGG".to_vec());
        sam.add(pair("p1").fr("chr1", 100..190));
        let genome = upper(&FastaBuilder::hg38_mini().bases("chr1", 100, 90));

        // R1 (forward, first read) → the builder's r1 adapter, as-is.
        let mut want_r1 = genome.clone();
        want_r1.extend_from_slice(b"TTTTTTTTTT");
        assert_eq!(*sam.records[0].sequence(), Sequence::from(want_r1));

        // R2 (reverse, second read) → the builder's r2 adapter, reverse-complemented.
        let mut want_r2 = reverse_complement(b"GGGGGGGGGG");
        want_r2.extend_from_slice(&genome);
        assert_eq!(*sam.records[1].sequence(), Sequence::from(want_r2));
    }

    #[test]
    fn next_id_increments_from_one() {
        let sam = SamBuilder::new();
        assert_eq!(sam.next_id(), "1");
        assert_eq!(sam.next_id(), "2");
        assert_eq!(sam.next_id(), "3");
    }

    #[test]
    fn next_id_uses_the_prefix_when_set() {
        let sam = SamBuilder::new().id_prefix("q");
        assert_eq!(sam.next_id(), "q1");
        assert_eq!(sam.next_id(), "q2");
    }

    #[test]
    fn next_id_is_callable_inline_inside_add() {
        // The `&self` + `Cell` design lets `next_id` be called in `add`'s argument
        // without fighting `add`'s `&mut self` borrow. Two ids handed out → next is 3.
        let mut sam = SamBuilder::new();
        sam.add(pair(sam.next_id()).at("chr1", 100, 300));
        sam.add(read().name(sam.next_id()).at("chr1", 500));
        assert_eq!(sam.next_id(), "3");
    }

    #[test]
    fn add_takes_reads_pairs_and_records_and_writes_a_bam() {
        let mut sam = SamBuilder::new();
        sam.add(read().at("chr1", 100)) // 1 record
            .add(pair("p").at("chr1", 100, 300)) // 2 records
            .add(read().at("chr1", 200).build()); // 1 raw RecordBuf
        let bam = sam.to_temp_bam().unwrap();

        let mut reader = bam::io::Reader::new(std::fs::File::open(bam.path()).unwrap());
        let _ = reader.read_header().unwrap();
        let mut record = bam::Record::default();
        let mut count = 0;
        while reader.read_record(&mut record).unwrap() != 0 {
            count += 1;
        }
        assert_eq!(count, 4, "1 read + 1 pair + 1 record = 4 records");
    }

    #[test]
    fn add_resolves_contig_id_against_the_builder_header_not_the_frozen_table() {
        // Header declares chr2 first (id 0), chr1 second (id 1) — the reverse of
        // DEFAULT_CONTIGS. A read added on chr2 must land on the header's chr2 slot (0),
        // not the frozen table's chr2 index (1).
        let mut sam =
            SamBuilder::with_contigs(&[("chr2".to_string(), 1000), ("chr1".to_string(), 1000)]);
        sam.add(pair("p").at("chr2", 100, 300)).add(read().name("r").at("chr1", 100));
        assert_eq!(sam.records[0].reference_sequence_id(), Some(0), "chr2 -> header id 0");
        assert_eq!(sam.records[0].mate_reference_sequence_id(), Some(0), "mate chr2 -> id 0");
        assert_eq!(sam.records[2].reference_sequence_id(), Some(1), "chr1 -> header id 1");
    }

    #[test]
    fn standalone_build_still_resolves_against_the_frozen_default_table() {
        // Without a SamBuilder there is no header, so chr2 resolves to the frozen index 1.
        assert_eq!(read().at("chr2", 100).build().reference_sequence_id(), Some(1));
    }

    #[test]
    #[should_panic(expected = "not declared in the SamBuilder header")]
    fn add_panics_when_a_read_names_a_contig_absent_from_the_header() {
        // Header has only chr1; adding a chr2 read fails loudly at add() rather than
        // stamping an out-of-range reference id into the record.
        let mut sam = SamBuilder::with_contigs(&[("chr1".to_string(), 1000)]);
        sam.add(read().at("chr2", 100));
    }

    #[test]
    #[should_panic(expected = "fragment span must be non-empty")]
    fn fr_rejects_an_empty_or_inverted_span() {
        let _ = pair("p").fr("chr1", 100..100);
    }

    #[test]
    #[should_panic(expected = "start must be >= 1")]
    fn a_mapped_read_with_start_zero_panics_clearly() {
        let _ = read().at("chr1", 0).build();
    }

    #[test]
    fn read_group_adds_an_rg_with_the_sample_tag() {
        use noodles::sam::header::record::value::map::read_group::tag;
        let sam = SamBuilder::new().read_group("rg0", "NA12878");
        let (id, rg) = sam.header().read_groups().iter().next().expect("one read group");
        assert_eq!(id.to_string(), "rg0");
        assert_eq!(
            rg.other_fields().get(&tag::SAMPLE).map(ToString::to_string),
            Some("NA12878".to_string())
        );
    }

    #[test]
    fn read_group_without_sample_omits_the_sample_tag() {
        use noodles::sam::header::record::value::map::read_group::tag;
        let sam = SamBuilder::new().read_group_without_sample("rg0");
        let (_, rg) = sam.header().read_groups().iter().next().expect("one read group");
        assert!(rg.other_fields().get(&tag::SAMPLE).is_none());
    }
}
