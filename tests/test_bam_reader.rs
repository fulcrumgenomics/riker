use std::io::Write;
use std::path::Path;

use anyhow::Result;
use noodles::core::Region;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::Op;
use riker_lib::sam::alignment_reader::{AlignmentReader, IndexedAlignmentReader};
use riker_lib::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use riker_lib::test_support::{
    FastaBuilder, SamBuilder, SortOrder, coord_builder, count_stored_cram_tags, pair, read,
};
use tempfile::NamedTempFile;

/// The comparable fields of one record: every field `RikerRecord` always
/// populates, plus the sequence and the `NM` and `MC` aux tags, which are
/// populated only when requested. Used to assert that two readers surface
/// identical records.
#[derive(Debug, PartialEq)]
struct Row {
    name: Option<Vec<u8>>,
    flags: u16,
    reference_sequence_id: Option<usize>,
    start: Option<usize>,
    end: Option<usize>,
    mapping_quality: Option<u8>,
    cigar: Vec<Op>,
    mate_reference_sequence_id: Option<usize>,
    mate_start: Option<usize>,
    template_length: i32,
    sequence_len: usize,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    nm: Option<i64>,
    mc: Option<Vec<u8>>,
}

/// Project a `RikerRecord` onto its comparable [`Row`].
fn record_row(record: &RikerRecord) -> Row {
    Row {
        name: record.name().map(|n| n.to_vec()),
        flags: record.flags().bits(),
        reference_sequence_id: record.reference_sequence_id(),
        start: record.alignment_start().map(|p| p.get()),
        end: record.alignment_end().map(|p| p.get()),
        mapping_quality: record.mapping_quality().map(|q| q.get()),
        cigar: record.cigar_ops().collect(),
        mate_reference_sequence_id: record.mate_reference_sequence_id(),
        mate_start: record.mate_alignment_start().map(|p| p.get()),
        template_length: record.template_length(),
        sequence_len: record.sequence_len(),
        sequence: record.sequence().to_vec(),
        quality_scores: record.quality_scores().to_vec(),
        nm: record.aux_tag(*b"NM").and_then(|v| v.as_int()),
        mc: record.aux_tag(*b"MC").and_then(|v| v.as_str().map(|s| s.to_vec())),
    }
}

/// Read every record from `path` sequentially with the given requirements and
/// decode-thread count, projecting each to a [`Row`]. `reference` is needed
/// only for CRAM.
fn collect_stream_rows(
    path: &Path,
    reference: Option<&Path>,
    requirements: &RikerRecordRequirements,
    decode_threads: usize,
) -> Result<Vec<Row>> {
    let mut reader = AlignmentReader::open(path, reference, decode_threads)?;
    reader.set_requirements(requirements.clone())?;
    let mut record = reader.empty_record();
    let mut rows = Vec::new();
    while reader.fill_record(&mut record)? {
        rows.push(record_row(&record));
    }
    Ok(rows)
}

/// A reference and records that between them exercise every field a CRAM
/// decode reconstructs: an FR pair (mate fields, template length, `MC`), a
/// pair split across chr1 and chr2 (mate fields, which CRAM stores per read
/// rather than linking the mates when they sit in different slices), a
/// reverse-strand read with varied qualities and planted mismatches (`NM`,
/// quality orientation), clipping and indels (CIGAR), and an unmapped read.
fn varied_records() -> (FastaBuilder, SamBuilder) {
    let reference = FastaBuilder::hg38_mini();
    let mut builder = SamBuilder::with_contigs(&reference.dict()).sort_order(SortOrder::Coordinate);
    builder.add(pair("pair1").at("chr1", 100, 300).len(50));
    builder.add(pair("split_pair").at("chr1", 5_000, 5_000).r2(|r2| r2.at("chr2", 5_000)));
    let reverse_quals: Vec<u8> = (0..60u8).map(|i| 10 + i % 30).collect();
    builder.add(
        read()
            .name("reverse_with_mismatches")
            .at("chr1", 1_000)
            .len(60)
            .mapq(30)
            .reverse()
            .quals(reverse_quals)
            .sub(10, b'A')
            .sub(20, b'C'),
    );
    builder.add(read().name("clipped_with_indels").at("chr2", 2_000).cigar("5S20M2I20M3D30M5S"));
    builder.add(read().name("unmapped").unmapped());
    (reference, builder)
}

/// [`varied_records`] written as `(fasta, bam, cram)`, the CRAM through htslib,
/// which drops each `NM` tag for the reader to regenerate.
fn varied_bam_and_cram() -> Result<(NamedTempFile, NamedTempFile, NamedTempFile)> {
    let (reference, builder) = varied_records();
    let fasta = reference.to_temp_fasta()?;
    let bam = builder.to_temp_bam()?;
    let cram = builder.to_temp_htslib_cram(fasta.path())?;
    Ok((fasta, bam, cram))
}

/// Read [`varied_bam_and_cram`]'s BAM and CRAM with `requirements` and return
/// the `(bam, cram)` rows.
fn bam_and_cram_rows(requirements: &RikerRecordRequirements) -> Result<(Vec<Row>, Vec<Row>)> {
    let (fasta, bam, cram) = varied_bam_and_cram()?;
    let bam_rows = collect_stream_rows(bam.path(), None, requirements, 0)?;
    let mut cram_rows = collect_stream_rows(cram.path(), Some(fasta.path()), requirements, 0)?;
    assert_eq!(bam_rows.len(), 7, "two pairs + two fragments + one unmapped read");
    clear_unmapped_mapq(&mut cram_rows);
    Ok((bam_rows, cram_rows))
}

/// htslib decodes an unmapped CRAM record's MAPQ as 0 whatever it is asked to decode,
/// where the BAM has 255 (missing). That is a CRAM-format difference, not one the
/// BAM-vs-CRAM comparisons are after, so clear it on the CRAM side.
fn clear_unmapped_mapq(cram_rows: &mut [Row]) {
    for row in cram_rows {
        if row.flags & Flags::UNMAPPED.bits() != 0 {
            row.mapping_quality = None;
        }
    }
}

#[test]
fn test_open_nonexistent_file() {
    match AlignmentReader::open(Path::new("/no/such/file.bam"), None, 0) {
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("/no/such/file.bam"), "error should contain path, got: {msg}");
        }
        Ok(_) => panic!("AlignmentReader::open should fail for nonexistent file"),
    }
}

#[test]
fn test_open_invalid_bam() -> Result<()> {
    let mut tmp = NamedTempFile::with_suffix(".bam")?;
    tmp.write_all(b"this is not a valid BAM file")?;
    tmp.flush()?;

    assert!(
        AlignmentReader::open(tmp.path(), None, 0).is_err(),
        "opening garbage data as BAM should fail"
    );
    Ok(())
}

/// End-to-end round-trip: build a BAM with `SamBuilder`, read it back
/// through `AlignmentReader::fill_record`, and check that `BamRec`'s
/// cached scalars (flags, positions, mapq, cigar, quality scores)
/// survive the encode/decode round trip. This exercises
/// `BamRec::read_from` + `refresh_cache` — the "fast path" that the
/// `RikerRecord` design hinges on.
#[test]
fn test_bamrec_roundtrip_preserves_scalars() -> Result<()> {
    let mut builder = SamBuilder::new();
    builder.add(pair("pair1").at("chr1", 100, 200).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).reverse().nm(3));
    let bam = builder.to_temp_bam()?;

    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM"))?;

    let mut record = reader.empty_record();
    let mut count = 0;
    while reader.fill_record(&mut record)? {
        count += 1;
        // Every record should expose its flags + start (the records we
        // wrote are all mapped) and have non-empty cigar/quality views.
        assert!(record.flags().bits() != 0, "flags should be populated");
        assert!(record.alignment_start().is_some(), "alignment_start should be Some");
        assert!(record.mapping_quality().is_some(), "mapping_quality should be Some");
        assert!(record.cigar_len() > 0, "cigar should be non-empty");
        assert!(!record.quality_scores().is_empty(), "quality scores should be present");
        assert!(!record.sequence().is_empty(), "sequence should be populated when requested");

        // The fragment carries an NM tag; the paired reads do not. Check
        // that aux scanning surfaces NM only when present.
        if let Some(name) = record.name()
            && &**name == b"frag1"
        {
            let nm = record.aux_tag(*b"NM").and_then(|v| v.as_int()).unwrap();
            assert_eq!(nm, 3, "NM tag value should round-trip");
        }
    }
    assert_eq!(count, 3, "expected one fragment + one pair (=3 records)");

    Ok(())
}

/// Sanity-check that the BAM fast path's cached `alignment_end` matches
/// what noodles computes on the raw record. Also exercises the
/// `RikerRecord::Bam` variant of the cigar iterator.
#[test]
fn test_bamrec_alignment_end_matches_cigar_span() -> Result<()> {
    let mut builder = SamBuilder::new();
    // 50bp match → alignment_end = start + 50 - 1 = start + 49.
    builder.add(read().name("r1").at("chr1", 100).len(50).no_nm());
    let bam = builder.to_temp_bam()?;

    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE)?;
    let mut record = reader.empty_record();
    assert!(reader.fill_record(&mut record)?);

    let RikerRecord::Bam(_) = &record else {
        panic!("expected RikerRecord::Bam variant for BAM input");
    };
    let start = record.alignment_start().unwrap().get();
    let end = record.alignment_end().unwrap().get();
    assert_eq!(start, 100);
    assert_eq!(end, 149, "alignment_end = start + cigar reference span - 1");

    Ok(())
}

/// CRAM round-trip via `fill_record`: write a small CRAM via
/// `SamBuilder::to_temp_cram`, open it with `AlignmentReader::open`,
/// fill records in place into a `RikerRecord::Htslib` slot, and check
/// they survive the encode/decode with their flags/positions/cigar
/// intact. Exercises the rust-htslib backend.
#[test]
fn test_cram_roundtrip_via_fill_record() -> Result<()> {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add(pair("pair1").at("chr1", 100, 200).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).no_nm());
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()), 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE)?;

    let mut record = reader.empty_record();
    // Variant check runs before the first `fill_record`, so what we're
    // locking in is that `empty_record()` allocated the right variant
    // for this format — not the type of the filled records (those reuse
    // this same slot).
    let RikerRecord::Htslib(_) = &record else {
        panic!("expected RikerRecord::Htslib for CRAM input");
    };

    let mut count = 0;
    while reader.fill_record(&mut record)? {
        assert!(record.alignment_start().is_some(), "alignment_start should be Some");
        assert!(record.cigar_len() > 0, "cigar should be non-empty");
        count += 1;
    }
    assert_eq!(count, 3, "expected 1 pair + 1 fragment = 3 records");

    Ok(())
}

/// CRAM via `drive_collector_single_threaded`: confirms the CRAM
/// single-collector flow round-trips through the
/// `fill_record`-driven loop in [`drive_records`].
#[test]
fn test_cram_drives_through_collector() -> Result<()> {
    use riker_lib::collector::{Collector, drive_collector_single_threaded};
    use riker_lib::progress::ProgressLogger;

    /// Trivial collector that just counts records seen — enough to
    /// verify the CRAM iterator path drives `accept` with the right
    /// number of records.
    struct Counter {
        n: u64,
    }
    impl Collector for Counter {
        fn initialize(&mut self, _h: &noodles::sam::Header) -> Result<()> {
            Ok(())
        }
        fn accept(&mut self, _r: &RikerRecord, _h: &noodles::sam::Header) -> Result<()> {
            self.n += 1;
            Ok(())
        }
        fn finish(&mut self) -> Result<()> {
            Ok(())
        }
        fn name(&self) -> &'static str {
            "counter"
        }
        fn field_needs(&self) -> RikerRecordRequirements {
            RikerRecordRequirements::NONE
        }
    }

    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add(pair("p1").at("chr1", 100, 200).len(50).no_nm());
    builder.add(pair("p2").at("chr1", 300, 400).len(50).no_nm());
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()), 0)?;
    let mut collector = Counter { n: 0 };
    let mut progress = ProgressLogger::new("test", "reads", 1_000_000);
    drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)?;
    assert_eq!(collector.n, 4, "expected 4 records (2 pairs)");
    Ok(())
}

/// Regression test: when a collector's `initialize` fails, `finish` must not
/// be called. Real collectors populate `Option` fields (e.g. a sequence
/// dictionary) in `initialize` and unwrap them in `finish`; calling finish
/// on a partially-constructed collector turns an honest error into a panic.
#[test]
fn test_drive_collector_does_not_finish_after_failed_initialize() -> Result<()> {
    use anyhow::anyhow;
    use riker_lib::collector::{Collector, drive_collector_single_threaded};
    use riker_lib::progress::ProgressLogger;

    struct FailingInit {
        finish_called: bool,
    }
    impl Collector for FailingInit {
        fn initialize(&mut self, _h: &noodles::sam::Header) -> Result<()> {
            Err(anyhow!("synthetic initialize failure"))
        }
        fn accept(&mut self, _r: &RikerRecord, _h: &noodles::sam::Header) -> Result<()> {
            Ok(())
        }
        fn finish(&mut self) -> Result<()> {
            self.finish_called = true;
            panic!("finish must not be called when initialize failed");
        }
        fn name(&self) -> &'static str {
            "failing_init"
        }
        fn field_needs(&self) -> RikerRecordRequirements {
            RikerRecordRequirements::NONE
        }
    }

    let mut builder = SamBuilder::new();
    builder.add(read().name("r1").at("chr1", 100).len(50).no_nm());
    let bam = builder.to_temp_bam()?;
    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    let mut collector = FailingInit { finish_called: false };
    let mut progress = ProgressLogger::new("test", "reads", 1_000_000);

    let result = drive_collector_single_threaded(&mut reader, &mut collector, &mut progress);
    let err = result.expect_err("expected initialize error to propagate");
    assert!(
        err.to_string().contains("synthetic initialize failure"),
        "expected init error to surface, got: {err}"
    );
    assert!(!collector.finish_called, "finish must be skipped after init failure");
    Ok(())
}

/// CRAM round-trip exercising the opt-in decoder paths: with `with_sequence()`
/// the SIMD nibble decoder runs against htslib's packed `seq().encoded`,
/// and with `with_aux_tag(NM)` the htslib aux walker populates the store.
#[test]
fn test_cram_fill_record_decodes_sequence_and_aux() -> Result<()> {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    // Fragment with NM:i:3 — the only record carrying an aux tag.
    builder.add(read().name("frag1").at("chr1", 100).len(50).nm(3));
    builder.add(pair("pair1").at("chr1", 200, 400).len(50).no_nm());
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()), 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM"))?;
    let mut record = reader.empty_record();

    let mut count = 0;
    let mut frag_seen = false;
    while reader.fill_record(&mut record)? {
        count += 1;
        assert!(!record.sequence().is_empty(), "sequence should be populated when requested");
        assert_eq!(record.sequence().len(), record.sequence_len(), "decoded len matches l_seq");
        for &b in record.sequence() {
            assert!(b.is_ascii(), "decoded base should be ASCII (got {b:#x})");
        }

        if let Some(name) = record.name()
            && &**name == b"frag1"
        {
            frag_seen = true;
            let nm = record.aux_tag(*b"NM").and_then(|v| v.as_int()).unwrap();
            assert_eq!(nm, 3, "NM tag value should round-trip via htslib aux walker");
            assert!(record.has_aux_tag(*b"NM"));
        }
    }
    assert_eq!(count, 3, "expected 1 fragment + 1 pair = 3 records");
    assert!(frag_seen, "fragment with NM tag should be present");
    Ok(())
}

/// Sanity-check that the CRAM path's cached `alignment_end` matches
/// what the cigar reference span implies. Mirrors the BAM equivalent
/// at [`test_bamrec_alignment_end_matches_cigar_span`].
#[test]
fn test_htslibrec_alignment_end_matches_cigar_span() -> Result<()> {
    let refa = FastaBuilder::new().contig("chr1", vec![b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add(read().name("r1").at("chr1", 100).len(50).no_nm());
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()), 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE)?;
    let mut record = reader.empty_record();
    assert!(reader.fill_record(&mut record)?);

    let RikerRecord::Htslib(_) = &record else {
        panic!("expected RikerRecord::Htslib variant for CRAM input");
    };
    let start = record.alignment_start().unwrap().get();
    let end = record.alignment_end().unwrap().get();
    assert_eq!(start, 100);
    assert_eq!(end, 149, "alignment_end = start + cigar reference span - 1");
    // A 50bp single-match read has exactly one CIGAR op. Catches a regression
    // where `cigar_len` could drift between op count and byte count.
    assert_eq!(record.cigar_len(), 1, "cigar_len must report op count, not bytes");
    Ok(())
}

/// SAM round-trip via `fill_record`. Mirrors the BAM scalar-roundtrip
/// test but exercises the plain-text SAM in-place fill path
/// (`fill_sam_slot`) which the BAM and CRAM tests don't cover.
#[test]
fn test_sam_fill_record_roundtrips_scalars() -> Result<()> {
    let mut builder = SamBuilder::new();
    builder.add(pair("pair1").at("chr1", 100, 200).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).reverse().nm(3));
    let sam = builder.to_temp_sam()?;

    let mut reader = AlignmentReader::open(sam.path(), None, 0)?;

    reader.set_requirements(RikerRecordRequirements::NONE.with_aux_tag(*b"NM"))?;
    let mut record = reader.empty_record();
    let mut count = 0;
    while reader.fill_record(&mut record)? {
        count += 1;
        // Fallback variant on SAM input.
        let RikerRecord::Fallback(_) = &record else {
            panic!("expected RikerRecord::Fallback for SAM input");
        };
        assert!(record.alignment_start().is_some(), "alignment_start should be Some");
        assert!(record.cigar_len() > 0, "cigar should be non-empty");
        assert!(record.flags().bits() != 0, "flags should be populated");

        if let Some(name) = record.name()
            && &**name == b"frag1"
        {
            let nm = record.aux_tag(*b"NM").and_then(|v| v.as_int()).unwrap();
            assert_eq!(nm, 3, "NM tag value should round-trip through SAM");
        }
    }
    assert_eq!(count, 3, "expected one fragment + one pair (=3 records)");
    Ok(())
}

/// Reading a BAM single-threaded (`decode_threads == 0`, plain noodles BGZF)
/// and multithreaded (`decode_threads > 0`, the `bgzf` crate's `MultithreadedReader`)
/// must surface byte-identical records. This is the core correctness guarantee
/// for the parallel BAM decode path.
#[test]
fn bam_records_identical_across_decode_thread_counts() -> Result<()> {
    let mut builder = SamBuilder::new();
    builder.add(pair("pair1").at("chr1", 100, 300).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).reverse().nm(3));
    builder.add(pair("pair2").at("chr1", 700, 900).len(60).mapq(55).no_nm());
    let bam = builder.to_temp_bam()?;

    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");
    let single = collect_stream_rows(bam.path(), None, &requirements, 0)?;
    let multithreaded = collect_stream_rows(bam.path(), None, &requirements, 4)?;

    assert_eq!(single.len(), 5, "two pairs + one fragment = 5 records");
    assert_eq!(single, multithreaded, "multithreaded BAM decode must match single-threaded");
    Ok(())
}

/// Indexed BAM region queries must return byte-identical records whether
/// htslib's decode pool is sized to 0 (single-threaded) or several threads.
#[test]
fn indexed_bam_query_identical_across_decode_thread_counts() -> Result<()> {
    let mut builder = coord_builder(&[("chr1", 2_000)]);
    builder.add(pair("pair1").at("chr1", 100, 300).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).nm(3));
    builder.add(pair("pair2").at("chr1", 1_200, 1_400).len(60).mapq(55).no_nm());
    let bam = builder.to_temp_indexed_bam()?;

    let region: Region = "chr1".parse().expect("valid region");
    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");

    let collect = |decode_threads: usize| -> Result<Vec<Row>> {
        let mut reader = IndexedAlignmentReader::open(bam.path(), None, decode_threads)?;
        let mut rows = Vec::new();
        reader.query_for_each(&region, &requirements, |record| {
            rows.push(record_row(record));
            Ok(())
        })?;
        Ok(rows)
    };

    let single = collect(0)?;
    let multithreaded = collect(4)?;
    assert_eq!(single.len(), 5, "two pairs + one fragment overlap the query");
    assert_eq!(single, multithreaded, "multithreaded indexed query must match single-threaded");
    Ok(())
}

/// A `.csi`-indexed BAM (no `.bai` present) must query correctly — the `.csi`
/// selected by [`select_bam_index`] is handed to htslib, which reads it — at
/// both a zero and non-zero decode-thread count.
#[test]
fn indexed_bam_query_works_with_csi_index() -> Result<()> {
    let mut builder = coord_builder(&[("chr1", 2_000)]);
    builder.add(pair("pair1").at("chr1", 100, 300).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).nm(3));
    builder.add(pair("pair2").at("chr1", 1_200, 1_400).len(60).mapq(55).no_nm());
    let bam = builder.to_temp_csi_indexed_bam()?;

    let region: Region = "chr1".parse().expect("valid region");
    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");

    let collect = |decode_threads: usize| -> Result<Vec<Row>> {
        let mut reader = IndexedAlignmentReader::open(bam.path(), None, decode_threads)?;
        let mut rows = Vec::new();
        reader.query_for_each(&region, &requirements, |record| {
            rows.push(record_row(record));
            Ok(())
        })?;
        Ok(rows)
    };

    let single = collect(0)?;
    let multithreaded = collect(4)?;
    assert_eq!(single.len(), 5, "two pairs + one fragment overlap the query");
    assert_eq!(single, multithreaded, "CSI query must match across decode thread counts");
    Ok(())
}

/// SAM has no BGZF layer, so `decode_threads` is meaningless and must be
/// silently ignored — the same records come back regardless of the value.
#[test]
fn sam_ignores_decode_threads() -> Result<()> {
    let mut builder = SamBuilder::new();
    builder.add(pair("pair1").at("chr1", 100, 300).len(50).no_nm());
    builder.add(read().name("frag1").at("chr1", 500).len(75).mapq(30).reverse().nm(3));
    let sam = builder.to_temp_sam()?;

    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");
    let unset = collect_stream_rows(sam.path(), None, &requirements, 0)?;
    let with_threads = collect_stream_rows(sam.path(), None, &requirements, 4)?;

    assert_eq!(unset.len(), 3, "one pair + one fragment = 3 records");
    assert_eq!(unset, with_threads, "decode_threads must not affect SAM reads");
    Ok(())
}

/// With no requirements, a CRAM yields the same always-populated fields as the
/// BAM it was made from (flags, positions, CIGAR, mate fields, qualities),
/// although htslib skips decoding its sequence and aux tags.
#[test]
fn cram_without_requirements_matches_bam() -> Result<()> {
    let (bam_rows, cram_rows) = bam_and_cram_rows(&RikerRecordRequirements::NONE)?;
    assert_eq!(cram_rows, bam_rows);
    Ok(())
}

/// Skipping the sequence decode must keep each read's length, which `gcbias`,
/// `hybcap` and `rna` read without declaring the sequence.
#[test]
fn cram_without_requirements_keeps_sequence_length() -> Result<()> {
    let (_, cram_rows) = bam_and_cram_rows(&RikerRecordRequirements::NONE)?;
    let lengths: Vec<usize> = cram_rows.iter().map(|row| row.sequence_len).collect();
    // In coordinate order: pair1 (x2), reverse_with_mismatches, split_pair read 1,
    // clipped_with_indels, split_pair read 2, unmapped.
    assert_eq!(lengths, [50, 50, 60, 100, 82, 100, 100]);
    assert!(cram_rows.iter().all(|row| row.sequence.is_empty()), "bases were not requested");
    Ok(())
}

/// Requesting the sequence makes a CRAM decode the bases, matching the BAM.
#[test]
fn cram_with_sequence_requirement_matches_bam() -> Result<()> {
    let (bam_rows, cram_rows) = bam_and_cram_rows(&RikerRecordRequirements::NONE.with_sequence())?;
    assert!(cram_rows.iter().all(|row| row.sequence.len() == row.sequence_len));
    assert_eq!(cram_rows, bam_rows);
    Ok(())
}

/// htslib drops `NM` from a CRAM when it matches the reference and regenerates
/// it on read, which takes the sequence decode. Requesting `NM` alone, as
/// `gcbias` does, must still surface the BAM's values.
#[test]
fn cram_with_only_nm_requirement_regenerates_dropped_nm() -> Result<()> {
    let (fasta, bam, cram) = varied_bam_and_cram()?;
    let stored_nm = count_stored_cram_tags(cram.path(), fasta.path(), *b"NM")?;
    assert_eq!(stored_nm, 0, "precondition: NM is not stored");

    let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
    let bam_rows = collect_stream_rows(bam.path(), None, &requirements, 0)?;
    let cram_rows = collect_stream_rows(cram.path(), Some(fasta.path()), &requirements, 0)?;
    let nm = |rows: &[Row]| rows.iter().map(|row| row.nm).collect::<Vec<_>>();
    let has_mismatches = nm(&bam_rows).iter().any(|n| n.is_some_and(|n| n > 0));
    assert!(has_mismatches, "precondition: a read has mismatches");
    assert_eq!(nm(&cram_rows), nm(&bam_rows));
    Ok(())
}

/// Requesting an aux tag the CRAM stores verbatim (`MC`, as `rna` does) decodes
/// the tags without the bases, matching the BAM.
#[test]
fn cram_with_only_mc_requirement_matches_bam() -> Result<()> {
    let (bam_rows, cram_rows) =
        bam_and_cram_rows(&RikerRecordRequirements::NONE.with_aux_tag(*b"MC"))?;
    assert!(cram_rows.iter().any(|row| row.mc.is_some()), "precondition: the pair carries MC");
    assert!(cram_rows.iter().all(|row| row.sequence.is_empty()), "bases were not requested");
    assert_eq!(cram_rows, bam_rows);
    Ok(())
}

/// The requirements reach htslib before its decode threads read ahead, so a
/// threaded CRAM read spanning several containers surfaces the same records as
/// a single-threaded one.
#[test]
fn cram_records_identical_across_decode_thread_counts() -> Result<()> {
    let reference = FastaBuilder::hg38_mini();
    let fasta = reference.to_temp_fasta()?;
    let mut builder = SamBuilder::with_contigs(&reference.dict()).sort_order(SortOrder::Coordinate);
    // htslib starts a new container every 10,000 records, so 25,000 span three.
    for i in 0..12_500 {
        builder.add(pair(format!("p{i}")).at("chr1", 100 + i % 4_000, 400 + i % 4_000).len(50));
    }
    let cram = builder.to_temp_htslib_cram(fasta.path())?;

    let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"MC");
    let single = collect_stream_rows(cram.path(), Some(fasta.path()), &requirements, 0)?;
    let multithreaded = collect_stream_rows(cram.path(), Some(fasta.path()), &requirements, 4)?;
    assert_eq!(single.len(), 25_000);
    assert_eq!(single, multithreaded);
    Ok(())
}

/// The owned-record iterator decodes per the reader's requirements on CRAM too.
#[test]
fn cram_riker_records_matches_bam() -> Result<()> {
    let (fasta, bam, cram) = varied_bam_and_cram()?;
    let requirements = RikerRecordRequirements::NONE.with_sequence();
    let expected = collect_stream_rows(bam.path(), None, &requirements, 0)?;

    let mut reader = AlignmentReader::open(cram.path(), Some(fasta.path()), 0)?;
    reader.set_requirements(requirements)?;
    let mut rows: Vec<Row> = reader
        .riker_records()?
        .map(|record| record.map(|r| record_row(&r)))
        .collect::<Result<_>>()?;
    clear_unmapped_mapq(&mut rows);
    assert_eq!(rows, expected);
    Ok(())
}

/// Requirements are fixed once reading starts: htslib may already have decoded
/// CRAM containers ahead of the reader with the old ones.
#[test]
fn set_requirements_after_first_read_is_an_error() -> Result<()> {
    let bam = varied_records().1.to_temp_bam()?;
    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    reader.set_requirements(RikerRecordRequirements::NONE)?;
    let mut record = reader.empty_record();
    assert!(reader.fill_record(&mut record)?);

    let err = reader
        .set_requirements(RikerRecordRequirements::NONE.with_sequence())
        .expect_err("changing requirements mid-read must fail");
    assert!(err.to_string().contains("before the first record"), "unexpected error: {err}");
    Ok(())
}

/// A reader never told what to decode refuses to read, rather than handing
/// back records with no bases or aux tags.
#[test]
fn fill_record_without_requirements_is_an_error() -> Result<()> {
    let bam = varied_records().1.to_temp_bam()?;
    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    let mut record = reader.empty_record();
    let err = reader.fill_record(&mut record).expect_err("reading without requirements must fail");
    assert!(err.to_string().contains("set_requirements"), "unexpected error: {err}");
    Ok(())
}

#[test]
fn riker_records_without_requirements_is_an_error() -> Result<()> {
    let bam = varied_records().1.to_temp_bam()?;
    let mut reader = AlignmentReader::open(bam.path(), None, 0)?;
    let err = reader.riker_records().err().expect("reading without requirements must fail");
    assert!(err.to_string().contains("set_requirements"), "unexpected error: {err}");
    Ok(())
}
