#[allow(dead_code)]
mod helpers;

use std::io::Write;
use std::path::Path;

use anyhow::Result;
use helpers::SamBuilder;
use riker_lib::sam::alignment_reader::AlignmentReader;
use riker_lib::sam::riker_record::{RikerRecord, RikerRecordRequirements};
use tempfile::NamedTempFile;

#[test]
fn test_open_nonexistent_file() {
    match AlignmentReader::open(Path::new("/no/such/file.bam"), None) {
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
        AlignmentReader::open(tmp.path(), None).is_err(),
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
    builder.add_pair("pair1", 0, 100, 200, 200, 60, 50, false, false);
    builder.add_unpaired("frag1", 0, 500, 30, 75, true, false, false, Some(3));
    let bam = builder.to_temp_bam()?;

    let mut reader = AlignmentReader::open(bam.path(), None)?;
    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");

    let mut record = reader.empty_record();
    let mut count = 0;
    while reader.fill_record(&requirements, &mut record)? {
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
    builder.add_unpaired("r1", 0, 100, 60, 50, false, false, false, None);
    let bam = builder.to_temp_bam()?;

    let mut reader = AlignmentReader::open(bam.path(), None)?;
    let mut record = reader.empty_record();
    assert!(reader.fill_record(&RikerRecordRequirements::NONE, &mut record)?);

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
    use helpers::{FastaBuilder, coord_builder};

    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add_pair("pair1", 0, 100, 200, 200, 60, 50, false, false);
    builder.add_unpaired("frag1", 0, 500, 30, 75, false, false, false, None);
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()))?;

    let requirements = RikerRecordRequirements::NONE;
    let mut record = reader.empty_record();
    // Variant check runs before the first `fill_record`, so what we're
    // locking in is that `empty_record()` allocated the right variant
    // for this format — not the type of the filled records (those reuse
    // this same slot).
    let RikerRecord::Htslib(_) = &record else {
        panic!("expected RikerRecord::Htslib for CRAM input");
    };

    let mut count = 0;
    while reader.fill_record(&requirements, &mut record)? {
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
    use helpers::{FastaBuilder, coord_builder};
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

    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add_pair("p1", 0, 100, 200, 200, 60, 50, false, false);
    builder.add_pair("p2", 0, 300, 400, 200, 60, 50, false, false);
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()))?;
    let mut collector = Counter { n: 0 };
    let mut progress = ProgressLogger::new("test", "reads", 1_000_000);
    drive_collector_single_threaded(&mut reader, &mut collector, &mut progress)?;
    assert_eq!(collector.n, 4, "expected 4 records (2 pairs)");
    Ok(())
}

/// CRAM round-trip exercising the opt-in decoder paths: with `with_sequence()`
/// the SIMD nibble decoder runs against htslib's packed `seq().encoded`,
/// and with `with_aux_tag(NM)` the htslib aux walker populates the store.
#[test]
fn test_cram_fill_record_decodes_sequence_and_aux() -> Result<()> {
    use helpers::{FastaBuilder, coord_builder};

    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    // Fragment with NM:i:3 — the only record carrying an aux tag.
    builder.add_unpaired("frag1", 0, 100, 60, 50, false, false, false, Some(3));
    builder.add_pair("pair1", 0, 200, 400, 200, 60, 50, false, false);
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()))?;
    let requirements = RikerRecordRequirements::NONE.with_sequence().with_aux_tag(*b"NM");
    let mut record = reader.empty_record();

    let mut count = 0;
    let mut frag_seen = false;
    while reader.fill_record(&requirements, &mut record)? {
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
    use helpers::{FastaBuilder, coord_builder};

    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 1_000]).to_temp_fasta()?;
    let mut builder = coord_builder(&[("chr1", 1_000)]);
    builder.add_unpaired("r1", 0, 100, 60, 50, false, false, false, None);
    let cram = builder.to_temp_cram(refa.path())?;

    let mut reader = AlignmentReader::open(cram.path(), Some(refa.path()))?;
    let mut record = reader.empty_record();
    assert!(reader.fill_record(&RikerRecordRequirements::NONE, &mut record)?);

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
    builder.add_pair("pair1", 0, 100, 200, 200, 60, 50, false, false);
    builder.add_unpaired("frag1", 0, 500, 30, 75, true, false, false, Some(3));
    let sam = builder.to_temp_sam()?;

    let mut reader = AlignmentReader::open(sam.path(), None)?;

    let requirements = RikerRecordRequirements::NONE.with_aux_tag(*b"NM");
    let mut record = reader.empty_record();
    let mut count = 0;
    while reader.fill_record(&requirements, &mut record)? {
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
