use std::path::Path;

use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;

/// Extract the sample name from the first read group's SM tag, or fall back to
/// the BAM file stem.
#[must_use]
pub fn derive_sample(input: &Path, header: &Header) -> String {
    use noodles::sam::header::record::value::map::read_group::tag;
    header
        .read_groups()
        .values()
        .find_map(|rg| rg.other_fields().get(&tag::SAMPLE).map(ToString::to_string))
        .unwrap_or_else(|| {
            input.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string()
        })
}

/// Extract an integer value from an auxiliary BAM tag, returning `None` if the
/// tag is absent.  Handles all signed/unsigned integer variants.
#[allow(clippy::cast_sign_loss)]
#[must_use]
pub fn get_integer_tag(record: &RecordBuf, tag: [u8; 2]) -> Option<u32> {
    match record.data().get(&tag) {
        Some(DataValue::UInt8(v)) => Some(u32::from(*v)),
        Some(DataValue::UInt16(v)) => Some(u32::from(*v)),
        Some(DataValue::UInt32(v)) => Some(*v),
        Some(DataValue::Int8(v)) => Some((*v).max(0) as u32),
        Some(DataValue::Int16(v)) => Some((*v).max(0) as u32),
        Some(DataValue::Int32(v)) => Some((*v).max(0) as u32),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use noodles::sam::header::record::value::map::read_group::tag;
    use noodles::sam::header::record::value::{Map, map::ReadGroup};

    // ── derive_sample ────────────────────────────────────────────────────────

    #[test]
    fn test_derive_sample_from_read_group() {
        let mut rg = Map::<ReadGroup>::default();
        rg.other_fields_mut().insert(tag::SAMPLE, "SampleA".into());
        let header = Header::builder().add_read_group("rg0", rg).build();

        let input = PathBuf::from("/data/sample.bam");
        assert_eq!(derive_sample(&input, &header), "SampleA");
    }

    #[test]
    fn test_derive_sample_fallback_to_filename() {
        // Read group present but no SM tag → falls back to file stem
        let rg = Map::<ReadGroup>::default();
        let header = Header::builder().add_read_group("rg0", rg).build();

        let input = PathBuf::from("/data/my_file.bam");
        assert_eq!(derive_sample(&input, &header), "my_file");
    }

    #[test]
    fn test_derive_sample_no_read_groups() {
        let header = Header::default();
        let input = PathBuf::from("/data/my_file.bam");
        assert_eq!(derive_sample(&input, &header), "my_file");
    }

    #[test]
    fn test_derive_sample_unknown_fallback() {
        // Path with no meaningful file stem (e.g. root path)
        let header = Header::default();
        let input = PathBuf::from("/");
        assert_eq!(derive_sample(&input, &header), "unknown");
    }

    // ── get_integer_tag ──────────────────────────────────────────────────────

    fn record_with_tag(tag_name: [u8; 2], value: DataValue) -> RecordBuf {
        let mut data = noodles::sam::alignment::record_buf::Data::default();
        data.insert(tag_name.into(), value);
        RecordBuf::builder().set_data(data).build()
    }

    #[test]
    fn test_get_integer_tag_uint8() {
        let rec = record_with_tag(*b"NM", DataValue::UInt8(42));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(42));
    }

    #[test]
    fn test_get_integer_tag_uint16() {
        let rec = record_with_tag(*b"NM", DataValue::UInt16(1000));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(1000));
    }

    #[test]
    fn test_get_integer_tag_uint32() {
        let rec = record_with_tag(*b"NM", DataValue::UInt32(100_000));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(100_000));
    }

    #[test]
    fn test_get_integer_tag_int8_positive() {
        let rec = record_with_tag(*b"NM", DataValue::Int8(50));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(50));
    }

    #[test]
    fn test_get_integer_tag_int8_negative_clamped() {
        let rec = record_with_tag(*b"NM", DataValue::Int8(-5));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(0));
    }

    #[test]
    fn test_get_integer_tag_int16_negative_clamped() {
        let rec = record_with_tag(*b"NM", DataValue::Int16(-100));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(0));
    }

    #[test]
    fn test_get_integer_tag_int32_negative_clamped() {
        let rec = record_with_tag(*b"NM", DataValue::Int32(-1));
        assert_eq!(get_integer_tag(&rec, *b"NM"), Some(0));
    }

    #[test]
    fn test_get_integer_tag_missing() {
        let rec = RecordBuf::default();
        assert_eq!(get_integer_tag(&rec, *b"NM"), None);
    }

    #[test]
    fn test_get_integer_tag_non_integer() {
        let rec = record_with_tag(*b"BC", DataValue::String("ACGT".into()));
        assert_eq!(get_integer_tag(&rec, *b"BC"), None);
    }
}
