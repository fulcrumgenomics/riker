pub mod alignment_reader;
pub mod mate_buffer;
pub mod pair_orientation;
pub mod riker_record;

use std::path::Path;

use noodles::sam::Header;
use noodles::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};

/// True if the header declares coordinate sort order (`@HD SO:coordinate`).
pub(crate) fn is_coordinate_sorted(header: &Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| AsRef::<[u8]>::as_ref(so) == COORDINATE)
}

/// Extract the sample name from the first read group's SM tag, or fall back to
/// the BAM file stem.
#[must_use]
pub(crate) fn derive_sample(input: &Path, header: &Header) -> String {
    use noodles::sam::header::record::value::map::read_group::tag;
    header
        .read_groups()
        .values()
        .find_map(|rg| rg.other_fields().get(&tag::SAMPLE).map(ToString::to_string))
        .unwrap_or_else(|| {
            input.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string()
        })
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use noodles::sam::Header;

    use super::derive_sample;
    use crate::test_support::SamBuilder;

    #[test]
    fn derive_sample_from_read_group() {
        let header = SamBuilder::new().read_group("rg0", "SampleA").header().clone();
        let input = PathBuf::from("/data/sample.bam");
        assert_eq!(derive_sample(&input, &header), "SampleA");
    }

    #[test]
    fn derive_sample_falls_back_to_filename_without_sm() {
        // Read group present but no SM tag → falls back to file stem.
        let header = SamBuilder::new().read_group_without_sample("rg0").header().clone();
        let input = PathBuf::from("/data/my_file.bam");
        assert_eq!(derive_sample(&input, &header), "my_file");
    }

    #[test]
    fn derive_sample_falls_back_to_filename_without_read_groups() {
        let header = Header::default();
        let input = PathBuf::from("/data/my_file.bam");
        assert_eq!(derive_sample(&input, &header), "my_file");
    }

    #[test]
    fn derive_sample_unknown_when_no_stem() {
        let header = Header::default();
        let input = PathBuf::from("/");
        assert_eq!(derive_sample(&input, &header), "unknown");
    }
}
