use std::collections::HashMap;
use std::ops::Index;

use noodles::fasta;
use noodles::sam::Header;

/// Metadata for a single reference sequence (contig).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceMetadata {
    /// 0-based positional index matching the BAM ref_id.
    index: usize,
    /// Contig name.
    name: String,
    /// Contig length in bases.
    length: usize,
}

impl SequenceMetadata {
    /// Return the 0-based index of this sequence.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }

    /// Return the name of this sequence.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the length of this sequence in bases.
    #[must_use]
    pub fn length(&self) -> usize {
        self.length
    }
}

/// A lookup table mapping reference sequence names to their indices and lengths.
///
/// Constructed from a BAM header or a FASTA index, this type consolidates
/// all name↔index↔length mapping into a single shared structure.
#[derive(Debug, Clone)]
pub struct SequenceDictionary {
    sequences: Vec<SequenceMetadata>,
    name_to_index: HashMap<String, usize>,
}

impl SequenceDictionary {
    /// Return the number of sequences.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Return `true` if the dictionary contains no sequences.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Look up a sequence by its 0-based index.
    #[must_use]
    pub fn get_by_index(&self, index: usize) -> Option<&SequenceMetadata> {
        self.sequences.get(index)
    }

    /// Look up a sequence by name.
    #[must_use]
    pub fn get_by_name(&self, name: &str) -> Option<&SequenceMetadata> {
        self.name_to_index.get(name).map(|&i| &self.sequences[i])
    }

    /// Iterate over all sequences in index order.
    pub fn iter(&self) -> impl Iterator<Item = &SequenceMetadata> {
        self.sequences.iter()
    }

    /// Return contig names in index order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.sequences.iter().map(|s| s.name.as_str()).collect()
    }
}

impl Index<usize> for SequenceDictionary {
    type Output = SequenceMetadata;

    fn index(&self, index: usize) -> &SequenceMetadata {
        &self.sequences[index]
    }
}

impl Index<&str> for SequenceDictionary {
    type Output = SequenceMetadata;

    fn index(&self, name: &str) -> &SequenceMetadata {
        let i = self.name_to_index[name];
        &self.sequences[i]
    }
}

impl From<&Header> for SequenceDictionary {
    fn from(header: &Header) -> Self {
        let mut sequences = Vec::new();
        let mut name_to_index = HashMap::new();

        for (i, (name, rs)) in header.reference_sequences().iter().enumerate() {
            let name = String::from_utf8_lossy(name.as_ref()).into_owned();
            let length = usize::from(rs.length());
            name_to_index.insert(name.clone(), i);
            sequences.push(SequenceMetadata { index: i, name, length });
        }

        Self { sequences, name_to_index }
    }
}

impl From<&fasta::fai::Index> for SequenceDictionary {
    fn from(index: &fasta::fai::Index) -> Self {
        let mut sequences = Vec::new();
        let mut name_to_index = HashMap::new();

        for (i, record) in index.as_ref().iter().enumerate() {
            let name = String::from_utf8_lossy(record.name().as_ref()).to_string();
            #[allow(clippy::cast_possible_truncation)]
            let length = record.length() as usize;
            name_to_index.insert(name.clone(), i);
            sequences.push(SequenceMetadata { index: i, name, length });
        }

        Self { sequences, name_to_index }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::num::NonZeroUsize;

    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    /// Build a header with the given contig names and lengths.
    fn make_header(contigs: &[(&str, usize)]) -> Header {
        let mut builder = Header::builder();
        for &(name, len) in contigs {
            let ref_seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(len).expect("len > 0"));
            builder = builder.add_reference_sequence(name.as_bytes(), ref_seq);
        }
        builder.build()
    }

    #[test]
    fn test_from_header() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000), ("chrX", 500)]);
        let dict = SequenceDictionary::from(&header);

        assert_eq!(dict.len(), 3);
        assert!(!dict.is_empty());

        assert_eq!(dict[0].name(), "chr1");
        assert_eq!(dict[0].length(), 1000);
        assert_eq!(dict[0].index(), 0);

        assert_eq!(dict[1].name(), "chr2");
        assert_eq!(dict[1].length(), 2000);

        assert_eq!(dict[2].name(), "chrX");
        assert_eq!(dict[2].length(), 500);
    }

    #[test]
    fn test_from_empty_header() {
        let header = Header::default();
        let dict = SequenceDictionary::from(&header);
        assert!(dict.is_empty());
        assert_eq!(dict.len(), 0);
    }

    #[test]
    fn test_get_by_index() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000)]);
        let dict = SequenceDictionary::from(&header);

        let meta = dict.get_by_index(0).unwrap();
        assert_eq!(meta.name(), "chr1");
        assert_eq!(meta.length(), 1000);

        let meta = dict.get_by_index(1).unwrap();
        assert_eq!(meta.name(), "chr2");
        assert_eq!(meta.length(), 2000);

        assert!(dict.get_by_index(2).is_none());
    }

    #[test]
    fn test_get_by_name() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000)]);
        let dict = SequenceDictionary::from(&header);

        let meta = dict.get_by_name("chr2").unwrap();
        assert_eq!(meta.index(), 1);
        assert_eq!(meta.length(), 2000);

        assert!(dict.get_by_name("chrZ").is_none());
    }

    #[test]
    fn test_index_by_usize() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000)]);
        let dict = SequenceDictionary::from(&header);
        assert_eq!(dict[1].name(), "chr2");
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn test_index_by_usize_out_of_bounds() {
        let header = make_header(&[("chr1", 1000)]);
        let dict = SequenceDictionary::from(&header);
        let _ = &dict[5];
    }

    #[test]
    fn test_index_by_str() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000)]);
        let dict = SequenceDictionary::from(&header);
        assert_eq!(dict["chr1"].length(), 1000);
    }

    #[test]
    #[should_panic(expected = "no entry found for key")]
    fn test_index_by_str_unknown() {
        let header = make_header(&[("chr1", 1000)]);
        let dict = SequenceDictionary::from(&header);
        let _ = &dict["nope"];
    }

    #[test]
    fn test_names() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000), ("chrX", 500)]);
        let dict = SequenceDictionary::from(&header);
        assert_eq!(dict.names(), vec!["chr1", "chr2", "chrX"]);
    }

    #[test]
    fn test_iter() {
        let header = make_header(&[("chr1", 1000), ("chr2", 2000)]);
        let dict = SequenceDictionary::from(&header);
        let names: Vec<&str> = dict.iter().map(SequenceMetadata::name).collect();
        assert_eq!(names, vec!["chr1", "chr2"]);
    }
}
