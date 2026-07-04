//! A fluent builder for small, sites-only VCFs used in tests.
//!
//! Mirrors fgbio's `VcfBuilder`, adapted to riker's Rust builder idiom: declare
//! the INFO/FILTER header lines up front, then add variants fluently. Genotypes
//! and FORMAT/sample columns are intentionally out of scope — the commands that
//! consume a VCF (currently `error`'s known-sites masking) read only the sites.
//!
//! The output is a bgzip-compressed `.vcf.gz` with a sibling tabix `.tbi`, which
//! is what [`crate::vcf::IndexedVcf`] requires for its region queries.
//!
//! ```ignore
//! let mut vcf = VcfBuilder::new()
//!     .info_header("DP", InfoType::Integer, "Total depth")
//!     .filter_header("LowAD", "Low allele depth");
//! vcf.add("chr1", 105, "C").alt("T"); // SNP: masks position 105
//! vcf.add("chr1", 200, "ACG").alt("A").filter("LowAD").info("DP", 30); // deletion: masks 200-202
//! let indexed = vcf.to_temp_vcf()?; // .vcf.gz + .tbi
//! ```

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::core::Position;
use noodles::tabix;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::header::record::value::map::{Contig, Filter, Info as InfoMap};
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::record_buf::info::field::Value;
use noodles::vcf::variant::record_buf::{AlternateBases, Filters, Ids, Info};
use noodles::vcf::{self, Header};
use tempfile::NamedTempFile;

use super::fasta_builder::FastaBuilder;

/// Builds a sites-only VCF and writes it as bgzip + tabix index. Construct via
/// [`VcfBuilder::new`] (frozen `chr1`/`chr2` contigs) or [`VcfBuilder::with_contigs`].
pub struct VcfBuilder {
    contigs: Vec<(String, usize)>,
    info_headers: Vec<(String, InfoType, String)>,
    filter_headers: Vec<(String, String)>,
    variants: Vec<VariantEntry>,
}

impl VcfBuilder {
    /// A builder whose header declares the frozen default `chr1`/`chr2` contigs
    /// (see [`FastaBuilder::hg38_mini`]), matching what `read()` derives bases from.
    #[must_use]
    pub fn new() -> Self {
        Self::from_contigs(FastaBuilder::hg38_mini().dict())
    }

    /// A builder whose header declares the given `(name, length)` contigs.
    #[must_use]
    pub fn with_contigs(contigs: &[(&str, usize)]) -> Self {
        Self::from_contigs(contigs.iter().map(|(name, len)| ((*name).to_string(), *len)).collect())
    }

    fn from_contigs(contigs: Vec<(String, usize)>) -> Self {
        Self { contigs, info_headers: Vec::new(), filter_headers: Vec::new(), variants: Vec::new() }
    }

    /// Declare an `##INFO` header line (`Number=1`). Any INFO key set on a variant
    /// must be declared here first, or [`Self::write_vcf`] panics.
    #[must_use]
    pub fn info_header(
        mut self,
        id: impl Into<String>,
        kind: InfoType,
        description: impl Into<String>,
    ) -> Self {
        self.info_headers.push((id.into(), kind, description.into()));
        self
    }

    /// Declare an `##FILTER` header line. Any filter set on a variant must be
    /// declared here first (`PASS` is always valid), or [`Self::write_vcf`] panics.
    #[must_use]
    pub fn filter_header(mut self, id: impl Into<String>, description: impl Into<String>) -> Self {
        self.filter_headers.push((id.into(), description.into()));
        self
    }

    /// Add a variant at 1-based `pos` with reference allele `reference`, returning
    /// a handle to set its alternate alleles, ID, QUAL, filters, and INFO.
    ///
    /// The masked reference span is `reference.len()` positions starting at `pos`,
    /// so a single-base `reference` is a SNP and a multi-base one is a deletion.
    #[expect(clippy::missing_panics_doc, reason = "last_mut on a just-pushed Vec is always Some")]
    pub fn add(
        &mut self,
        contig: impl Into<String>,
        pos: usize,
        reference: impl Into<String>,
    ) -> &mut VariantEntry {
        self.variants.push(VariantEntry {
            contig: contig.into(),
            pos,
            reference: reference.into(),
            alts: Vec::new(),
            id: None,
            qual: None,
            filters: Vec::new(),
            infos: Vec::new(),
        });
        self.variants.last_mut().expect("just pushed")
    }

    /// Write the VCF (bgzip) to `path` plus its tabix index at `<path>.tbi`.
    ///
    /// # Errors
    /// Returns an error if the file or index cannot be written.
    ///
    /// # Panics
    /// Panics if a variant references a FILTER or INFO key that was not declared
    /// via [`Self::filter_header`] / [`Self::info_header`], or an unknown contig.
    pub fn write_vcf(&self, path: &Path) -> io::Result<()> {
        self.validate();
        let header = self.build_header();

        let mut ordered: Vec<&VariantEntry> = self.variants.iter().collect();
        ordered.sort_by_key(|variant| (self.contig_index(&variant.contig), variant.pos));

        let file = File::create(path)?;
        let mut writer = vcf::io::Writer::new(bgzf::io::Writer::new(file));
        writer.write_header(&header)?;
        for variant in ordered {
            writer.write_variant_record(&header, &variant.to_record())?;
        }
        // Finalize the BGZF stream (EOF block) before it is re-read for indexing.
        writer.into_inner().finish()?;

        let index = vcf::fs::index(path)?;
        tabix::fs::write(with_suffix(path, ".tbi"), &index)?;
        Ok(())
    }

    /// Write to a temporary `.vcf.gz` (plus its `.tbi`) and return the handle.
    ///
    /// # Errors
    /// Returns an error if the temp file or its index cannot be written.
    pub fn to_temp_vcf(&self) -> io::Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".vcf.gz")?;
        self.write_vcf(tmp.path())?;
        Ok(tmp)
    }

    /// Build the noodles VCF header from the declared contigs and INFO/FILTER lines.
    fn build_header(&self) -> Header {
        let mut builder = Header::builder();
        for (name, length) in &self.contigs {
            let mut contig = Map::<Contig>::new();
            *contig.length_mut() = Some(*length);
            builder = builder.add_contig(name.as_str(), contig);
        }
        for (id, kind, description) in &self.info_headers {
            let map =
                Map::<InfoMap>::new(Number::Count(1), kind.to_noodles(), description.as_str());
            builder = builder.add_info(id.as_str(), map);
        }
        for (id, description) in &self.filter_headers {
            builder = builder.add_filter(id.as_str(), Map::<Filter>::new(description.as_str()));
        }
        builder.build()
    }

    /// Panic if any variant references an undeclared contig, FILTER, or INFO key.
    fn validate(&self) {
        for variant in &self.variants {
            assert!(
                self.contigs.iter().any(|(name, _)| name == &variant.contig),
                "variant on undeclared contig {:?}",
                variant.contig
            );
            for filter in &variant.filters {
                assert!(
                    filter == "PASS" || self.filter_headers.iter().any(|(id, _)| id == filter),
                    "filter {filter:?} used but no FILTER header declared"
                );
            }
            for (key, _) in &variant.infos {
                assert!(
                    self.info_headers.iter().any(|(id, _, _)| id == key),
                    "INFO key {key:?} used but no INFO header declared"
                );
            }
        }
    }

    /// The index of `contig` within the declared contig order (for coordinate sort).
    fn contig_index(&self, contig: &str) -> usize {
        self.contigs.iter().position(|(name, _)| name == contig).unwrap_or(usize::MAX)
    }
}

impl Default for VcfBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// A single variant being built. Obtained from [`VcfBuilder::add`]; every setter
/// returns `&mut Self` so calls chain within one statement.
pub struct VariantEntry {
    contig: String,
    pos: usize,
    reference: String,
    alts: Vec<String>,
    id: Option<String>,
    qual: Option<f32>,
    filters: Vec<String>,
    infos: Vec<(String, InfoValue)>,
}

impl VariantEntry {
    /// Add an alternate allele (repeatable for multi-allelic sites).
    pub fn alt(&mut self, alt: impl Into<String>) -> &mut Self {
        self.alts.push(alt.into());
        self
    }

    /// Set the variant ID (the `ID` column; default `.`).
    pub fn id(&mut self, id: impl Into<String>) -> &mut Self {
        self.id = Some(id.into());
        self
    }

    /// Set the QUAL column.
    pub fn qual(&mut self, qual: f32) -> &mut Self {
        self.qual = Some(qual);
        self
    }

    /// Add a FILTER value (repeatable). Must be declared via
    /// [`VcfBuilder::filter_header`] (or be `PASS`).
    pub fn filter(&mut self, filter: impl Into<String>) -> &mut Self {
        self.filters.push(filter.into());
        self
    }

    /// Set an INFO key/value (repeatable). The key must be declared via
    /// [`VcfBuilder::info_header`].
    pub fn info(&mut self, key: impl Into<String>, value: impl Into<InfoValue>) -> &mut Self {
        self.infos.push((key.into(), value.into()));
        self
    }

    /// Convert to a noodles record for writing.
    fn to_record(&self) -> RecordBuf {
        let mut builder = RecordBuf::builder()
            .set_reference_sequence_name(self.contig.clone())
            .set_variant_start(Position::new(self.pos).expect("position must be >= 1"))
            .set_reference_bases(self.reference.clone());
        if !self.alts.is_empty() {
            builder = builder.set_alternate_bases(AlternateBases::from(self.alts.clone()));
        }
        if let Some(id) = &self.id {
            builder = builder.set_ids([id.clone()].into_iter().collect::<Ids>());
        }
        if let Some(qual) = self.qual {
            builder = builder.set_quality_score(qual);
        }
        if !self.filters.is_empty() {
            builder = builder.set_filters(self.filters.iter().cloned().collect::<Filters>());
        }
        if !self.infos.is_empty() {
            let info: Info = self
                .infos
                .iter()
                .map(|(key, value)| (key.clone(), Some(value.to_noodles())))
                .collect();
            builder = builder.set_info(info);
        }
        builder.build()
    }
}

/// The type of an `##INFO` header line, mapped to the noodles equivalent.
#[derive(Clone, Copy)]
pub enum InfoType {
    Integer,
    Float,
    String,
}

impl InfoType {
    fn to_noodles(self) -> Type {
        match self {
            InfoType::Integer => Type::Integer,
            InfoType::Float => Type::Float,
            InfoType::String => Type::String,
        }
    }
}

/// A scalar INFO value. Construct via `From` (`30`, `0.5`, `"text"`).
pub enum InfoValue {
    Integer(i32),
    Float(f32),
    String(String),
}

impl InfoValue {
    fn to_noodles(&self) -> Value {
        match self {
            InfoValue::Integer(value) => Value::Integer(*value),
            InfoValue::Float(value) => Value::Float(*value),
            InfoValue::String(value) => Value::String(value.clone()),
        }
    }
}

impl From<i32> for InfoValue {
    fn from(value: i32) -> Self {
        InfoValue::Integer(value)
    }
}

impl From<f32> for InfoValue {
    fn from(value: f32) -> Self {
        InfoValue::Float(value)
    }
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "test INFO values are small literals; f32 is the VCF Float type"
)]
impl From<f64> for InfoValue {
    fn from(value: f64) -> Self {
        InfoValue::Float(value as f32)
    }
}

impl From<&str> for InfoValue {
    fn from(value: &str) -> Self {
        InfoValue::String(value.to_string())
    }
}

impl From<String> for InfoValue {
    fn from(value: String) -> Self {
        InfoValue::String(value)
    }
}

/// Append `suffix` to `path` (e.g. `foo.vcf.gz` + `.tbi` → `foo.vcf.gz.tbi`).
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut string = path.as_os_str().to_owned();
    string.push(suffix);
    PathBuf::from(string)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::IndexedVcf;

    /// A SNP sets exactly its single reference position in the mask.
    #[test]
    fn snp_masks_a_single_position() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
        vcf.add("chr1", 105, "C").alt("T");
        let tmp = vcf.to_temp_vcf().unwrap();

        let mut indexed = IndexedVcf::from_path(tmp.path()).unwrap();
        let mask = indexed.load_region("chr1", 0, 1000).unwrap();

        assert_eq!(mask.count_ones(), 1);
        assert!(mask[104], "0-based 104 (1-based 105) should be masked");
    }

    /// A deletion masks every reference position it spans (ref length positions).
    #[test]
    fn deletion_masks_its_whole_reference_span() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
        vcf.add("chr1", 200, "ACG").alt("A");
        let tmp = vcf.to_temp_vcf().unwrap();

        let mut indexed = IndexedVcf::from_path(tmp.path()).unwrap();
        let mask = indexed.load_region("chr1", 0, 1000).unwrap();

        // REF "ACG" at 1-based 200 spans 0-based [199, 202).
        assert_eq!(mask.count_ones(), 3);
        assert!(mask[199] && mask[200] && mask[201]);
    }

    /// Multiple variants across contigs each land at their own position.
    #[test]
    fn multiple_variants_are_written_and_indexed() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 500), ("chr2", 500)]);
        vcf.add("chr1", 10, "A").alt("G");
        vcf.add("chr1", 400, "T").alt("C");
        vcf.add("chr2", 50, "G").alt("A").alt("T"); // multi-allelic
        let tmp = vcf.to_temp_vcf().unwrap();

        let mut indexed = IndexedVcf::from_path(tmp.path()).unwrap();
        let chr1 = indexed.load_region("chr1", 0, 500).unwrap();
        let chr2 = indexed.load_region("chr2", 0, 500).unwrap();

        assert!(chr1[9] && chr1[399]);
        assert_eq!(chr1.count_ones(), 2);
        assert!(chr2[49]);
        assert_eq!(chr2.count_ones(), 1);
    }

    /// Declared FILTER and INFO fields write a valid, re-readable VCF.
    #[test]
    fn filters_and_info_round_trip_through_the_reader() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)])
            .info_header("DP", InfoType::Integer, "Total depth")
            .filter_header("LowAD", "Low allele depth");
        vcf.add("chr1", 300, "G").alt("A").filter("LowAD").info("DP", 42).qual(60.0);
        let tmp = vcf.to_temp_vcf().unwrap();

        let mut indexed = IndexedVcf::from_path(tmp.path()).unwrap();
        let mask = indexed.load_region("chr1", 0, 1000).unwrap();
        assert!(mask[299]);
    }

    /// Using a filter with no declared header is a builder misuse and panics.
    #[test]
    #[should_panic(expected = "no FILTER header declared")]
    fn undeclared_filter_panics() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
        vcf.add("chr1", 105, "C").alt("T").filter("LowAD");
        let _ = vcf.to_temp_vcf();
    }

    /// Using an INFO key with no declared header is a builder misuse and panics.
    #[test]
    #[should_panic(expected = "no INFO header declared")]
    fn undeclared_info_panics() {
        let mut vcf = VcfBuilder::with_contigs(&[("chr1", 1000)]);
        vcf.add("chr1", 105, "C").alt("T").info("DP", 10);
        let _ = vcf.to_temp_vcf();
    }
}
