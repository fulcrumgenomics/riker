use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, bail};
use bitvec::prelude::*;
use noodles::bcf;
use noodles::core::Region;
use noodles::vcf;
use noodles_bgzf as bgzf;

use crate::intervals::Intervals;
use crate::sequence_dict::SequenceDictionary;

/// An indexed variant reader supporting both VCF (bgzip-compressed) and BCF formats.
///
/// Automatically detects the format from the file extension (`.bcf` for BCF,
/// anything else for VCF) and validates the appropriate index exists.
pub struct IndexedVcf {
    inner: VariantReaderInner,
}

/// Internal enum dispatching between VCF and BCF indexed readers.
enum VariantReaderInner {
    Vcf { reader: vcf::io::IndexedReader<bgzf::io::Reader<File>>, header: vcf::Header },
    Bcf { reader: bcf::io::IndexedReader<bgzf::io::Reader<File>>, header: vcf::Header },
}

impl IndexedVcf {
    /// Open an indexed VCF or BCF file.
    ///
    /// Format is detected from the file extension: `.bcf` opens a BCF reader
    /// (requires `.csi` index), anything else opens a VCF reader (requires
    /// `.tbi` or `.csi` index).
    ///
    /// # Errors
    /// Returns an error if the file or its index cannot be opened.
    pub fn from_path(path: &Path) -> Result<Self> {
        let is_bcf = path.extension().is_some_and(|ext| ext.eq_ignore_ascii_case("bcf"));

        if is_bcf { Self::open_bcf(path) } else { Self::open_vcf(path) }
    }

    /// Open an indexed VCF file (bgzip-compressed, .tbi or .csi index).
    fn open_vcf(path: &Path) -> Result<Self> {
        let tbi_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".tbi");
            std::path::PathBuf::from(p)
        };
        let csi_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".csi");
            std::path::PathBuf::from(p)
        };

        if !tbi_path.exists() && !csi_path.exists() {
            bail!(
                "VCF index not found. Expected one of:\n  {}\n  {}\n\
                 Run `bcftools index -t {}` to create one.",
                tbi_path.display(),
                csi_path.display(),
                path.display(),
            );
        }

        let mut reader = vcf::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed VCF: {}", path.display()))?;

        let header = reader
            .read_header()
            .with_context(|| format!("Failed to read VCF header from: {}", path.display()))?;

        Ok(Self { inner: VariantReaderInner::Vcf { reader, header } })
    }

    /// Open an indexed BCF file (.csi index).
    fn open_bcf(path: &Path) -> Result<Self> {
        let csi_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".csi");
            std::path::PathBuf::from(p)
        };

        if !csi_path.exists() {
            bail!(
                "BCF index not found: expected {}\n\
                 Run `bcftools index {}` to create one.",
                csi_path.display(),
                path.display(),
            );
        }

        let mut reader = bcf::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed BCF: {}", path.display()))?;

        let header = reader
            .read_header()
            .with_context(|| format!("Failed to read BCF header from: {}", path.display()))?;

        Ok(Self { inner: VariantReaderInner::Bcf { reader, header } })
    }

    /// Load variant positions for a region into a `BitVec`.
    ///
    /// Returns a `BitVec` of length `end - start` where set bits indicate
    /// positions with known variants (relative to the region start).
    /// `start` and `end` are 0-based half-open coordinates.
    ///
    /// # Errors
    /// Returns an error if the query fails.
    pub fn load_region(&mut self, contig: &str, start: u32, end: u32) -> Result<BitVec> {
        let len = (end - start) as usize;
        let mut bits = bitvec![0; len];

        let region: Region =
            format!("{contig}:{}-{end}", start + 1).parse().with_context(|| {
                format!("Failed to parse query region: {contig}:{}-{end}", start + 1)
            })?;

        match &mut self.inner {
            VariantReaderInner::Vcf { reader, header } => {
                let query = reader
                    .query(header, &region)
                    .with_context(|| format!("Failed to query VCF for region: {region}"))?;
                for result in query.records() {
                    let record: vcf::Record =
                        result.context("Failed to read VCF record during query")?;
                    collect_variant_positions(&record, header, start, end, &mut bits)?;
                }
            }
            VariantReaderInner::Bcf { reader, header } => {
                let query = reader
                    .query(header, &region)
                    .with_context(|| format!("Failed to query BCF for region: {region}"))?;
                for result in query.records() {
                    let record: bcf::Record =
                        result.context("Failed to read BCF record during query")?;
                    collect_variant_positions(&record, header, start, end, &mut bits)?;
                }
            }
        }

        Ok(bits)
    }
}

/// Mask the full span of a variant record (start through end, inclusive).
/// For SNPs this is a single position; for deletions and multi-base variants
/// it covers all affected reference positions.
/// Works for both VCF and BCF records via the `vcf::variant::Record` trait.
#[expect(
    clippy::cast_possible_truncation,
    reason = "VCF/BCF positions fit in u32 for genomic coordinates"
)]
fn collect_variant_positions(
    record: &impl vcf::variant::Record,
    header: &vcf::Header,
    region_start: u32,
    region_end: u32,
    bits: &mut BitVec,
) -> Result<()> {
    let Some(pos_result) = record.variant_start() else { return Ok(()) };
    let var_start_0based =
        (pos_result.context("Failed to parse variant position")?.get() - 1) as u32;
    let var_end_0based =
        (record.variant_end(header).context("Failed to determine variant end")?.get() - 1) as u32;

    // Clamp to the region and set all bits in the variant's span
    let mask_start = var_start_0based.max(region_start);
    let mask_end = var_end_0based.min(region_end - 1);
    if mask_start <= mask_end {
        let lo = (mask_start - region_start) as usize;
        let hi = (mask_end - region_start) as usize;
        bits[lo..=hi].fill(true);
    }
    Ok(())
}

/// Pre-load variant masks for all relevant contigs from an indexed VCF or BCF.
///
/// Returns a map from contig name to a `BitVec` covering the full contig length,
/// where set bits indicate known variant positions. When `intervals` is provided,
/// contigs with no intervals are skipped.
///
/// # Panics
/// Panics if the sequence dictionary contains an index that cannot be resolved
/// (should never happen since iteration stays within `dict.len()`).
///
/// # Errors
/// Returns an error if any query fails.
#[expect(
    clippy::cast_possible_truncation,
    reason = "contig lengths fit in u32 for genomic coordinates"
)]
pub fn load_variant_masks(
    vcf: &mut IndexedVcf,
    dict: &SequenceDictionary,
    intervals: Option<&Intervals>,
) -> Result<std::collections::HashMap<String, BitVec>> {
    let mut masks = std::collections::HashMap::new();
    for ref_id in 0..dict.len() {
        // Skip contigs with no intervals when intervals are specified
        if let Some(ivs) = intervals
            && !ivs.has_contig(ref_id)
        {
            continue;
        }
        let meta = dict.get_by_index(ref_id).expect("ref_id in range");
        let bits = vcf.load_region(meta.name(), 0, meta.length() as u32)?;
        masks.insert(meta.name().to_string(), bits);
    }
    Ok(masks)
}
