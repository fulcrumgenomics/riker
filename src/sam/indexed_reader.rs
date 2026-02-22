use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles::core::Region;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, cram};
use noodles_bgzf as bgzf;

use super::alignment_reader::build_repository;

/// An indexed alignment reader supporting region-based queries for BAM and CRAM files.
pub enum IndexedAlignmentReader {
    Bam { reader: bam::io::IndexedReader<bgzf::io::Reader<File>>, header: Header },
    Cram { reader: cram::io::Reader<File>, index: cram::crai::Index, header: Header },
}

impl IndexedAlignmentReader {
    /// Open an indexed alignment file (BAM with `.bai`/`.csi`, or CRAM with `.crai`).
    ///
    /// A reference FASTA is required for CRAM files; it is ignored for BAM.
    ///
    /// # Errors
    /// Returns an error if the file or its index cannot be opened, or if a CRAM file
    /// is opened without a reference.
    pub fn new(path: &Path, reference: Option<&Path>) -> Result<Self> {
        Self::validate_index_exists(path)?;

        match (path.extension().is_some_and(|ext| ext.eq_ignore_ascii_case("cram")), reference) {
            (true, Some(ref_path)) => Self::open_indexed_cram(path, ref_path),
            (true, None) => {
                bail!("CRAM files require a reference FASTA (--reference): {}", path.display());
            }
            (false, _) => Self::open_indexed_bam(path),
        }
    }

    /// Returns the alignment header.
    #[must_use]
    pub fn header(&self) -> &Header {
        match self {
            Self::Bam { header, .. } | Self::Cram { header, .. } => header,
        }
    }

    /// Query records overlapping the given region, returning owned `RecordBuf`s.
    ///
    /// # Errors
    /// Returns an error if the query or record conversion fails.
    pub fn query(&mut self, region: &Region) -> Result<Vec<RecordBuf>> {
        match self {
            Self::Bam { reader, header } => {
                let query = reader
                    .query(header, region)
                    .with_context(|| format!("Failed to query BAM for region: {region}"))?;

                let mut records = Vec::new();
                for result in query.records() {
                    let record = result.context("Failed to read BAM record during query")?;
                    let record_buf = RecordBuf::try_from_alignment_record(header, &record)
                        .context("Failed to convert BAM record to RecordBuf")?;
                    records.push(record_buf);
                }
                Ok(records)
            }
            Self::Cram { reader, index, header } => {
                let query = reader
                    .query(header, index, region)
                    .with_context(|| format!("Failed to query CRAM for region: {region}"))?;

                let mut records = Vec::new();
                for result in query {
                    let record = result.context("Failed to read CRAM record during query")?;
                    let record_buf = RecordBuf::try_from_alignment_record(header, &record)
                        .context("Failed to convert CRAM record to RecordBuf")?;
                    records.push(record_buf);
                }
                Ok(records)
            }
        }
    }

    fn open_indexed_bam(path: &Path) -> Result<Self> {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed BAM: {}", path.display()))?;

        let header = reader
            .read_header()
            .with_context(|| format!("Failed to read header from: {}", path.display()))?;

        Ok(Self::Bam { reader, header })
    }

    fn open_indexed_cram(path: &Path, ref_path: &Path) -> Result<Self> {
        let repo = build_repository(ref_path)?;

        let crai_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".crai");
            std::path::PathBuf::from(p)
        };

        let index = cram::crai::fs::read(&crai_path)
            .with_context(|| format!("Failed to read CRAM index: {}", crai_path.display()))?;

        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed CRAM: {}", path.display()))?;

        let header = reader
            .read_header()
            .with_context(|| format!("Failed to read header from: {}", path.display()))?;

        Ok(Self::Cram { reader, index, header })
    }

    /// Validate that an index file exists for the given alignment file.
    fn validate_index_exists(path: &Path) -> Result<()> {
        // BAM index paths
        let bai_path = path.with_extension("bam.bai");
        let bai_path_alt = {
            let mut p = path.as_os_str().to_owned();
            p.push(".bai");
            std::path::PathBuf::from(p)
        };
        let csi_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".csi");
            std::path::PathBuf::from(p)
        };

        // CRAM index path
        let crai_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".crai");
            std::path::PathBuf::from(p)
        };

        if bai_path.exists() || bai_path_alt.exists() || csi_path.exists() || crai_path.exists() {
            return Ok(());
        }

        bail!(
            "Index not found for alignment file. Expected one of:\n  \
             {}\n  {}\n  {}\n  {}\n\
             Run `samtools index {}` to create one.",
            bai_path.display(),
            bai_path_alt.display(),
            csi_path.display(),
            crai_path.display(),
            path.display(),
        );
    }
}
