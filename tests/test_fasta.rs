#[allow(dead_code)]
mod helpers;

use std::path::Path;

use anyhow::Result;
use helpers::{FastaBuilder, SamBuilder};
use riker_lib::fasta::Fasta;

// ─── from_path error ─────────────────────────────────────────────────────────

#[test]
fn test_from_path_nonexistent() {
    match Fasta::from_path(Path::new("/no/such/reference.fa")) {
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("/no/such/reference.fa"), "error should contain path, got: {msg}");
        }
        Ok(_) => panic!("Fasta::from_path should fail for nonexistent file"),
    }
}

// ─── validate_bam_header ─────────────────────────────────────────────────────

#[test]
fn test_validate_bam_header_mismatch() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let fasta = Fasta::from_path(refa.path())?;

    // Build a BAM header with a contig not in the FASTA
    let builder = SamBuilder::with_contigs(&[("chr1".to_string(), 20), ("chrZ".to_string(), 100)]);

    match fasta.validate_bam_header(builder.header()) {
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("chrZ"), "error should mention missing contig, got: {msg}");
        }
        Ok(()) => panic!("validate_bam_header should fail when contig is missing from FASTA"),
    }
    Ok(())
}

#[test]
fn test_validate_bam_header_ok() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let fasta = Fasta::from_path(refa.path())?;

    let builder = SamBuilder::with_contigs(&[("chr1".to_string(), 20)]);
    assert!(fasta.validate_bam_header(builder.header()).is_ok());
    Ok(())
}

// ─── load_contig ─────────────────────────────────────────────────────────────

#[test]
fn test_load_contig_unknown() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", &[b'A'; 20]).to_temp_fasta()?;
    let mut fasta = Fasta::from_path(refa.path())?;

    match fasta.load_contig("chrZ", true) {
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("chrZ"), "error should mention the contig name, got: {msg}");
        }
        Ok(_) => panic!("load_contig should fail for unknown contig"),
    }
    Ok(())
}

#[test]
fn test_load_contig_ok() -> Result<()> {
    let refa = FastaBuilder::new().add_contig("chr1", b"ACGTacgt").to_temp_fasta()?;
    let mut fasta = Fasta::from_path(refa.path())?;

    let seq = fasta.load_contig("chr1", true)?;
    assert_eq!(seq, b"ACGTACGT"); // uppercase=true
    Ok(())
}

// ─── contig_length ───────────────────────────────────────────────────────────

#[test]
fn test_contig_length() -> Result<()> {
    let refa = FastaBuilder::new()
        .add_contig("chr1", &[b'A'; 100])
        .add_contig("chr2", &[b'G'; 50])
        .to_temp_fasta()?;
    let fasta = Fasta::from_path(refa.path())?;

    assert_eq!(fasta.contig_length("chr1"), Some(100));
    assert_eq!(fasta.contig_length("chr2"), Some(50));
    assert_eq!(fasta.contig_length("chrZ"), None);
    Ok(())
}

// ─── contig_names ────────────────────────────────────────────────────────────

#[test]
fn test_contig_names() -> Result<()> {
    let refa = FastaBuilder::new()
        .add_contig("chr1", &[b'A'; 20])
        .add_contig("chr2", &[b'G'; 30])
        .to_temp_fasta()?;
    let fasta = Fasta::from_path(refa.path())?;

    let names = fasta.contig_names();
    assert_eq!(names, vec!["chr1", "chr2"]);
    Ok(())
}
