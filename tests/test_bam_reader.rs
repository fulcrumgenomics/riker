#[allow(dead_code)]
mod helpers;

use std::io::Write;
use std::path::Path;

use anyhow::Result;
use riker_lib::sam::alignment_reader::AlignmentReader;
use tempfile::NamedTempFile;

#[test]
fn test_open_nonexistent_file() {
    match AlignmentReader::new(Path::new("/no/such/file.bam"), None) {
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("/no/such/file.bam"), "error should contain path, got: {msg}");
        }
        Ok(_) => panic!("AlignmentReader::new should fail for nonexistent file"),
    }
}

#[test]
fn test_open_invalid_bam() -> Result<()> {
    let mut tmp = NamedTempFile::with_suffix(".bam")?;
    tmp.write_all(b"this is not a valid BAM file")?;
    tmp.flush()?;

    assert!(
        AlignmentReader::new(tmp.path(), None).is_err(),
        "opening garbage data as BAM should fail"
    );
    Ok(())
}
