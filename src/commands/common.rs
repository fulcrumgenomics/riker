use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

pub const CORE_OPTIONS: &str = "Core Options";
pub const EXTENDED_OPTIONS: &str = "Extended Options";

/// Validate at startup that output files can be written for the given `-o`/`--output` prefix.
///
/// Commands derive every output path by appending a suffix to the prefix (see
/// [`super::command::output_path`]), so the prefix's parent directory must exist and be
/// writable. This is checked by creating and removing a probe file there. Calling this before
/// the (potentially long) processing pass means an unwritable destination fails immediately
/// rather than after the whole input has been read.
///
/// # Errors
/// Returns an error if the prefix's parent directory does not exist, is not a directory, or
/// cannot be written to.
pub fn validate_output_prefix(prefix: &Path) -> Result<()> {
    // `prefix.parent()` is `Some("")` for a bare filename (e.g. `out`) — treat that as the current
    // directory — and `None` when the prefix is itself a root (`/`, `C:\`), where the directory to
    // check is that root.
    let dir = match prefix.parent() {
        Some(p) if p.as_os_str().is_empty() => Path::new("."),
        Some(p) => p,
        None => prefix,
    };

    let meta = std::fs::metadata(dir).with_context(|| {
        format!("output directory '{}' does not exist or is inaccessible", dir.display())
    })?;
    if !meta.is_dir() {
        bail!("output prefix parent '{}' is not a directory", dir.display());
    }

    // Permission bits alone are unreliable across platforms; probe by actually creating a file.
    let probe = dir.join(format!(".riker-write-probe.{}", std::process::id()));
    std::fs::File::create(&probe)
        .with_context(|| format!("output directory '{}' is not writable", dir.display()))?;
    let _ = std::fs::remove_file(&probe);
    Ok(())
}

/// clap value parser that ensures a CLI argument refers to an existing,
/// non-directory path.
///
/// Used by every CLI option that names an input file (BAM, FASTA, intervals, VCF, etc.)
/// so that typos and missing files are reported with a clean error at argument-parse
/// time, before any tool starts opening readers or writing output.
///
/// `std::fs::metadata` follows symlinks, so symlinks to regular files are accepted
/// and broken symlinks fail the existence check. Directories are rejected so that
/// downstream code does not produce confusing "is a directory" errors. FIFOs and
/// character devices like `/dev/stdin` are intentionally allowed through — if the
/// downstream reader cannot consume them (e.g. because it needs to seek), it will
/// surface its own specific error.
///
/// # Errors
/// Returns `Err` if the path does not exist (including broken symlinks) or refers
/// to a directory.
pub fn parse_existing_file(s: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    let meta = std::fs::metadata(&path)
        .map_err(|e| format!("cannot access file '{}': {e}", path.display()))?;
    if meta.is_dir() {
        return Err(format!("'{}' is a directory, expected a file", path.display()));
    }
    Ok(path)
}

/// Options for specifying the input alignment file.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct InputOptions {
    /// Input SAM, BAM, or CRAM file.
    #[arg(short = 'i', long, value_name = "INPUT", value_parser = parse_existing_file)]
    pub input: PathBuf,
}

/// Options for specifying the output path prefix.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct OutputOptions {
    /// Output path prefix.
    #[arg(short = 'o', long, value_name = "PREFIX")]
    pub output: PathBuf,
}

/// Options for specifying a required reference FASTA.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct ReferenceOptions {
    /// Reference FASTA file (must be indexed with .fai). Required.
    #[arg(short = 'r', long, value_name = "FASTA", value_parser = parse_existing_file)]
    pub reference: PathBuf,
}

/// Options for specifying an optional reference FASTA.
#[derive(clap::Args, Debug, Clone)]
#[command()]
pub struct OptionalReferenceOptions {
    /// Reference FASTA file (must be indexed with .fai). Optional unless CRAM input supplied.
    #[arg(short = 'r', long, value_name = "FASTA", value_parser = parse_existing_file)]
    pub reference: Option<PathBuf>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_output_prefix_accepts_writable_dir() {
        let dir = tempfile::tempdir().unwrap();
        // Prefix is a stem inside an existing, writable directory.
        let prefix = dir.path().join("sample");
        assert!(validate_output_prefix(&prefix).is_ok());
        // The probe file must not be left behind.
        let leftovers: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
        assert!(leftovers.is_empty(), "probe file should be cleaned up");
    }

    #[test]
    fn validate_output_prefix_rejects_missing_dir() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("nonexistent_subdir").join("sample");
        let err = validate_output_prefix(&prefix).unwrap_err().to_string();
        assert!(err.contains("does not exist") || err.contains("inaccessible"), "got: {err}");
    }

    #[test]
    fn validate_output_prefix_rejects_file_as_parent() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("a_file");
        std::fs::write(&file, b"x").unwrap();
        // Using the file as if it were a directory for the prefix.
        let prefix = file.join("sample");
        assert!(validate_output_prefix(&prefix).is_err());
    }

    #[test]
    fn parse_existing_file_accepts_regular_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.txt");
        std::fs::write(&path, b"hello").unwrap();

        let parsed = parse_existing_file(path.to_str().unwrap()).unwrap();
        assert_eq!(parsed, path);
    }

    #[test]
    fn parse_existing_file_rejects_missing_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("does_not_exist.txt");

        let err = parse_existing_file(path.to_str().unwrap()).unwrap_err();
        assert!(err.contains("cannot access file"), "unexpected error: {err}");
        assert!(err.contains("does_not_exist.txt"), "unexpected error: {err}");
    }

    #[test]
    fn parse_existing_file_rejects_directory() {
        let dir = tempfile::tempdir().unwrap();

        let err = parse_existing_file(dir.path().to_str().unwrap()).unwrap_err();
        assert!(err.contains("is a directory"), "unexpected error: {err}");
    }

    #[test]
    fn parse_existing_file_follows_symlinks_to_regular_files() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("target.txt");
        std::fs::write(&target, b"data").unwrap();
        let link = dir.path().join("link");
        std::os::unix::fs::symlink(&target, &link).unwrap();

        let parsed = parse_existing_file(link.to_str().unwrap()).unwrap();
        assert_eq!(parsed, link);
    }

    #[test]
    fn parse_existing_file_rejects_broken_symlink() {
        let dir = tempfile::tempdir().unwrap();
        let link = dir.path().join("dangling");
        std::os::unix::fs::symlink(dir.path().join("nope"), &link).unwrap();

        let err = parse_existing_file(link.to_str().unwrap()).unwrap_err();
        assert!(err.contains("cannot access file"), "unexpected error: {err}");
    }

    #[test]
    fn parse_existing_file_accepts_fifo() {
        // FIFOs are not regular files but they are not directories either, so
        // we let them through for tools that can consume streams.
        let dir = tempfile::tempdir().unwrap();
        let fifo = dir.path().join("p");
        let status = std::process::Command::new("mkfifo")
            .arg(&fifo)
            .status()
            .expect("mkfifo should be available on the test host");
        assert!(status.success());

        let parsed = parse_existing_file(fifo.to_str().unwrap()).unwrap();
        assert_eq!(parsed, fifo);
    }
}
