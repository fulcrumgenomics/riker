use std::path::PathBuf;

pub const CORE_OPTIONS: &str = "Core Options";
pub const EXTENDED_OPTIONS: &str = "Extended Options";

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
