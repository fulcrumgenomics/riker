//! Small filesystem-path helpers shared across the I/O modules.

use std::path::{Path, PathBuf};

/// Returns `path` with `suffix` (e.g. `".bai"`) appended *after* the existing
/// extension — used to derive sibling index-file paths (`.bai` / `.crai` /
/// `.fai` / `.tbi` / `.csi`). Unlike [`Path::with_extension`], this appends
/// rather than replacing the extension.
#[must_use]
pub(crate) fn append_extension(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn appends_after_the_existing_extension() {
        assert_eq!(append_extension(Path::new("a/b.bam"), ".bai"), PathBuf::from("a/b.bam.bai"));
    }

    #[test]
    fn appends_when_there_is_no_extension() {
        assert_eq!(append_extension(Path::new("reads"), ".fai"), PathBuf::from("reads.fai"));
    }
}
