use std::path::{Path, PathBuf};

use anyhow::Result;

/// Trait implemented by every riker subcommand.
#[enum_dispatch::enum_dispatch]
pub trait Command {
    /// Execute the command.
    ///
    /// # Errors
    /// Returns an error if the command fails.
    fn execute(&self) -> Result<()>;
}

/// Build an output path by appending a suffix to a prefix path.
///
/// E.g. `output_path(Path::new("out/sample"), ".metrics.txt")` → `out/sample.metrics.txt`
#[must_use]
pub fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.to_string_lossy()))
}
