use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::Result;

/// How a command spends its total thread budget: some threads decode input,
/// the rest do parallel compute. Produced by [`Command::plan_threads`].
#[derive(Debug, Clone, Copy)]
pub struct ThreadPlan {
    /// BGZF/CRAM decode workers handed to the reader, in addition to the
    /// calling thread. `0` reads single-threaded.
    pub decode_threads: usize,
    /// Parallel compute workers the command runs. `1` for the single-pass
    /// tools (their compute is serial on the main thread); `multi` uses more.
    pub compute_workers: usize,
}

/// Trait implemented by every riker subcommand.
#[enum_dispatch::enum_dispatch]
pub trait Command {
    /// Execute the command.
    ///
    /// `threads` is the toolkit-wide `--threads` value (`None` when unset) —
    /// the total OS-thread budget, counting the main thread. Each command
    /// resolves it against [`default_threads`](Command::default_threads) and
    /// divides it via [`plan_threads`](Command::plan_threads).
    ///
    /// # Errors
    /// Returns an error if the command fails.
    fn execute(&self, threads: Option<u8>) -> Result<()>;

    /// Total thread budget to use when `--threads` is unset. The single-pass
    /// tools stay serial (`1`); `multi` overrides this to parallelize by
    /// default.
    fn default_threads(&self) -> NonZero<usize> {
        NonZero::<usize>::MIN
    }

    /// Divide the resolved total budget between input decoding and compute.
    /// The default hands the whole remainder to the reader because a
    /// single-pass tool's own work runs on the calling thread; `multi`
    /// overrides this with a format-aware split.
    fn plan_threads(&self, total: NonZero<usize>, _input: &Path) -> ThreadPlan {
        ThreadPlan { decode_threads: total.get() - 1, compute_workers: 1 }
    }
}

/// Resolve the toolkit-wide `--threads` value against a command's
/// [`default_threads`](Command::default_threads). The CLI value parser
/// guarantees any provided value is `>= 1`.
#[must_use]
pub fn resolve_threads(threads: Option<u8>, default: NonZero<usize>) -> NonZero<usize> {
    threads.and_then(|n| NonZero::new(usize::from(n))).unwrap_or(default)
}

/// Build an output path by appending a suffix to a prefix path.
///
/// E.g. `output_path(Path::new("out/sample"), ".metrics.txt")` → `out/sample.metrics.txt`
#[must_use]
pub fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.to_string_lossy()))
}
