use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::Result;

/// How a command spends its total thread budget: some threads decode input,
/// the rest do parallel compute. Produced by [`Command::thread_plan`].
#[derive(Debug, Clone, Copy)]
pub struct ThreadPlan {
    /// BGZF/CRAM decode workers handed to the reader, in addition to the
    /// calling thread. `0` reads single-threaded.
    pub decode_threads: usize,
    /// Parallel compute workers the command runs. `multi` interprets `0` as
    /// serial-main mode (collectors run on the main thread with decode
    /// offloaded to the reader pool) and `>= 1` as its dispatch+worker
    /// pipeline. The single-pass tools ignore this (their compute is always
    /// serial on the main thread) and only read `decode_threads`.
    pub compute_workers: usize,
}

/// Trait implemented by every riker subcommand.
pub trait Command {
    /// Execute the command.
    ///
    /// `threads` is the toolkit-wide `--threads` value (`None` when unset) —
    /// the total OS-thread budget, counting the main thread. Each command
    /// resolves it against [`default_threads`](Command::default_threads) and
    /// divides it via [`thread_plan`](Command::thread_plan).
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

    /// Resolve `--threads` against [`default_threads`](Self::default_threads)
    /// and split it into decode and compute threads. The default gives a
    /// single-pass tool the whole remainder for input decoding — its own work
    /// runs on the calling thread, so `compute_workers` stays `1`. A tool that
    /// wants to parallelize its own work (or `multi`) overrides this to claim
    /// `compute_workers`; an overriding tool can read whatever it needs
    /// (input path, options) from `self`.
    fn thread_plan(&self, threads: Option<u8>) -> ThreadPlan {
        let total = resolve_threads(threads, self.default_threads()).get();
        ThreadPlan { decode_threads: total - 1, compute_workers: 1 }
    }
}

/// Resolve the toolkit-wide `--threads` value against a command's
/// [`default_threads`](Command::default_threads). The CLI value parser
/// guarantees any provided value is `>= 1`.
#[must_use]
pub(crate) fn resolve_threads(threads: Option<u8>, default: NonZero<usize>) -> NonZero<usize> {
    threads.and_then(|n| NonZero::new(usize::from(n))).unwrap_or(default)
}

/// Build an output path by appending a suffix to a prefix path.
///
/// E.g. `output_path(Path::new("out/sample"), ".metrics.txt")` → `out/sample.metrics.txt`
#[must_use]
pub fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.to_string_lossy()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nz(n: usize) -> NonZero<usize> {
        NonZero::new(n).expect("test uses non-zero")
    }

    #[test]
    fn resolve_threads_unset_uses_the_default() {
        assert_eq!(resolve_threads(None, nz(1)).get(), 1);
        assert_eq!(resolve_threads(None, nz(7)).get(), 7);
    }

    #[test]
    fn resolve_threads_uses_the_provided_value() {
        assert_eq!(resolve_threads(Some(1), nz(9)).get(), 1);
        assert_eq!(resolve_threads(Some(6), nz(9)).get(), 6);
    }

    /// A minimal command exercising the trait defaults.
    struct Dummy;
    impl Command for Dummy {
        fn execute(&self, _threads: Option<u8>) -> Result<()> {
            Ok(())
        }
    }

    #[test]
    fn default_command_is_single_threaded() {
        assert_eq!(Dummy.default_threads().get(), 1);
    }

    #[test]
    fn default_thread_plan_hands_the_whole_remainder_to_decoding() {
        // Single-pass tools compute on the calling thread and decode with the
        // rest, so decode_threads = total - 1 and compute stays at 1.
        // Dummy's default budget is 1, so an unset `--threads` decodes serially.
        let plan = Dummy.thread_plan(None);
        assert_eq!((plan.decode_threads, plan.compute_workers), (0, 1));
        let plan = Dummy.thread_plan(Some(4));
        assert_eq!((plan.decode_threads, plan.compute_workers), (3, 1));
    }
}
