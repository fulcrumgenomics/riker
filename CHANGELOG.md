# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Toolkit-wide `--threads` option** for multithreaded input decoding — the
  number of cores riker tries to saturate (it may spin a few more threads
  internally to do so). It is a single, global budget (accepted before or after
  the subcommand) that each command divides between decoding the input and its
  own work. The single-pass tools (`wgs`, `alignment`, …) run their own work on
  one thread and hand the rest to input decoding; `multi` lays the budget out
  across decode threads, a dispatch thread, and parallel collector workers,
  **tuned to the input format** (from benchmarking): BAM inflate is cheap so it
  favors collector workers, while CRAM decode is expensive and parallelizes
  well so it favors decode threads. BAM uses `noodles-bgzf`'s multithreaded BGZF
  reader; CRAM sizes htslib's decode pool. Gains are sub-linear and flatten
  past ~6 threads for BAM / ~8 for CRAM, but the sweet spot is platform- and
  input-dependent. Left unset, the single-pass tools run single-threaded and
  `multi` uses a small default pool.

### Changed

- **`multi --threads` is now the toolkit `--threads`.** It is reinterpreted as
  the total thread budget (which `multi` divides between input decoding and
  collector workers) rather than a bare worker count, and with no value given
  `multi` now uses up to four threads (clamped to the core count) instead of a
  fixed two. Existing `riker multi --threads N` invocations keep working.

### Fixed

- **`multi` output is now deterministic across thread counts.** Each collector
  is pinned to a single worker and receives batches in file order, so a
  parallel run is byte-identical to a single-threaded one; previously the
  per-record output order of some collectors could vary with the thread count.

## [0.3.0] - 2026-06-24

### Added

- **`gcbias --exclude-intervals`** (BED or IntervalList, auto-detected) to mask
  artifact regions — e.g. poly-G stretches or adapter constructs — from the
  GC-bias signal. Both the read and its reference window at an excluded position
  are dropped, keeping the numerator/denominator normalization consistent.
  Exclusion is keyed on each read's computed start position and applied per read.
  ([#31](https://github.com/fulcrumgenomics/riker/pull/31))

### Changed

- **`gcbias` summary schema:** read accounting is now an explicit, auditable
  funnel — `total_reads` → `aligned_reads` → `filtered_reads` /
  `frac_filtered_reads` — and secondary/QC-fail records are excluded outright
  (counted nowhere) rather than inflating the totals.
  ([#31](https://github.com/fulcrumgenomics/riker/pull/31))
- **`wgs`: replaced `--exclude-unpaired-reads` with `--include-unpaired-reads`.**
  The old flag was a boolean that defaulted to `true` with no way to negate it
  from the command line, so unpaired reads (and reads with an unmapped mate)
  were always excluded and the option could never actually be toggled. The new
  `--include-unpaired-reads` flag defaults to `false` (behavior unchanged) and,
  when passed, counts unpaired reads toward coverage. This is a breaking change
  to the CLI, but the removed flag had no working behavior to preserve.
  ([#34](https://github.com/fulcrumgenomics/riker/issues/34),
  [#35](https://github.com/fulcrumgenomics/riker/issues/35),
  [#36](https://github.com/fulcrumgenomics/riker/pull/36))

### Fixed

- Missing input files and other invalid inputs now produce a clear error
  instead of a panic. ([#29](https://github.com/fulcrumgenomics/riker/pull/29))
- The x86_64 multivers launcher exited 255 under Docker amd64 emulation (e.g.
  the biocontainer on Apple Silicon). Pinned `cargo-multivers >=0.12.0`, which
  carries the fexecve/memfd fix, in CI and the benchmark install script.
  ([#32](https://github.com/fulcrumgenomics/riker/issues/32),
  [#33](https://github.com/fulcrumgenomics/riker/pull/33))
- The usage line could read `Usage: 11 [OPTIONS] <COMMAND>` instead of
  `Usage: riker ...`. clap derives its displayed program name from `argv[0]`,
  and under the multivers fexecve/memfd launcher `argv[0]` is the descriptor
  path (`/proc/self/fd/11`), whose basename is the fd number. Pinned the
  displayed name with an explicit `bin_name = "riker"`.
  ([#38](https://github.com/fulcrumgenomics/riker/pull/38))

### Internal

- Added `cargo release` configuration for cutting workspace releases, plus this
  changelog and automated `[Unreleased]` → version stamping at release time.
  ([#30](https://github.com/fulcrumgenomics/riker/pull/30))

## [0.2.0] - 2026-05-08

A performance-focused release: several tools got substantially faster and CRAM
processing is now several times quicker than in 0.1.0.

### Performance

- `riker basic` is ~2.4× faster.
  ([#18](https://github.com/fulcrumgenomics/riker/pull/18))
- `riker wgs` is ~2× faster; the depth pipeline was restructured around a shared
  mate buffer. ([#10](https://github.com/fulcrumgenomics/riker/pull/10))
- CRAM decoding is several times faster, now via rust-htslib.
  ([#24](https://github.com/fulcrumgenomics/riker/pull/24))
- Added reusable byte-level SIMD kernels, adopted in `alignment`, `hybcap`, and
  `error`. ([#9](https://github.com/fulcrumgenomics/riker/pull/9))
- BAM/SAM reads reuse `RecordBuf` allocations; `multi` was rewritten around a new
  `RikerRecord` with pooled records and a crossbeam work queue.
  ([#4](https://github.com/fulcrumgenomics/riker/pull/4),
  [#12](https://github.com/fulcrumgenomics/riker/pull/12),
  [#13](https://github.com/fulcrumgenomics/riker/pull/13),
  [#15](https://github.com/fulcrumgenomics/riker/pull/15))
- Switched the flate2 backend to zlib-ng.
  ([#2](https://github.com/fulcrumgenomics/riker/pull/2))

### Added

- A `sample` column on all per-row metric outputs; `hybcap` per-base output is
  now gzipped. ([#21](https://github.com/fulcrumgenomics/riker/pull/21))
- `cargo-multivers` x86_64 release artifacts that dispatch to the best CPU
  variant at startup. ([#19](https://github.com/fulcrumgenomics/riker/pull/19))
- A Snakemake performance benchmark pipeline.
  ([#26](https://github.com/fulcrumgenomics/riker/pull/26))

### Changed

- `riker isize` trims the histogram TSV to the median ± `DEVIATIONS` × MAD.
  ([#22](https://github.com/fulcrumgenomics/riker/pull/22))
- Upgraded kuva to 0.2.0 and tightened plot tick/label layout.
  ([#25](https://github.com/fulcrumgenomics/riker/pull/25))
- Dropped the release `overflow-checks` override after a WGS-scale audit.
  ([#17](https://github.com/fulcrumgenomics/riker/pull/17))
- Bumped noodles (0.107 → 0.110), strum (0.27 → 0.28), and other dependencies.
  ([#3](https://github.com/fulcrumgenomics/riker/pull/3),
  [#27](https://github.com/fulcrumgenomics/riker/pull/27))

### Fixed

- Plot text failed to render when running in a bioconda biocontainer.
  ([#18](https://github.com/fulcrumgenomics/riker/pull/18))
- Some CRAM files triggered errors or panics.
  ([#24](https://github.com/fulcrumgenomics/riker/pull/24))
- Fixed cross-contig orphan processing in `riker error`.
  ([#11](https://github.com/fulcrumgenomics/riker/pull/11))
- Read `.sam.gz` files with `MultiGzDecoder` so multi-member streams decode
  fully. ([#2](https://github.com/fulcrumgenomics/riker/pull/2))

## [0.1.0] - 2026-04-08

First (alpha) release. A fast Rust CLI toolkit that ports key QC-metrics tools
from Picard/fgbio, processing SAM/BAM/CRAM files and emitting clean TSV.

### Added

- `alignment` — alignment summary metrics (cf. Picard CollectAlignmentSummaryMetrics).
- `basic` — base distribution, mean quality, and quality score distribution
  (cf. Picard CollectBaseDistributionByCycle / MeanQualityByCycle /
  QualityScoreDistribution).
- `error` — base-level error metrics: mismatch, overlap, indel
  (cf. Picard CollectSamErrorMetrics).
- `isize` — insert size distribution metrics (cf. Picard CollectInsertSizeMetrics).
- `wgs` — whole-genome coverage metrics (cf. Picard CollectWgsMetrics).
- `gcbias` — GC bias metrics (cf. Picard CollectGcBiasMetrics).
- `hybcap` — hybrid capture (bait/target) metrics (cf. Picard CollectHsMetrics).
- `multi` — run multiple collectors in a single BAM pass
  (cf. Picard CollectMultipleMetrics).
- `docs` — print metric field documentation.
