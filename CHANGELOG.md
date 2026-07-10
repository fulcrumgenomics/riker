# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.1] - 2026-07-10

### Changed

- **Upgraded the `kuva` plotting library to 0.4.0 and removed the plot-layout workarounds its fixes make unnecessary.** The insert-size charts (`isize` and the `rna` transcript-space insert size) no longer round the x-axis up to a whole tick — the axis now fits tight to the data while keeping the readable ~100 bp / ~50 bp gridlines — and the `gcbias` chart drops a manual y2-axis-label offset that was compensating for label mis-positioning in the old release. The upgrade additionally tightens categorical bar spacing, text measurement, and legend-box sizing across every chart with no code change on our side.

### Fixed

- **`basic` base-distribution-by-cycle now reverse-complements reverse-strand reads.** Reverse-strand reads had their per-cycle index reversed but their bases left un-complemented, so each reverse read contributed the *complement* of the base that was sequenced. This flattened the per-cycle curves toward near-equal A≈T and C≈G (most visible on strand-asymmetric libraries such as EM-seq/bisulfite). The distribution is now computed in sequencing order, so an aligned BAM agrees with its Picard `RevertSam`-reverted counterpart. Mean-quality-by-cycle and quality-score-distribution were unaffected.
- **The mean-quality-by-cycle (`basic`) and transcript-coverage (`rna`) charts now anchor their y-axis at 0.** Both are filled-area charts whose data sits well above zero, so the y-axis previously auto-floored just below the data and the fill "floated", visually exaggerating the falloff (e.g. transcript coverage looked like it crashed to zero at the 5'/3' ends when it was really ~0.4). Anchoring at 0 makes the fill height read as true magnitude.

## [0.4.0] - 2026-07-06

### Added

- **Toolkit-wide `--threads` option** to spread work across multiple cores. It
  is a single global budget — accepted before or after the subcommand — that
  each command splits between decoding the input and its own work. The
  single-pass tools (`wgs`, `alignment`, …) do their own work on one thread and
  hand the rest to input decoding; `multi` also spreads the budget across
  parallel collector workers. Left unset, the single-pass tools stay
  single-threaded and `multi` uses a small default pool. Gains are sub-linear
  and taper off past a handful of threads, and the sweet spot depends on the
  machine and whether the input is BAM or CRAM.
  ([#44](https://github.com/fulcrumgenomics/riker/pull/44))
- **`rna` command** — RNA-seq QC metrics in a single pass. It ports Picard
  `CollectRnaSeqMetrics` (base composition, strand specificity, 5'/3' transcript
  coverage bias) and fgbio `EstimateRnaSeqInsertSize` (insert size in transcript
  space), and adds read-level genomic origin, gene detection, splice-junction
  annotation, transcript integrity (TIN), and a per-biotype read-count file. The
  gene model can be UCSC refFlat, GTF, or GFF3 (GENCODE / Ensembl / RefSeq),
  auto-detected from the file, with contig-name reconciliation (`chr` add/strip,
  `MT`↔`chrM`, RefSeq accessions). Strandedness auto-detects by default. Requires
  a coordinate-sorted input. Several Picard/fgbio behaviours are corrected along
  the way — see the "Differences in rna" section of the README and the new
  [ERRATA.md](./ERRATA.md) file..
  ([#40](https://github.com/fulcrumgenomics/riker/pull/40))

### Changed

- **`multi`'s `--threads` is now the total thread budget** — `multi` divides it
  between input decoding and collector workers rather than treating it as a bare
  worker count, and with no value given it now uses up to four threads (clamped
  to the core count) instead of two. Existing `riker multi --threads N`
  invocations keep working.
  ([#44](https://github.com/fulcrumgenomics/riker/pull/44))
- **`hybcap` now requires a coordinate-sorted BAM** and aborts on an
  out-of-order record; it was previously order-agnostic. (This is what lets the
  faster coordinate-order sweep below work.)
  ([#45](https://github.com/fulcrumgenomics/riker/pull/45))
- **All commands validate the output path at startup.** A misspelled or
  unwritable output directory now fails immediately, before any input is read,
  instead of after a full pass over the data.

### Performance

- **`wgs` is roughly 30% faster and uses about a third less peak memory.**
  Several changes to the depth pipeline combine to do it: a SIMD pileup kernel
  for the common no-overlapping-mate read, a slice-based depth buffer, and
  precomputing the reference's non-N regions once up front instead of reloading
  each contig. Measured on a 12× WGS BAM at four threads; output is
  byte-identical.
- **`hybcap` is ~1.5× faster** (about a quarter less CPU), from a
  coordinate-order coverage and bait/target sweep in place of the old per-base
  and interval-tree lookups. Output is byte-identical; note it now requires a
  sorted BAM (see Changed).
  ([#45](https://github.com/fulcrumgenomics/riker/pull/45))
- **`multi` threads more efficiently**, splitting it's thread budget between
  input decoding and worker threads with the split being optimized separately
  for BAM vs. CRAM input.
  ([#44](https://github.com/fulcrumgenomics/riker/pull/44),
  [#49](https://github.com/fulcrumgenomics/riker/pull/49))
- **Hash-heavy paths run faster** — frequency counters and the `error` model's
  covariate lookups moved from SipHash to fxhash.
  ([#48](https://github.com/fulcrumgenomics/riker/pull/48))

### Fixed

- **A failed single-pass run no longer leaves partial output on disk.** The
  single-threaded collector driver now skips finalization when the read loop
  errors (matching `multi`'s all-or-nothing behavior); previously a mid-stream
  error could still write a complete-looking but incomplete metrics file.
- **`multi` output is now deterministic across thread counts.** Each collector
  is pinned to a single worker and receives batches in file order, so a
  parallel run is byte-identical to a single-threaded one; previously the
  per-record output order of some collectors could, in rare circumstances,
  vary with the thread count.

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
