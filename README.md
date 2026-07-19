![Build](https://github.com/fulcrumgenomics/riker/actions/workflows/check.yml/badge.svg)
[![Version at crates.io](https://img.shields.io/crates/v/riker-ngs)](https://crates.io/crates/riker-ngs)
[![Documentation at docs.rs](https://img.shields.io/docsrs/riker-ngs)](https://docs.rs/riker-ngs)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/riker.svg?label=bioconda)](https://bioconda.github.io/recipes/riker/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/riker/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21445767.svg)](https://doi.org/10.5281/zenodo.21445767)

# riker

Fast Rust CLI toolkit for sequencing QC metrics -- ports key QC metrics tools from [Picard](https://github.com/broadinstitute/picard) with cleaner output and better performance.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/riker/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/riker/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/riker/main/.github/logos/fulcrumgenomics-light.svg" height="125">
</picture>
</a><img src="https://raw.githubusercontent.com/fulcrumgenomics/riker/main/.github/logos/riker.png" height="125" style="padding-left: 30px;"/>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with riker and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Overview

Riker is the spiritual successor to [picard](https://github.com/broadinstitute/picard) for next-generation sequencing QC.  It aims to provide modern, fast implementations of NGS QC metrics that tell you both "how does your sample look?" and often "why does it look bad?".  It is **not** intended to be a drop-in replacement for Picard.  Notably:

* Command line structure is different
* Output file formats are simplified and column names changed to resolve inconsistencies and confusion (and to use lower-case)
* Bugs that are difficult to fix in Picard due to implementation choices have been fixed in Riker, yielding slightly different outputs

Instead, Riker provides broadly similar tools to many of Picard's most widely used QC tools, with a goal of leading you to the same conclusions as before, in a lot less time.  Riker will remain focused on QC tools, and not attempt to replicate other tools from Picard.  Though other tools have been re-implemented for the better elsewhere too (e.g. [fgumi](https://github.com/fulcrumgenomics/fgumi)'s , `fastq`, `sort`, `dedup`, and `zipper` tools).

See the [Available Tools](#available-tools) section for a list of current tools.

## Motivation

The obvious question is: why not just fix up Picard?  Riker exists for a number of reasons:

* **Speed**: Riker's tools run 12–38× faster than their Picard counterparts, measured across WGS, exome, and RNA-seq metrics (see [Performance](#performance)). Starting fresh also made it easy to compute _all_ the metrics in a single pass over the input via one `multi` command — cutting both complexity and compute, since in Picard `CollectWgsMetrics`, `CollectHsMetrics`, `CollectRnaSeqMetrics`, and `CollectMultipleMetrics` each require a separate invocation and a separate read of the data.
* **Cleaner Output**: Picard's tools output an inconvenient mostly-TSV format with a variable number of `#`-commented header lines and sometimes multiple independent tables in the same file, making them harder than they should be to parse programmatically _and_ annoying to review manually.  Riker's outputs are simple TSVs with one table per output file that can easily be routed into `cat file | column -t` or python's `csv.DictReader` with no fuss.  In addition Riker standardizes on having the `sample` be the first column in every file - so you can easily concatenate files for many samples.
* **Lightweight Distribution**: Picard carries a lot of baggage. It needs a JVM. It needs R and much of the tidyverse in order to produce plots.  The bioconda distribution requires a python interpreter to run it's wrapper script.  Running `pixi init && pixi add picard` results in a 1.2GB environment!  In contrast Riker is distributed as a single executable of < 10MB with no external dependencies.
* **Maintenance**: Maintenance of Picard has been [minimal for some time](https://github.com/broadinstitute/picard/commits/master/), with what little activity there is coming mostly from the community.  Picard is owned by, but no longer actively led by, the Broad Institute, making its path forward unclear.  Picard also suffers from some now-unnecessary complexity and years of less-than-necessary maintenance.  All of this makes a fresh start more appealing.

## Performance

Numbers below are from a [reproducible benchmark pipeline](benchmark-pipeline/) run on **July 5th 2026** on a single AWS `m8a.2xlarge` instance (8 vCPU AMD EPYC 9R45, 32 GB RAM, gp3 EBS provisioned at 1000 MiB/s), against publicly-available 1000 Genomes WGS BAMs (transcoded from CRAM and downsampled to 4–30×), the GIAB Ashkenazi exome trio, and two ENCODE RNA-seq BAMs. The benchmark was performed using **riker 0.4.0** vs. the latest version of each tool available from `bioconda` as of the benchmark date.

### Riker vs. Picard

Picard is the tool every one of Riker's tools was written to match, so we'll keep this part short: tool for tool, Riker computes the same metrics **12–38× faster**.

![Riker vs. Picard: single-tool speedup for WGS, exome, and RNA-seq metrics, from 12.5× to 38.3×](https://raw.githubusercontent.com/fulcrumgenomics/riker/main/docs/images/fig1_vs_picard.png)

Whether single-threaded tool vs. tool, or multi-threaded via `multi`, Riker is massively faster than the aging Picard.

### WGS depth: Riker vs. mosdepth

[mosdepth](https://github.com/brentp/mosdepth) is well known for being the fastest WGS coverage tool... until now. With both tools configured as closely to each other as possible, Riker is faster at every coverage *and* every thread count we measured — by **1.2–1.5×**, with the biggest multipliers at low coverage.

![Riker vs. mosdepth wall time across 4×–30× coverage, single-threaded; Riker 1.2–1.4× faster at every depth](https://raw.githubusercontent.com/fulcrumgenomics/riker/main/docs/images/fig2a_mosdepth_coverage.png)

![Riker vs. mosdepth wall time on a 30× genome across 1–4 threads; Riker 1.2–1.3× faster at every thread count](https://raw.githubusercontent.com/fulcrumgenomics/riker/main/docs/images/fig2b_mosdepth_threads.png)

Riker wins while doing *more* per base: it applies base quality score filtering and quality-aware mate-overlap correction (matching Picard's WGS semantics), while mosdepth is not base quality aware. Base quality filtering is the one difference we couldn't harmonize away — mosdepth has no switch for it — so Riker is doing strictly more work and still finishes first.

### RNA-seq: Riker vs. Picard & RustQC

For RNA-seq the field also includes Picard `CollectRnaSeqMetrics` and **RustQC** (Seqera's Rust RNA-seq QC toolkit, which is advertised as ~60× faster than the tools it reimplements), scope-matched so all three compute a comparable set of metrics. Single-threaded, Riker runs **31–38× faster than Picard** and **18–21× faster than RustQC**. RustQC threads well, but even at four threads it is still **~10× behind** Riker.

![Riker vs. Picard and RustQC on single- and paired-end RNA-seq; Riker 10–38× faster, with Picard shown as single-threaded only](https://raw.githubusercontent.com/fulcrumgenomics/riker/main/docs/images/fig3_rna.png)

## Installation

### Install from bioconda

Using pixi, after adding the `bioconda` channel to your configuration:

```
pixi add riker
```

Or using your favorite conda client:

```
conda install -c bioconda riker
```

### Installing with Cargo

Riker can be installed from [crates.io](https://crates.io/).  If you're unfamiliar with cargo (the Rust build tool) but want to go this route, start by installing [rustup](https://rustup.rs/).

```
cargo install riker-ngs
```

### Building from source

Similarly, to build from source you'll also need cargo.  Once installed you can:

```shell
# clone the repo
git clone https://github.com/fulcrumgenomics/riker

# build the release version:
cd riker
cargo build --release
```

## Available Tools

| Command | Description | Equivalent Tool(s) |
|---------|-------------|------------------|
| `alignment` | Collect alignment summary metrics | `picard CollectAlignmentSummaryMetrics` |
| `basic` | Collect base distribution, mean quality, and quality score distribution | `picard CollectBaseDistributionByCycle`, `MeanQualityByCycle`, `QualityScoreDistribution` |
| `error` | Collect base-level error metrics (mismatch, overlap, indel) | `picard CollectSamErrorMetrics` |
| `isize` | Collect insert size distribution metrics | `picard CollectInsertSizeMetrics` |
| `wgs` | Collect whole-genome coverage metrics | `picard CollectWgsMetrics` |
| `gcbias` | Collect GC bias metrics | `picard CollectGcBiasMetrics` |
| `hybcap` | Collect hybrid capture (bait/target) metrics | `picard CollectHsMetrics` |
| `rna` | Collect RNA-seq metrics (read origin, base composition, gene detection, splice junctions, strand, coverage bias, transcript integrity, transcript-space insert size) | `picard CollectRnaSeqMetrics`, `fgbio EstimateRnaSeqInsertSize`, RSeQC, RNA-SeQC, Qualimap |
| `multi` | Run multiple collectors in a single BAM pass | `picard CollectMultipleMetrics` |
| `docs` | Print metric field documentation | -- |

## Usage

For detailed usage of each command, run:
```bash
riker <command> --help
```

### Examples

Collect basic QC metrics (base distribution, mean quality, quality score distribution):

```bash
riker basic -i sample.bam -o out_prefix
```

Collect alignment summary metrics with a reference:

```bash
riker alignment -i sample.bam -r ref.fa -o out_prefix
```

Collect base-level error metrics with stratification:

```bash
riker error -i sample.bam -r ref.fa -o out_prefix
riker error -i sample.bam -r ref.fa -o out_prefix --vcf known.vcf.gz --stratify-by read_num,cycle bq
```

Collect insert size metrics:

```bash
riker isize -i sample.bam -o out_prefix
```

Collect whole-genome coverage metrics:

```bash
riker wgs -i sample.bam -r ref.fa -o out_prefix
riker wgs -i sample.bam -r ref.fa -o out_prefix -L intervals.bed
```

Collect GC bias metrics:

```bash
riker gcbias -i sample.bam -r ref.fa -o out_prefix
riker gcbias -i sample.bam -r ref.fa -o out_prefix --exclude-intervals artifacts.bed
```

Collect hybrid capture metrics:

```bash
riker hybcap -i sample.bam -o out_prefix --baits baits.bed --targets targets.bed
riker hybcap -i sample.bam -o out_prefix --baits baits.bed --targets targets.bed -r ref.fa
```

The input BAM must be coordinate-sorted.

Collect RNA-seq metrics (the gene-model format — UCSC refFlat, GTF, or GFF3 — is auto-detected):

```bash
riker rna -i sample.bam -o out_prefix --gene-model gencode.gtf.gz
riker rna -i sample.bam -o out_prefix --gene-model refFlat.txt.gz --ribosomal-intervals rRNA.interval_list
riker rna -i sample.bam -o out_prefix --gene-model genes.gff3.gz --strand reverse
```

Run multiple collectors in a single pass:

```bash
riker multi \
  -i sample.bam \
  -r ref.fa \
  -o out_prefix \
  --tools alignment isize basic hybcap \
  --hybcap::baits baits.bed \
  --hybcap::targets targets.bed
```

### Threads

Riker takes a single, whole-toolkit `--threads` budget. It's a global option, so
it works either before or after the subcommand:

```bash
riker --threads 4 wgs -i sample.bam -r ref.fa -o out_prefix
riker multi --threads 4 -i sample.bam -r ref.fa -o out_prefix --tools wgs gcbias basic alignment isize
```

`--threads` is the number of cores riker will try to saturate; it may spin a
few more threads internally (e.g. a dispatch thread) to keep them busy. Each
command spends the budget as it sees fit. The single-pass tools (`wgs`,
`alignment`, …) run their own work on one thread and hand the remainder to
input decoding; `multi` lays it out across input-decode threads, a dispatch
thread, and parallel collector workers, tuned to the input format. Left unset,
the single-pass tools run single-threaded and `multi` uses a small default pool.

The payoff is format-dependent: **BAM** inflate is cheap, so `multi` is
compute-bound and leans on collector workers; **CRAM** decode is expensive and
parallelizes well, so it leans on decode threads (and CRAM peak memory rises
with thread count, especially for the heavier `archive`/`small` codecs, while
BAM stays essentially flat). Gains are sub-linear — as a rough guide returns
flatten past ~6 threads for BAM and ~8 for CRAM — but the sweet spot is
platform- and input-dependent, and over-supplying threads can slow a run down,
so it's worth measuring on your own data.

## Output Format

Riker produces PDF plots using [kuva](https://github.com/Psy-Fer/kuva), and plain TSV output designed for easy downstream consumption:

- **Lowercase snake_case headers** (e.g., `total_reads`, `mean_insert_size`)
- **Tab-separated** with no metadata or comment lines
- **Fractions use `frac_` prefix** instead of `pct_` (e.g., `frac_aligned` not `pct_aligned`)
- **No per read-group or library breakdown** -- all reads in the file are combined

To see detailed documentation on the columns output in each file type run `riker docs`.

## Differences vs. Picard and Other Tools

Riker aims to reproduce the metrics of Picard — and the other tools it ports — as closely as
possible. Where behavior intentionally differs (corrections of known bugs, or deliberate design
choices), every difference is documented per command in **[ERRATA.md](ERRATA.md)**:

- **[`alignment`](ERRATA.md#alignment)** — `mean_aligned_read_length` is over aligned (not all PF) reads; improper-pair counting requires a mapped mate.
- **[`hybcap`](ERRATA.md#hybcap)** — chimeric-pair overlap clipping is fixed; per-target GC is N-aware.
- **[`gcbias`](ERRATA.md#gcbias)** — explicit read-accounting schema, a forward-strand window-binning fix, and optional `--min-mapq` / `--exclude-intervals`.
- **[`error`](ERRATA.md#error)** — reference-N positions skipped, mismatch denominator excludes insertions, raw-rate Q-scores, and 17 of Picard's 32 stratifiers ported.
- **[`rna`](ERRATA.md#rna)** — a single-pass superset of Picard `CollectRnaSeqMetrics` + fgbio `EstimateRnaSeqInsertSize`, plus RSeQC / RNA-SeQC / Qualimap metrics (read origin, gene detection, splice junctions, transcript integrity (TIN), biotype); several Picard coverage/bias off-by-ones are corrected; coordinate-sorted input is required. Transcript-space insert size uses the `MC` (mate CIGAR) tag when present and otherwise estimates FR pairs from a spec-compliant `TLEN` (add `MC` with `samtools fixmate -m` for exact results). The biggest absolute difference from the reference tools is **TIN**, which riker reports ~14–16 points higher than RSeQC/RustQC (same shape, different scale — riker scores one well-covered transcript per gene); both track degradation equally well and riker is notably more depth-robust — see [docs/rna-integrity.md](docs/rna-integrity.md).

## Authors

- [Tim Fennell](https://github.com/tfenne)

## Sponsors

Development of riker is supported by [Fulcrum Genomics](https://www.fulcrumgenomics.com).

[Become a sponsor](https://github.com/sponsors/fulcrumgenomics)

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fulcrumgenomics/riker/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/riker/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/riker/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.
