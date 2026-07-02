![Build](https://github.com/fulcrumgenomics/riker/actions/workflows/check.yml/badge.svg)
[![Version at crates.io](https://img.shields.io/crates/v/riker-ngs)](https://crates.io/crates/riker-ngs)
[![Documentation at docs.rs](https://img.shields.io/docsrs/riker-ngs)](https://docs.rs/riker-ngs)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/riker.svg?label=bioconda)](https://bioconda.github.io/recipes/riker/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/riker/blob/main/LICENSE)

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

**Warning: ALPHA SOFTWARE - USE AT YOUR OWN RISK**

This software is currently in **ALPHA**. While we have extensively tested these
tools across a wide variety of data, **no guarantees are made** regarding
correctness or stability.

## Overview

Riker is the spiritual successor to [picard](https://github.com/broadinstitute/picard) for next-generation sequencing QC.  It aims to provide modern, fast implementations of NGS QC metrics that tell you both "how does your sample look?" and often "why does it look bad?".  It is **not** intended to be a drop-in replacement for Picard.  Notably:

* Command line structure is different
* Output file formats are simplified and column names changed to resolve inconsistencies and confusion (and to use lower-case)
* Bugs that are difficult to fix in Picard due to implementation choices have been fixed in Riker, yielding slightly different outputs

Instead, Riker provides broadly similar tools to many of Picard's most widely used QC tools, with a goal of leading you to the same conclusions as before, in a lot less time.  Riker will remain focused on QC tools, and not attempt to replicate other tools from Picard.  Though other tools have been re-implemented for the better elsewhere too (e.g. [fgumi](https://github.com/fulcrumgenomics/fgumi)'s , `fastq`, `sort`, `dedup`, and `zipper` tools).

See the [Available Tools](#available-tools) section for a list of current tools.

## Motivation

The obvious question is: why not just fix up Picard?  Riker exists for a number of reasons:

* **Speed**: Most tools in riker are 4-6x faster than their counterparts in Picard; some are _much_ faster. With a fresh start it was much easier to make _all_ tools runnable via a single `multi` command, reducing complexity and saving more compute time (e.g. in Picard the `wgs` `hybcap`, and `error` tools all require separate invocations).
* **Cleaner Output**: Picard's tools output an inconvenient mostly-TSV format with a variable number of `#`-commented header lines and sometimes multiple independent tables in the same file, making them harder than they should be to parse programmatically _and_ annoying to review manually.  Riker's outputs are simple TSVs with one table per output file that can easily be routed into `cat file | column -t` or python's `csv.DictReader` with no fuss.  In addition Riker standardizes on having the `sample` be the first column in every file - so you can easily concatenate files for many samples.
* **Lightweight Distribution**: Picard carries a lot of baggage. It needs a JVM. It needs R and much of the tidyverse in order to produce plots.  The bioconda distribution requires a python interpreter to run it's wrapper script.  Running `pixi init && pixi add picard` results in a 1.2GB environment!  In contrast Riker is distributed as a single executable of < 10MB with no external dependencies.
* **Maintenance**: Maintenance of Picard has been [minimal for some time](https://github.com/broadinstitute/picard/commits/master/), with what little activity there is coming mostly from the community.  Picard is owned by, but no longer actively led by, the Broad Institute, making its path forward unclear.  Picard also suffers from some now-unnecessary complexity and years of less-than-necessary maintenance.  All of this makes a fresh start more appealing.

## Performance

Numbers below are from a [reproducible benchmark pipeline](benchmark-pipeline/) run on **2026-05-06** on a single AWS `r8id.xlarge` instance (4 vCPU, 32 GB RAM, local NVMe), against publicly-available 1000 Genomes 30× WGS BAMs (transcoded from CRAM) and the GIAB Ashkenazi exome trio.

Tool versions:

- **Riker**: 0.2.0 release candidate (built from source on the host)
- **Picard**: 3.4.0 on OpenJDK 25.0.2 (bioconda)
- **mosdepth**: 0.3.14 (bioconda)
- **samtools**: 1.23.1 (bioconda; used for staging — CRAM→BAM transcoding, downsampling with `view --subsample`, indexing)

### WGS — Riker wgs vs. Picard CollectWgsMetrics

| Sample    | BAM    | Cov | Riker Wall | Picard Wall | Speedup   | Riker RSS | Picard RSS |
|---        |---     |--- |---         |---          |---        |---        |---         |
| HG02675   | 7.1 GB | 4×  | 0:56       | 12:24       | **13.2×** | 0.74 GB   | 1.56 GB    |
| HG02675   | 21.7 GB| 15× | 2:56       | 32:20       | **11.0×** | 0.74 GB   | 5.97 GB    |
| HG02675   | 28.1 GB| 20× | 3:47       | 40:24       | **10.7×** | 0.74 GB   | 5.58 GB    |
| HG00188   | 37.5 GB| 30× | 5:13       | 1:02:30     | **12.0×** | 0.74 GB   | 5.24 GB    |
| HG02675   | 41.1 GB| 30× | 5:32       | 1:01:20     | **11.1×** | 0.74 GB   | 5.20 GB    |

### WGS — Riker multi vs. Picard Collect*Metrics

| Sample    | BAM    | Cov | Riker Wall | Picard Wall (CMM + CWM) | Speedup   | Riker RSS | Picard peak RSS |
|---        |---     |--- |---         |---                      |---        |---        |---              |
| HG02675   | 7.1 GB | 4×  | 1:25       | 21:16   (8:44 + 12:32)  | **15.1×** | 1.44 GB   | 1.60 GB         |
| HG02675   | 21.7 GB| 15× | 4:34       | 58:38   (26:32 + 32:06) | **12.8×** | 1.47 GB   | 5.97 GB         |
| HG02675   | 28.1 GB| 20× | 5:58       | 1:15:12 (34:21 + 40:51) | **12.6×** | 1.48 GB   | 5.58 GB         |
| HG00188   | 37.5 GB| 30× | 8:14       | 1:46:16 (46:36 + 59:40) | **12.9×** | 1.48 GB   | 5.24 GB         |
| HG02675   | 41.1 GB| 30× | 8:45       | 1:51:55 (51:13 + 60:42) | **12.8×** | 1.47 GB   | 5.20 GB         |

Riker was tested with a single invocation of `riker --threads 4 multi --tools wgs gcbias alignment basic isize`.
Picard was run twice, once for `CollectWgsMetrics` and once for `CollectMultipleMetrics` to generate a matching set of outputs.

"Picard peak RSS" is the larger of the two sequential JVM runs — typically dominated by `CollectWgsMetrics`, which scales with genome size + coverage.

### WGS — Riker wgs vs. mosdepth

| Sample        | BAM    | Cov | Riker Wall | mosdepth Wall | Δ                    | Riker RSS | mosdepth RSS |
|---            |---     |--- |---         |---            |---                   |---        |---           |
| HG02675_4x    | 7.1 GB | 4×  | 0:56       | 1:13          | Riker 23 % faster    | 0.74 GB   | 3.07 GB      |
| HG02675_15x   | 21.7 GB| 15× | 2:56       | 2:29          | mosdepth 15 % faster | 0.74 GB   | 2.42 GB      |
| HG02675_20x   | 28.1 GB| 20× | 3:47       | 2:59          | mosdepth 21 % faster | 0.74 GB   | 2.48 GB      |
| HG00188_30x   | 37.5 GB| 30× | 5:13       | 3:48          | mosdepth 27 % faster | 0.74 GB   | 2.52 GB      |
| HG02675_30x   | 41.1 GB| 30× | 5:32       | 4:05          | mosdepth 26 % faster | 0.74 GB   | 2.59 GB      |

Both tools running pure single-thread. `mosdepth` was run with its default `-t 0` for zero _extra_ decompression threads, and with `--no-per-base`.

Surprisingly Riker outperforms mosdepth on the low-coverage (4×) sample, while mosdepth wins at higher depths.  The 15-27% delta is largely explainable by the fact that Riker is performing per-base quality score filtering, and quality-score aware mate-overlap computations whereas mosdepth does not examine quality scores.

### Hybcap — Riker hybcap vs. Picard CollectHsMetrics

Hybcap measurements are the mean of three GIAB Ashkenazi trio samples (HG002, HG003, HG004), each a ~9.8 GB exome BAM aligned to hs37d5 with the Agilent SureSelect Human All Exon V5 capture kit.

| Sample          | BAM     | Kit         | Riker Wall | Picard Wall | Speedup   | Riker RSS | Picard RSS |
|---              |---      |---          |---         |---          |---        |---        |---         |
| AJ trio (mean)  | ~9.8 GB | Agilent v5  | 1:26       | 14:12       | **9.9×**  | 0.98 GB   | 2.45 GB    |

### Hybcap — Riker multi vs. Picard Collect*Metrics

| Sample          | BAM     | Kit         | Riker Wall | Picard Wall (CMM + CHsM) | Speedup   | Riker RSS | Picard peak RSS |
|---              |---      |---          |---         |---                       |---        |---        |---              |
| AJ trio (mean)  | ~9.8 GB | Agilent v5  | 1:45       | 22:07 (7:57 + 14:09)     | **12.6×** | 0.99 GB   | 3.23 GB         |

Riker was tested with a single invocation of `riker --threads 4 multi --tools hybcap alignment basic isize`.
Picard was run twice, once for `CollectHsMetrics` and once for `CollectMultipleMetrics` to generate a matching set of outputs.

"Picard peak RSS" is the larger of the two sequential JVM runs — dominated by `CollectHsMetrics` on the hybcap trio.

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
