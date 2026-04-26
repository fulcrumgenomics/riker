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
  <source media="(prefers-color-scheme: dark)" srcset=".github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset=".github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src=".github/logos/fulcrumgenomics-light.svg" height="125">
</picture>
</a><img src=".github/logos/riker.png" height="125" style="padding-left: 30px;"/>
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
```

Collect hybrid capture metrics:

```bash
riker hybcap -i sample.bam -o out_prefix --baits baits.bed --targets targets.bed
riker hybcap -i sample.bam -o out_prefix --baits baits.bed --targets targets.bed -r ref.fa
```

Run multiple collectors in a single pass:

```bash
riker multi \
  -i sample.bam \
  -r ref.fa \
  -o out_prefix \
  --tools alignment isize basic hypcap \
  --hybcap::baits baits.bed \
  --hybcap::targets targets.bed
```

## Output Format

Riker produces PDF plots using [kuva](https://github.com/Psy-Fer/kuva), and plain TSV output designed for easy downstream consumption:

- **Lowercase snake_case headers** (e.g., `total_reads`, `mean_insert_size`)
- **Tab-separated** with no metadata or comment lines
- **Fractions use `frac_` prefix** instead of `pct_` (e.g., `frac_aligned` not `pct_aligned`)
- **No per read-group or library breakdown** -- all reads in the file are combined

To see detailed documentation on the columns output in each file type run `riker docs`.

## Key Differences vs. Picard

Riker aims to produce metrics identical to Picard's equivalent tools. This section documents where there are known and expected functional differences in Riker vs. Picard.

### Differences in alignment vs. CollectAlignmentSummaryMetrics

#### `mean_aligned_read_length`

**Picard** computes mean aligned read length over all PF reads, including unmapped reads which contribute zero to the sum. The denominator is `PF_READS`, so unmapped reads dilute the average toward zero.

**riker** uses `aligned_reads` as the denominator, giving the mean length of reads that actually aligned.

**Impact:** riker produces slightly higher values than Picard — typically ~0-2bp for WGS data with high alignment rates. The difference grows with the fraction of unmapped reads.

#### `reads_improperly_paired` / `frac_reads_improperly_paired`

**Picard** counts all mapped, paired, non-proper reads as improperly paired, including reads whose mate is unmapped (no `is_mate_unmapped()` guard).

**riker** requires the mate to also be mapped before counting a read toward `aligned_reads_in_pairs` or `reads_improperly_paired`. This avoids inflating the improper-pair count with reads whose mate simply failed to align.

**Impact:** riker reports fewer improperly paired reads than Picard. On typical WGS data the difference is usually small, in the single digit percent range. `aligned_reads_in_pairs` and `frac_reads_improperly_paired` are also affected.

### Differences in hybcap vs. CollectHsMetrics.

#### `frac_exc_overlap` (Picard: `PCT_EXC_OVERLAP`)

riker reports a slightly lower `frac_exc_overlap` than Picard — typically below 1% relative difference (e.g. 0.032542 vs 0.032646 on a 1000G exome dataset).

**Cause:** Picard's overlap-clipping function (`SAMUtils.getNumOverlappingAlignedBasesToClip()` in htsjdk) does not verify that the read and its mate are mapped to the same contig before computing the overlap. It compares `getMateAlignmentStart()` against `getAlignmentStart()` as raw integers regardless of reference sequence. For chimeric read pairs — where one end maps to a different chromosome — Picard may compute a spurious overlap when the mate's coordinate on a different contig happens to be within the range of the reads aligned coordinates.  See https://github.com/broadinstitute/picard/issues/2039.

riker's equivalent function checks `reference_sequence_id != mate_reference_sequence_id` and returns 0 for chimeric pairs, correctly skipping overlap clipping when the reads are on different contigs.

**Impact:** `frac_exc_overlap` and `frac_exc_off_target` are affected; the latter because bases now correctly *not* categorized as overlapping are usually excluded subsequently for being off-target. The magnitude depends on the fraction of chimeric read pairs in the data, but is usually very small. 

#### Per-target GC content and GC dropout

riker and Picard compute per-target GC fraction differently when a target's reference sequence contains ambiguous (N) bases.

**Picard** counts N bases in the denominator but not the numerator: `gc_fraction = (G + C) / (A + T + G + C + N)`. This dilutes the GC fraction toward zero in proportion to the number of N bases in the target.

**riker** treats N bases as maximally uncertain, contributing 0.5 to the GC numerator and 1.0 to the denominator: `gc_fraction = (G + C + 0.5 * N) / (A + T + G + C + N)`. A target that is entirely N bases reports a GC fraction of 0.5 (maximum uncertainty) rather than 0.0.

Neither approach is fully correct — the true GC content of N bases is unknown — but riker's treatment avoids the bias toward low GC that Picard introduces. In practice the difference is only visible for targets where N bases make up a meaningful fraction of the reference sequence.

**Impact:** Per-target GC values can differ for N-containing targets, which in turn affects the GC bias curve and the `at_dropout` and `gc_dropout` summary metrics. All other per-target and summary metrics are unaffected.

### Differences in error vs. CollectSamErrorMetrics

#### Reference N bases

**Picard** counts bases at reference-N positions toward `TOTAL_BASES` (and may count them as errors since the read base does not match N).

**riker** skips positions where the reference base is N, since errors cannot be meaningfully assessed without a known reference base.

**Impact:** riker reports slightly fewer `total_bases` and `error_bases` than Picard. The difference is proportional to the number of reference-N positions covered by reads — typically negligible for well-assembled references.

#### Mismatch `total_bases` includes insertion bases in Picard

**Picard's** mismatch metric (`ERROR`) includes insertion bases in `TOTAL_BASES`. This happens because `BaseErrorCalculator.addBase()` increments `nBases` for both aligned (Match) bases and insertion bases, and `SimpleErrorCalculator` inherits this count as the denominator.

**riker** counts only aligned bases in the mismatch `total_bases`. Insertion bases are counted separately in the indel metric.

**Impact:** Picard's mismatch `TOTAL_BASES` is higher than riker's by the number of inserted bases passing filters. The mismatch `error_bases` (numerator) is identical between the two tools — only the denominator differs. riker's separation is arguably cleaner since insertion bases are not relevant to the substitution error rate.

#### Insert size stratification

**Picard** caps insert size at `readLength * 10` for the `INSERT_LENGTH` stratifier. Reads with larger absolute insert sizes (e.g., chimeric pairs) are binned at the cap value.

**riker** excludes reads with absolute insert size above `--max-isize` (default 1000) from the `isize` stratifier entirely. These reads are still counted by all other stratifiers and in the `all` group. The threshold can be adjusted via `--max-isize`.

**Impact:** Picard may have a large bin at the cap value containing chimeric reads, while riker omits them. Overall error rates are unaffected.

#### Insertions at read start

**Picard** counts insertions that occur before any aligned base in the read (e.g., CIGARs starting with `nI` or `nSnI`). The locus iterator attaches these to the preceding reference position.

**riker** skips insertions that have no preceding aligned base (no anchor position), since there is no reference context for stratification.

**Impact:** Picard reports slightly more insertions and inserted bases than riker. The difference is small — typically a few hundred events out of tens of thousands. Deletion counts are unaffected and match exactly between the tools.

#### Q-score computation

**Picard** computes Q-scores using a Bayesian prior: `Q = -10 * log10((errors + 0.001) / (total + 1))`, rounded to the nearest integer. The prior (configurable via `PRIOR_Q`, default 30) prevents infinite Q-scores when there are zero errors.

**riker** computes Q-scores from the raw error rate: `Q = -10 * log10(errors / total)`, reported to two decimal places. When there are zero errors, a cap of Q99 is used.

**Impact:** Q-scores differ slightly due to the prior and rounding. The underlying counts (numerator and denominator) are comparable; only the derived Q-score differs.

#### Stratifiers not ported from Picard

Picard's `CollectSamErrorMetrics` defines 32 stratifiers. riker ports 15 of them (`all`, `bq`, `mapq`, `cycle`, `read_num`, `strand`, `pair_orientation`, `isize`, `gc`, `read_base`, `ref_base`, `hp_len`, `pre_dinuc`, `post_dinuc`, `context_3bp`). The following Picard stratifiers are **not** available in riker:

- `PAIR_PROPERNESS` — whether the read is in a proper pair
- `HOMOPOLYMER` — the homopolymer base (A/C/G/T) at the current position
- `BINNED_HOMOPOLYMER` — homopolymer length bucketed into bins
- `BINNED_CYCLE` — machine cycle bucketed into bins
- `SOFT_CLIPS` — number of soft-clipped bases in the read
- `MISMATCHES_IN_READ` — total mismatches in the read (`nm` is similar)
- `TWO_BASE_PADDED_CONTEXT` — 5-base context (2bp each side) (`Context3bp` is provided)
- `CONSENSUS` — whether the read is a consensus/duplex read
- `NS_IN_READ` — number of N bases in the read
- `INSERTIONS_IN_READ` — number of insertion events in the read
- `DELETIONS_IN_READ` — number of deletion events in the read
- `INDELS_IN_READ` — number of indel events in the read
- `FLOWCELL_TILE` — flowcell tile from the read name
- `FLOWCELL_X` — flowcell X coordinate from the read name
- `FLOWCELL_Y` — flowcell Y coordinate from the read name
- `READ_GROUP` — the read group of the read

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
