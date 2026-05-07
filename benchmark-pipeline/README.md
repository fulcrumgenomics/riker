# riker benchmark-pipeline

Reproducible Snakemake pipeline that benchmarks **riker** against
**Picard** and **mosdepth** on real public WGS and exome data
(**Qualimap** is implemented as an optional comparator; see
[Comparator matrix](#profiles-benchmarked) below). Runs end-to-end on
a fresh AWS EC2 instance (or any clean Linux box) with public-internet
access only — no laptop-local files required.

## Profiles benchmarked

Four riker invocations, each with its corresponding fair comparators:

| Profile | Riker | Picard | mosdepth |
|---|---|---|---|
| `wgs-only` | `riker wgs` | `CollectWgsMetrics` | `mosdepth -x` |
| `wgs-bundle` | `riker multi --tools wgs alignment basic isize gcbias` | `CollectMultipleMetrics+CollectGcBiasMetrics` + `CollectWgsMetrics` | — |
| `hybcap-only` | `riker hybcap` | `CollectHsMetrics` | `mosdepth --by targets.bed` |
| `hybcap-bundle` | `riker multi --tools hybcap alignment basic isize` | `CollectMultipleMetrics` + `CollectHsMetrics` | — |

Threading: single-tool profiles run **1 thread** for every tool. Bundle
profiles run **riker `--threads 4`**; Picard `CollectMultipleMetrics`
runs with `samjdk.async_io_*=true` JVM properties to be at least
somewhat fair vs riker's multithreading. Supplementary Picard runs
(`CollectWgsMetrics`, `CollectHsMetrics`) stay at 1 thread.

**Qualimap** is implemented in `workflow/rules/run_qualimap.smk` but
not enabled in the default `COMPARATORS` matrix (Snakefile). It OOMs
on 30× WGS at both 16 GB and 24 GB heap on a 32 GB host even with the
retry loop; re-enable by adding `("qualimap", "main")` entries to
`COMPARATORS` if running on a larger box.

## Data (locked picks, no per-run randomness)

**WGS** — three NYGC ~30x CRAMs picked by file-size percentile across
all 2,504 samples in the public 1000 Genomes 30x bucket, plus a
downsample series anchored on the median-size sample:

BAM sizes are after CRAM→BAM transcoding; the empirical CRAM→BAM
ratio for 30× Illumina NovaSeq WGS is ~2.55× (lower coverage inflates
further — the 4× downsample is +29 % above linear extrapolation).

| Sample | Source | Approx BAM | Notes |
|---|---|---|---|
| `HG00188_30x` | `s3://1000genomes/.../ERR3240174/HG00188.final.cram` | ~37 GB | small (5th pct by CRAM size) |
| `HG02675_30x` | `s3://1000genomes/.../ERR3242389/HG02675.final.cram` | ~41 GB | medium (50th pct); downsample anchor |
| `HG02675_20x` | downsampled from `HG02675_30x` | ~28 GB | seed=1 |
| `HG02675_15x` | downsampled from `HG02675_30x` | ~22 GB | seed=2 |
| `HG02675_4x`  | downsampled from `HG02675_30x` | ~7 GB | seed=3, smoke fixture |
| `HG04131_30x` | `s3://1000genomes/.../ERR3243060/HG04131.final.cram` | ~52 GB | large (95th pct); **disabled by default** in `samples.wgs.tsv` (re-enable on hosts with ≥400 GB local NVMe) |

The depth series uses `samtools view --subsample <frac> --subsample-seed <N>`
with explicit seeds (1, 2, 3 — recorded in `config/samples.wgs.tsv`) so
the downsamples are deterministic and reproducible.

**Capture** — GIAB Ashkenazi trio with the Oslo University Hospital
Agilent SureSelect Human All Exon V5 prep (~9 GB BAMs):

| Sample | Source | Approx BAM |
|---|---|---|
| `HG002_av5` | GIAB FTP / Oslo / Agilent v5 | 9.6 GB |
| `HG003_av5` | GIAB FTP / Oslo / Agilent v5 | 8.4 GB |
| `HG004_av5` | GIAB FTP / Oslo / Agilent v5 | 9.5 GB |

Caveat: Agilent v5 is from 2015. It's the only GIAB-trio exome
dataset that fits the 5–15 GB BAM size target. The pipeline schema
accepts additional kits/samples without code changes — extend
`config/samples.hybcap.tsv` and `config/kits.yaml`.

`config/samples.hybcap.tsv` declares `reference: grch37_b37` for these
rows on the assumption the Oslo BAMs are b37-aligned. **Verify after
first staging:**

    samtools view -H stage/HG002_av5/input.bam | grep '^@SQ' | head -3

`@SQ` headers should start with `SN:1`, `SN:2`, ... (no `chr` prefix).
If the BAMs turn out to be hs37d5 instead, switch the `grch37_b37`
entry's `fasta_url` in `config/references.yaml` to point at hs37d5 and
re-run staging.

## Setup

### Pre-`./install.sh` checklist (fresh EC2 / Linux box)

These steps aren't done by `install.sh` itself:

1. **Mount the local NVMe instance store** (c5d / i4i / m6id / r8id-class
   instances ship with an unmounted ephemeral disk):

   ```bash
   sudo mkfs.xfs /dev/nvme1n1                 # check `lsblk` for the device
   sudo mkdir -p /mnt/scratch
   sudo mount /dev/nvme1n1 /mnt/scratch
   sudo chown $(whoami):$(whoami) /mnt/scratch
   ```

2. **Install `git` and clone the repo** (AL2023-minimal does not ship
   git):

   ```bash
   sudo dnf install -y -q git
   git clone https://github.com/fulcrumgenomics/riker.git
   cd riker/benchmark-pipeline
   ```

3. **`sudo` access**: `install.sh` shells out to `sudo dnf` / `sudo
   apt-get` to install the C/C++ toolchain if missing. The script will
   hang on a password prompt if your user lacks passwordless sudo.

4. **Clean AWS environment**: the staging rules use `aws s3 cp
   --no-sign-request` to read public 1KG bucket objects. If you have
   stale or expired AWS credentials in your shell environment
   (`AWS_PROFILE`, `AWS_ACCESS_KEY_ID`, etc.), some `aws-cli` versions
   ignore `--no-sign-request` and try to sign with the broken creds.
   `unset AWS_PROFILE AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY
   AWS_SESSION_TOKEN` before running the pipeline if you suspect this.

### `./install.sh`

```bash
./install.sh              # idempotent; fetches pixi, builds isolated rust, builds riker
```

Flags:
- `--system-rust`: use the cargo on PATH instead of installing rustup into `.rust/`.
- `--skip-build`: set up envs but don't (re)build riker.

### Pointing `stage_dir` at the NVMe

The default `stage_dir` in `config/performance.config.yaml` is the
relative path `stage/` (i.e. inside the repo). To stage on the NVMe
instead, override the value at runtime — Snakemake's `--config` flag is
forwarded after `--`:

```bash
./run.sh config/performance.config.yaml --cores 4 -- \
    --config stage_dir=/mnt/scratch/riker-bench/stage
```

Don't edit `performance.config.yaml` in place — that change tracks back
into your fork's git history.

## Smoke test (validates the pipeline on the smallest fixture)

```bash
./run.sh config/smoke.config.yaml --cores 4
```

Smoke runs the full DAG on **one sample** (HG02675 downsampled to ~4x,
~3.5 GB BAM) with a single replicate. Total wall time should be well
under 30 minutes on a modern EC2 instance with local NVMe.

After the smoke run:

- `results/bench.tsv` — one row per `(sample, profile, tool, rep)` cell
- `results/bench_summary.tsv` — replicates collapsed
- `results/plots/0[1-6]*.pdf` — six diagnostic plots
- `results/host.json` — host metadata (CPU, mem, EC2 instance type, etc.)

## Full performance run

```bash
./run.sh config/performance.config.yaml --cores $(nproc)
```

Storage budget: the staged BAMs total ~180 GB with the default
five-WGS sample set (HG04131 disabled), or ~230 GB with all six. Use a
c5d / i4i / m6id / r8id-class EC2 instance with local NVMe and
configure `stage_dir:` (in `config/performance.config.yaml`) to point
there. **Do not use EBS** for staged BAMs — random-read variance
contaminates measurements. The pipeline was first validated on
`r8id.xlarge` (4 vCPU / 32 GB / 220 GB NVMe) which fits the default
sample set with ~45 GB headroom.

Replicates default to 1; bundle profiles run riker at 4 threads.
With ~45 timed cells × ~2-50 minutes each, expect ~12-20 hours of
wall-clock on a 4-vCPU host. Bump `replicates:` in the config to 3
once you have a multi-day window and want min/max variance bars.

### Resuming after a failure

`./run.sh` already passes `--rerun-incomplete`. Re-running the same
command will resume — Snakemake skips rules whose declared `output:`
files already exist. Per-tool diagnostics (`cmdline.txt` and
`tool.log`) are now declared as `log:` rather than `output:` so they
survive Snakemake's auto-cleanup on rule failure: when a tool errors
out, look at
`results/run/<sample>/<profile>/<tool>/rep<N>/tool.log`.

If you rebuild riker or edit a rule's shell text and want Snakemake to
re-fire only the file-mtime-affected rules (rather than everything
whose params hash changed), pass `--rerun-triggers mtime` after `--`:

```bash
./run.sh config/performance.config.yaml --cores 4 -- \
    --config stage_dir=/mnt/scratch/riker-bench/stage \
    --rerun-triggers mtime
```

### Kit interval-list `@SQ` rewrite

Picard `CollectHsMetrics` strictly requires the bait/target
interval_list's `@SQ` headers to match the BAM's exactly (any size
mismatch is a fatal error). Real-world kit interval-lists are often
shipped with a slightly different `@SQ` set than the BAM's reference.
The `fetch_kit_intervals` rule (`workflow/rules/stage_inputs.smk`)
strips the upstream `@SQ` headers and concatenates the kit's data rows
under the reference's full `.dict`. Riker's `hybcap` is permissive
about this and would work with the original headers — only Picard
needs the rewrite.

## Output schema

`results/bench.tsv` — wide TSV, one row per `(sample, profile, tool, rep)`:

| Column | Meaning |
|---|---|
| `sample`, `input_type`, `profile`, `tool`, `tool_family` | identifiers |
| `tool_version` | first-line of the tool's `--version` output |
| `picard_bundle_role` | `main` / `base` / `supp` / `""` — for summing into "fair Picard bundle" |
| `threads` | requested thread count |
| `effective_threads` | actual (always 1 for picard) |
| `rep` | 1..N |
| `input_format` | always `bam` (CRAMs are transcoded during staging) |
| `sample_size_gb` | on-disk BAM size after staging |
| `wall_s`, `user_s`, `sys_s`, `cpu_percent` | from `/usr/bin/time -v` |
| `max_rss_kb`, `max_rss_gb` | peak resident set size |
| `throughput_mb_per_s` | sample_size_gb × 1024 / wall_s |
| `exit_status` | 0 = success |
| `host_*` | hostname, arch, cpu_model, cpu_count_logical, mem_bytes, os, os_release, aws_instance_type |

`results/bench_summary.tsv` — same schema with replicates collapsed
(median + min + max per numeric column; `rep_count` records how many
replicates contributed to each row).

`results/plots/` — what to look at first:
**`01_wall_time_leaderboard.pdf`** is the canonical "who's fastest"
bar chart and is the recommended starting point. **`04_riker_bundle_vs_picard.pdf`**
is the comparative claim of the project (riker single-pass vs sum of
Picard parts). The full list:

1. `01_wall_time_leaderboard.pdf` — bars per `(tool, sample)` faceted by profile (headline)
2. `02_throughput_vs_size.pdf` — log-scale x = sample_size_gb (carried by the depth series)
3. `03_max_rss_vs_size.pdf` — same axes, y = max_rss_gb (catches OOM-prone tools)
4. `04_riker_bundle_vs_picard.pdf` — `wall(riker bundle) / sum(wall(picard parts))` ratio per `(sample, profile)`
5. `05_riker_vs_picard_cmm.pdf` — riker multi vs Picard CollectMultipleMetrics head-to-head
6. `06_per_comparator_ratio.pdf` — wall ratio relative to riker (riker = 1.0)

## Resource lock semantics

Every timed rule reserves `bench=100` (the full pool, set by `run.sh`).
Every other rule (staging, aggregation, plotting) reserves `bench=1`.
Snakemake therefore never dispatches a non-timed rule **while** a timed
rule is running, but staging and DAG fill-in run in parallel during
benchmark idle time.

**Do not raise `bench` above 100.** It would silently allow concurrent
timed rules and contaminate every measurement.

## Repository conventions

- `config/` — TSV sample sheets and YAML configs. Edit these to extend
  the benchmark; no code changes required.
- `workflow/rules/*.smk` — one snakemake rule file per tool family. Each
  rule has its argv inline as a `shell:` block (no `tools/<tool>/render.py`
  indirection — for ~9 distinct tool invocations the inline form is
  simpler than the chelae-style render dispatch).
- `workflow/scripts/*.py`, `*.R` — aggregation, summarization, plotting.
  Ported from chelae's pipeline where possible (`parse_gnu_time.py`,
  `host_info.py`).
- `select_nygc.py` is **documentation-only**: it was used once to pick
  the three NYGC samples currently locked in `samples.wgs.tsv`. Re-run
  it to verify the picks are still small/medium/large, but the picks
  themselves are baked in for reproducibility.

## Out of scope (v1)

- Accuracy correlation between riker and Picard outputs. Not planned.
- Thread-scaling sweep (one thread setting per profile per spec).
- Multi-kit hybcap comparison.
- samtools stats.
- Multi-node / cluster execution.
- riker `error` subcommand (not in any of the four profiles).
- Cold-cache vs warm-cache comparison (warm only).
