# Contributing to Riker

Riker is a fast Rust CLI toolkit that ports key QC metrics tools from
[Picard](https://broadinstitute.github.io/picard/) and
[fgbio](https://fulcrumgenomics.github.io/fgbio/). When porting from these
tools, the goal is to match correctness while improving the interface and output
format.

## Getting Started

**Prerequisites:** Rust stable (see `Cargo.toml` for minimum version) and
[cargo-nextest](https://nexte.st/).

```bash
cargo build              # debug build
cargo build --release    # optimized build (portable, runs anywhere)
```

`cargo build --release` deliberately does **not** set `target-cpu=native`
so that `cargo install riker-ngs` from crates.io produces a binary that
runs on any reasonable hardware. Distribution channels handle per-CPU
tuning separately (`cargo multivers --profile dist` for x86_64 release
artifacts and bioconda; `cargo build --profile dist` for aarch64). For a
locally-tuned benchmarking build, opt in:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Verification Checklist

Run all four before submitting changes:

```bash
cargo ci-fmt             # check formatting
cargo ci-lint            # clippy, pedantic, warnings-as-errors
cargo clippy --all-targets
cargo ci-test            # all tests via nextest
```

The `ci-*` aliases are defined in `.cargo/config.toml`.

## Crate Structure

| Crate | Purpose |
|-------|---------|
| `riker` | Binary and library (published as `riker_lib`) |
| `riker_derive` | Proc-macro crate: `#[derive(MetricDocs)]` and `#[multi_options]` |

## Architecture

### Commands

Each subcommand is a struct that implements the `Command` trait
(`src/commands/command.rs`) and is dispatched via the `Subcommand` enum in
`src/main.rs`.

### Collectors

Each command's core logic lives in a `Collector` struct implementing the
`Collector` trait (`src/collector.rs`). Collectors store all configuration
(output paths, reference, thresholds) as fields set at construction. The trait
methods only receive the BAM header and records — this design enables the `multi`
command to feed a single BAM pass to multiple collectors in parallel.

### Metrics

Per-command metric structs live alongside their collector (e.g.
`InsertSizeMetric` in `src/commands/isize.rs`). Metric structs derive
`Serialize`, `Deserialize`, and `MetricDocs`. Shared serialization utilities
live in `src/metrics.rs`.

### Shared CLI Options

Reusable option groups in `src/commands/common.rs` (`InputBamOptions`,
`ReferenceOptions`, `OptionalReferenceOptions`, `IntervalOptions`,
`DuplicateOptions`) are composed into commands via `#[command(flatten)]`.

## Adding a New Command

This section walks through every step needed to add a new metric collector to
riker, using `isize` and `hybcap` as reference examples.

### 1. Create the command module

Create `src/commands/<name>.rs`. The file should contain four things, in order:

**File suffixes** — constants for each output file:

```rust
pub const METRICS_SUFFIX: &str = ".<name>-metrics.txt";
```

**Command struct** — the CLI entry point. Flatten shared option groups and your
options struct:

```rust
#[derive(Args, Debug, Clone)]
#[command(long_about, after_long_help = "Examples:\n  riker <name> -i in.bam -o out")]
pub struct MyCommand {
    #[command(flatten)]
    pub input: InputBamOptions,

    #[arg(short = 'o', long, value_name = "PREFIX")]
    pub output: PathBuf,

    #[command(flatten)]
    pub reference: ReferenceOptions,  // if needed

    #[command(flatten)]
    pub options: MyOptions,
}
```

**Options struct** — tool-specific tuning parameters. This struct serves double
duty: it configures the standalone command *and* (via the `#[multi_options]`
macro) generates the prefixed options for the `multi` command. See
[The `#[multi_options]` macro](#the-multi_options-macro) below for details.

```rust
#[riker_derive::multi_options("<name>", "<Help Heading>")]
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Tuning")]
pub struct MyOptions {
    /// Include duplicate reads in metric calculations.
    #[arg(long, default_value_t = false)]
    pub include_duplicates: bool,
    // ...
}
```

**Metric struct(s)** — one struct per output file, deriving `Serialize`,
`Deserialize`, and `MetricDocs`:

```rust
/// Summary metrics for the <name> collector.
#[derive(Debug, Serialize, Deserialize, MetricDocs)]
pub struct MyMetric {
    /// Description of this field (becomes the `riker docs` entry).
    pub field_name: u64,

    /// A fractional field with 5 decimal places.
    #[serde(serialize_with = "serialize_f64_5dp")]
    pub frac_something: f64,
}
```

**Collector struct** — implements the `Collector` trait:

```rust
pub struct MyCollector { /* fields */ }

impl MyCollector {
    #[must_use]
    pub fn new(prefix: &Path, options: &MyOptions) -> Self { /* ... */ }
}

impl Collector for MyCollector {
    fn initialize(&mut self, header: &Header) -> Result<()> { /* ... */ }
    fn accept(&mut self, record: &RecordBuf, header: &Header) -> Result<()> { /* ... */ }
    fn finish(&mut self) -> Result<()> { /* ... */ }
    fn name(&self) -> &'static str { "<name>" }
}
```

**`Command` impl** — the standalone execution path. Use `for_each_record` to
iterate records with buffer reuse (avoids per-record allocation for BAM/SAM):

```rust
impl Command for MyCommand {
    fn execute(&self) -> Result<()> {
        let (mut reader, header) = AlignmentReader::new(&self.input.input, ...)?;
        let mut collector = MyCollector::new(&self.output, &self.options);
        collector.initialize(&header)?;
        let mut progress = ProgressLogger::new("<name>", "reads", 5_000_000);
        reader.for_each_record(&header, |record| {
            progress.record_with(record, &header);
            collector.accept(record, &header)
        })?;
        progress.finish();
        collector.finish()
    }
}
```

### 2. Register the command

- Add `pub mod <name>;` to `src/commands/mod.rs`
- Add a variant to the `Subcommand` enum in `src/main.rs`:
  ```rust
  MyCmd(MyCommand),
  ```
- Add a match arm in `impl Command for Subcommand::execute()`:
  ```rust
  Subcommand::MyCmd(c) => c.execute(),
  ```

### 3. Integrate with the multi command

The `multi` command (`src/commands/multi.rs`) runs multiple collectors in a
single BAM pass. Every new collector should be wired in:

1. **Add a `CollectorKind` variant** and update its `Display` impl
2. **Import** `MyCollector` and `MultiMyOptions` (generated by the macro)
3. **Add a `#[command(flatten)]`** field to the `Multi` struct:
   ```rust
   #[command(flatten)]
   pub my_opts: MultiMyOptions,
   ```
4. **Add a match arm** in `build_collectors()`:
   ```rust
   CollectorKind::MyColl => {
       let opts = self.my_opts.clone().validate()?;
       collectors.push(Box::new(MyCollector::new(&self.output, &opts)));
   }
   ```
5. **Decide the default** — add the new kind to the `default_values_t` list on
   `Multi::tools` if the tool should run by default, or leave it opt-in (like
   hybcap) if it requires extra inputs.

### 4. Register metric docs

Add your metric structs to `src/commands/docs.rs` — add
`render_metric_docs_text::<MyMetric>` and `render_metric_docs_markdown::<MyMetric>`
calls in both the `"text"` and `"markdown"` arms.

### 5. Add tests

- **Integration tests** in `tests/test_<name>.rs` — build BAM data
  programmatically using `SamBuilder`, run the command, deserialize output with
  `read_metrics_tsv`, and assert on results.
- **Multi integration tests** in `tests/test_multi.rs` — add at least one test
  that runs your collector through the `multi` command.
- **Unit tests** in an inline `#[cfg(test)]` module for any non-trivial internal
  logic.

### 6. Verify

Run the full verification checklist.

### The `#[multi_options]` macro

The `#[multi_options("prefix", "Heading")]` attribute macro
(in `riker_derive`) generates a `Multi<StructName>` companion struct for use
in the `multi` command. It prefixes every CLI flag with `prefix::` (e.g.
`--isize::include-duplicates`) and groups them under the given help heading.

The macro classifies each field into one of three categories:

| Category | Original type | Has `default_value_t`? | Multi struct type | In `validate()` |
|----------|--------------|------------------------|-------------------|----------------|
| **Optional** | `Option<T>` | n/a | `Option<T>` | Passed through |
| **Defaulted** | `T` | Yes | `T` (with default) | Passed through |
| **Required** | `T` | No | `Option<T>` | `.ok_or_else(...)` with error message |

**What it generates:**

- `Multi<Name>` struct with `clap::Args`, all fields prefixed.
- `validate(self) -> anyhow::Result<Name>` — converts `Multi<Name>` back to the
  original struct, checking that required fields are `Some`. Required fields get
  an error message like `"--hybcap::baits is required when hybcap is selected"`.
- `From<Name> for Multi<Name>` — wraps required fields in `Some()`, passes
  others through. Useful for constructing multi options from standalone options
  in tests.

**Rules for Options structs:**

- The struct must have `#[derive(Args, Debug, Clone)]` and implement `Default`.
- All fields must be flat — no `#[command(flatten)]`.
- Every non-`Option` field that isn't inherently required at the CLI level
  should have `#[arg(long, default_value_t = ...)]`.
- Fields that *are* required (like file paths) should be plain types without
  `default_value_t`. The macro wraps them in `Option` for the multi struct and
  generates validation. The `Default` impl can use a placeholder
  (e.g. `PathBuf::new()`) since standalone clap validation ensures it's always
  provided, and `validate()` checks it in the multi path.
- Doc comments on fields become the `--help` text; for required fields the macro
  appends `" Required when <prefix> is selected."`

## Output Format Conventions

- **TSV output**: lowercase `snake_case` headers, tab-separated, no metadata or
  comment lines
- **Fractions**: use `frac_` prefix (not `pct_`)
- **Float precision**: use the serializers in `src/metrics.rs`
  (`serialize_f64_2dp`, `serialize_f64_5dp`, `serialize_f64_6dp`)
- **Scope**: all reads in file combined — no per-read-group, library, or sample
  breakdown
- **PF filtering**: statistics are computed on PF reads by default; field names
  do not carry a `pf_` prefix

## Code Style

**Priorities:** correctness > readability > performance.

- Write idiomatic Rust — don't fight the language
- Prefer meaningful names even if longer (`aligned_base_count` over `cnt`)
- Extract small/medium functions with clear inputs and outputs
- Keep related code together (metric struct next to its collector)
- Use `pub(crate)` for shared internal utilities
- Doc comments on all public items; code comments explain *why*, not *what*
- Don't over-generalize — solve the problem at hand

## Testing

**No checked-in test data.** All BAM, FASTA, and interval data is built
programmatically using the helpers in `tests/helpers/mod.rs`:

| Helper | Purpose |
|--------|---------|
| `SamBuilder` | Build `RecordBuf` instances and write to a temp BAM |
| `FastaBuilder` | Build a temp FASTA with `.fai` index |
| `read_metrics_tsv::<T>(path)` | Deserialize TSV output for assertions |
| `assert_float_eq!(a, b, eps)` | Float comparison with epsilon |

**Guidelines:**

- Integration tests go in `tests/test_<command>.rs`; unit tests in inline
  `#[cfg(test)]` modules
- Prefer many small, focused tests over parameterized/table-driven tests
- Test function, not implementation — tests should survive a significant refactor
- Cover expected results, error conditions, and boundary cases

## Reporting Bugs

Good bug reports are the fastest path to a fix. Please include:

1. **Version**: output of `riker --version`
2. **Command**: the exact command line that failed
3. **Error output**: full stderr/stdout
4. **OS**: operating system and version
5. **Expected vs. actual behavior**: what should happen vs. what did happen

### Creating Minimal Test Data

Full-size BAM files are often tens of gigabytes and can't be shared in an issue.
Instead, create a small file that still reproduces the problem:

**Extract a genomic region** (preferred — keeps all reads in a specific area):

```bash
samtools view -h -b input.bam chr1:1000-2000 > minimal.bam
samtools index minimal.bam
```

**Downsample** (when the bug isn't region-specific):

```bash
# Keeps ~0.1% of reads; seed=42 for reproducibility; preserves read pairs
samtools view -h -b -s 42.001 input.bam > minimal.bam
samtools index minimal.bam
```

**For reference FASTA**: extract just the relevant contigs, or describe the
genome build (e.g. GRCh38) and source so we can obtain it.

**For BED/interval files**: trim to only the intervals that trigger the bug.

Before submitting, **verify the bug still reproduces** on your minimal files:

```bash
riker <your-command> --input minimal.bam [other args]
```

**Common gotcha:** chromosome naming mismatches (`chr1` vs `1`) between your
BAM, reference, and interval files are a frequent source of errors — note which
style your files use.

## Performance

- Get it correct and well-tested first, then optimize
- Profile before tuning — use [samply](https://github.com/mstange/samply) or
  cargo-flamegraph to identify actual hot paths
- Verify correctness after any optimization (diff outputs against a baseline and
  run the full test suite)
