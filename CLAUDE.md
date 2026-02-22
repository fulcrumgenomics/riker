# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Riker is a fast Rust CLI toolkit that ports key QC metrics tools from Picard/fgbio. It processes SAM, BAM, and CRAM files and produces sequencing QC metrics in TSV format — with lowercase headers, no metadata lines, and `frac_` instead of `pct_` naming conventions.

## Commands

```bash
# Build
cargo build
cargo build --release

# Test (uses nextest)
cargo ci-test                          # all tests
cargo nextest run <test_name>          # single test by name
cargo nextest run --test test_isize    # all tests in one integration test file

# Lint and format
cargo ci-lint                          # clippy with pedantic warnings as errors
cargo ci-fmt                           # check formatting (--check mode)
cargo fmt --package riker              # apply formatting
```

The `ci-*` aliases are defined in `.cargo/config.toml`. CI runs all three.

## Architecture

### Crate Structure

- `riker` — binary + library (`riker_lib`)
- `riker_derive` — proc-macro crate for `#[derive(MetricDocs)]` and `#[multi_options]`

The library is named `riker_lib` (see `[lib]` in `Cargo.toml`). External code and tests reference it as `riker_lib::...`.

### Command Pattern

Each subcommand is a struct implementing the `Command` trait (`src/commands/command.rs`). The `Subcommand` enum in `main.rs` dispatches to them via a `match` in `impl Command for Subcommand`.

To add a new command:
1. Create `src/commands/<name>.rs` with the command struct + options struct + collector + metric structs
2. Add `pub mod <name>;` to `src/commands/mod.rs`
3. Add a variant to `Subcommand` in `src/main.rs` and a match arm in `Command::execute()`
4. Integrate with the `multi` command (add `CollectorKind` variant, flatten `MultiMyOptions`, add `build_collectors` arm)
5. Register metric structs in `src/commands/docs.rs`

See the "Adding a New Command" section in `CONTRIBUTING.md` for the full walkthrough.

### Command Module Ordering

Command modules must follow this ordering convention:

1. **`use` statements**
2. **Constants and type aliases**
3. **Structs/enums, each immediately followed by all its impl blocks:**
   - Options struct (e.g. `ErrorOptions`) + `impl` + `impl Default`
   - Command struct (e.g. `Error`) + `impl` + `impl Command`
   - Collector struct (e.g. `ErrorCollector`) + `impl` + `impl Collector`
   - Metric structs (output row types)
   - Helper structs/enums — ordered higher-to-lower level (if A uses B, A comes first), then by importance to the overall implementation
4. **Module-level functions** (functions operating on primitives or external types)
5. **`#[cfg(test)] mod tests`**

**Impl block rules:**
- Consolidate all inherent methods into one `impl` block per struct; keep each trait impl as a separate block
- Trait impls go immediately after the struct's own `impl` block
- Within an impl block, order methods callers-before-callees (higher in the call stack first); constructors (`new`) always come first
- Functions that naturally belong to a type we own should be methods, not standalone functions

### Collector Pattern

Each command's core logic lives in a `Collector` struct implementing the `Collector` trait (`src/collector.rs`):

```rust
pub trait Collector: Send {
    fn initialize(&mut self, header: &Header) -> Result<()>;
    fn accept(&mut self, record: &RecordBuf, header: &Header) -> Result<()>;
    fn finish(&mut self) -> Result<()>;
    fn name(&self) -> &'static str;
}
```

Collectors store all configuration (output paths, reference handle, thresholds) as fields set at construction. The trait methods only receive the BAM header and records. This design enables the `multi` command to feed a single BAM pass to multiple collectors in parallel.

### Metrics Serialization

- All serialization utilities live in `src/metrics.rs`
- Per-tool metric structs live in their command file (e.g., `InsertSizeMetric` in `src/commands/isize.rs`)
- Float serializers: `serialize_f64_2dp`, `serialize_f64_5dp`, `serialize_f64_6dp`, `serialize_opt_f64_5dp`
- `write_tsv<T: Serialize>(path, rows)` writes tab-separated output with no metadata lines
- The `#[derive(MetricDocs)]` macro (from `riker_derive`) extracts doc comments from metric structs at compile time for the `riker docs` subcommand
- The `#[multi_options("prefix", "Heading")]` macro (from `riker_derive`) generates prefixed `Multi<Name>` companion structs for per-tool options in the `multi` command, with `validate()` for required-field checking

### Shared CLI Option Groups

Reusable `#[derive(clap::Args)]` structs in `src/commands/common.rs`:
- `InputOptions` — `--input`
- `ReferenceOptions` — `--reference` (required FASTA with .fai)
- `OptionalReferenceOptions` — `--reference` (optional FASTA with .fai)
- `IntervalOptions` — `--intervals` (IntervalList or BED)
- `DuplicateOptions` — `--include-duplicates`

### Key Source Modules

| Module | Purpose |
|--------|---------|
| `src/sam/alignment_reader.rs` | `open_alignment()` — opens SAM/BAM/CRAM via noodles |
| `src/sam/indexed_reader.rs` | `IndexedAlignmentReader` — indexed BAM/CRAM with region queries |
| `src/sam/record_filter.rs` | Filtering predicates (secondary, duplicate, PF, MAPQ) |
| `src/sam/pair_orientation.rs` | `PairOrientation` enum (FR/RF/Tandem) |
| `src/fasta.rs` | FASTA loading via noodles `IndexedReader`, wrapped in `Fasta` |
| `src/vcf.rs` | Indexed VCF/BCF reading, variant mask pre-loading |
| `src/intervals.rs` | Auto-detect IntervalList vs BED, parse to 0-based half-open |
| `src/progress.rs` | Progress logger (logs every N records) |
| `src/plotting.rs` | Plotting utilities and color palette |

### Error Command Architecture

The `error` command (`src/commands/error.rs`) has a unique dual execution path:

- **Standalone** (`riker error`): Uses `IndexedAlignmentReader` for per-region queries. Faster for solo runs.
- **Multi** (`riker multi --tools error`): Implements the `Collector` trait, receiving records sequentially. Loads full-contig `RegionContext` on contig transitions.

Performance-critical design choices (benchmarked):
- `CovariateValue` enum is deliberately 8 bytes — do not add larger variants without benchmarking
- Uses `rustc-hash` (fxhash) for HashMap keys — benchmarked 2.8× faster than std SipHash for covariate key types
- `StringInterner` avoids per-base string allocation for covariate values
- `BaseCovariateCache` computes each stratifier value at most once per base
- VCF variant masks are pre-loaded at initialization (not stored as `IndexedVcf` on the collector) to avoid `unsafe impl Send`

## Testing Conventions

- **No checked-in test data** — all BAM/FASTA/interval data is built programmatically
- Integration tests live in `tests/test_<command>.rs`; unit tests are inline `#[cfg(test)]` modules
- Test helpers in `tests/helpers/mod.rs`:
  - `SamBuilder` — builds `RecordBuf` instances and writes to a temp BAM (`to_temp_bam()` or `to_temp_indexed_bam()` with BAI index)
  - `FastaBuilder` — builds a temp FASTA with `.fai` index
  - `read_metrics_tsv::<T>(path)` — deserializes TSV output for assertions
  - `assert_float_eq!(a, b, eps)` macro for float comparisons

## Output Format Conventions

- TSV files: lowercase snake_case headers, tab-separated, no metadata/comment lines
- Fractions use `frac_` prefix (not `pct_`)
- PF-filtered reads don't use `pf_` field prefix
- No per-RG/library/sample breakdown (all reads in file combined)

## Code Safety

- `#![deny(unsafe_code)]` is enforced in both `src/lib.rs` and `src/main.rs` — avoid `unsafe` and find safe alternatives (e.g. pre-loading data rather than storing non-`Send` types on collectors)
