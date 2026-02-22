# Picard IntervalList File Format Specification

This document specifies the Picard IntervalList file format as defined by its
reference implementation in [htsjdk](https://github.com/samtools/htsjdk).  No
formal specification has been published elsewhere; this specification was derived
by reading the htsjdk source code and validated by testing edge cases with
Picard v3.4.0.

The key htsjdk source files that define the format are:

- `IntervalListCodec.java` — parsing logic
- `IntervalList.java` — data model and file I/O
- `IntervalListWriter.java` — output format
- `Interval.java` — interval data structure
- `Strand.java` — strand encoding


## Overview

An IntervalList file defines a set of genomic intervals against a reference
genome.  It is a plain-text, tab-delimited format consisting of a SAM-style
header section followed by interval data lines.

- **File extension**: `.interval_list` (by convention; not enforced)
- **Compression**: gzip and bgzip are supported when the file has a `.gz`
  extension (e.g., `.interval_list.gz`).  Picard/htsjdk uses the `.gz` extension
  to decide whether to decompress — compressed files without a `.gz` extension
  will fail.  Both gzip and bgzip content work since bgzip is gzip-compatible.
- **Character encoding**: UTF-8 (inherited from the SAM specification)
- **Line terminator**: `\n`; `\r\n` is tolerated by the SAM header parser
- **Created by**: The Picard / htsjdk project (Broad Institute)


## Header Section

The file begins with a **required** SAM-style header.  A file with no header
lines is invalid and will be rejected.

Header lines start with `@` and are parsed by htsjdk's `SAMTextHeaderCodec`,
which follows the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

### Required header line types

| Line type | Purpose | Required tags |
|-----------|---------|---------------|
| `@SQ` | Defines a reference sequence (contig) | `SN` (name), `LN` (length) |

At least one `@SQ` line must be present to define the sequences that intervals
may reference.  Additional `@SQ` tags (`AS`, `UR`, `M5`, `SP`) are optional.

### Optional header line types

| Line type | Purpose |
|-----------|---------|
| `@HD` | Header metadata (SAM version, sort order). Conventional but not required. |
| `@RG` | Read group definitions. Carried through but not used by interval processing. |
| `@PG` | Program records. Picard appends one when writing output. |
| `@CO` | Free-text comment lines. |

### Header termination

The header ends at the first line that does not start with `@`.  That line
becomes the first data line (or the file ends if there are no more lines).

**Blank lines must not appear before or within the header.**  A blank line
before the first `@` line causes Picard to report "Interval list file must
contain header."  A blank line between header lines (e.g., between `@HD` and
`@SQ`) terminates header parsing early — any `@SQ` lines after the blank are
never seen, causing sequence lookups to fail silently.

### Example header

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@CO	Exome capture targets, hg38
```


## Data Lines

Each data line represents one genomic interval.  Fields are separated by
**exactly one TAB character** (`\t`).  Spaces are not valid delimiters — a
space-separated line will be treated as a single field and rejected.

### Field definitions

Every data line must contain **exactly 5 fields**.  Lines with fewer or more
fields are rejected with a `TribbleException`.

| # | Field | Type | Description |
|---|-------|------|-------------|
| 1 | Sequence | String | Reference sequence name. Must match an `@SQ SN:` value in the header. |
| 2 | Start | Integer | Start position, **1-based, inclusive**. Must be ≥ 1. |
| 3 | End | Integer | End position, **1-based, inclusive**. See validation rules below. |
| 4 | Strand | Character | `+` (forward) or `-` (reverse). No other values are accepted. |
| 5 | Name | String | An identifier for the interval. May contain spaces. |

### Example data lines

```
chr1	1000	2000	+	target_001
chr1	5000	6000	-	target_002
chr2	100	500	+	target_003
```


## Coordinate System

Coordinates are **1-based and closed** (inclusive on both ends).

- The interval `chr1  1000  2000  +  foo` covers bases 1000 through 2000 on
  chr1, inclusive — a total of **1001 bases**.
- The minimum valid start position is **1** (not 0).
- A single-base interval has `start == end` (e.g., `chr1  500  500  +  snp`
  covers 1 base).

### Zero-length intervals

An interval where `start = end + 1` is valid and represents a zero-length
interval.  For example, `chr1  101  100  +  empty` is accepted by Picard and
reports 0 bases.  This is the *only* case where start may exceed end.

### Interval length formula

```
length = end - start + 1
```

For zero-length intervals, this yields 0 (since `start = end + 1`).


## Validation Rules

The following constraints are enforced by the htsjdk reference implementation.

1. **Header is required.**
   A file must begin with at least one `@`-prefixed header line.

2. **Exactly 5 fields per data line.**
   Each data line must contain exactly 5 tab-separated fields.

3. **Start ≥ 1.**
   The start position must be a positive integer (0 and negative values are
   invalid).

4. **Start ≤ End + 1.**
   Start may exceed end by at most 1 (a zero-length interval).  Any larger
   gap is invalid.

5. **End ≤ sequence length.**
   The end position must not exceed the sequence length declared in the
   corresponding `@SQ` header line.

6. **Strand must be `+` or `-`.**
   No other values are accepted, including `.` (period) and `*` (asterisk).

7. **Sequence must exist in header.**
   If the sequence name does not match any `@SQ SN:` entry, the interval is
   **silently skipped** with a warning.  This is notably *not* an error — the
   file continues to parse.

8. **Blank lines in the data section are skipped.**
   Empty or whitespace-only lines *after* the header are silently ignored.
   However, blank lines *before* or *within* the header are **not safe** — they
   truncate header parsing, causing the file to either be rejected as headerless
   or to lose `@SQ` definitions (see "Header termination" above).

9. **`#` is not a comment character.**
   Lines beginning with `#` are treated as data lines and will fail parsing
   (typically as a 1-field record).


## Output Format

The canonical writer (`IntervalListWriter`) produces output as follows:

1. The SAM header, including any `@PG` records added by the writing tool
2. One data line per interval: `{contig}\t{start}\t{end}\t{strand}\t{name}`
3. If an interval has a null/absent name, it is written as `.` (period)
4. Strand is always `+` or `-`

### Example output

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@PG	ID:1	CL:IntervalListTools ...	PN:IntervalListTools
chr1	1000	2000	+	target_001
chr1	5000	6000	-	target_002
chr2	100	500	+	.
```


## Complete Annotated Example

```
@HD	VN:1.6	SO:coordinate              ← SAM header line (optional)
@SQ	SN:chr1	LN:248956422               ← Sequence dictionary (required)
@SQ	SN:chr2	LN:242193529               ← One @SQ per reference contig
@CO	Exome capture targets, hg38        ← Comment (optional)
chr1	1000	2000	+	exon_BRCA1_1       ← 1001 bases on forward strand
chr1	5000	5000	+	snp_rs12345        ← Single base (start == end)
chr1	8001	8000	+	empty_interval     ← Zero-length (start == end+1)
chr2	100	500	-	exon_TP53_3            ← 401 bases on reverse strand
```


## Comparison with BED Format

| Property | IntervalList | BED |
|----------|-------------|-----|
| Coordinate system | 1-based, closed (inclusive both ends) | 0-based, half-open `[start, end)` |
| Header | SAM-style header, **required** | None required (track lines optional) |
| Fields per line | Exactly 5 | 3 minimum, up to 12+ |
| Delimiter | TAB only | TAB only |
| Strand field | Required (`+` or `-`) | Optional (field 6); allows `.` |
| Sequence validation | Against header dictionary | None inherent |
| Comment lines | `@CO` in header only; `#` not supported | `#` lines are comments |
| File extension | `.interval_list` | `.bed` |

### Coordinate conversion

| Direction | Start | End |
|-----------|-------|-----|
| IntervalList → 0-based half-open | `start_0 = start_1 - 1` | `end_halfopen = end_1` |
| 0-based half-open → IntervalList | `start_1 = start_0 + 1` | `end_1 = end_halfopen` |

For example, the IntervalList interval `chr1  1000  2000` (1-based closed)
corresponds to `chr1  999  2000` in BED (0-based half-open).  Both represent the
same 1001 bases.


## Riker Implementation Notes

Riker reads but does not write the IntervalList format.  While matching htsjdk
on many fronts, it is somewhat more lenient.

- **Format auto-detection.**  Riker checks the first non-empty line: if it
  starts with `@`, the file is treated as IntervalList; otherwise as BED.
  File extension is not considered.

- **Transparent decompression.**  Gzip and bgzip files are detected by magic
  bytes (`0x1f 0x8b`) and decompressed automatically, regardless of file
  extension.

- **Header is parsed and validated.**  The SAM header is parsed using noodles
  and the `@SQ` lines are validated as a **prefix** of the BAM file's sequence
  dictionary: same names, same lengths, same order.  The BAM may have additional
  contigs beyond those in the interval list, but not vice versa.

- **Individual interval validation.**  Each interval is validated against the
  parsed header: start must be ≥ 1, end must not exceed the sequence length,
  and the contig must exist in the interval list's own header.  All errors are
  collected and reported together.

- **Strand is ignored.**  The strand field (column 4) is consumed but discarded.
  Riker treats all intervals as unstranded.

- **Flexible field count.**  Riker does not require exactly 5 fields.  The name
  (field 5) and strand (field 4) are optional.  When absent, a name is
  synthesized from the coordinates (e.g., `chr1:999-2000`).

## References

- **htsjdk source**: https://github.com/samtools/htsjdk
  - `src/main/java/htsjdk/tribble/IntervalList/IntervalListCodec.java`
  - `src/main/java/htsjdk/samtools/util/IntervalList.java`
  - `src/main/java/htsjdk/samtools/util/IntervalListWriter.java`
  - `src/main/java/htsjdk/samtools/util/Interval.java`
  - `src/main/java/htsjdk/tribble/annotation/Strand.java`
- **SAM specification**: https://samtools.github.io/hts-specs/SAMv1.pdf
