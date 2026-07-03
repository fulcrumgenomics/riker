# Riker Errata — Known and Intentional Differences

Riker aims to reproduce the metrics of Picard — and the other reference tools it ports (fgbio,
RSeQC, RNA-SeQC, Qualimap) — as closely as possible. This document records every known and
intentional functional difference, organized by command, so results can be reconciled precisely.
Each command's entry in the [README](README.md#differences-vs-picard-and-other-tools) links here.

Unless noted otherwise, a difference is a deliberate correction or design choice, not a regression.

---

## alignment

*Reproduces Picard `CollectAlignmentSummaryMetrics`.*

### `mean_aligned_read_length`

**Picard** computes mean aligned read length over all PF reads, including unmapped reads which contribute zero to the sum. The denominator is `PF_READS`, so unmapped reads dilute the average toward zero.

**riker** uses `aligned_reads` as the denominator, giving the mean length of reads that actually aligned.

**Impact:** riker produces slightly higher values than Picard — typically ~0-2bp for WGS data with high alignment rates. The difference grows with the fraction of unmapped reads.

### `reads_improperly_paired` / `frac_reads_improperly_paired`

**Picard** counts all mapped, paired, non-proper reads as improperly paired, including reads whose mate is unmapped (no `is_mate_unmapped()` guard).

**riker** requires the mate to also be mapped before counting a read toward `aligned_reads_in_pairs` or `reads_improperly_paired`. This avoids inflating the improper-pair count with reads whose mate simply failed to align.

**Impact:** riker reports fewer improperly paired reads than Picard. On typical WGS data the difference is usually small, in the single digit percent range. `aligned_reads_in_pairs` and `frac_reads_improperly_paired` are also affected.

---

## hybcap

*Reproduces Picard `CollectHsMetrics`.*

### `frac_exc_overlap` (Picard: `PCT_EXC_OVERLAP`)

riker reports a slightly lower `frac_exc_overlap` than Picard — typically below 1% relative difference (e.g. 0.032542 vs 0.032646 on a 1000G exome dataset).

**Cause:** Picard's overlap-clipping function (`SAMUtils.getNumOverlappingAlignedBasesToClip()` in htsjdk) does not verify that the read and its mate are mapped to the same contig before computing the overlap. It compares `getMateAlignmentStart()` against `getAlignmentStart()` as raw integers regardless of reference sequence. For chimeric read pairs — where one end maps to a different chromosome — Picard may compute a spurious overlap when the mate's coordinate on a different contig happens to be within the range of the reads aligned coordinates.  See https://github.com/broadinstitute/picard/issues/2039.

riker's equivalent function checks `reference_sequence_id != mate_reference_sequence_id` and returns 0 for chimeric pairs, correctly skipping overlap clipping when the reads are on different contigs.

**Impact:** `frac_exc_overlap` and `frac_exc_off_target` are affected; the latter because bases now correctly *not* categorized as overlapping are usually excluded subsequently for being off-target. The magnitude depends on the fraction of chimeric read pairs in the data, but is usually very small.

### Per-target GC content and GC dropout

riker and Picard compute per-target GC fraction differently when a target's reference sequence contains ambiguous (N) bases.

**Picard** counts N bases in the denominator but not the numerator: `gc_fraction = (G + C) / (A + T + G + C + N)`. This dilutes the GC fraction toward zero in proportion to the number of N bases in the target.

**riker** treats N bases as maximally uncertain, contributing 0.5 to the GC numerator and 1.0 to the denominator: `gc_fraction = (G + C + 0.5 * N) / (A + T + G + C + N)`. A target that is entirely N bases reports a GC fraction of 0.5 (maximum uncertainty) rather than 0.0.

Neither approach is fully correct — the true GC content of N bases is unknown — but riker's treatment avoids the bias toward low GC that Picard introduces. In practice the difference is only visible for targets where N bases make up a meaningful fraction of the reference sequence.

**Impact:** Per-target GC values can differ for N-containing targets, which in turn affects the GC bias curve and the `at_dropout` and `gc_dropout` summary metrics. All other per-target and summary metrics are unaffected.

---

## gcbias

*Reproduces Picard `CollectGcBiasMetrics`.*

### Summary metrics schema and read accounting

riker reworks the summary schema to make read accounting explicit and auditable. Picard's `GcBiasSummaryMetrics` exposes `ACCUMULATION_LEVEL`, `READS_USED`, `WINDOW_SIZE`, `TOTAL_CLUSTERS`, `ALIGNED_READS`, `AT_DROPOUT`, `GC_DROPOUT`, and `GC_NC_0_19 … GC_NC_80_100`, with `SAMPLE`/`LIBRARY`/`READ_GROUP` columns.

**Dropped fields.** riker emits a single summary row for the whole file, so `ACCUMULATION_LEVEL`, `LIBRARY`, and `READ_GROUP` are omitted (consistent with riker's no-per-read-group convention). `READS_USED` is also dropped: where Picard emits *two* rows — `READS_USED=ALL` and `READS_USED=UNIQUE` — when handling duplicates, riker selects duplicate inclusion with `--exclude-duplicates` and emits one row.

**Added fields.** riker adds an explicit read-accounting funnel:

- `total_reads` — every record except secondary and QC-fail (includes unmapped, duplicate, supplementary, and low-MAPQ reads).
- `aligned_reads` — the mapped subset of `total_reads`.
- `filtered_reads` — aligned reads not used in the GC computation (excluded as duplicate/supplementary, low MAPQ, in an excluded interval, or lacking a valid GC window), derived as `aligned_reads` minus the total binned read starts.
- `frac_filtered_reads` — `filtered_reads / aligned_reads`.

So `total_reads ≥ aligned_reads ≥ (aligned_reads − filtered_reads)`, the last being the reads actually binned. Note Picard's `ALIGNED_READS` in its `UNIQUE` row *excludes* duplicates, whereas riker's `aligned_reads` always counts all mapped reads and reports the excluded portion via `filtered_reads`.

**Renames** follow riker's global conventions: `GC_NC_x_y` → `gc_x_y_normcov`, and lowercase `frac_`-prefixed headers throughout.

### Read filtering

**Picard** processes every record (`SinglePassSamProgram` applies no record filter), so `TOTAL_CLUSTERS` and `ALIGNED_READS` include secondary, supplementary, and QC-fail reads.

**riker** discards secondary and QC-fail records entirely (counted nowhere), which avoids double-counting molecules represented by secondary alignments. Supplementary reads are included by default — they are real molecules — but can be dropped from the GC bins with `--exclude-supplementary`. Both tools count unmapped reads toward `total_clusters`/`TOTAL_CLUSTERS`, so riker matches Picard there.

**Impact:** for files containing secondary or QC-fail reads, riker's `total_clusters` and `aligned_reads` are lower than Picard's by the count of those reads.

### Read-to-window GC binning (forward-strand offset)

Both tools precompute, for every reference position, the %GC of the window starting there, then assign each read to a single %GC bin using the window at one computed position: the alignment start for forward-strand reads, and `alignment_end − window_size` for reverse-strand reads.

Picard stores these per-window GC values in a 0-based array (`gc[i]` is the window covering reference bases `[i, i + windowSize)`) but, for forward reads, indexes it with `SAMRecord.getAlignmentStart()`, which is **1-based**. A forward read therefore lands on `gc[start]` — the window beginning one base 3′ of the read's actual leftmost aligned base. Reverse reads use `getAlignmentEnd() − scanWindowSize`, which resolves to the intended window, so they are unaffected.

riker converts the alignment start to 0-based before indexing, so a forward read is binned by the window anchored exactly at its leftmost aligned base. riker's reverse-strand handling matches Picard's.

**Impact:** forward-strand reads are binned by windows that sit one base apart between the two tools; reverse-strand reads agree. This affects every forward read — not only N-containing or flagged regions — so it is the one difference present even with default flags and an N-free reference. Near a local GC transition a read can fall into an adjacent integer %GC bin, producing a small, systematic shift in the detail GC curve and in the `at_dropout` / `gc_dropout` and normalized-coverage values derived from it.

### `--min-mapq` and `--exclude-intervals` (vs. Picard PR #2030)

riker supports a `--min-mapq` threshold and an `--exclude-intervals` mask (BED or IntervalList) for skipping artifact regions such as poly-G stretches and adapter constructs. Released Picard's `CollectGcBiasMetrics` has neither; both correspond to the unreleased Picard PR [#2030](https://github.com/broadinstitute/picard/pull/2030), with two deliberate improvements over that PR:

- **Exclusion fixes the denominator.** The PR removes excluded reads from the numerator but leaves the reference-window distribution (the denominator) untouched, creating a numerator/denominator asymmetry. riker drops both the read *and* the reference window at an excluded position, keeping normalization consistent.
- **Exclusion is keyed on the read's computed start position**, not full-read overlap. This matches how reads are binned (by start position), is simpler, and lets users pad intervals when they want span semantics. (The PR uses `overlapsAny`.)

Both `--min-mapq` and `--exclude-intervals` are applied **per read, not per template**: each mate is evaluated independently, mirroring gcbias's per-read-start accounting, so one mate of a pair may be filtered while the other is still counted. riker does not drop a pair as a unit when only one mate fails a filter; widen the intervals to cover both mates, or pre-filter the BAM, if you need pair-level semantics.

**Impact:** these are opt-in; with neither flag set, riker's behavior matches released Picard (modulo the schema, read-filtering, and forward-strand binning differences above).

---

## error

*Reproduces Picard `CollectSamErrorMetrics`.*

### Reference N bases

**Picard** counts bases at reference-N positions toward `TOTAL_BASES` (and may count them as errors since the read base does not match N).

**riker** skips positions where the reference base is N, since errors cannot be meaningfully assessed without a known reference base.

**Impact:** riker reports slightly fewer `total_bases` and `error_bases` than Picard. The difference is proportional to the number of reference-N positions covered by reads — typically negligible for well-assembled references.

### Mismatch `total_bases` includes insertion bases in Picard

**Picard's** mismatch metric (`ERROR`) includes insertion bases in `TOTAL_BASES`. This happens because `BaseErrorCalculator.addBase()` increments `nBases` for both aligned (Match) bases and insertion bases, and `SimpleErrorCalculator` inherits this count as the denominator.

**riker** counts only aligned bases in the mismatch `total_bases`. Insertion bases are counted separately in the indel metric.

**Impact:** Picard's mismatch `TOTAL_BASES` is higher than riker's by the number of inserted bases passing filters. The mismatch `error_bases` (numerator) is identical between the two tools — only the denominator differs. riker's separation is arguably cleaner since insertion bases are not relevant to the substitution error rate.

### Insert size stratification

**Picard** caps insert size at `readLength * 10` for the `INSERT_LENGTH` stratifier. Reads with larger absolute insert sizes (e.g., chimeric pairs) are binned at the cap value.

**riker** excludes reads with absolute insert size above `--max-isize` (default 1000) from the `isize` stratifier entirely. These reads are still counted by all other stratifiers and in the `all` group. The threshold can be adjusted via `--max-isize`.

**Impact:** Picard may have a large bin at the cap value containing chimeric reads, while riker omits them. Overall error rates are unaffected.

### Insertions at read start

**Picard** counts insertions that occur before any aligned base in the read (e.g., CIGARs starting with `nI` or `nSnI`). The locus iterator attaches these to the preceding reference position.

**riker** skips insertions that have no preceding aligned base (no anchor position), since there is no reference context for stratification.

**Impact:** Picard reports slightly more insertions and inserted bases than riker. The difference is small — typically a few hundred events out of tens of thousands. Deletion counts are unaffected and match exactly between the tools.

### Q-score computation

**Picard** computes Q-scores using a Bayesian prior: `Q = -10 * log10((errors + 0.001) / (total + 1))`, rounded to the nearest integer. The prior (configurable via `PRIOR_Q`, default 30) prevents infinite Q-scores when there are zero errors.

**riker** computes Q-scores from the raw error rate: `Q = -10 * log10(errors / total)`, reported to two decimal places. When there are zero errors, a cap of Q99 is used.

**Impact:** Q-scores differ slightly due to the prior and rounding. The underlying counts (numerator and denominator) are comparable; only the derived Q-score differs.

### Stratifiers not ported from Picard

Picard's `CollectSamErrorMetrics` defines 32 stratifiers. riker ports 17 of them (`all`, `bq`, `mapq`, `cycle`, `read_num`, `strand`, `pair_orientation`, `isize`, `gc`, `read_base`, `ref_base`, `hp_len`, `pre_dinuc`, `post_dinuc`, `context_3bp`, `nm`, `indel_len`). The following Picard stratifiers are **not** available in riker:

- `PAIR_PROPERNESS` — whether the read is in a proper pair
- `HOMOPOLYMER` — the homopolymer base (A/C/G/T) at the current position
- `BINNED_HOMOPOLYMER` — homopolymer length bucketed into bins
- `BINNED_CYCLE` — machine cycle bucketed into bins
- `SOFT_CLIPS` — number of soft-clipped bases in the read
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

---

## rna

*Reproduces Picard `CollectRnaSeqMetrics` and fgbio `EstimateRnaSeqInsertSize`, and extends them in a
single pass with metrics drawn from RSeQC, RNA-SeQC, and Qualimap.*

riker's `rna` combines Picard `CollectRnaSeqMetrics` (base composition, strand specificity, 5'/3'
coverage bias) and fgbio `EstimateRnaSeqInsertSize` (insert size in transcript space) into one pass
over a single gene-model load, then adds the read-origin, gene-detection, splice-junction,
transcript-integrity, and per-biotype metrics described below. On identical inputs (same refFlat,
`--strand unstranded`, no ribosomal intervals) riker reproduces Picard's read/base accounting
essentially exactly (`PF_BASES` and `PF_ALIGNED_BASES` match to the base; base composition agrees
within ~0.5%); the differences below are intentional.

### Additional QC metrics beyond Picard/fgbio

These are folded into the same single pass and written to the wider `rna-metrics.txt` (plus the
separate biotype file). They are informed by RSeQC, RNA-SeQC, and Qualimap; where a gate or
definition differs from those tools it is noted.

- **Read accounting & genomic origin** (read-level): `total_reads`, `mapped_reads`, `unmapped_reads`,
  and the origin breakdown `exonic_reads` / `intronic_reads` / `intergenic_reads` / `ambiguous_reads`
  (exonic to >1 gene) / `ribosomal_reads`, with `assigned_reads` = reads uniquely assigned to one
  gene. Under default filters the four origin categories plus ribosomal **partition** the mapped
  reads, and `exonic_reads = assigned_reads + ambiguous_reads`. Fractions (`frac_*_reads`) are over
  `mapped_reads`.
- **Gene detection**: `genes_detected` — genes (by name, so a gene mapping to disjoint loci counts
  once) with at least `--genes-detected-min-reads` (default 5) uniquely-assigned reads.
- **Splice-junction annotation** (RSeQC `junction_annotation`): observed introns (CIGAR `N` gaps
  ≥ `--junction-min-intron`, default 50) classified against the annotated junctions built from the
  gene model. A junction is **known** if both splice sites match an annotated intron, **partial-novel**
  if exactly one site is annotated, **novel** otherwise. Reported both per-observation (`*_splice_obs`,
  counting read multiplicity) and as distinct junctions (`known_juncs` / `partial_juncs` /
  `novel_juncs`), plus `spliced_reads` / `frac_spliced_reads` and `frac_known_juncs`.
- **Transcript integrity (TIN)**: `median_tin` / `mean_tin` / `stddev_tin` over `tin_transcripts`
  transcripts — see the dedicated section below.
- **Per-biotype read counts**: a separate `<prefix>.rna-biotype.txt` (one row per gene biotype:
  `reads`, `frac_reads` over `assigned_reads`, most-assigned first), for spotting rRNA / globin /
  pseudogene contamination. refFlat carries no biotype column, so its rows are the inferred
  `protein_coding` / `noncoding` split (from CDS presence) — the only distinction recoverable without
  an annotated biotype. Fine biotypes (lincRNA, pseudogene, **rRNA**, miRNA …) cannot be recovered
  from refFlat, so ribosomal-from-biotype still requires `--ribosomal-intervals` for refFlat input.

### Transcript Integrity Number (TIN) — formula, interpretation, and gating

TIN is a per-transcript measure of how *evenly* a transcript is covered (0–100), an analog of the
Bioanalyzer RIN computed from the aligned reads. It originates in RSeQC `tin.py`. riker computes:

```
TIN = 100 · e^H / n,   where  H = −Σ pᵢ·ln pᵢ,   pᵢ = cᵢ / Σⱼ cⱼ
```

over the `n` exonic positions of the transcript (`cᵢ` = read coverage at position `i`). `e^H` is the
*perplexity* of the coverage distribution — the effective number of equally-covered positions —
so dividing by `n` gives "what fraction of the transcript is covered as evenly as if it were
uniform." Equivalently:

```
TIN = 100 · exp( −D_KL(normalized coverage ‖ uniform) )
```

i.e. **TIN measures how far the coverage *shape* is from flat**, and is **scale-invariant** — it
depends only on the relative coverage profile, not on absolute depth. Perfectly uniform coverage
gives `e^H = n` and TIN = 100; coverage concentrated on a sub-region drives it down.

**What does and does not move TIN.** Two behaviors follow directly and are worth internalizing:

- **Smooth spikes and dips, with every base covered → TIN stays high.** Ordinary multiplicative
  wobble (say 0.5×–1.5× the mean everywhere) is a small KL divergence, so TIN remains in the high
  90s. TIN does **not** punish a transcript for being merely bumpy.
- **Zero / near-zero gaps dominate.** Uncovered (or relatively near-zero) stretches act like missing
  positions: `e^H` falls toward the number of well-covered positions, so `TIN ≈ 100 · (well-covered
  fraction)`. A transcript with ~10% of its length in a deep valley scores ≈ 90; ~20% ≈ 80; ~5% ≈ 95.
  This is the 3′-degradation signal TIN is designed to capture — a degraded mRNA loses its 5′ end,
  creating exactly such a gap.

**Scale-invariance caveat (important for interpretation).** Because only the *relative* shape
matters, a region at, say, 1–5× coverage inside a transcript whose mean is 500× is ~1% of the mean —
which the entropy treats as **essentially empty**, even though 1–5× is perfectly adequate to call
expression. Such a transcript's TIN drops roughly in proportion to the length of that low region
(10% of the length → TIN ≈ 90), regardless of the high absolute depth elsewhere. So a "low" TIN on a
deeply-sequenced gene means **relative unevenness**, *not* inadequate coverage. TIN is a
relative-uniformity metric, not an absolute-adequacy one — read it accordingly.

**Behavior with multiple isoforms.** TIN is computed on a single transcript model, so co-expressed
isoforms are a potential confounder, but a bounded one:

- If two isoforms share most exons and each has a few unique exons, and **both are expressed**, the
  shared exons receive reads from both isoforms while the unique exons receive reads from one — a
  coverage **step** (e.g. ~2× on shared vs. unique exons), but **no gaps**. The TIN penalty is mild
  (≈ 2–3 points for a 2× step).
- If **only one isoform is expressed**, that isoform's model is fully covered (TIN ≈ 100), while the
  other isoform's model has its unique exons at zero (a gap → low TIN). riker computes TIN on **one
  representative transcript per gene — the most-expressed (highest mean coverage) one** — which
  naturally selects the fully-covered isoform and avoids the gappy unexpressed one. This makes riker
  *more* isoform-robust than computing TIN over all transcripts (as RSeQC does), since the latter's
  median is dragged down by unexpressed-isoform artifacts.

**Transcript selection and gating.** A transcript contributes to the TIN statistics only if:

1. it is the **highest-mean-coverage** transcript of its gene (one representative per gene), and
2. its transcribed length is at least `--minimum-length` (default 500), and
3. its **mean** per-base coverage exceeds `--tin-min-coverage` (default 10).

The mean-coverage gate (rather than peak depth) is deliberate: a narrow pileup barely moves the
whole-transcript mean, so it cannot sneak a noisy, sparsely-covered transcript past the gate, while a
genuinely deeply-but-unevenly covered transcript is kept. Length and best-per-gene selection mirror
how riker already chooses transcripts for the coverage/bias metrics, keeping the two consistent.

**Divergence from RSeQC (the field is not standardized here).** RSeQC `tin.py` gates on the number of
**unique read-start positions** (> `minCov`, default 10) and computes entropy over a *strided sample*
of ~100 exonic positions; it scores **every** transcript. riker instead uses a mean-coverage +
length gate, computes entropy over the **full** transcript (no sampling), and scores **one transcript
per gene**. These choices place riker's TIN in the same family as the coverage/bias metrics of
Picard, RNA-SeQC, and riker itself (mean-coverage + length + best-per-gene selection), rather than
RSeQC's read-start gate. Consequently absolute TIN values will not match RSeQC to the decimal, though
the metric is the same in spirit (100 = perfectly even, lower = more degraded/uneven). The number of
contributing transcripts (`tin_transcripts`) will typically be far smaller than RSeQC's, because the
mean > 10 gate restricts to well-expressed genes — by design, since entropy estimates from sparse
coverage are noisy.

**Scale difference (validated, and deliberate).** In practice riker's median TIN runs roughly
**14–16 points higher** than RSeQC's / RustQC's on the same data — because riker scores one
well-covered transcript per gene while RSeQC scores the full isoform set (whose marginal members pull
the median down). This is the single largest absolute numeric difference between riker `rna` and the
reference tools, so it is worth stating plainly: **the two report the same shape on a different
scale.** We validated this on a controlled RNA-degradation ladder (Sigurgeirsson et al. 2014, RIN
10→2 from one RNA source) and a depth-downsampling sweep: both implementations track degradation
monotonically with very similar sensitivity (riker 88→41, RSeQC 74→26 across the ladder), and both
are robust to sequencing depth, with riker holding markedly steadier (≈0.5% vs ≈4% TIN change across
a 5× depth reduction). The upshot for users: read TIN on each tool's own scale and compare
like-for-like; do not expect `riker` and `RSeQC` to agree on the absolute number. Full methodology,
data, and plots are in [docs/rna-integrity.md](docs/rna-integrity.md).

### Gene-model formats and contig naming

**Picard** accepts only UCSC refFlat (`REF_FLAT`).

**riker** accepts refFlat, GTF, and GFF3 (GENCODE / Ensembl / RefSeq) and **auto-detects** the format
from file contents. Provider conventions are handled without a flag: attribute keys are tried in
order (`gene_type`/`gene_biotype`/`biotype`, `gene_name`/`Name`/`gene_id`), and contig names are
reconciled to the BAM header by adding/stripping a `chr` prefix, aliasing `MT`↔`chrM`, and mapping
RefSeq accessions (e.g. `NC_000001.11`) to common names via the GFF `region` features.

### Ribosomal regions: biotype-derived plus explicit intervals

**Picard** marks ribosomal bases only from a separate `RIBOSOMAL_INTERVALS` interval_list; without
it, `RIBOSOMAL_BASES` is null.

**riker** derives ribosomal territory from the gene-model biotype (`rRNA` / `Mt_rRNA` /
`rRNA_pseudogene`) when present, takes the **union** with an optional `--ribosomal-intervals` file
(BED or IntervalList), and leaves the ribosomal metrics blank only when neither source exists. Note
that biotype/annotation-derived rRNA is a **lower bound** — GRCh38 collapses the rDNA tandem arrays,
so for an accurate rRNA fraction supply curated intervals (or align to a reference with a dedicated
rDNA contig). refFlat carries no biotype, so refFlat input needs `--ribosomal-intervals`.

### Strand-specificity auto-detection

**Picard** requires `STRAND_SPECIFICITY` (`NONE` / `FIRST_READ_TRANSCRIPTION_STRAND` /
`SECOND_READ_TRANSCRIPTION_STRAND`).

**riker** defaults to `--strand auto`, inferring `forward` / `reverse` / `unstranded` from the
observed read-1 vs read-2 transcription-strand counts, and always reports the
`r1_tx_strand_reads`/`r2_tx_strand_reads` and `frac_r1_tx_strand_reads`/`frac_r2_tx_strand_reads`
regardless of the setting. `--strand {unstranded,forward,reverse}` forces a specific model.

### Transcript-space insert size (from fgbio EstimateRnaSeqInsertSize)

Not part of Picard's tool. riker computes the fragment insert size in **transcript space** — the
5'-to-5' distance summed over the overlapping exons, so introns are collapsed — per pair orientation
(FR / RF / tandem). It needs the mate's reference span, obtained one of two ways:

- **`MC` (mate CIGAR) tag — exact, any orientation.** When present, the mate's exact blocks are
  used; this is the recommended input (add it with `samtools fixmate -m`).
- **TLEN fallback — FR pairs only, when `MC` is absent.** Evaluated at the forward (leftmost) read,
  the mate's right end is `start + |TLEN| - 1`, so the pair span and both 5' ends are known even
  without the mate's CIGAR. Because the mate's internal splicing is then unknown, the pair is kept
  only when both mate ends fall in an exon and the mate's exonic footprint clears the same overlap
  threshold (a properly spliced mate still passes — the intron drops out of the exonic length). This
  is a best guess, not exact: it assumes a **spec-compliant TLEN** and roughly equal mate read
  length. A BAM with neither `MC` nor a valid `TLEN` produces no (or incorrect) insert sizes — so
  `samtools fixmate -m` remains the way to guarantee correctness.

The size is never computed from read length alone.

riker accepts a pair when **exactly one gene encloses the whole pair span** (and the read and mate
clear the per-transcript exon-overlap threshold, with all transcripts agreeing on the size). fgbio
instead requires that exactly one gene *overlaps* the pair span, so a non-enclosing bystander gene
near the pair causes fgbio to drop it. riker's rule keeps those pairs — they are still unambiguously
inside a single gene — which on real data uses a few percent more pairs (≈ +7% on our test set) with
a negligible shift in the distribution (mean and median move < 1–2 bp).

### Coordinate-sorted input required

**Picard / fgbio** query an interval tree per read and do not require a particular record order.

**riker** walks the gene model with a coordinate sweep, so it **requires a coordinate-sorted**
BAM/CRAM (`@HD SO:coordinate`) and fails fast otherwise (sort with `samtools sort`). This is what
makes the single-pass locus lookup allocation-free and is satisfied by essentially all aligned
RNA-seq data in practice.

### Coverage and bias fixes

Picard's coverage-based metrics (`MEDIAN_CV_COVERAGE`, `MEDIAN_5PRIME_BIAS`, `MEDIAN_3PRIME_BIAS`,
`MEDIAN_5PRIME_TO_3PRIME_BIAS`) contain several off-by-one issues that riker corrects, so these
values differ by ~1% in our testing (always in the corrected direction):

- **Top-N selection**: Picard keeps the top *1001* most-expressed transcripts
  (`coverages[length - 1001]`); riker keeps the top **1000**.
- **Short transcripts**: Picard's 5'/3' windows are each `end_bias_bases` wide but it only requires
  transcripts ≥ `max(minimum_length, end_bias_bases)`, so the windows overlap for transcripts
  shorter than `2 × end_bias_bases`. riker **excludes** those transcripts from the CV/bias set
  instead of changing the formula.
- **Per-base coverage**: Picard's coverage accumulator iterates a half-open range against an
  inclusive end, dropping the last base of every alignment block; riker counts **every** aligned base.
- **Normalized coverage-by-position**: Picard's per-percent windows overlap at the boundaries,
  double-counting boundary bases; riker samples 101 **non-overlapping** normalized positions.

### Ribosomal fragment overlap

**Picard** takes the single largest rRNA-interval intersection per fragment (`Math.max`) and derives
the fragment end from the inferred insert size (`TLEN`), which is degenerate when `TLEN` is 0.

**riker** sums the **union** of overlaps with the (merged, disjoint) ribosomal intervals, and derives
the fragment end from the mate's `MC` tag, falling back to `TLEN` (never read length).

### Gene → locus → transcript model

**Picard** keys genes by name, treating all transcripts of a name as one gene spanning their extent.

**riker** groups a gene's transcripts into **loci** (same contig/strand, overlapping spans), so a gene
mapping to disjoint locations becomes independent loci. The strand/template-strand metrics count a
read only when it overlaps exactly **one locus** (vs Picard's one *gene*), which moves a small number
of reads near multi-locus genes to `unexplained_reads`; the reported strand *fractions* are
unaffected.

### Filtering, duplicates, and accumulation levels

- **Duplicates** are counted by default (matching Picard, which never filters them); riker adds
  `--exclude-duplicates` and always reports an observed `duplicate_rate`.
- **`--min-mapq`** (default 0, matching Picard's lack of a MAPQ filter) gates the alignment-based
  metrics but not the `bases` (PF) total. Supplementary reads are counted in base metrics but
  excluded from the strand/template and read-origin metrics, as in Picard.
- **`frac_usable_bases`** keeps Picard's intentional denominator — total sequenced bases (`bases`),
  not aligned bases — so it remains directly comparable to Picard's `PCT_USABLE_BASES`.
- riker reports **all reads combined**; Picard's `METRIC_ACCUMULATION_LEVEL` (per sample / library /
  read group) is not replicated.
