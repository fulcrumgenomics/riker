# Staging rules — fetch raw inputs from the public internet and produce
# a uniform stage/{sample}/input.bam for every sample, transcoding CRAM
# to BAM (and deleting the CRAM after) so all benchmarks read BAM. All
# rules here reserve bench=1 — they can run during the DAG fill-in but
# are excluded by the bench=100 lock during timed runs.

# ---- Reference fetch -----------------------------------------------------

rule fetch_reference:
    """Download a reference FASTA from its public URL, gunzip if needed,
    and produce both a .fai (for samtools/riker) and a .dict (for picard
    CollectHsMetrics' bait/target interval-list conversion)."""
    output:
        fasta = f"{STAGE_DIR}/refs/{{ref}}.fa",
        fai   = f"{STAGE_DIR}/refs/{{ref}}.fa.fai",
        dict  = f"{STAGE_DIR}/refs/{{ref}}.dict",
    wildcard_constraints:
        ref = "|".join(REFERENCES.keys()),
    params:
        url = lambda w: REFERENCES[w.ref]["fasta_url"],
        md5 = lambda w: REFERENCES[w.ref].get("md5") or "",
    threads: 4
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        url={params.url:q}
        out_fa={output.fasta:q}
        mkdir -p "$(dirname "$out_fa")"

        # Download to a temp file (or directly to .fa for non-.gz URLs).
        # Keeping curl/aws and decompress separate so curl's --retry can
        # actually re-pull on a mid-stream error rather than silently
        # corrupting a piped gunzip.
        if [[ "$url" == *.gz ]]; then
            tmp="${{out_fa}}.gz.partial"
        else
            tmp="${{out_fa}}.partial"
        fi
        rm -f "$tmp"

        case "$url" in
            s3://*)
                aws s3 cp --no-sign-request "$url" "$tmp"
                ;;
            http://*|https://*|ftp://*)
                curl -sSfL --retry 10 --retry-all-errors --retry-delay 10 \
                     "$url" -o "$tmp"
                ;;
            *)
                echo "ERROR: unsupported reference URL scheme: $url" >&2
                exit 1
                ;;
        esac

        if [[ "$url" == *.gz ]]; then
            gunzip -c "$tmp" > "$out_fa"
            rm -f "$tmp"
        else
            mv "$tmp" "$out_fa"
        fi

        # Optional md5 check. Empty md5 means "not pinned" — the check is
        # skipped. To pin, edit references.yaml manually with the value of
        # `md5sum stage_dir/refs/<ref>.fa` after a successful first fetch;
        # subsequent runs will then fail loudly on any mismatch.
        #
        # Bind to a bash variable first: Snakemake's :q formatter renders
        # an empty string as nothing at all (not as ''), so inlining
        # {{params.md5:q}} into a [[ ... ]] test would yield a bash syntax
        # error when md5 is unset. Assigning to a quoted variable avoids
        # the trap and lets us test cleanly.
        md5_expected={params.md5:q}
        if [[ -n "$md5_expected" ]]; then
            actual=$(md5sum "$out_fa" | awk '{{print $1}}')
            if [[ "$actual" != "$md5_expected" ]]; then
                echo "ERROR: md5 mismatch on reference {wildcards.ref}" >&2
                echo "  expected: $md5_expected" >&2
                echo "  actual:   $actual" >&2
                exit 1
            fi
        fi

        samtools faidx "$out_fa"
        samtools dict "$out_fa" -o {output.dict:q}
        """


# ---- Capture-kit interval_list fetch + BED derivation -------------------

rule fetch_kit_intervals:
    """Download a capture kit's intervals file and rewrite its @SQ
    header from the kit's declared reference's .dict, so Picard
    CollectHsMetrics' strict sequence-dictionary equality check
    against the BAM passes. The upstream interval_list (e.g. the
    GIAB-shared b37 list) often has a slightly different @SQ set than
    the BAM's hs37d5 alignment dictionary — picard rejects any size
    mismatch even when the canonical contigs are identical.

    Riker hybcap is permissive about this and works with the original
    headers, but we need the rewritten file for picard."""
    input:
        ref_dict = lambda w: f"{STAGE_DIR}/refs/{KITS[w.kit]['reference']}.dict",
    output:
        intervals = f"{STAGE_DIR}/kits/{{kit}}/intervals.list",
    wildcard_constraints:
        kit = "|".join(KITS.keys()) if KITS else "DOESNOTMATCH",
    params:
        url = lambda w: KITS[w.kit]["intervals_url"],
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.intervals:q})"
        raw="$(dirname {output.intervals:q})/raw.list"
        url={params.url:q}
        case "$url" in
            s3://*)        aws s3 cp --no-sign-request "$url" "$raw" ;;
            http*|ftp://*) curl -sSfL --retry 10 --retry-all-errors --retry-delay 10 "$url" -o "$raw" ;;
            *) echo "ERROR: unsupported intervals URL scheme: $url" >&2; exit 1 ;;
        esac
        # Reject empty / @-only files so a malformed upstream doesn't produce
        # a silent header-only intervals.list (downstream tools would then
        # benchmark zero capture regions). Count data rows with `grep -vc`
        # rather than `grep -v ... | grep -q .`: under `set -o pipefail`, the
        # `grep -q` short-circuits and SIGPIPEs the upstream grep (exit 141)
        # on a LARGE kit file, so the negated pipeline falsely fires the guard.
        if [[ ! -s "$raw" ]] || [[ "$(grep -vc '^@' "$raw")" -eq 0 ]]; then
            echo "ERROR: kit intervals file empty or header-only: $raw" >&2
            exit 1
        fi
        # Concatenate the reference's .dict (provides @HD + matching @SQ
        # headers) with the data rows from the upstream interval_list.
        cat {input.ref_dict} <(grep -v '^@' "$raw") > {output.intervals:q}
        rm -f "$raw"
        """


rule kit_intervals_to_bed:
    """Convert a Picard interval_list to BED (0-based half-open) for
    mosdepth. Picard interval_list rows are: chr<TAB>1-based start<TAB>
    1-based end<TAB>strand<TAB>name. BED rows are: chr<TAB>0-based
    start<TAB>0-based-exclusive end."""
    input:  intervals = f"{STAGE_DIR}/kits/{{kit}}/intervals.list"
    output: bed       = f"{STAGE_DIR}/kits/{{kit}}/intervals.bed"
    wildcard_constraints:
        kit = "|".join(KITS.keys()) if KITS else "DOESNOTMATCH",
    resources: bench=1
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}} /^@/ {{ next }} {{ print $1, $2 - 1, $3 }}' \
            {input.intervals:q} > {output.bed:q}
        """


# ---- Source download (CRAM or BAM) --------------------------------------

# Fetch a CRAM source. Runs once per sample whose source_url ends in .cram.
# Output is intentionally a "scratch" file — the transcode rule deletes it
# after producing input.bam.
rule fetch_source_cram:
    output: cram = f"{STAGE_DIR}/{{sample}}/source.cram"
    wildcard_constraints:
        sample = "|".join(CRAM_SOURCE_SAMPLES) if CRAM_SOURCE_SAMPLES else "DOESNOTMATCH",
    params:
        url = lambda w: SAMPLES.loc[w.sample, "source_url"],
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.cram:q})"
        url={params.url:q}
        case "$url" in
            s3://*) aws s3 cp --no-sign-request "$url" {output.cram:q} ;;
            http*|ftp://*) curl -sSfL --retry 5 --retry-delay 10 "$url" -o {output.cram:q} ;;
            *) echo "ERROR: unsupported source URL scheme: $url" >&2; exit 1 ;;
        esac
        """


# Fetch a BAM source directly to input.bam. Used for samples whose
# source_url ends in .bam (the GIAB Oslo exomes). Index produced via
# samtools index in a separate rule that's reused by the transcode and
# downsample paths too.
rule fetch_source_bam:
    output: bam = f"{STAGE_DIR}/{{sample}}/input.bam"
    wildcard_constraints:
        sample = "|".join(BAM_SOURCE_SAMPLES) if BAM_SOURCE_SAMPLES else "DOESNOTMATCH",
    params:
        url = lambda w: SAMPLES.loc[w.sample, "source_url"],
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam:q})"
        url={params.url:q}
        case "$url" in
            s3://*) aws s3 cp --no-sign-request "$url" {output.bam:q} ;;
            http*|ftp://*) curl -sSfL --retry 5 --retry-delay 10 "$url" -o {output.bam:q} ;;
            *) echo "ERROR: unsupported source URL scheme: $url" >&2; exit 1 ;;
        esac
        """


# ---- CRAM -> BAM transcode ----------------------------------------------

rule transcode_cram_to_bam:
    """Convert a CRAM to BAM using the sample's reference, then delete
    the CRAM. Output is the canonical stage/{sample}/input.bam consumed
    by every benchmark rule. Per the design spec we always benchmark BAM
    so that mosdepth/picard/qualimap aren't penalized by htsjdk's slower
    CRAM decoder."""
    input:
        cram = f"{STAGE_DIR}/{{sample}}/source.cram",
        ref  = lambda w: ref_for_sample(w.sample),
    output:
        bam  = f"{STAGE_DIR}/{{sample}}/input.bam",
    wildcard_constraints:
        sample = "|".join(CRAM_SOURCE_SAMPLES) if CRAM_SOURCE_SAMPLES else "DOESNOTMATCH",
    threads: lambda w: int(config["staging_threads"])
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        samtools view -@ {threads} -b -T {input.ref} -o {output.bam} {input.cram}
        # Disk hygiene per spec: remove the CRAM so we don't double-store.
        rm -f {input.cram}
        """


# ---- Downsample one BAM into another ------------------------------------

rule downsample_bam:
    """Subsample a parent BAM by `--subsample <frac> --subsample-seed <N>`.
    Used for the HG02675 30x/20x/15x/4x depth series."""
    input:
        bam = lambda w: f"{STAGE_DIR}/{SAMPLES.loc[w.sample, 'downsample_from']}/input.bam",
    output:
        bam = f"{STAGE_DIR}/{{sample}}/input.bam",
    wildcard_constraints:
        sample = "|".join(DOWNSAMPLE_SAMPLES) if DOWNSAMPLE_SAMPLES else "DOESNOTMATCH",
    params:
        frac = lambda w: SAMPLES.loc[w.sample, "downsample_frac"],
        seed = lambda w: SAMPLES.loc[w.sample, "downsample_seed"],
    threads: lambda w: int(config["staging_threads"])
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        samtools view -@ {threads} -b \
            --subsample {params.frac} --subsample-seed {params.seed} \
            {input.bam} -o {output.bam}
        """


# ---- BAM index -----------------------------------------------------------

rule index_bam:
    """Produce a .bai sidecar for every staged BAM. samtools index is
    cheap and the same code path serves the transcoded, downsampled, and
    direct-fetch BAMs."""
    input:  bam = f"{STAGE_DIR}/{{sample}}/input.bam"
    output: bai = f"{STAGE_DIR}/{{sample}}/input.bam.bai"
    threads: lambda w: int(config["staging_threads"])
    resources: bench=1
    shell:
        "samtools index -@ {threads} {input.bam}"


# ---- RNA gene-model (GTF) fetch + refFlat derivation --------------------
# The RNA BAMs themselves reuse fetch_source_bam (their source_url ends in
# .bam) and index_bam — the ENCODE STAR alignments are already coordinate-
# sorted, which is all riker rna + rustqc require. No markdup/fixmate: the
# RustQC subset drops dupRadar (its only marked-dup consumer) and riker uses
# its non-MC insert-size fallback.

rule fetch_gene_model_gtf:
    """Download a gene-model GTF (gzip). Read directly by riker rna and
    rustqc; also the source for the Picard refFlat + ribosomal interval_list
    derived below."""
    output:
        gtf = f"{STAGE_DIR}/annotations/{{annotation}}/genes.gtf.gz",
    wildcard_constraints:
        annotation = "|".join(ANNOTATIONS.keys()) if ANNOTATIONS else "DOESNOTMATCH",
    params:
        url = lambda w: ANNOTATIONS[w.annotation]["gtf_url"],
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gtf:q})"
        url={params.url:q}
        case "$url" in
            s3://*)        aws s3 cp --no-sign-request "$url" {output.gtf:q} ;;
            http*|ftp://*) curl -sSfL --retry 10 --retry-all-errors --retry-delay 10 "$url" -o {output.gtf:q} ;;
            *) echo "ERROR: unsupported annotation URL scheme: $url" >&2; exit 1 ;;
        esac
        # Sanity: must be a valid non-empty gzip (guards against a mirror
        # serving an HTML error page with a 200).
        gzip -t {output.gtf:q}
        """


rule gtf_to_refflat:
    """Convert a GTF to Picard refFlat via UCSC gtfToGenePred. Picard
    CollectRnaSeqMetrics needs refFlat; riker + rustqc read the GTF directly,
    so all three share one annotation source."""
    input:
        gtf = f"{STAGE_DIR}/annotations/{{annotation}}/genes.gtf.gz",
    output:
        refflat = f"{STAGE_DIR}/annotations/{{annotation}}/refFlat.txt",
    wildcard_constraints:
        annotation = "|".join(ANNOTATIONS.keys()) if ANNOTATIONS else "DOESNOTMATCH",
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        d="$(dirname {output.refflat:q})"
        # gunzip -c (not zcat): zcat on macOS only handles .Z, so this keeps
        # the recipe portable to a non-Linux dev box.
        gunzip -c {input.gtf:q} > "$d/genes.gtf"
        # genePredExt columns 1-10 are name,chrom,strand,txStart,txEnd,
        # cdsStart,cdsEnd,exonCount,exonStarts,exonEnds; column 12 is name2
        # (the gene name, via -geneNameAsName2). refFlat = gene name followed
        # by those first 10 columns.
        gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons \
            "$d/genes.gtf" "$d/genes.genePred"
        awk 'BEGIN{{OFS="\t"}} {{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' \
            "$d/genes.genePred" > {output.refflat:q}
        rm -f "$d/genes.gtf" "$d/genes.genePred"
        """


rule build_ribosomal_intervals:
    """Build a Picard ribosomal interval_list for CollectRnaSeqMetrics from
    the sample BAM's sequence dictionary (@HD/@SQ) plus the rRNA / Mt_rRNA
    gene loci in the annotation GTF. Per-sample because the @SQ set comes
    from the BAM header. riker rna auto-derives the same rRNA loci from the
    GTF's biotypes, so this gives Picard the equivalent ribosomal input."""
    input:
        bam = f"{STAGE_DIR}/{{sample}}/input.bam",
        gtf = lambda w: gene_model_gtf_for_sample(w.sample),
    output:
        intervals = f"{STAGE_DIR}/{{sample}}/ribosomal.interval_list",
    wildcard_constraints:
        sample = "|".join(RNA_SAMPLES) if RNA_SAMPLES else "DOESNOTMATCH",
    resources: bench=1
    shell:
        r"""
        set -euo pipefail
        # The SAM sequence dictionary (@HD + @SQ) becomes the interval_list header.
        samtools view -H {input.bam:q} | grep -E '^@HD|^@SQ' > {output.intervals:q}
        # Append rRNA / Mt_rRNA gene loci as interval_list data rows:
        #   chrom<TAB>start(1-based)<TAB>end(1-based)<TAB>strand<TAB>gene_id
        # GTF is 1-based inclusive, matching Picard interval_list coordinates.
        # gene_id is extracted portably (POSIX split, not gawk match(...,arr)).
        # gunzip -c, never zcat (zcat is .Z-only on macOS / non-portable).
        gunzip -c {input.gtf:q} | awk -F'\t' '
            $3 == "gene" && $9 ~ /gene_type "(rRNA|Mt_rRNA)"/ {{
                split($9, a, "gene_id \"");
                split(a[2], b, "\"");
                print $1"\t"$4"\t"$5"\t"$7"\t"b[1]
            }}' >> {output.intervals:q}
        # Fail loudly if the annotation yielded zero rRNA rows (Picard would
        # otherwise silently report all-zero ribosomal).
        if ! grep -qv '^@' {output.intervals:q}; then
            echo "ERROR: no rRNA/Mt_rRNA gene loci found in {input.gtf}" >&2
            exit 1
        fi
        """
