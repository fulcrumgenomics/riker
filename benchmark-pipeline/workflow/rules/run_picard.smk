# Picard benchmark rules — one rule per Picard tool variant. Each rule is
# disjoint by `tool` wildcard; the path schema results/run/{sample}/{profile}/
# {tool}/rep{rep}/ keeps outputs separated.
#
# Picard rules:
#   picard-cwm   CollectWgsMetrics — used for both wgs-only AND wgs-bundle
#                (in wgs-bundle it's the supplementary tool not covered by
#                CollectMultipleMetrics).
#   picard-chsm  CollectHsMetrics — same dual role for hybcap.
#   picard-cmm   CollectMultipleMetrics — only the bundle profiles. For
#                wgs-bundle, adds PROGRAM=CollectGcBiasMetrics so it
#                covers everything riker multi does except CollectWgsMetrics.
#                For hybcap-bundle, the default 5 PROGRAMs cover what
#                riker multi adds beyond hybcap itself.
#
# Threading: per spec, picard runs at 1 thread for the single-tool
# profiles AND for supplementary picard-cwm/-chsm in bundle profiles.
# picard-cmm in bundle profiles gets USE_ASYNC_IO_*=true to be at least
# somewhat fair against riker --threads 4.

# ---- CollectWgsMetrics ---------------------------------------------------

rule run_picard_cwm:
    """CollectWgsMetrics. The single-threaded WGS comparator in the fair
    three-way (wgs-t1), plus the supplementary 'CWM on top of CMM' datum for
    wgs-bundle."""
    input:
        bam = f"{STAGE_DIR}/{{sample}}/input.bam",
        bai = f"{STAGE_DIR}/{{sample}}/input.bam.bai",
        ref = lambda w: ref_for_sample(w.sample),
        fai = lambda w: ref_for_sample(w.sample) + ".fai",
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cwm/rep{{rep}}/time.txt",
    log:
        # See note on log: vs output: in run_riker.smk — these survive
        # Snakemake's auto-cleanup on failure so we keep the JVM stderr.
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cwm/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cwm/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "wgs-t1|wgs-bundle",
    params:
        xmx = config["picard_xmx_gb"],
        # Harmonize to riker/mosdepth in the fair three-way (wgs-t1): count
        # orphan/unpaired reads. In wgs-bundle, CWM is a Picard-faithful
        # supplement to riker multi (default wgs), so leave it off there.
        count_unpaired = lambda w: "true" if w.profile == "wgs-t1" else "false",
    threads: 1
    resources: bench=100
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        cmd=(picard CollectWgsMetrics
             "I={input.bam}"
             "R={input.ref}"
             "COUNT_UNPAIRED={params.count_unpaired}"
             "O=$outdir/picard.wgs_metrics.txt")
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        _JAVA_OPTIONS="-Xmx{params.xmx}g" command time -v -o {output.time:q} \
            "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """


# ---- CollectHsMetrics ----------------------------------------------------

rule run_picard_chsm:
    """CollectHsMetrics. Headline tool for hybcap-only AND supplementary
    in hybcap-bundle. The kit's interval_list is used directly for both
    BAIT_INTERVALS and TARGET_INTERVALS — no BedToIntervalList step."""
    input:
        bam       = f"{STAGE_DIR}/{{sample}}/input.bam",
        bai       = f"{STAGE_DIR}/{{sample}}/input.bam.bai",
        intervals = lambda w: kit_intervals_for_sample(w.sample),
        ref       = lambda w: ref_for_sample(w.sample),
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-chsm/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-chsm/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-chsm/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "hybcap-only|hybcap-bundle",
    params:
        xmx = config["picard_xmx_gb"],
    threads: 1
    resources: bench=100
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"

        cmd=(picard CollectHsMetrics
             "I={input.bam}"
             "R={input.ref}"
             "BAIT_INTERVALS={input.intervals}"
             "TARGET_INTERVALS={input.intervals}"
             "O=$outdir/picard.hs_metrics.txt")
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        _JAVA_OPTIONS="-Xmx{params.xmx}g" command time -v -o {output.time:q} \
            "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """


# ---- CollectMultipleMetrics (bundle profiles only) ----------------------

# Programs to bundle for each profile. wgs-bundle adds CollectGcBiasMetrics
# on top of the default-5 so it matches riker multi's wgs-bundle (minus
# the WGS coverage tool, which has its own CWM run as a supplementary).
_CMM_PROGRAMS = {
    "wgs-bundle": [
        "CollectAlignmentSummaryMetrics",
        "CollectInsertSizeMetrics",
        "QualityScoreDistribution",
        "MeanQualityByCycle",
        "CollectBaseDistributionByCycle",
        "CollectGcBiasMetrics",
    ],
    "hybcap-bundle": [
        "CollectAlignmentSummaryMetrics",
        "CollectInsertSizeMetrics",
        "QualityScoreDistribution",
        "MeanQualityByCycle",
        "CollectBaseDistributionByCycle",
    ],
}


rule run_picard_cmm:
    """CollectMultipleMetrics. Single JVM that produces several metrics in
    one BAM pass — Picard's analog to riker `multi`. We enable async I/O
    so it's at least somewhat fair vs riker --threads 4."""
    input:
        bam = f"{STAGE_DIR}/{{sample}}/input.bam",
        bai = f"{STAGE_DIR}/{{sample}}/input.bam.bai",
        ref = lambda w: ref_for_sample(w.sample),
        fai = lambda w: ref_for_sample(w.sample) + ".fai",
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cmm/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cmm/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-cmm/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "wgs-bundle|hybcap-bundle",
    params:
        xmx       = config["picard_xmx_gb"],
        # Pre-build the PROGRAM= flag list as a single shell-quoted string.
        programs  = lambda w: " ".join(f'"PROGRAM={p}"' for p in _CMM_PROGRAMS[w.profile]),
    threads: 1
    resources: bench=100
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        cmd=(picard CollectMultipleMetrics
             "I={input.bam}"
             "R={input.ref}"
             "O=$outdir/picard.cmm")
        # Append PROGRAM= entries (already shell-quoted by params.programs).
        eval "cmd+=( {params.programs} )"
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        # Async-IO is a JVM system property in picard 3.x, not a CLI flag.
        # Enable it so picard's BAM read/write threads at least somewhat
        # match riker --threads 4 in the bundle profiles.
        _JAVA_OPTIONS="-Xmx{params.xmx}g -Dsamjdk.async_io_read_samtools=true -Dsamjdk.async_io_write_samtools=true" \
            command time -v -o {output.time:q} \
                "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """


# ---- CollectRnaSeqMetrics (rna profiles) --------------------------------

rule run_picard_crsm:
    """CollectRnaSeqMetrics — the direct analog of riker rna's base-level
    RNA metrics. Single-threaded, so it only appears in the rna-t1 profile.
    REF_FLAT (derived from the same GTF riker/rustqc use), RIBOSOMAL_INTERVALS
    (per-sample rRNA loci), and the sample's STRAND_SPECIFICITY."""
    input:
        bam       = f"{STAGE_DIR}/{{sample}}/input.bam",
        bai       = f"{STAGE_DIR}/{{sample}}/input.bam.bai",
        refflat   = lambda w: refflat_for_sample(w.sample),
        ribosomal = lambda w: ribosomal_intervals_for_sample(w.sample),
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-crsm/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-crsm/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/picard-crsm/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "rna-t1",
    params:
        xmx    = config["picard_xmx_gb"],
        strand = lambda w: picard_strand_for_sample(w.sample),
    threads: 1
    resources: bench=100
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        cmd=(picard CollectRnaSeqMetrics
             "I={input.bam}"
             "REF_FLAT={input.refflat}"
             "RIBOSOMAL_INTERVALS={input.ribosomal}"
             "STRAND_SPECIFICITY={params.strand}"
             "O=$outdir/picard.rna_metrics.txt")
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        _JAVA_OPTIONS="-Xmx{params.xmx}g" command time -v -o {output.time:q} \
            "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
