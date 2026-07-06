# RustQC (Seqera) benchmark rule — single-pass RNA-seq QC, the modern
# Rust peer to `riker rna`. We run a scope-matched SUBSET (see
# config/rustqc.rna.subset.yaml) so the comparison is fair: the config
# disables the modules riker doesn't compute (dupRadar, preseq,
# read_duplication, junction_saturation, the samtools passthroughs) and
# keeps the ones it does (read distribution, strand, junction annotation,
# inner distance, TIN, gene-body coverage via qualimap, gene assignment
# via featureCounts). `--skip-dup-check` lets it run on the ENCODE BAMs,
# which aren't duplicate-marked (dupRadar is off, so nothing needs it).
#
# RustQC is multithreaded, so it runs across the rna-t1/-t2/-t4 thread
# sweep at the profile's budget, alongside riker. Installed from Seqera's
# own bioconda distribution (pixi), so this is riker-as-shipped vs
# RustQC-as-shipped.

rule run_rustqc:
    input:
        bam    = f"{STAGE_DIR}/{{sample}}/input.bam",
        bai    = f"{STAGE_DIR}/{{sample}}/input.bam.bai",
        gtf    = lambda w: gene_model_gtf_for_sample(w.sample),
        config = str(RUSTQC_RNA_CONFIG),
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/rustqc/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/rustqc/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/rustqc/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "rna-t1|rna-t2|rna-t4",
    threads: lambda w: thread_count(w.profile)
    resources: bench=100
    params:
        stranded = lambda w: strand_for_sample(w.sample),
        layout   = lambda w: SAMPLES.loc[w.sample, "layout"],
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        # --paired only for PE. -Q 0 matches riker (--min-mapq 0) and Picard
        # (no MAPQ filter) — RustQC otherwise defaults to MAPQ>=30, which
        # would process fewer reads and be unfair. -s pins strandedness to
        # match riker --strand / Picard STRAND_SPECIFICITY. The scope match
        # (which modules run) lives entirely in the -c config file.
        paired=()
        if [[ "{params.layout}" == "PE" ]]; then paired=(--paired); fi
        cmd=(rustqc rna {input.bam}
             --gtf {input.gtf}
             --config {input.config}
             --skip-dup-check
             -Q 0
             -s {params.stranded}
             "${{paired[@]}}"
             --threads {threads}
             --outdir "$outdir/rustqc")
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        command time -v -o {output.time:q} "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
