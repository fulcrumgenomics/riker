# Aggregation: collect every per-run time.txt and emit results/bench.tsv.
# Then collapse replicates into bench_summary.tsv. Then plots.
#
# All rules here reserve bench=1 — they can run during DAG fill but
# never DURING a timed rule.

rule host_info:
    output: f"{RESULTS_DIR}/host.json"
    resources: bench=1
    shell:
        "python {WORKDIR}/workflow/scripts/host_info.py > {output}"


rule aggregate_bench:
    """Walk results/run/**/time.txt and emit a wide TSV with one row per
    (sample, profile, tool, rep)."""
    input:
        time_files = all_bench_targets(),
        host       = f"{RESULTS_DIR}/host.json",
    output: f"{RESULTS_DIR}/bench.tsv"
    params:
        samples_wgs    = str(SAMPLES_WGS_TSV),
        samples_hybcap = str(SAMPLES_HYBCAP_TSV),
        samples_rna    = str(SAMPLES_RNA_TSV) if SAMPLES_RNA_TSV else "",
        run_root       = f"{RESULTS_DIR}/run",
        riker_bin      = str(RIKER_BIN),
        stage_dir      = str(STAGE_DIR),
    resources: bench=1
    shell:
        r"""
        # --samples-rna is optional; bind to a shell var first because
        # Snakemake's :q renders an empty string as nothing (not ''), which
        # would let argparse swallow the next flag as its value.
        samples_rna={params.samples_rna:q}
        rna_arg=()
        [[ -n "$samples_rna" ]] && rna_arg=(--samples-rna "$samples_rna")
        python {WORKDIR}/workflow/scripts/merge_results.py \
            --host {input.host:q} \
            --root {params.run_root:q} \
            --samples-wgs {params.samples_wgs:q} \
            --samples-hybcap {params.samples_hybcap:q} \
            "${{rna_arg[@]}}" \
            --riker-bin {params.riker_bin:q} \
            --stage-dir {params.stage_dir:q} \
            --out {output:q}
        """


rule summarize_bench:
    """Collapse replicates: median wall/user/sys/RSS plus min/max bars."""
    input: f"{RESULTS_DIR}/bench.tsv"
    output: f"{RESULTS_DIR}/bench_summary.tsv"
    resources: bench=1
    shell:
        "python {WORKDIR}/workflow/scripts/summarize.py --in {input:q} --out {output:q}"


rule plots:
    input:
        bench   = f"{RESULTS_DIR}/bench.tsv",
        summary = f"{RESULTS_DIR}/bench_summary.tsv",
    output:
        leaderboard = f"{RESULTS_DIR}/plots/01_wall_time_leaderboard.pdf",
        throughput  = f"{RESULTS_DIR}/plots/02_throughput_vs_size.pdf",
        rss         = f"{RESULTS_DIR}/plots/03_max_rss_vs_size.pdf",
        bundle      = f"{RESULTS_DIR}/plots/04_riker_bundle_vs_picard.pdf",
        head2head   = f"{RESULTS_DIR}/plots/05_riker_vs_picard_cmm.pdf",
        ratio       = f"{RESULTS_DIR}/plots/06_per_comparator_ratio.pdf",
    resources: bench=1
    shell:
        r"""
        mkdir -p "$(dirname {output.leaderboard:q})"
        Rscript {WORKDIR}/workflow/scripts/plot.R \
            {input.bench:q} {input.summary:q} {RESULTS_DIR}/plots
        """
