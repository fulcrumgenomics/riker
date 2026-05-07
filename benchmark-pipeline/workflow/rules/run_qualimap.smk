# Qualimap bamqc benchmark rule. One rule covers all four profiles —
# threading: single-tool profiles run -nt 1, bundle profiles run -nt 4
# to roughly mirror riker --threads 4.
#
# Qualimap heap behavior on 30x WGS is documented as memory-hungry;
# default qualimap_xmx_gb=16 with attempts:3 doubling on retry. If
# qualimap fails repeatedly even at 64 GB, drop it from COMPARATORS in
# the Snakefile and document.
#
# We do NOT parse qualimap's output in v1; we only need the timing.

def _qualimap_inputs(wildcards):
    """For hybcap profiles, qualimap reads the kit BED via -gff. We use
    the BED form (qualimap's -gff accepts BED despite the flag name)."""
    inputs = {
        "bam": f"{STAGE_DIR}/{wildcards.sample}/input.bam",
        "bai": f"{STAGE_DIR}/{wildcards.sample}/input.bam.bai",
    }
    if wildcards.profile.startswith("hybcap-"):
        inputs["target_bed"] = kit_bed_for_sample(wildcards.sample)
    return inputs


rule run_qualimap:
    input: unpack(_qualimap_inputs)
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/qualimap/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/qualimap/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/qualimap/rep{{rep}}/tool.log",
    threads: lambda w: thread_count(w.profile)
    # On Java OOM, retry once with doubled heap (capped by
    # qualimap_xmx_gb_max). If the doubled heap still OOMs the rule fails
    # and the pipeline aborts — qualimap is the most OOM-prone comparator
    # on 30x WGS and getting past that limit is best-effort.
    retries: 1
    resources:
        bench = 100,
        # `attempt` is 1 on first try, 2 on retry. Multiply the configured
        # initial heap by `attempt` and cap at `qualimap_xmx_gb_max`.
        qualimap_xmx = lambda w, attempt: min(
            int(config["qualimap_xmx_gb"]) * attempt,
            int(config.get("qualimap_xmx_gb_max", int(config["qualimap_xmx_gb"]) * 2)),
        ),
    params:
        # Always a string — empty for wgs profiles, so the shell `if`
        # branch never references it.
        target_bed = lambda w: kit_bed_for_sample(w.sample) if w.profile.startswith("hybcap-") else "",
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir/qmreport"
        cmd=(qualimap bamqc
             -bam {input.bam:q}
             -outdir "$outdir/qmreport"
             -outformat HTML
             -nt {threads}
             "--java-mem-size={resources.qualimap_xmx}G")
        if [[ "{wildcards.profile}" == hybcap-* ]]; then
            cmd+=(-gff {params.target_bed:q})
        fi
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # qualimap's wrapper would otherwise try to use an X display;
        # force headless explicitly.
        JAVA_OPTS="-Djava.awt.headless=true" \
            command time -v -o {output.time:q} "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
