# mosdepth benchmark rule. Mosdepth is a single-purpose coverage tool; it
# only makes sense as a comparator for the wgs-only and hybcap-only
# profiles. Bundle profiles don't include mosdepth (no analog of riker
# multi).
#
# We always run with `--no-per-base -x` to skip the heavy per-base output
# riker doesn't produce either; this matches "what users actually run for
# QC". For hybcap, `--by <target.bed>` constrains to capture regions.
#
# Threading: deliberately omit `-t` so mosdepth uses its default
# (`-t 0` = no extra decompression threads, htslib runs single-threaded).
# Riker's standalone wgs/hybcap subcommands are also single-threaded, so
# the comparison is fair only without mosdepth's extra decompression
# thread. We measured a 124-149% CPU footprint with `-t 1` because
# htslib spins up an extra block-decompress thread on top of the main
# thread; on 30x WGS that buys mosdepth ~30-40% wall-time over riker
# even though riker uses ~10% LESS total CPU time.

def _mosdepth_inputs(wildcards):
    inputs = {
        "bam": f"{STAGE_DIR}/{wildcards.sample}/input.bam",
        "bai": f"{STAGE_DIR}/{wildcards.sample}/input.bam.bai",
    }
    if wildcards.profile == "hybcap-only":
        # mosdepth needs BED, not interval_list, so use the derived BED.
        inputs["target_bed"] = kit_bed_for_sample(wildcards.sample)
    return inputs


rule run_mosdepth:
    input: unpack(_mosdepth_inputs)
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/mosdepth/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/mosdepth/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/mosdepth/rep{{rep}}/tool.log",
    wildcard_constraints:
        # mosdepth doesn't run for bundle profiles.
        profile = "wgs-only|hybcap-only",
    threads: 1
    resources: bench=100
    params:
        target_bed = lambda w: kit_bed_for_sample(w.sample) if w.profile == "hybcap-only" else "",
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        if [[ "{wildcards.profile}" == "hybcap-only" ]]; then
            cmd=(mosdepth --by {params.target_bed:q} --no-per-base \
                          "$outdir/mosdepth" {input.bam})
        else
            cmd=(mosdepth -x --no-per-base \
                          "$outdir/mosdepth" {input.bam})
        fi
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        command time -v -o {output.time:q} "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
