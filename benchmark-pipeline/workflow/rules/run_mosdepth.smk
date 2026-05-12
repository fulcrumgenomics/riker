# mosdepth benchmark rule. Mosdepth is a single-purpose coverage tool;
# it only makes sense as a comparator for the wgs-only / wgs-mosdepth /
# hybcap-only profiles. Bundle profiles don't include mosdepth (no
# analog of `riker multi`).
#
# Two tool tokens, mapping to two invocations:
#
#   * `mosdepth`           runs with `-x` in the wgs profiles, without
#                          `-x` in hybcap-only. `-x` is mosdepth's
#                          documented "fast mode": skip the per-record
#                          CIGAR walk (each fragment counts as one
#                          contiguous block) AND skip mate-overlap
#                          deduplication. That is faster but does NOT
#                          match riker's accuracy.
#   * `mosdepth-accurate`  always runs WITHOUT `-x`, so per-record
#                          CIGAR walking and mate-overlap correction
#                          are enabled. This is the apples-to-apples
#                          comparator against riker's wgs.
#
# Both invocations also pass `--no-per-base`, since neither riker `wgs`
# nor Picard `CollectWgsMetrics` writes a per-base output file. The
# hybcap-only branch passes `--by <target.bed>` to constrain mosdepth
# to the capture regions.
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


def _mosdepth_use_fast_mode(profile: str, tool: str) -> bool:
    """`-x` (fast mode) selection:
    - `mosdepth-accurate` token: never (the whole point of the token).
    - hybcap-only profile:       never (accuracy parity with Picard CHsM).
    - everything else (wgs-only, wgs-mosdepth + mosdepth token): yes.
    """
    if tool == "mosdepth-accurate":
        return False
    if profile == "hybcap-only":
        return False
    return True


rule run_mosdepth:
    input: unpack(_mosdepth_inputs)
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/tool.log",
    wildcard_constraints:
        # mosdepth doesn't run for bundle profiles.
        profile = "wgs-only|wgs-mosdepth|hybcap-only",
        tool    = "mosdepth|mosdepth-accurate",
    threads: 1
    resources: bench=100
    params:
        target_bed = lambda w: kit_bed_for_sample(w.sample) if w.profile == "hybcap-only" else "",
        fast_mode  = lambda w: _mosdepth_use_fast_mode(w.profile, w.tool),
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        if [[ "{wildcards.profile}" == "hybcap-only" ]]; then
            cmd=(mosdepth --by {params.target_bed:q} --no-per-base \
                          "$outdir/mosdepth" {input.bam})
        else
            if [[ "{params.fast_mode}" == "True" ]]; then
                cmd=(mosdepth -x --no-per-base \
                              "$outdir/mosdepth" {input.bam})
            else
                cmd=(mosdepth --no-per-base \
                              "$outdir/mosdepth" {input.bam})
            fi
        fi
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        command time -v -o {output.time:q} "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
