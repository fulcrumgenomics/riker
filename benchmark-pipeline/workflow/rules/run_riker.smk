# Riker benchmark rule. One rule covers all four profiles by reading
# PROFILES[wildcards.profile] for the argv suffix and thread count.
#
# Output schema:
#   results/run/{sample}/{profile}/riker/rep{rep}/
#     time.txt     /usr/bin/time -v output (parsed by parse_gnu_time.py)
#     cmdline.txt  the exact command line that was run
#     tool.log     stdout + stderr from riker
#     out.*        riker's own metric TSVs / PDFs (not parsed in v1)

import shlex


def _riker_argv(wildcards) -> list[str]:
    """Build the riker argv. Returned as a list so we can shlex.join it
    upstream of the shell block.

    Standalone subcommands (wgs, hybcap, ...) are single-threaded and
    take their kit BEDs as plain `--baits`/`--targets`. Only `riker
    multi` has `--threads` and uses the `--hybcap::` prefix for kit
    BEDs."""
    profile = PROFILES[wildcards.profile]
    sample = wildcards.sample
    rep = wildcards.rep
    is_multi = profile["riker_args"][0] == "multi"

    bam = f"{STAGE_DIR}/{sample}/input.bam"
    out_prefix = f"{RESULTS_DIR}/run/{sample}/{wildcards.profile}/riker/rep{rep}/out"

    argv = [str(RIKER_BIN), *profile["riker_args"]]
    argv += ["--input", bam]
    argv += ["--output", out_prefix]
    if is_multi:
        argv += ["--threads", str(thread_count(wildcards.profile))]

    needs_ref = profile["needs_ref"] or (is_multi and (
        "wgs" in profile["riker_args"] or "gcbias" in profile["riker_args"]
    ))
    if needs_ref:
        argv += ["--reference", ref_for_sample(sample)]

    if is_multi and "hybcap" in profile["riker_args"]:
        argv += [
            "--hybcap::baits",   bait_bed_for_sample(sample),
            "--hybcap::targets", target_bed_for_sample(sample),
        ]
    elif not is_multi and profile["riker_args"][0] == "hybcap":
        argv += [
            "--baits",   bait_bed_for_sample(sample),
            "--targets", target_bed_for_sample(sample),
        ]
    return argv


def _riker_inputs(wildcards):
    """Files that must exist before riker runs. Reference / BEDs are
    optional depending on profile."""
    sample = wildcards.sample
    profile = PROFILES[wildcards.profile]
    inputs = {
        "bam": f"{STAGE_DIR}/{sample}/input.bam",
        "bai": f"{STAGE_DIR}/{sample}/input.bam.bai",
    }
    needs_ref = profile["needs_ref"] or "wgs" in profile["riker_args"] or "gcbias" in profile["riker_args"]
    if needs_ref:
        inputs["ref"] = ref_for_sample(sample)
        inputs["fai"] = ref_for_sample(sample) + ".fai"
    if "hybcap" in profile["riker_args"]:
        inputs["baits"]   = bait_bed_for_sample(sample)
        inputs["targets"] = target_bed_for_sample(sample)
    return inputs


rule run_riker:
    input: unpack(_riker_inputs)
    output:
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/riker/rep{{rep}}/time.txt",
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/riker/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/riker/rep{{rep}}/tool.log",
    params:
        # Pre-shell-quote the argv so paths with spaces are safe.
        argv_str = lambda w: shlex.join(_riker_argv(w)),
    threads: lambda w: thread_count(w.profile)
    resources: bench=100
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.time:q})"
        printf '%s\n' {params.argv_str:q} > {output.cmdline:q}
        # `command time` bypasses the bash builtin so we get the conda-forge
        # GNU time binary (which supports -v / -o).
        command time -v -o {output.time:q} {params.argv_str} > {output.tool_log:q} 2>&1
        """
