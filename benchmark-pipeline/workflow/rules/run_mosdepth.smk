# mosdepth benchmark rule. Comparator for the wgs-t{1..4} (fair, harmonized
# three-way WGS) and hybcap-only profiles. Bundle profiles have no mosdepth
# analog (no `riker multi` equivalent).
#
# WGS (wgs-t{N}): mosdepth runs in its ACCURATE mode (no -x fast-mode, no -a
# fragment-mode), harmonized to riker/Picard's read selection so the three
# tools count the same reads:
#   -Q 20          MAPQ floor (matches riker --min-mapq default / Picard default)
#   -F 3332        exclude unmapped+secondary+supplementary+dup — matches
#                  riker's whole-read filters; orphan/unpaired reads are still
#                  counted (as riker --include-unpaired-reads / Picard
#                  COUNT_UNPAIRED=true do), and qcfail is counted (as riker does)
#   --no-per-base  neither riker wgs nor Picard CollectWgsMetrics writes per-base
# mosdepth cannot base-quality filter, so riker/Picard's BQ>=20 is the single
# accepted asymmetry (an intentional riker/Picard feature).
#
# Threads: a TOTAL-thread sweep. mosdepth's -t counts EXTRA BAM-decompression
# threads on top of the main thread, so N total threads = `-t (N-1)`.
#
# hybcap-only: `--by <targets.bed>` constrains coverage to the capture
# regions; single-threaded (matches riker hybcap), also accurate (no -x).

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
        time     = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/time.txt",
    log:
        cmdline  = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/cmdline.txt",
        tool_log = f"{RESULTS_DIR}/run/{{sample}}/{{profile}}/{{tool}}/rep{{rep}}/tool.log",
    wildcard_constraints:
        profile = "wgs-t1|wgs-t2|wgs-t3|wgs-t4|hybcap-only",
        tool    = "mosdepth",
    threads: lambda w: thread_count(w.profile)
    resources: bench=100
    params:
        target_bed = lambda w: kit_bed_for_sample(w.sample) if w.profile == "hybcap-only" else "",
        # mosdepth -t is EXTRA decompression threads, so total N => -t (N-1).
        decomp_threads = lambda w: max(thread_count(w.profile) - 1, 0),
    shell:
        r"""
        set -euo pipefail
        outdir="$(dirname {output.time:q})"
        mkdir -p "$outdir"
        if [[ "{wildcards.profile}" == "hybcap-only" ]]; then
            cmd=(mosdepth --by {params.target_bed:q} --no-per-base \
                          "$outdir/mosdepth" {input.bam})
        else
            # wgs-t{{N}}: accurate + harmonized; total N threads = -t (N-1).
            cmd=(mosdepth -t {params.decomp_threads} -Q 20 -F 3332 --no-per-base \
                          "$outdir/mosdepth" {input.bam})
        fi
        printf '%s ' "${{cmd[@]}}" > {log.cmdline:q}; echo >> {log.cmdline:q}
        # Cold cache: drop the page cache so the tool reads its inputs cold
        # from disk — deterministic + size-uniform across the coverage ladder
        # (no BAM fits-in-RAM advantage). Fails loudly without sudo/root.
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
        command time -v -o {output.time:q} "${{cmd[@]}}" > {log.tool_log:q} 2>&1
        """
