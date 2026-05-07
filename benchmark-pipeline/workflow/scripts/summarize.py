#!/usr/bin/env python3
"""Collapse replicates in bench.tsv to per-cell medians + min/max bars.

Group by (sample, profile, tool); compute median wall_s, user_s, sys_s,
max_rss_gb, throughput_mb_per_s; carry through metadata columns
unchanged (they're constant within a group)."""

import argparse
from pathlib import Path

import pandas as pd


GROUP_KEYS = ["sample", "input_type", "profile", "tool", "tool_family",
              "picard_bundle_role", "threads", "effective_threads",
              "input_format"]

NUMERIC_COLS = ["wall_s", "user_s", "sys_s", "max_rss_kb", "max_rss_gb",
                "cpu_percent", "throughput_mb_per_s", "sample_size_gb"]

CARRY_COLS = ["tool_version", "approx_bam_gb",
              "host_hostname", "host_arch", "host_cpu_model",
              "host_cpu_count_logical", "host_mem_bytes",
              "host_os", "host_os_release", "host_aws_instance_type"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.inp, sep="\t")
    if len(df) == 0:
        Path(args.out).write_text("")
        return

    # Per-cell aggregation: median + min + max for each numeric column.
    agg_dict: dict[str, list[str]] = {c: ["median", "min", "max"] for c in NUMERIC_COLS if c in df.columns}
    grouped = df.groupby(GROUP_KEYS, dropna=False).agg(agg_dict)
    grouped.columns = [f"{col}_{stat}" for col, stat in grouped.columns]
    grouped = grouped.reset_index()

    # Carry through metadata that's constant within a group.
    for c in CARRY_COLS:
        if c in df.columns:
            first = df.groupby(GROUP_KEYS, dropna=False)[c].first().reset_index()
            grouped = grouped.merge(first, on=GROUP_KEYS, how="left")

    # Replicate count per group (sanity check column).
    counts = df.groupby(GROUP_KEYS, dropna=False).size().reset_index(name="rep_count")
    grouped = grouped.merge(counts, on=GROUP_KEYS, how="left")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    grouped.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
