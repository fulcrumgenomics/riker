#!/usr/bin/env python3
"""Walk results/run/**/time.txt, parse GNU time output, join with sample
metadata + host info + tool versions, emit a single wide TSV.

Output columns are documented in benchmark-pipeline/README.md §"Output schema"
and in the plan file. One row per (sample, profile, tool, rep)."""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# parse_gnu_time sits next to this script.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from parse_gnu_time import parse as parse_time


# results/run/{sample}/{profile}/{tool}/rep{rep}/time.txt
PATH_RE = re.compile(
    r"results?/run/"
    r"(?P<sample>[^/]+)/"
    r"(?P<profile>[^/]+)/"
    r"(?P<tool>[^/]+)/"
    r"rep(?P<rep>\d+)/"
    r"time\.txt$"
)


# Map profile + tool token to the picard_bundle_role column. None for
# non-picard tools.
def picard_bundle_role(profile: str, tool: str) -> str | None:
    if not tool.startswith("picard-"):
        return None
    if profile.endswith("-only"):
        return "main"
    # bundle profiles
    if tool == "picard-cmm":
        return "base"
    return "supp"


def thread_count(profile: str, tool: str) -> int:
    """Snakefile sets threads via the rule's threads: directive; we mirror
    that here so bench.tsv accurately reflects what was requested."""
    if profile.endswith("-bundle") and tool in ("riker", "qualimap"):
        return 4
    return 1


def effective_threads(tool: str, requested: int) -> int:
    """Picard is single-threaded per tool — even with USE_ASYNC_IO_*=true
    the compute is mostly serial. Mosdepth uses 1 here; qualimap uses
    requested; riker uses requested."""
    if tool.startswith("picard-"):
        return 1
    return requested


def tool_family(tool: str) -> str:
    if tool == "riker": return "riker"
    if tool.startswith("picard-"): return "picard"
    return tool


def safe_version(cmd: list[str]) -> str:
    """Run a `--version` command and return the first non-empty line.
    Best-effort — empty string on failure."""
    try:
        out = subprocess.check_output(cmd, text=True, timeout=15,
                                      stderr=subprocess.STDOUT)
        for line in out.splitlines():
            line = line.strip()
            if line:
                return line
    except Exception:
        return ""
    return ""


def collect_versions(riker_bin: str) -> dict[str, str]:
    """Capture tool versions once. The riker binary path comes in as a
    config value because it lives outside the pixi env."""
    return {
        "riker":    safe_version([riker_bin, "--version"]),
        "picard":   safe_version(["picard", "MarkDuplicates", "--version"]),
        "mosdepth": safe_version(["mosdepth", "--version"]),
        # qualimap doesn't have a --version. `qualimap --help` greets with
        # "Java memory size is set to ..." which is useless; grab the
        # second non-empty line which actually carries the version.
        "qualimap": _qualimap_version(),
    }


def _qualimap_version() -> str:
    try:
        out = subprocess.check_output(["qualimap", "--help"], text=True,
                                      timeout=15, stderr=subprocess.STDOUT)
        lines = [l.strip() for l in out.splitlines() if l.strip()]
        # Look for a line containing "QualiMap" with a version-y word.
        for line in lines:
            if "qualimap" in line.lower() and any(c.isdigit() for c in line):
                return line
        # Fallback: any line that has "v" + digits, common for version strings.
        for line in lines:
            if any(seg.startswith("v") and seg[1:2].isdigit() for seg in line.split()):
                return line
        return lines[0] if lines else "unknown"
    except Exception as e:
        return f"error: {e}"


def bam_size_gb(stage_dir: Path, sample: str) -> float:
    bam = stage_dir / sample / "input.bam"
    if not bam.exists():
        return 0.0
    return bam.stat().st_size / 1e9


def load_samples(wgs_path: Path, hyb_path: Path) -> pd.DataFrame:
    frames = []
    if wgs_path.exists():
        wgs = pd.read_csv(wgs_path, sep="\t", dtype=str).fillna("-")
        if len(wgs) > 0:
            wgs["input_type"] = "wgs"
            frames.append(wgs)
    if hyb_path.exists():
        hyb = pd.read_csv(hyb_path, sep="\t", dtype=str).fillna("-")
        if len(hyb) > 0:
            hyb["input_type"] = "hybcap"
            frames.append(hyb)
    if not frames:
        return pd.DataFrame(columns=["name", "input_type"])
    df = pd.concat(frames, ignore_index=True)
    return df.set_index("name", drop=False)


def parse_path(path: Path) -> dict:
    m = PATH_RE.search(str(path))
    if not m:
        raise ValueError(f"path doesn't match expected schema: {path}")
    d = m.groupdict()
    d["rep"] = int(d["rep"])
    return d


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--host", required=True, help="host.json")
    p.add_argument("--root", required=True, help="results/run dir")
    p.add_argument("--samples-wgs", required=True)
    p.add_argument("--samples-hybcap", required=True)
    p.add_argument("--riker-bin", required=True)
    p.add_argument("--stage-dir", required=True,
                   help="Where stage/<sample>/input.bam lives — passed through "
                        "from the Snakefile config so this script doesn't have "
                        "to guess when stage_dir is overridden via --config")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with open(args.host) as fh:
        host = json.load(fh)

    samples = load_samples(Path(args.samples_wgs), Path(args.samples_hybcap))
    versions = collect_versions(args.riker_bin)

    run_root = Path(args.root)
    stage_dir = Path(args.stage_dir)

    rows = []
    for time_txt in sorted(run_root.glob("*/*/*/rep*/time.txt")):
        meta = parse_path(time_txt)
        timings = parse_time(time_txt)
        sample = meta["sample"]
        profile = meta["profile"]
        tool = meta["tool"]

        tf = tool_family(tool)
        family_version = versions.get(tf, "")

        sample_row = samples.loc[sample] if sample in samples.index else None
        input_type = sample_row["input_type"] if sample_row is not None else "-"
        approx_bam_gb = sample_row.get("approx_bam_gb", "-") if sample_row is not None else "-"

        size_gb = bam_size_gb(stage_dir, sample)
        wall = float(timings.get("wall_s") or 0.0)
        threads = thread_count(profile, tool)
        max_rss_kb = int(timings.get("max_rss_kb") or 0)

        row = {
            "sample": sample,
            "input_type": input_type,
            "profile": profile,
            "tool": tool,
            "tool_family": tf,
            "tool_version": family_version,
            "picard_bundle_role": picard_bundle_role(profile, tool) or "",
            "threads": threads,
            "effective_threads": effective_threads(tool, threads),
            "rep": meta["rep"],
            "input_format": "bam",
            "sample_size_gb": round(size_gb, 3),
            "approx_bam_gb": approx_bam_gb,
            "wall_s": wall,
            "user_s": float(timings.get("user_s") or 0.0),
            "sys_s":  float(timings.get("sys_s")  or 0.0),
            "cpu_percent": int(timings.get("cpu_percent") or 0),
            "max_rss_kb": max_rss_kb,
            "max_rss_gb": round(max_rss_kb / (1024 * 1024), 3),
            "throughput_mb_per_s": round((size_gb * 1024) / wall, 2) if wall > 0 else 0.0,
            "exit_status": int(timings.get("exit_status") or 0),
            "host_hostname": host.get("hostname"),
            "host_arch": host.get("arch"),
            "host_cpu_model": host.get("cpu_model"),
            "host_cpu_count_logical": host.get("cpu_count_logical"),
            "host_mem_bytes": host.get("total_mem_bytes"),
            "host_os": host.get("os"),
            "host_os_release": host.get("os_release"),
            "host_aws_instance_type": host.get("aws_instance_type"),
        }
        rows.append(row)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        print(f"WARNING: no time.txt files under {run_root}", file=sys.stderr)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
