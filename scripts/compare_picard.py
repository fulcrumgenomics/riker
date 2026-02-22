#!/usr/bin/env python3
"""Compare riker error output against Picard CollectSamErrorMetrics.

Usage:
    python3 scripts/compare_picard.py [--riker-prefix /tmp/riker_r1] [--picard-prefix /tmp/picard_r1]

Expects both tools to have already been run with the given output prefixes.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# ── Stratifier mapping: riker_name -> picard_file_suffix ──────────────────────
RIKER_TO_PICARD_SUFFIX = {
    "all": "all",
    "bq": "base_quality",
    "mapq": "mapping_quality",
    "cycle": "cycle",
    "read_num": "read_ordinality",
    "strand": "read_direction",
    "pair_orientation": "pair_orientation",
    "isize": "insert_length",
    "gc": "gc",
    "read_base": "read_base",
    "ref_base": "ref_base",
    "hp_len": "homopolymer_length",
    "pre_dinuc": "pre_dinuc",
    "post_dinuc": "post_dinuc",
    "context_3bp": "one_base_padded_context",
    "nm": "mismatches_in_read",
    "indel_len": "indel_length",
}

# Error type mapping
ERROR_TYPES = {
    "mismatch": {"riker_suffix": "error-mismatch", "picard_suffix": "error"},
    "overlap": {"riker_suffix": "error-overlap", "picard_suffix": "overlapping_error"},
    "indel": {"riker_suffix": "error-indel", "picard_suffix": "indel_error"},
}

# Columns to compare per error type
COMPARE_COLUMNS = {
    "mismatch": {
        "riker": ["total_bases", "error_bases"],
        "picard": ["TOTAL_BASES", "ERROR_BASES"],
    },
    "overlap": {
        "riker": ["overlapping_read_bases", "bases_matching_mate_but_not_ref", "bases_mismatching_ref_and_mate", "bases_in_three_way_disagreement"],
        "picard": ["NUM_BASES_WITH_OVERLAPPING_READS", "NUM_DISAGREES_WITH_REFERENCE_ONLY", "NUM_DISAGREES_WITH_REF_AND_MATE", "NUM_THREE_WAYS_DISAGREEMENT"],
    },
    "indel": {
        "riker": ["total_bases", "num_insertions", "num_inserted_bases", "num_deletions", "num_deleted_bases"],
        "picard": ["TOTAL_BASES", "NUM_INSERTIONS", "NUM_INSERTED_BASES", "NUM_DELETIONS", "NUM_DELETED_BASES"],
    },
}


def parse_picard(path: Path) -> list[dict]:
    """Parse a Picard metrics file, returning list of dicts keyed by column name."""
    rows = []
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith("##"):
                continue
            if not line.strip():
                if header is not None:
                    break  # blank after data = end of metrics section
                continue
            if header is None:
                header = line.split("\t")
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            rows.append(row)
    return rows


def parse_riker(path: Path, stratifier: str) -> list[dict]:
    """Parse a riker TSV file, filtering to rows matching the given stratifier."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["stratifier"] == stratifier:
                rows.append(row)
    return rows


def compare_one(strat: str, error_key: str, riker_prefix: str, picard_prefix: str, diff_dir: Path) -> str:
    """Compare one stratifier x error_type. Returns 'PASS', 'FAIL', or 'SKIP'."""
    riker_file = Path(f"{riker_prefix}.{ERROR_TYPES[error_key]['riker_suffix']}.txt")
    picard_suffix = RIKER_TO_PICARD_SUFFIX[strat]
    picard_file = Path(f"{picard_prefix}.{ERROR_TYPES[error_key]['picard_suffix']}_by_{picard_suffix}")

    diff_file = diff_dir / f"{strat}_{error_key}.diff"

    if not picard_file.exists():
        diff_file.write_text(f"SKIP: Picard file not found: {picard_file}\n")
        return "SKIP"

    if not riker_file.exists():
        diff_file.write_text(f"SKIP: Riker file not found: {riker_file}\n")
        return "SKIP"

    riker_rows = parse_riker(riker_file, strat)
    picard_rows = parse_picard(picard_file)

    # Build maps keyed by covariate
    r_cols = COMPARE_COLUMNS[error_key]["riker"]
    p_cols = COMPARE_COLUMNS[error_key]["picard"]
    col_labels = [rc for rc in r_cols]  # use riker names as labels

    riker_map = {}
    for row in riker_rows:
        cov = row["covariate"]
        riker_map[cov] = [int(row[c]) for c in r_cols]

    picard_map = {}
    for row in picard_rows:
        cov = row["COVARIATE"]
        picard_map[cov] = [int(row[c]) for c in p_cols]

    # Merge covariates
    all_covs = list(dict.fromkeys(list(riker_map.keys()) + list(picard_map.keys())))

    lines = []
    n_match = 0
    n_diff = 0
    n_riker_only = 0
    n_picard_only = 0
    total_diffs = [0] * len(col_labels)

    for cov in all_covs:
        r_vals = riker_map.get(cov)
        p_vals = picard_map.get(cov)

        # Treat all-zero rows as absent — Picard emits empty bins (e.g., overlap
        # rows with TOTAL_BASES > 0 but all overlap counts zero) that riker skips.
        if r_vals is not None and all(v == 0 for v in r_vals):
            r_vals = None
        if p_vals is not None and all(v == 0 for v in p_vals):
            p_vals = None

        if r_vals is None and p_vals is None:
            continue
        if r_vals is None:
            n_picard_only += 1
            lines.append(f"  PICARD_ONLY  {cov}: picard={p_vals}")
            continue
        if p_vals is None:
            n_riker_only += 1
            lines.append(f"  RIKER_ONLY   {cov}: riker={r_vals}")
            continue

        diffs = [r - p for r, p in zip(r_vals, p_vals)]
        for i, d in enumerate(diffs):
            total_diffs[i] += abs(d)

        if all(d == 0 for d in diffs):
            n_match += 1
        else:
            n_diff += 1
            parts = []
            for label, rv, pv, d in zip(col_labels, r_vals, p_vals, diffs):
                if d != 0:
                    pct = f" ({d/pv*100:+.4f}%)" if pv != 0 else ""
                    parts.append(f"{label}: riker={rv} picard={pv} diff={d:+d}{pct}")
            lines.append(f"  DIFF  {cov}: {'; '.join(parts)}")

    total = len(all_covs)
    summary = f"Total: {total} covariates | Match: {n_match} | Differ: {n_diff} | Riker-only: {n_riker_only} | Picard-only: {n_picard_only}"

    if total_diffs:
        agg_parts = [f"{label}: {td}" for label, td in zip(col_labels, total_diffs)]
        summary += f"\nAggregate absolute diffs: {', '.join(agg_parts)}"

    with open(diff_file, "w") as f:
        f.write(f"Comparison: {strat} / {error_key}\n")
        f.write(f"Riker: {riker_file}\n")
        f.write(f"Picard: {picard_file}\n")
        f.write(f"\n{summary}\n\n")
        if lines:
            f.write("Details (only differences shown):\n")
            for line in lines:
                f.write(line + "\n")

    if n_diff == 0 and n_riker_only == 0 and n_picard_only == 0:
        return "PASS"
    return "FAIL"


def main():
    parser = argparse.ArgumentParser(description="Compare riker vs Picard error metrics")
    parser.add_argument("--riker-prefix", default="/tmp/riker_r1")
    parser.add_argument("--picard-prefix", default="/tmp/picard_r1")
    parser.add_argument("--diff-dir", default="/tmp/error_comparison")
    args = parser.parse_args()

    diff_dir = Path(args.diff_dir)
    diff_dir.mkdir(parents=True, exist_ok=True)

    stratifiers = [s for s in RIKER_TO_PICARD_SUFFIX.keys() if s != "indel_len"]  # Picard crashes on INDEL_LENGTH
    error_keys = ["mismatch", "overlap", "indel"]

    total = 0
    passed = 0
    failed = 0
    skipped = 0
    fail_list = []

    print()
    print("Comparing riker vs Picard CollectSamErrorMetrics")
    print("=" * 55)
    print()

    for strat in stratifiers:
        for ek in error_keys:
            total += 1
            result = compare_one(strat, ek, args.riker_prefix, args.picard_prefix, diff_dir)

            if result == "PASS":
                passed += 1
                icon = "PASS"
            elif result == "SKIP":
                skipped += 1
                icon = "SKIP"
            else:
                failed += 1
                fail_list.append(f"{strat}/{ek}")
                icon = "FAIL"

            print(f"  {strat:<25s} {ek:<15s}  {icon}")

    print()
    print("=" * 55)
    print(f"  TOTAL: {total}  |  PASS: {passed}  |  FAIL: {failed}  |  SKIP: {skipped}")
    print("=" * 55)

    if fail_list:
        print()
        print(f"Failed comparisons (see {diff_dir}/ for details):")
        for f in fail_list:
            strat, ek = f.split("/")
            diff_file = diff_dir / f"{strat}_{ek}.diff"
            if diff_file.exists():
                # Print summary line
                for line in diff_file.read_text().splitlines():
                    if line.startswith("Total:") or line.startswith("Aggregate"):
                        print(f"  {f}: {line}")

    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
