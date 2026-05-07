#!/usr/bin/env python3
"""One-time helper that selected the NYGC samples currently locked into
config/samples.wgs.tsv. Kept under workflow/scripts/ for reproducibility,
NOT invoked by the pipeline.

Lists the public NYGC 1000 Genomes 30x bucket via `aws s3`, sorts CRAMs
by size, and prints recommended small/medium/large picks (5th / 50th /
95th percentile by file size).

Usage:
    python select_nygc.py            # uses --no-sign-request (anonymous)
    python select_nygc.py --bucket   # alternate bucket if NYGC moves
"""

import argparse
import subprocess
import sys


DEFAULT_BUCKET = "1000genomes"
DEFAULT_PREFIX = "1000G_2504_high_coverage/data/"


def list_crams(bucket: str, prefix: str) -> list[tuple[int, str]]:
    """Return [(size_bytes, key), ...] sorted ascending by size."""
    proc = subprocess.run(
        ["aws", "s3api", "list-objects-v2",
         "--no-sign-request",
         "--bucket", bucket,
         "--prefix", prefix,
         "--query", "Contents[?ends_with(Key, `.final.cram`)].[Size, Key]",
         "--output", "text"],
        check=True, capture_output=True, text=True,
    )
    out = []
    for line in proc.stdout.splitlines():
        size_str, _, key = line.partition("\t")
        out.append((int(size_str), key.strip()))
    out.sort()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bucket", default=DEFAULT_BUCKET)
    p.add_argument("--prefix", default=DEFAULT_PREFIX)
    args = p.parse_args()

    crams = list_crams(args.bucket, args.prefix)
    if not crams:
        sys.exit("No CRAMs found at s3://%s/%s" % (args.bucket, args.prefix))

    n = len(crams)
    pct5  = crams[int(n * 0.05)]
    pct50 = crams[int(n * 0.50)]
    pct95 = crams[int(n * 0.95)]

    print(f"Total NYGC samples: {n}")
    print(f"Smallest: {crams[0][0] / 1e9:6.2f} GB  {crams[0][1]}")
    print(f"Largest:  {crams[-1][0] / 1e9:6.2f} GB  {crams[-1][1]}")
    print()
    print("Recommended picks (5th / 50th / 95th percentile by size):")
    for label, (size, key) in [(" 5th", pct5), ("50th", pct50), ("95th", pct95)]:
        sample = key.split("/")[-1].replace(".final.cram", "")
        print(f"  {label} pct: {size / 1e9:6.2f} GB  {sample}  ({key})")


if __name__ == "__main__":
    main()
