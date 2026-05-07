#!/usr/bin/env python3
"""Parse `/usr/bin/time -v` output into a dict.

GNU time's verbose format is one `key: value` per line with some keys
carrying units (e.g. `Maximum resident set size (kbytes): 12345`). We
keep a curated subset — the fields most consumers want — and ignore
the rest."""

import re
from pathlib import Path


_HMS_RE = re.compile(r"^(?:(\d+):)?(\d+):(\d+(?:\.\d+)?)$")


def _wall_to_seconds(value: str) -> float:
    m = _HMS_RE.match(value.strip())
    if not m:
        # Fallback: seconds directly.
        return float(value)
    h, mm, ss = m.groups()
    return (int(h) if h else 0) * 3600 + int(mm) * 60 + float(ss)


def parse(path: str | Path) -> dict:
    out: dict = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            # Split on ": " (colon-space) from the RIGHT — some keys
            # contain their own colons.
            if ": " not in line:
                continue
            key, _, value = line.strip().rpartition(": ")
            key = key.strip()
            value = value.strip()
            if key == "Elapsed (wall clock) time (h:mm:ss or m:ss)":
                out["wall_s"] = _wall_to_seconds(value)
            elif key == "User time (seconds)":
                out["user_s"] = float(value)
            elif key == "System time (seconds)":
                out["sys_s"] = float(value)
            elif key == "Maximum resident set size (kbytes)":
                out["max_rss_kb"] = int(value)
            elif key == "Percent of CPU this job got":
                out["cpu_percent"] = int(value.rstrip("%"))
            elif key == "Exit status":
                out["exit_status"] = int(value)
    return out


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(parse(sys.argv[1]), indent=2))
