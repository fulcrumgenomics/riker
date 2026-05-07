#!/usr/bin/env python3
"""Emit a small JSON blob describing the host: CPU model, arch, cores,
memory, kernel, OS. Written once per pipeline run and joined onto every
row of bench.tsv at aggregation time."""

import json
import os
import platform
import re
import subprocess
import sys


def cpu_model() -> str:
    if sys.platform == "darwin":
        try:
            return subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.brand_string"], text=True
            ).strip()
        except Exception:
            return "unknown"
    # Linux
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith(("model name", "Model", "cpu model")):
                    return line.split(":", 1)[1].strip()
            # aarch64 /proc/cpuinfo often only has implementer/part — probe.
            fh.seek(0)
            for line in fh:
                if line.startswith("CPU part"):
                    return f"aarch64 CPU part {line.split(':')[1].strip()}"
    except FileNotFoundError:
        pass
    return "unknown"


def total_mem_bytes() -> int:
    if sys.platform == "darwin":
        try:
            return int(subprocess.check_output(["sysctl", "-n", "hw.memsize"]).strip())
        except Exception:
            return 0
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(re.search(r"(\d+)", line).group(1)) * 1024
    except FileNotFoundError:
        pass
    return 0


def aws_instance_type() -> str | None:
    """Best-effort AWS IMDSv2 probe. Returns None off-AWS or on timeout."""
    try:
        token = subprocess.check_output(
            ["curl", "-s", "-X", "PUT",
             "-H", "X-aws-ec2-metadata-token-ttl-seconds: 60",
             "--max-time", "1",
             "http://169.254.169.254/latest/api/token"],
            text=True,
        ).strip()
        if not token:
            return None
        out = subprocess.check_output(
            ["curl", "-s", "-H", f"X-aws-ec2-metadata-token: {token}",
             "--max-time", "1",
             "http://169.254.169.254/latest/meta-data/instance-type"],
            text=True,
        ).strip()
        return out or None
    except Exception:
        return None


def main():
    info = {
        "hostname": platform.node(),
        "os": platform.system(),
        "os_release": platform.release(),
        "arch": platform.machine(),
        "cpu_model": cpu_model(),
        "cpu_count_logical": os.cpu_count(),
        "total_mem_bytes": total_mem_bytes(),
        "python": platform.python_version(),
        "aws_instance_type": aws_instance_type(),
    }
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
