#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import csv
import glob
import os
from typing import Optional, Dict

LOGDIR = "build/logs"
OUTCSV = "results/benchmark_fm.csv"

re_cmd = re.compile(r'Command being timed:\s+"([^"]+)"')
re_search_time = re.compile(r"Search time:\s*([0-9.]+)\s*seconds\.")
re_hits = re.compile(r"queries=(\d+)\s+errors=(\d+)\s+hits=(\d+)")
re_maxrss = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
re_wall = re.compile(r"Elapsed \(wall clock\) time.*:\s*([0-9:.]+)")
re_exit = re.compile(r"Exit status:\s*(\d+)")

def parse_wall_to_seconds(s: str) -> Optional[float]:
    """
    Convert time -v wall string:
    - "m:ss.xx"
    - "h:mm:ss"
    - "h:mm:ss.xx"
    into seconds (float).
    """
    s = s.strip()
    if not s:
        return None
    parts = s.split(":")
    try:
        if len(parts) == 2:
            m = int(parts[0])
            sec = float(parts[1])
            return m * 60 + sec
        if len(parts) == 3:
            h = int(parts[0])
            m = int(parts[1])
            sec = float(parts[2])
            return h * 3600 + m * 60 + sec
    except ValueError:
        return None
    return None

def get_arg(cmd: str, name: str) -> Optional[str]:
    """
    Extract command arg value: --name VALUE
    Works even if multiple spaces.
    """
    m = re.search(rf"(?:^|\s)--{re.escape(name)}\s+(\S+)", cmd)
    return m.group(1) if m else None

def infer_query_len(query_path: Optional[str]) -> Optional[int]:
    """
    From '...illumina_reads_100.fasta.gz' -> 100
    Works for .fasta, .fa, .fastq, with/without .gz.
    """
    if not query_path:
        return None
    base = os.path.basename(query_path)
    m = re.search(r"illumina_reads_(\d+)\.(?:fasta|fa|fastq)(?:\.gz)?$", base)
    return int(m.group(1)) if m else None

def parse_one_err(path: str) -> Optional[Dict]:
    txt = open(path, "r", errors="ignore").read()

    cmd_m = re_cmd.search(txt)
    if not cmd_m:
        return None
    cmd = cmd_m.group(1)

    prog = os.path.basename(cmd.split()[0])
    # Normalize program name
    if "pigeon" in prog:
        program = "fmindex_pigeon_search"
    elif "construct" in prog:
        program = "fmindex_construct"
    else:
        program = "fmindex_search"

    query_ct = get_arg(cmd, "query_ct")
    errors = get_arg(cmd, "errors")
    query_file = get_arg(cmd, "query")
    query_len = infer_query_len(query_file)

    st_m = re_search_time.search(txt)
    search_time_s = float(st_m.group(1)) if st_m else None

    hits_m = re_hits.search(txt)
    hits = int(hits_m.group(3)) if hits_m else None

    rss_m = re_maxrss.search(txt)
    max_rss_kb = int(rss_m.group(1)) if rss_m else None

    wall_m = re_wall.search(txt)
    wall_raw = wall_m.group(1).strip() if wall_m else None
    wall_time_s = parse_wall_to_seconds(wall_raw) if wall_raw else None

    exit_m = re_exit.search(txt)
    exit_code = int(exit_m.group(1)) if exit_m else None

    # Some runs may not have --query (e.g., construct)
    return {
        "log_file": os.path.basename(path),
        "program": program,
        "query_len": query_len,
        "query_ct": int(query_ct) if query_ct is not None else None,
        "errors": int(errors) if errors is not None else None,
        "search_time_s": search_time_s,
        "wall_time_s": wall_time_s,
        "max_rss_kb": max_rss_kb,
        "hits": hits,
        "exit_code": exit_code,
        "cmd": cmd,
    }

def main():
    os.makedirs(os.path.dirname(OUTCSV), exist_ok=True)

    rows = []
    for path in sorted(glob.glob(os.path.join(LOGDIR, "*.err"))):
        r = parse_one_err(path)
        if r is not None:
            rows.append(r)

    # Write CSV (keep cmd/log_file to debug; you can drop them later)
    fields = [
        "program","query_len","query_ct","errors",
        "search_time_s","wall_time_s","max_rss_kb","hits",
        "exit_code","log_file","cmd"
    ]
    with open(OUTCSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})

    print(f"Parsed {len(rows)} logs -> {OUTCSV}")

if __name__ == "__main__":
    main()

