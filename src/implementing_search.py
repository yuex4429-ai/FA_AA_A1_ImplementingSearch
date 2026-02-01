#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gc
import gzip
import os
import struct
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Iterable, List, Tuple, Optional, Set

# ---------- optional libs ----------
try:
    import iv2py as iv  # fmindex
    HAS_IV2PY = True
except Exception:
    HAS_IV2PY = False

try:
    from pydivsufsort import divsufsort  # suffix array
    HAS_DIVSUFSORT = True
except Exception:
    HAS_DIVSUFSORT = False


# =========================================================
# Globals for process workers (avoid pickling huge objects per task)
# =========================================================
_G_TEXT: Optional[str] = None
_G_SA_BYTES: Optional[bytes] = None
_G_SA_N: int = 0
_G_QUERIES: Optional[List[str]] = None


def _init_text_queries_worker(text: str, queries: List[str]) -> None:
    global _G_TEXT, _G_QUERIES
    _G_TEXT = text
    _G_QUERIES = queries


def _init_sa_worker(text: str, sa_bytes: bytes, sa_n: int, queries: List[str]) -> None:
    global _G_TEXT, _G_SA_BYTES, _G_SA_N, _G_QUERIES
    _G_TEXT = text
    _G_SA_BYTES = sa_bytes
    _G_SA_N = sa_n
    _G_QUERIES = queries


# =========================================================
# dna5 normalization (match SeqAn3 dna5): map non-ACGTN -> N
# =========================================================
_DNA5_KEEP = set("ACGTN")


def to_dna5_str(seq: str) -> str:
    s = seq.strip().upper()
    if not s:
        return ""
    return "".join(ch if ch in _DNA5_KEEP else "N" for ch in s)


# =========================================================
# FASTA reader (supports .gz)
# =========================================================
def iter_fasta_sequences(path: str) -> Iterable[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        seq_chunks: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_chunks:
                    yield to_dna5_str("".join(seq_chunks))
                    seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_chunks:
            yield to_dna5_str("".join(seq_chunks))


def load_reference(path: str) -> List[str]:
    ref = list(iter_fasta_sequences(path))
    if not ref:
        raise RuntimeError(f"Empty reference: {path}")
    return ref


def load_queries(path: str) -> List[str]:
    qs = list(iter_fasta_sequences(path))
    if not qs:
        raise RuntimeError(f"Empty query file: {path}")
    return qs


# =========================================================
# Duplicate queries (match your C++ doubling logic)
# while size < n: resize(2*old); copy old block once
# =========================================================
def duplicate_to_n_cppstyle(items: List[str], n: int) -> List[str]:
    if n <= 0 or not items:
        return []
    out = list(items)
    while len(out) < n:
        old = len(out)
        out.extend(out[:old])
    return out[:n]


# =========================================================
# chunk ranges (FIXED): never collapse to a single block for small n
# =========================================================
def chunk_ranges(n: int, workers: int, min_block: int = 64) -> List[Tuple[int, int]]:
    """
    Split [0, n) into blocks.

    - Always create at least min(workers, n) blocks (so parallelism is observable even for small n).
    - For large n, try to keep blocks >= min_block to reduce overhead.
    """
    if n <= 0:
        return []
    workers = max(1, workers)
    min_block = max(1, min_block)

    min_blocks = min(workers, n)
    blocks_by_min_block = (n + min_block - 1) // min_block
    blocks = max(min_blocks, blocks_by_min_block)
    blocks = min(blocks, n)

    block_size = (n + blocks - 1) // blocks

    rs: List[Tuple[int, int]] = []
    for t in range(blocks):
        b = t * block_size
        e = min(n, b + block_size)
        if b < e:
            rs.append((b, e))
    return rs


# =========================================================
# NAIVE (count all occurrences)
# =========================================================
def count_occurrences_naive(text: str, query: str) -> int:
    """
    Count all occurrences (including overlaps) using str.find loop.
    """
    if not query or len(query) > len(text):
        return 0
    total = 0
    start = 0
    while True:
        pos = text.find(query, start)
        if pos == -1:
            break
        total += 1
        start = pos + 1
    return total


def _proc_worker_naive_range(b: int, e: int) -> int:
    assert _G_TEXT is not None and _G_QUERIES is not None
    s = 0
    for i in range(b, e):
        s += count_occurrences_naive(_G_TEXT, _G_QUERIES[i])
    return s


# =========================================================
# SUFFIX ARRAY (divsufsort), packed uint32
# =========================================================
def build_sa_packed(text: str) -> Tuple[bytes, int]:
    if not HAS_DIVSUFSORT:
        raise RuntimeError("suffixarray mode requires pydivsufsort. Install: pip install pydivsufsort")
    sa = divsufsort(text.encode("ascii"))
    n = len(sa)
    buf = bytearray(n * 4)
    for i, v in enumerate(sa):
        struct.pack_into("<I", buf, i * 4, int(v))
    return bytes(buf), n


def _sa_u32(sa_bytes: bytes, idx: int) -> int:
    return struct.unpack_from("<I", sa_bytes, idx * 4)[0]


def _lower_bound(text: str, sa_bytes: bytes, sa_n: int, pat: str) -> int:
    lo, hi = 0, sa_n
    m = len(pat)
    while lo < hi:
        mid = (lo + hi) // 2
        pos = _sa_u32(sa_bytes, mid)
        s = text[pos:pos + m]
        if s < pat:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _upper_bound(text: str, sa_bytes: bytes, sa_n: int, pat: str, start: int) -> int:
    lo, hi = start, sa_n
    m = len(pat)
    while lo < hi:
        mid = (lo + hi) // 2
        pos = _sa_u32(sa_bytes, mid)
        s = text[pos:pos + m]
        if s <= pat:
            lo = mid + 1
        else:
            hi = mid
    return lo


def count_occurrences_sa(text: str, sa_bytes: bytes, sa_n: int, pat: str) -> int:
    if not pat:
        return 0
    lb = _lower_bound(text, sa_bytes, sa_n, pat)
    ub = _upper_bound(text, sa_bytes, sa_n, pat, lb)
    return ub - lb


def _proc_worker_sa_range(b: int, e: int) -> int:
    assert _G_TEXT is not None and _G_SA_BYTES is not None and _G_QUERIES is not None
    s = 0
    for i in range(b, e):
        s += count_occurrences_sa(_G_TEXT, _G_SA_BYTES, _G_SA_N, _G_QUERIES[i])
    return s


# =========================================================
# FM-index via IV2py
# =========================================================
class IVFM:
    def __init__(self, reference_records: List[str], sampling_rate: int = 16):
        if not HAS_IV2PY:
            raise RuntimeError("fmindex modes require iv2py. Install: pip install iv2py")
        self.reference = reference_records
        self.index = iv.fmindex(reference=reference_records, samplingRate=int(sampling_rate))

    def search_all(self, pattern: str, k: int = 0):
        return self.index.search(pattern, k=int(k))


# =========================================================
# Pigeon helpers
# =========================================================
@dataclass(frozen=True)
class CandidateKey:
    ref_id: int
    start: int


def verify_hamming(text: str, query: str, start_pos: int, max_errors: int) -> bool:
    if start_pos < 0:
        return False
    end = start_pos + len(query)
    if end > len(text):
        return False
    errors = 0
    for i, ch in enumerate(query):
        if text[start_pos + i] != ch:
            errors += 1
            if errors > max_errors:
                return False
    return True


def count_verified_hits_pigeon(reference: List[str], fm: IVFM, query: str, errors: int) -> int:
    if not query:
        return 0
    if errors <= 0:
        return len(fm.search_all(query, k=0))

    m = len(query)
    k = min(errors + 1, m)
    b = [(i * m) // k for i in range(k + 1)]

    seen: Set[CandidateKey] = set()
    hits = 0

    for part_id in range(k):
        pb, pe = b[part_id], b[part_id + 1]
        if pe <= pb:
            continue
        part = query[pb:pe]
        if not part:
            continue

        for res in fm.search_all(part, k=0):
            try:
                ref_id, hit_pos = int(res[0]), int(res[1])
            except Exception:
                continue

            start = hit_pos - pb
            if start < 0:
                continue

            key = CandidateKey(ref_id, start)
            if key in seen:
                continue
            seen.add(key)

            if 0 <= ref_id < len(reference) and verify_hamming(reference[ref_id], query, start, errors):
                hits += 1

    return hits


def _thread_worker_fm_count(fm: IVFM, queries: List[str], k: int, b: int, e: int) -> int:
    s = 0
    for i in range(b, e):
        s += len(fm.search_all(queries[i], k=k))
    return s


def _thread_worker_pigeon(reference: List[str], fm: IVFM, queries: List[str], errors: int, b: int, e: int) -> int:
    s = 0
    for i in range(b, e):
        s += count_verified_hits_pigeon(reference, fm, queries[i], errors)
    return s


# =========================================================
# Main
# =========================================================
def main():
    ap = argparse.ArgumentParser(prog="implementing_search.py")
    ap.add_argument("--mode", choices=["naive", "suffixarray", "fmindex", "fmindex_pigeon"], required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--query_ct", type=int, default=100)
    ap.add_argument("--errors", type=int, default=0)
    ap.add_argument("--threads", type=int, default=0, help="0=cpu_count(); used for parallel search")
    ap.add_argument("--sampling_rate", type=int, default=16, help="IV2py FM samplingRate")
    ap.add_argument("--min_block", type=int, default=64, help="min queries per task explaining process pool overhead")
    ap.add_argument("--verbose", action="store_true", help="print partition info")
    args = ap.parse_args()

    # ---- load queries first (small), then reference ----
    queries = duplicate_to_n_cppstyle(load_queries(args.query), args.query_ct)
    if not queries:
        raise RuntimeError("No queries after duplication.")

    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)
    threads = max(1, min(threads, len(queries)))

    # ---- load reference ----
    reference_records = load_reference(args.reference)

    # IMPORTANT memory fix:
    # Do NOT keep both:
    #   - reference_records (list of strings) AND
    #   - concat_text (one giant string)
    # at the same time unless the mode needs it.
    concat_text: Optional[str] = None

    if args.mode in ("naive", "suffixarray"):
        sep = "N" * 50
        concat_text = sep.join(reference_records) if len(reference_records) > 1 else reference_records[0]
        # free list to reduce baseline RSS
        del reference_records
        reference_records = []  # type: ignore
        gc.collect()

    # ---- index construction ----
    sa_bytes: Optional[bytes] = None
    sa_n: int = 0
    fm: Optional[IVFM] = None

    if args.mode == "suffixarray":
        assert concat_text is not None
        t0 = time.perf_counter()
        sa_bytes, sa_n = build_sa_packed(concat_text)
        t1 = time.perf_counter()
        print(f"Index Construction time: {t1 - t0} seconds.", flush=True)

    elif args.mode in ("fmindex", "fmindex_pigeon"):
        t0 = time.perf_counter()
        fm = IVFM(reference_records, sampling_rate=args.sampling_rate)
        t1 = time.perf_counter()
        print(f"Index Construction time: {t1 - t0} seconds.", flush=True)
        # fmindex modes don't need concat_text at all
        # (keep reference_records; pigeon needs it for verification)

    # ---- search phase ----
    t0 = time.perf_counter()
    total_hits = 0

    ranges = chunk_ranges(len(queries), threads, min_block=max(1, int(args.min_block)))
    if args.verbose:
        print(f"[DEBUG] queries={len(queries)} threads={threads} min_block={args.min_block} blocks={len(ranges)}",
              flush=True)

    if args.mode == "naive":
        assert concat_text is not None
        with ProcessPoolExecutor(
            max_workers=threads,
            initializer=_init_text_queries_worker,
            initargs=(concat_text, queries),
        ) as ex:
            futs = [ex.submit(_proc_worker_naive_range, b, e) for (b, e) in ranges]
            for fu in as_completed(futs):
                total_hits += fu.result()

    elif args.mode == "suffixarray":
        assert concat_text is not None and sa_bytes is not None
        with ProcessPoolExecutor(
            max_workers=threads,
            initializer=_init_sa_worker,
            initargs=(concat_text, sa_bytes, sa_n, queries),
        ) as ex:
            futs = [ex.submit(_proc_worker_sa_range, b, e) for (b, e) in ranges]
            for fu in as_completed(futs):
                total_hits += fu.result()

    elif args.mode == "fmindex":
        assert fm is not None
        k = max(0, int(args.errors))
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(_thread_worker_fm_count, fm, queries, k, b, e) for (b, e) in ranges]
            for fu in as_completed(futs):
                total_hits += fu.result()

    else:  # fmindex_pigeon
        assert fm is not None
        emax = max(0, int(args.errors))
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(_thread_worker_pigeon, reference_records, fm, queries, emax, b, e) for (b, e) in ranges]
            for fu in as_completed(futs):
                total_hits += fu.result()

    t1 = time.perf_counter()
    print(f"Search time: {t1 - t0} seconds.", flush=True)

    if args.mode == "fmindex_pigeon":
        print(f"queries={len(queries)} errors={int(args.errors)} threads={threads} verified_hits={total_hits}",
              flush=True)
    else:
        print(f"queries={len(queries)} errors={int(args.errors)} threads={threads} hits={total_hits}",
              flush=True)


if __name__ == "__main__":
    main()
