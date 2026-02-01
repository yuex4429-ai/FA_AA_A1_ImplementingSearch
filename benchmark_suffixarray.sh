#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Suffix Array benchmark (Assignment-compatible, simple)
# - build index once
# - benchmark search with varying query counts / lengths
# - use /usr/bin/time -v
# ==========================================================

# ---------- paths ----------
REF="data/hg38_partial.fasta.gz"
Q100="data/illumina_reads_100.fasta.gz"
Q40="data/illumina_reads_40.fasta.gz"
Q60="data/illumina_reads_60.fasta.gz"
Q80="data/illumina_reads_80.fasta.gz"

BIN_CONSTRUCT="build/bin/suffixarray_construct"
BIN_SEARCH="build/bin/suffixarray_search"

INDEX="build/indexes/hg38_partial.sa.bin"
LOGDIR="build/logs"

mkdir -p build/indexes "$LOGDIR"

# ==========================================================
# 1) Build suffix array index
# ==========================================================
echo "[BUILD] suffix array index"

OUT="$LOGDIR/suffixarray_construct_hg38_partial.out"
ERR="$LOGDIR/suffixarray_construct_hg38_partial.err"

(/usr/bin/time -v \
  "$BIN_CONSTRUCT" \
    --reference "$REF" \
    --index "$INDEX" \
) >"$OUT" 2>"$ERR"

# ==========================================================
# 2) Benchmark search: varying query count (len=100)
# ==========================================================
COUNTS=(1000 10000 100000)

for CT in "${COUNTS[@]}"; do
  echo "[SEARCH] suffix array | len=100 | queries=$CT"

  OUT="$LOGDIR/suffixarray_search_len100_ct${CT}.out"
  ERR="$LOGDIR/suffixarray_search_len100_ct${CT}.err"

  (/usr/bin/time -v \
    "$BIN_SEARCH" \
      --reference "$REF" \
      --index "$INDEX" \
      --query "$Q100" \
      --query_ct "$CT" \
  ) >"$OUT" 2>"$ERR"
done

# ==========================================================
# 3) Benchmark search: varying query length (ct=1000)
# ==========================================================
FIXED_CT=1000

for Q in "$Q40" "$Q60" "$Q80" "$Q100"; do
  LEN=$(echo "$Q" | grep -oP '[0-9]+(?=\.fasta)')
  echo "[SEARCH] suffix array | len=$LEN | queries=$FIXED_CT"

  OUT="$LOGDIR/suffixarray_search_len${LEN}_ct${FIXED_CT}.out"
  ERR="$LOGDIR/suffixarray_search_len${LEN}_ct${FIXED_CT}.err"

  (/usr/bin/time -v \
    "$BIN_SEARCH" \
      --reference "$REF" \
      --index "$INDEX" \
      --query "$Q" \
      --query_ct "$FIXED_CT" \
  ) >"$OUT" 2>"$ERR"
done

echo "[DONE] suffix array benchmark finished."

