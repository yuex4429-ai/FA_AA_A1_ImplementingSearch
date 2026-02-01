#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Benchmark for FM-index + Pigeonhole (SeqAn3)
# - Uses /usr/bin/time -v for memory + wall time
# - Separates index construction time (fmindex_construct stdout)
# - Logs:
#   stdout -> .out
#   stderr -> .err   (contains time -v output, and SeqAn debug_stream)
#   exit code -> .code
# ==========================================================

# --------- paths (edit if needed) ----------
BIN_DIR="build/bin"          # where your compiled binaries live
LOG_DIR="build/logs"
IDX_DIR="build/indexes"
RESULTS_DIR="results"

# small reference (assignment 1/2 partial) and full reference (RefSeq .fna.gz)
REFERENCE_SMALL="data/hg38_partial.fasta.gz"
REFERENCE_FULL="data/hg38_partial.fasta.gz"
#REFERENCE_FULL="data/hg38_refseq_GCF_000001405.26_GRCh38_genomic.fna.gz"

# queries (you can add more)
Q40="data/illumina_reads_40.fasta.gz"
Q60="data/illumina_reads_60.fasta.gz"
Q80="data/illumina_reads_80.fasta.gz"
Q100="data/illumina_reads_100.fasta.gz"

# query counts required by assignment 2
COUNTS=(1000 10000 100000 1000000)

# lengths required by assignment 1/2
LENS=(40 60 80 100)

# for pigeonhole (assignment 3)
KS=(1 2 3)

# Choose which reference to use for each part:
# - For assignment 2, they want full humane reference genome.
# - For quick testing, swap REFERENCE_FULL -> REFERENCE_SMALL.
REF_FOR_A2="${REFERENCE_FULL}"
REF_FOR_A3="${REFERENCE_FULL}"

# If 1, rebuild index every run; if 0, build once per reference and reuse.
REBUILD_INDEX_EACH_RUN=0

# ------------------------------------------
mkdir -p "${LOG_DIR}" "${IDX_DIR}" "${RESULTS_DIR}"

# helper: run command with time -v and log files
run_one () {
  local tag="$1"; shift
  local out="${LOG_DIR}/${tag}.out"
  local err="${LOG_DIR}/${tag}.err"
  local codef="${LOG_DIR}/${tag}.code"

  echo "----------------------------------------------------------"
  echo "[RUN] ${tag}"
  echo "  stdout: ${out}"
  echo "  stderr: ${err}"
  echo "  cmd   : $*"

  set +e
  /usr/bin/time -v "$@" 1> "${out}" 2> "${err}"
  local rc=$?
  set -e
  echo "${rc}" > "${codef}"

  echo "[DONE] exit=${rc}"
  return 0
}

# helper: safe filename slug
slug () {
  echo "$1" | sed 's/[\/ ]/_/g'
}

# build FM index (serialized)
build_index () {
  local ref="$1"
  local idx="$2"

  local tag="fmindex_construct_$(slug "${ref}")"
  run_one "${tag}" "${BIN_DIR}/fmindex_construct" --reference "${ref}" --index "${idx}"
}

# -------------- A2: exact FM-index benchmark --------------
# Required: runtime + memory for counts, query length 100
benchmark_fm_exact_counts () {
  local ref="$1"
  local qfile="$2"
  local idx="$3"

  echo "==================== A2: FM-index exact | counts ===================="

  for ct in "${COUNTS[@]}"; do
    local tag="fmindex_search_$(slug "${ref}")_$(slug "${qfile}")_ct${ct}_e0"
    run_one "${tag}" "${BIN_DIR}/fmindex_search" --index "${idx}" --query "${qfile}" --query_ct "${ct}" --errors 0
  done
}

# A1/A2: runtime for varying query lengths (choose a suitable number of queries)
benchmark_fm_exact_lengths () {
  local ref="$1"
  local idx="$2"
  local fixed_ct=1000   # "suitable number of queries" – adjust if too slow/fast

  echo "==================== A1/A2: FM-index exact | lengths ===================="
  for len in "${LENS[@]}"; do
    local qvar="data/illumina_reads_${len}.fasta.gz"
    if [[ ! -f "${qvar}" ]]; then
      echo "[SKIP] missing query file ${qvar}"
      continue
    fi
    local tag="fmindex_search_$(slug "${ref}")_illumina_reads_${len}_ct${fixed_ct}_e0"
    run_one "${tag}" "${BIN_DIR}/fmindex_search" --index "${idx}" --query "${qvar}" --query_ct "${fixed_ct}" --errors 0
  done
}

# -------------- A3: pigeonhole benchmark (k=1,2,3) --------------
benchmark_pigeon_vs_plain () {
  local ref="$1"
  local qfile="$2"
  local idx="$3"
  local fixed_ct=1000   # keep stable across k comparisons

  echo "==================== A3: FM-index vs pigeonhole | k=1..3 ===================="

  for k in "${KS[@]}"; do
    # plain fmindex_search with "errors=k" (your code uses Hamming-only mismatches)
    local tag_plain="fmindex_search_$(slug "${ref}")_$(slug "${qfile}")_ct${fixed_ct}_e${k}"
    run_one "${tag_plain}" "${BIN_DIR}/fmindex_search" --index "${idx}" --query "${qfile}" --query_ct "${fixed_ct}" --errors "${k}"

    # pigeonhole search (exact pieces + verify in reference)
    local tag_pig="fmindex_pigeon_$(slug "${ref}")_$(slug "${qfile}")_ct${fixed_ct}_e${k}"
    run_one "${tag_pig}" "${BIN_DIR}/fmindex_pigeon_search" --index "${idx}" --reference "${ref}" --query "${qfile}" --query_ct "${fixed_ct}" --errors "${k}"
  done
}

# -------------- main --------------
main () {
  # sanity checks
  for exe in fmindex_construct fmindex_search fmindex_pigeon_search; do
    if [[ ! -x "${BIN_DIR}/${exe}" ]]; then
      echo "ERROR: missing executable ${BIN_DIR}/${exe}"
      exit 1
    fi
  done

  # -------- index for A2 (full ref) --------
  local idx_a2="${IDX_DIR}/hg38_refseq.fmidx.bin"
  if [[ "${REBUILD_INDEX_EACH_RUN}" -eq 1 || ! -f "${idx_a2}" ]]; then
    build_index "${REF_FOR_A2}" "${idx_a2}"
  else
    echo "[INFO] reuse existing index: ${idx_a2}"
  fi

  # A2 counts for length 100
  if [[ -f "${Q100}" ]]; then
    benchmark_fm_exact_counts "${REF_FOR_A2}" "${Q100}" "${idx_a2}"
  else
    echo "[SKIP] missing ${Q100}"
  fi

  # A1/A2 lengths
  benchmark_fm_exact_lengths "${REF_FOR_A2}" "${idx_a2}"

  # -------- index for A3 (same ref unless you want different) --------
  local idx_a3="${idx_a2}"
  if [[ "${REF_FOR_A3}" != "${REF_FOR_A2}" ]]; then
    idx_a3="${IDX_DIR}/a3_$(slug "${REF_FOR_A3}").fmidx.bin"
    if [[ "${REBUILD_INDEX_EACH_RUN}" -eq 1 || ! -f "${idx_a3}" ]]; then
      build_index "${REF_FOR_A3}" "${idx_a3}"
    else
      echo "[INFO] reuse existing index: ${idx_a3}"
    fi
  fi

  # A3 compare fmindex vs pigeonhole using length 100 (can change to 40/60/80 too)
  if [[ -f "${Q100}" ]]; then
    benchmark_pigeon_vs_plain "${REF_FOR_A3}" "${Q100}" "${idx_a3}"
  else
    echo "[SKIP] missing ${Q100}"
  fi

  echo "==================== DONE ===================="
  echo "Logs are in: ${LOG_DIR}"
  echo "Indexes are in: ${IDX_DIR}"
}

main "$@"

