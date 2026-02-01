#!/usr/bin/env bash
set -euo pipefail

# =========================
# 0) Paths (edit if needed)
# =========================
REF_PARTIAL="data/hg38_partial.fasta.gz"
REF_REFSEQ="data/hg38_partial.fasta.gz"
# REF_REFSEQ="data/hg38_refseq_GCF_000001405.26_GRCh38_genomic.fna.gz"

Q40="data/illumina_reads_40.fasta.gz"
Q60="data/illumina_reads_60.fasta.gz"
Q80="data/illumina_reads_80.fasta.gz"
Q100="data/illumina_reads_100.fasta.gz"

BIN_NAIVE="./build/bin/naive_search"
BIN_SA_CONSTRUCT="./build/bin/suffixarray_construct"
BIN_SA_SEARCH="./build/bin/suffixarray_search"
BIN_FM_CONSTRUCT="./build/bin/fmindex_construct"
BIN_FM_SEARCH="./build/bin/fmindex_search"
BIN_PIGEON="./build/bin/fmindex_pigeon_search"

IDX_DIR="build/indexes"
LOG_DIR="build/logs"
OUT_CSV="results/benchmark_cpp_total.csv"

IDX_SA_PARTIAL="${IDX_DIR}/hg38_partial.sa.bin"
IDX_FM_PARTIAL="${IDX_DIR}/hg38_partial.fm.bin"

IDX_FM_REFSEQ="${IDX_DIR}/hg38_refseq.fm.bin"

mkdir -p "${IDX_DIR}" "${LOG_DIR}" "results"

# =========================
# 1) helpers
# =========================
get_maxrss_kb() {
  awk -F': *' '/Maximum resident set size \(kbytes\)/{print $2}' "$1" | tail -n 1
}


get_first_number_after_prefix_any() {
  # usage: get_first_number_after_prefix_any "Search time:" out_file err_file
  local prefix="$1"
  local f1="$2"
  local f2="$3"

  awk -v p="$prefix" '
    index($0,p)==1 {
      for(i=1;i<=NF;i++){
        if($i ~ /^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/){print $i; exit}
      }
    }' "$f1" "$f2" 2>/dev/null || true
}

get_hits_any() {
  # usage: get_hits_any out_file err_file
  local f1="$1"
  local f2="$2"
  local h
  h="$(grep -Eo 'hits=[0-9]+' "$f1" "$f2" 2>/dev/null | tail -n 1 | cut -d= -f2 || true)"
  if [[ -z "${h}" ]]; then
    h="$(awk -F'[: ]+' '/^hits[: ]/{print $2}' "$f1" "$f2" 2>/dev/null | tail -n 1 || true)"
  fi
  echo "${h:-NA}"
}

run_and_record() {
  # args:
  #   program, reference, query, qlen, qct, errors, threads, cmd...
  local program="$1"; shift
  local reference="$1"; shift
  local query="$1"; shift
  local qlen="$1"; shift
  local qct="$1"; shift
  local errors="$1"; shift
  local threads="$1"; shift

  local logbase="${program}_len${qlen}_ct${qct}_e${errors}_t${threads}_$(date +%s%N)"
  local out="${LOG_DIR}/${logbase}.out"
  local err="${LOG_DIR}/${logbase}.err"
  local timef="${LOG_DIR}/${logbase}.time"

  set +e
  /usr/bin/time -v -o "${timef}" "$@" >"${out}" 2>"${err}"
  local ec=$?
  set -e

  local idx_s="NA"
  local search_s="NA"

  idx_s="$(get_first_number_after_prefix_any "Index Construction time:" "${out}" "${err}")"
  [[ -z "${idx_s}" ]] && idx_s="NA"

  search_s="$(get_first_number_after_prefix_any "Search time:" "${out}" "${err}")"
  [[ -z "${search_s}" ]] && search_s="NA"

  local maxrss
  maxrss="$(get_maxrss_kb "${timef}")"
  maxrss="${maxrss:-NA}"

  local hits
  hits="$(get_hits_any "${out}" "${err}")"

  local cmd_str
  cmd_str="$(printf "%q " "$@")"

  echo "${program},${reference},${query},${qlen},${qct},${errors},${threads},${idx_s},${search_s},${maxrss},${hits},${ec},${logbase},\"${cmd_str}\"" >> "${OUT_CSV}"
}

run_construct_and_record() {
  # args: program, reference, index_path, cmd...
  local program="$1"; shift
  local reference="$1"; shift
  local index_path="$1"; shift

  local logbase="${program}_$(basename "${reference}")_$(date +%s%N)"
  local out="${LOG_DIR}/${logbase}.out"
  local err="${LOG_DIR}/${logbase}.err"
  local timef="${LOG_DIR}/${logbase}.time"

  set +e
  /usr/bin/time -v -o "${timef}" "$@" >"${out}" 2>"${err}"
  local ec=$?
  set -e

  local idx_s="NA"
  idx_s="$(get_first_number_after_prefix_any "Index Construction time:" "${out}" "${err}")"
  [[ -z "${idx_s}" ]] && idx_s="NA"

  local maxrss
  maxrss="$(get_maxrss_kb "${timef}")"
  maxrss="${maxrss:-NA}"

  local cmd_str
  cmd_str="$(printf "%q " "$@")"

  echo "${program},${reference},NA,NA,NA,NA,NA,${idx_s},NA,${maxrss},NA,${ec},${logbase},\"${cmd_str}\"" >> "${OUT_CSV}"
}

# =========================
# 2) CSV header
# =========================
echo "program,reference,query,query_len,query_ct,errors,threads,index_time_s,search_time_s,max_rss_kb,hits,exit_code,log_base,cmd" > "${OUT_CSV}"

# =========================
# 3) Build indices (once)
# =========================
echo "[1/4] Build indices (if missing)..."

if [[ ! -s "${IDX_SA_PARTIAL}" ]]; then
  echo "  - building SA index: ${IDX_SA_PARTIAL}"
  run_construct_and_record "suffixarray_construct" "${REF_PARTIAL}" "${IDX_SA_PARTIAL}" \
    "${BIN_SA_CONSTRUCT}" --reference "${REF_PARTIAL}" --index "${IDX_SA_PARTIAL}"
else
  echo "  - SA index exists: ${IDX_SA_PARTIAL}"
fi

if [[ ! -s "${IDX_FM_PARTIAL}" ]]; then
  echo "  - building FM index: ${IDX_FM_PARTIAL}"
  run_construct_and_record "fmindex_construct" "${REF_PARTIAL}" "${IDX_FM_PARTIAL}" \
    "${BIN_FM_CONSTRUCT}" --reference "${REF_PARTIAL}" --index "${IDX_FM_PARTIAL}"
else
  echo "  - FM index exists: ${IDX_FM_PARTIAL}"
fi

# whole genome FM-index (optional)
if [[ -f "${REF_REFSEQ}" ]]; then
  if [[ ! -s "${IDX_FM_REFSEQ}" ]]; then
    echo "  - building FM index (refseq): ${IDX_FM_REFSEQ}"
    run_construct_and_record "fmindex_construct_refseq" "${REF_REFSEQ}" "${IDX_FM_REFSEQ}" \
      "${BIN_FM_CONSTRUCT}" --reference "${REF_REFSEQ}" --index "${IDX_FM_REFSEQ}"
  else
    echo "  - FM index (refseq) exists: ${IDX_FM_REFSEQ}"
  fi
else
  echo "  - skip refseq (file not found): ${REF_REFSEQ}"
fi

# =========================
# 4) Run searches
# =========================
echo "[2/4] Run searches..."


COUNTS=(1000)
LENS=(40 60 80 100)

THREADS_NAIVE=14
THREADS_INDEXED=1

# ---- A1 Task1: varying counts (len=100)
LEN=100
Q="${Q100}"

for CT in "${COUNTS[@]}"; do
  # ✅ FIX: naive_search 不要传 --query_len
  run_and_record "naive_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_NAIVE}" \
    "${BIN_NAIVE}" --reference "${REF_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors 0 --threads "${THREADS_NAIVE}"

  run_and_record "suffixarray_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_INDEXED}" \
    "${BIN_SA_SEARCH}" --reference "${REF_PARTIAL}" --index "${IDX_SA_PARTIAL}" --query "${Q}" --query_ct "${CT}"

  run_and_record "fmindex_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_INDEXED}" \
    "${BIN_FM_SEARCH}" --index "${IDX_FM_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors 0
done

# ---- A1 Task2: varying lengths (ct=1000)
CT=1000
for LEN in "${LENS[@]}"; do
  case "${LEN}" in
    40) Q="${Q40}" ;;
    60) Q="${Q60}" ;;
    80) Q="${Q80}" ;;
    100) Q="${Q100}" ;;
  esac

  run_and_record "naive_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_NAIVE}" \
    "${BIN_NAIVE}" --reference "${REF_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors 0 --threads "${THREADS_NAIVE}"

  run_and_record "suffixarray_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_INDEXED}" \
    "${BIN_SA_SEARCH}" --reference "${REF_PARTIAL}" --index "${IDX_SA_PARTIAL}" --query "${Q}" --query_ct "${CT}"

  run_and_record "fmindex_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" 0 "${THREADS_INDEXED}" \
    "${BIN_FM_SEARCH}" --index "${IDX_FM_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors 0
done

# ---- A3: approximate (len=100, ct=1000), k=1,2,3
LEN=100
Q="${Q100}"
CT=1000

for K in 1 2 3; do
  run_and_record "fmindex_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" "${K}" "${THREADS_INDEXED}" \
    "${BIN_FM_SEARCH}" --index "${IDX_FM_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors "${K}"

  run_and_record "fmindex_pigeon_search" "${REF_PARTIAL}" "${Q}" "${LEN}" "${CT}" "${K}" "${THREADS_INDEXED}" \
    "${BIN_PIGEON}" --index "${IDX_FM_PARTIAL}" --reference "${REF_PARTIAL}" --query "${Q}" --query_ct "${CT}" --errors "${K}"
done

echo "[3/4] Done."
echo "CSV: ${OUT_CSV}"
echo "Logs: ${LOG_DIR}/"
echo "[4/4] Tip: check naive rows quickly:"
echo "      grep naive_search ${OUT_CSV} | tail"
