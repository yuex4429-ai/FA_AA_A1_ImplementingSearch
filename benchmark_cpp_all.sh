#!/bin/bash
set -euo pipefail

# ==========================================================
# Unified C++ benchmark (Assignment 2 compliant)
# Python-style logging:
#   - stdout -> .out
#   - stderr -> .err  (also contains RSS peak markers)
#   - /usr/bin/time -v -> .time (clean)
#   - Samples RSS of the *real algorithm process* (child),
#     and (optionally) its children/grandchildren if any.
# Threads policy:
#   - ONLY naive_search uses THREADS
#   - others forced to 1
# ==========================================================

REFERENCE_SMALL="data/hg38_partial.fasta.gz"
REFERENCE_FULL="data/hg38_partial.fasta.gz"

QFILE_100="data/illumina_reads_100.fasta.gz"

COUNTS=(1000 10000 100000 1000000)
LENGTHS=(40 60 80 100)
FIXED_COUNT=1000

ERRORS=(0)
THREADS=14
MIN_BLOCK=256
REPEATS=1

RUN_NAIVE_FULL=0
RUN_SA_FULL=0

BIN_NAIVE="./build/bin/naive_search"
BIN_SA="./build/bin/suffixarray_search"
BIN_FM="./build/bin/fmindex_search"
BIN_PIGEON="./build/bin/fmindex_pigeon_search"

CSV_FILE="results/benchmark_cpp_all.csv"

RSS_SAMPLE_INTERVAL="0.05"

mkdir -p results build/logs

for b in "$BIN_NAIVE" "$BIN_SA" "$BIN_FM" "$BIN_PIGEON"; do
  [ -x "$b" ] || { echo "Missing binary: $b"; exit 1; }
done
[ -f "$REFERENCE_SMALL" ] || { echo "Missing $REFERENCE_SMALL"; exit 1; }
[ -f "$REFERENCE_FULL" ]  || { echo "Missing $REFERENCE_FULL";  exit 1; }
[ -f "$QFILE_100" ]       || { echo "Missing $QFILE_100";       exit 1; }

echo "Algorithm,Reference,Query_File,Query_Len,Query_Count,Errors,Threads,Repeat,Index_Time_s,Search_Time_s,Wall_Time_s,Max_RSS_KB,Peak_Total_RSS_KB,Peak_Main_RSS_KB,Exit_Code" > "$CSV_FILE"

sanitize() { echo "$1" | tr '/: ' '___' | tr -cd 'A-Za-z0-9_.-'; }
infer_query_len() { [[ "$1" =~ _([0-9]+)\.fasta(\.gz)?$ ]] && echo "${BASH_REMATCH[1]}" || echo "N/A"; }

extract_seconds_value() { sed -n 's/.*: \([0-9]\+\(\.[0-9]\+\)\?\) seconds\.*/\1/p' | tail -n 1; }
extract_time_from_logs() {
  local pattern="$1" out="$2" err="$3"
  local v=""
  if grep -q "$pattern" "$out" 2>/dev/null; then v=$(grep "$pattern" "$out" | extract_seconds_value || true); fi
  if [ -z "$v" ] && grep -q "$pattern" "$err" 2>/dev/null; then v=$(grep "$pattern" "$err" | extract_seconds_value || true); fi
  echo "$v"
}

# Parse from /usr/bin/time -v output file
parse_wall_str_timefile() { grep "Elapsed (wall clock) time" "$1" | awk '{print $NF}' | tail -n 1 || echo ""; }
parse_max_rss_kb_timefile() { grep "Maximum resident set size" "$1" | awk '{print $NF}' | tail -n 1 || echo ""; }

# Markers written into .err by wrapper
parse_peak_total_rss_kb() { grep -E '^__PEAK_RSS_TOTAL_KB__=' "$1" | tail -n 1 | sed 's/.*=//' || echo ""; }
parse_peak_main_rss_kb()  { grep -E '^__PEAK_RSS_MAIN_KB__='  "$1" | tail -n 1 | sed 's/.*=//' || echo ""; }
parse_inner_rc()          { grep -E '^__INNER_RC__='          "$1" | tail -n 1 | sed 's/.*=//' || echo ""; }

to_seconds() {
  local t="$1"
  if [ -z "$t" ]; then echo ""; return; fi
  if [[ "$t" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then echo "$t"; return; fi
  IFS=':' read -r a b c <<< "$t"
  if [ -z "${b:-}" ]; then echo ""; return; fi
  if [ -z "${c:-}" ]; then awk -v m="$a" -v s="$b" 'BEGIN{print (m*60 + s)}'
  else awk -v h="$a" -v m="$b" -v s="$c" 'BEGIN{print (h*3600 + m*60 + s)}'
  fi
}

# ----------------------------------------------------------
# Python-style timed runner:
#  - time -v -> timef
#  - program stdout -> out
#  - program stderr -> err
#  - wrapper prints RSS markers to err
# ----------------------------------------------------------
run_cmd_timed() {
  local out="$1" err="$2" timef="$3" codefile="$4"; shift 4

  set +e
  /usr/bin/time -v -o "$timef" bash -c '
    set -euo pipefail

    "$@" &           # run program in background
    pid=$!

    peak_total=0
    peak_main=0

    get_main_rss_kb() {
      local p="$1"
      if [ -r "/proc/${p}/status" ]; then
        awk "/^VmRSS:/ {print \$2+0}" "/proc/${p}/status" 2>/dev/null || echo 0
      else
        echo 0
      fi
    }

    # Sum RSS (KB) for pid + children + grandchildren
    get_total_rss_kb() {
      local p="$1"
      local c g
      c=$(pgrep -P "$p" 2>/dev/null || true)
      g=""
      if [ -n "$c" ]; then
        g=$(echo "$c" | xargs -r -n1 pgrep -P 2>/dev/null || true)
      fi
      ps -o rss= -p "$p" $c $g 2>/dev/null | awk "{s+=\$1} END{print s+0}"
    }

    while kill -0 "$pid" 2>/dev/null; do
      cur_main=$(get_main_rss_kb "$pid")
      cur_total=$(get_total_rss_kb "$pid")
      if [ "$cur_main"  -gt "$peak_main"  ]; then peak_main="$cur_main"; fi
      if [ "$cur_total" -gt "$peak_total" ]; then peak_total="$cur_total"; fi
      sleep "'"$RSS_SAMPLE_INTERVAL"'"
    done

    set +e
    wait "$pid"
    rc=$?
    set -e

    echo "__PEAK_RSS_TOTAL_KB__=${peak_total}" 1>&2
    echo "__PEAK_RSS_MAIN_KB__=${peak_main}"   1>&2
    echo "__INNER_RC__=${rc}"                  1>&2

    exit 0
  ' _ "$@" \
    > >(tee "$out") \
    2> >(tee "$err" >&2)

  local wrapper_rc=${PIPESTATUS[0]}
  set -e

  local inner_rc
  inner_rc=$(parse_inner_rc "$err")
  if [ -z "$inner_rc" ]; then
    inner_rc="$wrapper_rc"
  fi
  echo "$inner_rc" > "$codefile"
}

run_algo() {
  local algo="$1" bin="$2" ref="$3" qfile="$4" qct="$5" errc="$6" threads_in="$7" rep="$8"

  # Only naive_search uses multi-threading; all others forced to 1
  local threads="$threads_in"
  if [[ "$algo" != "naive_search" ]]; then
    threads=1
  fi

  local qlen qbase tag out errf timef codefile code
  qlen="$(infer_query_len "$qfile")"
  qbase="$(basename "$qfile")"
  tag="$(sanitize "${algo}_$(basename "$ref")_${qbase}_ct${qct}_e${errc}_t${threads}_r${rep}")"
  out="build/logs/${tag}.out"
  errf="build/logs/${tag}.err"
  timef="build/logs/${tag}.time"
  codefile="build/logs/${tag}.code"

  echo
  echo "----------------------------------------------------------"
  echo "[RUN] $algo | ref=$(basename "$ref") | q=$qbase | len=$qlen | ct=$qct | e=$errc | t=$threads | rep=$rep"
  echo "  stdout: $out"
  echo "  stderr: $errf"
  echo "  time  : $timef"

  case "$algo" in
    naive_search)
      run_cmd_timed "$out" "$errf" "$timef" "$codefile" "$bin" \
        --reference "$ref" \
        --query "$qfile" \
        --query_ct "$qct" \
        --threads "$threads" \
        --min_block "$MIN_BLOCK"
      ;;
    suffixarray_search)
      run_cmd_timed "$out" "$errf" "$timef" "$codefile" "$bin" \
        --reference "$ref" \
        --query "$qfile" \
        --query_ct "$qct"
      ;;
    fmindex_search|fmindex_pigeon_search)
      run_cmd_timed "$out" "$errf" "$timef" "$codefile" "$bin" \
        --reference "$ref" \
        --query "$qfile" \
        --query_ct "$qct" \
        --errors "$errc" \
        --threads "$threads"
      ;;
    *)
      echo "Unknown algo: $algo" >&2
      return 1
      ;;
  esac

  code="$(cat "$codefile")"

  local idx_s search_s wall_s rss_kb peak_total_kb peak_main_kb
  idx_s=$(extract_time_from_logs "Index Construction time" "$out" "$errf")
  search_s=$(extract_time_from_logs "Search time" "$out" "$errf")

  wall_s=$(to_seconds "$(parse_wall_str_timefile "$timef")")
  rss_kb=$(parse_max_rss_kb_timefile "$timef")

  peak_total_kb=$(parse_peak_total_rss_kb "$errf")
  peak_main_kb=$(parse_peak_main_rss_kb "$errf")

  [ -z "$search_s" ] && search_s="$wall_s"

  echo "[DONE] exit=$code | index_s=${idx_s:-NA} | search_s=${search_s:-NA} | wall_s=${wall_s:-NA} | max_rss_kb=${rss_kb:-NA} | peak_total_kb=${peak_total_kb:-NA} | peak_main_kb=${peak_main_kb:-NA}"
  echo "${algo},${ref},${qbase},${qlen},${qct},${errc},${threads},${rep},${idx_s:-},${search_s:-},${wall_s:-},${rss_kb:-},${peak_total_kb:-},${peak_main_kb:-},${code}" >> "$CSV_FILE"
}

echo "--- Task 1: varying query counts ---"
for ct in "${COUNTS[@]}"; do
  for r in $(seq 1 "$REPEATS"); do
    run_algo naive_search          "$BIN_NAIVE"  "$REFERENCE_SMALL" "$QFILE_100" "$ct" 0 "$THREADS" "$r"
    run_algo suffixarray_search    "$BIN_SA"     "$REFERENCE_SMALL" "$QFILE_100" "$ct" 0 "$THREADS" "$r"
    run_algo fmindex_search        "$BIN_FM"     "$REFERENCE_SMALL" "$QFILE_100" "$ct" 0 "$THREADS" "$r"
    run_algo fmindex_pigeon_search "$BIN_PIGEON" "$REFERENCE_SMALL" "$QFILE_100" "$ct" 0 "$THREADS" "$r"
  done
done

echo "--- Task 2: varying query length (full genome) ---"
for len in "${LENGTHS[@]}"; do
  Q_FILE="data/illumina_reads_${len}.fasta.gz"
  [ -f "$Q_FILE" ] || { echo "Missing $Q_FILE, skipping."; continue; }
  for r in $(seq 1 "$REPEATS"); do
    (( RUN_NAIVE_FULL )) && run_algo naive_search "$BIN_NAIVE" "$REFERENCE_FULL" "$Q_FILE" "$FIXED_COUNT" 0 "$THREADS" "$r"
    (( RUN_SA_FULL ))    && run_algo suffixarray_search "$BIN_SA" "$REFERENCE_FULL" "$Q_FILE" "$FIXED_COUNT" 0 "$THREADS" "$r"
    run_algo fmindex_search        "$BIN_FM"     "$REFERENCE_FULL" "$Q_FILE" "$FIXED_COUNT" 0 "$THREADS" "$r"
    run_algo fmindex_pigeon_search "$BIN_PIGEON" "$REFERENCE_FULL" "$Q_FILE" "$FIXED_COUNT" 0 "$THREADS" "$r"
  done
done

echo "Done. CSV written to $CSV_FILE"
echo "Logs: build/logs/"
