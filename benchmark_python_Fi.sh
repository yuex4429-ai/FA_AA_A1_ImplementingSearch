#!/bin/bash
set -euo pipefail

# ==========================================================
# Python Benchmarks (Assignment 2 compliant)
# - Four modes: naive, suffixarray, fmindex, fmindex_pigeon
# - Captures:
#     * "Index Construction time: X seconds." (if printed)
#     * "Search time: X seconds." (if printed)
#     * /usr/bin/time -v written to per-run .time file (clean logs)
#     * PEAK Total RSS (python pid + worker descendants) via sampling
# - Writes numeric seconds CSV.
# ==========================================================

REFERENCE_SMALL="data/hg38_partial.fasta.gz"
REFERENCE_FULL="data/hg38_partial.fasta.gz"
# REFERENCE_FULL="data/GCF_000001405.26_GRCh38_genomic.fna"
# REFERENCE_SMALL="data/GCF_000001405.26_GRCh38_genomic.fasta.gz"
# REFERENCE_FULL="data/GCF_000001405.26_GRCh38_genomic.fasta.gz"

QFILE_100="data/illumina_reads_100.fasta.gz"

COUNTS=(1000 10000 100000 1000000)      
LENGTHS=(40 60 80 100)        
FIXED_COUNT=1000

ERRORS=(0)

PROCS=14
REPEATS=1

PY="python3"
BIN="src/implementing_search.py"
CSV="results/benchmark_python_all.csv"

MODES=(naive suffixarray fmindex)
# MODES=(naive suffixarray fmindex fmindex_pigeon)

SAMPLING_RATE=16
MIN_BLOCK=256
RSS_SAMPLE_INTERVAL="0.2"

mkdir -p results build/logs

echo "Mode,Reference,Query_File,Query_Len,Query_Count,Errors,Procs,Repeat,Index_Time_s,Search_Time_s,Wall_Time_s,Max_RSS_KB,Peak_Total_RSS_KB,Peak_Main_RSS_KB,Exit_Code" > "$CSV"

sanitize() { echo "$1" | tr '/: ' '___' | tr -cd 'A-Za-z0-9_.-'; }
infer_len() { [[ "$1" =~ _([0-9]+)\.fasta(\.gz)?$ ]] && echo "${BASH_REMATCH[1]}" || echo "N/A"; }

extract_seconds_value() { sed -n 's/.*: \([0-9]\+\(\.[0-9]\+\)\?\) seconds\.*/\1/p' | tail -n 1; }
extract_time_from_logs() {
  local pattern="$1" out="$2" err="$3"
  local v=""
  if grep -q "$pattern" "$out" 2>/dev/null; then v=$(grep "$pattern" "$out" | extract_seconds_value || true); fi
  if [ -z "$v" ] && grep -q "$pattern" "$err" 2>/dev/null; then v=$(grep "$pattern" "$err" | extract_seconds_value || true); fi
  echo "$v"
}

# Parse from /usr/bin/time -v output file (NOT from .err, so logs are clean)
parse_wall_str_timefile() { grep "Elapsed (wall clock) time" "$1" | awk '{print $NF}' | tail -n 1 || echo ""; }
parse_max_rss_kb_timefile() { grep "Maximum resident set size" "$1" | awk '{print $NF}' | tail -n 1 || echo ""; }

# Markers: match only at line start
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

run_cmd_timed() {
  local out="$1" err="$2" timef="$3" codefile="$4"; shift 4

  # time -v output goes to $timef (clean), wrapper stderr (markers) goes to $err
  set +e
  /usr/bin/time -v -o "$timef" bash -c '
    set -euo pipefail

    "$@" &
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

run_one() {
  local mode="$1" ref="$2" qfile="$3" qct="$4" errc="$5" procs="$6" rep="$7"
  local qlen qbase tag out errf timef codefile code
  qlen="$(infer_len "$qfile")"
  qbase="$(basename "$qfile")"

  if [ "$mode" = "naive" ] || [ "$mode" = "suffixarray" ]; then
    errc=0
  fi

  tag="$(sanitize "${mode}_$(basename "$ref")_${qbase}_ct${qct}_e${errc}_p${procs}_r${rep}")"
  out="build/logs/${tag}.out"
  errf="build/logs/${tag}.err"
  timef="build/logs/${tag}.time"
  codefile="build/logs/${tag}.code"

  echo
  echo "----------------------------------------------------------"
  echo "[RUN] mode=$mode | ref=$(basename "$ref") | q=$qbase | len=$qlen | ct=$qct | e=$errc | procs=$procs | rep=$rep"
  echo "  stdout: $out"
  echo "  stderr: $errf"
  echo "  time  : $timef"

  args=(
    "$PY" "$BIN"
    --mode "$mode"
    --reference "$ref"
    --query "$qfile"
    --query_ct "$qct"
    --errors "$errc"
    --threads "$procs"
    --min_block "$MIN_BLOCK"
  )
  if [ "$mode" = "fmindex" ] || [ "$mode" = "fmindex_pigeon" ]; then
    args+=( --sampling_rate "$SAMPLING_RATE" )
  fi

  run_cmd_timed "$out" "$errf" "$timef" "$codefile" "${args[@]}"
  code="$(cat "$codefile")"

  local idx_s search_s wall_s rss_kb peak_total_kb peak_main_kb
  idx_s=$(extract_time_from_logs "Index Construction time" "$out" "$errf")
  search_s=$(extract_time_from_logs "Search time" "$out" "$errf")

  wall_s=$(to_seconds "$(parse_wall_str_timefile "$timef")")
  rss_kb=$(parse_max_rss_kb_timefile "$timef")

  peak_total_kb=$(parse_peak_total_rss_kb "$errf")
  peak_main_kb=$(parse_peak_main_rss_kb "$errf")

  if [ -z "$search_s" ]; then search_s="$wall_s"; fi

  echo "[DONE] exit=$code | index_s=${idx_s:-NA} | search_s=${search_s:-NA} | wall_s=${wall_s:-NA} | max_rss_kb=${rss_kb:-NA} | peak_total_kb=${peak_total_kb:-NA} | peak_main_kb=${peak_main_kb:-NA}"
  echo "${mode},${ref},${qbase},${qlen},${qct},${errc},${procs},${rep},${idx_s:-},${search_s:-},${wall_s:-},${rss_kb:-},${peak_total_kb:-},${peak_main_kb:-},${code}" >> "$CSV"
}

# ---- sanity checks ----
[ -f "$REFERENCE_SMALL" ] || { echo "Error: missing $REFERENCE_SMALL" >&2; exit 1; }
[ -f "$REFERENCE_FULL" ]  || { echo "Error: missing $REFERENCE_FULL"  >&2; exit 1; }
[ -f "$QFILE_100" ]       || { echo "Error: missing $QFILE_100"       >&2; exit 1; }

task1_runs=$(( ${#MODES[@]} * ${#COUNTS[@]} * ${#ERRORS[@]} * REPEATS ))

exist_len=0
for len in "${LENGTHS[@]}"; do
  q="data/illumina_reads_${len}.fasta.gz"
  if [ -f "$q" ]; then exist_len=$((exist_len+1)); fi
done
task2_runs=$(( ${#MODES[@]} * exist_len * ${#ERRORS[@]} * REPEATS ))
total_runs=$(( task1_runs + task2_runs ))

echo "=========================================================="
echo "Python Benchmarks (Assignment 2 compliant)"
echo "Task1 ref: $REFERENCE_SMALL"
echo "Task2 ref: $REFERENCE_FULL"
echo "CSV: $CSV"
echo "BIN: $BIN"
echo "MODES=${MODES[*]}"
echo "PROCS=$PROCS REPEATS=$REPEATS ERRORS=${ERRORS[*]}"
echo "FM sampling_rate=$SAMPLING_RATE"
echo "ProcessPool MIN_BLOCK=$MIN_BLOCK"
echo "RSS_SAMPLE_INTERVAL=$RSS_SAMPLE_INTERVAL"
echo "Task1 counts: ${COUNTS[*]} (len=100)"
echo "Task2 lengths: ${LENGTHS[*]} (ct=$FIXED_COUNT)"
echo "Planned runs: TOTAL_RUNS=$total_runs (Task1=$task1_runs, Task2=$task2_runs)"
echo "=========================================================="

echo
echo "--- Task 1: varying query counts (len=100) on hg38_partial ---"
for e in "${ERRORS[@]}"; do
  for ct in "${COUNTS[@]}"; do
    for r in $(seq 1 "$REPEATS"); do
      for mode in "${MODES[@]}"; do
        run_one "$mode" "$REFERENCE_SMALL" "$QFILE_100" "$ct" "$e" "$PROCS" "$r"
      done
    done
  done
done

echo
echo "--- Task 2: varying query lengths (ct=$FIXED_COUNT) on full genome ---"
for e in "${ERRORS[@]}"; do
  for len in "${LENGTHS[@]}"; do
    q="data/illumina_reads_${len}.fasta.gz"
    if [ ! -f "$q" ]; then
      echo "Warning: missing $q, skipping."
      continue
    fi
    for r in $(seq 1 "$REPEATS"); do
      for mode in "${MODES[@]}"; do
        run_one "$mode" "$REFERENCE_FULL" "$q" "$FIXED_COUNT" "$e" "$PROCS" "$r"
      done
    done
  done
done

echo "=========================================================="
echo "Done. CSV: $CSV"
echo "Logs: build/logs/"
echo "=========================================================="
