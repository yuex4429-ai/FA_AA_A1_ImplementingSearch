#!/usr/bin/env bash
set -euo pipefail

LOGDIR="build/logs"
OUT="results/benchmark_suffixarray.csv"

mkdir -p results
echo "program,query_len,query_ct,search_time_s,wall_time,max_rss_kb,hits,log_base,cmd" > "$OUT"

for ERR in "$LOGDIR"/*.err; do
    BASE=$(basename "$ERR" .err)
    OUTFILE="$LOGDIR/$BASE.out"

    # safer command extraction (ignore leading spaces)
    CMD_LINE=$(grep -m1 -E "Command being timed:" "$ERR" || true)
    CMD=$(echo "$CMD_LINE" | sed 's/.*"//;s/"//')

    # fallback: decide by filename if grep failed
    if [[ -z "$CMD" ]]; then
        if [[ "$BASE" == suffixarray_* ]]; then
            CMD="(from filename) $BASE"
        else
            continue
        fi
    fi

    # only suffixarray
    if [[ "$CMD" != *suffixarray* && "$BASE" != suffixarray_* ]]; then
        continue
    fi

    # program
    if [[ "$CMD" == *construct* || "$BASE" == *construct* ]]; then
        PROGRAM="suffixarray_construct"
    else
        PROGRAM="suffixarray_search"
    fi

    # params from CMD
    QUERY_CT=$(echo "$CMD" | grep -oP -- '--query_ct\s+\K[0-9]+' || echo "")
    QUERY_FILE=$(echo "$CMD" | grep -oP -- '--query\s+\K\S+' || echo "")
    QUERY_LEN=$(echo "$QUERY_FILE" | grep -oP '[0-9]+(?=\.fasta)' || echo "")

    # stdout
    SEARCH_TIME=""
    HITS=""
    if [[ -f "$OUTFILE" ]]; then
        SEARCH_TIME=$(grep -oP 'Search time:\s+\K[0-9\.]+' "$OUTFILE" || echo "")
        HITS=$(grep -oP 'hits=\K[0-9]+' "$OUTFILE" || echo "")
    fi

    # stderr
    WALL_TIME=$(grep -oP 'Elapsed .*:\s+\K.*' "$ERR" || echo "")
    MAX_RSS=$(grep -oP 'Maximum resident set size \(kbytes\):\s+\K[0-9]+' "$ERR" || echo "")

    echo "$PROGRAM,$QUERY_LEN,$QUERY_CT,$SEARCH_TIME,$WALL_TIME,$MAX_RSS,$HITS,$BASE,\"$CMD\"" >> "$OUT"
done

echo "Written to $OUT"

