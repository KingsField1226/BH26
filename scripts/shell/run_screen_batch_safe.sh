#!/usr/bin/env bash
# Same as run_screen_batch.sh, but with a hard per-window wall-clock timeout.
# One pathological window (e.g. a repeat-dense locus that makes odgi extract
# balloon or hang) must not be allowed to stall the whole scan the way a
# previous vg surject run once did for ~21 hours -- see README.md "Technical
# limitations" for that incident. A timed-out window is recorded as a
# documented failure (odgi_stats.tsv / odgi_degree.tsv simply won't exist for
# it) rather than silently retried or allowed to block everything after it.
#
# Usage: run_screen_batch_safe.sh <manifest.tsv> <input.og> <path_name> <out_dir> [threads] [timeout_sec]
set -uo pipefail

MANIFEST="$1"
INPUT_OG="$2"
PATH_NAME="$3"
OUT_DIR="$4"
THREADS="${5:-12}"
TIMEOUT_SEC="${6:-300}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${OUT_DIR}"

TIMEOUT_LOG="${OUT_DIR}/../timed_out_windows.tsv"
if [ ! -f "${TIMEOUT_LOG}" ]; then
    printf "region_id\tstart\tend\ttimeout_sec\n" > "${TIMEOUT_LOG}"
fi

HEADER=$(head -n1 "${MANIFEST}")
IFS=$'\t' read -r -a COLS <<< "${HEADER}"
region_id_idx=-1; start_idx=-1; end_idx=-1
for i in "${!COLS[@]}"; do
    case "${COLS[$i]}" in
        region_id) region_id_idx=$i ;;
        start) start_idx=$i ;;
        end) end_idx=$i ;;
    esac
done
if [ "${region_id_idx}" -lt 0 ] || [ "${start_idx}" -lt 0 ] || [ "${end_idx}" -lt 0 ]; then
    echo "ERROR: manifest must have region_id, start, end columns" >&2
    exit 1
fi

TOTAL=$(($(wc -l < "${MANIFEST}") - 1))
N=0

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r -a ROW; do
    N=$((N + 1))
    region_id="${ROW[$region_id_idx]}"
    start="${ROW[$start_idx]}"
    end="${ROW[$end_idx]}"
    echo "=== [${N}/${TOTAL}] screening ${region_id} (${PATH_NAME}:${start}-${end}) at $(date -u +%H:%M:%S) ==="

    bash "${SCRIPT_DIR}/extract_region_screen.sh" "${INPUT_OG}" "${PATH_NAME}" "${start}" "${end}" "${region_id}" "${OUT_DIR}" "${THREADS}" &
    child_pid=$!

    elapsed=0
    interval=5
    while kill -0 "${child_pid}" 2>/dev/null; do
        sleep "${interval}"
        elapsed=$((elapsed + interval))
        if [ "${elapsed}" -ge "${TIMEOUT_SEC}" ]; then
            echo "TIMEOUT: ${region_id} exceeded ${TIMEOUT_SEC}s, killing" >&2
            # Kill the wrapper script and any odgi process working on this region's output path.
            kill -9 "${child_pid}" 2>/dev/null
            pkill -9 -f "${OUT_DIR}/${region_id}/${region_id}.og" 2>/dev/null
            pkill -9 -f "odgi extract.*${start}-${end}" 2>/dev/null
            printf "%s\t%s\t%s\t%s\n" "${region_id}" "${start}" "${end}" "${TIMEOUT_SEC}" >> "${TIMEOUT_LOG}"
            break
        fi
    done
    wait "${child_pid}" 2>/dev/null
done

echo "=== safe screening batch complete ==="
