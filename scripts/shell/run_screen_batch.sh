#!/usr/bin/env bash
# Run extract_region_screen.sh (fast, no -E) for every row in a candidate-region manifest.
#
# Usage: run_screen_batch.sh <manifest.tsv> <input.og> <path_name> <out_dir> [threads]
set -euo pipefail

MANIFEST="$1"
INPUT_OG="$2"
PATH_NAME="$3"
OUT_DIR="$4"
THREADS="${5:-12}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${OUT_DIR}"

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

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r -a ROW; do
    region_id="${ROW[$region_id_idx]}"
    start="${ROW[$start_idx]}"
    end="${ROW[$end_idx]}"
    echo "=== screening ${region_id} (${PATH_NAME}:${start}-${end}) ==="
    bash "${SCRIPT_DIR}/extract_region_screen.sh" "${INPUT_OG}" "${PATH_NAME}" "${start}" "${end}" "${region_id}" "${OUT_DIR}" "${THREADS}"
done

echo "=== screening batch complete ==="
