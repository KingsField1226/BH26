#!/usr/bin/env bash
# Fast candidate-window extraction for complexity screening (spec sections 31-36).
#
# Unlike extract_region.sh (used for the final representative regions before
# Giraffe indexing), this does NOT use `odgi extract -E/--full-range`. -E
# collects all nodes in the sorted min/max span touched by the path range,
# which is correct/important for building a clean, un-fragmented GFA for
# indexing, but on highly repetitive/complex loci that span can balloon far
# beyond the requested window and make extraction extremely slow (observed:
# a single 100kb candidate window near the chr21 centromere took ~2 hours
# with -E). For screening/ranking candidate windows by complexity, the
# unlaced extraction is a reasonable, much faster proxy -- odgi stats/degree
# on it still reflect real local topology. The 4 finally-selected
# representative regions are re-extracted properly (with -E) by
# extract_region.sh before indexing.
#
# Usage:
#   extract_region_screen.sh <input.og> <path_name> <start_0based> <end> <region_id> <out_dir> [threads]
set -euo pipefail

INPUT_OG="$1"
PATH_NAME="$2"
START="$3"
END="$4"
REGION_ID="$5"
OUT_ROOT="$6"
THREADS="${7:-12}"

REGION_DIR="${OUT_ROOT}/${REGION_ID}"
mkdir -p "${REGION_DIR}"

LOG="${REGION_DIR}/extract.log"
exec > >(tee "${LOG}") 2>&1

echo "=== extract_region_screen.sh (fast, no -E) ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "input_og: ${INPUT_OG}"
echo "path_name: ${PATH_NAME}"
echo "range: ${START}-${END} (0-based)"
echo "region_id: ${REGION_ID}"
echo "threads: ${THREADS}"

OG_OUT="${REGION_DIR}/${REGION_ID}.og"

echo "--- odgi extract (screening, no -E, no subpath merging) ---"
/usr/bin/time -l odgi extract \
    -i "${INPUT_OG}" \
    -r "${PATH_NAME}:${START}-${END}" \
    -d 0 \
    -O \
    -t "${THREADS}" \
    -P \
    -o "${OG_OUT}"

echo "--- odgi stats -S ---"
odgi stats -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_stats.tsv"

echo "--- odgi degree -S ---"
odgi degree -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_degree.tsv"

echo "=== done ==="
