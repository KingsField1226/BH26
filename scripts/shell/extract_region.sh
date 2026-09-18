#!/usr/bin/env bash
# NOT USED in the final pipeline -- kept for the record. This uses
# `odgi extract -E/--full-range`, which was found to be impractical on the
# chr21 ODGI graph built for this project (see README.md "Technical
# limitations" and logs/vg_build_summary.txt-adjacent notes): it expands to
# nearly the whole chromosome graph regardless of locus. Use
# extract_region_screen.sh (fast candidate screening) or
# extract_region_final.sh (final representative regions, default subpath
# merging, no -E) instead.
#
# Extract a local subgraph from a chromosome-level ODGI graph for a given
# GRCh38 path range, and compute its graph complexity metrics.
#
# Usage:
#   extract_region.sh <input.og> <path_name> <start_0based> <end> <region_id> <out_dir> [threads]
#
# Outputs (under <out_dir>/<region_id>/):
#   <region_id>.og
#   <region_id>.gfa
#   odgi_stats.tsv
#   odgi_degree.tsv
#   extract.log
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

echo "=== extract_region.sh ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "input_og: ${INPUT_OG}"
echo "path_name: ${PATH_NAME}"
echo "range: ${START}-${END} (0-based)"
echo "region_id: ${REGION_ID}"
echo "threads: ${THREADS}"

OG_OUT="${REGION_DIR}/${REGION_ID}.og"
GFA_OUT="${REGION_DIR}/${REGION_ID}.gfa"

echo "--- odgi extract ---"
/usr/bin/time -l odgi extract \
    -i "${INPUT_OG}" \
    -r "${PATH_NAME}:${START}-${END}" \
    -E \
    -O \
    -t "${THREADS}" \
    -P \
    -o "${OG_OUT}"

echo "--- odgi view (GFA export) ---"
odgi view -i "${OG_OUT}" -g -t "${THREADS}" > "${GFA_OUT}"

echo "--- odgi stats -S ---"
odgi stats -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_stats.tsv"

echo "--- odgi degree -S ---"
odgi degree -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_degree.tsv"

echo "=== done ==="
