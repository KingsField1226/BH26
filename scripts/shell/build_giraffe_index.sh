#!/usr/bin/env bash
# Build a Giraffe (GBZ + distance + minimizer/zipcode) index for one region's GFA.
# Spec sections 39-42: record wall time, peak RSS, exit code, index file sizes.
#
# Usage: build_giraffe_index.sh <region.gfa> <region_id> <out_dir> [threads] [target_mem]
set -uo pipefail  # NOT -e: we want to capture and record autoindex failures, not abort the script

VG_BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/software/vg/bin/vg"

REGION_GFA="$1"
REGION_ID="$2"
OUT_DIR="$3"
THREADS="${4:-12}"
TARGET_MEM="${5:-85G}"

INDEX_DIR="${OUT_DIR}/${REGION_ID}"
mkdir -p "${INDEX_DIR}"

LOG="${INDEX_DIR}/index.log"
exec > >(tee "${LOG}") 2>&1

echo "=== build_giraffe_index.sh ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "region_id=${REGION_ID} gfa=${REGION_GFA} threads=${THREADS} target_mem=${TARGET_MEM}"
echo "vg: $("${VG_BIN}" version | head -1)"

/usr/bin/time -l "${VG_BIN}" autoindex \
    --workflow sr-giraffe \
    --gfa "${REGION_GFA}" \
    --prefix "${INDEX_DIR}/${REGION_ID}" \
    --target-mem "${TARGET_MEM}" \
    --threads "${THREADS}" \
    --verbosity 2
EXIT_CODE=$?

echo "INDEX_EXIT_CODE=${EXIT_CODE}"
echo "=== index files ==="
ls -la "${INDEX_DIR}/" || true

echo "=== done ==="
date -u +%Y-%m-%dT%H:%M:%SZ
exit "${EXIT_CODE}"
