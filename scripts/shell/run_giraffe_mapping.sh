#!/usr/bin/env bash
# Map simulated paired-end reads with vg giraffe, then surject to GRCh38.
# Spec sections 47-49: keep indexing and mapping resource benchmarks separate.
#
# NOTE: <path_prefix> is a PREFIX (e.g. "GRCh38#0#chr21"), not necessarily the
# exact full path name -- odgi extract renames a reference path that only
# partially covers the original chromosome path into a PanSN "subpath" name
# with a coordinate suffix (e.g. "GRCh38#0#chr21:38134180-38234292"), and that
# exact suffix differs per region since extraction snaps to node boundaries.
# The actual path name is looked up dynamically from the GBZ via `vg paths -L`.
#
# Usage:
#   run_giraffe_mapping.sh <gbz> <r1.fastq.gz> <r2.fastq.gz> <path_prefix> <region_id> <out_dir> [threads]
set -uo pipefail

VG_BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/software/vg/bin/vg"

GBZ="$1"
R1="$2"
R2="$3"
PATH_PREFIX="$4"
REGION_ID="$5"
OUT_DIR="$6"
THREADS="${7:-12}"

MAP_DIR="${OUT_DIR}/${REGION_ID}"
mkdir -p "${MAP_DIR}"

LOG="${MAP_DIR}/mapping.log"
exec > >(tee "${LOG}") 2>&1

echo "=== run_giraffe_mapping.sh ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "region_id=${REGION_ID} gbz=${GBZ} r1=${R1} r2=${R2} path_prefix=${PATH_PREFIX} threads=${THREADS}"

PATH_NAME=$("${VG_BIN}" paths -x "${GBZ}" -L | grep "^${PATH_PREFIX}" | head -1)
if [ -z "${PATH_NAME}" ]; then
    echo "ERROR: no path in ${GBZ} matches prefix '${PATH_PREFIX}'"
    "${VG_BIN}" paths -x "${GBZ}" -L
    exit 1
fi
echo "resolved_path_name=${PATH_NAME}"

GAM_OUT="${MAP_DIR}/${REGION_ID}.gam"

echo "--- vg giraffe ---"
/usr/bin/time -l "${VG_BIN}" giraffe \
    -Z "${GBZ}" \
    -f "${R1}" -f "${R2}" \
    -t "${THREADS}" \
    -p \
    > "${GAM_OUT}"
GIRAFFE_EXIT=$?
echo "GIRAFFE_EXIT_CODE=${GIRAFFE_EXIT}"

if [ "${GIRAFFE_EXIT}" -ne 0 ]; then
    echo "=== giraffe failed, skipping surject ==="
    exit "${GIRAFFE_EXIT}"
fi

BAM_OUT="${MAP_DIR}/${REGION_ID}.surjected.bam"

echo "--- vg surject ---"
/usr/bin/time -l "${VG_BIN}" surject \
    -x "${GBZ}" \
    -D short \
    -p "${PATH_NAME}" \
    -i \
    -b \
    -t "${THREADS}" \
    "${GAM_OUT}" \
    > "${BAM_OUT}"
SURJECT_EXIT=$?
echo "SURJECT_EXIT_CODE=${SURJECT_EXIT}"

echo "=== output files ==="
ls -la "${MAP_DIR}/"

echo "=== done ==="
date -u +%Y-%m-%dT%H:%M:%SZ
exit "${SURJECT_EXIT}"
