#!/usr/bin/env bash
# Orchestrates the full single-region 100kb pipeline: extract -> giraffe index
# -> wgsim read simulation -> vg giraffe mapping -> vg surject, writing clear
# stage-completion markers to the log so a wrapping watchdog (or a human) can
# tell exactly how far it got. Intended to be launched inside a tmux session
# together with an external wall-clock watchdog (see README "Technical
# limitations" for why: vg giraffe/surject have shown pathological
# multi-minute-per-read hangs on specific reads in specific local subgraphs,
# unrelated to true memory pressure, that a plain foreground run would block
# on indefinitely).
#
# Usage:
#   run_full_pipeline_100kb.sh <region_id> <chrom> <start_0based> <end> <run_dir> [threads] [merge_distance]
set -uo pipefail

source "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
conda activate bh26

REGION_ID="$1"
CHROM="$2"
START="$3"
END="$4"
RUN_DIR="$5"
THREADS="${6:-12}"
MERGE_DISTANCE="${7:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CHR_OG="${PROJECT_ROOT}/data/graphs/chr21.og"
CHR_FASTA="${PROJECT_ROOT}/data/reference/chr21.fa"
TRUTH_TSV="${RUN_DIR}/metadata/truth.tsv"

echo "=== run_full_pipeline_100kb.sh ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "region_id=${REGION_ID} range=${CHROM}:${START}-${END} threads=${THREADS} merge_distance=${MERGE_DISTANCE:-<default>}"

echo "--- stage: extract ---"
if [ -n "${MERGE_DISTANCE}" ]; then
    bash "${SCRIPT_DIR}/extract_region_final.sh" "${CHR_OG}" "GRCh38#0#${CHROM}" "${START}" "${END}" "${REGION_ID}" "${RUN_DIR}/work/extracted_final" "${THREADS}" "${MERGE_DISTANCE}"
else
    bash "${SCRIPT_DIR}/extract_region_final.sh" "${CHR_OG}" "GRCh38#0#${CHROM}" "${START}" "${END}" "${REGION_ID}" "${RUN_DIR}/work/extracted_final" "${THREADS}"
fi
EXTRACT_EXIT=$?
if [ "${EXTRACT_EXIT}" -ne 0 ]; then
    echo "EXTRACT_EXIT_CODE=${EXTRACT_EXIT}"
    exit "${EXTRACT_EXIT}"
fi
echo "EXTRACT_DONE"

GFA="${RUN_DIR}/work/extracted_final/${REGION_ID}/${REGION_ID}.gfa"

echo "--- stage: index ---"
bash "${SCRIPT_DIR}/build_giraffe_index.sh" "${GFA}" "${REGION_ID}" "${RUN_DIR}/work/indexes" "${THREADS}"
INDEX_EXIT=$?
echo "INDEX_EXIT_CODE=${INDEX_EXIT}"
if [ "${INDEX_EXIT}" -ne 0 ]; then
    exit "${INDEX_EXIT}"
fi
echo "INDEX_DONE"

GBZ="${RUN_DIR}/work/indexes/${REGION_ID}/${REGION_ID}.giraffe.gbz"

echo "--- stage: simulate reads ---"
bash "${SCRIPT_DIR}/simulate_region_reads.sh" "${CHR_FASTA}" "${CHROM}" "${START}" "${END}" "${REGION_ID}" "${RUN_DIR}/work/reads" "${TRUTH_TSV}"
SIM_EXIT=$?
if [ "${SIM_EXIT}" -ne 0 ]; then
    echo "SIM_EXIT_CODE=${SIM_EXIT}"
    exit "${SIM_EXIT}"
fi
echo "SIM_DONE"

R1="${RUN_DIR}/work/reads/${REGION_ID}_R1.fastq.gz"
R2="${RUN_DIR}/work/reads/${REGION_ID}_R2.fastq.gz"

echo "--- stage: map + surject ---"
bash "${SCRIPT_DIR}/run_giraffe_mapping.sh" "${GBZ}" "${R1}" "${R2}" "GRCh38#0#${CHROM}" "${REGION_ID}" "${RUN_DIR}/work/mapping" "${THREADS}"
MAP_EXIT=$?
echo "MAP_EXIT_CODE=${MAP_EXIT}"
if [ "${MAP_EXIT}" -eq 0 ]; then
    echo "MAP_DONE"
fi

echo "=== pipeline finished for ${REGION_ID} at $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
exit "${MAP_EXIT}"
