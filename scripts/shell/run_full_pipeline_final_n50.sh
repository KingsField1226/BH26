#!/usr/bin/env bash
# The single, clean, one-pass pipeline for the final n=50 stratified analysis.
# Unlike the earlier ad-hoc n=20 -> n=30 -> n=50 growth (kept only in
# logs/hazard_stratified_robustness_analysis.md), this runs all 50
# pre-selected regions end to end in one batch, including mapping evaluation
# (the earlier run_full_pipeline_batch_safe.sh stopped after surject and left
# evaluate_mapping.py to a separate, unscripted step).
#
# Learns from every prior incident in this project:
#   - a hard per-region wall-clock kill (background PID + timeout loop, since
#     macOS lacks GNU timeout), at 30 minutes (not the original 10-20 minute
#     watchdogs, which produced 2 false failures during the earlier growth).
#   - automatic detection + retry of the "default subpath merging balloons
#     the window past its requested size" pathology (recurred 4 times
#     independently across this project) -- merge_distance=0 retry if the
#     first extraction exceeds 3x the requested window size.
#   - every region's outcome is recorded in results/n50_final_progress.tsv
#     (ephemeral -- deleted once combined_results_n50_final.tsv is built)
#     rather than silently skipped.
#
# Usage:
#   run_full_pipeline_final_n50.sh <manifest.tsv> <run_dir> [threads] [per_region_timeout_sec]
#
# <manifest.tsv> needs columns: region_id, chrom, start, end (extra columns ignored).
set -uo pipefail

source "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
conda activate bh26

MANIFEST="$1"
RUN_DIR="$2"
THREADS="${3:-12}"
TIMEOUT_SEC="${4:-1800}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${RUN_DIR}/../.." && pwd)"
WORK_DIR="${RUN_DIR}/work"
TRUTH_TSV="${RUN_DIR}/metadata/truth_n50_final.tsv"
MAPPING_STATS="${RUN_DIR}/results/mapping_stats_n50_final.tsv"
PROGRESS="${RUN_DIR}/results/n50_final_progress.tsv"

mkdir -p "${WORK_DIR}/extracted_final" "${WORK_DIR}/indexes" "${WORK_DIR}/reads" "${WORK_DIR}/mapping" \
    "${RUN_DIR}/logs/n50_final" "${RUN_DIR}/results"

if [ ! -f "${PROGRESS}" ]; then
    printf "region_id\tstart\tend\textract_bp\tblowup_fixed\tindex_exit\tmap_exit\tsurject_exit\teval_exit\tstatus\telapsed_sec\n" > "${PROGRESS}"
fi

run_one_region() {
    local region_id="$1" chrom="$2" start="$3" end="$4"
    local region_log="${RUN_DIR}/logs/n50_final/${region_id}.log"
    local extract_dir="${WORK_DIR}/extracted_final/${region_id}"
    local og="${extract_dir}/${region_id}.og"
    local gfa="${extract_dir}/${region_id}.gfa"
    local window_size=$((end - start))
    local blowup_fixed="False"

    {
        echo "=== n50_final: ${region_id} (${chrom}:${start}-${end}) ==="
        date -u +%Y-%m-%dT%H:%M:%SZ

        echo "--- stage: extract (default merge) ---"
        bash "${SCRIPT_DIR}/extract_region_final.sh" "${PROJECT_ROOT}/data/graphs/chr21.og" "GRCh38#0#${chrom}" "${start}" "${end}" "${region_id}" "${WORK_DIR}/extracted_final" "${THREADS}"
        extract_exit=$?
        if [ "${extract_exit}" -ne 0 ]; then
            echo "EXTRACT_EXIT_CODE=${extract_exit}"
            echo "STATUS=extract_failed"
            exit 1
        fi

        extract_bp=$(awk 'NR==2{print $1}' "${extract_dir}/odgi_stats.tsv" 2>/dev/null)
        echo "extract_bp=${extract_bp} window_size=${window_size}"
        if [ -n "${extract_bp}" ] && [ "${extract_bp}" -gt $((window_size * 3)) ]; then
            echo "BLOWUP DETECTED (${extract_bp} > 3x ${window_size}) -- re-extracting with merge_distance=0"
            rm -rf "${extract_dir}"
            bash "${SCRIPT_DIR}/extract_region_final.sh" "${PROJECT_ROOT}/data/graphs/chr21.og" "GRCh38#0#${chrom}" "${start}" "${end}" "${region_id}" "${WORK_DIR}/extracted_final" "${THREADS}" 0
            extract_exit=$?
            blowup_fixed="True"
            extract_bp=$(awk 'NR==2{print $1}' "${extract_dir}/odgi_stats.tsv" 2>/dev/null)
            echo "corrected extract_bp=${extract_bp}"
            if [ "${extract_exit}" -ne 0 ]; then
                echo "EXTRACT_EXIT_CODE=${extract_exit} (after blowup retry)"
                echo "STATUS=extract_failed_after_retry"
                exit 1
            fi
        fi
        echo "EXTRACT_DONE"

        echo "--- stage: index ---"
        bash "${SCRIPT_DIR}/build_giraffe_index.sh" "${gfa}" "${region_id}" "${WORK_DIR}/indexes" "${THREADS}"
        index_exit=$?
        echo "INDEX_EXIT_CODE=${index_exit}"
        if [ "${index_exit}" -ne 0 ]; then
            echo "STATUS=index_failed"
            exit 1
        fi
        echo "INDEX_DONE"

        gbz="${WORK_DIR}/indexes/${region_id}/${region_id}.giraffe.gbz"

        echo "--- stage: simulate ---"
        bash "${SCRIPT_DIR}/simulate_region_reads.sh" "${PROJECT_ROOT}/data/reference/chr21.fa" "${chrom}" "${start}" "${end}" "${region_id}" "${WORK_DIR}/reads" "${TRUTH_TSV}"
        sim_exit=$?
        if [ "${sim_exit}" -ne 0 ]; then
            echo "SIM_EXIT_CODE=${sim_exit}"
            echo "STATUS=sim_failed"
            exit 1
        fi
        echo "SIM_DONE"

        r1="${WORK_DIR}/reads/${region_id}_R1.fastq.gz"
        r2="${WORK_DIR}/reads/${region_id}_R2.fastq.gz"

        echo "--- stage: map + surject ---"
        bash "${SCRIPT_DIR}/run_giraffe_mapping.sh" "${gbz}" "${r1}" "${r2}" "GRCh38#0#${chrom}" "${region_id}" "${WORK_DIR}/mapping" "${THREADS}"
        map_exit=$?
        echo "MAP_EXIT_CODE=${map_exit}"
        if [ "${map_exit}" -ne 0 ]; then
            echo "STATUS=map_or_surject_failed"
            exit 1
        fi
        echo "MAP_DONE"

        bam="${WORK_DIR}/mapping/${region_id}/${region_id}.surjected.bam"

        echo "--- stage: evaluate ---"
        python "${SCRIPT_DIR}/../python/evaluate_mapping.py" \
            --bam "${bam}" --truth "${TRUTH_TSV}" \
            --region-id "${region_id}" --window-size "${window_size}" \
            --out "${MAPPING_STATS}"
        eval_exit=$?
        echo "EVAL_EXIT_CODE=${eval_exit}"
        if [ "${eval_exit}" -ne 0 ]; then
            echo "STATUS=eval_failed"
            exit 1
        fi
        echo "EVAL_DONE"
        echo "STATUS=success"
    } > "${region_log}" 2>&1
    return 0
}

TOTAL=$(($(wc -l < "${MANIFEST}") - 1))
N=0

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r region_id chrom start end _rest; do
    N=$((N + 1))
    echo "=== [${N}/${TOTAL}] ${region_id} (${chrom}:${start}-${end}) at $(date -u +%H:%M:%S) ==="

    t0=$(date +%s)
    run_one_region "${region_id}" "${chrom}" "${start}" "${end}" &
    child_pid=$!

    elapsed=0
    interval=15
    while kill -0 "${child_pid}" 2>/dev/null; do
        sleep "${interval}"
        elapsed=$((elapsed + interval))
        if [ "${elapsed}" -ge "${TIMEOUT_SEC}" ]; then
            echo "TIMEOUT: ${region_id} exceeded ${TIMEOUT_SEC}s, killing" >&2
            kill -9 "${child_pid}" 2>/dev/null
            pkill -9 -f "${region_id}" 2>/dev/null
            break
        fi
    done
    wait "${child_pid}" 2>/dev/null
    t1=$(date +%s)
    elapsed_total=$((t1 - t0))

    region_log="${RUN_DIR}/logs/n50_final/${region_id}.log"
    status=$(grep -oE "STATUS=[a-z_]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    status="${status:-timeout_or_crash}"
    extract_bp=$(grep -oE "^extract_bp=[0-9]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    blowup_fixed=$(grep -q "BLOWUP DETECTED" "${region_log}" 2>/dev/null && echo "True" || echo "False")
    index_exit=$(grep -oE "^INDEX_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    map_exit=$(grep -oE "^MAP_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    surject_exit=$(grep -oE "^SURJECT_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    eval_exit=$(grep -oE "^EVAL_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${region_id}" "${start}" "${end}" "${extract_bp:-NA}" "${blowup_fixed}" \
        "${index_exit:-NA}" "${map_exit:-NA}" "${surject_exit:-NA}" "${eval_exit:-NA}" "${status}" "${elapsed_total}" \
        >> "${PROGRESS}"

    echo "  -> status=${status} elapsed=${elapsed_total}s"
done

echo "=== n50_final batch complete ==="
