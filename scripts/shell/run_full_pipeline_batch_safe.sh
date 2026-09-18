#!/usr/bin/env bash
# Run the full single-region pipeline (extract -> index -> simulate -> map)
# across MANY regions in one protected batch, for the n=100 complexity-vs-
# outcome scan. Learns from every prior incident in this project:
#   - a hard per-region wall-clock kill (background PID + timeout loop, since
#     macOS lacks GNU timeout) so one pathologically-hanging vg giraffe/surject
#     read (seen repeatedly: implausible multi-TB "memory growth" warnings)
#     cannot stall the whole 100-region batch the way an early unmonitored
#     job once stalled for 21 hours.
#   - automatic detection + retry of the "default subpath merging balloons
#     the window past its requested size" pathology (seen 3 times
#     independently: cand_chr21_010 46x, cand_chr21_354 2.5x, cand_chr21_093
#     5.4x) -- if the first (default-merge) extraction comes out more than
#     3x the requested window size, it is discarded and redone with
#     merge_distance=0.
#   - every region's outcome (success or the stage it failed at) is recorded
#     in a summary TSV rather than silently skipped, so a partial batch is
#     still a usable, honest result.
#
# Usage:
#   run_full_pipeline_batch_safe.sh <manifest.tsv> <run_dir> [threads] [per_region_timeout_sec]
#
# <manifest.tsv> needs columns: region_id, chrom, start, end (extra columns ignored).
set -uo pipefail

source "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
conda activate bh26

MANIFEST="$1"
RUN_DIR="$2"
THREADS="${3:-12}"
TIMEOUT_SEC="${4:-2400}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMMARY="${RUN_DIR}/results/batch100_summary.tsv"
mkdir -p "${RUN_DIR}/logs/batch100"

if [ ! -f "${SUMMARY}" ]; then
    printf "region_id\tstart\tend\textract_bp\tblowup_fixed\tindex_exit\tmap_exit\tsurject_exit\tstatus\telapsed_sec\n" > "${SUMMARY}"
fi

run_one_region() {
    local region_id="$1" chrom="$2" start="$3" end="$4"
    local region_log="${RUN_DIR}/logs/batch100/${region_id}.log"
    local extract_dir="${RUN_DIR}/work/extracted_final/${region_id}"
    local og="${extract_dir}/${region_id}.og"
    local gfa="${extract_dir}/${region_id}.gfa"
    local window_size=$((end - start))
    local blowup_fixed="False"

    {
        echo "=== batch100: ${region_id} (${chrom}:${start}-${end}) ==="
        date -u +%Y-%m-%dT%H:%M:%SZ

        echo "--- stage: extract (default merge) ---"
        bash "${SCRIPT_DIR}/extract_region_final.sh" "${RUN_DIR}/../../data/graphs/chr21.og" "GRCh38#0#${chrom}" "${start}" "${end}" "${region_id}" "${RUN_DIR}/work/extracted_final" "${THREADS}"
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
            bash "${SCRIPT_DIR}/extract_region_final.sh" "${RUN_DIR}/../../data/graphs/chr21.og" "GRCh38#0#${chrom}" "${start}" "${end}" "${region_id}" "${RUN_DIR}/work/extracted_final" "${THREADS}" 0
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
        bash "${SCRIPT_DIR}/build_giraffe_index.sh" "${gfa}" "${region_id}" "${RUN_DIR}/work/indexes" "${THREADS}"
        index_exit=$?
        echo "INDEX_EXIT_CODE=${index_exit}"
        if [ "${index_exit}" -ne 0 ]; then
            echo "STATUS=index_failed"
            exit 1
        fi
        echo "INDEX_DONE"

        gbz="${RUN_DIR}/work/indexes/${region_id}/${region_id}.giraffe.gbz"

        echo "--- stage: simulate ---"
        bash "${SCRIPT_DIR}/simulate_region_reads.sh" "$(cd "${RUN_DIR}/../.." && pwd)/data/reference/chr21.fa" "${chrom}" "${start}" "${end}" "${region_id}" "${RUN_DIR}/work/reads" "${RUN_DIR}/metadata/truth.tsv"
        sim_exit=$?
        if [ "${sim_exit}" -ne 0 ]; then
            echo "SIM_EXIT_CODE=${sim_exit}"
            echo "STATUS=sim_failed"
            exit 1
        fi
        echo "SIM_DONE"

        r1="${RUN_DIR}/work/reads/${region_id}_R1.fastq.gz"
        r2="${RUN_DIR}/work/reads/${region_id}_R2.fastq.gz"

        echo "--- stage: map + surject ---"
        bash "${SCRIPT_DIR}/run_giraffe_mapping.sh" "${gbz}" "${r1}" "${r2}" "GRCh38#0#${chrom}" "${region_id}" "${RUN_DIR}/work/mapping" "${THREADS}"
        map_exit=$?
        echo "MAP_EXIT_CODE=${map_exit}"
        if [ "${map_exit}" -eq 0 ]; then
            echo "MAP_DONE"
            echo "STATUS=success"
        else
            echo "STATUS=map_or_surject_failed"
        fi
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

    region_log="${RUN_DIR}/logs/batch100/${region_id}.log"
    status=$(grep -oE "STATUS=[a-z_]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    status="${status:-timeout_or_crash}"
    extract_bp=$(grep -oE "^extract_bp=[0-9]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    blowup_fixed=$(grep -q "BLOWUP DETECTED" "${region_log}" 2>/dev/null && echo "True" || echo "False")
    index_exit=$(grep -oE "^INDEX_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    map_exit=$(grep -oE "^MAP_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)
    surject_exit=$(grep -oE "^SURJECT_EXIT_CODE=[0-9-]+" "${region_log}" 2>/dev/null | tail -1 | cut -d= -f2)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${region_id}" "${start}" "${end}" "${extract_bp:-NA}" "${blowup_fixed}" \
        "${index_exit:-NA}" "${map_exit:-NA}" "${surject_exit:-NA}" "${status}" "${elapsed_total}" \
        >> "${SUMMARY}"

    echo "  -> status=${status} elapsed=${elapsed_total}s"
done

echo "=== batch100 complete ==="
