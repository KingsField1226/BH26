#!/usr/bin/env bash
# Final regional extraction for Giraffe indexing.
#
# NOTE ON METHOD CHANGE: odgi extract's -E/--full-range flag (used in an
# earlier version of this script and in the now-unused extract_region.sh)
# was found to be impractical on this particular chr21 ODGI graph: on all
# tested loci (both a near-centromere window and four mid-chromosome
# representative windows spanning 16-42 Mb), -E expanded the "sorted order
# of the graph" span to cover nearly the *entire* chromosome graph (~19 GB
# .og / ~10 GB GFA for a requested 100 kb window), and once actually crashed
# with a Bus error under sustained near-total memory pressure. This is
# consistent with odgi's own documentation warning ("Be careful to use it
# with very complex graphs") and likely reflects node IDs in this graph not
# being strictly coordinate-sorted (chr21.og was built without odgi build -s).
# This script extracts WITHOUT -E. It DOES use odgi extract's default
# subpath-merging distance (-d, default 300000bp) -- unlike the fast
# candidate-screening script (extract_region_screen.sh), which disables
# merging entirely (-d 0) for speed across many candidates. For only 4 final
# regions the extra ~1-4 minutes/region for merging is affordable, and it
# matters: with -d 0, a visual check (odgi viz) showed most non-GRCh38
# haplotype paths were chopped into many small disconnected fragments rather
# than one contiguous local path per haplotype, which is a poor basis for
# both the graph-topology figure and Giraffe indexing. With default merging,
# haplotype paths are stitched back into (close to) one piece each.
# Consequence still documented as a limitation in README.md: haplotype paths
# may be very slightly fragmented at the window boundary rather than fully
# "laced" the way -E would (impractically) attempt -- the region itself is
# still a faithful ~100 kb multi-haplotype subgraph.
#
# Usage:
#   extract_region_final.sh <input.og> <path_name> <start_0based> <end> <region_id> <out_dir> [threads] [merge_distance]
#
# [merge_distance] defaults to odgi's own default (300000bp merging, i.e. the
# flag is simply omitted) if not given. Pass 0 to disable subpath merging
# entirely (-d 0) -- needed for windows that sit at the boundary of a very
# large collapsed repeat structure, where the default merge distance can
# chase fragments across megabases and balloon the "100 kb" region to many
# times its requested size (observed directly: one such window ballooned to
# 4.6 Mb, 46x the requested 100 kb, before this override was added).
set -euo pipefail

INPUT_OG="$1"
PATH_NAME="$2"
START="$3"
END="$4"
REGION_ID="$5"
OUT_ROOT="$6"
THREADS="${7:-12}"
MERGE_DISTANCE="${8:-}"

REGION_DIR="${OUT_ROOT}/${REGION_ID}"
mkdir -p "${REGION_DIR}"

LOG="${REGION_DIR}/extract.log"
exec > >(tee "${LOG}") 2>&1

echo "=== extract_region_final.sh (no -E; see script header) ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "input_og: ${INPUT_OG}"
echo "path_name: ${PATH_NAME}"
echo "range: ${START}-${END} (0-based)"
echo "region_id: ${REGION_ID}"
echo "threads: ${THREADS}"
echo "merge_distance: ${MERGE_DISTANCE:-<default 300000>}"

OG_OUT="${REGION_DIR}/${REGION_ID}.og"
GFA_OUT="${REGION_DIR}/${REGION_ID}.gfa"

echo "--- odgi extract ---"
if [ -n "${MERGE_DISTANCE}" ]; then
    /usr/bin/time -l odgi extract \
        -i "${INPUT_OG}" \
        -r "${PATH_NAME}:${START}-${END}" \
        -d "${MERGE_DISTANCE}" \
        -O \
        -t "${THREADS}" \
        -P \
        -o "${OG_OUT}"
else
    /usr/bin/time -l odgi extract \
        -i "${INPUT_OG}" \
        -r "${PATH_NAME}:${START}-${END}" \
        -O \
        -t "${THREADS}" \
        -P \
        -o "${OG_OUT}"
fi

echo "--- odgi view (GFA export) ---"
odgi view -i "${OG_OUT}" -g -t "${THREADS}" > "${GFA_OUT}"

echo "--- odgi stats -S ---"
odgi stats -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_stats.tsv"

echo "--- odgi degree -S ---"
odgi degree -i "${OG_OUT}" -S -t "${THREADS}" | tee "${REGION_DIR}/odgi_degree.tsv"

echo "=== done ==="
ls -la "${REGION_DIR}/"
