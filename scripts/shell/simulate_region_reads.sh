#!/usr/bin/env bash
# Simulate error-free paired-end reads for one region from the GRCh38 reference
# FASTA (spec sections 43-45), rename reads, and build the truth table.
#
# Usage:
#   simulate_region_reads.sh <chrom_fasta> <chrom> <start_0based> <end> <region_id> <out_dir> <truth_tsv> \
#       [n_pairs] [read_len] [frag_mean] [frag_sd] [seed]
set -euo pipefail

source "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
conda activate bh26

CHROM_FASTA="$1"
CHROM="$2"
START0="$3"
END="$4"
REGION_ID="$5"
OUT_DIR="$6"
TRUTH_TSV="$7"
N_PAIRS="${8:-30000}"
READ_LEN="${9:-150}"
FRAG_MEAN="${10:-350}"
FRAG_SD="${11:-35}"
SEED="${12:-42}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${OUT_DIR}"

LOG="${OUT_DIR}/${REGION_ID}_simulate.log"
exec > >(tee "${LOG}") 2>&1

echo "=== simulate_region_reads.sh ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "region_id=${REGION_ID} chrom=${CHROM} start0=${START0} end=${END} n_pairs=${N_PAIRS} read_len=${READ_LEN} frag=${FRAG_MEAN}+-${FRAG_SD} seed=${SEED}"

START1=$((START0 + 1))
REGION_FASTA="${OUT_DIR}/${REGION_ID}_reference.fa"
samtools faidx "${CHROM_FASTA}" "${CHROM}:${START1}-${END}" > "${REGION_FASTA}"

RAW_R1="${OUT_DIR}/${REGION_ID}_raw_R1.fastq"
RAW_R2="${OUT_DIR}/${REGION_ID}_raw_R2.fastq"

# Error-free, mutation-free, indel-free primary simulation (spec section 43/46).
wgsim -1 "${READ_LEN}" -2 "${READ_LEN}" -d "${FRAG_MEAN}" -s "${FRAG_SD}" \
      -N "${N_PAIRS}" -e 0 -r 0 -R 0 -S "${SEED}" \
      "${REGION_FASTA}" "${RAW_R1}" "${RAW_R2}"

FINAL_R1="${OUT_DIR}/${REGION_ID}_R1.fastq.gz"
FINAL_R2="${OUT_DIR}/${REGION_ID}_R2.fastq.gz"

python "${SCRIPT_DIR}/../python/process_wgsim_output.py" \
    --r1-in "${RAW_R1}" --r2-in "${RAW_R2}" \
    --r1-out "${FINAL_R1}" --r2-out "${FINAL_R2}" \
    --region-id "${REGION_ID}" --chrom "${CHROM}" \
    --window-start-0based "${START0}" \
    --truth-out "${TRUTH_TSV}"

rm -f "${RAW_R1}" "${RAW_R2}" "${REGION_FASTA}.fai"

echo "=== done ==="
ls -la "${FINAL_R1}" "${FINAL_R2}"
