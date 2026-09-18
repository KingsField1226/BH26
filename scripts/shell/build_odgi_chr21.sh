#!/usr/bin/env bash
# Build the ODGI graph for the chr21 chromosome-level PGGB GFA.
# Run inside tmux per project policy (long-running, may take a while on a 16GB GFA).
set -euo pipefail

source "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
conda activate bh26

cd /Users/kazumichi/Downloads/BioHack26/PGGB_short_reads_mapping

echo "=== build_odgi_chr21.sh ==="
date -u +%Y-%m-%dT%H:%M:%SZ
echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"

/usr/bin/time -l odgi build \
    -g data/graphs/chr21.gfa \
    -o data/graphs/chr21.og \
    -t 12 \
    -P

echo "=== done ==="
date -u +%Y-%m-%dT%H:%M:%SZ
ls -la data/graphs/chr21.og
