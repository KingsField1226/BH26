#!/usr/bin/env bash
# Local-execution environment verification, per analysis spec section 4.1.
# Run this at the start of any analysis session and before long-running jobs.
set -euo pipefail

echo "=== pwd ==="
pwd

echo "=== uname -m ==="
uname -m

echo "=== sw_vers ==="
sw_vers

echo "=== CONDA_DEFAULT_ENV ==="
echo "${CONDA_DEFAULT_ENV:-<not set>}"

echo "=== which python ==="
which python

echo "=== python --version ==="
python --version
