#!/usr/bin/env bash
# Monitor system memory while a named tmux session runs; log periodic snapshots,
# and if resident memory of the tracked process pattern exceeds the stop
# threshold (spec section 63: RSS > ~100-105 GB and still climbing), kill the
# tmux session and record why.
#
# Usage: monitor_resources.sh <tmux_session_name> <proc_pattern> <log_file> [interval_sec] [rss_limit_gb]
set -uo pipefail

SESSION="$1"
PATTERN="$2"
LOG="$3"
INTERVAL="${4:-30}"
RSS_LIMIT_GB="${5:-100}"

echo "=== monitor_resources.sh started $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" >> "${LOG}"
echo "tracking tmux session '${SESSION}', process pattern '${PATTERN}', RSS stop limit ${RSS_LIMIT_GB}GB" >> "${LOG}"

LOW_FREE_STREAK=0

while tmux has-session -t "${SESSION}" 2>/dev/null; do
    # Sum RSS (KB, from ps -o rss=) across all matching processes.
    RSS_KB=$(ps -axo rss=,command= | grep -i "${PATTERN}" | grep -v grep | awk '{sum+=$1} END {print sum+0}')
    RSS_GB=$(( RSS_KB / 1024 / 1024 ))
    TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    PHYSMEM=$(top -l 1 -n 0 | grep "PhysMem")
    # Free-percentage safety net (spec section 63: critical macOS memory pressure).
    # top's "unused" figure (MB) is parsed out; sustained (3 consecutive checks)
    # near-zero unused memory, even if the tracked process's own RSS looks fine,
    # indicates system-wide pressure worth stopping for (this is how a prior
    # `odgi extract` run crashed with a Bus error without ever tripping the
    # per-process RSS limit below).
    UNUSED_MB=$(echo "${PHYSMEM}" | grep -oE '[0-9]+M unused' | grep -oE '^[0-9]+')
    UNUSED_MB="${UNUSED_MB:-99999}"
    echo "${TS} tracked_rss_gb=${RSS_GB} unused_mb=${UNUSED_MB} ${PHYSMEM}" >> "${LOG}"

    if [ "${RSS_GB}" -ge "${RSS_LIMIT_GB}" ]; then
        echo "${TS} STOP CONDITION: tracked RSS ${RSS_GB}GB >= limit ${RSS_LIMIT_GB}GB. Killing tmux session ${SESSION}." >> "${LOG}"
        tmux kill-session -t "${SESSION}" 2>/dev/null
        break
    fi

    if [ "${UNUSED_MB}" -lt 1000 ]; then
        LOW_FREE_STREAK=$((LOW_FREE_STREAK + 1))
    else
        LOW_FREE_STREAK=0
    fi
    if [ "${LOW_FREE_STREAK}" -ge 5 ]; then
        echo "${TS} STOP CONDITION: system unused memory < 1GB for ${LOW_FREE_STREAK} consecutive checks. Killing tmux session ${SESSION}." >> "${LOG}"
        tmux kill-session -t "${SESSION}" 2>/dev/null
        break
    fi
    sleep "${INTERVAL}"
done

echo "=== monitor_resources.sh exiting $(date -u +%Y-%m-%dT%H:%M:%SZ) (session ${SESSION} no longer running) ===" >> "${LOG}"
