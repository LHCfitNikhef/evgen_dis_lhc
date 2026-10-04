#!/usr/bin/env bash
# The POWHEG-V2mc lane: 1 TeV first (the high-statistics point), then the
# other energies, proton then neutron, V2_PAR batches at a time.
# Usage: powheg/production/production_mc.sh [max_parallel=2]
# Re-invocable: run_v2.sh skips generated batches and verified showers.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
PAR=${1:-2}
LOG=$BENCH_REPO/logs
# keep the machine awake for the lane (config.sh decides; empty on Linux)
case "$BENCH_NOSLEEP" in
  caffeinate*) $BENCH_NOSLEEP -s -w $$ > /dev/null 2>&1 & ;;
esac
for e in 1000 400 700 2000 4000; do
    for t in p n; do
        echo "[$(date '+%F %T')] V2mc $e $t: start"
        if ! V2_PAR=$PAR "$HERE/run_v2_mc.sh" "$e" "$t" \
                >> "$LOG/powheg_v2mc_${e}_${t}.log" 2>&1; then
            echo "[$(date '+%F %T')] V2mc $e $t: FAILED -- logs/powheg_v2mc_${e}_${t}.log; lane stopped"
            exit 1
        fi
        echo "[$(date '+%F %T')] V2mc $e $t: done"
    done
done
echo "[$(date '+%F %T')] V2mc lane: ALL DONE"
