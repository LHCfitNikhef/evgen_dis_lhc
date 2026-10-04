#!/usr/bin/env bash
# The POWHEG part of the paper-plots production (PRODUCTION.md): both
# codes, both nucleons, all five energies.
#
# Usage: nohup powheg/production/production.sh [all|res|] > logs/powheg_v2_production.log 2>&1 &
#
# TWO LANES, 4 CORES IN TOTAL (other productions share the machine):
#   RES lane  POWHEG-RES mu NC, 3 processes at a time  (run_res.sh, V2_PAR=3)
#   V2  lane  POWHEG-V2  nu CC, 1 process at a time    (run_v2.sh,  V2_PAR=1)
# Each lane walks (energy, nucleon) in the order below and STOPS at the first
# failure -- including a disk-floor refusal -- rather than skipping ahead.
# Both drivers skip finished work, so after fixing the cause the same command
# picks up where the lane stopped.
#
# Per-(energy, nucleon) logs: logs/powheg_v2_<code>_<TAG>_<t>.log.
#
# ORDER: 400 GeV first (a small point, so the first showered job of every
# code x nucleon can be inspected early), then the 1 TeV anchor, then the rest.
#
# Validation (powheg/production/validate.sh) is NOT run from here: it must have been
# run and its numbers recorded before this is started.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"

WHAT=${1:-all}
ENERGIES="400 1000 700 2000 4000"
LOGDIR=$BENCH_REPO/logs
mkdir -p "$LOGDIR"

lane () {   # $1 = res|v2
    local code=$1 e t tag cur drv par
    case "$code" in
        res) cur=mu; drv=$HERE/run_res.sh; par=3 ;;
        v2)  cur=nu; drv=$HERE/run_v2.sh;  par=1 ;;
    esac
    for e in $ENERGIES; do
        for t in p n; do
            tag=$(v2_tag "$cur" "$e")
            echo "[$(date '+%F %T')] lane $code: start $tag $t"
            if ! V2_PAR=$par "$drv" "$e" "$t" \
                    >> "$LOGDIR/powheg_v2_${code}_${tag}_${t}.log" 2>&1; then
                echo "[$(date '+%F %T')] lane $code: FAILED at $tag $t --" \
                     "see $LOGDIR/powheg_v2_${code}_${tag}_${t}.log; lane stopped"
                return 1
            fi
            echo "[$(date '+%F %T')] lane $code: done $tag $t"
        done
    done
    echo "[$(date '+%F %T')] lane $code: ALL DONE"
}

pids=()
case "$WHAT" in
    all) lane res & pids+=($!); lane v2 & pids+=($!) ;;
    res) lane res & pids+=($!) ;;
    v2)  lane v2 & pids+=($!) ;;
    *) echo "usage: production.sh [all|res|v2]" >&2; exit 1 ;;
esac
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ "$fail" = 0 ] && echo "POWHEG V2 PRODUCTION DONE" || { echo "a lane failed" >&2; exit 1; }
