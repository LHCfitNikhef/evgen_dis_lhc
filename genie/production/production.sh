#!/usr/bin/env bash
# THE  GENIE PRODUCTION (pp03/pp04): one lane per configuration.
#
# Usage:  genie/production/production.sh <cfg> [<cfg> ...]
#         (each cfg runs as its own lane; lanes run in parallel, and inside a
#          lane the (target, energy) points run one after another, each with
#          its jobs in parallel)
#
# STATISTICS (user, 2026-09-14: "if statistics are an issue, feel free to
# increase sample size (for 1 TeV of course)").  1 TeV carries the
# differential figure (pp04), so it is the high-statistics point; the other
# energies only need sigma_fid (pp03).  Events per NUCLEON:
#
#   G18_02a (mu_grv, mu_nnpdf, nu_grv, nu_nnpdf):  1 TeV 500k (8 jobs),
#                                                  400 / 700 GeV 100k (4 jobs)
#   HEDIS (nu_hedis):                              1 TeV 400k (8 jobs),
#                                                  400/700/2000/4000 100k (4 jobs)
#
# the earlier production's 1 TeV samples were 200k (G18_02a) and 300k (HEDIS) on the proton alone.
# G18_02a stops at 1 TeV (its declared validity; genie_job.sh refuses above).
# Measured throughput per core at 1 TeV: G18_02a mu ~285 ev/s, nu ~670 ev/s,
# HEDIS ~30 ev/s -- HEDIS is the long lane.
#
# Every point is skipped once its jobs carry V2_OK, so a lane can be killed
# and re-invoked (CONVENTIONS.md 1c).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
LOG=$BENCH_REPO/logs/genie_job
mkdir -p "$LOG"

lane () {
    local c=$1 e t
    case "$c" in
        nu_hedis) points="1000:50000:8 400:25000:4 700:25000:4 2000:25000:4 4000:25000:4" ;;
        *)        points="1000:62500:8 400:25000:4 700:25000:4" ;;
    esac
    for pt in $points; do
        e=${pt%%:*}; rest=${pt#*:}
        # HEDIS runs ~30 events/s per core, so its two nucleons run side by
        # side; the fast G18_02a lanes keep them in sequence
        for t in p n; do
            ( "$HERE/genie_job.sh" run "$c" "$t" "$e" "${rest%%:*}" "${rest#*:}" \
                >> "$LOG/${c}_${t}.log" 2>&1 \
                || echo "  FAILED $c $t $e GeV -- see logs/genie_job/${c}_${t}.log" ) &
            [ "$c" = nu_hedis ] || wait
        done
        wait
        echo "  $(date +%H:%M) $c $e GeV done"
    done
}

for c in "$@"; do lane "$c" & done
wait
echo "== GENIE lanes finished: $* =="
