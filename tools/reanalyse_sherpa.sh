#!/usr/bin/env bash
# Re-derive EVERY Sherpa result file after the shower-charm subtraction
# (analyze.SHOWER_CHARM_SUBTRACTED_PREFIXES, user decision 2026-09-03).
#
# Usage: tools/reanalyse_sherpa.sh [--dry-run] [--only <substring>]   (MAXJOBS=8)
#
# The job list is DERIVED FROM THE RESULT FILES ON DISK, not typed: every
# results{,_nu}/histos_sherpa*.json names a (key, selection, energy, Q2 floor)
# that was published, and each is regenerated through the same command that
# made it.  Hard-process-tag diagnostics (`_charm`, `_charmany`) and the
# shower-less ME-level samples (`_me`) are untouched by the subtraction and
# are skipped.  After it: tools/make_all_plots.sh,
# tools/sync_paper_figures.sh, tools/check_paper_claims.py.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
MAXJOBS=${MAXJOBS:-8}
DRY=0; ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in --dry-run) DRY=1; shift ;; --only) ONLY=$2; shift 2 ;; *) shift ;; esac
done
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
n=0
for cur in mu nu; do
  case $cur in mu) d=results;    script=analyze.py    ;;
               nu) d=results_nu; script=analyze_nu.py ;; esac
  for f in "$BENCH_REPO"/$d/histos_sherpa*.json; do
    name=$(basename "$f" .json); name=${name#histos_}
    case "$name" in *_charm|*_charmany|*_me|*_me_*) continue ;; esac
    [ -n "$ONLY" ] && case "$name" in *"$ONLY"*) ;; *) continue ;; esac
    q2=4;  case "$name" in *_q2min11) q2=11; name=${name%_q2min11} ;; esac
    e=1000; case "$name" in *_400GeV) e=400; name=${name%_400GeV} ;; *_4TeV) e=4000; name=${name%_4TeV} ;; esac
    sel=inclusive
    case "$name" in *_faser_dimuon) sel=faser_dimuon; name=${name%_faser_dimuon} ;;
                    *_faser_e)      sel=faser_e;      name=${name%_faser_e} ;;
                    *_faser_s)      sel=faser_s;      name=${name%_faser_s} ;; esac
    n=$((n+1))
    log=$BENCH_REPO/logs/reanalyse_sherpa/${cur}_${name}_${sel}_E${e}_q$q2.log
    mkdir -p "$(dirname "$log")"
    echo "[$n] $cur $name sel=$sel E=$e Q2min=$q2"
    [ "$DRY" = 1 ] && continue
    limit
    ( cd "$BENCH_REPO" && BENCH_SELECTION=$sel BENCH_ENERGY=$e BENCH_Q2MIN=$q2 \
        python3 analysis/$script "$name" > "$log" 2>&1 \
      && echo "  done: $cur $name $sel E=$e q$q2" \
      || echo "  FAILED: $cur $name $sel E=$e q$q2 -- see $log" ) &
  done
done
wait
echo "reanalysed $n Sherpa result(s)"
