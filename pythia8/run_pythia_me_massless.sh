#!/usr/bin/env bash
# MASSLESS-QUARK DIAGNOSTIC production (ME level, dis_mu1TeV_me_massless.cmnd):
# same as run_pythia_me.sh but with all five quark masses set to zero, so that
# Pythia's phase space matches the benchmark's ZM (massless-charm) convention
# and the strictly massless YADISM reference. Output in me_ml_job_*/.
# NOT a benchmark sample -- a diagnostic; see the header of the card.
# Usage: ./run_pythia_me_massless.sh <events_per_job> <njobs>
set -e
NEV=${1:-400000}
NJOBS=${2:-3}
HERE=$(cd "$(dirname "$0")" && pwd)
# the shell profile points PYTHIA8DATA at conda's 8.312 xmldoc; use our 8.311
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo
for i in $(seq 1 $NJOBS); do
  d=$HERE/me_ml_job_$i
  mkdir -p $d
  (cd $d && bench_run $HERE/main_dis $HERE/dis_mu1TeV_me_massless.cmnd \
      events $NEV $((4451 + i)) > pythia.log 2>&1) &
done
wait
echo "all $NJOBS Pythia massless-ME jobs done"
grep -h "sigma_gen_mb" $HERE/me_ml_job_*/events_xsec.json
