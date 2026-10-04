#!/usr/bin/env bash
# Matrix-element-level Pythia8 production (PartonLevel:all = off): N seeded
# jobs with the dis_mu1TeV_me.cmnd card, output in me_job_*/.
# Usage: ./run_pythia_me.sh <events_per_job> <njobs>
set -e
NEV=${1:-50000}
NJOBS=${2:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
# the shell profile points PYTHIA8DATA at conda's 8.312 xmldoc; use our 8.311
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo
for i in $(seq 1 $NJOBS); do
  d=$HERE/me_job_$i
  mkdir -p $d
  (cd $d && bench_run $HERE/main_dis $HERE/dis_mu1TeV_me.cmnd events \
      $NEV $((7321 + i)) > pythia.log 2>&1) &
done
wait
echo "all $NJOBS Pythia ME jobs done"
grep -h "sigma_gen_mb" $HERE/me_job_*/events_xsec.json
