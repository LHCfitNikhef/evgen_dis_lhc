#!/usr/bin/env bash
# Parallel fixed-order NLO (parton-level) event production: N independent
# Sherpa jobs with distinct random seeds sharing the cached integration grids
# (each gets a private copy of Results.zip so nothing writes to the shared one).
#
# The output MUST be HepMC3_Short: only the "short" writer serialises the
# correlated NLO sub-event lists of the real-emission (RS) part
# (SHERPA/Tools/HepMC3_Interface.C: SubEvtList2ShortHepMC).  With the full
# HepMC3_GenEvent writer the sub-events are silently dropped and the RS
# contribution would be misrepresented.
#
# Usage: ./run_parallel.sh <events_per_job> <njobs> [extra Sherpa args...]
set -e

# This script also runs from the Sherpa install tree, where it is a symlink
# back into the repo (tools/link_cards.sh).  Resolve the link to find
# config.sh; $0's own directory is the install tree, not the repo.
SELF=$0
if [ -L "$SELF" ]; then SELF=$(readlink "$SELF"); fi
BENCH_REPO=${BENCH_REPO:-$(cd "$(dirname "$SELF")/../../.." && pwd)}
. "$BENCH_REPO/config.sh"          # bench_run
NEV=${1:-25000}
NJOBS=${2:-8}
shift 2 2>/dev/null || true
for i in $(seq 1 $NJOBS); do
  d=job_$i
  mkdir -p $d
  ln -sf ../Sherpa.yaml $d/
  [ -e $d/Process ] || ln -s ../Process $d/Process
  cp -f Results.zip $d/
  (cd $d && bench_run ../../../install/bin/Sherpa -e $NEV \
     "RANDOM_SEED: $((7100 + i))" \
     'EVENT_OUTPUT: HepMC3_Short[evtme]' "$@" \
     > sherpa.log 2>&1) &
done
wait
echo "all $NJOBS jobs done"
grep -h "Generated events" job_*/sherpa.log
