#!/usr/bin/env bash
# Parallel Pythia8 production for the CC neutrino-DIS benchmark: N seeded jobs.
# Usage: ./run_pythia_nu.sh <events_per_job> <njobs> [beam_energy_gev]
#
# The muon-side twin is run_pythia.sh; the two are kept in step deliberately
# (CONVENTIONS.md rule 2b).  The anchor keeps dis_nu1TeV.cmnd and nu_job_N; other
# energies use dis_nu<tag>.cmnd and nu_job_<tag>_N, both named by beams.py.
set -e
NEV=${1:-50000}
NJOBS=${2:-8}
EBEAM=${3:-1000}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag nu "$EBEAM")
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name nu_job "$EBEAM")
CARD=$HERE/dis_nu$TAG.cmnd
[ -f "$CARD" ] || {
    echo "no card $CARD -- run tools/make_energy_cards.py first" >&2; exit 1; }

if [ "$EBEAM" = "1000" ]; then SEED0=7321; else SEED0=$((7321 + EBEAM)); fi

for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "Pythia8 nu: $NJOBS x $NEV events, nu_mu $EBEAM GeV [$TAG] -> ${JOBBASE}_N"
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  mkdir -p "$d"
  (cd "$d" && bench_run "$HERE/main_dis" "$CARD" events \
      "$NEV" $((SEED0 + i)) > pythia.log 2>&1) &
done
wait
echo "all $NJOBS Pythia neutrino jobs done"
grep -h "sigma_gen_mb" "$HERE"/${JOBBASE}_*/events_xsec.json
