#!/usr/bin/env bash
# Parallel Herwig 7.3.0 LO production for the CC neutrino-DIS benchmark.
# Usage: ./run_herwig_nu.sh <events_per_job> <njobs> [beam_energy_gev]
#
# The muon-side twin is run_herwig.sh; kept in step deliberately (rule 2b).
# The anchor keeps DIS-nu.in and nu_job_N; other energies use DIS-nu-<tag>.in
# and nu_job_<tag>_N.  The .run is COMPILED here when missing or stale.
#
# NOTE the run name: Herwig names its output after the card's `saverun`, so
# DIS-nu.in produces DIS-nu.run and DIS-nu-400GeV.in produces
# DIS-nu-400GeV.run.  Leaving saverun alone in a generated card would have it
# overwrite the 1 TeV run file -- which is exactly what happened on the muon
# side before tools/make_energy_cards.py learned to rename it.
set -e
NEV=${1:-25000}
NJOBS=${2:-4}
EBEAM=${3:-1000}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
HW=$HERWIG_BIN

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag nu "$EBEAM")
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name nu_job "$EBEAM")
# THE POLARIZED CARD, at every energy.  DIS-nu.in leaves the neutrino
# beam unpolarised and Herwig then spin-averages over two helicities
# when a massless neutrino has one, HALVING the cross-section.  The
# published 1 TeV sample is DIS-nu-POL, which is why
# HERWIG_NU_SPIN_FACTOR is 1.
if [ "$EBEAM" = "1000" ]; then CARD=DIS-nu-POL.in; else CARD=DIS-nu-POL-$TAG.in; fi
RUNF=${CARD%.in}.run
[ -f "$HERE/$CARD" ] || {
    echo "no card $HERE/$CARD -- run tools/make_energy_cards.py first" >&2
    exit 1; }

hw_env() { bench_run env -i HOME="$HOME" \
    PATH="$HW:$BENCH_EXTRA_BIN:/usr/bin:/bin" \
    LHAPDF_DATA_PATH="$LHAPDF_DATA_PATH" "$@"; }

if [ ! -f "$HERE/$RUNF" ] || [ "$HERE/$CARD" -nt "$HERE/$RUNF" ]; then
    echo "compiling $CARD -> $RUNF"
    (cd "$HERE" && hw_env Herwig read "$CARD" > "read_nu_$TAG.log" 2>&1) || {
        echo "Herwig read failed; see herwig7/read_nu_$TAG.log" >&2; exit 1; }
    [ -f "$HERE/$RUNF" ] || {
        echo "Herwig read produced no $RUNF -- check the card's saverun" >&2
        exit 1; }
fi

for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "Herwig nu: $NJOBS x $NEV events, nu_mu $EBEAM GeV [$TAG] -> ${JOBBASE}_N"
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  mkdir -p "$d"
  cp "$HERE/$RUNF" "$d/"
  (cd "$d" && hw_env Herwig run "$RUNF" -N "$NEV" -s $((7500 + i)) -d 0 \
      > herwig_run.log 2>&1) &
done
wait
echo "all $NJOBS Herwig neutrino jobs done"
for i in $(seq 1 "$NJOBS"); do
  grep -A7 "Statistics for event handler" "$HERE/${JOBBASE}_$i"/*-S*.out | tail -3
done
