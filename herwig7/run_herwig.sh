#!/usr/bin/env bash
# Parallel Herwig 7.3.0 LO production for the muon-DIS benchmark.
# Usage: ./run_herwig.sh <events_per_job> <njobs> [beam_energy_gev]
#
# The energy defaults to the 1 TeV anchor, which keeps its original card
# (DIS-mu.in) and untagged job dirs so nothing published moves.  Other
# energies use DIS-mu-<tag>.in and job_<tag>_N, both named by
# analysis/beams.py -- the one place that knows the rule.
#
#   ./run_herwig.sh 25000 4          # 1 TeV, the published sample
#   ./run_herwig.sh 25000 8 400      # 400 GeV -> job_400GeV_1..8
#
# The .run file is COMPILED from the .in here if missing or stale, so a new
# energy needs no separate `Herwig read` step.
set -e
NEV=${1:-25000}
NJOBS=${2:-4}
EBEAM=${3:-1000}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo
HW=$HERWIG_BIN

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag mu "$EBEAM")
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name job "$EBEAM")
if [ "$EBEAM" = "1000" ]; then CARD=DIS-mu.in; else CARD=DIS-mu-$TAG.in; fi
RUNF=${CARD%.in}.run
[ -f "$HERE/$CARD" ] || {
    echo "no card $HERE/$CARD -- run tools/make_energy_cards.py first" >&2
    exit 1; }

# Herwig's env: `env -i` on purpose, so a stray LHAPDF/PATH in the shell
# cannot change what the generator sees.
hw_env() { bench_run env -i HOME="$HOME" \
    PATH="$HW:$BENCH_EXTRA_BIN:/usr/bin:/bin" \
    LHAPDF_DATA_PATH="$LHAPDF_DATA_PATH" "$@"; }

# Compile the card if the .run is missing or older than its source -- OR than
# the LeptonicDISCut plugin.  A .run file carries the cut object SERIALISED,
# field by field in the order the class declares them, so rebuilding the
# plugin with a new field and then reusing an old .run reads the new field out
# of the old bytes: no error, a cut that is not the one the card asks for.
# The plugin gained MinW2 on 2026-09-11, which is what made this reachable.
if [ ! -f "$HERE/$RUNF" ] || [ "$HERE/$CARD" -nt "$HERE/$RUNF" ] \
   || { [ -f "$HERE/LeptonicDISCut.so" ] \
        && [ "$HERE/LeptonicDISCut.so" -nt "$HERE/$RUNF" ]; }; then
    echo "compiling $CARD -> $RUNF"
    (cd "$HERE" && hw_env Herwig read "$CARD" > "read_$TAG.log" 2>&1) || {
        echo "Herwig read failed; see herwig7/read_$TAG.log" >&2; exit 1; }
fi

# A job dir holding events is a production sample -- refuse rather than
# silently overwrite (FORCE=1 to override), as run_pythia.sh does.
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "Herwig: $NJOBS x $NEV events, mu- $EBEAM GeV [$TAG] -> ${JOBBASE}_N"
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  mkdir -p "$d"
  cp "$HERE/$RUNF" "$d/"
  (cd "$d" && hw_env Herwig run "$RUNF" -N "$NEV" -s $((7100 + i)) -d 0 \
      > herwig_run.log 2>&1) &
done
wait
echo "all $NJOBS Herwig jobs done"
for i in $(seq 1 "$NJOBS"); do
  grep -A7 "Statistics for event handler" "$HERE/${JOBBASE}_$i"/*-S*.out | tail -3
done
