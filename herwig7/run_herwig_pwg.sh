#!/usr/bin/env bash
# Parallel Herwig 7.3.0 POWHEG-matched NLO production, either current, either
# half of the NLO weight.
#
# Usage: ./run_herwig_pwg.sh <mu|nu> <pos|neg> <nev_per_job> <njobs> [beam_gev]
#
# BOTH HALVES ARE REQUIRED FOR A PHYSICAL RESULT.  Herwig's DISBase splits the
# NLO weight into its positive and negative parts and generates them as two
# separate POSITIVE-weight samples, selected by the Contribution switch
# (DISBase.cc::NLOWeight returns max(0,w) for 1 and max(0,-w) for 2).  Herwig's
# defaults put Contribution 1 on PowhegMEDIS, so a single run out of the box is
# only HALF the answer and is too LARGE -- it looks perfectly reasonable, which
# is why this is trap 5 in the Herwig notes.  herwig7/combine_nlo.py assembles
# pos - neg.
#
# Replaces run_herwig_mu_pwg.sh and run_herwig_nu_pwg.sh, which differed only
# by the current, were hardcoded to 1 TeV job directories, had no negative-half
# runner at all, and would happily overwrite the anchor's samples.
set -e
CUR=${1:-mu}
HALF=${2:-pos}
NEV=${3:-25000}
NJOBS=${4:-6}
EBEAM=${5:-1000}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
HW=$HERWIG_BIN

case "$CUR" in mu|nu) ;; *) echo "current must be mu or nu" >&2; exit 1;; esac
# STEM is the card's `saverun` name and hence the .run file; ANCHORCARD is the
# 1 TeV card's FILENAME, which is not the same thing for the positive half --
# DIS-mu-POWHEG.in declares `saverun DIS-mu-PWG`.  The generated per-energy
# cards are named after the saverun stem, so they need no such exception.
case "$HALF" in
  pos) STEM="DIS-$CUR-PWG";    ANCHORCARD="DIS-$CUR-POWHEG.in"
       JOBSTEM="${CUR}pwg_job";;
  neg) STEM="DIS-$CUR-PWGNEG"; ANCHORCARD="DIS-$CUR-PWGNEG.in"
       JOBSTEM="${CUR}pwgneg_job";;
  *) echo "half must be pos or neg" >&2; exit 1;;
esac

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$CUR" "$EBEAM")
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name "$JOBSTEM" "$EBEAM")
# the anchor keeps its original card and job names, so re-running it cannot
# move a published number
if [ "$EBEAM" = "1000" ]; then RUNSTEM=$STEM; CARD=$ANCHORCARD
else RUNSTEM="$STEM-$TAG"; CARD=$RUNSTEM.in; fi
RUNF=$RUNSTEM.run
[ -f "$HERE/$CARD" ] || {
    echo "no card $HERE/$CARD -- run tools/make_energy_cards.py first" >&2
    exit 1; }

hw_env() { bench_run env -i HOME="$HOME" \
    PATH="$HW:$BENCH_EXTRA_BIN:/usr/bin:/bin" \
    LHAPDF_DATA_PATH="$LHAPDF_DATA_PATH" "$@"; }

if [ ! -f "$HERE/$RUNF" ] || [ "$HERE/$CARD" -nt "$HERE/$RUNF" ]; then
    echo "compiling $CARD -> $RUNF"
    (cd "$HERE" && hw_env Herwig read "$CARD" > "read_${RUNSTEM}.log" 2>&1) || {
        echo "Herwig read failed; see herwig7/read_${RUNSTEM}.log" >&2; exit 1; }
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

echo "Herwig POWHEG $CUR $HALF: $NJOBS x $NEV events, $EBEAM GeV [$TAG] -> ${JOBBASE}_N"
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  mkdir -p "$d"
  cp "$HERE/$RUNF" "$d/"
  [ -f "$HERE/LeptonicDISCut.so" ] && cp "$HERE/LeptonicDISCut.so" "$d/"
  [ -d "$HERE/Herwig-cache/$RUNSTEM" ] && cp -R "$HERE/Herwig-cache" "$d/"
  (cd "$d" && hw_env Herwig run "$RUNF" -N "$NEV" -s $((7500 + i)) -d 0 \
      > herwig_run.log 2>&1) &
done
wait
echo "all $NJOBS Herwig POWHEG $CUR $HALF jobs done"
for i in $(seq 1 "$NJOBS"); do
  grep -A7 "Statistics for event handler" "$HERE/${JOBBASE}_$i"/${RUNSTEM}*.out | tail -3
done
