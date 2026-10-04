#!/usr/bin/env bash
# Diagnostic: the SAME Sherpa LO card with charm and bottom MASSIVE IN THE
# PARTON SHOWER ONLY (MASSIVE_PS: [4, 5]), the matrix element untouched.
#
# Usage: ./run_sherpa_massive_ps.sh <mu|nu> [lo|nlo] [nev_per_job] [njobs]
#        -> $SHERPA_RUNS/<MuonDIS_LO|NuDIS_LO|MuonDIS_NLO|NuDIS_NLO_ckm3>_massivePS/job_N/evtfull
#
# AT NLO the same override on the MC@NLO card tests the second half of the
# question: whether Sherpa's hard-process charm at NLO, 20.63% of the CC
# rate against 18.87% for POWHEG-V2, is the massless first emission of the
# MC@NLO S-events being counted as hard.
#
# WHY (N4, 2026-09-03).  The benchmark's RESPECT_MASSIVE_FLAG: true makes
# every flavour that is massless in the ME massless in the shower too, so
# Sherpa's CSS splits g -> c cbar with no charm-mass threshold, where Pythia
# and Herwig keep m_c in their showers.  The charm-tag anatomy finds the
# shower piece of Sherpa's charm at 1.06% of the CC rate against 0.14%
# (POWHEG-V2) and 0.03% (Herwig), and 0.44% against 0.03% on the NC, and it
# is the piece that grows with beam energy.  This run tests that reading:
# if the shower charm collapses to the others' level with the mass back in,
# the mechanism is established.
#
# THE KNOWN COST.  With charm massive in the shower the ME -> shower
# interface must put an incoming cbar on shell, which at these energies it
# often cannot (memory: sherpa-dis-traps, trap 1): those events are DISCARDED
# and regenerated.  So the delivered sample is depleted in charm-initiated
# events and its event-level cross-section falls below the integrator's.
# That is measured and reported by analysis/charm_shower_mass_check.py; the
# quantity of interest, shower charm in events whose hard process has none,
# is not affected by it.  This is a diagnostic, not a benchmark entry.
set -e
CUR=${1:-nu}
ORDER=${2:-lo}
NEV=${3:-50000}
NJOBS=${4:-4}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
case "$CUR:$ORDER" in
  mu:lo) BASE=MuonDIS_LO ;;  nu:lo) BASE=NuDIS_LO ;;
  mu:nlo) BASE=MuonDIS_NLO ;; nu:nlo) BASE=NuDIS_NLO_ckm3 ;;
  *) echo "usage: run_sherpa_massive_ps.sh <mu|nu> [lo|nlo] [nev] [njobs]" >&2; exit 1 ;; esac
VARIANT=${BASE}_massivePS
SRC=$SHERPA_RUNS/$BASE
DST=$SHERPA_RUNS/$VARIANT
REPOCARD=$BENCH_REPO/sherpa/Runs/$BASE/Sherpa.yaml
[ -d "$SRC/Process" ] || { echo "no compiled Process/ in $SRC" >&2; exit 1; }
mkdir -p "$DST"
ln -sfn "$REPOCARD" "$DST/Sherpa.yaml"
[ -e "$DST/Process" ] || ln -s "$SRC/Process" "$DST/Process"
SHERPA_BIN=$SHERPA_INSTALL/bin/Sherpa
OVERRIDE='MASSIVE_PS: [4, 5]'
if [ ! -f "$DST/Results.zip" ]; then
    echo "$VARIANT: integrating"
    ( cd "$DST" && bench_run "$SHERPA_BIN" -e 0 -R 4321 "$OVERRIDE" \
        > integ.log 2>&1 ) || { echo "integration failed; see $DST/integ.log" >&2; exit 1; }
fi
for i in $(seq 1 "$NJOBS"); do
  d=$DST/job_$i
  if [ -s "$d/evtfull" ] && [ "${FORCE:-0}" != "1" ]; then
      echo "refusing to overwrite $d/evtfull -- FORCE=1" >&2; exit 1; fi
  mkdir -p "$d"; ln -sfn "$REPOCARD" "$d/Sherpa.yaml"
  [ -e "$d/Process" ] || ln -s "$DST/Process" "$d/Process"
  cp -f "$DST/Results.zip" "$d/"
  ( cd "$d" && bench_run "$SHERPA_BIN" -e "$NEV" "RANDOM_SEED: $((1234 + i))" \
      "$OVERRIDE" 'EVENT_OUTPUT: ["HepMC3_GenEvent[evtfull]"]' > sherpa.log 2>&1 ) &
done
wait
for i in $(seq 1 "$NJOBS"); do
  echo "  job_$i: $(grep -c '^E ' "$DST/job_$i/evtfull") events; $(grep -m1 'Massive PS flavours' "$DST/job_$i/sherpa.log")"
done
echo "$VARIANT done"
