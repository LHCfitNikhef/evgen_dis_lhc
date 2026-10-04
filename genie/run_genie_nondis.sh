#!/usr/bin/env bash
# GENIE default tune, FULLY INCLUSIVE generation for the non-DIS study
# (layer 2: the differential distributions, user request 2026-09-04).  No
# selection, no acceptance cut, no Q2 or y window at generation: every
# channel of the default event-generator list, on a free proton.
#
# Two runs per (current, energy):
#   full   the stock list (CC for neutrinos, EM for muons): the inclusive
#          sample, every channel at its own rate.  This is the closure --
#          its per-channel event fractions must reproduce the spline ones.
#   nondis the same list MINUS the DIS generators (CCNONDIS / EMNONDIS from
#          make_nondis_overlay.py): the non-DIS shapes with full statistics,
#          to be normalised to the spline cross-sections in the analysis.
#
# Output per job: events.gst.root, GENIE's summary ntuple (gntpc -f gst),
# which carries the per-event channel flags (qel, res, dis, dfr, charm) and
# the lepton kinematics (Q2, x, y, W, El, cthl).  The GHEP file is deleted
# once the ntuple is written: it is 4x larger and nothing reads it
# (CONVENTIONS.md rule 1, intermediates do not outlive their results).
#
# Usage: ./run_genie_nondis.sh <mu|nu> <beam_gev> [nev_full] [nev_nondis] [njobs] [p|n]
#        defaults 1000000 200000 8 (per-job counts are the totals / njobs)
#
# THE TARGET ARGUMENT IS THE  PRODUCTION (user, 2026-09-14: "regenerate
# GENIE's non-DIS sample on p+n, so corresponding to a tungsten nucleus").
# Given, it runs on that free nucleon and names the jobs
# v2nd{full,only}_<cur>_<p|n>_job[_TAG]_N, combined per nucleon in the
# analysis as (74 p + 110 n) / 184.  Omitted, it is the earlier proton run and its
# earlier names (ndfull_/ndonly_), whose event files were deleted 2026-09-14.
# The neutron seeds are offset by 500 so the two nucleons draw disjoint
# random streams (their fluctuations would otherwise correlate in the sum).
#
# THE MUON SIDE IS Q2 > 4 ONLY, by this build's EM floor (genie/README.md,
# patches/genie-em-q2-floor-and-pythia8-teardown.diff): "fully inclusive"
# for muons means every EM channel above 4 GeV2.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
source "$BENCH_REPO/config.sh" > /dev/null
source $HERE/genie_setup.sh > /dev/null

CUR=${1:?usage: run_genie_nondis.sh <mu|nu> <beam_gev> [nev_full] [nev_nondis] [njobs]}
EBEAM=${2:?beam energy in GeV}
NFULL=${3:-1000000}
NNONDIS=${4:-200000}
JOBS=${5:-8}
TGT=${6:-}
TUNE=G18_02a_00_000
case "$TGT" in
  "") NUC=p; TCODE=1000010010; SOFF=0;   STEMF=ndfull_${CUR}_job;         STEMO=ndonly_${CUR}_job;;
  p)  NUC=p; TCODE=1000010010; SOFF=0;   STEMF=v2ndfull_${CUR}_p_job;     STEMO=v2ndonly_${CUR}_p_job;;
  n)  NUC=n; TCODE=1000000010; SOFF=500; STEMF=v2ndfull_${CUR}_n_job;     STEMO=v2ndonly_${CUR}_n_job;;
  *) echo "target must be p or n" >&2; exit 2;;
esac
case "$CUR" in
  nu) PID=14; FULL=CC; NONDIS=CCNONDIS; SPLINE=$HERE/splines/faser/14_${NUC}_cc_e8000.xml;;
  mu) PID=13; FULL=EM; NONDIS=EMNONDIS; SPLINE=$HERE/splines/faser/13_${NUC}_em_e8000.xml;;
  *) echo "current must be mu or nu" >&2; exit 2;;
esac
# the full-list splines reach 8 TeV (tools/genie_faser_splines.sh,
# tools/genie_nondis_splines.sh); a missing spline HANGS gevgen rather than
# failing it, so check first
[ -s "$SPLINE" ] || { echo "missing spline $SPLINE" >&2; exit 1; }
python3 $HERE/make_nondis_overlay.py > /dev/null
# the overlay goes FIRST (its lists and its GVLD-Emax), the p8 overlay that
# genie_setup.sh prepared still supplies everything else
export GXMLPATH=$HERE/config/nondis:$GXMLPATH

run_set() {   # <list> <nev_total> <jobbase-stem>
  local LIST=$1 NEV=$2 STEM=$3
  local PER=$((NEV / JOBS))
  local JOBBASE
  JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name "$STEM" "$EBEAM")
  for i in $(seq 1 $JOBS); do
    d=$HERE/${JOBBASE}_$i
    if [ -s "$d/events.gst.root" ] && [ "${FORCE:-0}" != "1" ]; then
      echo "exists: $d/events.gst.root (FORCE=1 to redo)"; continue
    fi
    mkdir -p $d
    # the GENIE environment is re-sourced INSIDE the wrapper: macOS strips
    # DYLD_LIBRARY_PATH across caffeinate (as every run_genie*.sh does)
    (cd $d && bench_run "$BENCH_SHELL" -c "
        source $HERE/genie_setup.sh > /dev/null
        export GXMLPATH=$HERE/config/nondis:\$GXMLPATH
        gevgen -n $PER -p $PID -t $TCODE -e $EBEAM \
          --tune $TUNE --event-generator-list $LIST \
          --cross-sections $SPLINE --seed $((2468 + SOFF + 100 * ${#LIST} + i)) \
          -o events.ghep.root > genie.log 2>&1 \
        && gntpc -i events.ghep.root -f gst -o events.gst.root >> genie.log 2>&1 \
        && rm -f events.ghep.root \
        && echo 'list=$LIST nev=$PER beam=$EBEAM pid=$PID target=$TCODE spline=$SPLINE' > events.meta") &
  done
  wait
  # gevgen prints every event record: 1 GB of log per job.  Reduce each to
  # its rejection counts (the mechanism behind the delivered channel mix,
  # see nondis_log_digest.py) plus a 200-line tail, then delete it.
  python3 $HERE/nondis_log_digest.py $HERE/${JOBBASE}_[0-9]*
  echo "$LIST: $JOBS x $PER events -> ${JOBBASE}_N"
}

echo "GENIE non-DIS study: $CUR on $NUC ($TCODE) at $EBEAM GeV, full list $FULL ($NFULL) and $NONDIS ($NNONDIS)"
run_set $FULL   $NFULL   $STEMF
run_set $NONDIS $NNONDIS $STEMO
echo "done: $CUR $EBEAM GeV"
