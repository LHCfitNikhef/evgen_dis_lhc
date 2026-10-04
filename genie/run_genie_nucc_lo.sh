#!/usr/bin/env bash
# GENIE nu_mu CC DIS production at 1 TeV on a free proton, CLASSIC LO path:
# tune G18_02a_00_000 = Bodek-Yang DIS built on GRV98LO, AGKY/Pythia8
# hadronisation through the benchmark config overlay.  This is the LO
# counterpart of run_genie_hedis.sh, whose GHE19_00a tune is BGR18 NLO.
#
# The event-generator list is CCDISCHARM, written by make_nucc_lo_overlay.py into a
# throwaway overlay directory: GENIE's inclusive DIS-CC cross-section has the
# charm piece SUBTRACTED internally (QPMDISPXSec.cxx), so DIS-CC-CHARM has to
# run alongside it -- see that script's docstring.
#
# Usage: ./run_genie_nucc_lo.sh [events_per_job] [njobs|job list "1 4"] [beam_gev]
#        ./run_genie_nucc_lo.sh spline [beam_gev]     # build the spline only
#
# ENERGY.  The ANCHOR (1000) keeps its original nucclo_job_N directories and
# the _e2000 spline, so re-running it cannot move a published number; other
# energies get nucclo_job_<tag>_N from beams.at_energy() and the spline range
# genie_spline_emax() prescribes.  Kept in step with run_genie.sh and
# run_genie_hedis.sh (CONVENTIONS.md rule 2b).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
LIST=CCDISCHARM
TUNE=G18_02a_00_000

source $HERE/genie_setup.sh > /dev/null

# the beam energy sits in a different argument for the two entry points
if [ "$1" = "spline" ]; then EBEAM=${2:-1000}; else EBEAM=${3:-1000}; fi
EMAX=$(genie_spline_emax "$EBEAM")
SPLINE=$HERE/splines/numu_ccdis_p_e$EMAX.xml
if [ "$EMAX" -gt 2000 ]; then
  export GXMLPATH=$HERE/config/hienergy:$GXMLPATH
  echo "note: GVLD-Emax raised via config/hienergy -- $TUNE is declared valid"
  echo "      to 1000 GeV, so $EBEAM GeV is an EXTRAPOLATION of the tune."
fi
python3 $HERE/make_nucc_lo_overlay.py > /dev/null
# the extra list overlay goes FIRST; everything else still comes from the p8
# overlay that genie_setup.sh set up
export GXMLPATH=$HERE/config/nu_cc_lo:$GXMLPATH
# ABOVE 1 TeV, raise the declared validity range so a long enough spline can be
# built at all.  PREPENDED BEFORE nu_cc_lo: first match wins on GXMLPATH, and
# nu_cc_lo's own CommonParam.xml carries GVLD-Emax = 5000, short of the 8000
# the 4 TeV spline needs.  nu_cc_lo's EventGeneratorListAssembler.xml is still
# found there, since this overlay holds no such file.  Same tune and the same
# caveat as the muon side: G18_02a is declared valid to 1000 GeV, so 4 TeV
# EXTRAPOLATES past it.  See genie/config/hienergy/README.md.

build_spline() {
  mkdir -p $HERE/splines
  if [ -s "$SPLINE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "$SPLINE already exists; FORCE=1 to rebuild."
    echo "  NB rebuilding the _e2000 spline would change the ANCHOR's"
    echo "  sigma_gen through different knot spacing.  Do not."
    return 0
  fi
  echo "building $SPLINE ($LIST, nu_mu on H1, up to $EMAX GeV) -- this is slow"
  # the spline must extend WELL past the beam energy: GENIE evaluates the last
  # knot by reading past the knot array and returns 0 (see run_genie.sh)
  (cd $HERE/splines && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$HERE/config/nu_cc_lo:\$GXMLPATH
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gmkspl -p 14 -t 1000010010 -e $EMAX --tune $TUNE \
             --event-generator-list $LIST -o $SPLINE" \
      > $HERE/splines/gmkspl_ccdis_e$EMAX.log 2>&1)
  echo "spline done: $SPLINE"
}

if [ "$1" = "spline" ]; then
  build_spline
  exit 0
fi

NEV=${1:-25000}
JOBS=${2:-4}
[ -f "$SPLINE" ] || build_spline
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name nucclo_job "$EBEAM")

# NOTE "${JOBLIST[@]}" below, not $JOBLIST: bash expands the latter to the
# FIRST element only, which would silently run one job instead of N.
case "$JOBS" in (*' '*) JOBLIST=($JOBS);; (*) JOBLIST=($(seq 1 $JOBS));; esac

for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "GENIE CC LO: ${#JOBLIST[@]} x $NEV events, nu_mu $EBEAM GeV -> ${JOBBASE}_N (spline e$EMAX)"
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  mkdir -p $d
  # re-source INSIDE the wrapper shell (macOS (SIP strips DYLD_LIBRARY_PATH); a no-op on Linux
  (cd $d && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$HERE/config/nu_cc_lo:\$GXMLPATH
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gevgen -n $NEV -p 14 -t 1000010010 -e $EBEAM \
        --tune $TUNE --event-generator-list $LIST \
        --cross-sections $SPLINE --seed $((4321 + i)) \
        -o events.ghep.root > genie.log 2>&1 \
      && $HERE/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && python3 $HERE/spline_to_json.py $SPLINE $EBEAM $NEV > events_xsec.json") &
done
wait
echo "all CC LO jobs done"
grep -h "sigma_gen_mb" $HERE/${JOBBASE}_*/events_xsec.json
