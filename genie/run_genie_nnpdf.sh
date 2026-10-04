#!/usr/bin/env bash
# GENIE muon-DIS production with NNPDF40_nnlo_as_01180 replacing GRV98 LO
# inside Bodek-Yang (overlay genie/config/nnpdf). Mirrors run_genie.sh;
# jobs land in nnpdf_job_<i>/ (or nnpdf_job_<tag>_<i> away from the anchor).
#
# Usage: ./run_genie_nnpdf.sh <events_per_job> <njobs|job list "1 4"> [beam_gev]
#        ./run_genie_nnpdf.sh spline [beam_gev]     # build the spline only
#
# ENERGY-AWARE SINCE 2026-08-27.  It used to hardcode the spline name AND
# `-e 1000`, so it could only ever make the anchor -- which is why the muon
# NNPDF row existed at 1 TeV alone while every other GENIE row spanned the
# scan.  The ANCHOR keeps its original nnpdf_job_N directories so nothing
# published moves; other energies get nnpdf_job_<tag>_N from beams.at_energy().
#
# NB genie_spline_emax maps BOTH 400 and 1000 GeV to the e2000 spline, so
# 400 GeV needs no new spline -- only 4 TeV does, and there G18_02a is being
# EXTRAPOLATED past its declared 1000 GeV validity (see genie/config/hienergy).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
source $HERE/genie_setup.sh > /dev/null

if [ "$1" = "spline" ]; then EBEAM=${2:-1000}; else EBEAM=${3:-1000}; fi
EMAX=$(genie_spline_emax "$EBEAM")
NEV=${1:-25000}
JOBS=${2:-8}
SPLINE=$HERE/splines/mu_emdis_p_nnpdf40_e$EMAX.xml

build_spline() {
  mkdir -p $HERE/splines
  if [ -s "$SPLINE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "$SPLINE already exists; FORCE=1 to rebuild."
    echo "  NB rebuilding the e2000 spline would change the ANCHOR's"
    echo "  sigma_gen through different knot spacing.  Do not."
    return 0
  fi
  echo "building $SPLINE (EMDIS, mu- on H1, NNPDF4.0, up to $EMAX GeV)"
  echo "  this is slow"
  (cd $HERE/splines && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$HERE/config/nnpdf:\$GXMLPATH
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gmkspl -p 13 -t 1000010010 -e $EMAX --tune G18_02a_00_000 \
             --event-generator-list EMDIS -o $SPLINE" \
      > $HERE/splines/gmkspl_nnpdf40_e$EMAX.log 2>&1)
  [ -s "$SPLINE" ] || { echo "gmkspl produced nothing -- see" \
      "splines/gmkspl_nnpdf40_e$EMAX.log" >&2; exit 1; }
  echo "spline done: $SPLINE"
}

if [ "$1" = "spline" ]; then build_spline; exit 0; fi

# A missing spline does not fail fast: gevgen evaluates it as 0 and spins
# forever in "Could not select interaction".
[ -f "$SPLINE" ] || build_spline
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name nnpdf_job "$EBEAM")

# NOTE "${JOBLIST[@]}" below, not $JOBLIST: bash expands the latter to the
# FIRST element only, which would silently run one job instead of N.
case "$JOBS" in (*' '*) JOBLIST=($JOBS);; (*) JOBLIST=($(seq 1 $JOBS));; esac
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  mkdir -p $d
  # re-source INSIDE the wrapper shell (macOS (SIP strips DYLD_LIBRARY_PATH); a no-op on Linux,
  # then select the NNPDF overlay
  (cd $d && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$HERE/config/nnpdf
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gevgen -n $NEV -p 13 -t 1000010010 -e $EBEAM \
        --tune G18_02a_00_000 --event-generator-list EMDIS \
        --cross-sections $SPLINE --seed $((7654 + i)) \
        -o events.ghep.root > genie.log 2>&1 \
      && $HERE/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && python3 $HERE/spline_to_json.py $SPLINE $EBEAM $NEV > events_xsec.json") &
done
wait
echo "all NNPDF jobs done"
grep -h "sigma_gen_mb" $HERE/${JOBBASE}_*/events_xsec.json