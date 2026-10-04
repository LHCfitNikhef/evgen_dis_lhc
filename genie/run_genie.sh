#!/usr/bin/env bash
# Parallel GENIE production for the muon-DIS benchmark: N seeded jobs.
# mu- on a free proton, EMDIS list, tune G18_02a (Pythia8 overlay via
# GXMLPATH from genie/genie_setup.sh). Requires the benchmark GENIE build
# (EM Q2 floor raised to 4 GeV^2) and the precomputed EMDIS spline.
#
# Usage: ./run_genie.sh <events_per_job> <njobs|job list "1 4 7"> [beam_gev]
#        ./run_genie.sh spline [beam_gev]        # build the spline only
#
# The neutrino twins are run_genie_hedis.sh (NLO, BGR18) and
# run_genie_nucc_lo.sh (LO); kept in step deliberately (CONVENTIONS.md rule 2b).
#
# ENERGY.  The ANCHOR (1000) keeps its original job_N directories and the
# _e2000 spline, so re-running it cannot move a published number; other
# energies get job_<tag>_N from beams.at_energy() and the spline range
# genie_spline_emax() prescribes.  Before this was parameterised the script
# hardcoded `-e 1000` AND `job_$i`, so running it for the scan would have
# overwritten the 1 TeV anchor samples with 4 TeV events -- the anchor
# overwrite trap for the sixth time in this project.
#
# The spline range must extend WELL past the beam energy: evaluating a GENIE
# spline exactly at its last knot reads past the knot array
# (Spline::FindClosestKnot GetKnot(iknot+1)) and returns 0 -> an infinite
# "Could not select interaction" loop.  The rule lives in genie_setup.sh.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
source $HERE/genie_setup.sh > /dev/null

TUNE=G18_02a_00_000
LIST=EMDIS

# the beam energy sits in a different argument for the two entry points
if [ "$1" = "spline" ]; then EBEAM=${2:-1000}; else EBEAM=${3:-1000}; fi
EMAX=$(genie_spline_emax "$EBEAM")
SPLINE=$HERE/splines/mu_emdis_p_q2gt4_e$EMAX.xml

# ABOVE 1 TeV, raise GENIE's declared validity range so gmkspl will build a
# spline that long at all (see genie/config/hienergy/README.md).  G18_02a is
# Bodek-Yang, declared valid to 1000 GeV, so the 4 TeV point EXTRAPOLATES past
# it -- deliberately, because Bodek-Yang is what FASER uses and the benchmark
# needs the reference even where the model is outside its stated range.
# Applied ONLY above 1 TeV, so 400 GeV and the 1 TeV anchor keep exactly the
# configuration their published numbers came from.
if [ "$EMAX" -gt 2000 ]; then
  export GXMLPATH=$HERE/config/hienergy:$GXMLPATH
  echo "note: GVLD-Emax raised via config/hienergy -- $TUNE is declared valid"
  echo "      to 1000 GeV, so $EBEAM GeV is an EXTRAPOLATION of the tune."
fi

if [ "$1" = "spline" ]; then
  mkdir -p $HERE/splines
  if [ -s "$SPLINE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "$SPLINE already exists; FORCE=1 to rebuild."
    echo "  NB rebuilding the _e2000 spline would change the ANCHOR's"
    echo "  sigma_gen through different knot spacing.  Do not."
    exit 0
  fi
  echo "building $SPLINE (EMDIS, mu- on H1, up to $EMAX GeV) -- this is slow"
  (cd $HERE/splines && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gmkspl -p 13 -t 1000010010 -e $EMAX --tune $TUNE \
             --event-generator-list $LIST -o $SPLINE" \
      > $HERE/splines/gmkspl_emdis_e$EMAX.log 2>&1)
  echo "spline done: $SPLINE"
  exit 0
fi

# A missing spline does not fail fast: gevgen evaluates it as 0 and spins
# forever in "Could not select interaction".
[ -f "$SPLINE" ] || {
  echo "missing spline: $SPLINE"
  echo "  build it with: $0 spline $EBEAM"
  exit 1; }

NEV=${1:-25000}
JOBS=${2:-8}
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name job "$EBEAM")

# NOTE "${JOBLIST[@]}" below, not $JOBLIST: bash expands the latter to the
# FIRST element only, which would silently run one job instead of N.
case "$JOBS" in (*' '*) JOBLIST=($JOBS);; (*) JOBLIST=($(seq 1 $JOBS));; esac

# genie/job_1 and friends hold real production samples; refuse to clobber them.
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "GENIE mu: ${#JOBLIST[@]} x $NEV events, mu- $EBEAM GeV -> ${JOBBASE}_N (spline e$EMAX)"
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  mkdir -p $d
  # macOS SIP strips DYLD_LIBRARY_PATH when exec'ing a protected binary
  # (caffeinate, /bin/bash), so the GENIE environment must be re-sourced
  # INSIDE the wrapper.  Harmless on Linux, where the re-source is a no-op.
  (cd $d && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gevgen -n $NEV -p 13 -t 1000010010 -e $EBEAM \
        --tune $TUNE --event-generator-list $LIST \
        --cross-sections $SPLINE --seed $((8765 + i)) \
        -o events.ghep.root > genie.log 2>&1 \
      && $HERE/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && python3 $HERE/spline_to_json.py $SPLINE $EBEAM $NEV > events_xsec.json") &
done
wait
echo "all ${#JOBLIST[@]} GENIE jobs done"
grep -h "sigma_gen_mb" $HERE/${JOBBASE}_*/events_xsec.json
