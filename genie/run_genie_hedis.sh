#!/usr/bin/env bash
# GENIE nu_mu DIS production on a free proton via HEDIS
# (GHE19_00a = BGR18 NLO structure functions through APFEL).
#
# Usage: ./run_genie_hedis.sh <CC|NC> <events_per_job> <njobs|job list> [beam_gev]
#        ./run_genie_hedis.sh spline [beam_gev] [CC|NC]   # build the spline only
#
# The muon twin is run_genie.sh; kept in step deliberately (CONVENTIONS.md rule 2b).
#
# ENERGY.  The ANCHOR (1000) keeps its original nucc_job_N directories and the
# _e2000 spline, so re-running it cannot move a published number; other
# energies get nucc_job_<tag>_N from beams.at_energy().  The spline range comes
# from genie_spline_emax() in genie_setup.sh -- too short a spline does not
# fail, it hangs gevgen forever in "Could not select interaction".
#
# FIXED 2026-08-25: the job directory was built with "nu${MODE:l}_job_$i", a
# ZSH lowercase expansion.  Under bash -- which the shebang asks for, and which
# CONVENTIONS.md fixes as the target -- ${MODE:l} is SUBSTRING expansion with the
# unset variable l as the offset, i.e. ${MODE:0}, so it silently yields
# "nuCC_job_1" while analyze_nu.py globs "nucc_job_*".  It never errored; the
# committed samples simply date from a run under zsh.  Same family as the
# `for i in $ARR` trap in CONVENTIONS.md.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
source $HERE/genie_setup.sh > /dev/null

TUNE=GHE19_00a_00_000

if [ "$1" = "spline" ]; then EBEAM=${2:-1000}; MODE=${3:-CC}
else MODE=${1:-CC}; NEV=${2:-25000}; JOBS=${3:-4}; EBEAM=${4:-1000}; fi

# explicit lowercase tag -- see the FIXED note above; bash 3.2 has no ${v,,}
case $MODE in
  CC) LIST=CCHEDIS; SUB=cc;;
  NC) LIST=NCHEDIS; SUB=nc;;
  *) echo "mode must be CC or NC"; exit 1;;
esac
EMAX=$(genie_spline_emax "$EBEAM")
SPLINE=$HERE/splines/numu_${SUB}hedis_p_e$EMAX.xml

if [ "$1" = "spline" ]; then
  mkdir -p $HERE/splines
  if [ -s "$SPLINE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "$SPLINE already exists; FORCE=1 to rebuild."
    echo "  NB rebuilding the _e2000 spline would change the ANCHOR's"
    echo "  sigma_gen through different knot spacing.  Do not."
    exit 0
  fi
  echo "building $SPLINE ($LIST, nu_mu on H1, up to $EMAX GeV) -- this is slow"
  (cd $HERE/splines && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      gmkspl -p 14 -t 1000010010 -e $EMAX --tune $TUNE \
             --event-generator-list $LIST -o $SPLINE" \
      > $HERE/splines/gmkspl_${SUB}hedis_e$EMAX.log 2>&1)
  echo "spline done: $SPLINE"
  exit 0
fi

# A missing spline does not fail fast: gevgen evaluates it as 0 and spins
# forever in "Could not select interaction".  NOTE: the NC spline has never
# been produced (gmkspl_nchedis.log is a stub), so NC mode stops here.
[ -f "$SPLINE" ] || {
  echo "missing spline: $SPLINE"
  echo "  build it with: $0 spline $EBEAM $MODE"
  exit 1; }

JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name "nu${SUB}_job" "$EBEAM")

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

echo "GENIE HEDIS $MODE: ${#JOBLIST[@]} x $NEV events, nu_mu $EBEAM GeV -> ${JOBBASE}_N (spline e$EMAX)"
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  mkdir -p $d
  # re-source INSIDE the wrapper shell (macOS SIP strips DYLD_LIBRARY_PATH); a no-op on Linux
  (cd $d && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      gevgen -n $NEV -p 14 -t 1000010010 -e $EBEAM \
        --tune $TUNE --event-generator-list $LIST \
        --cross-sections $SPLINE --seed $((6543 + i)) \
        -o events.ghep.root > genie.log 2>&1 \
      && $HERE/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && python3 $HERE/spline_to_json.py $SPLINE $EBEAM $NEV > events_xsec.json") &
done
wait
echo "all HEDIS $MODE jobs done"
grep -h "sigma_gen_mb" $HERE/${JOBBASE}_*/events_xsec.json
