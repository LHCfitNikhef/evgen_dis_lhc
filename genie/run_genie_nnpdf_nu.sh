#!/usr/bin/env bash
# GENIE nu_mu CC DIS with NNPDF40_nnlo_as_01180 replacing GRV98 LO inside
# Bodek-Yang.  The neutrino counterpart of run_genie_nnpdf.sh (CONVENTIONS.md rule
# 2b), and the CC counterpart of run_genie_nucc_lo.sh.
#
# WHY THIS SAMPLE EXISTS (user request, 2026-08-27).  Every other LO generator
# in the benchmark now runs NNPDF4.0 NNLO, and GENIE does not: its G18_02a
# tune builds Bodek-Yang on GRV98 LO.  So the GENIE row differs from the rest
# in TWO ways at once -- the parton densities AND everything else GENIE does
# differently (Bodek-Yang low-Q2 modelling, AGKY hadronisation, its own
# structure-function treatment).  With the PDF held common, the residual
# difference isolates "everything else", which is the quantity of interest.
#
# The muon side already has this (nnpdf_job_*, results/histos_genie_nnpdf*);
# this closes the gap on the neutrino side.
#
# HOW THE TWO OVERLAYS COMPOSE, AND WHY THE ORDER IS WHAT IT IS.
# GXMLPATH is FIRST MATCH WINS, PER FILE.  Two overlays are needed here:
#
#   config/nu_cc_lo   EventGeneratorListAssembler.xml -> the CCDISCHARM list.
#                     GENIE's inclusive DIS-CC cross-section has the charm
#                     piece SUBTRACTED internally (QPMDISPXSec.cxx), so
#                     DIS-CC-CHARM must run alongside it.  Written by
#                     make_nucc_lo_overlay.py.
#   config/nnpdf      BYPDF.xml -> Uncorr-PDF-Set = LHAPDF6/NNPDF40_nnlo_0,
#                     PDF-Q2min = 2.7225; plus LHAPDF6.xml naming the set.
#
# nu_cc_lo goes FIRST so its generator list wins.  nnpdf holds no
# EventGeneratorListAssembler.xml, so nothing shadows it; nu_cc_lo holds no
# BYPDF.xml, so the PDF swap comes through.  The one file BOTH provide is
# CommonParam.xml, and the two copies are BYTE IDENTICAL (both raise
# GVLD-Emax to 5000 from the same source file), so the order cannot change
# which parameters are used -- verified 2026-08-27 before this was written.
#
# TRAP THAT MADE THIS WORTH CHECKING: GXMLPATH FAILS OPEN.  A nonexistent
# directory does not error -- GENIE falls through to the built-in tune and
# silently uses GRV98 after all, which is exactly the sample this script
# exists to avoid producing.  The spline and the events must be built with the
# SAME GXMLPATH or the cross-section and the events disagree without a word.
#
# Usage: ./run_genie_nnpdf_nu.sh [events_per_job] [njobs|job list "1 4"] [beam_gev]
#        ./run_genie_nnpdf_nu.sh spline [beam_gev]     # build the spline only
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
LIST=CCDISCHARM
TUNE=G18_02a_00_000

source $HERE/genie_setup.sh > /dev/null

if [ "$1" = "spline" ]; then EBEAM=${2:-1000}; else EBEAM=${3:-1000}; fi
EMAX=$(genie_spline_emax "$EBEAM")
SPLINE=$HERE/splines/numu_ccdis_p_nnpdf40_e$EMAX.xml

python3 $HERE/make_nucc_lo_overlay.py > /dev/null
# nu_cc_lo first (generator list), nnpdf second (the PDF swap).  See above.
OVERLAY=$HERE/config/nu_cc_lo:$HERE/config/nnpdf
for d in $HERE/config/nu_cc_lo $HERE/config/nnpdf; do
    [ -d "$d" ] || { echo "missing overlay $d -- GXMLPATH fails OPEN, so this" \
                          "would silently run GRV98" >&2; exit 1; }
done

build_spline() {
  mkdir -p $HERE/splines
  if [ -s "$SPLINE" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "$SPLINE already exists; FORCE=1 to rebuild."
    return 0
  fi
  echo "building $SPLINE ($LIST, nu_mu on H1, NNPDF4.0, up to $EMAX GeV)"
  echo "  this is slow -- the muon NNPDF spline took hours"
  (cd $HERE/splines && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$OVERLAY:\$GXMLPATH
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gmkspl -p 14 -t 1000010010 -e $EMAX --tune $TUNE \
             --event-generator-list $LIST -o $SPLINE" \
      > $HERE/splines/gmkspl_ccdis_nnpdf40_e$EMAX.log 2>&1)
  [ -s "$SPLINE" ] || { echo "gmkspl produced nothing -- see" \
      "splines/gmkspl_ccdis_nnpdf40_e$EMAX.log" >&2; exit 1; }
  echo "spline done: $SPLINE"
}

if [ "$1" = "spline" ]; then
  build_spline
  exit 0
fi

NEV=${1:-25000}
JOBS=${2:-8}
# A missing spline does not fail fast: gevgen evaluates it as 0 and spins
# forever in "Could not select interaction".
[ -f "$SPLINE" ] || build_spline
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name nnpdfnu_job "$EBEAM")

# The guard used to sweep jobs 1..20 unconditionally, and it ran BEFORE
# JOBLIST was built.  That made the job-list argument useless for the one thing
# it exists for -- EXTENDING a sample with new job numbers -- because any
# existing job in that range refused the run.  It now checks exactly the jobs
# about to be written, so `... "5 6 7"` adds to a sample rather than being
# blocked by job 1.  See the JOBLIST construction below.

# NOTE "${JOBLIST[@]}" below, not $JOBLIST: bash expands the latter to the
# FIRST element only, which would silently run one job instead of N.
case "$JOBS" in (*' '*) JOBLIST=($JOBS);; (*) JOBLIST=($(seq 1 $JOBS));; esac

for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc -- FORCE=1 to replace" >&2
    exit 1
  fi
done

echo "GENIE CC LO + NNPDF4.0: ${#JOBLIST[@]} x $NEV events, nu_mu $EBEAM GeV" \
     "-> ${JOBBASE}_N (spline e$EMAX)"
for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  mkdir -p $d
  # re-source INSIDE the wrapper shell (macOS SIP strips DYLD_LIBRARY_PATH;
  # a no-op on Linux)
  (cd $d && bench_run "$BENCH_SHELL" -c "
      source $HERE/genie_setup.sh > /dev/null
      export GXMLPATH=$OVERLAY:\$GXMLPATH
      [ "$EMAX" -gt 2000 ] && export GXMLPATH=$HERE/config/hienergy:\$GXMLPATH
      gevgen -n $NEV -p 14 -t 1000010010 -e $EBEAM \
        --tune $TUNE --event-generator-list $LIST \
        --cross-sections $SPLINE --seed $((9876 + i)) \
        -o events.ghep.root > genie.log 2>&1 \
      && $HERE/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && python3 $HERE/spline_to_json.py $SPLINE $EBEAM $NEV > events_xsec.json") &
done
wait
echo "all CC LO NNPDF jobs done"
grep -h "sigma_gen_mb" $HERE/${JOBBASE}_*/events_xsec.json
echo
echo "next: BENCH_ENERGY=$EBEAM analysis/analyze_nu.py genie_nnpdf"
echo "      (the analyze_nu entry still has to be added)"
