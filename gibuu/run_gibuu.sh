#!/usr/bin/env bash
# Generate one GiBUU sample and convert it to HepMC3.
#
# Usage: gibuu/run_gibuu.sh [--current nu|nubar] [--target p|n|W] \
#                           <beam energy in GeV> [ensembles] [runs] [job]
#
# >>> READ gibuu/README.md FIRST.  Its neutrino module is documented for
# 1-50 GeV and this benchmark runs at 400 GeV - 4 TeV. <<<  It does run there
# and the energy dependence is right, but every number it produces above
# 50 GeV is an extrapolation and the figures say so.
#
# DIS ONLY, as for NuWro: every other channel is switched off, because this is
# a DIS benchmark and the non-DIS channels are per mille of the fiducial rate.
#
# THE TWO TARGETS DIFFER IN MORE THAN Z AND A.
#   p, n   a free nucleon, numTimeSteps = 0 -- no transport at all, so the
#          final state is the primary vertex and this is like for like with
#          every other generator in this benchmark.
#   W      tungsten (Z = 74, A = 184) with the transport switched on, which
#          is what GiBUU is FOR: the hadrons are propagated out through the
#          nucleus by the Boltzmann-Uehling-Uhlenbeck equation rather than
#          having a cascade bolted on.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"

CURRENT=nu
TARGET=p
while [ $# -gt 0 ]; do
  case "$1" in
    --current) CURRENT=$2; shift 2 ;;
    --target)  TARGET=$2;  shift 2 ;;
    *) break ;;
  esac
done
EBEAM=${1:?usage: run_gibuu.sh [--current nu|nubar] [--target p|n|W] <E_GeV> [ensembles] [runs] [job]}
ENS=${2:-4000}
RUNS=${3:-25}
JOB=${4:-1}

case "$CURRENT" in
  nu)    PROC=2  ;;   # 2 = CC, -2 = antiCC
  nubar) PROC=-2 ;;
  *) echo "--current takes nu or nubar" >&2; exit 1 ;;
esac
# >>> length_perturbative HAS TO GROW WITH THE TARGET, AND IT IS NOT OBVIOUS.
# <<<  It is the per-ensemble room for perturbative particles.  On a free
# nucleon 200 is ample; on tungsten with the transport on, the first attempt
# (200) died with
#     setIntoVector, step 3 reached: no Hole found!
#     error setIntoVector in initNeutrino
# after 400 ensembles -- a 1 TeV DIS event makes tens of hadrons and the
# cascade multiplies them, and every one needs a slot for the whole run.
case "$TARGET" in
  p) Z=1;  A=1;   STEPS=0   ; NUC=p     ; LENP=200  ; DT=0.2  ;;
  n) Z=0;  A=1;   STEPS=0   ; NUC=n     ; LENP=200  ; DT=0.2  ;;
  # >>> AND A SMALLER TIME STEP, because the transport's own collision
  # criterion breaks at TeV energies. <<<  With GiBUU's default delta_T = 0.2
  # fm the run died with
  #     Error in localCollisionCriteria: probability>1
  # -- the per-step collision probability exceeding unity, which means the
  # step is too long for the densities and energies the cascade is being
  # handed.  The distance travelled, numTimeSteps x delta_T, is kept at 24 fm
  # so the hadrons still leave the nucleus.
  W) Z=74; A=184; STEPS=480 ; NUC=mixed ; LENP=6000 ; DT=0.05 ;;
  *) echo "--target takes p, n or W" >&2; exit 1 ;;
esac
# THE ROOM GROWS WITH THE ENERGY TOO.  Every target nucleon of every ensemble
# carries one perturbative event, so the per-ensemble count is set by the
# multiplicity, not by numEnsembles: 6000 suffices at 50 GeV, but at 100 GeV
# (2026-10-01) the run died with "Perturbative particle vector too small!"
# after every ensemble had filled.  BENCH_GIBUU_LENP overrides it.
LENP=${BENCH_GIBUU_LENP:-$LENP}

[ -x "$GIBUU_DIR/testRun/GiBUU.x" ] || {
    echo "no $GIBUU_DIR/testRun/GiBUU.x -- build it, see gibuu/README.md" >&2
    exit 1; }
[ -d "$BUUINPUT" ] || { echo "no \$BUUINPUT at $BUUINPUT" >&2; exit 1; }

BASE_NAME=$(python3 "$BENCH_REPO/analysis/beams.py" name \
                "gibuu_job_${CURRENT}_${TARGET}" "$EBEAM")
d=$HERE/${BASE_NAME}_${JOB}
mkdir -p "$d"
echo "GiBUU: $CURRENT on $TARGET at $EBEAM GeV, $ENS x $RUNS -> $(basename "$d")"

# THE JOBCARD IS WRITTEN HERE AND NOT KEPT AS A TEMPLATE WITH HOLES.  Every
# value in it depends on the run, and a template with six substitutions is a
# template nobody can read.  The variables that are NOT obvious:
#   numTimeSteps = 0   no transport; the free-nucleon runs want the primary
#                      vertex and nothing else
#   nuXsectionMode = 6 dSigmaMC, i.e. a monoenergetic beam whose energy is
#                      &nl_SigmaMC enu (in GeV)
#   nuExp = 0          no experimental flux folded in
#   length_perturbative  the per-ensemble room for perturbative particles; at
#                      TeV energies the final states are long and the default
#                      is not enough
cat > "$d/job.card" <<EOF
! Written by gibuu/run_gibuu.sh -- do not edit, edit the driver.
&input
      eventtype           = 5
      numEnsembles        = $ENS
      numTimeSteps        = $STEPS
      num_runs_SameEnergy = $RUNS
      num_Energies        = 1
      delta_T             = $DT
      path_To_Input       = '$BUUINPUT'
      length_perturbative = $LENP
      localEnsemble       = .true.
/
&initRandom
      SEED = $((1000 + JOB))
/
&target
      target_Z = $Z
      target_A = $A
/
&neutrino_induced
      process_ID     = $PROC
      flavor_ID      = 2
      nuXsectionMode = 6
      nuExp          = 0
      includeQE      = .false.
      includeDELTA   = .false.
      includeRES     = .false.
      include1pi     = .false.
      include2p2hQE  = .false.
      include2pi     = .false.
      includeDIS     = .true.
/
&nl_SigmaMC
      enu = $EBEAM
/
&EventOutput
      WritePerturbativeParticles = .true.
      EventFormat = 1
/
EOF

( cd "$d" && bench_run "$GIBUU_DIR/testRun/GiBUU.x" < job.card > run.log 2>&1 )

# >>> ONE LES HOUCHES FILE PER RUN, NOT ONE PER JOB. <<<  GiBUU writes
# EventOutput.Pert.<run>.lhe for each of num_runs_SameEnergy, so taking the
# first would have silently thrown away (RUNS - 1)/RUNS of the sample -- and
# the cross-section would still have looked right, because it is a sum of
# weights over whatever was read.  They are concatenated instead, and the
# count is checked against the number of runs.
lhes=$(ls "$d"/EventOutput.Pert.*.lhe 2>/dev/null)
nlhe=$(echo "$lhes" | grep -c . || true)
[ "$nlhe" -gt 0 ] || { echo "  !! GiBUU wrote no event file; see $d/run.log" >&2; exit 1; }
[ "$nlhe" -eq "$RUNS" ] || {
    echo "  !! $nlhe Les Houches files for $RUNS runs -- GiBUU stopped early;" >&2
    echo "     see $d/run.log" >&2; exit 1; }
python3 "$HERE/gibuu2hepmc.py" --current "$CURRENT" --nucleon "$NUC" \
        --out "$d/events.hepmc" $lhes
# The Les Houches file is an intermediate (CONVENTIONS.md rule 1) and is several
# times the size of what it converts to.
rm -f "$d"/EventOutput.Pert.*.lhe
echo "  $(grep -c '^E ' "$d/events.hepmc") events in $(basename "$d")/events.hepmc"
