#!/usr/bin/env bash
# The analytic sigma_fid(E) ladder for ONE fiducial region, on one or both
# currents.
#
# WHY A LADDER PER REGION.  faser_rates.py rides every generator on an
# analytic sigma(E) SHAPE, because the FASER flux is broad and handing a whole
# bin its representative energy overestimated the rate by 68% once already.
# The shape is a property of the REGION and cannot be borrowed: the
# 0.2 < y < 0.9 window CLOSES as the beam energy falls (sigma_fid at 10 GeV is
# a tenth of its 400 GeV value, where a power law says a half), while a region
# with no y window does not close the same way.  Reusing the benchmark's shape
# for a foreign region would bury that difference inside R_g, where it is then
# extrapolated below the lowest computed point.
#
# THE BENCHMARK ENERGIES ARE INCLUDED HERE, unlike tools/yadism_energy_ladder.sh
# which skips them because the benchmark proper computes them: for a foreign
# region nothing else does.
#
# Usage:  tools/yadism_ladder_region.sh <selection> [nu|mu|both] [energies...]
#
# Regions this is used for, and what to read before quoting either:
#
#   dis2506  Q > 1.65 GeV, W > 2 GeV -- arXiv:2506.13889's region, muon only.
#   ally     Q2 > 4 GeV2, no y window -- for the arXiv:2402.13318 comparison,
#            both currents.
#
# BOTH REACH x -> 0.999, WHERE THIS CALCULATION HAS NO TARGET-MASS
# CORRECTIONS.  Measured at 1 TeV: the three independent integrators inside
# the yadism drivers agree to six figures on the benchmark region, spread by
# 0.72% on dis2506 and by 0.16% on ally, and yadism prints "Some NaNs are
# encountered and set to zero" on both and not on the benchmark region.  The
# absolute normalisation largely CANCELS in faser_rates, because R_g is formed
# against this same curve at the benchmark energies and only the shape between
# them survives -- but a number taken from one of these ladders directly is
# worth about a per cent less than one from the benchmark ladder, and the
# missing TMCs are a further, unmeasured error on top.
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
. "$HERE/config.sh"
cd "$HERE"

SEL=${1:?usage: yadism_ladder_region.sh <selection> [nu|mu|both] [energies...]}
WHICH=${2:-both}
shift 2 2>/dev/null || shift $#

LADDER="10 15 25 40 60 100 200 300 400 632 1000 1500 2000 3000 4000 6000"
[ $# -gt 0 ] && LADDER="$*"

case "$WHICH" in
  nu)   SCRIPTS="analysis/yadism_cc_calc.py" ;;
  mu)   SCRIPTS="analysis/yadism_calc.py" ;;
  both) SCRIPTS="analysis/yadism_cc_calc.py analysis/yadism_calc.py" ;;
  *) echo "unknown current: $WHICH" >&2; exit 2 ;;
esac

for e in $LADDER; do
  for s in $SCRIPTS; do
    echo "=== $s nlo, $SEL, $e GeV"
    BENCH_SELECTION="$SEL" BENCH_ENERGY="$e" \
        python3 "$s" nlo 2>&1 \
        | grep -E "sigma_fid|Traceback|Error" || true
  done
done
echo "LADDERDONE $SEL $WHICH"
