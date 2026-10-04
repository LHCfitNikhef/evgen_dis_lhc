#!/usr/bin/env bash
# POWHEG-V2mc FOR PAPER PLOTS V2: massive-charm neutrino CC on ONE nucleon.
#
# Usage: powheg/production/run_v2_mc.sh <energy_gev> <p|n>
#
# User, 2026-09-14: "massive-charm run is still very much part of the
# plan ... (and launch the corresponding generation!)".  The massless 
# driver run_v2.sh, pointed at a SEPARATE code entry through the knobs it
# already has, so nothing of the massless sample can be touched:
#   card     powheg/cards/production/POWHEG-V2mc/nu<TAG>-<t>   (make_cards.py: qmass 1.51,
#                                                      numflav 4, no iupperfsr)
#   run dir  $POWHEG_V2/v2-numc<TAG>-<t>[-b<j>]
#   samples  powheg/v2numc_<t>_job[_TAG]_<j>
#   seeds    the massless stream + 50000, disjoint from it
# Showered with powheg_nu_v2.cmnd (LesHouches:matchInOut = off), like the
# massless sample.  The sample is CHARM PRODUCTION ONLY: every event carries
# the massive charm, so its cross-section sits beside charm, never inclusive.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
E=${1:?usage: run_v2_mc.sh <energy_gev> <p|n>}
T=${2:?usage: run_v2_mc.sh <energy_gev> <p|n>}
TAG=$(v2_tag nu "$E")
export V2_CARD=$BENCH_REPO/powheg/cards/production/POWHEG-V2mc/nu$TAG-$T/powheg.input
export V2_RUNDIR=$POWHEG_V2/v2-numc$TAG-$T
export V2_JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name "v2numc_${T}_job" "$E")
export V2_SEEDBASE=$(( $(v2_seed_base "$T" "$E") + 50000 ))
grep -Eq '^qmass[[:space:]]+1\.51d0' "$V2_CARD" || {
    echo "$V2_CARD is not the massive-charm card -- refusing" >&2; exit 1; }
exec "$HERE/run_v2.sh" "$E" "$T"
