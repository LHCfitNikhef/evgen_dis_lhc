#!/usr/bin/env bash
# POWHEG-BOX-V2 DIS at NLO: generate and shower, one energy.
#
# Usage: ./run_powheg_v2.sh [--current mu|nu] [--mc] <beam_energy_gev>
#
#   ./run_powheg_v2.sh 400                     # neutrino CC, as before
#   ./run_powheg_v2.sh 4000
#   ./run_powheg_v2.sh --current mu 1000       # the muon NC cross-variant
#   ./run_powheg_v2.sh --mc 1000               # POWHEG-V2mc: MASSIVE charm, CC
#
# --mc SELECTS A DIFFERENT CODE ENTRY, not a switch on the same one.  It reads
# the POWHEG-V2mc cards, writes to prod-numc<TAG>, showers into numc_job_<TAG>
# and is labelled POWHEG-V2mc throughout, so it can never be confused with the
# massless POWHEG-V2 results -- which stay exactly as they are.  CC only, by
# user decision (2026-08-28); the muon NC card forbids a non-zero qmass.
#
# THE CURRENT IS A PARAMETER, THE CODE IS NOT.  This driver runs POWHEG-V2,
# whose native side is neutrino CC -- hence the default -- but the same single
# invocation generates the muon NC cross-variant, a third independent matching
# on the muon side.  Card, run directory and Pythia steering come from
# powheg_variants.sh, which is also where the NC variant's mandatory q2cut is
# written down.
#
# The other driver is run_powheg_res.sh, and the two are deliberately NOT the
# same shape (CONVENTIONS.md rule 2b asks for parity of coverage, not of mechanism).
# POWHEG-RES runs five parallel stages with manyseeds; POWHEG-V2's cards have
# no parallelstage at all and run ONCE, caching their grids through
# use-old-grid / use-old-ubound.  Hence a single invocation here and a single
# pwgevents.lhe rather than ten.  That is a property of the CODE, so it holds
# for both currents.
#
# The 1 TeV run was done by hand and never committed, exactly as the
# POWHEG-RES one was; this is the missing script.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/powheg_variants.sh"     # powheg_variant / powheg_parse_current

# --mc anywhere in the arguments picks the massive-charm entry
PV_CODE=V2
NEWARGS=()
for a in "$@"; do
    if [ "$a" = "--mc" ]; then PV_CODE=V2MC; else NEWARGS+=("$a"); fi
done
set -- "${NEWARGS[@]}"

powheg_parse_current nu "$@"     # default: POWHEG-V2's native current
set -- "${PV_ARGS[@]}"
if [ "$PV_CODE" = V2MC ] && [ "$PV_CURRENT" != nu ]; then
    echo "POWHEG-V2mc is charged-current only: the muon NC card requires" >&2
    echo "  qmass = 0 and the code enforces it." >&2
    exit 1
fi
EBEAM=${1:?usage: run_powheg_v2.sh [--current mu|nu] <beam_energy_gev>}

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$PV_CURRENT" "$EBEAM")
powheg_variant "$PV_CODE" "$PV_CURRENT" "$EBEAM"
CARD=$PV_CARD
RUNDIR=$PV_RUNDIR
JOBBASE=$PV_JOBBASE
PWHG=$POWHEG_V2/pwhg_main
[ -x "$PWHG" ] || { echo "no POWHEG-V2 binary at $PWHG" >&2; exit 1; }

if [ -s "$RUNDIR/pwgevents.lhe" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $RUNDIR/pwgevents.lhe -- FORCE=1 to replace" >&2
    exit 1
fi
if [ -s "$HERE/${JOBBASE}_1/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $HERE/${JOBBASE}_1/events.hepmc -- FORCE=1" >&2
    exit 1
fi

mkdir -p "$RUNDIR"
powheg_install_card "$CARD" "$RUNDIR/powheg.input"
# POWHEG keeps its own refusal to overwrite events, separate from the guard
# above, and it RETURNS 0 when it fires -- so the stage looks like it ran.
[ "${FORCE:-0}" = "1" ] && rm -f "$RUNDIR"/pwgevents.lhe

echo "$PV_LABEL at NLO: $EBEAM GeV [$TAG] -> $RUNDIR"
( cd "$RUNDIR" && bench_run "$PWHG" > run_v2.log 2>&1 ) || {
    echo "pwhg_main failed; see $RUNDIR/run_v2.log" >&2; exit 1; }
[ -s "$RUNDIR/pwgevents.lhe" ] || {
    echo "NO pwgevents.lhe produced -- check $RUNDIR/run_v2.log" >&2; exit 1; }
echo "  LHE: $(du -h "$RUNDIR/pwgevents.lhe" | cut -f1)"

# The same guard the RES shower step applies -- see run_powheg_shower.sh and
# powheg/strip_lhe_nan.py.  No POWHEG-V2 sample has ever carried a NaN, but a
# silent 14.5% truncation is not something to check for on only one side
# (CONVENTIONS.md rule 2b).
"$BENCH_REPO/powheg/strip_lhe_nan.py" "$RUNDIR/pwgevents.lhe"

echo "showering -> ${JOBBASE}_1"
d=$HERE/${JOBBASE}_1
mkdir -p "$d"
# $PV_CMND is used at EVERY energy on purpose: it sets no beams,
# main_powheg.cc takes them from the LHE, so the beam travels with the events.
# It IS keyed by current, though -- the two cards differ in the lepton.
( cd "$HERE" && bench_run ./main_powheg "$PV_CMND" \
    "$RUNDIR/pwgevents.lhe" "$d/events" > "$d/shower.log" 2>&1 )
echo "  ${JOBBASE}_1: $(grep -c '^E ' "$d/events.hepmc" 2>/dev/null) events, $(cat "$d/events_xsec.json" 2>/dev/null)"
if [ "$PV_CURRENT" = nu ]; then
    echo "next: BENCH_ENERGY=$EBEAM analysis/analyze_nu.py powheg_nu"
else
    echo "next: BENCH_ENERGY=$EBEAM analysis/analyze.py powheg_v2_mu   (once the entry exists)"
fi
