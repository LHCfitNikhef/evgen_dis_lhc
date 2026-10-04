#!/usr/bin/env bash
# Re-shower a POWHEG-V2 LHE so its HepMC carries `lhe_index`.
#
# Usage: ./reshower_v2_indexed.sh [--current mu|nu] <beam_energy_gev> [lhe_path]
#
# WHY.  A PDF band on a HADRON-LEVEL selection -- the dimuon tier above all --
# needs the member weights (which live in the LHE, harvested to weights.npz)
# joined to the showered event (which is what the selection cuts on).  The
# join is by `lhe_index`, an attribute `main_powheg` stamps on every HepMC
# event, and it CANNOT be replaced by position: Pythia aborts some events
# (940 of 50000 at 1 TeV), so the showered file is shorter than the LHE and
# every weight after the first abort would belong to the wrong event.
#
# THE PUBLISHED PRODUCTION SAMPLES PREDATE THAT ATTRIBUTE.  nu_job_1,
# nu_job_400GeV_1 and nu_job_4TeV_1 carry no lhe_index, so they cannot be
# joined to the weights at all -- only the 1 TeV rwgtnu_job_1 sample can.
# This re-showers the SAME LHE with the current main_powheg to produce an
# indexed twin.  The physics is unchanged; what is added is the bookkeeping.
#
# THE MUON SIDE NEEDS THE SAME THING, for the SCALE band rather than the PDF
# one (user, 2026-09-08: "can you also please show MHOUs in the POWHEG-
# prediction?").  powheg/v2mu_job_1 was showered on 2026-08-27, before the
# attribute existed, so it carries no index either; and the muon band has to
# come from POWHEG-V2 in any case, because POWHEG-RES runs with `manyseeds`
# over a hundred-odd Les Houches files and is not reweightable in one pass.
#
# Output goes to powheg/rwgt<base>_job_<TAG>_1 so it never collides with the
# published sample, which stays the source of every published number.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$BENCH_REPO/powheg/powheg_variants.sh"

CURRENT=nu
if [ "$1" = "--current" ]; then CURRENT=$2; shift 2; fi
case "$CURRENT" in mu|nu) ;; *) echo "--current takes mu or nu" >&2; exit 1;; esac

EBEAM=${1:?usage: reshower_v2_indexed.sh [--current mu|nu] <beam_energy_gev> [lhe_path]}
powheg_variant V2 "$CURRENT" "$EBEAM"
TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$CURRENT" "$EBEAM")
LHE=${2:-$PV_RUNDIR/pwgevents.lhe}
[ -s "$LHE" ] || { echo "no LHE at $LHE" >&2; exit 1; }

case "$CURRENT" in mu) BASENAME=rwgtv2mu_job ;; nu) BASENAME=rwgtnu_job ;; esac
d=$HERE/${BASENAME}_${TAG}_1
mkdir -p "$d"
echo "re-showering $LHE ($(du -h "$LHE" | cut -f1)) -> $(basename "$d")"
( cd "$HERE" && bench_run ./main_powheg "$PV_CMND" "$LHE" "$d/events" \
      > "$d/shower.log" 2>&1 )
n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null); n=${n:-0}
i=$(grep -c 'lhe_index' "$d/events.hepmc" 2>/dev/null); i=${i:-0}
echo "  $n events, $i carrying lhe_index"
[ "$n" -gt 0 ] && [ "$i" -eq "$n" ] || {
    echo "  !! every event must carry lhe_index -- the join depends on it" >&2
    exit 1; }
