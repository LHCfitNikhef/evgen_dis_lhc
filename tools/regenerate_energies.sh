#!/usr/bin/env bash
# Regenerate the event files away from 1 TeV (REGENERATE_ENERGIES.md).
#
# The 400/700/2000/4000 GeV samples were deleted on 2026-09-14 once the
# integrated results had been analysed from them.  A NEW PER-EVENT PASS away
# from 1 TeV -- the first was the single-inclusive pion study moved to the 
# region (user, 2026-09-18) -- needs them back, and each production driver
# rebuilds only what is missing, so this is the list of commands and nothing
# else.  Every driver is re-invocable; the GENIE one skips on its V2_OK marker
# alone, which the deleted directories still carry, so the marker is cleared
# there first (the seed is a function of cfg, nucleon, energy and job number,
# so the regenerated events are the ones the marker described).
#
# Usage:  tools/regenerate_energies.sh powheg|herwig|genie|sherpa|all
#         (each under nohup with its own log; `all` runs them in sequence)
#
# Do not re-run tools/analyse_production.py on the regenerated samples: the tracked
# histos_*_<tag>.json were analysed from the ORIGINAL samples, and a
# regenerated POWHEG shower or Sherpa run reproduces cross-sections only to
# their statistical error.  A per-event pass (faser_pions.py) is normalised
# by those tracked ladders and reads only shapes from the events.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=$(cd "$HERE/.." && pwd)
. "$BENCH_REPO/config.sh" > /dev/null
LOG=$BENCH_REPO/logs
mkdir -p "$LOG"
STAMP=$(date +%F)

powheg () {   # both codes, both nucleons, all five energies (1 TeV is skipped)
    "$BENCH_REPO/powheg/production/production.sh" all > "$LOG/powheg_v2_regen_$STAMP.log" 2>&1
}

herwig () {
    /bin/bash "$BENCH_REPO/herwig7/production/production.sh" > "$LOG/herwig_v2_regen_$STAMP.log" 2>&1
}

genie () {
    # the pion-study rows: FASER's tune (<= 1 TeV) and HEDIS (every energy)
    # on both nucleons, 300 GeV included (beams.SIDIS_ENERGIES); 8 x 25k
    # per point, of which jobs 1-4 existed before at the energies
    local G=$BENCH_REPO/genie/production/genie_job.sh cfg t e base d
    local specs=${GENIE_SPECS:-"nu_grv 300 400 700;mu_grv 300 400 700;nu_hedis 300 400 700 2000 4000"}
    local IFS_=$IFS
    IFS=';'; set -f; set -- $specs; set +f; IFS=$IFS_
    for spec in "$@"; do
        set -- $spec; cfg=$1; shift
        for e in "$@"; do
            for t in p n; do
                base=$(python3 "$BENCH_REPO/analysis/beams.py" name "v2g_${cfg}_${t}_job" "$e")
                for d in "$BENCH_REPO"/genie/${base}_[0-9]*; do
                    [ -d "$d" ] && [ ! -s "$d/events.hepmc" ] && rm -f "$d/V2_OK"
                done
                # four at a time, twice, rather than eight at once (shared machine)
                "$G" run "$cfg" "$t" "$e" 25000 4
                "$G" run "$cfg" "$t" "$e" 25000 8
            done
        done
    done > "$LOG/genie_v2_regen_$STAMP.log" 2>&1
}

sherpa () {
    ENERGIES="400 700 2000 4000" "$BENCH_REPO/tools/sherpa_production.sh" run \
        > "$LOG/sherpa_v2_regen_$STAMP.log" 2>&1
}

# THE INCLUSIVE-HADRON PASS (user, 2026-09-30): the SIDIS yields figure now
# shows every stable charged hadron, which the pion passes never recorded, so
# the neutrino samples of the whole SIDIS ladder come back -- the NEUTRINO
# current only (the study draws only it, CONVENTIONS.md 2b carve-out), POWHEG-V2's
# points below 300 GeV included.  Then tools/faser_pions_passes.sh with
# MU_KEYS=" " FORCE=1, and tools/prune_energies.py as before.
sidis_nu () {
    local e t
    for e in 300 400 700 2000 4000 20 50 100 200; do
        for t in p n; do
            V2_PAR=2 "$BENCH_REPO/powheg/production/run_v2.sh" "$e" "$t" \
                >> "$LOG/powheg_v2_sidisnu_regen_$STAMP.log" 2>&1 \
                || echo "POWHEG-V2 $e $t FAILED" >> "$LOG/powheg_v2_sidisnu_regen_$STAMP.log"
        done
    done
    V2_CURRENTS=nu /bin/bash "$BENCH_REPO/herwig7/production/production.sh" \
        > "$LOG/herwig_v2_sidisnu_regen_$STAMP.log" 2>&1
    GENIE_SPECS="nu_grv 300 400 700;nu_hedis 300 400 700 2000 4000" genie
    CURRENTS=nu ENERGIES="300 400 700 2000 4000" "$BENCH_REPO/tools/sherpa_production.sh" run \
        > "$LOG/sherpa_v2_sidisnu_regen_$STAMP.log" 2>&1
}

case "${1:-}" in
    powheg|herwig|genie|sherpa) "$1" ;;
    sidis-nu) sidis_nu ;;
    all) powheg; herwig; genie; sherpa ;;
    *) echo "usage: $0 powheg|herwig|genie|sherpa|all|sidis-nu" >&2; exit 1 ;;
esac
