#!/usr/bin/env bash
# GENIE ANTINEUTRINO samples for the SIDIS yields (paper plot 15), 2026-10-04:
# the default tune (G18_02a, to its 1 TeV validity) and HEDIS, on p and n,
# at the SIDIS ladder energies, 4 x 25k events per point.  Splines first
# (genie_job.sh spline), then events (genie_job.sh run).  events.hepmc is kept
# until analysis/faser_pions.py has read it (tools/nubar_sidis.py prunes).
# Usage: nohup /bin/bash tools/genie_nubar_sidis.sh > logs/genie_nubar_sidis.log 2>&1 &
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
J=$HERE/../genie/production/genie_job.sh
LOG=$HERE/../logs/genie_job; mkdir -p "$LOG"
echo "GENIE nubar start $(date)"
# splines: all six in parallel (each is one gmkspl)
for t in p n; do
    ( "$J" spline nubar_grv "$t" 2000 >> "$LOG/nubar_grv_$t.log" 2>&1 ) &
    ( "$J" spline nubar_hedis "$t" 2000 >> "$LOG/nubar_hedis_$t.log" 2>&1 ) &
    ( "$J" spline nubar_hedis "$t" 8000 >> "$LOG/nubar_hedis_$t.log" 2>&1 ) &
done
wait
echo "splines done $(date)"
lane () {   # cfg energies
    local c=$1 e t; shift
    for e in "$@"; do
        for t in p n; do
            ( "$J" run "$c" "$t" "$e" 25000 4 >> "$LOG/${c}_$t.log" 2>&1 \
              || echo "  FAILED $c $t $e" ) &
        done
        wait
        echo "  $(date +%H:%M) $c $e GeV done"
    done
}
lane nubar_grv 1000 400 700 300 &
lane nubar_hedis 1000 400 4000 700 2000 300 &
wait
echo "GENIE nubar end $(date)"
