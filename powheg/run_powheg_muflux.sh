#!/usr/bin/env bash
# POWHEG-RES muon NC DIS on the FASERnu MUON FLUX: the arXiv:2506.13889 setup.
#
# Usage: ./run_powheg_muflux.sh [nseeds]        (default 8; FORCE=1 to redo)
#
# WHY THIS IS A SEPARATE DRIVER.  run_powheg_res.sh runs the benchmark at a
# FIXED beam energy through the variant table; this run draws the incoming
# muon from the flux grid instead (fixed_lepton_beam 0 + nupdf in the card),
# so nothing in it is a benchmark point and it must not be confused with one.
# It exists to check that, under the paper's own cuts, the generator chain
# reproduces the paper's Table 2.1 (user, 2026-09-03: "an important check").
#
# THE FLUX GRID IS FOUND THROUGH LHAPDF_DATA_PATH, which is PREPENDED with
# the in-repo layout data/faser_muon_flux/lhapdf/ (a set directory plus a
# pdfsets.index naming its LHAPDF id).  Prepended, not appended: an index
# collision resolves to whichever pdfsets.index is read first, and the
# vendored id 111222333 was chosen by the flux authors, not registered.
#
# THE BINARY IS A SEPARATE COPY OF THE PROCESS, DIS_v_muflux, carrying the
# paper's lepton_flux.f (patches/powheg-res-dis-v-muflux-lepton-flux.diff):
# the stock DIS_v's flux route zeroes the beam-1 quark entries and POWHEG's
# limit check then loops forever.  The four stages are as in
# run_powheg_res.sh; the output lands in $POWHEG_RES_MUFLUX/parallel-muflux-2506.13889/; shower and analyse with
#     powheg/run_powheg_muflux.sh --shower   ->  powheg/muflux_job_N/events.hepmc
#     analysis/faser_muflux_check.py         ->  results/faser_muflux_check.json
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/powheg_variants.sh"        # powheg_install_card

export LHAPDF_DATA_PATH=$BENCH_REPO/data/faser_muon_flux/lhapdf:$LHAPDF_DATA_PATH
CARD=$HERE/cards/POWHEG-RES/muflux-2506.13889/powheg.input
RUNDIR=$POWHEG_RES_MUFLUX/parallel-muflux-2506.13889
PWHG=$POWHEG_RES_MUFLUX/pwhg_main   # the lepton-flux copy, see the card header
[ -x "$PWHG" ] || { echo "no POWHEG binary at $PWHG" >&2; exit 1; }

if [ "${1:-}" = "--shower" ]; then
    # Pythia8 on each LHE, exactly as run_powheg_shower.sh does for the
    # benchmark, but with the lepton "PDF" switched on so that Pythia accepts
    # an incoming muon carrying x < 1 of the 7 TeV beam declared in the LHE.
    CMND=$HERE/powheg_muflux.cmnd
    for f in "$RUNDIR"/pwgevents-*.lhe; do
        s=$(basename "$f" .lhe); s=${s#pwgevents-}; s=$((10#$s))
        job=$HERE/muflux_job_$s
        if [ -s "$job/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
            echo "  $job/events.hepmc exists -- FORCE=1 to redo"; continue
        fi
        mkdir -p "$job"
        ( cd "$job" && bench_run "$HERE/main_powheg" "$CMND" "$f" events \
              > shower.log 2>&1 ) &
    done
    wait
    echo "showered: $(ls -d "$HERE"/muflux_job_* | wc -l | tr -d ' ') job dir(s)"
    exit 0
fi

NSEEDS=${1:-8}
SEEDLIST=($(seq 1 "$NSEEDS"))
for s in "${SEEDLIST[@]}"; do
    f=$RUNDIR/pwgevents-$(printf "%04d" "$s").lhe
    if [ -s "$f" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "refusing to overwrite $f -- FORCE=1 to replace" >&2; exit 1
    fi
done
mkdir -p "$RUNDIR"
powheg_install_card "$CARD" "$RUNDIR/powheg.input"
if [ "${FORCE:-0}" = "1" ]; then
    for s in "${SEEDLIST[@]}"; do
        rm -f "$RUNDIR/pwgevents-$(printf "%04d" "$s").lhe"
    done
fi
seq 1 "$NSEEDS" > "$RUNDIR/pwgseeds.dat"

set_stage () {   # $1 = parallelstage, $2 = xgriditeration
    python3 - "$RUNDIR/powheg.input" "$1" "$2" <<'PY'
import re, sys
path, stage, xg = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
s = re.sub(r"^parallelstage\s+\d+", f"parallelstage {stage}", s, flags=re.M)
s = re.sub(r"^xgriditeration\s+\d+", f"xgriditeration {xg}", s, flags=re.M)
open(path, "w").write(s)
PY
}
run_pass () {    # $1 = label, $2 = parallelstage, $3 = xgriditeration
    local label=$1
    set_stage "$2" "$3"
    echo "  $label: seeds ${SEEDLIST[*]} ..."
    echo "beg $label $(date)" >> "$RUNDIR/Timings.txt"
    for s in "${SEEDLIST[@]}"; do
        ( cd "$RUNDIR" && echo "$s" | bench_run "$PWHG" \
            > "run_${label}_$s.log" 2>&1 ) &
    done
    wait
    echo "end $label $(date)" >> "$RUNDIR/Timings.txt"
    if grep -qil "error\|not found\|abort" "$RUNDIR"/run_${label}_${SEEDLIST[0]}.log 2>/dev/null; then
        echo "  WARNING: $label log mentions an error -- check" \
             "$RUNDIR/run_${label}_${SEEDLIST[0]}.log" >&2
    fi
}

echo "POWHEG-RES mu NC on the FASERnu muon flux (arXiv:2506.13889 setup), seeds ${SEEDLIST[*]} -> $RUNDIR"
run_pass st1xg1 1 1
run_pass st1xg2 1 2
run_pass st2    2 1
run_pass st3    3 1
run_pass st4    4 1
n=$(ls "$RUNDIR"/pwgevents-*.lhe 2>/dev/null | wc -l | tr -d ' ')
echo "done: $n LHE file(s) in $RUNDIR"
[ "$n" -gt 0 ] || { echo "NO LHE FILES PRODUCED -- check the stage logs" >&2; exit 1; }
echo "next: powheg/run_powheg_muflux.sh --shower"
