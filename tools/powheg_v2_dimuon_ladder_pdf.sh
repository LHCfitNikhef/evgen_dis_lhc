#!/usr/bin/env bash
# PDF-member (and scale) weights on every point of the POWHEG-V2 dimuon
# ladder (user, 2026-09-07: the dimuon cut flow for NNPDF4.0, CT18, MSHT20,
# ATLASpdf21, ABMP16 with their PDF bands, and GRV98 without one).
#
# The exact-subset trick of powheg/reweight/run_reweight_selected.sh: only
# the events that can enter a cut-flow row are reweighted -- the charm events
# for every row from "with a charm hadron" down, and a fixed random sample
# for the CC row (analysis/faser_dimuon_cutflow.py subset) -- which is the
# same calculation as reweighting the whole file, since an event outside a
# row contributes nothing to that row's member spread.  ~22k of 100k events
# per point, 316 weights, a few minutes per point.
#
# Proton points read powheg/reweight/rwl_dimuon_p.xml, neutron points the
# _n twin, whose identical weight ids map to the isospin-mirrored sets of
# tools/make_neutron_pdf.py --all-members (make_rwl_dimuon.py).  The pristine
# pwgevents.lhe is never touched: everything happens in <point>/pdf/.
#
# Usage: tools/powheg_v2_dimuon_ladder_pdf.sh [E1 E2 ...]     (MAXJOBS=4)
# Output: $POWHEG_V2/ladder-dimuon/<point>/pdf/{weights.npz,
#         pwgevents-subset.index.json, subset.json}
#
# COST, MEASURED: 0.7 ms per (event, weight) whatever the batching, so
# ~5k events x 316 weights is ~18 min a point and the ladder is 1.5 h on
# 8 cores.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-dimuon
RW=$BENCH_REPO/powheg/reweight
ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
MAXJOBS=${MAXJOBS:-8}
PERBATCH=${PERBATCH:-6}      # ~15 members in one pass aborts pwhg_main (SIGABRT); 6 is known to survive
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
for t in p n; do [ -s "$RW/rwl_dimuon_$t.xml" ] || { echo "no $RW/rwl_dimuon_$t.xml -- run make_rwl_dimuon.py" >&2; exit 1; }; done

one () {   # <point dir> <p|n>
  local R=$1 t=$2 D=$1/pdf
  if [ -s "$D/weights.npz" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $D/weights.npz"; return 0; fi
  [ -s "$R/pwgevents.lhe" ] && [ -s "$R/dimuon_events.npz" ] || { echo "incomplete point: $R"; return 0; }
  rm -rf "$D"; mkdir -p "$D/wgt"
  cp "$R"/pwggrid*.dat "$R"/pwgubound*.dat "$R"/pwgxgrid.dat "$D/" 2>/dev/null || true
  python3 "$BENCH_REPO/analysis/faser_dimuon_cutflow.py" subset "$R/dimuon_events.npz" "$D/subset.json" > "$D/subset.log" 2>&1
  python3 "$RW/subset_lhe.py" "$R/pwgevents.lhe" "$D/pwgevents-subset.lhe" "$D/subset.json" >> "$D/subset.log" 2>&1
  local NB
  NB=$(python3 "$RW/split_rwl.py" "$RW/rwl_dimuon_$t.xml" "$D" "$PERBATCH")
  local b x npz
  for b in $(seq 1 "$NB"); do
    x=$(printf "rwl_b%03d.xml" "$b")
    npz=$D/wgt/$(printf "b%03d.npz" "$b")
    [ -f "$npz" ] && continue
    cat "$R/powheg.input" > "$D/powheg.input"
    cat >> "$D/powheg.input" <<KEYS
! ---- a-posteriori reweighting of the dimuon subset, pass $b of $NB ----
rwl_file '$x'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
    rm -f "$D/pwgevents-rwgt.lhe"
    cp "$D/pwgevents-subset.lhe" "$D/pwgevents.lhe"
    ( cd "$D" && bench_run "$PWHG" > "reweight_b$b.log" 2>&1 ) || { echo "FAILED pass $b: $D (see reweight_b$b.log)"; return 0; }
    [ -s "$D/pwgevents-rwgt.lhe" ] || { echo "FAILED pass $b produced nothing: $D"; return 0; }
    python3 "$RW/harvest_weights.py" "$D/pwgevents-rwgt.lhe" "$npz" > /dev/null
    rm -f "$D/pwgevents-rwgt.lhe"
  done
  python3 "$RW/harvest_weights.py" --merge "$D/wgt" "$D/weights.npz" > /dev/null
  rm -f "$D/pwgevents.lhe" "$D"/rwl_b*.xml
  rm -rf "$D/wgt"
  echo "done: $(basename "$R")  $(python3 -c "import numpy as np; z=np.load('$D/weights.npz'); print(z['weights'].shape)")"
}

# NOT `... | while`: a piped loop's background jobs escape `wait`
while read -r c t E; do
  R=$BASEDIR/${c}_${t}_E$E
  limit
  one "$R" "$t" &
done < <(for c in nu nubar; do for t in p n; do for E in $ES; do echo "$c $t $E"; done; done; done)
wait
echo "dimuon ladder PDF reweighting finished"
