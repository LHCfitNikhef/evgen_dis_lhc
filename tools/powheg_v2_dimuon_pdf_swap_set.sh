#!/usr/bin/env bash
# Re-reweight ONE PDF set on every point of the POWHEG-V2 dimuon ladder and
# splice its columns into the existing <point>/pdf/weights.npz, leaving the
# other ~280 weights alone.  Written 2026-10-05 when ABMP16 was swapped for
# its alpha_s = 0.118 variant (user): the set takes the old one's weight ids
# in powheg/reweight/make_rwl_dimuon.py, so the analysis reads it unchanged.
#
# Needs what tools/powheg_v2_dimuon_ladder_pdf.sh left in <point>/pdf/
# (pwgevents-subset.lhe, grids, weights.npz).  The nominal weight 1001 is
# recomputed with the set and compared with the stored one before anything
# is overwritten (harvest_weights.py --replace).
#
# Usage: tools/powheg_v2_dimuon_pdf_swap_set.sh <SET as in rwl_dimuon_p.json> [E1 E2 ...]
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
SET=${1:?usage: $0 <SET> [E ...]}; shift
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-dimuon
RW=$BENCH_REPO/powheg/reweight
ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
MAXJOBS=${MAXJOBS:-8}
PERBATCH=${PERBATCH:-6}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }

# the set's weight lines (plus the nominal 1001) as a small rwl file per target
for t in p n; do
  $BENCH_PYTHON - "$RW/rwl_dimuon_$t.xml" "$RW/rwl_dimuon_$t.json" "$SET" > "$TMPDIR/swap_rwl_$t.xml" <<'PY'
import json, re, sys
xml, js, name = sys.argv[1:4]
ids = set(json.load(open(js))["sets"][name]["weight_ids"]) | {"1001"}
w = [l for l in open(xml) if (m := re.search(r"<weight id='(\d+)'>", l)) and m.group(1) in ids]
print("<initrwgt>\n<weightgroup name='swap' combine='none'>\n" + "".join(w) + "</weightgroup>\n</initrwgt>")
PY
done

one () {   # <point dir> <p|n>
  local R=$1 t=$2 D=$1/pdf S=$1/pdf/swap
  [ -s "$D/weights.npz" ] && [ -s "$D/pwgevents-subset.lhe" ] || { echo "incomplete point: $R"; return 0; }
  rm -rf "$S"; mkdir -p "$S/wgt"
  cp "$D"/pwggrid*.dat "$D"/pwgubound*.dat "$D"/pwgxgrid.dat "$S/" 2>/dev/null || true
  local NB b x
  NB=$($BENCH_PYTHON "$RW/split_rwl.py" "$TMPDIR/swap_rwl_$t.xml" "$S" "$PERBATCH")
  for b in $(seq 1 "$NB"); do
    x=$(printf "rwl_b%03d.xml" "$b")
    cat "$R/powheg.input" > "$S/powheg.input"
    printf "rwl_file '%s'\nrwl_add 1\nrwl_group_events 2000\nrwl_format_rwgt 1\n" "$x" >> "$S/powheg.input"
    rm -f "$S/pwgevents-rwgt.lhe"
    cp "$D/pwgevents-subset.lhe" "$S/pwgevents.lhe"
    ( cd "$S" && bench_run "$PWHG" > "reweight_b$b.log" 2>&1 ) || { echo "FAILED pass $b: $S"; return 0; }
    [ -s "$S/pwgevents-rwgt.lhe" ] || { echo "FAILED pass $b produced nothing: $S"; return 0; }
    $BENCH_PYTHON "$RW/harvest_weights.py" "$S/pwgevents-rwgt.lhe" "$S/wgt/$(printf b%03d.npz "$b")" > /dev/null
  done
  $BENCH_PYTHON "$RW/harvest_weights.py" --merge "$S/wgt" "$S/new.npz" > /dev/null
  if $BENCH_PYTHON "$RW/harvest_weights.py" --replace "$D/weights.npz" "$S/new.npz" > "$S/replace.log" 2>&1; then
    echo "done: $(basename "$R") $(cat "$S/replace.log")"; rm -rf "$S"
  else
    echo "REFUSED: $(basename "$R") $(cat "$S/replace.log")"
  fi
}

while read -r c t E; do
  limit
  one "$BASEDIR/${c}_${t}_E$E" "$t" &
done < <(for c in nu nubar; do for t in p n; do for E in $ES; do echo "$c $t $E"; done; done; done)
wait
echo "swap of $SET finished"
