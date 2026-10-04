#!/usr/bin/env bash
# POWHEG-V2 ANTINEUTRINO samples for the SIDIS yields (paper plot 15), 2026-10-04.
#
# Reuses the pp12 FASER ladder ($POWHEG_V2/ladder-faser-v2/nubar_<t>_E<E>,
# ubexcess_correct 1, 100k events, genuine neutron for t = n) at 400 GeV-4 TeV
# and generates the missing 300 GeV point from the 400 GeV card (beam energy
# and seed changed, nothing else).  Each point is then showered with the
# production card (powheg_nu_v2.cmnd, matchInOut off) into
#   <point>/sidis_events.hepmc
# which analysis/faser_pions.py reads (key powheg_nubar) and tools/nubar_sidis.py
# deletes after the extraction.
# Usage: nohup /bin/bash tools/powheg_nubar_sidis.sh > logs/powheg_nubar_sidis.log 2>&1 &
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PY=${BENCH_PYTHON:-python3}
V2L=$POWHEG_V2/ladder-faser-v2
CMND=$BENCH_REPO/powheg/powheg_nu_v2.cmnd
PWHG=$POWHEG_V2/pwhg_main
MAXJOBS=${MAXJOBS:-4}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
grep -q '^LesHouches:matchInOut = off$' "$CMND" || { echo "$CMND lost matchInOut = off" >&2; exit 1; }

gen300 () {   # t
  local t=$1
  local R=$V2L/nubar_${t}_E300 src=$V2L/nubar_${t}_E400 want
  want="-14 2212"; [ "$t" = n ] && want="-14 2112"
  if [ -s "$R/pwgevents.lhe.nevents" ]; then echo "exists: nubar_${t}_E300"; return 0; fi
  mkdir -p "$R"
  "$PY" - "$src/powheg.input" "$R/powheg.input" <<'PY'
import re, sys, os, zlib
src, dst = sys.argv[1:3]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "analysis"))
import beams
eb = beams.Beams("nu", 300.0).sqrt_s / 2.0
s = open(src).read()
for pat, rep in ((r"^ebeam1 .*", f"ebeam1 {eb:.14f}d0"), (r"^ebeam2 .*", f"ebeam2 {eb:.14f}d0"),
                 (r"^iseed .*", f"iseed {zlib.crc32(dst.encode()) % 900000 + 100}")):
    s, n = re.subn(pat, rep, s, count=1, flags=re.M)
    assert n == 1, pat
open(dst, "w").write(s)
PY
  ( cd "$R" && bench_run "$PWHG" > run_v2.log 2>&1 ) || true
  tail -c 200 "$R/pwgevents.lhe" 2>/dev/null | grep -q "</LesHouchesEvents>" || { echo "gen FAILED nubar_${t}_E300" >&2; return 1; }
  "$BENCH_REPO/powheg/strip_lhe_nan.py" "$R/pwgevents.lhe" > "$R/strip_nan.log" 2>&1
  [ "$(awk '/<init>/ {getline; print $1, $2; exit}' "$R/pwgevents.lhe")" = "$want" ] || {
      echo "nubar_${t}_E300: LHE beams are not '$want'" >&2; return 1; }
  echo "gen done nubar_${t}_E300: $(grep -m1 'total (btilde+remnants) cross section' "$R/pwg-stat.dat")"
}

shower () {   # t E
  local t=$1 E=$2
  local R=$V2L/nubar_${t}_E$E
  if [ -s "$R/SIDIS_SHOWER_OK" ]; then echo "skip shower nubar_${t}_E$E"; return 0; fi
  [ -s "$R/pwgevents.lhe" ] || { echo "no LHE in $R" >&2; return 1; }
  ( cd "$BENCH_REPO/powheg" \
      && bench_run ./main_powheg "$CMND" "$R/pwgevents.lhe" "$R/sidis_events" > "$R/shower_sidis.log" 2>&1 \
      && tail -c 200 "$R/sidis_events.hepmc" | grep -q "END_EVENT_LISTING" \
      && date '+%F %T' > "$R/SIDIS_SHOWER_OK" \
      && echo "shower done nubar_${t}_E$E" ) || { echo "shower FAILED nubar_${t}_E$E" >&2; return 1; }
}

echo "POWHEG-V2 nubar SIDIS start $(date)"
gen300 p & gen300 n & wait
for E in 1000 400 700 2000 4000 300; do for t in p n; do limit; shower "$t" "$E" & done; done
wait
echo "POWHEG-V2 nubar SIDIS end $(date)"
