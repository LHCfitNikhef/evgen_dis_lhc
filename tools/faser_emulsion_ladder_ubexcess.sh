#!/usr/bin/env bash
# THE pp12 POWHEG-V2 LADDER REGENERATED WITH ubexcess_correct 1 (2026-09-19,
# user: regenerate POWHEG-V2 with it -- see powheg/production/regen_ubexcess.sh and
# PRODUCTION.md, "POWHEG-V2 low-x deficit").
#
#   tools/faser_emulsion_ladder_ubexcess.sh [max_parallel=8]
#   then: python3 analysis/faser_emulsion_shapes.py combine
#
# The ladder of tools/faser_emulsion_ladder.sh has 40 POWHEG points: nu and
# nubar, p and n, 10 energies.  Its NEUTRON points were generated in
# $POWHEG_V2/ladder-faser-v2/<c>_n_E<E>; its PROTON points were SYMLINKS to the
# earlier LHE in $POWHEG_V2/ladder-dimuon/<c>_p_E<E>, which the FASER dimuon cut
# flow also reads.  Here, per point, with the point's OWN integration grids
# and upper bounds (never re-integrated -- refused if they change) and its OWN
# iseed:
#   n: regenerated in place with `ubexcess_correct 1`;
#   p: the earlier point's grids, bounds, pwg-stat.dat and card are COPIED into
#      ladder-faser-v2/<c>_p_E<E> (replacing the symlinks) and regenerated
#      there -- ladder-dimuon, and with it the dimuon cut flow, is untouched;
# then showered with the card (matchInOut off) and extracted exactly as
# faser_emulsion_ladder.sh powheg-shower does.  A done point carries
# UBEXCESS_OK next to a fresh faserdata_events.npz; safe to re-invoke.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PY=${BENCH_PYTHON:-python3}
MAXJOBS=${1:-8}
V1=$POWHEG_V2/ladder-dimuon
V2L=$POWHEG_V2/ladder-faser-v2
CMND=$BENCH_REPO/powheg/powheg_nu_v2.cmnd
PWHG=$POWHEG_V2/pwhg_main
ES="30 60 100 200 400 700 1000 2000 4000 6800"
GRIDS="pwggrid.dat pwgubound.dat pwgborngrid.top pwgxgrid.dat FlavRegList bornequiv virtequiv realequivregions-btl realequivregions-rad pwg-btlgrid.top pwg-rmngrid.top pwg-stat.dat"
grep -q '^LesHouches:matchInOut = off$' "$CMND" || { echo "$CMND lost matchInOut = off" >&2; exit 1; }
[ "$(df -Pk "$PHYSICS24" | awk 'NR==2 {print int($4/1048576)}')" -ge 158 ] || { echo "disk floor" >&2; exit 1; }
md5f () { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }

point () {   # c t E   (foreground)
  local c=$1 t=$2 E=$3 ih1=14 R want g0 u0 f
  [ "$c" = nubar ] && ih1=-14
  R=$V2L/${c}_${t}_E$E
  want="$ih1 2212"; [ "$t" = n ] && want="$ih1 2112"
  if [ -f "$R/UBEXCESS_OK" ] && [ -s "$R/faserdata_events.npz" ] && [ "$R/faserdata_events.npz" -nt "$R/UBEXCESS_OK" ]; then
    echo "skip ${c}_${t}_E$E"; return 0; fi
  if [ "$t" = p ] && [ ! -f "$R/UBEXCESS_OK" ]; then
    [ -s "$V1/${c}_p_E$E/pwggrid.dat" ] || { echo "no dimuon-ladder grids for ${c}_p_E$E" >&2; return 1; }
    rm -f "$R/pwgevents.lhe" "$R/pwg-stat.dat"          # the symlinks to the dimuon ladder
    mkdir -p "$R"
    for f in $GRIDS; do [ -f "$V1/${c}_p_E$E/$f" ] && cp -f "$V1/${c}_p_E$E/$f" "$R/$f"; done
    cp -f "$V1/${c}_p_E$E/powheg.input" "$R/powheg.input"
  fi
  if [ ! -f "$R/UBEXCESS_OK" ]; then
    [ -s "$R/pwggrid.dat" ] && [ -s "$R/pwgubound.dat" ] || { echo "no grids in $R" >&2; return 1; }
    grep -q '^use-old-grid *1' "$R/powheg.input" && grep -q '^use-old-ubound *1' "$R/powheg.input" || {
      echo "$R/powheg.input does not reuse grid and bound" >&2; return 1; }
    g0=$(md5f "$R/pwggrid.dat"); u0=$(md5f "$R/pwgubound.dat")
    grep -q '^ubexcess_correct' "$R/powheg.input" || \
      printf 'ubexcess_correct 1  ! faser_emulsion_ladder_ubexcess.sh, 2026-09-19\n' >> "$R/powheg.input"
    rm -f "$R/pwgevents.lhe" "$R/pwgevents.lhe.nevents" "$R/faserdata_events.npz"
    ( cd "$R" && bench_run "$PWHG" > run_ubexcess.log 2>&1 ) || true
    tail -c 200 "$R/pwgevents.lhe" 2>/dev/null | grep -q "</LesHouchesEvents>" || { echo "gen FAILED ${c}_${t}_E$E" >&2; return 1; }
    [ "$(md5f "$R/pwggrid.dat")" = "$g0" ] && [ "$(md5f "$R/pwgubound.dat")" = "$u0" ] || {
      echo "${c}_${t}_E$E: grids CHANGED (re-integrated) -- refusing" >&2; return 1; }
    "$BENCH_REPO/powheg/strip_lhe_nan.py" "$R/pwgevents.lhe" > "$R/strip_nan.log" 2>&1
    [ "$(awk '/<init>/ {getline; print $1, $2; exit}' "$R/pwgevents.lhe")" = "$want" ] || {
      echo "${c}_${t}_E$E: LHE beams are not '$want'" >&2; return 1; }
    date '+%F %T' > "$R/UBEXCESS_OK"
  fi
  ( cd "$BENCH_REPO/powheg" \
      && bench_run ./main_powheg "$CMND" "$R/pwgevents.lhe" "$R/fd_events" > "$R/shower_faserdata.log" 2>&1 \
      && "$PY" "$BENCH_REPO/analysis/faser_emulsion_shapes.py" extract \
           --beam="$ih1" --target="$t" --energy="$E" --generator=powheg_v2 \
           "$R/fd_events.hepmc" "$R/faserdata_events.npz" > "$R/faserdata_extract.log" 2>&1 \
      && [ -s "$R/faserdata_events.npz" ] \
      && rm -f "$R/fd_events.hepmc" ) || { echo "shower/extract FAILED ${c}_${t}_E$E" >&2; return 1; }
  echo "done ${c}_${t}_E$E: $(tail -1 "$R/faserdata_extract.log")"
}

if [ "${2:-}" = "--one" ]; then point "$3" "$4" "$5"; exit $?; fi
for c in nu nubar; do for t in p n; do for E in $ES; do
  limit
  point "$c" "$t" "$E" &
done; done; done
wait
echo "== ladder regeneration finished $(date '+%F %T')"
