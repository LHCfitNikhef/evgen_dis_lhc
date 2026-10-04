#!/usr/bin/env bash
# THE  LADDERS FOR PAPER PLOT 12, the comparison with FASER's emulsion data.
#
# Usage: tools/faser_emulsion_ladder.sh <stage> [E1 E2 ...]
#   powheg-gen     POWHEG-V2 NEUTRON points, genuine neutron beam (ih2 = 2)
#   powheg-shower  shower + extract ALL POWHEG points (p from the kept earlier LHE,
#                  n from powheg-gen) with the card, matchInOut off
#   sherpa-n       Sherpa MC@NLO NEUTRON points, genuine 2112 beam, patched copy
#   status
# Knobs: MAXJOBS (default 3), NEV (POWHEG, 100000), NEV_SHERPA (20000).
#
# WHAT CHANGES FROM the earlier production, AND WHY (user, 2026-09-14).  pp12 stays at FASER's own
# cuts -- no Q2 window beyond the ladders' generation floor, no y, the lepton
# among the tracks -- so its region is not the benchmark region and nothing
# is re-cut.  Two things were earlier in the samples themselves:
#
#   1. POWHEG-V2 was showered before the charm fix.  Pythia dropped ~2% of the
#      events, nearly all charm ("setting mass failed"); powheg_nu_v2.cmnd sets
#      LesHouches:matchInOut = off.  Every point is re-showered with it.
#   2. THE NEUTRON POINTS WERE PROTONS WITH A NEUTRON PDF.  earlier ran POWHEG-V2
#      and Sherpa on a 2212 beam carrying NNPDF40_nnlo_as_01180_n: the right
#      cross-section and ISR, but a PROTON remnant (uud), one unit of charge
#      off in the hadronic system -- and the track multiplicity is one of the
#      four observables.  User: "Regenerate n points".  Here the neutron is a
#      genuine neutron everywhere, as in the production:
#        POWHEG-V2  ih2 = 2 with the PROTON set (it swaps u<->d itself), LHE
#                   beam 2112; Pythia swaps the proton pSet for 2112 itself.
#        Sherpa     2112 beam with the neutron set and PDF_LIBRARY, run from
#                   the patched copy $SHERPA_NPATCH, library and remnant checked
#                   per point (tools/sherpa_check_events.py).
#      The proton points' LHE and Sherpa runs are unchanged (they were right).
#      GENIE's ladder always used the genuine neutron target code.
#
# OUTPUT.  POWHEG: $POWHEG_V2/ladder-faser-v2/<nu|nubar>_<p|n>_E<E>/ with
# faserdata_events.npz and pwg-stat.dat (p: a link to the earlier point's).  Sherpa:
# $SHERPA_RUNS/V2_NuDIS_NLO_FASER_<nu|nubar>_n_E<E>/.  The earlier directories are
# not written.  analysis/faser_emulsion_shapes.py reads these for .
#
# Safe to kill and re-invoke: finished points are skipped.  CONVENTIONS.md 1c: do
# not edit while running.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PY=${BENCH_PYTHON:-python3}
STAGE=${1:-status}
shift || true
MAXJOBS=${MAXJOBS:-3}
NEV=${NEV:-100000}
NEV_SHERPA=${NEV_SHERPA:-20000}
V1=$POWHEG_V2/ladder-dimuon
V2L=$POWHEG_V2/ladder-faser-v2
CARD=$BENCH_REPO/powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input
CMND=$BENCH_REPO/powheg/powheg_nu_v2.cmnd
PWHG=$POWHEG_V2/pwhg_main
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 10; done; }
free_gb () { df -Pk "$PHYSICS24" | awk 'NR==2 {print int($4/1048576)}'; }
disk_guard () {
  if [ "$(free_gb)" -lt 158 ]; then
    echo "only $(free_gb) GB free; the floor is 150 GB (CONVENTIONS.md rule 1b) -- stopping" >&2
    exit 1
  fi
}
grep -q '^LesHouches:matchInOut = off$' "$CMND" || { echo "$CMND lost matchInOut = off" >&2; exit 1; }
grep -q '^PDF:pSet = LHAPDF6:NNPDF40_nnlo_as_01180$' "$CMND" || { echo "$CMND lost the proton pSet" >&2; exit 1; }

PW_ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
SH_ES=${*:-"200 400 700 1000 2000 4000"}

powheg_gen () {
  disk_guard
  mkdir -p "$V2L"
  for c in nu nubar; do for E in $PW_ES; do
    local ih1=14; [ "$c" = nubar ] && ih1=-14
    local R=$V2L/${c}_n_E$E
    if [ -s "$R/pwgevents.lhe.nevents" ]; then echo "exists: ${c}_n_E$E"; continue; fi
    mkdir -p "$R"
    "$PY" - "$CARD" "$R/powheg.input" "$E" "$ih1" "$NEV" <<'PY'
import re, sys, os, zlib
src, dst, e, ih1, nev = sys.argv[1:6]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "analysis"))
import beams
eb = beams.Beams("nu", float(e)).sqrt_s / 2.0     # s on M_P for both nucleons
s = open(src).read()
def sub1(pat, rep):
    global s
    s, n = re.subn(pat, rep, s, count=1, flags=re.M)
    assert n == 1, pat
sub1(r"^numevts .*", f"numevts {nev}")
sub1(r"^ih1 .*", f"ih1   {ih1}")
sub1(r"^ih2 .*", "ih2   2            ! GENUINE NEUTRON (PROTON set; POWHEG swaps u<->d)")
assert re.search(r"^lhans1 331100", s, re.M) and re.search(r"^lhans2 331100", s, re.M)
sub1(r"^ebeam1 .*", f"ebeam1 {eb:.14f}d0")
sub1(r"^ebeam2 .*", f"ebeam2 {eb:.14f}d0")
seed = zlib.crc32(f"v2n:{float(e):g}:{ih1}".encode()) % 900000 + 100
sub1(r"^iseed .*", f"iseed {seed}")
open(dst, "w").write(s)
PY
    limit
    ( cd "$R" && bench_run "$PWHG" > run_v2.log 2>&1 \
        && "$BENCH_REPO/powheg/strip_lhe_nan.py" pwgevents.lhe > strip_nan.log 2>&1 \
        && b=$(awk '/<init>/ {getline; print $1, $2; exit}' pwgevents.lhe) \
        && { [ "$b" = "$ih1 2112" ] || { echo "BAD BEAMS '$b' in $R" >&2; rm -f pwgevents.lhe.nevents; false; }; } \
        && echo "gen done: ${c}_n_E$E  $(grep -m1 'total (btilde+remnants) cross section' pwg-stat.dat)" \
        || echo "gen FAILED: $R" >&2 ) &
  done; done
  wait
}

powheg_shower () {
  for c in nu nubar; do for t in p n; do for E in $PW_ES; do
    local ih1=14; [ "$c" = nubar ] && ih1=-14
    local R=$V2L/${c}_${t}_E$E lhe want
    mkdir -p "$R"
    if [ "$t" = p ]; then
      [ -s "$V1/${c}_p_E$E/pwgevents.lhe" ] || { echo "no dimuon-ladder proton LHE for ${c}_p_E$E" >&2; continue; }
      ln -sfn "$V1/${c}_p_E$E/pwgevents.lhe" "$R/pwgevents.lhe"
      ln -sfn "$V1/${c}_p_E$E/pwg-stat.dat" "$R/pwg-stat.dat"
      want="$ih1 2212"
    else
      want="$ih1 2112"
    fi
    lhe=$R/pwgevents.lhe
    [ -s "$lhe" ] || { echo "no LHE yet: ${c}_${t}_E$E"; continue; }
    if [ -s "$R/faserdata_events.npz" ]; then echo "done: ${c}_${t}_E$E"; continue; fi
    [ "$(awk '/<init>/ {getline; print $1, $2; exit}' "$lhe")" = "$want" ] || {
        echo "LHE beams of ${c}_${t}_E$E are not '$want' -- refusing" >&2; continue; }
    disk_guard
    limit
    ( cd "$BENCH_REPO/powheg" \
        && bench_run ./main_powheg "$CMND" "$lhe" "$R/fd_events" > "$R/shower_faserdata.log" 2>&1 \
        && "$PY" "$BENCH_REPO/analysis/faser_emulsion_shapes.py" extract \
             --beam="$ih1" --target="$t" --energy="$E" --generator=powheg_v2 \
             "$R/fd_events.hepmc" "$R/faserdata_events.npz" > "$R/faserdata_extract.log" 2>&1 \
        && [ -s "$R/faserdata_events.npz" ] \
        && rm -f "$R/fd_events.hepmc" \
        && echo "powheg ${c}_${t}_E$E: $(tail -1 "$R/faserdata_extract.log")" \
        || echo "powheg shower FAILED: ${c}_${t}_E$E" >&2 ) &
  done; done; done
  wait
}

# ---------------------------------------------------------------- Sherpa
SRC=$SHERPA_RUNS/NuDIS_NLO_ckm3
NUBAR_SEED=$SHERPA_RUNS/NuDIS_NLO_FASER_nubar_p_E1000

bvi_rs () {                       # <integ.log> -> "BVI RS"
  sed 's/\x1b\[[0-9;]*m//g' "$1" | sed -nE 's/.*\((BVI|RS)\) : ([-0-9.e+]+) pb.*/\1 \2/p' \
      | awk '{v[$1]=$2} END {print v["BVI"], v["RS"]}'
}
same_to_1e5 () {                  # "<BVI> <RS>" "<BVI> <RS>"
  awk -v a="$1" -v b="$2" 'BEGIN {split(a, x, " "); split(b, y, " ");
      for (i = 1; i <= 2; i++) { if (x[i] == "" || y[i] == "") exit 1;
      r = (x[i] - y[i]) / x[i]; if (r < 0) r = -r; if (r > 1e-5) exit 1 } }'
}

sherpa_point () {    # <nu|nubar> <E>   (foreground)
  local c=$1 E=$2 pid=14 lep=13
  [ "$c" = nubar ] && { pid=-14; lep=-13; }
  local d=$SHERPA_RUNS/V2_NuDIS_NLO_FASER_${c}_n_E$E
  if [ -s "$d/faserdata_events.npz" ] && grep -q "^OK" "$d/remnant_check.txt" 2>/dev/null; then
    echo "done: sherpa ${c}_n_E$E"; return 0; fi
  mkdir -p "$d"
  local libs=$SRC/Process
  [ "$c" = nubar ] && libs=$NUBAR_SEED/Process
  [ -d "$libs/Amegic/lib" ] || [ -d "$libs" ] || { echo "no Process libraries at $libs" >&2; return 1; }
  [ -e "$d/Process" ] || ln -s "$libs" "$d/Process"
  # the earlier ladder's card, with the neutron beam of tools/make_sherpa_cards.py
  "$PY" - "$SHERPA_RUNS/NuDIS_NLO_FASER_${c}_p_E$E/Sherpa.yaml" "$d/Sherpa.yaml" "$E" <<'PY'
import os, re, sys
src, dst, e = sys.argv[1:4]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "tools"))
import make_sherpa_cards as mk
ea, eb = mk.cm_energies("nu", "n", float(e), neutron_beam=True)
t = open(src).read()
def sub1(pat, rep):
    global t
    t, n = re.subn(pat, rep, t, count=1, flags=re.M)
    assert n == 1, pat
m = re.search(r"^BEAMS: \[(-?\d+), 2212\]$", t, re.M)
assert m, "the proton card has no [lepton, 2212] beam line"
sub1(r"^BEAMS: .*$", f"BEAMS: [{m.group(1)}, 2112]")
sub1(r"^BEAM_ENERGIES: .*$", f"BEAM_ENERGIES: [{ea:.7f}, {eb:.7f}]")
sub1(r"^PDF_SET: .*$", "PDF_LIBRARY: [None, LHAPDFSherpa]\nPDF_SET: [None, NNPDF40_nnlo_as_01180_n]")
sub1(r"^MPI_PDF_SET: .*$", "MPI_PDF_SET: NNPDF40_nnlo_as_01180_n")
assert "INEL" not in t
t = ("# v2 NEUTRON POINT of the FASER emulsion ladder -- GENERATED by\n"
     "# tools/faser_emulsion_ladder.sh from the v1 proton point's card: a\n"
     "# GENUINE 2112 beam (run only from the patched copy $SHERPA_NPATCH).\n" + t)
open(dst, "w").write(t)
PY
  local SH=$SHERPA_NPATCH/bin/Sherpa v1n=$SHERPA_RUNS/NuDIS_NLO_FASER_${c}_n_E$E
  # REUSE the earlier point's grid (2212 beam + neutron set, same region): the
  # 2112 card accepts it and gives the same cross-section, as the 
  # production measured (tools/sherpa_production.sh).  Accepted only if the 2112 run
  # does not re-integrate and repeats BVI and RS to 1e-5; else integrate afresh.
  if [ ! -s "$d/Results.zip" ] && [ -s "$v1n/Results.zip" ]; then
    cp -f "$v1n/Results.zip" "$d/Results.zip"; touch "$d/REUSE_ATTEMPT"
  fi
  ( cd "$d" \
    && export SHERPA_LIBRARY_PATH="$SHERPA_NPATCH/lib/SHERPA-MC/" SHERPA_SHARE_PATH="$SHERPA_NPATCH/share/SHERPA-MC/" \
    && { grep -q "^Time: " integ.log 2>/dev/null || bench_run "$SH" -e 0 > integ.log 2>&1; } \
    && { [ ! -f REUSE_ATTEMPT ] || {
           a=$(bvi_rs "$v1n/integ.log"); b=$(bvi_rs integ.log)
           if ! grep -q "Starting the calculation" integ.log && same_to_1e5 "$a" "$b"; then
             echo "reused the dimuon-ladder grid: BVI RS $a -> $b" > RESULTS_REUSED
           else
             echo "reuse refused ($a -> $b), integrating afresh"
             mv integ.log integ_reuse_refused.log; rm -f Results.zip Results.zip~
             bench_run "$SH" -e 0 > integ.log 2>&1
           fi
           rm -f REUSE_ATTEMPT; }; } \
    && grep -q "Beam 2: n " <(sed 's/\x1b\[[0-9;]*m//g' integ.log) \
    && bench_run "$SH" -e "$NEV_SHERPA" \
         'EVENT_OUTPUT: ["HepMC3_GenEvent[evtfull]"]' > gen.log 2>&1 \
    && "$PY" "$BENCH_REPO/tools/sherpa_check_events.py" evtfull --target n --beam "$pid" \
         > remnant_check.txt 2>&1 \
    && grep -q "^OK" remnant_check.txt \
    && "$PY" "$BENCH_REPO/analysis/faser_emulsion_shapes.py" extract \
         --beam="$pid" --target=n --energy="$E" --generator=sherpa_nlo \
         evtfull "$d/faserdata_events.npz" > extract.log 2>&1 \
    && [ -s faserdata_events.npz ] \
    && rm -f evtfull \
    && echo "sherpa ${c}_n_E$E: $(head -1 remnant_check.txt | cut -c1-80); $(tail -1 extract.log)" \
    || echo "sherpa FAILED: ${c}_n_E$E (see $d)" >&2 )
}

sherpa_n () {
  [ -x "$SHERPA_NPATCH/bin/Sherpa" ] && [ -f "$SHERPA_NPATCH/NPATCH_STAMP" ] || {
      echo "no patched Sherpa copy at $SHERPA_NPATCH" >&2; exit 1; }
  [ -d "$NUBAR_SEED/Process" ] || { echo "no antineutrino libraries in $NUBAR_SEED" >&2; exit 1; }
  for c in nu nubar; do for E in $SH_ES; do
    disk_guard
    limit
    sherpa_point "$c" "$E" &
  done; done
  wait
}

status () {
  local c t E n=0 m=0
  for c in nu nubar; do for t in p n; do for E in $PW_ES; do
    [ -s "$V2L/${c}_${t}_E$E/faserdata_events.npz" ] && n=$((n + 1))
  done; done; done
  for c in nu nubar; do for E in $SH_ES; do
    [ -s "$SHERPA_RUNS/V2_NuDIS_NLO_FASER_${c}_n_E$E/faserdata_events.npz" ] && m=$((m + 1))
  done; done
  echo "POWHEG points extracted: $n / 40;  Sherpa neutron points: $m / 12"
}

case "$STAGE" in
  powheg-gen)    powheg_gen ;;
  powheg-shower) powheg_shower ;;
  sherpa-n)      sherpa_n ;;
  status)        status ;;
  *) sed -n 2,12p "$0"; exit 2 ;;
esac
