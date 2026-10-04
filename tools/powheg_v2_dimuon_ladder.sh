#!/usr/bin/env bash
# POWHEG-V2 ladder for the FASER DIMUON CUT FLOW (user, 2026-09-07): the
# opposite-sign dimuon signal of charm production in nu_mu CC DIS, as a
# FASER collaboration slide estimates it with GENIE at 300 fb^-1, reproduced
# with our own machinery -- see analysis/faser_dimuon_cutflow.py for the cut
# flow itself and the physics.
#
# WHY A SECOND LADDER.  tools/powheg_v2_faser_ladder.sh stops at the Les
# Houches events, which is exact for a cut on Q2 and y but says nothing about
# a muon from a charm-hadron decay.  Every step of the cut flow after the first
# needs the shower, the hadronisation and the decays, so here each ladder
# point is SHOWERED (Pythia 8, the production steering card, so the charm
# fragmentation and the ctau > 10 mm stability convention are exactly the
# benchmark's), and the showered events are reduced to a compact per-event
# table (dimuon_events.npz) that carries every muon's energy, angle and charge.
# The HepMC is then deleted: at 6.4 kB/event, 40 points of 100k events would
# be 26 GB backing a few MB of tables, and the tables hold everything a
# re-cut on the muons can ask for (CONVENTIONS.md rule 1, intermediates do not
# outlive their results).  The Les Houches file is kept, so a re-shower costs
# seconds and a re-extraction never needs pwhg_main again.
#
# Usage: tools/powheg_v2_dimuon_ladder.sh [gen|shower|extract|all] [E1 E2 ...]
#        NEV=100000 MAXJOBS=4 by default.
# Output: $POWHEG_V2/ladder-dimuon/<nu|nubar>_<p|n>_E<GeV>/
#             pwg-stat.dat pwgevents.lhe dimuon_events.npz events_xsec.json
#
# THE NEUTRON POINTS SHOWER WITH THE NEUTRON PDF SET.  The Les Houches init
# block says "proton" (POWHEG's ih2 = 1) with lhans = 339900, the isospin-
# mirrored NNPDF4.0 of tools/make_neutron_pdf.py; the Pythia steering card
# names the proton set for its initial-state shower, so a per-point copy of
# the card swaps in the mirrored set.  A small effect on a decay-muon
# efficiency, but the matrix element and the shower should not disagree
# about which nucleon they are on.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
CARD=$BENCH_REPO/powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input
CMND=$BENCH_REPO/powheg/powheg_nu1TeV.cmnd
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-dimuon
STAGE=${1:-all}
shift || true
ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
MAXJOBS=${MAXJOBS:-4}
NEV=${NEV:-100000}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }

# CONVENTIONS.md rule 1b: 150 GB stay free.  The showers are transient (one
# HepMC per running job, deleted after extraction) but the Les Houches files
# stay: ~0.9 kB per event, 40 points.
free_gb () { df -Pk "$PHYSICS24" | awk 'NR==2 {print int($4/1048576)}'; }   # -k: macOS and Linux alike
if [ "$(free_gb)" -lt 160 ]; then
  echo "only $(free_gb) GB free on the $PHYSICS24 volume; the floor is 150 GB" >&2
  echo "and this ladder writes ~4 GB of Les Houches files plus transient showers" >&2
  exit 1
fi

points () {
  for cur in nu:14 nubar:-14; do
    for tgt in p:331100 n:339900; do
      for E in $ES; do
        echo "${cur%%:*} ${cur#*:} ${tgt%%:*} ${tgt#*:} $E"
      done
    done
  done
}

do_gen () {
  # NOT `points | while`: a piped loop runs in a subshell, whose background
  # jobs the `wait` below would never see.  Process substitution keeps the
  # loop in this shell (bash 3.2 has it).
  while read -r c ih1 t set2 E; do
    R=$BASEDIR/${c}_${t}_E$E
    if [ -s "$R/pwgevents.lhe" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $R/pwgevents.lhe"; continue; fi
    mkdir -p "$R"
    python3 - "$CARD" "$R/powheg.input" "$E" "$ih1" "$set2" "$NEV" <<'PY'
import re, sys, os
src, dst, e, ih1, set2, nev = sys.argv[1:7]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "analysis"))
import beams
eb = beams.Beams("nu", float(e)).sqrt_s / 2.0
s = open(src).read()
s = re.sub(r"^numevts .*", f"numevts {nev}", s, flags=re.M)
s = re.sub(r"^ih1 .*", f"ih1   {ih1}", s, flags=re.M)
# BOTH ids: POWHEG-V2 demands its own QCDLambda5 whenever lhans1 differs
# from lhans2 (init_phys.f), and a neutrino beam reads no PDF anyway
s = re.sub(r"^lhans1 .*", f"lhans1 {set2}", s, flags=re.M)
s = re.sub(r"^lhans2 .*", f"lhans2 {set2}", s, flags=re.M)
s = re.sub(r"^ebeam1 .*", f"ebeam1 {eb:.14f}d0", s, flags=re.M)
s = re.sub(r"^ebeam2 .*", f"ebeam2 {eb:.14f}d0", s, flags=re.M)
# a different, REPRODUCIBLE seed per point (zlib, not hash(): Python's
# string hash is salted per process), so two points never share a stream
import zlib
seed = zlib.crc32(f"{float(e):g}:{ih1}:{set2}".encode()) % 900000 + 100
s = re.sub(r"^iseed .*", f"iseed {seed}", s, flags=re.M)
open(dst, "w").write(s)
PY
    limit
    ( cd "$R" && bench_run "$PWHG" > run_v2.log 2>&1 \
        && "$BENCH_REPO/powheg/strip_lhe_nan.py" pwgevents.lhe > strip_nan.log 2>&1 \
        && echo "gen done: ${c}_${t}_E$E  $(grep -m1 'total (btilde+remnants) cross section' pwg-stat.dat)" \
        || echo "gen FAILED: $R" ) &
  done < <(points)
  wait
}

do_shower () {
  # NOT `points | while`: a piped loop runs in a subshell, whose background
  # jobs the `wait` below would never see.  Process substitution keeps the
  # loop in this shell (bash 3.2 has it).
  while read -r c ih1 t set2 E; do
    R=$BASEDIR/${c}_${t}_E$E
    [ -s "$R/pwgevents.lhe" ] || { echo "no LHE yet: $R"; continue; }
    if [ -s "$R/dimuon_events.npz" ] && [ "${FORCE:-0}" != "1" ]; then echo "extracted already: $R"; continue; fi
    if [ -s "$R/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then echo "showered already: $R"; continue; fi
    card=$CMND
    if [ "$t" = n ]; then
      card=$R/powheg_shower_n.cmnd
      sed 's/^PDF:pSet = LHAPDF6:NNPDF40_nnlo_as_01180$/PDF:pSet = LHAPDF6:NNPDF40_nnlo_as_01180_n/' "$CMND" > "$card"
      grep -q "as_01180_n" "$card" || { echo "neutron PDF substitution failed for $card" >&2; exit 1; }
    fi
    limit
    ( cd "$BENCH_REPO/powheg" && bench_run ./main_powheg "$card" "$R/pwgevents.lhe" "$R/events" > "$R/shower.log" 2>&1 \
        && echo "shower done: ${c}_${t}_E$E  $(grep -c '^E ' "$R/events.hepmc") events" \
        || echo "shower FAILED: $R" ) &
  done < <(points)
  wait
}

do_extract () {
  # NOT `points | while`: a piped loop runs in a subshell, whose background
  # jobs the `wait` below would never see.  Process substitution keeps the
  # loop in this shell (bash 3.2 has it).
  while read -r c ih1 t set2 E; do
    R=$BASEDIR/${c}_${t}_E$E
    if [ -s "$R/dimuon_events.npz" ] && [ "${FORCE:-0}" != "1" ]; then echo "extracted already: $R"; continue; fi
    [ -s "$R/events.hepmc" ] || { echo "no shower to extract: $R"; continue; }
    limit
    ( python3 "$BENCH_REPO/analysis/faser_dimuon_cutflow.py" extract \
          --beam="$ih1" --target="$t" --energy="$E" --generator=powheg_v2 \
          "$R/events.hepmc" "$R/dimuon_events.npz" > "$R/extract.log" 2>&1 \
        && [ -s "$R/dimuon_events.npz" ] \
        && rm -f "$R/events.hepmc" \
        && echo "extracted: ${c}_${t}_E$E  $(tail -1 "$R/extract.log")" \
        || echo "extract FAILED: $R (HepMC kept)" ) &
  done < <(points)
  wait
}

case "$STAGE" in
  gen)     do_gen ;;
  shower)  do_shower ;;
  extract) do_extract ;;
  all)     do_gen; do_shower; do_extract ;;
  *) echo "stage must be gen, shower, extract or all" >&2; exit 2 ;;
esac
echo "stage $STAGE done in $BASEDIR"
