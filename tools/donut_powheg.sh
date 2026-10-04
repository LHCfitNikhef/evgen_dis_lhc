#!/usr/bin/env bash
# POWHEG-V2 + Pythia 8 at DONUT energies, for the multiplicity comparison
# against arXiv:0711.0728 FIG. 9 (user, 2026-09-10).
#
# A LADDER, NOT A FLUX RUN.  POWHEG has no flux driver, so the beam energy is
# scanned and the multiplicity is folded afterwards with DONUT's own FIG. 2
# spectrum of INTERACTING neutrinos (analysis/donut_nch.py).  Each point
# generates the same number of events, so the ladder needs no cross-section:
# FIG. 2 already carries sigma(E), and multiplying by POWHEG's own sigma as
# well would count it twice.
#
# >>> THE TARGET IS THE FREE-NUCLEON AVERAGE OF IRON, IN ONE PDF SET. <<<
# DONUT's Emulsion Cloud Chambers are 1 mm stainless-steel sheets interleaved
# with emulsion, so about 86% of the mass is iron.  POWHEG scatters off a
# nucleon, so the target is (26 p + 30 n) / 56, built as an LHAPDF set by
# tools/make_isoscalar_pdf.py --Z 26 --N 30.  A set rather than two runs
# combined afterwards, because DIS is linear in the parton density and the
# average nucleon's cross-section IS the average of the two -- exactly the
# argument the nuclear-PDF study rests on, and it halves the ladder.
#
# WHAT POWHEG CANNOT SAY HERE, and it has to be said in place rather than left
# to look like a discrepancy:
#   * DIS ONLY.  The card generates above q2cut = 2.25 GeV2.  DONUT's sample
#     is every channel the beam produces, and at these energies quasi-elastic
#     and resonance production are a large part of it -- they are what makes
#     the data's n_ch = 1 bin.
#   * CHARGED CURRENT ONLY.  POWHEG-V2's nu-DIS process is W exchange; the
#     neutral current is about a quarter of the sample and has one track fewer.
#   * A FREE NUCLEON, with no final-state interaction.  Iron re-scatters the
#     hadrons on the way out and adds knocked-out protons; GENIE models that
#     and POWHEG does not.
# All three push the POWHEG multiplicity in known directions, which is the
# point of showing it beside GENIE rather than instead of it.
#
# Usage: tools/donut_powheg.sh [E1 E2 ...]
# Output: $POWHEG_V2/ladder-donut/E<GeV>/{pwgevents.lhe,events.hepmc}
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
CARD=$BENCH_REPO/powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-donut
CMND=powheg_nu1TeV.cmnd          # sets no beams: main_powheg reads the LHE
# Eight points from below the spectrum's rise to above its tail, log-spaced.
# The fold interpolates the SHAPE in log E between them, so the ladder has to
# be dense where the shape moves fastest, which is at the bottom.
ES=${*:-"12 20 32 50 75 110 160 230"}
MAXJOBS=${MAXJOBS:-4}
NEV=${NEV:-20000}
SET=341400                       # NNPDF40_nnlo_as_01180_Fe56free

limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }

for E in $ES; do
  R=$BASEDIR/E$E
  if [ -s "$R/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "exists: $R"; continue
  fi
  mkdir -p "$R"
  BENCH_REPO=$BENCH_REPO python3 - "$CARD" "$R/powheg.input" "$E" "$SET" "$NEV" <<'PY'
import os, re, sys
src, dst, e, pdfset, nev = sys.argv[1:6]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "analysis"))
import beams
eb = beams.Beams("nu", float(e)).sqrt_s / 2.0
s = open(src).read()
s = re.sub(r"^numevts .*", f"numevts {nev}", s, flags=re.M)
# BOTH ids: POWHEG-V2 demands its own QCDLambda5 whenever lhans1 differs from
# lhans2 (init_phys.f), and a neutrino beam reads no PDF anyway.  The trap is
# recorded in memory note faser-neutrino-reproduction: it exits SILENTLY.
s = re.sub(r"^lhans1 .*", f"lhans1 {pdfset}", s, flags=re.M)
s = re.sub(r"^lhans2 .*", f"lhans2 {pdfset}", s, flags=re.M)
s = re.sub(r"^ebeam1 .*", f"ebeam1 {eb:.14f}d0", s, flags=re.M)
s = re.sub(r"^ebeam2 .*", f"ebeam2 {eb:.14f}d0", s, flags=re.M)
s = re.sub(r"^iseed .*", f"iseed {1700 + int(float(e))}", s, flags=re.M)
open(dst, "w").write(s)
PY
  limit
  ( cd "$R" && bench_run "$PWHG" > run_v2.log 2>&1 \
      && ( cd "$BENCH_REPO/powheg" && bench_run ./main_powheg "$CMND" \
             "$R/pwgevents.lhe" "$R/events" > "$R/shower.log" 2>&1 ) \
      && echo "done: E$E  $(grep -c '^E ' "$R/events.hepmc" 2>/dev/null) showered" \
      || echo "FAILED: $R" ) &
done
wait
echo "ladder done in $BASEDIR"
