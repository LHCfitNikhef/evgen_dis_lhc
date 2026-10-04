#!/usr/bin/env bash
# Sherpa MC@NLO for the comparison with FASER's emulsion data, at THEIR cuts.
#
# Usage: tools/sherpa_faser_ladder.sh [E1 E2 ...]      (default: the 5-point ladder)
#        NEV=20000 MAXJOBS=4 by default; ONE=nu:p:1000 runs a single point.
#
# WHY THIS IS A REGENERATION AND NOT A RE-ANALYSIS (user, 2026-09-08: "for the
# paper version of the plot, let's also add SHERPA NLO, I understand this
# requires regeneration but this is an important plot").  FASER imposes NO y
# window, and Sherpa cuts 0.2 < y < 0.9 AT GENERATION through its INEL
# selector, so its samples contain no event outside that window at all.  The
# card here is NuDIS_NLO_ckm3 with that selector removed and the beam energy
# moved; everything else -- the full-CKM patch, the massless-charm convention,
# the ctau > 10 mm stability rule -- is untouched.
#
# >>> THE y WINDOW WAS AN AHADIC WORKAROUND, SO REMOVING IT IS THE RISK. <<<
# It was never physics (see the memory note on the y cut's provenance): it was
# introduced because Sherpa's cluster hadronisation struggles at the edges of
# the range.  Whether it still does is what the first point tests, and the
# driver reports the hadronisation-failure rate per point rather than letting
# a quietly depleted sample through.
#
# SIX ENERGIES, NOT TEN.  The selection needs a 200 GeV lepton, so nothing
# below that can contribute at all -- measured, every POWHEG ladder point up
# to and including 200 GeV passes 0.00% -- and the flux-weighted selected
# spectrum is spent by 4 TeV.  Each point costs about forty minutes of
# integration, so the ladder is put where the rate is.
#
# >>> BUT 200 GeV IS IN IT, AND IT IS NOT WASTED. <<<  analysis/faser_emulsion_shapes.py
# ZEROES the cross-section below a ladder's FIRST point and holds the selected
# fraction flat at its ends, so a ladder starting at 400 GeV contributes
# nothing at all below 400 -- while the POWHEG and GENIE ladders, which have a
# 200 GeV point at exactly 0.00%, interpolate up from it and put 6.5% of the
# selected rate between 200 and 400 GeV.  Without the 200 GeV anchor the
# Sherpa curve would be the same physics folded over a 6.5%-narrower spectrum,
# in a figure whose entire content is a comparison of shapes.  The point costs
# an integration and delivers no selected events, which is the correct answer
# and the reason it is needed.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
# THE NEUTRON IS A LOCAL PDF SET.  tools/make_neutron_pdf.py writes the
# isospin-swapped NNPDF40_nnlo_as_01180_n into data/pdfs/lhapdf/ (id 339900),
# which is not on LHAPDF's own search path -- Sherpa dies with "PDF ... does
# not exist in any of the loaded libraries", which is how every neutron point
# of the first pass failed.  Same line as the POWHEG ladders.
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
SRC=$SHERPA_RUNS/NuDIS_NLO_ckm3
ES=${*:-"200 400 700 1000 2000 4000"}
NEV=${NEV:-20000}
MAXJOBS=${MAXJOBS:-4}
PY=${BENCH_PYTHON:-python3}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 10; done; }

[ -d "$SRC/Process" ] || { echo "no $SRC/Process -- the NLO libraries must be built first" >&2; exit 1; }

# The antineutrino seed: the one point that COMPILES the -14 93 -> -13 93
# libraries.  Every other antineutrino point links to its Process/.
NUBAR_SEED=$SHERPA_RUNS/NuDIS_NLO_FASER_nubar_p_E1000

# Amegic writes C++ for a process it has not seen and then EXITS, asking for
# ./makelibs.  A driver that does not expect that reads the exit as a failed
# integration; one that does not RERUN afterwards integrates nothing.  Only
# the seed ever takes this branch, and only once.
seed_libs () {                   # <dir> <nu|nubar>
  local d=$1 cur=$2
  [ "$cur" = nubar ] || return 0
  [ "$d" = "$NUBAR_SEED" ] || return 0
  [ -d "$d/Process/Amegic/lib" ] && return 0
  cp -f "$SRC/makelibs" "$d/makelibs"
  ( cd "$d" && bench_run "$SHERPA_INSTALL/bin/Sherpa" -e 0 > libs.log 2>&1 \
      || true )
  ( cd "$d" && ./makelibs -j 4 >> libs.log 2>&1 ) || {
      echo "makelibs failed in $d -- see libs.log" >&2; return 1; }
  echo "  built the antineutrino process libraries in $(basename "$d")"
}

point () {                       # <nu|nubar> <p|n> <E>
  local cur=$1 tgt=$2 E=$3
  local d=$SHERPA_RUNS/NuDIS_NLO_FASER_${cur}_${tgt}_E${E}
  if [ -s "$d/faserdata_events.npz" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "done already: $(basename "$d")"; return 0; fi
  mkdir -p "$d"
  # THE COMPILED LIBRARIES ARE SHARED WITHIN A BEAM AND NOT ACROSS ONE.
  # Only the beam energy and the target PDF change between points of the same
  # beam, so Process/ is linked rather than rebuilt -- but Amegic writes and
  # compiles process-specific code, and `-14 93 -> -13 93` is a DIFFERENT
  # process from `14 93 -> 13 93`.  The antineutrino points therefore share a
  # second set, built once by the seed point below.
  local libs=$SRC/Process
  [ "$cur" = nubar ] && libs=$NUBAR_SEED/Process
  if [ "$d" != "$NUBAR_SEED" ]; then
      [ -e "$d/Process" ] || ln -s "$libs" "$d/Process"
  fi
  local pid=14; [ "$cur" = nubar ] && pid=-14
  local lep=13; [ "$cur" = nubar ] && lep=-13
  # SHERPA WRITES `evtfull`, NOT `evtfull.hepmc`.  EVENT_OUTPUT names the
  # file exactly; the production runs have always written it that way and the
  # analysis has always read it that way.  Getting this wrong cost one full
  # 37-minute integration whose events were on disk the whole time.
  local pdf=NNPDF40_nnlo_as_01180
  [ "$tgt" = n ] && pdf=NNPDF40_nnlo_as_01180_n
  "$PY" - "$SRC/Sherpa.yaml" "$d/Sherpa.yaml" "$E" "$pid" "$lep" "$pdf" <<'PY'
import re, os, sys
src, dst, e, pid, lep, pdf = sys.argv[1:7]
sys.path.insert(0, os.path.join(os.environ["BENCH_REPO"], "analysis"))
import beams
b = beams.Beams("nu", float(e))
# the EXACT c.m. energies of a massless neutrino of energy E on a proton at
# rest, the convention the production card documents and uses
s = b.sqrt_s ** 2
mp = beams.M_P if hasattr(beams, "M_P") else 0.9382720813
enu = (s - mp * mp) / (2.0 * b.sqrt_s)
ep = (s + mp * mp) / (2.0 * b.sqrt_s)
t = open(src).read()
t = re.sub(r"^BEAMS: .*$", f"BEAMS: [{pid}, 2212]", t, flags=re.M)
t = re.sub(r"^BEAM_ENERGIES: .*$", f"BEAM_ENERGIES: [{enu:.6f}, {ep:.6f}]",
           t, flags=re.M)
t = re.sub(r"^PDF_SET: .*$", f"PDF_SET: [None, {pdf}]", t, flags=re.M)
t = re.sub(r"^MPI_PDF_SET: .*$", f"MPI_PDF_SET: {pdf}", t, flags=re.M)
# >>> THE y WINDOW GOES, AND NOTHING ELSE DOES. <<<  FASER imposes no y cut;
# the Q2 > 4 floor stays, because it is the perturbative floor this whole
# benchmark is defined above and the POWHEG ladder carries its own.
t = re.sub(r"^- \[INEL, .*\]\n", "", t, flags=re.M)
# ...AND THE COMMENT ABOVE IT, which otherwise survives to describe a cut that
# is no longer there.  A generated card whose comment says "0.2<y<0.9" while
# its selector list does not is exactly the kind of thing that gets believed
# six months later.
t = t.replace("# same fiducial region as the rest of the benchmark: "
              "Q2 > 4 GeV2, 0.2<y<0.9",
              "# FASER's own region: Q2 > 4 GeV2 and NOTHING else.  The y\n"
              "# window of the benchmark is removed here on purpose -- FASER\n"
              "# imposes none -- by tools/sherpa_faser_ladder.sh.")
t = re.sub(r"^- \[Q2, \d+, \d+,", f"- [Q2, {pid}, {lep},", t, flags=re.M)
# >>> AND THE PROCESS ITSELF, which the first pass forgot. <<<  Changing
# BEAMS to an antineutrino while the PROCESSES block still asks for
# `14 93 -> 13 93` leaves Sherpa with a beam it cannot use: it reports
# "Error in initialising ISR (vmub -> vmu)", then "No hard process found",
# and exits with status 0 -- a silent failure of exactly the kind CONVENTIONS.md
# rule 2 is about.  Every antineutrino point of the first pass died here.
t = re.sub(r"^- 14 93 -> 13 93:", f"- {pid} 93 -> {lep} 93:", t, flags=re.M)
open(dst, "w").write(t)
PY
  grep -q "INEL" "$d/Sherpa.yaml" && { echo "the y selector survived in $d" >&2; return 1; }
  # THE CARD MUST NAME THIS POINT'S BEAM, in the process as well as in BEAMS.
  grep -q "^- $pid 93 -> $lep 93:" "$d/Sherpa.yaml" || {
      echo "the process block in $d is not ${pid} 93 -> ${lep} 93" >&2; return 1; }
  (
    cd "$d" \
      && seed_libs "$d" "$cur" \
      && bench_run "$SHERPA_INSTALL/bin/Sherpa" -e 0 > integ.log 2>&1 \
      && bench_run "$SHERPA_INSTALL/bin/Sherpa" -e "$NEV" \
           'EVENT_OUTPUT: ["HepMC3_GenEvent[evtfull]"]' > gen.log 2>&1 \
      && "$PY" "$BENCH_REPO/analysis/faser_emulsion_shapes.py" extract \
           --beam="$pid" --target="$tgt" --energy="$E" --generator=sherpa_nlo \
           evtfull "$d/faserdata_events.npz" > extract.log 2>&1 \
      && rm -f evtfull \
      && echo "sherpa ${cur}_${tgt}_E${E}: $(tail -1 extract.log)" \
      || echo "sherpa FAILED: ${cur}_${tgt}_E${E}" >&2
  )
}

if [ -n "${ONE:-}" ]; then
  IFS=: read -r c t e <<< "$ONE"
  point "$c" "$t" "$e"
  exit 0
fi
# NEUTRINO POINTS FIRST AND IN PARALLEL; they share libraries that exist.
for tgt in p n; do
  for E in $ES; do limit; point nu "$tgt" "$E" & done
done
wait
# THEN THE ANTINEUTRINO SEED, ALONE.  It compiles the libraries the other
# nine link to, so nothing may start before it has finished.
point nubar p 1000
[ -d "$NUBAR_SEED/Process/Amegic/lib" ] || {
    echo "the antineutrino libraries were not built -- stopping rather than" >&2
    echo "starting nine points that would each try to build their own." >&2
    exit 1; }
for tgt in p n; do
  for E in $ES; do
    [ "$tgt" = p ] && [ "$E" = 1000 ] && continue
    limit; point nubar "$tgt" "$E" &
  done
done
wait
echo "sherpa FASER ladder done"
