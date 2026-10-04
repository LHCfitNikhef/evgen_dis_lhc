#!/usr/bin/env bash
# POWHEG-V2 sigma_fid(E) ladder for the FASER neutrino-rate reproduction
# (Step B, user 2026-09-03: "sufficient to use POWHEG-V2 ... we still need
# Q2 > 4 GeV2"): nu_mu and nubar_mu on a proton and on a neutron, at a ladder
# of beam energies, each a short POWHEG-V2 run whose Les Houches events give
# the fiducial fraction exactly -- Q2 and y are leptonic and the shower does
# not touch them with QED off -- so no showering is needed.
#
# Usage: tools/powheg_v2_faser_ladder.sh [E1 E2 ...]      (MAXJOBS=4 parallel)
# Output: $POWHEG_V2/ladder-faser/<nu|nubar>_<p|n>_E<GeV>/{pwg-stat.dat,pwgevents.lhe}
#         read by analysis/faser_powheg_rates.py
#
# THE NEUTRON is the isospin-swapped PDF set data/pdfs/lhapdf/ (id 339900,
# tools/make_neutron_pdf.py), reached by PREPENDING LHAPDF_DATA_PATH.  The
# card is the 1 TeV production card with the beam, the lepton and the target
# set changed and 20000 events instead of 150000.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
CARD=$BENCH_REPO/powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-faser
ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
MAXJOBS=${MAXJOBS:-4}
NEV=${NEV:-20000}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
for cur in nu:14 nubar:-14; do
  for tgt in p:331100 n:339900; do
    for E in $ES; do
      c=${cur%%:*}; ih1=${cur#*:}; t=${tgt%%:*}; set2=${tgt#*:}
      R=$BASEDIR/${c}_${t}_E$E
      if [ -s "$R/pwgevents.lhe" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $R"; continue; fi
      mkdir -p "$R"
      python3 - "$CARD" "$R/powheg.input" "$E" "$ih1" "$set2" "$NEV" <<'PY'
import re, sys
src, dst, e, ih1, set2, nev = sys.argv[1:7]
import os
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
open(dst, "w").write(s)
PY
      limit
      ( cd "$R" && bench_run "$PWHG" > run_v2.log 2>&1 \
          && echo "done: ${c}_${t}_E$E  $(grep -m1 'total (btilde+remnants) cross section' pwg-stat.dat)" \
          || echo "FAILED: $R" ) &
    done
  done
done
wait
echo "ladder done in $BASEDIR"
