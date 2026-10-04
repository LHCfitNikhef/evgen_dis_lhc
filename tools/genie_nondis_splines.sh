#!/usr/bin/env bash
# GENIE default-tune TOTAL electromagnetic cross-sections for the muon side
# of the "Impact of non-DIS processes" study: mu- on a free proton and a free
# neutron, every EM process of the default list (QEL-EM, RES-EM, DIS-EM,
# MEC-EM), up to 8 TeV.  The NEUTRINO side of the same study reads the CC
# splines tools/genie_faser_splines.sh already builds for the FASER rate
# reproduction, which are the full default CC list on both nucleons; this
# script is the muon counterpart that CONVENTIONS.md rule 2b requires, and it
# follows that one line for line.
#
# Usage: tools/genie_nondis_splines.sh            (2 gmkspl runs in parallel)
# Output: genie/splines/faser/13_<p|n>_em_e8000.xml, read by
#         analysis/genie_nondis.py, which splits them process by process.
#
# THE VALIDITY CAP.  G18_02a declares itself valid to 1 TeV and gmkspl
# silently returns a 1 TeV spline beyond it (memory: genie-energy-validity);
# config/hienergy raises GVLD-Emax to 10 TeV.  Above 1 TeV this is an
# extrapolation of the tune, as it is for every 4 TeV GENIE number here.
# ALWAYS check the last <E> knot of a new spline (it is echoed below).
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
G=$BENCH_REPO/genie
TUNE=G18_02a_00_000
EMAX=8000
OUT=$G/splines/faser
mkdir -p "$OUT"
pid=13
for tgt in p:1000010010 n:1000000010; do
  lab=${tgt%%:*}; code=${tgt#*:}
  f=$OUT/${pid}_${lab}_em_e$EMAX.xml
  if [ -s "$f" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $f"; continue; fi
  # the environment is set INSIDE the wrapped shell: macOS strips
  # DYLD_LIBRARY_PATH across the caffeinate wrapper (as run_genie*.sh do)
  ( . "$BENCH_REPO/config.sh" > /dev/null
    cd "$OUT" && bench_run "$BENCH_SHELL" -c "
        source $G/genie_setup.sh > /dev/null
        export GXMLPATH=$G/config/hienergy:$G/config/p8
        gmkspl -p $pid -t $code -e $EMAX --tune $TUNE \
               --event-generator-list EM -o $f" \
        > "$OUT/gmkspl_${pid}_${lab}_em.log" 2>&1 \
    && echo "done: $f  (last knot $(grep -o '<E> *[0-9.]*' "$f" | tail -1))" ) &
done
wait
echo "all EM splines built in $OUT"
