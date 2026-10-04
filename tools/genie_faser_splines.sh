#!/usr/bin/env bash
# GENIE default-tune TOTAL charged-current cross-sections for the FASER rate
# reproduction (arXiv:2402.13318 Table I uses "the default GENIE cross
# section"): nu_e, nubar_e, nu_mu, nubar_mu on a free proton and a free
# neutron, every CC process of the default list, up to 8 TeV.  The TAU
# flavours are here too (2026-09-04): arXiv:2402.13318 Table I quotes a
# nu_tau row and the POWHEG charm flux of arXiv:2309.12793 supplies the
# neutrinos for it, so leaving them out would be the one flavour of that
# table we could not answer.  Their cross-section is NOT the muon one -- the
# tau mass suppresses it strongly at FASER energies -- so it has to be built
# rather than borrowed, which is exactly why it is built here.
#
# Usage: tools/genie_faser_splines.sh            (8 gmkspl runs in parallel)
# Output: genie/splines/faser/<pid>_<p|n>_cc_e8000.xml, read by
#         analysis/faser_genie_rates.py, which combines p and n as tungsten.
#
# WHY FREE NUCLEONS.  Tungsten is 74 p + 110 n; at the TeV energies that
# dominate the FASERnu rate the nuclear model is a per-cent matter, and a
# nuclear-target spline costs hours where a nucleon one costs minutes.
#
# THE VALIDITY CAP.  G18_02a declares itself valid to 1 TeV and gmkspl
# silently returns a 1 TeV spline beyond it (memory: genie-energy-validity);
# config/hienergy raises GVLD-Emax to 10 TeV.  The flux reaches 5.7 TeV, so
# this is an extrapolation of the tune above 1 TeV, exactly as the published
# prediction must also have been.  The p8 overlay is the benchmark's Pythia8
# decay/hadronisation substitution (this build has no Pythia6); it does not
# touch the cross-section model.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
G=$BENCH_REPO/genie
TUNE=G18_02a_00_000
EMAX=8000
OUT=$G/splines/faser
mkdir -p "$OUT"
for pid in 12 -12 14 -14 16 -16; do
  for tgt in p:1000010010 n:1000000010; do
    lab=${tgt%%:*}; code=${tgt#*:}
    f=$OUT/${pid}_${lab}_cc_e$EMAX.xml
    if [ -s "$f" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $f"; continue; fi
    # the environment is set INSIDE the wrapped shell: macOS strips
    # DYLD_LIBRARY_PATH across the caffeinate wrapper (as run_genie*.sh do)
    ( . "$BENCH_REPO/config.sh" > /dev/null
      cd "$OUT" && bench_run "$BENCH_SHELL" -c "
          source $G/genie_setup.sh > /dev/null
          export GXMLPATH=$G/config/hienergy:$G/config/p8
          gmkspl -p $pid -t $code -e $EMAX --tune $TUNE \
                 --event-generator-list CC -o $f" \
          > "$OUT/gmkspl_${pid}_${lab}.log" 2>&1 \
      && echo "done: $f  (last knot $(grep -o '<E> *[0-9.]*' "$f" | tail -1))" ) &
  done
done
wait
echo "all splines built in $OUT"
