#!/usr/bin/env bash
# The SHOWER-VARIATION arms under Tier E, with the emulsion observables.
#
# Usage: tools/tier_e_shower_passes.sh
#
# paper plot 9 (user, 2026-09-08) compares four matched samples under FASER
# Tier E: POWHEG + Pythia 8, POWHEG + Herwig 7, Sherpa MC@NLO with its default
# CSS shower and Sherpa MC@NLO with Dire.  The first is done by
# tools/tier_e_hadron_passes.sh, which paper plot 8 also needs; these are the
# other three, on both currents.
#
# THEY NEED THE SAME TWO OBSERVABLES paper plot 8 introduced -- nch05, the
# multiplicity inside the emulsion's tan(theta) < 0.5 acceptance, and dphix,
# the azimuth against the summed hadron system -- and the Sherpa Dire results
# predate both while the Herwig ones predate nch05's integer binning.  A
# re-parse, not a regeneration.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
PY=${BENCH_PYTHON:-python3}

for K in powheg_hw sherpa_dire; do
  echo "===== muon NC, Tier E, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
for K in powheg_nu_hw sherpa_nlo_dire; do
  echo "===== neutrino CC, Tier E, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze_nu.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
