#!/usr/bin/env bash
# Re-analyse the 1 TeV Tier E samples so they carry the emulsion observables.
#
# Usage: tools/tier_e_hadron_passes.sh
#
# WHY (2026-09-08).  paper plot 8 draws what FASER's emulsion detector
# actually reconstructs, and two of its four observables did not exist as
# histograms until that day:
#
#   nch05  charged hadrons above 1 GeV INSIDE the emulsion's tan(theta) < 0.5
#          track acceptance.  `nch1` counts every charged hadron above 1 GeV,
#          including tracks that leave the acceptance.
#   dphix  the azimuth between the lepton and the VECTOR SUM of the charged
#          hadrons -- the back-to-back variable arXiv:2403.12520 cuts on, and
#          the one Tier E has always cut on.  `dphi` is the minimum over
#          INDIVIDUAL hadrons: a different quantity with the same everyday
#          name, and the neutrino parser did not have even that one.
#
# Both were added to analyze.py and analyze_nu.py; this re-parses the six
# samples the figure needs.  A RE-PARSE, NOT A REGENERATION: the samples are
# untouched and every number already published from them is reproduced, with
# two more histograms beside it.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
PY=${BENCH_PYTHON:-python3}

for K in powheg sherpa genie; do
  echo "===== muon NC, Tier E, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
for K in powheg_nu sherpa_nlo genie_lo; do
  echo "===== neutrino CC, Tier E, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze_nu.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
