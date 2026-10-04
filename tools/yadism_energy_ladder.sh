#!/usr/bin/env bash
# The ANALYTIC sigma_fid(E) ladder that the FASER flux convolution rides on.
#
# WHY THIS EXISTS.  The FASER flux is a broad spectrum: two thirds of the
# neutrino flux and 70% of the muon flux sit BELOW 632 GeV, and their
# flux-averaged energy is 171 GeV -- not the 400 GeV point the benchmark
# computed.  Handing a whole bin the cross-section of its representative
# energy therefore overestimates the rate by 68% on the neutrino side, where
# sigma_fid is very nearly linear in E.  The cure is to convolute against a
# CONTINUOUS sigma_fid(E), and this script measures it: YADISM NLO (ZM-VFNS)
# on the benchmark's own fiducial region, at a ladder of energies spanning the
# flux, for BOTH currents.
#
# Each generator then rides this shape, normalised to its own computed points
# at 400 / 1000 / 4000 GeV -- see analysis/faser_rates.py.  A point costs
# about a minute, so the whole ladder is a few minutes and is cheap to redo.
#
# Usage:  tools/yadism_energy_ladder.sh [energies...]     (default: the ladder)
HERE=$(cd "$(dirname "$0")/.." && pwd)
. "$HERE/config.sh"

# Spans the vendored fluxes: the neutrino files start at 11.3 GeV and the muon
# files at 8.1 GeV; both die out above ~5 TeV.  400 / 1000 / 4000 are already
# computed by the benchmark proper and are NOT repeated here.
LADDER="10 15 25 40 60 100 200 300 632 1500 2000 3000 6000"
[ $# -gt 0 ] && LADDER="$*"

for e in $LADDER; do
  for s in analysis/yadism_calc.py analysis/yadism_cc_calc.py; do
    echo "=== $s at $e GeV"
    BENCH_ENERGY=$e python3 "$HERE/$s" nlo 2>&1 | grep -E "sigma_fid|CLOSURE|wrote"
  done
done
