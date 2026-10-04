#!/usr/bin/env bash
# Analyse the POWHEG-LHE-showered-by-Herwig arms under every selection.
#
# Usage: tools/powheg_herwig_passes.sh [energies...]     (default: 1000)
#
# The hadronisation study's second half (user, 2026-09-08): one POWHEG matrix
# element showered by Pythia 8 and by Herwig 7, so that the shower and the
# hadronisation model are the only things that differ.  The samples are made
# by powheg/run_powheg_herwig.sh; this turns them into results.
#
# BOTH CURRENTS, EVERY SELECTION, AND THE CHARM TAG (CONVENTIONS.md rule 2b and the
# FASER tiers the rest of the benchmark is quoted under).  The Pythia arms
# `powheg` and `powheg_nu` already have all of these, so the comparison is
# available wherever the benchmark makes one.
#
# READ THE CHARM ROWS WITH THE DECLARED DEFICIT IN HAND.  Herwig discards
# 3.6% (muon) and 5.2% (neutrino) of these events in the shower, and the
# discards are flavour dependent -- charm hardest of all.  The rate is
# restored by normalising to POWHEG's integrator; the charm SHAPE is not.
# analyze.DECLARED_DEFICITS carries the measurement.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
ES=${*:-"1000"}
PY=${BENCH_PYTHON:-python3}

for E in $ES; do
  for S in inclusive faser_s faser_e faser_dimuon; do
    echo "===== POWHEG + Herwig7, $E GeV, $S ====="
    for K in powheg_hw powheg_hw_charmfinal; do
      BENCH_ENERGY=$E BENCH_SELECTION=$S $PY analysis/analyze.py "$K" \
        2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
    done
    for K in powheg_nu_hw powheg_nu_hw_charmfinal; do
      BENCH_ENERGY=$E BENCH_SELECTION=$S $PY analysis/analyze_nu.py "$K" \
        2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
    done
  done
done
