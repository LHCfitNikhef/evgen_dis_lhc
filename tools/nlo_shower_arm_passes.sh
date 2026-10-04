#!/usr/bin/env bash
# The NLO-MATCHED shower comparisons, under every selection, on both currents.
#
# Usage: tools/nlo_shower_arm_passes.sh
#
# WHY THIS IS NOT `pythia8/run_pythia.sh --shower` (user, 2026-09-08: "it is
# better to test the impact of PS variations for NLO matched MCs, rather than
# LO ones").  Switching PartonShowers:model inside Pythia does not work on a
# POWHEG sample: the matching veto is implemented as SimpleShower UserHooks,
# which Vincia and Dire do not expose, and the two deliver 0 and 225 events of
# 59775 while exiting 0 (pythia8/shower/README.md, logs/shower_probe_powheg.log).
# So the LO samples carry the model switch, and the NLO-matched shower
# comparison is made in the two ways that DO work:
#
#   1. ONE POWHEG LES HOUCHES FILE, TWO GENERATORS.  powheg vs powheg_hw on
#      the muon side, powheg_nu vs powheg_nu_hw on the neutrino side -- the
#      same events showered by Pythia 8 and by Herwig 7.  This varies the
#      shower AND the hadronisation model together, which is what makes it
#      the hadronisation study as well; tools/powheg_herwig_passes.sh runs it.
#
#   2. SHERPA MC@NLO WITH ITS TWO OWN SHOWERS, CSS and Dire.  Native
#      matching on both sides, so nothing is worked around: the ME, the
#      matching and the hadronisation (AHADIC) are identical and only the
#      shower changes.  The samples have existed since the shower study but
#      were only ever analysed INCLUSIVELY -- this script adds the FASER
#      tiers, which is a re-analysis and not a regeneration.
#
# This script is (2).  Both currents, CONVENTIONS.md rule 2b.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
PY=${BENCH_PYTHON:-python3}

for S in faser_s faser_e faser_dimuon; do
  echo "===== Sherpa MC@NLO Dire vs CSS, 1 TeV, $S ====="
  for K in sherpa_dire sherpa_dire_charmfinal; do
    BENCH_SELECTION=$S $PY analysis/analyze.py "$K" \
      2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
  done
  for K in sherpa_nlo_dire sherpa_nlo_dire_charmfinal; do
    BENCH_SELECTION=$S $PY analysis/analyze_nu.py "$K" \
      2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
  done
done
