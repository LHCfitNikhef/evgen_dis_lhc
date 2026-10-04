#!/usr/bin/env bash
# The QED arms under the FASERnu selection, with the emulsion observables.
#
# Usage: tools/tier_e_qed_passes.sh
#
# Paper plot 11 (user, 2026-09-08) shows the impact of QED corrections on the
# HADRON-LEVEL distributions inside the FASERnu selection, as paper plot 10
# shows it on the DIS observables in the inclusive region.  The hadron-level
# half needs nch05 and dphix -- the multiplicity inside the emulsion's track
# acceptance and the azimuth against the summed hadron system -- and the QED
# tier results predate both.  A re-parse, not a regeneration.
#
# >>> THE BASELINE OF THE QED STUDY IS POWHEG-RES ON BOTH CURRENTS, and that
# is not the benchmark's neutrino row. <<<  Every QED arm re-showers the SAME
# POWHEG-RES Les Houches events with radiation switched on in stages, so the
# comparison is exact -- same matrix element, same matching, same events --
# and the denominator has to be the POWHEG-RES sample, not POWHEG-V2.  On the
# neutrino side that means `powheg_res_nu`, whose own caveat (diagonal CKM,
# inclusive only) does not touch a RATIO in which it cancels.
#
# `fsrisr` EXISTS ON THE MUON SIDE ONLY, and it is physics: initial-state
# radiation off the lepton line needs a CHARGED incoming lepton, and the
# neutrino beam has none.  Said in the figure rather than left as a gap.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
PY=${BENCH_PYTHON:-python3}

for K in powheg_qed_fsr powheg_qed_fsrisr powheg_qed_full; do
  echo "===== muon NC, FASERnu selection, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
for K in powheg_res_nu powheg_res_nu_qed_fsr powheg_res_nu_qed_full; do
  echo "===== neutrino CC, FASERnu selection, $K ====="
  BENCH_SELECTION=faser_e $PY analysis/analyze_nu.py "$K" \
    2>&1 | grep -E "sigma_fid|closure|wrote|FAILED|refus" | sed 's/^/  /'
done
