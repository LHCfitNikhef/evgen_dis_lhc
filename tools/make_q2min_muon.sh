#!/usr/bin/env bash
# RE-DERIVE THE WHOLE MUON SIDE AT A RAISED Q2 FLOOR.
#
# Usage:  tools/make_q2min_muon.sh [Q2MIN] [energies...]
#         tools/make_q2min_muon.sh 11 1000          # the anchor only
#         tools/make_q2min_muon.sh 11               # all three energies
#         tools/make_q2min_muon.sh 4                # back to the default
#
# WHY (user, 2026-09-01): "for muon DIS we should use Q2 > 11 GeV2 in all
# paper plots.  This breaks the symmetry with neutrino DIS but it is fine (it
# reflects different physics) and this way we ensure that scale variations are
# well defined and also that NNLO corrections make sense.  But let's always
# keep the option of going back to muon DIS with Q2min = 4 GeV2."
#
# That option is what this script is: ONE command re-derives every muon input
# a paper figure reads, at whatever floor it is given, and `selection.Q2_FLOOR`
# writes the floor into every filename -- so the Q2 > 4 results are never
# touched and going back is running it again with 4 (or simply not passing the
# knob, since the default results are already on disk).
#
# THE NEUTRINO SIDE IS DELIBERATELY NOT TOUCHED.  CONVENTIONS.md rule 2b asks that
# whatever is done for one current be done for the other, and here the
# asymmetry IS the physics: the neutral current's scale band and its NNLO
# correction are dominated by the bottom of the Q2 range (the 1/Q^4 photon
# propagator pins <Q2> near 15 GeV2), while the W propagator does not hold the
# charged current there and its band is +-1% at Q2 > 4 already.  Each figure
# states its own floor in its panel title rather than leaving the reader to
# assume they match.
#
# ORDER MATTERS: the FONLL inclusive reference is ASSEMBLED from the ZM
# inclusive and the two charm references, so those come first; the Herwig
# halves come before their combination; the bands come last, since they read
# the analysis's bins.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
cd "$REPO"
PY=${BENCH_PYTHON:-python3}

Q=${1:-11}
shift || true
ES=${*:-"1000 400 4000"}
SUF=""
[ "$Q" = "4" ] || SUF="_q2min${Q}"
export BENCH_Q2MIN=$Q

echo "== muon side at Q2 > $Q GeV2, energies: $ES  (suffix '${SUF:-none}') =="

for E in $ES; do
  export BENCH_ENERGY=$E
  echo "-- $E GeV: analytic references"
  $PY analysis/yadism_calc.py nlo
  $PY analysis/yadism_calc.py nnlo
  for O in nlo nnlo; do
    $PY analysis/yadism_charm_calc.py charm $O zm
    $PY analysis/yadism_charm_calc.py charm $O fonll_damp
  done
  # needs the three above, and writes the muon FONLL inclusive from them
  $PY analysis/yadism_fonll_inclusive.py nlo  || true
  $PY analysis/yadism_fonll_inclusive.py nnlo || true

  echo "-- $E GeV: generators"
  for K in powheg powheg_charmfinal sherpa sherpa_charmfinal; do
    $PY analysis/analyze.py "$K"
  done
  $PY herwig7/make_histos_nlo.py powheg
  $PY herwig7/make_histos_nlo.py powheg charmfinal
  $PY herwig7/make_histos_nlo.py powhegneg
  $PY herwig7/make_histos_nlo.py powhegneg charmfinal
  # the energy tag sits before the Q2 suffix, as analyze.result_path builds it
  ETAG=$($PY analysis/beams.py tag mu "$E")
  T=""; [ "$E" = "1000" ] || T="_$ETAG"
  R=$REPO/results
  $PY herwig7/combine_nlo.py \
      "$R/histos_herwig_nlo_powheg$T$SUF.json" \
      "$R/histos_herwig_nlo_powheg_neg$T$SUF.json" \
      "$R/histos_herwig_nlo_powheg_full$T$SUF.json"
  $PY herwig7/combine_nlo.py \
      "$R/histos_herwig_nlo_powheg_charmfinal$T$SUF.json" \
      "$R/histos_herwig_nlo_powheg_neg_charmfinal$T$SUF.json" \
      "$R/histos_herwig_nlo_powheg_full_charmfinal$T$SUF.json" \
      "Herwig 7.3.0 NLO (POWHEG, charm)"
done

# THE DIFFERENTIAL BANDS, at the anchor only -- paper plot 2 and 4 are 1 TeV
# figures and nothing else reads them.  The INTEGRATED bands that plots 1 and
# 2 use come from analysis/mhou_sigma_fonll.py, which has carried its own
# BENCH_Q2MIN suffix since 2026-08-30 and already covers all three energies.
unset BENCH_ENERGY
echo "-- 1 TeV: differential scale bands"
$PY analysis/mhou_diff.py --current mu
$PY analysis/mhou_diff.py --current mu --pto 2 --scheme fonll
$PY analysis/mhou_diff.py --current mu --flavour charm
$PY analysis/mhou_diff.py --current mu --pto 2 --scheme fonll --flavour charm

echo
echo "== done.  Now: tools/make_all_plots.sh && analysis/make_report.py =="
