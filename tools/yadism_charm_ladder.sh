#!/usr/bin/env bash
# The FONLL charm-production ladder: YADISM NLO charm in FONLL at the same
# energies as the inclusive ladder (tools/yadism_energy_ladder.sh), on both
# currents, for the "Charm production in NC DIS" application of the
# Predictions for FASER tab (user, 2026-09-07) and its CC partner.
#
# WHY A LADDER OF ITS OWN.  faser_rates.py convolves the flux with
# sigma_yadism(E) x R_g(E), R_g held constant outside the three computed
# energies.  For charm the ratio to the INCLUSIVE ladder is the charm
# fraction, which rises from a few per cent at 100 GeV to 14% at 4 TeV and
# keeps falling below 400 GeV towards the threshold; holding it at its
# 400 GeV value under the low-energy half of the muon flux would overstate
# the charm rate there.  A charm ladder carries that shape itself, and R_g
# is then the O(1) generator-to-FONLL charm ratio the tab is about.
#
# Usage: tools/yadism_charm_ladder.sh [mu|nu|both] [E1 E2 ...]   (MAXJOBS=4)
# Output: results/histos_yadism_charm_nlo_fonll_damp_<E>.json (mu),
#         results_nu/... (nu); the three benchmark energies already exist.
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
. "$HERE/config.sh" > /dev/null
WHICH=${1:-both}
shift || true
LADDER=${*:-"10 15 25 40 60 100 200 300 632 1500 2000 3000 6000"}
MAXJOBS=${MAXJOBS:-4}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
mkdir -p "$HERE/logs"
run_one () {   # <current> <energy>
  local cur=$1 e=$2 script resdir tag out
  case "$cur" in mu) script=analysis/yadism_charm_calc.py; resdir=results;;
                 nu) script=analysis/yadism_cc_charm_calc.py; resdir=results_nu;; esac
  tag=$(python3 "$HERE/analysis/beams.py" tag "$cur" "$e")
  out=$HERE/$resdir/histos_yadism_charm_nlo_fonll_damp_$tag.json
  if [ -s "$out" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $out"; return 0; fi
  limit
  # the two calculators take different arguments: the NC one names the
  # flavour first, the CC one is charm-only
  local args="charm nlo fonll_damp"; [ "$cur" = nu ] && args="nlo fonll_damp"
  ( BENCH_ENERGY=$e python3 "$HERE/$script" $args \
      > "$HERE/logs/charm_ladder_${cur}_${e}.log" 2>&1 \
    && [ -s "$out" ] && echo "done: $cur $e GeV  $(grep -m1 -E 'sigma_fid' "$HERE/logs/charm_ladder_${cur}_${e}.log" | cut -c1-90)" \
    || echo "FAILED: $cur $e GeV (logs/charm_ladder_${cur}_${e}.log)" ) &
}
for e in $LADDER; do
  case "$WHICH" in
    mu)   run_one mu "$e" ;;
    nu)   run_one nu "$e" ;;
    both) run_one mu "$e"; run_one nu "$e" ;;
    *) echo "first argument must be mu, nu or both" >&2; exit 2 ;;
  esac
done
wait
echo "charm ladder ($WHICH) finished"
