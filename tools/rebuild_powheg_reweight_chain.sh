#!/usr/bin/env bash
# Re-derive everything that depends on the POWHEG-V2 member reweighting.
#
# Usage: tools/rebuild_powheg_reweight_chain.sh [energies...]   (default: all)
#
# WHY THIS IS A SCRIPT.  Two bugs found on 2026-08-28 invalidated every
# reweighting-derived number at 1 TeV, and the repair touches six outputs in a
# fixed order.  Doing that by hand is how a stale file survives a repair --
# so the order is written down once, here, and re-running it is the way to
# reproduce the whole chain.
#
# THE TWO BUGS, both silent:
#   1. `powheg/rwgtnu_job_1` was showered TWO MINUTES before the lhe_index fix
#      landed, so its indices counted Pythia's retries and the LHE->HepMC join
#      was scrambled.  It still produced plausible bands: the dimuon band came
#      out at the INCLUSIVE size because a scrambled join makes any subset
#      behave like the whole sample.  Caught by asking whether dimuon-selected
#      events are charm -- 28% against the ~100% physics demands.
#   2. `powheg_variants.sh` resolved 1 TeV nu to `prod-nu1TeV`, the superseded
#      pre-CKM-fix tree, while the published sample comes from
#      `prod-nu1TeV-ckm`.  So the reweighting ran on a sample the published
#      results are not from.  Caught by reading nu_job_1/shower.log.
#
# 400 GeV and 4 TeV were produced after the CKM fix and were never affected.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
cd "$BENCH_REPO"

ES=${*:-"400 1000 4000"}
for E in $ES; do
  TAG=$(python3 analysis/beams.py tag nu "$E")
  W=$POWHEG_V2/rwgt-nu$TAG/weights.npz
  if [ ! -f "$W" ]; then
    echo "== nu $E GeV: no weights.npz yet, skipping =="
    continue
  fi
  echo "== nu $E GeV [$TAG] =="
  # 1. inclusive charm fraction with PDF bands and the scale envelope
  BENCH_ENERGY=$E python3 analysis/powheg_pdf_members.py
  # 2. the hadron-level selections: central sets + scale, and the member bands
  for SEL in inclusive faser_s faser_dimuon; do
    if [ -f "$BENCH_REPO/powheg/rwgtnu_job_${TAG}_1/events.hepmc" ]; then
      BENCH_ENERGY=$E BENCH_SELECTION=$SEL \
        python3 analysis/powheg_selection_members.py || \
        echo "   ($SEL: skipped -- see the message above)"
    fi
  done
done

echo "== figures =="
python3 analysis/plot_unc_methods.py
python3 analysis/plot_pdf_dependence_mc.py
python3 analysis/plot_pdf_crosscheck.py || echo "  (crosscheck skipped)"
echo "== done; rebuild the report next =="
