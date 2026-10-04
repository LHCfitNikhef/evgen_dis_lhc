#!/usr/bin/env bash
# Regenerate the report figures that are drawn from result JSONs outside the
# paper plots: the NOMAD comparison, the FASER dimuon cut flow, the
# single-inclusive pion spectra and the DONUT multiplicity.  Cheap --
# matplotlib over result JSONs, no event files -- so it is always safe to
# re-run in full.  The paper plots have their own scripts, one per figure, in
# analysis/paper_plots/.
#
# ONLY THE FINAL REGION SURVIVES (user, 2026-09-19).  The per-(selection,
# energy) comparison views, the cross-energy figures of the 0.2 < y < 0.9
# scan, the FASER inclusive and charm rates, the arXiv:2402.13318 and
# arXiv:2506.13889 reproductions and the nuclear-PDF study were removed from
# this script together with the results they drew.
#
# Usage: tools/make_all_plots.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
# the analysis interpreter (config.sh: $BENCH_PYTHON, not a bare python3)
. "$BENCH_REPO/config.sh"

# the NOMAD dimuon comparison (Comparison with data tab)
echo "===== NOMAD ====="
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_nomad.py" 2>&1 \
  | grep -E 'wrote|skipp|Traceback|Error' | sed 's/^/  /' || true

echo "===== FASER dimuons, pions, DONUT ====="
# the dimuon cut flow (2026-09-07): its rates come from the per-event tables
# of tools/genie_dimuon_ladder.sh and tools/powheg_v2_dimuon_ladder.sh,
# recomputed here so the figure and the report's table read one JSON
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_dimuon_cutflow.py" rates 2>&1 \
  | grep -E 'written|Traceback|Error|incomplete' | sed 's/^/  /'
# (the two figures compare with a FASER-internal estimate and are not part of
#  the public release; skipped where absent)
for s in plot_faser_dimuon_cutflow plot_faser_dimuon_rates; do
  [ -f "$BENCH_REPO/analysis/$s.py" ] || continue
  "$BENCH_PYTHON" "$BENCH_REPO/analysis/$s.py" 2>&1 \
    | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
done
# the PDF dependence of the POWHEG-V2 cut flow: from the member weights of
# tools/powheg_v2_dimuon_ladder_pdf.sh, skipped (with a message) until they exist
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_dimuon_cutflow.py" pdf 2>&1 \
  | grep -E 'written|Traceback|Error|fewer' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_faser_dimuon_pdf.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
# the single-inclusive pion spectra (region): rates first, figures after
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_pions.py" rates 2>&1 \
  | grep -E 'written|Traceback|Error' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_faser_pions.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
# and the NNLO SIDIS comparison and yields, which read the same JSON and were
# missing from this script until 2026-09-21 -- a figure nothing regenerates
# goes stale under a correct pipeline, which is the failure mode of rule 2
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_sidis.py" 2>&1 \
  | grep -E 'wrote|written|Traceback|Error|missing' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_faser_sidis.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
# and the Tier E efficiency factors (user, 2026-09-21), the bridge between a
# calculation in the region alone and a FASER yield
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_sidis_efficiency.py" 2>&1 \
  | grep -E 'written|Traceback|Error' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_faser_sidis_efficiency.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
# the region-only POWHEG-V2 benchmark for an analytic calculation, and the
# package that shares it with the fluxes (user, 2026-09-21)
"$BENCH_PYTHON" "$BENCH_REPO/analysis/faser_sidis_region_yields.py" 2>&1 \
  | grep -E 'written|Traceback|Error' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_faser_sidis_region_yields.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/tools/export_sidis_share.py" 2>&1 \
  | grep -E 'Traceback|Error' | sed 's/^/  /'
# the fragmentation functions a calculation of those yields could use
# (2026-09-21); the LHAPDF sets come from tools/install_ff_sets.sh
"$BENCH_PYTHON" "$BENCH_REPO/analysis/ff_comparison.py" 2>&1 \
  | grep -E 'written|Traceback|Error' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_ff_comparison.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'

# the DONUT multiplicity comparison (user, 2026-09-10): the per-generator
# extractions are done by their own drivers, so this only re-merges and
# re-draws.  Skipped with a message when no generator has run yet.
"$BENCH_PYTHON" "$BENCH_REPO/analysis/donut_nch.py" --merge 2>&1 \
  | grep -E 'wrote|Traceback|Error|no per-generator' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/plot_donut_nch.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'
"$BENCH_PYTHON" "$BENCH_REPO/analysis/donut_xsec.py" 2>&1 \
  | grep -E 'wrote|Traceback|Error|missing' | sed 's/^/  /'

echo "plots done."
