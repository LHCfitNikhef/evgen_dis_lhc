#!/usr/bin/env bash
# Copy the paper plots into paper/figures/, where the LaTeX picks them up.
# The paper includes the vector PDF of each figure (user decision, 2026-10-09);
# the PNG twin stays in results_nu/ for a quick look only.
#
# WHY A COPY AND NOT A REFERENCE.  results_nu/ is gitignored -- everything a
# generator or an analysis writes is regenerable and stays out of git -- but a
# paper figure is a deliverable: a clone must build the document without first
# re-running the benchmark.  So the paper's copies are TRACKED, and this script
# is what keeps them in step with the report (user rule, 2026-08-29: the two
# are always in sync).
#
# Run it after regenerating the paper plots, before committing.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
mkdir -p "$REPO/paper/figures"
# >>> THE PAPER PLOTS V2, BY NAME. <<<  Listing them by name rather than by
# glob keeps a new report plot from slipping into the paper's tracked figures
# unnoticed; tools/check_paper_figures_used.py checks the list against each
# script's IN_PAPER flag.  The earlier paper_*.png glob that used to follow was
# removed with the 0.2 < y < 0.9 results (user, 2026-09-19): only survives.
EXTRA="results_nu/pp_sigma_vs_energy.pdf results_nu/pp_dis_distributions.pdf results_nu/pp_genie_vs_energy.pdf results_nu/pp_genie_dis_distributions.pdf results_nu/pp_genie_nondis.pdf results_nu/pp_genie_nondis_100GeV.pdf results_nu/pp_npdf_impact.pdf results_nu/pp_npdf_hadron.pdf results_nu/pp_genie_charm_vs_energy.pdf results_nu/pp_genie_charm_distributions.pdf results_nu/pp_charm_vs_energy.pdf results_nu/pp_charm_distributions.pdf results_nu/pp_dmeson_ratios.pdf results_nu/pp_hadron_level_faser_e.pdf results_nu/pp_shower_faser_e.pdf results_nu/pp_qed_dis.pdf results_nu/pp_qed_hadron.pdf results_nu/pp_faser_emulsion.pdf results_nu/pp_faser_electronic.pdf results_nu/pp_faser_dimuon_pdf.pdf results_nu/pp_sidis_yields.pdf results_nu/pp_nu_generators.pdf"

n=0
for e in $EXTRA; do
    [ -f "$REPO/$e" ] || { echo "  MISSING $e" >&2; continue; }
    cp "$REPO/$e" "$REPO/paper/figures/"
    echo "  $(basename "$e")  $(du -h "$REPO/$e" | cut -f1)"
    n=$((n + 1))
done
echo "$n paper figure(s) synced to paper/figures/"
