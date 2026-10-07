# The paper figures

One script per figure.  Each reads only tracked result files (`results/` for
muon DIS, `results_nu/` for neutrino DIS) and writes its PNG next to them, so
every figure rebuilds from a clone without any event sample:

    for f in analysis/paper_plots/pp*.py; do "$BENCH_PYTHON" "$f"; done
    tools/sync_paper_figures.sh        # copy the PNGs into paper/figures/
    tools/check_paper_plots.py         # each figure's quoted numbers vs its inputs
    tools/check_paper_claims.py        # the numbers quoted in the paper text

Every script exposes the same small contract (`tools/check_paper_plots.py`
enforces it): `TITLE`, `OUTPUT`, `RESULTS`, `IN_PAPER`, a `MESSAGE` (the
figure's physics statement, as shown in the report) and `CLAIMS`, a list of
checks that re-derive every number the message quotes from the result files.

| Fig. | Script | Content |
|---|---|---|
| 4.1 | `pp01_sigma_vs_energy.py` | fiducial NLO cross-section vs beam energy, both currents |
| 4.2 | `pp02_dis_distributions.py` | NLO differential DIS distributions at 1 TeV |
| 4.3 | `pp03_genie_vs_energy.py` | GENIE against the NLO matchings vs energy |
| 4.4 | `pp04_genie_dis_distributions.py` | GENIE differential distributions at 1 TeV |
| 5.1 | `pp06_charm_vs_energy.py` | fiducial charm cross-section vs energy |
| 5.2 | `pp05_genie_charm_vs_energy.py` | inclusive charm production, GENIE and references |
| 5.3 | `pp07b_dmeson_ratios.py` | D-meson ratios (charm fragmentation) vs energy |
| 5.4 | `pp07_charm_distributions.py` | differential charm distributions at 1 TeV |
| 5.5 | `pp05b_genie_charm_distributions.py` | GENIE differential charm distributions |
| 6.1 | `pp08_hadron_level_faser_e.py` | hadron-level observables under the FASERν selection |
| 6.2 | `pp09_shower_faser_e.py` | parton-shower dependence |
| 6.3 | `pp10_qed_dis.py` | QED corrections, DIS observables |
| 6.4 | `pp11_qed_hadron.py` | QED corrections, hadron-level observables |
| 7.1 | `pp12_faser_emulsion.py` | comparison with the FASERν emulsion measurement |
| 7.2 | `pp12b_faser_electronic.py` | comparison with the FASER electronic-detector measurement |
| 7.3 | `pp13_faser_dimuon_pdf.py` | dimuon yields under six PDF sets |
| 7.4 | `pp15_sidis_yields.py` | charged-hadron SIDIS yields in z (ν + ν̄) |
| B.1 | `ppA1b_genie_nondis_100GeV.py` | non-DIS processes in GENIE at 100 GeV |
| B.2 | `ppA1_genie_nondis.py` | non-DIS processes in GENIE |
| C.1 | `ppA2_npdf_impact.py` | impact of nuclear PDFs |
| D.1 | `ppA3_nu_generators.py` | GENIE, NuWro and GiBUU on a proton and on tungsten |

`ppA2b_npdf_hadron.py` (`IN_PAPER = False` pending the authors' decision) is
the hadron-level companion of Fig. C.1: nPDF bands on the four observables of
Fig. 6.1, neutrino DIS only (user decision 2026-10-07, an illustration), from
POWHEG-V2 reweighting (`analysis/npdf_hadron.py`).

`pp16_sidis_yields_kaon.py` (`IN_PAPER = False`) is the kaon companion of
Fig. 7.4, shown in the HTML report only.  Shared style lives in
`analysis/plotstyle.py` (LaTeX text, no silent fallback) and the legend labels
in `analysis/labels.py`.
