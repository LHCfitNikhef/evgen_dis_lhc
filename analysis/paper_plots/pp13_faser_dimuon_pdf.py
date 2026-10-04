#!/usr/bin/env python3
"""Paper plot 13: the POWHEG-V2 NLO dimuon signal at FASER under six
PDF sets.

Promoted from the report's "Predictions for FASER -> Dimuon production in
CC DIS" section (analysis/plot_faser_dimuon_pdf.py) at the user's request,
2026-09-21: "the POWHEG-V2 NLO dimuons at FASER with 6 different PDF sets
must be added to the list of paper plots, before the SIDIS ones".  It takes
the number 13, freed the same day when the old pp13 (pions against the NNLO
SIDIS calculation) and pp14 (pions in the emulsion, generator against
generator) were dropped from the paper plots -- pp15 and pp16 carry the
SIDIS message.

Upper panel: events at 300 fb^-1 for every row of the FASER dimuon cut flow,
nu_mu + nubar_mu, one marker per set with its own 68% CL PDF uncertainty
from LHAPDF's prescription for that set (replicas for NNPDF, Hessian for the
others); GRV98 is a single member and is drawn without a bar.  Lower panel:
the ratio to NNPDF4.0, with the NNPDF4.0 band and the 7-point scale envelope.
The momentum rows carry the spectrometer acceptance, both muons inside
25 mrad.

THE REGION IS FASER'S CUT FLOW, NOT THE BENCHMARK'S, as for paper plots 12
and 12b: POWHEG-V2 above Q2 > 4 GeV2 on tungsten by isospin (genuine proton
and neutron ladders), folded over the FASER flux on a ten-energy ladder,
with the collaboration's dimuon rows.  There is no W cut.

THE LADDER PREDATES THE ubexcess_correct REGENERATION of 2026-09-19
(tools/powheg_v2_dimuon_ladder.sh ran on 2026-09-07).  That fix moved the
POWHEG-V2 integrals by 0.3% at most and acted at x < 0.011, against a 17%
spread between sets here; regenerating the ladder and its 309-member
reweighting is the cost of making it exact.

GENIE IS NOT A ROW (rule 1b's paper-plot carve-out): this figure is about the
PDF dependence of one NLO calculation.  GENIE's default tune carries GRV98,
which IS a row, and the GENIE-against-POWHEG-V2 dimuon comparison is the cut
flow table of the same section of the paper.

Inputs: results_nu/faser_dimuon_cutflow_pdf.json (analysis/faser_dimuon_cutflow.py pdf).
Usage: analysis/paper_plots/pp13_faser_dimuon_pdf.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SLUG = "pp13_faser_dimuon_pdf"
V1_NUMBER = None
IN_PAPER = True
TITLE = "Dimuon production at FASER: POWHEG-V2 NLO under six PDF sets"
OUTPUT = "pp_faser_dimuon_pdf.png"
RESULTS = "results_nu"

SRC = f"{BASE}/results_nu/faser_dimuon_cutflow_pdf.json"
SURFACE = "white"
WHO = "sum"                     # nu_mu + nubar_mu
NOM = "NNPDF40_nnlo_as_01180"
# the colours of plot_faser_dimuon_pdf.py
COLOUR = {"NNPDF40_nnlo_as_01180": ("NNPDF4.0 NNLO", "#2a78d6", "o"),
          "CT18NNLO": ("CT18 NNLO", "#eb6834", "s"),
          "MSHT20nnlo_as118": ("MSHT20 NNLO", "#3fa66a", "D"),
          # the T = 1 tolerance variant, NNLO (its LHAPDF SetDesc)
          "ATLASpdf21_T1": ("ATLASpdf21 (T=1) NNLO", "#8e5bd0", "^"),
          "ABMP16_5_nnlo": ("ABMP16 NNLO", "#d4a017", "v"),
          "GRV98lo": ("GRV98 LO (no band)", "#777777", "x")}
NNLO_SETS = [k for k in COLOUR if k != "GRV98lo"]
# user, 2026-09-21: no Q2 cut on the first row and no nu_mu + nubar_mu on the
# y axis -- both belong to the dimuon tier's definition, stated in the
# caption; the momentum is that of the SUBLEADING muon.
ROWS = [("cc", r"CC Fiducial"), ("charm", r"charm"), ("charm_mu", r"charm$\to\mu$"),
        ("p20", r"$p_{\mu,\rm subl}>20$"), ("p50", r"$p_{\mu,\rm subl}>50$"),
        ("p100", r"$p_{\mu,\rm subl}>100$")]

CAPTION = (
    "The FASER dimuon cut flow for &nu;<sub>&mu;</sub> + "
    "&nu;&#772;<sub>&mu;</sub> charged-current scattering on tungsten at 300 "
    "fb<sup>&minus;1</sup>, predicted by POWHEG-V2 at NLO showered with "
    "Pythia 8 (Q&sup2; &gt; 4 GeV&sup2;) under six PDF sets: the "
    "charged-current total, events with a charm hadron, with a muon from its "
    "decay, and with both muons inside 25 mrad and the subleading one above "
    "20, 50 and 100 GeV. "
    "Each set carries its own PDF uncertainty at 68% confidence level; "
    "GRV98 has no error members. Lower panel: the ratio to NNPDF4.0, with "
    "the NNPDF4.0 PDF uncertainty (shaded) and the seven-point scale envelope "
    "(hatched).")

MESSAGE = """The dimuon rate at FASER is a measurement of the strange sea. On the charged-current total the six PDF sets agree to 6%, and every NNLO set's own band there is below 1.5%. Once a charm hadron is required they spread over 22%: CT18 at 0.78 of NNPDF4.0, ATLASpdf21 at 0.80, ABMP16 at 0.86 and MSHT20 at 0.91, with the leading-order GRV98 that GENIE's default tune carries at 0.56. The spread survives to the spectrometer rows: with both muons above 100 GeV and inside 25 mrad the NNLO sets predict 14.3 to 17.3 events, a 17% spread, against 4% for NNPDF4.0's own band and a scale envelope of -3.6% to +4.0%. CT18's band is the widest, 15%, and it still does not reach NNPDF4.0: the sets disagree by more than their uncertainties, where the charm-free total does not tell them apart."""

_J = {}


def _d():
    if "d" not in _J:
        with open(SRC) as f:
            _J["d"] = json.load(f)
    return _J["d"]


def _ev(s, k):
    return _d()["sets"][s]["rows"][WHO][k]["events"]


def _ratio(s, k):
    return _ev(s, k) / _ev(NOM, k)


def _rel(s, k, side="err_plus"):
    r = _d()["sets"][s]["rows"][WHO][k]
    return r[side] / r["events"]


def _spread(k, sets):
    v = [_ev(s, k) for s in sets]
    return 1.0 - min(v) / max(v)


def _scale(k):
    sc = _d()["scale"][WHO][k]
    return sc["lo"] / sc["central"] - 1.0, sc["hi"] / sc["central"] - 1.0


CLAIMS = [
    {"what": "the nominal set is NNPDF4.0 and all six sets are present",
     "check": lambda: _d()["nominal_set"] == NOM
     and all(s in _d()["sets"] for s in COLOUR),
     "detail": lambda: ", ".join(_d()["sets"])},
    {"what": "on the charged-current total the six sets agree to 6% and "
             "every NNLO band is below 1.5%",
     "check": lambda: _spread("cc", COLOUR) < 0.06
     and all(_rel(s, "cc") < 0.015 and _rel(s, "cc", "err_minus") < 0.015
             for s in NNLO_SETS),
     "detail": lambda: f"spread {100*_spread('cc', COLOUR):.1f}%; bands " + ", ".join(
         f"{s} {100*_rel(s, 'cc'):.2f}%" for s in NNLO_SETS)},
    {"what": "on the charm row the sets spread over 22%: CT18 0.78, ATLASpdf21 "
             "0.80, ABMP16 0.86, MSHT20 0.91, GRV98 0.56 of NNPDF4.0",
     "check": lambda: abs(_spread("charm", NNLO_SETS) - 0.22) < 0.01
     and all(abs(_ratio(s, "charm") - v) < 0.006 for s, v in
             (("CT18NNLO", 0.78), ("ATLASpdf21_T1", 0.80), ("ABMP16_5_nnlo", 0.86),
              ("MSHT20nnlo_as118", 0.91), ("GRV98lo", 0.56))),
     "detail": lambda: f"NNLO spread {100*_spread('charm', NNLO_SETS):.1f}%; " + ", ".join(
         f"{s} {_ratio(s, 'charm'):.3f}" for s in COLOUR)},
    {"what": "the p > 100 GeV spectrometer row: NNLO sets 14.3 to 17.3 events, "
             "a 17% spread, against a 4% NNPDF4.0 band",
     "check": lambda: abs(min(_ev(s, "p100") for s in NNLO_SETS) - 14.3) < 0.05
     and abs(max(_ev(s, "p100") for s in NNLO_SETS) - 17.3) < 0.05
     and abs(_spread("p100", NNLO_SETS) - 0.17) < 0.005
     and abs(_rel(NOM, "p100") - 0.04) < 0.005,
     "detail": lambda: ", ".join(f"{s} {_ev(s, 'p100'):.2f}" for s in COLOUR)
     + f"; spread {100*_spread('p100', NNLO_SETS):.1f}%, NNPDF band "
       f"{100*_rel(NOM, 'p100'):.2f}%"},
    {"what": "the scale envelope on the p > 100 GeV row is -3.6% / +4.0%",
     "check": lambda: abs(_scale("p100")[0] + 0.036) < 0.0015
     and abs(_scale("p100")[1] - 0.040) < 0.0015,
     "detail": lambda: "%+.2f%% / %+.2f%%" % tuple(100 * v for v in _scale("p100"))},
    {"what": "CT18's band is the widest (about 15%) and does not reach NNPDF4.0 "
             "on the p > 100 GeV row",
     "check": lambda: max(NNLO_SETS, key=lambda s: _rel(s, "p100")) == "CT18NNLO"
     and abs(_rel("CT18NNLO", "p100") - 0.15) < 0.01
     and _ev("CT18NNLO", "p100") * (1 + _rel("CT18NNLO", "p100"))
     < _ev(NOM, "p100") * (1 - _rel(NOM, "p100")),
     "detail": lambda: f"CT18 {_ev('CT18NNLO', 'p100'):.2f} +{100*_rel('CT18NNLO', 'p100'):.1f}%"
                       f" vs NNPDF {_ev(NOM, 'p100'):.2f} -{100*_rel(NOM, 'p100', 'err_minus'):.1f}%"},
    {"what": "the reweighted nominal closes on the unweighted cut flow",
     "check": lambda: all(abs(v - 1.0) < 1e-6
                          for v in _d()["closure_vs_cutflow"].values()),
     "detail": lambda: ", ".join(f"{k} {v:.6f}" for k, v in _d()["closure_vs_cutflow"].items())},
]


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    d = _d()
    nom = d["sets"][NOM]["rows"][WHO]
    yn = np.array([nom[k]["events"] for k, _ in ROWS])
    x = np.arange(len(ROWS))
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.8, 7.6), sharex=True,
                                  gridspec_kw={"height_ratios": (2.2, 1.0), "hspace": 0.06})
    n = len(COLOUR)
    for i, (name, (lab, col, mk)) in enumerate(COLOUR.items()):
        if name not in d["sets"]:
            continue
        r = d["sets"][name]["rows"][WHO]
        off = (i - (n - 1) / 2.0) * 0.11
        y = np.array([r[k]["events"] for k, _ in ROWS])
        ep = np.array([r[k]["err_plus"] for k, _ in ROWS])
        em = np.array([r[k]["err_minus"] for k, _ in ROWS])
        mfc = "none" if name == "GRV98lo" else col
        ax.errorbar(x + off, y, yerr=[em, ep], color=col, marker=mk, ms=6.5, lw=0,
                    elinewidth=1.6, capsize=3, label=tex(lab), mfc=mfc, mew=1.5)
        axr.errorbar(x + off, y / yn, yerr=[em / yn, ep / yn], color=col, marker=mk,
                     ms=6.5, lw=0, elinewidth=1.6, capsize=3, mfc=mfc, mew=1.5)
    bp = np.array([nom[k]["err_plus"] for k, _ in ROWS]) / yn
    bm = np.array([nom[k]["err_minus"] for k, _ in ROWS]) / yn
    for i in range(len(ROWS)):
        axr.fill_between([i - 0.45, i + 0.45], [1 - bm[i]] * 2, [1 + bp[i]] * 2,
                         color="#2a78d6", alpha=0.14, lw=0,
                         label=tex("NNPDF4.0 PDF uncertainty") if i == 0 else None)
        sc = d["scale"][WHO][ROWS[i][0]]
        axr.fill_between([i - 0.45, i + 0.45], [sc["lo"] / sc["central"]] * 2,
                         [sc["hi"] / sc["central"]] * 2, facecolor="none",
                         edgecolor="#555555", hatch="////", lw=0,
                         label=tex("7-point scale envelope") if i == 0 else None)
    axr.axhline(1.0, color="#111111", lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel(tex(r"events at 300 fb$^{-1}$"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"Dimuon Event Yields at FASER"),
                 fontsize=plotstyle.FS_PANEL_TITLE + 4)
    ax.legend(fontsize=plotstyle.FS_LEGEND + 1, frameon=False, loc="upper right", ncol=2)
    axr.set_xticks(x)
    axr.set_xticklabels([tex(l) for _, l in ROWS], fontsize=plotstyle.FS_TICKS)
    axr.set_ylabel(tex("ratio to NNPDF4.0"), fontsize=plotstyle.FS_YLABEL - 1)
    axr.legend(fontsize=plotstyle.FS_LEGEND, frameon=False, loc="upper left", ncol=2)
    axr.set_ylim(0.45, 1.45)
    plotstyle.ticks(ax, labelbottom=False)
    plotstyle.ticks(axr)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
