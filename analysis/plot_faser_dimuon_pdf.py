#!/usr/bin/env python3
"""The POWHEG-V2 dimuon cut flow at FASER under six PDF sets, with bands.

User, 2026-09-07: "the same event counts but now only for POWHEG-V2 for
different PDF sets, including their PDF uncertainty band: NNPDF4.0, CT18,
MSHT20, ATLASpdf21, ABMP16, and GRV98 (the latter without PDF error band)."

Upper panel: events at 300 fb^-1 for every row of the cut flow,
nu_mu + nubar_mu, one marker per set with its own 68% CL PDF uncertainty
from LHAPDF's prescription for that set (replicas for NNPDF, Hessian for
the others); GRV98 is a single member and is drawn without a bar.  Lower
panel: the ratio to NNPDF4.0, with the NNPDF4.0 band as the grey reference
and the 7-point scale envelope as a hatched band, so the reader sees at
once whether the sets disagree by more than their own uncertainties -- on
the dimuon rows they do, as the PDF-dependence study of the 1 TeV sample
found (memory: dimuon-pdf-amplification), because the dimuon rate is the
strange-quark PDF.

The rows carry the spectrometer acceptance (both muons inside 25 mrad).

Usage: analysis/plot_faser_dimuon_pdf.py
Writes results_nu/cmp_faser_dimuon_pdf.png
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plotstyle  # noqa: E402

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

tex = plotstyle.tex
SURFACE = "white"
RES = f"{BASE}/results_nu"

# the same colours as plot_dimuon_pdf_sets.py, GRV98 added
COLOUR = {"NNPDF40_nnlo_as_01180": ("NNPDF4.0 NNLO", "#2a78d6", "o"),
          "CT18NNLO": ("CT18 NNLO", "#eb6834", "s"),
          "MSHT20nnlo_as118": ("MSHT20 NNLO", "#3fa66a", "D"),
          "ATLASpdf21_T1": ("ATLASpdf21", "#8e5bd0", "^"),
          "ABMP16_5_nnlo": ("ABMP16 NNLO", "#d4a017", "v"),
          "GRV98lo": ("GRV98 LO (no band)", "#777777", "x")}
ROWS = [("cc", r"CC, $Q^2>4$"), ("charm", r"charm"), ("charm_mu", r"charm$\to\mu$"),
        ("p20", r"$p>20$"), ("p50", r"$p>50$"), ("p100", r"$p>100$")]
WHO = "sum"


def main():
    p = f"{RES}/faser_dimuon_cutflow_pdf.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run analysis/faser_dimuon_cutflow.py pdf")
    with open(p) as f:
        d = json.load(f)
    nom = d["sets"][d["nominal_set"]]["rows"][WHO]
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
        ax.errorbar(x + off, y, yerr=[em, ep], color=col, marker=mk, ms=6.5, lw=0,
                    elinewidth=1.6, capsize=3, label=tex(lab),
                    mfc="none" if name == "GRV98lo" else col, mew=1.5)
        yn = np.array([nom[k]["events"] for k, _ in ROWS])
        axr.errorbar(x + off, y / yn, yerr=[em / yn, ep / yn], color=col, marker=mk,
                     ms=6.5, lw=0, elinewidth=1.6, capsize=3,
                     mfc="none" if name == "GRV98lo" else col, mew=1.5)
    # the reference band and the scale envelope, on the ratio panel
    yn = np.array([nom[k]["events"] for k, _ in ROWS])
    bp = np.array([nom[k]["err_plus"] for k, _ in ROWS]) / yn
    bm = np.array([nom[k]["err_minus"] for k, _ in ROWS]) / yn
    for i in range(len(ROWS)):
        axr.fill_between([i - 0.45, i + 0.45], [1 - bm[i]] * 2, [1 + bp[i]] * 2,
                         color="#2a78d6", alpha=0.14, lw=0,
                         label=tex("NNPDF4.0 PDF band") if i == 0 else None)
        sc = d["scale"][WHO][ROWS[i][0]]
        axr.fill_between([i - 0.45, i + 0.45], [sc["lo"] / sc["central"]] * 2,
                         [sc["hi"] / sc["central"]] * 2, facecolor="none",
                         edgecolor="#555555", hatch="////", lw=0,
                         label=tex("7-point scale envelope") if i == 0 else None)
    axr.axhline(1.0, color="#111111", lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel(tex(r"events at 300 fb$^{-1}$, $\nu_\mu+\bar\nu_\mu$"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"POWHEG-V2 NLO dimuons at FASER under six PDF sets, "
                     r"both $\mu$ inside 25 mrad"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=False, loc="upper right", ncol=2)
    axr.set_xticks(x)
    axr.set_xticklabels([tex(l) for _, l in ROWS], fontsize=plotstyle.FS_TICKS)
    axr.set_ylabel(tex("ratio to NNPDF4.0"), fontsize=plotstyle.FS_YLABEL - 1)
    axr.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False, loc="upper left", ncol=2)
    axr.set_ylim(0.45, 1.45)
    plotstyle.ticks(ax, labelbottom=False)
    plotstyle.ticks(axr)
    out = f"{RES}/cmp_faser_dimuon_pdf.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
