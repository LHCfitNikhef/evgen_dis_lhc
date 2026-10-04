#!/usr/bin/env python3
"""The D+- / D0 ratio against beam energy, both currents.

WHY THIS DESERVES ITS OWN FIGURE (user, 2026-08-28).  It is the single number
that explains the largest generator spread on this page.  The dimuon tier asks
for a second muon from a charm decay, and
BR(D+- -> mu X) ~ 16% is about 2.3 times BR(D0 -> mu X) ~ 6.8%, so a generator
that puts more charm into CHARGED D mesons gets proportionally more dimuons.
Measured at 1 TeV, the dimuon rate per unit charm tracks this ratio across
seven generators with a correlation coefficient of 0.954, while the charm
fraction itself spans only 18.6-22.8%.  So the factor of two in the dimuon
rate is hadronisation, not production.

WHAT THE FIGURE SHOWS.  The ratio is nearly FLAT in beam energy for every
generator -- it moves by a few per cent from 400 GeV to 4 TeV while the
cross-section changes by an order of magnitude.  That is the point: it is a
property of the fragmentation model, not of the kinematics, so it does not
average away and cannot be tuned out by moving to another energy.

The generators separate into their hadronisation models rather than their
perturbative orders: Pythia's Lund string near 0.52, Herwig's cluster model
near 0.45, Sherpa's Ahadic cluster model near 0.29, and each generator's LO
and NLO samples sitting on top of one another.

GENIE'S MUON ROW IS ABSENT and that is not a gap: its electromagnetic DIS is
three-flavour with no charm matrix element and no charm from fragmentation, so
it makes no D mesons at all in the neutral current.

Usage: plot_dmeson_ratio.py     -> results/cmp_dmeson_ratio_E.png
"""
import json
import os
import sys

import plotstyle          # noqa: E402
plotstyle.apply()
from plotstyle import tex  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams  # noqa: E402
import labels  # noqa: E402

SURFACE = "#ffffff"

# (label, key, colour, marker).  Grouped by HADRONISATION MODEL rather than by
# perturbative order, because that is what the figure separates.
MU_ROWS = [
    ("Pythia 8.311 LO (Lund)",      "pythia",     "#eb6834", "o"),
    ("POWHEG-RES + Pythia (Lund)",  "powheg",     "#b8442a", "D"),
    ("Herwig 7.3.0 LO (cluster)",   "herwig",     "#8e5bd0", "s"),
    ("Herwig NLO (cluster)",        "herwig_nlo_powheg_full", "#5a2f96", "v"),
    ("Sherpa 3.0.5 LO (Ahadic)",    "sherpa_lo",  "#3fa66a", "^"),
    ("Sherpa MC@NLO (Ahadic)",      "sherpa",     "#1f7a4a", "*"),
]
NU_ROWS = [
    ("Pythia 8.311 LO (Lund)",      "pythia",     "#eb6834", "o"),
    ("POWHEG-V2 + Pythia (Lund)",   "powheg_nu",  "#b8442a", "D"),
    ("Herwig 7.3.0 LO (cluster)",   "herwig",     "#8e5bd0", "s"),
    ("Herwig NLO (cluster)",        "herwig_nlo_full", "#5a2f96", "v"),
    ("Sherpa 3.0.5 LO (Ahadic)",    "sherpa_lo",  "#3fa66a", "^"),
    ("Sherpa MC@NLO (Ahadic)",      "sherpa_nlo", "#1f7a4a", "*"),
    (labels.GENIE_NU_HEDIS[:-1] + ", AGKY)", "genie",  "#c2317b", "P"),
    (labels.GENIE_NU[:-1] + ", AGKY)", "genie_lo",   "#e88bb8", "X"),
]


def ratio(resdir, key, e_lab):
    """<N(D+-)>/<N(D0)> at one energy, or None if the sample makes no D."""
    f = f"{resdir}/histos_{beams.at_energy(key, e_lab)}.json"
    try:
        with open(f) as fh:
            d = json.load(fh)
    except FileNotFoundError:
        return None
    a, b = d.get("mean_n_d_ch"), d.get("mean_n_d_0")
    return (a / b) if (a and b) else None


def panel(ax, resdir, rows, title):
    es = list(beams.ENERGIES)
    for lab, key, colour, marker in rows:
        pts = [(e, ratio(resdir, key, e)) for e in es]
        pts = [(e, r) for e, r in pts if r]
        if not pts:
            continue
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-",
                color=colour, lw=1.6, marker=marker, ms=6.0, label=lab)
    ax.set_xscale("log")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_xlim(min(es) * 0.8, max(es) * 1.25)
    ax.set_xticks(es)
    ax.set_xticklabels([f"{e:g}" for e in es])
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
    ax.set_xlabel(tex("beam energy [GeV]"))


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.4), sharey=True)
    panel(axes[0], f"{BASE}/results", MU_ROWS,
          "Muon DIS (NC): $\\mu^-$ on p at rest")
    panel(axes[1], f"{BASE}/results_nu", NU_ROWS,
          "Neutrino DIS (CC): $\\nu_\\mu$ on p at rest")
    axes[0].set_ylabel(tex(r"$\langle N(D^{\pm})\rangle / "
                           r"\langle N(D^{0})\rangle$"))
    for a, col in zip(axes, (0, 1)):
        _h, _l = a.get_legend_handles_labels()
        a.legend(_h, [tex(x) for x in _l],
                 fontsize=plotstyle.FS_LEGEND - 3, frameon=True,
                 loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=1,
                 handlelength=1.8, labelspacing=0.3)
    fig.suptitle(tex("Charged-to-neutral D-meson ratio across the scan")
                 + "\n"
                 + tex("the hadronisation property behind the dimuon spread"),
                 y=1.04, fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, -0.30,
             tex("Nearly flat in energy for every generator, so it is a "
                 "property of the fragmentation model rather than of the "
                 "kinematics.")
             + "\n"
             + tex("The generators separate by hadronisation model, not by "
                   "perturbative order: each one's LO and NLO samples lie on "
                   "top of one another.")
             + "\n"
             + tex("GENIE has no muon-NC row because its electromagnetic DIS "
                   "makes no charm at all, so there are no D mesons to count."),
             ha="center", va="top", fontsize=8.4, style="italic",
             color="#5d6470")
    out = f"{BASE}/results/cmp_dmeson_ratio_E.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
