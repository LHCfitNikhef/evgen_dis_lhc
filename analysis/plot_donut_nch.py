#!/usr/bin/env python3
"""DONUT's charged-particle multiplicity against the benchmark's generators.

The data and DONUT's own LEPTO-based simulation come from FIG. 9 of
arXiv:0711.0728, read off the paper's vector graphics by
tools/digitise_donut.py.  Our generators are folded with DONUT's FIG. 2
spectrum of interacting neutrinos.

EVERY CURVE IS NORMALISED TO THE SAME NUMBER OF EVENTS, so this is a shape
comparison and nothing else.  It has to be: DONUT's own Monte Carlo carries a
trigger, a scan and a location efficiency that depends on the topology, and
ours carries none of them, so the two normalisations are not comparable
objects.  DONUT normalised its simulation to its data for the same reason.

THE BAND AROUND EACH GENERATOR IS THE EMULSION ACCEPTANCE, not a theory
uncertainty.  DONUT counts tracks reconstructed in emulsion and states no
acceptance, so the prediction is computed with three of them and the envelope
is drawn -- a momentum threshold of 0.1 or 0.3 GeV/c and an angular one of
tan(theta) < 0.5 or < 1.0.  It is worth about one track in the mean, which is
larger than the difference between the two generators, and a comparison that
ignored it would read that as physics.

Usage: analysis/plot_donut_nch.py
Writes results_nu/cmp_donut_nch.png
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plotstyle                                              # noqa: E402

plotstyle.apply()
import matplotlib.pyplot as plt                               # noqa: E402

tex = plotstyle.tex
SURFACE = "white"
RES = f"{BASE}/results_nu"

LABEL = {"genie": "GENIE (GRV98LO), Fe (all channels, FSI)",
         "powheg": "POWHEG-V2 + Pythia 8, Fe (CC DIS)"}
COLOUR = {"genie": "#c0392b", "powheg": "#1f5fa9"}
# the acceptance variants whose envelope is drawn as the band; `raw` is left
# out of it deliberately -- no emulsion sees every track, so it is not a
# candidate for the measurement, only a reference point in the prose
ENVELOPE = ("default", "soft", "wide")


def main():
    p = f"{RES}/donut_nch.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run analysis/donut_nch.py --merge")
    d = json.load(open(p))
    n = np.array(d["nch"], dtype=float)
    data = np.array(d["data"])
    err = np.array(d["data_err"])
    mc = np.array(d["donut_mc"])

    fig = plt.figure(figsize=(8.6, 7.4))
    gs = plotstyle.fused_gridspec(fig, 2, (2.3, 1.0))
    ax, axr = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    edges = np.append(n - 0.5, n[-1] + 0.5)
    ax.errorbar(n, data, yerr=err, fmt="o", color="#111111", ms=5.5,
                lw=1.2, capsize=0, label=tex("DONUT data"), zorder=5)
    ax.step(edges, np.append(mc, mc[-1]), where="post", color="#777777",
            lw=1.5, ls="--", label=tex("DONUT simulation (LEPTO)"))
    for key, g in d["generators"].items():
        v = np.array(g["variants"]["default"]["events"])
        lo = np.minimum.reduce([np.array(g["variants"][k]["events"])
                                for k in ENVELOPE if k in g["variants"]])
        hi = np.maximum.reduce([np.array(g["variants"][k]["events"])
                                for k in ENVELOPE if k in g["variants"]])
        c = COLOUR.get(key, "#444444")
        ax.step(edges, np.append(v, v[-1]), where="post", color=c, lw=1.8,
                label=tex(LABEL.get(key, key)))
        ax.fill_between(edges, np.append(lo, lo[-1]), np.append(hi, hi[-1]),
                        step="post", color=c, alpha=0.20, lw=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(data > 0, v / data, np.nan)
            rl = np.where(data > 0, lo / data, np.nan)
            rh = np.where(data > 0, hi / data, np.nan)
        axr.step(edges, np.append(r, r[-1]), where="post", color=c, lw=1.6)
        axr.fill_between(edges, np.append(rl, rl[-1]), np.append(rh, rh[-1]),
                         step="post", color=c, alpha=0.20, lw=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rmc = np.where(data > 0, mc / data, np.nan)
        rerr = np.where(data > 0, err / data, np.nan)
    axr.step(edges, np.append(rmc, rmc[-1]), where="post", color="#777777",
             lw=1.4, ls="--")
    axr.errorbar(n, np.ones_like(n), yerr=rerr, fmt="o", color="#111111",
                 ms=4.0, lw=1.0, capsize=0, zorder=5)

    ax.set_ylabel(tex("Number of events"), fontsize=plotstyle.FS_YLABEL)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0.0, 1.35 * max(data.max(), mc.max()))
    ax.set_title(tex("DONUT: charged particles at the primary vertex, "
                     "all located events"),
                 fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    ax.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=False, loc="upper right")
    axr.axhline(1.0, color="#111111", lw=0.9, ls="--")
    axr.set_ylim(0.0, 2.4)
    axr.set_xlim(edges[0], edges[-1])
    axr.set_ylabel(tex("ratio to data"), fontsize=plotstyle.FS_YLABEL - 3)
    axr.set_xlabel(tex(r"$n_{\rm ch}$"), fontsize=plotstyle.FS_XLABEL)
    plotstyle.ticks(ax, labelbottom=False)
    plotstyle.ticks(axr)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.94, bottom=0.10)

    out = f"{RES}/cmp_donut_nch.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(out, BASE))


if __name__ == "__main__":
    main()
