#!/usr/bin/env python3
"""The Tier E efficiency factors in z, for the SIDIS tab.

Two figures (user, 2026-09-21: "move the rightmost plot to a separate plot,
adding also the ratio to Tier E selection as reference"):

  cmp_faser_sidis_efficiency.png         the flux-averaged eps(z), nu_e left
      and nu_mu right, pions and kaons.  The line is the MEAN of the five
      generators and the band their envelope -- that spread IS the model
      uncertainty on the factor, because eps is a ratio of two rates
      computed on the same sample and almost everything else cancels in it.
  cmp_faser_sidis_efficiency_energy.png  the energy dependence, POWHEG-V2
      alone (GENIE's default tune stops at 1 TeV, so a five-generator mean
      would change composition along the axis): the Tier E event efficiency
      and eps at two z values, with a lower panel dividing each by the Tier
      E event efficiency.  A calculation done at one beam energy needs eps
      at THAT energy rather than the flux average.

Inputs: results_nu/faser_sidis_efficiency.json
Output: results_nu/cmp_faser_sidis_efficiency.png,
        results_nu/cmp_faser_sidis_efficiency_energy.png
Usage:  analysis/plot_faser_sidis_efficiency.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

SRC = f"{BASE}/results_nu/faser_sidis_efficiency.json"
OUT = f"{BASE}/results_nu/cmp_faser_sidis_efficiency.png"
OUT_E = f"{BASE}/results_nu/cmp_faser_sidis_efficiency_energy.png"
# the energy dependence only: the z panels draw the five-generator mean
REF = "powheg_nu"
# the z range the study quotes: above the z > 0.1 cut, and stopping where the
# region's own yield runs out
Z_MAX = 0.8
COLS = [("nu_e", r"$\nu_e$ CC"), ("nu_mu", r"$\nu_\mu$ CC")]
SPECIES = [("pi", r"$\pi^\pm$", "#1f5fa9"), ("K", r"$K^\pm$", "#c0392b")]


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    with open(SRC) as f:
        d = json.load(f)
    edges = np.array(d["z_edges"])
    hi = int(np.argmax(edges[:-1] >= Z_MAX - 1e-12))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3),
                             gridspec_kw={"wspace": 0.26})

    for col, (fl, beam) in enumerate(COLS):
        ax = axes[col]
        mean = d["flux_averaged_mean"][fl]
        for h, hlab, colour in SPECIES:
            y = np.array(mean[h]["eps_total"])[:hi]
            lo = np.array(mean[h]["eps_total_min"])[:hi]
            up = np.array(mean[h]["eps_total_max"])[:hi]
            ax.fill_between(edges[:hi + 1],
                            np.append(lo, lo[-1]), np.append(up, up[-1]),
                            step="post", color=colour, alpha=0.18, lw=0)
            ax.step(edges[:hi + 1], np.append(y, y[-1]), where="post",
                    color=colour, lw=1.8, label=tex(hlab))
        ax.set_xlim(0.0, Z_MAX)
        ax.set_ylim(0.0, 0.35)
        ax.axvspan(0.0, d.get("z_min", 0.1), color="#000000", alpha=0.07, lw=0)
        ax.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        ax.set_ylabel(tex(r"$\varepsilon(z)$"), fontsize=plotstyle.FS_YLABEL)
        ax.set_title(tex("%s, FASER selection" % beam),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=False,
                  loc="lower left")
        plotstyle.ticks(ax)

    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")
    energy_figure(d, edges)


def energy_figure(d, edges):
    """The energy dependence of the event selection, which the flux average
    hides, with each curve divided by the Tier E event efficiency below."""
    import plotstyle
    from plotstyle import tex
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    pe = d["per_energy"][REF]
    # from the Tier E zero up: POWHEG-V2's ladder reaches 20 GeV for the
    # region-only yield (beams.SIDIS_ENERGIES_LOW), but below 200 GeV the
    # factor is identically zero and the ratio undefined
    es = [e for e in sorted(pe, key=float) if float(e) >= 200.0]
    x = [float(e) for e in es]
    ev = np.array([pe[e]["tier_e_event_efficiency"] for e in es])
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(5.6, 5.6), sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1],
                                              "hspace": 0.06})
    ax.plot(x, ev, "o-", color="#111111", lw=1.7, ms=5.5,
            label=tex("Tier E, events"))
    rx.axhline(1.0, color="#111111", lw=1.2)
    for h, hlab, colour in SPECIES:
        for i, ls in ((2, "-"), (10, "--")):
            y = np.array([pe[e][h]["eps_total"][i] for e in es])
            lab = tex(r"%s, $z \simeq %.3f$" % (hlab, 0.5 * (edges[i] + edges[i + 1])))
            ax.plot(x, y, ls, color=colour, lw=1.5, label=lab)
            rx.plot(x, y / ev, ls, color=colour, lw=1.5)
    ax.set_ylabel(tex(r"efficiency"), fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex("energy dependence, POWHEG-V2 + Pythia 8"),
                 fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    ax.set_ylim(0.0, 0.72)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 4, frameon=False, ncol=1,
              loc="upper left")
    rx.set_ylabel(tex(r"ratio to Tier E"), fontsize=plotstyle.FS_YLABEL - 2)
    rx.set_xlabel(tex(r"$E_\nu$  [GeV]"), fontsize=plotstyle.FS_XLABEL)
    for a in (ax, rx):
        a.set_xscale("log")
        a.set_xticks(x)
        a.set_xticklabels([f"{v:g}" for v in x])
        # the log axis draws its own minor labels on top of ours otherwise
        a.xaxis.set_minor_formatter(ticker.NullFormatter())
        a.xaxis.set_minor_locator(ticker.NullLocator())
        plotstyle.ticks(a)
    ax.tick_params(labelbottom=False)
    fig.savefig(OUT_E, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_E}")


if __name__ == "__main__":
    main()
