#!/usr/bin/env python3
"""POWHEG-V2 SIDIS in the benchmark region alone, for the SIDIS tab's
"Region yields" view (user, 2026-09-21): the numbers an analytic NLO
calculation with fragmentation functions is benchmarked against.

  cmp_faser_sidis_region_yields.png   flux-folded counts per z bin at the
      study's luminosity, nu_e left and nu_mu right, pions and kaons: the
      region (Q2 > 4, W > 3, nothing else) as solid lines, and the same times
      the five-generator eps(z) -- the FASERnu yield -- dashed, with the
      eps envelope as its band.
  cmp_faser_sidis_region_dsigma.png   dsigma/dz per nucleon at every energy
      of the SIDIS ladder, pions left and kaons right: the fixed-energy
      comparison, before any flux.

Inputs: results_nu/faser_sidis_region_yields.json,
        results_nu/faser_sidis_efficiency.json
Usage:  analysis/plot_faser_sidis_region_yields.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

SRC = f"{BASE}/results_nu/faser_sidis_region_yields.json"
EFF = f"{BASE}/results_nu/faser_sidis_efficiency.json"
OUT = f"{BASE}/results_nu/cmp_faser_sidis_region_yields.png"
OUT_E = f"{BASE}/results_nu/cmp_faser_sidis_region_dsigma.png"
COLS = [("nu_e", r"$\nu_e$ CC"), ("nu_mu", r"$\nu_\mu$ CC")]
SPECIES = [("pi", r"$\pi^\pm$", "#1f5fa9"), ("K", r"$K^\pm$", "#c0392b")]


def _st(y):
    return np.append(y, y[-1])


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    with open(SRC) as f:
        d = json.load(f)
    with open(EFF) as f:
        eff = json.load(f)
    z = np.array(d["z_edges"])

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3),
                             gridspec_kw={"wspace": 0.26})
    for ax, (fl, beam) in zip(axes, COLS):
        r = d["flux_folded"][fl]
        for h, hlab, colour in SPECIES:
            y = np.array(r[h])
            m = eff["flux_averaged_mean"][fl][h]
            ax.step(z, _st(y), where="post", color=colour, lw=1.8,
                    label=tex(f"{hlab}, region"))
            ax.fill_between(z, _st(y * np.array(m["eps_total_min"])),
                            _st(y * np.array(m["eps_total_max"])),
                            step="post", color=colour, alpha=0.18, lw=0)
            ax.step(z, _st(y * np.array(m["eps_total"])), where="post",
                    color=colour, lw=1.5, ls="--",
                    label=tex(rf"{hlab}, $\times\,\varepsilon(z)$"))
        ax.set_yscale("log")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(1e-1, 1e5)
        ax.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        ax.set_ylabel(tex(r"hadrons per $z$ bin, %g fb$^{-1}$" % d["lumi_fb"]),
                      fontsize=plotstyle.FS_YLABEL)
        ax.set_title(tex(f"{beam}, POWHEG-V2"),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False,
                  loc="upper right")
        plotstyle.ticks(ax)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")

    fe = d["fixed_energy"]
    es = sorted(fe, key=float)
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3),
                             gridspec_kw={"wspace": 0.26})
    for ax, (h, hlab, _c) in zip(axes, SPECIES):
        for k, e in enumerate(es):
            y = np.array(fe[e][h]["dsigma_pb"]) / np.diff(z)
            ax.step(z, _st(y), where="post", lw=1.5,
                    color=cmap(k / max(len(es) - 1, 1)),
                    label=tex(f"{float(e):g} GeV"))
        ax.set_yscale("log")
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        ax.set_ylabel(tex(r"$d\sigma/dz$  [pb per nucleon]"),
                      fontsize=plotstyle.FS_YLABEL)
        ax.set_title(tex(rf"$\nu_\mu$ CC $\to$ {hlab}, POWHEG-V2"),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.legend(fontsize=plotstyle.FS_LEGEND - 4, frameon=False,
                  loc="lower left", ncol=2)
        plotstyle.ticks(ax)
    fig.savefig(OUT_E, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_E}")


if __name__ == "__main__":
    main()
