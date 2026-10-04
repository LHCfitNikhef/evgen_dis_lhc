#!/usr/bin/env python3
"""The NOMAD dimuon comparison figure: R = sigma_mumu/sigma_CC vs E_nu.

Reads results_nu/nomad_dimuon_fonll.json (written by analysis/nomad_dimuon.py)
and draws the measurement against the calculation, with the ratio beneath.

THE FIGURE IS DRAWN EVEN THOUGH THE COMPARISON DOES NOT WORK YET, and that is
deliberate: the calculation sits 1.4 to 2.2 times above the data with the
right shape, and a disagreement of that size is a result to look at rather
than something to withhold until it is fixed.  The report says so in place.

Usage: analysis/plot_nomad.py  -> results_nu/cmp_nomad_dimuon.png
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import plotstyle                                            # noqa: E402
plotstyle.apply()
from plotstyle import tex                                   # noqa: E402
import matplotlib.pyplot as plt                             # noqa: E402

SRC = f"{BASE}/results_nu/nomad_dimuon_fonll.json"
OUT = f"{BASE}/results_nu/cmp_nomad_dimuon.png"
SURFACE = "#fcfcfb"
TH_COLOUR = "#0072b2"


def main():
    if not os.path.exists(SRC):
        sys.exit(f"{SRC} not found -- run analysis/nomad_dimuon.py first")
    with open(SRC) as f:
        d = json.load(f)
    p = d["points"]
    e = np.array([r["e"] for r in p])
    lo = np.array([r["e_lo"] for r in p])
    hi = np.array([r["e_hi"] for r in p])
    dat = np.array([r["r"] for r in p]) * 1e3
    err = np.array([np.hypot(r["stat"], r["sys"]) for r in p]) * 1e3
    th = np.array([r["r_th"] for r in p]) * 1e3
    thi = np.array([r["r_th_hi"] for r in p]) * 1e3
    tlo = np.array([r["r_th_lo"] for r in p]) * 1e3
    zm = np.array([r.get("r_zm", np.nan) for r in p]) * 1e3
    ff = np.array([r.get("r_ffns", np.nan) for r in p]) * 1e3

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.4, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.0})

    # the bin WIDTH is drawn as the horizontal bar, because the theory is
    # evaluated at the bin centre and a reader should see how wide the bin
    # that centre stands for actually is -- the first one spans 6 to 22 GeV
    ax.errorbar(e, dat, yerr=err, xerr=[e - lo, hi - e], fmt="o", ms=5,
                color="#111111", lw=1.2, capsize=0, zorder=4,
                label=tex("NOMAD dimuon data"))
    ax.plot(e, th, "-", color=TH_COLOUR, lw=1.9, zorder=3,
            label=tex("YADISM NLO (FONLL)"))
    ax.fill_between(e, tlo, thi, color=TH_COLOUR, alpha=0.20, lw=0, zorder=2,
                    label=tex("MHOU (7-point scale)"))
    # >>> THE SAME RATIO WITH TWO OTHER CHARM NUMERATORS, AND THIS IS WHERE
    # THE COMPARISON WAS FAILING. <<<  FONLL is the massless result plus a
    # DAMPED mass correction, and over NOMAD's kinematics the damping factor
    # is 0.03 at Q2 = 2.7 GeV2 -- so the benchmark's own convention is very
    # nearly the massless calculation there, with no charm threshold
    # suppression at all.  Drawing all three makes that visible instead of
    # leaving a factor of two unexplained.
    if np.isfinite(zm).all():
        ax.plot(e, zm, ":", color="#8e5bd0", lw=1.7, zorder=3,
                label=tex("massless charm"))
    if np.isfinite(ff).all():
        ax.plot(e, ff, "--", color="#009e73", lw=1.9, zorder=3,
                label=tex(r"massive charm ($n_f = 3$)"))

    axr.errorbar(e, dat / th, yerr=err / th, xerr=[e - lo, hi - e], fmt="o",
                 ms=5, color="#111111", lw=1.2, capsize=0, zorder=4)
    axr.fill_between(e, tlo / th, thi / th, color=TH_COLOUR, alpha=0.20,
                     lw=0, zorder=2)
    axr.axhline(1.0, color=TH_COLOUR, lw=1.4, zorder=3)
    if np.isfinite(ff).all():
        axr.plot(e, dat / ff, "s", ms=4.5, color="#009e73", zorder=4,
                 label=tex("data / massive"))
    if np.isfinite(zm).all():
        axr.plot(e, dat / zm, "^", ms=4, color="#8e5bd0", zorder=4)

    ax.set_xscale("log")
    ax.set_ylabel(tex(r"$\sigma_{\mu\mu}/\sigma_{\rm CC}$  [$10^{-3}$]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"Dimuon fraction in $\nu_{\mu}$ charged-current DIS"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 1, frameon=True,
              loc="upper left")
    axr.set_ylabel(tex("data / theory"), fontsize=plotstyle.FS_YLABEL - 2)
    axr.set_xlabel(tex(r"$E_{\nu}$  [GeV]"), fontsize=plotstyle.FS_XLABEL)
    for a in (ax, axr):
        a.grid(alpha=0.22, lw=0.6)
        a.set_xlim(lo[0] * 0.9, hi[-1] * 1.05)
    axr.set_ylim(0.3, 1.05)
    axr.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=True,
               loc="lower right")
    fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.10)
    fig.savefig(OUT, dpi=160, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
