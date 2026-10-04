#!/usr/bin/env python3
"""Figures for the report's "Fragmentation functions" tab: charge-summed
pion and kaon FFs, z D_i^{h+ + h-}(z, Q2), from every LHAPDF family, in the
z range of the FASER SIDIS study.

One figure per (hadron, order): columns are the LO nu_mu CC combination and
the flavours u, d, s, c, g; rows are z D at Q2 = 10 GeV2, its ratio to
MAPFF1.0 of the same order, and the same pair at Q2 = 100 GeV2.  Bands are
68% CL.  The z range outside 0.1 < z < 0.8 is shaded.

Inputs: results_nu/ff_comparison.json   (analysis/ff_comparison.py)
Output: results_nu/cmp_ff_<pi|K>_<nlo|nnlo>.png
Usage:  analysis/plot_ff_comparison.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

SRC = f"{BASE}/results_nu/ff_comparison.json"
COLOURS = {"NNFF1.0": "#1f5fa9", "MAPFF1.0": "#c0392b", "NPC23": "#2e8b57",
           "JAM19": "#e08e0b", "JAM20-SIDIS": "#7d3c98", "JAM24": "#6e4b2a"}
TITLES = {"nucc": r"$\nu_\mu$ CC (LO weights)", "u": r"$u$", "d": r"$d$",
          "s": r"$s$", "c": r"$c$", "g": r"$g$"}
HADRON = {"pi": r"\pi^+ + \pi^-", "K": r"K^+ + K^-"}


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    with open(SRC) as f:
        d = json.load(f)
    z = np.array(d["z"])
    zlo, zhi = d["z_window_faser"]
    fls = d["flavours"]
    for h in ("pi", "K"):
        for order in ("nlo", "nnlo"):
            sets = [(n, r) for n, r in d["sets"].items()
                    if r["hadron"] == h and r["order"] == order]
            ref = next(r for n, r in sets if r["family"] == d["reference"][order])
            fig, axes = plt.subplots(4, len(fls), figsize=(2.9 * len(fls), 10.5),
                                     sharex=True,
                                     gridspec_kw={"height_ratios": [2, 1, 2, 1],
                                                  "hspace": 0.08, "wspace": 0.30})
            for k, q2 in enumerate(d["q2"]):
                qk = f"{q2:g}"
                for c, fl in enumerate(fls):
                    ax, rx = axes[2 * k, c], axes[2 * k + 1, c]
                    rc = np.array(ref["values"][qk][fl]["central"])
                    for name, r in sets:
                        v = r["values"][qk][fl]
                        cen, lo, hi = (np.array(v[x]) for x in ("central", "lo", "hi"))
                        col = COLOURS[r["family"]]
                        ax.fill_between(z, lo, hi, color=col, alpha=0.15, lw=0)
                        ax.plot(z, cen, color=col, lw=1.4, label=r["family"])
                        with np.errstate(divide="ignore", invalid="ignore"):
                            rx.fill_between(z, lo / rc, hi / rc, color=col,
                                            alpha=0.15, lw=0)
                            rx.plot(z, cen / rc, color=col, lw=1.4)
                    for a in (ax, rx):
                        a.axvspan(0.0, zlo, color="#000000", alpha=0.06, lw=0)
                        a.axvspan(zhi, 1.0, color="#000000", alpha=0.06, lw=0)
                        a.set_xlim(z[0], z[-1])
                        plotstyle.ticks(a)
                    ax.set_yscale("log")
                    # the window a reader looks at: the FASER z range, down
                    # to 1e-3 of the largest value (a set that vanishes at
                    # large z would otherwise stretch the axis to 1e-6)
                    win = (z >= zlo) & (z <= zhi)
                    top = max(np.max(np.array(r["values"][qk][fl]["central"])[win])
                              for _n, r in sets)
                    bot = min(np.min(np.array(r["values"][qk][fl]["central"])[win])
                              for _n, r in sets)
                    ax.set_ylim(max(bot * 0.5, top * 1e-3), top * 3.0)
                    rx.set_ylim(0.0, 2.0)
                    ax.tick_params(labelbottom=False)
                    if k == 0:
                        rx.tick_params(labelbottom=False)
                    rx.axhline(1.0, color="#111111", lw=0.8)
                    if k == 0:
                        ax.set_title(tex(TITLES[fl]),
                                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
                    if c == 0:
                        ax.set_ylabel(tex(r"$z D(z)$, $Q^2 = %g$ GeV$^2$" % q2),
                                      fontsize=plotstyle.FS_YLABEL - 4)
                        rx.set_ylabel(tex("ratio to %s" % d["reference"][order]),
                                      fontsize=plotstyle.FS_YLABEL - 6)
            for c in range(len(fls)):
                axes[-1, c].set_xlabel(tex(r"$z$"), fontsize=plotstyle.FS_XLABEL)
            axes[0, 0].legend(fontsize=plotstyle.FS_LEGEND - 5, frameon=False,
                              loc="lower left")
            fig.suptitle(tex(r"$%s$ fragmentation functions, %s" % (HADRON[h], order.upper())),
                         fontsize=plotstyle.FS_PANEL_TITLE + 2, x=0.1, ha="left", y=0.93)
            out = f"{BASE}/results_nu/cmp_ff_{h}_{order}.png"
            fig.savefig(out, dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {out}")


if __name__ == "__main__":
    main()
