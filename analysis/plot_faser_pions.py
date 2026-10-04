#!/usr/bin/env python3
"""Single-inclusive charged pions at FASERnu: the predicted spectra at
300 fb^-1 under the emulsion tier, one figure per hadron.

>>> NEUTRINO ONLY, BY USER DECISION (2026-09-07). <<<  "to simplify the
discussion (I know is a rule violation, but accepted) show only neutrino DIS
here, remove the muon DIS in this specific study".  So this study is the one
documented exception to CONVENTIONS.md rule 2b, which otherwise requires whatever
is done for one current to be done for the other; the muon-side spectra are
still EXTRACTED and still in results/, they are simply not drawn.  Every
other study on the tab keeps both currents.

Three panels, generator against generator (no analytic reference exists for
a hadron spectrum; the spread IS the message):
  left    pions per bin of z = E_h / nu, the fragmentation variable;
  centre  pions per bin of lab energy E_h, on a log axis;
  right   pions per bin of transverse momentum to the beam axis.
The vertical axis is a COUNT IN THE BIN, not a density (user, 2026-09-07:
"replace 'pions per unit z' by 'Events per bin' ... this is easier to
interpret"), so the bins are not divided by their widths.  The one place a
density survives is the multiplicity figure of plot_faser_sidis.py, which is
compared against a published dM/dz and has to stay one.
Only pions inside the emulsion track acceptance (tan theta < 0.5) are
counted, in events passing FASER Tier E nested in the region (Q2 > 4 GeV2,
W > 3 GeV, NO y window; user, 2026-09-18), on tungsten by isospin from the
samples.  Lower panels: the ratio to POWHEG-V2, the NLO-matched entry
that anchors the rest of the tab.

Usage: analysis/plot_faser_pions.py [pi|K]
Writes results_nu/cmp_faser_pions_nu_<hadron>.png
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

COLOURS = {"powheg_nu": "#1f5fa9", "powheg": "#1f5fa9",
           "herwig_nlo_full": "#3fa66a", "herwig_nlo_powheg_full": "#3fa66a",
           "sherpa_nlo": "#8e5bd0", "sherpa": "#8e5bd0",
           "genie_lo": "#c0392b", "genie": "#c0392b",
           "genie_nnpdf": "#eb6834"}
REF = {"nu": "powheg_nu", "mu": "powheg"}
HAD_LABEL = {"pi": r"\pi^\pm", "K": r"K^\pm"}
# Tier E in the region (== q4w3_faser_e), the SIDIS study's name for it
SEL = "sidis_e"


def load():
    p = f"{RES}/faser_pions.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run analysis/faser_pions.py rates")
    with open(p) as f:
        return json.load(f)


def panel(ax, axr, d, cur, had, obs, edges, xlabel, logx=False):  # noqa: C901
    gens = d["currents"][cur]["generators"]
    ref = REF[cur]
    ref_y = None
    if ref in gens and SEL in gens[ref]["selections"]:
        ref_y = np.array(gens[ref]["selections"][SEL][f"{had}_emul_{obs}"])
    x = np.array(edges)
    for key, g in gens.items():
        if SEL not in g["selections"]:
            continue
        y = np.array(g["selections"][SEL][f"{had}_emul_{obs}"])
        # GENIE HEDIS is the non-default tune; the default tune is FASER's
        ls = "--" if key == "genie" and cur == "nu" else "-"
        ax.step(x, np.append(y, y[-1]), where="post",
                color=COLOURS.get(key, "#444444"), lw=1.7, ls=ls,
                label=tex(g["label"]))
        if ref_y is not None and axr is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(ref_y > 0, y / ref_y, np.nan)
            axr.step(x, np.append(r, r[-1]), where="post",
                     color=COLOURS.get(key, "#444444"), lw=1.4, ls=ls)
    if logx:
        ax.set_xscale("log")
        if axr is not None:
            axr.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x[0] if not logx else max(x[0], 0.3), x[-1])
    ax.set_ylabel(tex(r"Events per bin, 300 fb$^{-1}$"),
                  fontsize=plotstyle.FS_YLABEL - 1)
    if axr is not None:
        axr.axhline(1.0, color="#111111", lw=1.0, ls="--")
        axr.set_ylim(0.4, 1.8)
        axr.set_ylabel(tex("ratio to POWHEG-V2"),
                       fontsize=plotstyle.FS_YLABEL - 3)
        axr.set_xlabel(tex(xlabel), fontsize=plotstyle.FS_XLABEL)
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(axr)
    else:
        ax.set_xlabel(tex(xlabel), fontsize=plotstyle.FS_XLABEL)
        plotstyle.ticks(ax)


def main():
    # NEUTRINO ONLY (user, 2026-09-07) -- see the module docstring.
    curs = ["nu"]
    hads = [a for a in sys.argv[1:] if a in ("pi", "K")] or ["pi", "K"]
    d = load()
    for cur in curs:
        for had in hads:
            gens = d["currents"][cur]["generators"]
            if not gens:
                print(f"no spectra for {cur}")
                continue
            any_g = next(iter(gens.values()))
            fig, axes = plt.subplots(2, 3, figsize=(18.6, 6.6), sharex="col",
                                     gridspec_kw={"height_ratios": (2.4, 1.0),
                                                  "hspace": 0.06, "wspace": 0.26})
            e_edges = any_g["selections"][SEL]["E_edges_common"]
            panel(axes[0, 0], axes[1, 0], d, cur, had, "z", d["z_edges"], r"$z = E_h/\nu$")
            panel(axes[0, 1], axes[1, 1], d, cur, had, "E", e_edges, r"$E_h$ [GeV]", logx=True)
            panel(axes[0, 2], axes[1, 2], d, cur, had, "pt", d["pt_edges"], r"$p_T$ [GeV]")
            beam = r"$\nu_\mu$ CC" if cur == "nu" else r"$\mu^\pm$ NC"
            axes[0, 0].set_title(tex(r"%s: single-inclusive $%s$ in the emulsion, "
                                     r"FASER Tier E, $Q^2 > 4$ GeV$^2$, $W > 3$ GeV, no $y$ cut"
                                     % (beam, HAD_LABEL[had])),
                                 fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
            axes[0, 0].legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False,
                              loc="lower left")
            out = f"{RES}/cmp_faser_pions_{cur}_{had}.png"
            fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
            plt.close(fig)
            print("wrote", out)


if __name__ == "__main__":
    main()
