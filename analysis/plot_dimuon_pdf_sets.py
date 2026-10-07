#!/usr/bin/env python3
"""PDF dependence of the DIMUON rate, all sets, from POWHEG-V2 reweighting.

WHY THIS FIGURE (user, 2026-08-28).  The dimuon-vs-energy figure in the energy
scan carries a PDF band on its POWHEG-V2 row, but only NNPDF4.0's -- the
benchmark's own set.  The question this figure answers is the other one: how
much do the SETS disagree with each other on the dimuon rate, and how does
that compare with any single set's own uncertainty?

It is the dimuon counterpart of "PDF dependence of the CC charm fraction",
computed the same way -- every member of every set, reweighted -- and drawn to
the same recipe so the two can be read side by side.

WHY THE DIMUON TIER IS THE INTERESTING ONE.  It is a charm tag built on a
decay muon, so it weighs the strange density directly, where the inclusive
rate is carried by valence quarks.  Measured here, NNPDF4.0's own band is
4.2-4.8% across the scan against 0.85% on the inclusive cross-section from the
SAME sample and the SAME members -- a fivefold amplification -- while the sets
disagree with one another by about 20%.

THE BAND NEEDED THE SHOWERED EVENTS.  A dimuon does not exist in a Les Houches
file, so the member weights are joined to the showered events on `lhe_index`
(analysis/powheg_selection_members.py).  At 400 GeV the sample is 1M events
and the tier keeps 84 of them, so the weights come from the exact selected
subset -- identical to the full reweighting, proven to 3e-13 by
tools/check_subset_reweight.py, and the difference between seconds and 32
hours.

Usage: plot_dimuon_pdf_sets.py -> results_nu/cmp_dimuon_pdf_sets_E.png
"""
import json
import os
import sys

import numpy as np

import plotstyle                                            # noqa: E402
plotstyle.apply()
from plotstyle import tex                                   # noqa: E402
import matplotlib.pyplot as plt                             # noqa: E402
import matplotlib.ticker as ticker                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams                                                # noqa: E402

SURFACE = "#ffffff"
REFERENCE = "NNPDF40_nnlo_as_01180"
FILES = {400.0: "powheg_selmembers_faser_dimuon_400GeV.json",
         1000.0: "powheg_selmembers_faser_dimuon_1TeV.json",
         4000.0: "powheg_selmembers_faser_dimuon_4TeV.json"}
# same colours as the analytic PDF study, so a reader moving between the two
# figures does not have to relearn them
COLOUR = {"NNPDF40_nnlo_as_01180": ("NNPDF4.0 NNLO", "#2a78d6"),
          "CT18NNLO": ("CT18 NNLO", "#eb6834"),
          "MSHT20nnlo_as118": ("MSHT20 NNLO", "#3fa66a"),
          "ATLASpdf21_T1": ("ATLASpdf21", "#8e5bd0"),
          "ABMP16als118_5_nnlo": ("ABMP16 NNLO", "#d4a017")}   # alpha_s = 0.118 (user, 2026-10-05)


def main():
    data = {}
    for e, fn in FILES.items():
        p = f"{BASE}/results_nu/{fn}"
        if os.path.exists(p):
            with open(p) as f:
                data[e] = json.load(f)
    if not data:
        sys.exit("no dimuon member results -- run "
                 "analysis/powheg_selection_members.py")
    es = sorted(data)
    have_frac = all(
        any("fraction" in r for r in data[e]["sets"].values()) for e in es)
    if not have_frac:
        sys.exit("the member files carry no `fraction` -- rerun "
                 "analysis/powheg_pdf_members.py (for the inclusive "
                 "denominator) and then "
                 "analysis/powheg_selection_members.py")

    # DRAWN LIKE "PDF dependence of the CC charm fraction" (user, 2026-09-01):
    # "upper panel is the dimuon cross-section normalised to the inclusive
    # cross-section, bottom panel the ratio to NNPDF4.0, in both cases showing
    # for each PDF set the corresponding PDF uncertainties".  The bands are
    # SHADED rather than error bars, which is what lets five sets share a
    # panel without their caps overlapping, and it is the recipe the figure
    # above uses -- the two are meant to be read one after the other.
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8.4, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.0})

    def series(key):
        x, y, ep, em = [], [], [], []
        for e in es:
            r = data[e]["sets"].get(key)
            if not r or r.get("fraction") is None:
                continue
            x.append(e)
            y.append(r["fraction"])
            ep.append(r["fraction_err_plus"])
            em.append(r["fraction_err_minus"])
        return (np.array(x), np.array(y), np.array(ep), np.array(em))

    rx, ry, rep, rem = series(REFERENCE)
    if not len(rx):
        sys.exit("the reference set has no fraction points")

    for key, (lab, col) in COLOUR.items():
        x, y, ep, em = series(key)
        if not len(x):
            continue
        ax.plot(x, 100 * y, "-o", color=col, ms=6, lw=1.5, label=tex(lab))
        ax.fill_between(x, 100 * (y - em), 100 * (y + ep),
                        color=col, alpha=0.20, lw=0)
        r = y / ry
        axr.plot(x, r, "-o", color=col, ms=6, lw=1.5)
        axr.fill_between(x, (y - em) / ry, (y + ep) / ry,
                         color=col, alpha=0.20, lw=0)
    # the reference's own band, once, about unity -- as on the charm-fraction
    # figure, so the reader can see whether a set's offset exceeds it
    axr.fill_between(rx, 1 - rem / ry, 1 + rep / ry, color="#888888",
                     alpha=0.25, lw=0, zorder=0)
    axr.axhline(1.0, color="#111111", lw=1.0, ls="--")

    for a in (ax, axr):
        a.set_xscale("log")
        a.set_xlim(min(es) * 0.8, max(es) * 1.25)
        a.grid(alpha=0.25, which="both")
        a.xaxis.set_minor_formatter(ticker.NullFormatter())
    # LOG on the top panel: the dimuon fraction runs from 0.013% at 400 GeV
    # to 1.5% at 4 TeV, because the tier needs a SECOND muon above 100 GeV
    # and at 400 GeV there is barely the energy for one.  On a linear axis the
    # two lower energies are a flat line at zero.
    ax.set_yscale("log")
    ax.set_ylabel(tex(r"$\sigma_{\rm dimuon}/\sigma_{\rm inclusive}$  [%]"),
                  fontsize=plotstyle.FS_YLABEL)
    _h, _l = ax.get_legend_handles_labels()
    ax.legend(_h, _l, fontsize=plotstyle.FS_LEGEND, frameon=True,
              loc="upper left", ncol=2)
    ax.tick_params(labelbottom=False)
    axr.set_ylabel(tex("ratio to") + "\n" + tex("NNPDF4.0 NNLO"),
                   fontsize=plotstyle.FS_YLABEL)
    axr.set_xlabel(tex(r"$E_{\nu}$  [GeV]"), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])

    fig.suptitle(tex("PDF dependence of the dimuon rate") + "\n"
                 + tex("POWHEG-V2 NLO+PS (ZM-VFNS) by LHE reweighting; "
                       "68% CL PDF errors"),
                 fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, 0.015,
             "The dimuon rate as a fraction of the inclusive one, formed "
             "member by member so the correlation between numerator and\n"
             "denominator is kept; the inclusive denominator comes from the "
             "full-sample member weights. The tier amplifies the PDF\n"
             "uncertainty about fivefold over the inclusive rate, and the "
             "sets disagree by more than any one of them claims.",
             ha="center", va="bottom", fontsize=8.4, style="italic",
             color="#5d6470")
    fig.subplots_adjust(top=0.89, bottom=0.145, left=0.125, right=0.97)
    out = f"{BASE}/results_nu/cmp_dimuon_pdf_sets_E.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
