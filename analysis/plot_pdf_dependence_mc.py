#!/usr/bin/env python3
"""PDF dependence of the CC charm fraction, computed with POWHEG-V2.

THE MONTE CARLO TWIN of `pdf_dependence.py`'s figure (user, 2026-08-28): the
same quantity, the same PDF sets, the same energies, the same ratio-to-NNPDF4.0
lower panel -- but the numbers come from reweighting POWHEG-V2's Les Houches
events over every PDF member instead of from convoluting YADISM structure
functions.  Drawn to the same recipe on purpose, so the two figures can be put
side by side and read as one comparison.

WHAT DIFFERS FROM THE ANALYTIC FIGURE, and neither difference is cosmetic:

  * THE SCHEME.  YADISM is in FONLL; these samples are ZM-VFNS, because every
    POWHEG-V2 production card sets `qmass 0d0` to meet the benchmark's
    massless-charm convention.  That is a run setting, not a limitation --
    POWHEG-V2 is the massive-charm code.  It shifts the absolute fractions;
    the ratio panel, being a ratio of ratios, very largely cancels it.
  * GRV98 IS ABSENT.  The installed set declares an LHAPDF id that the index
    resolves to a different, absent set, so POWHEG-V2 cannot be asked for it.
    One member, no band, so only a row is lost -- but it is SAID here rather
    than left as a gap the reader has to notice.

THE BAND IS FORMED ON THE FRACTION, member by member, for the same reason the
analytic study does it: numerator and denominator move together under a PDF
variation, and combining their separate uncertainties roughly doubles it.

Usage: plot_pdf_dependence_mc.py -> results_nu/cmp_pdf_dependence_mc_E.png
"""
import json
import os
import sys

import plotstyle                                            # noqa: E402
plotstyle.apply()
from plotstyle import tex                                   # noqa: E402
import matplotlib.pyplot as plt                             # noqa: E402
import matplotlib.ticker as ticker                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SURFACE = "#ffffff"
REFERENCE = "NNPDF40_nnlo_as_01180"
FILES = {400.0: "powheg_pdf_members_400GeV.json",
         1000.0: "powheg_pdf_members.json",
         4000.0: "powheg_pdf_members_4TeV.json"}


def main():
    with open(f"{BASE}/results_nu/pdf_dependence_charm_fraction.json") as f:
        an = json.load(f)          # borrowed only for labels and colours
    data = {}
    for e, fn in FILES.items():
        p = f"{BASE}/results_nu/{fn}"
        if os.path.exists(p):
            with open(p) as f:
                data[e] = json.load(f)
    if not data:
        sys.exit("no POWHEG-V2 member results -- run powheg_pdf_members.py")
    es = sorted(data)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8.4, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.0})

    def series(key):
        x, y, ep, em = [], [], [], []
        for e in es:
            r = data[e]["sets"].get(key)
            if not r:
                continue
            x.append(e); y.append(100*r["fraction"])
            ep.append(100*r["err_plus"]); em.append(100*r["err_minus"])
        return x, y, ep, em

    ref = {e: data[e]["sets"].get(REFERENCE, {}).get("fraction")
           for e in es}
    for key, rec in an["sets"].items():
        x, y, ep, em = series(key)
        if not x:
            continue
        c = rec["colour"]
        ax.errorbar(x, y, yerr=[em, ep], marker="o", ms=6.5, lw=1.6,
                    color=c, capsize=3, elinewidth=1.1,
                    label=tex(rec["label"]))
        rr = [(e, data[e]["sets"][key]["fraction"]/ref[e],
               100*data[e]["sets"][key]["err_plus"]/100/ref[e],
               100*data[e]["sets"][key]["err_minus"]/100/ref[e])
              for e in x if ref.get(e)]
        if rr:
            axr.errorbar([p[0] for p in rr], [p[1] for p in rr],
                         yerr=[[p[3] for p in rr], [p[2] for p in rr]],
                         marker="o", ms=6.5, lw=1.6, color=c, capsize=3,
                         elinewidth=1.1)
    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="--")
    for a in (ax, axr):
        a.set_xscale("log")
        a.grid(alpha=0.25, lw=0.6)
        a.set_xlim(min(es)*0.8, max(es)*1.25)
    ax.set_ylabel(tex(r"$\sigma_{\rm charm}/\sigma_{\rm inclusive}$  [%]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.tick_params(labelbottom=False)
    _h, _l = ax.get_legend_handles_labels()
    ax.legend(_h, _l, fontsize=plotstyle.FS_LEGEND, frameon=True,
              loc="upper left")
    axr.set_ylabel(tex("ratio to") + "\n" + tex("NNPDF4.0 NNLO"),
                   fontsize=plotstyle.FS_YLABEL)
    axr.set_xlabel(tex(r"$E_{\nu}$ [GeV]"), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])
    axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    fig.suptitle(tex("PDF dependence of the CC charm fraction")
                 + "\n"
                 + tex("POWHEG-V2 NLO+PS (ZM-VFNS) by LHE reweighting, "
                       "$Q^2>4$ GeV$^2$, $0.2<y<0.9$; 68% CL PDF errors"),
                 fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, 0.015,
             tex("The Monte Carlo twin of the analytic figure: same quantity, "
                 "same sets, same energies, computed by reweighting POWHEG-V2 "
                 "events") + "\n"
             + tex("over every PDF member rather than by convoluting YADISM "
                   "structure functions. The band is formed on the FRACTION, "
                   "member by member.") + "\n"
             + tex("GRV98 has no arm here: the installed set declares an "
                   "LHAPDF id the index resolves to a different, absent set. "
                   "It has one member and no band."),
             ha="center", va="bottom", fontsize=8.4, style="italic",
             color="#5d6470")
    fig.subplots_adjust(top=0.89, bottom=0.135, left=0.115, right=0.97)
    out = f"{BASE}/results_nu/cmp_pdf_dependence_mc_E.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
