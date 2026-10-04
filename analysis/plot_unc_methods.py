#!/usr/bin/env python3
"""The same uncertainties, two ways: YADISM against POWHEG-V2 reweighting.

WHY THIS FIGURE (user, 2026-08-28).  The PDF and missing-higher-order bands on
this page are ANALYTIC -- YADISM structure functions convoluted with each PDF
member and each scale point.  POWHEG-V2 answers the same two questions on the
same process by reweighting the Les Houches events it already generated.
Nothing is shared: different codes, different perturbative machinery,
different definition of what "the calculation" is.

So the figure is not about the central values -- those are different objects,
and the schemes differ (see below).  It is about whether the RELATIVE
uncertainties agree.  If they do, then "the PDF choice beats the missing
higher orders" rests on two independent calculations rather than one.

WHAT THE PANELS SHOW.
  Left   the relative PDF band per set, at each energy, one method against
         the other.  Points on the diagonal agree.
  Right  the same for the 7-point scale envelope, which is one number per
         energy rather than one per set.

ON THE SCHEMES, which are NOT the same and must not be silently equated.  The
analytic study is FONLL; the POWHEG-V2 samples are ZM-VFNS, because every
production card sets `qmass 0d0` to meet the benchmark's massless-charm
convention -- a RUN SETTING, not a limitation of the code, which is the
massive-charm one.  A relative band is far less sensitive to that than a
central value, which is exactly why this figure compares bands.

GRV98 IS ABSENT FROM THE MONTE CARLO ARM and that is said on the figure: the
installed set declares an LHAPDF id that the index resolves to a different,
absent set, so POWHEG-V2 cannot be asked for it.  It has one member and no
band, so nothing is lost here beyond a row.

Usage: plot_unc_methods.py   -> results_nu/cmp_unc_methods.png
"""
import json
import os
import sys

import plotstyle                                            # noqa: E402
plotstyle.apply()
from plotstyle import tex                                   # noqa: E402
import matplotlib.pyplot as plt                             # noqa: E402
import numpy as np                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams                                                # noqa: E402

SURFACE = "#ffffff"
MC_FILES = {400.0: "powheg_pdf_members_400GeV.json",
            1000.0: "powheg_pdf_members.json",
            4000.0: "powheg_pdf_members_4TeV.json"}
# one marker per energy, one colour per PDF set -- so a reader can see at once
# whether a disagreement belongs to a set or to an energy
EMARK = {400.0: "o", 1000.0: "s", 4000.0: "^"}
ANALYTIC = f"{BASE}/results_nu/pdf_dependence_charm_fraction.json"


def load():
    with open(ANALYTIC) as f:
        an = json.load(f)
    mc = {}
    for e, fn in MC_FILES.items():
        p = f"{BASE}/results_nu/{fn}"
        if os.path.exists(p):
            with open(p) as f:
                mc[e] = json.load(f)
    return an, mc


def main():
    an, mc = load()
    if not mc:
        sys.exit("no POWHEG-V2 member results -- run powheg_pdf_members.py")

    fig, (axp, axs) = plt.subplots(1, 2, figsize=(11.6, 5.6))

    # ---------------- left: the PDF band, set by set, energy by energy ------
    pts = []
    for e, d in sorted(mc.items()):
        ek = f"{e:g}"
        for key, rec in d["sets"].items():
            arec = an["sets"].get(key, {}).get("points", {}).get(ek)
            if not arec or not arec.get("fraction"):
                continue
            a = 100.0 * arec["err_plus"] / arec["fraction"]
            m = 100.0 * rec["rel_plus"]
            col = an["sets"][key]["colour"]
            lab = an["sets"][key]["label"]
            axp.plot(a, m, EMARK[e], color=col, ms=8.5, mfc=col,
                     mec="#33383f", mew=0.8, zorder=3)
            pts.append((a, m, lab, e))
    if not pts:
        sys.exit("nothing in common between the analytic and MC results")
    hi = max(max(p[0] for p in pts), max(p[1] for p in pts)) * 1.12
    axp.plot([0, hi], [0, hi], "-", color="#9aa1a9", lw=1.1, zorder=1)
    axp.set_xlim(0, hi)
    axp.set_ylim(0, hi)
    axp.set_aspect("equal")
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_xlabel(tex("YADISM (FONLL), relative PDF band [%]"))
    axp.set_ylabel(tex("POWHEG-V2 (ZM-VFNS), relative PDF band [%]"))
    axp.set_title(tex("PDF uncertainty on the charm fraction"),
                  fontsize=plotstyle.FS_PANEL_TITLE)

    # legends: colour = set, marker = energy, built by hand so the two
    # dimensions read separately rather than as one list of 15 entries
    from matplotlib.lines import Line2D
    seen, hs = set(), []
    for key, rec in an["sets"].items():
        if key in ("GRV98lo",):
            continue
        if rec["label"] in seen:
            continue
        seen.add(rec["label"])
        hs.append(Line2D([], [], marker="o", ls="none", color=rec["colour"],
                         ms=8, mec="#33383f", mew=0.8, label=tex(rec["label"])))
    l1 = axp.legend(handles=hs, fontsize=plotstyle.FS_LEGEND - 2,
                    loc="upper left", frameon=True)
    axp.add_artist(l1)
    hs2 = [Line2D([], [], marker=EMARK[e], ls="none", color="#5d6470", ms=8,
                  label=tex(f"{e:g} GeV" if e < 1000 else f"{e/1000:g} TeV"))
           for e in sorted(mc)]
    axp.legend(handles=hs2, fontsize=plotstyle.FS_LEGEND - 2,
               loc="lower right", frameon=True)

    d = [abs(p[1] - p[0]) for p in pts]
    axp.text(0.03, 0.60,
             tex(f"{len(pts)} comparisons") + "\n"
             + tex(f"max |diff| {max(d):.2f} pp") + "\n"
             + tex(f"mean {sum(d)/len(d):.2f} pp"),
             transform=axp.transAxes, fontsize=9, color="#33383f",
             va="top", bbox=dict(fc="#f5f6f7", ec="#d0d4d9", lw=0.7))

    # ---------------- right: the scale envelope, one per energy ------------
    es = sorted(mc)
    with open(f"{BASE}/results_nu/mhou_charm_fraction.json") as f:
        mh = json.load(f)
    a_half, m_half = [], []
    for e in es:
        rows = (mh.get("points") or {}).get(f"{e:g}", [])
        fr = [r["fraction"] for r in rows if r.get("fraction")]
        cen = [r["fraction"] for r in rows
               if r.get("xiR") == 1 and r.get("xiF") == 1]
        a_half.append(100.0 * 0.5 * (max(fr) - min(fr)) / cen[0]
                      if fr and cen else np.nan)
        s = mc[e].get("scale") or {}
        m_half.append(100.0 * s.get("rel_half_width", np.nan))
    axs.plot(es, a_half, "-o", color="#2a78d6", lw=1.7, ms=8,
             label=tex("YADISM (FONLL), 7-point"))
    axs.plot(es, m_half, "--s", color="#eb6834", lw=1.7, ms=8,
             label=tex("POWHEG-V2 (ZM-VFNS), 7-point"))
    axs.set_xscale("log")
    axs.set_xticks(es)
    axs.set_xticklabels([f"{e:g}" for e in es])
    axs.set_xlim(min(es)*0.8, max(es)*1.25)
    axs.set_ylim(0, max([v for v in a_half + m_half
                         if np.isfinite(v)] or [1]) * 1.45)
    axs.grid(alpha=0.25, lw=0.6)
    axs.set_xlabel(tex(r"$E_{\nu}$ [GeV]"))
    axs.set_ylabel(tex("scale envelope, half-width [%]"))
    axs.set_title(tex("Missing higher orders"),
                  fontsize=plotstyle.FS_PANEL_TITLE)
    _h, _l = axs.get_legend_handles_labels()
    axs.legend(_h, _l, fontsize=plotstyle.FS_LEGEND - 1, frameon=True,
               loc="upper right")

    fig.suptitle(tex("The same uncertainties computed two ways")
                 + "\n"
                 + tex("analytic structure functions against Monte Carlo "
                       "reweighting, CC charm fraction"),
                 y=1.03, fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, -0.09,
             tex("Nothing is shared between the two: different codes, "
                 "different machinery. Points on the diagonal agree.")
             + "\n"
             + tex("The schemes are NOT the same -- YADISM in FONLL, "
                   "POWHEG-V2 in ZM-VFNS by its qmass = 0 run setting -- "
                   "which a relative band is far less sensitive to than a "
                   "central value.")
             + "\n"
             + tex("GRV98 has no POWHEG-V2 arm: the installed set declares an "
                   "LHAPDF id the index resolves to a different, absent set. "
                   "It has one member and no band."),
             ha="center", va="top", fontsize=8.4, style="italic",
             color="#5d6470")
    out = f"{BASE}/results_nu/cmp_unc_methods.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  {len(pts)} comparisons, max |diff| {max(d):.2f} pp, "
          f"mean {sum(d)/len(d):.2f} pp")


if __name__ == "__main__":
    main()
