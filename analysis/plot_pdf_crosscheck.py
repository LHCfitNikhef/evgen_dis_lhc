#!/usr/bin/env python3
"""Analytic against Monte Carlo: the same PDF and scale study, two ways.

WHY THIS EXISTS (user, 2026-08-27).  The PDF-dependence and MHOU studies on
this page are ANALYTIC: YADISM structure functions convoluted with each PDF
member and each scale point.  POWHEG-V2 can answer the same two questions on
the SAME process by a completely different route -- reweighting the Les
Houches events it already generated -- so putting the two side by side is a
real cross-check rather than a presentation choice.  Nothing is shared between
them: different codes, different perturbative machinery, different definitions
of "the calculation".

WHAT IS AND IS NOT LIKE-FOR-LIKE.

  * The CHARM TAG.  YADISM's numerator is F2_charm, a parton-level object.
    The POWHEG tag is >= 1 outgoing charm QUARK in the LHE, also parton level.
    These match, and deliberately so (user, 2026-08-27): the hadron-level
    `_charmfinal` tag used elsewhere on the page would NOT be comparable to a
    structure function.
  * The SCHEME.  Both are ZM-VFNS here, but for POWHEG-V2 that is a RUN
    SETTING and not a property of the code, and the distinction matters
    because it is the opposite of what the code is for.  POWHEG-V2 is
    `powheg-cmass`: massive charm in the matrix element is its whole purpose,
    and it is the neutrino entry partly for that reason.  Every production
    card here nonetheless sets `qmass 0d0` -- deliberately, to meet the
    benchmark's massless-charm convention -- so THESE SAMPLES are ZM while
    the code is not.  `testrun-nu1TeV-cmass` carries qmass = 1.5 GeV and is
    what a massive-charm comparison would be built from.
    Saying "POWHEG-V2 generates with massless charm" as if it were a
    limitation of the code inverts why POWHEG-V2 is here at all.
  * The ORDER.  YADISM NLO against POWHEG-V2 NLO+PS.  The parton shower is not
    involved in either number, since both are evaluated before showering.
  * NOT like-for-like: POWHEG-V2's reweighting carries five of the analytic
    study's six sets -- GRV98 has no arm, because the installed set declares
    an LHAPDF id that the index resolves to a different, absent set.  Only
    the sets common to both are drawn.
  * ALSO NOT like-for-like as of 2026-08-28: the analytic study moved to
    FONLL, while the POWHEG-V2 samples remain ZM by the qmass = 0 setting
    above.  The RELATIVE uncertainties compared here are far less sensitive
    to that than the central values would be, but a direct comparison of
    absolute fractions now crosses schemes and is not drawn as such.

Usage: plot_pdf_crosscheck.py
Writes results_nu/cmp_pdf_crosscheck_E.png and cmp_mhou_crosscheck_E.png
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

RESULTS_NU = f"{BASE}/results_nu"
SURFACE = "#ffffff"

# the four sets POWHEG's reweighting carries, with the analytic study's colours
COMMON = [("NNPDF40_nnlo_as_01180", "NNPDF4.0 NNLO", "#2a78d6"),
          ("CT18NNLO",              "CT18 NNLO",     "#eb6834"),
          ("MSHT20nnlo_as118",      "MSHT20 NNLO",   "#3fa66a"),
          ("ATLASpdf21_T1",         "ATLASpdf21",    "#8e5bd0")]


def powheg_file(e):
    """The reweighting result for one beam energy, or None."""
    tag = beams.Beams("nu", e).tag
    name = ("nlo_unc_powheg_nu.json" if e == beams.ANCHOR_ENERGY
            else f"nlo_unc_powheg_nu_{tag}.json")
    p = f"{RESULTS_NU}/{name}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def yadism_pdf():
    """{set: {energy: fraction}} from the analytic PDF study."""
    p = f"{RESULTS_NU}/pdf_dependence_charm_fraction.json"
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        d = json.load(f)
    out = {}
    for name, rec in d.get("sets", {}).items():
        pts = rec.get("points") or {}
        out[name] = {float(k): v.get("fraction") for k, v in pts.items()
                     if v.get("fraction")}
    return out


def yadism_mhou():
    """{energy: [ (label, fraction) ]} from the analytic scale study."""
    p = f"{RESULTS_NU}/mhou_charm_fraction.json"
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        d = json.load(f)
    out = {}
    for ek, rows in (d.get("points") or {}).items():
        out[float(ek)] = [(f"({r['xiR']:g}, {r['xiF']:g})", r.get("fraction"))
                          for r in rows]
    return out


def pdf_figure():
    es = list(beams.ENERGIES)
    yad = yadism_pdf()
    pw = {e: powheg_file(e) for e in es}
    if not any(pw.values()):
        print("(no POWHEG reweighting results -- run "
              "analysis/powheg_nlo_uncertainty.py first)")
        return
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.4, 7.0), sharex=True,
                                  gridspec_kw={"height_ratios": [2.1, 1],
                                               "hspace": 0.0})
    for name, lab, col in COMMON:
        ya = [(e, yad.get(name, {}).get(e)) for e in es]
        ya = [(e, v) for e, v in ya if v]
        if len(ya) >= 2:
            ax.plot([p[0] for p in ya], [100*p[1] for p in ya], "-",
                    color=col, lw=1.8, marker="o", ms=5.0,
                    label=f"{lab} (YADISM)")
        po = []
        for e in es:
            d = pw.get(e)
            if not d or not d.get("charm_fraction_by_pdf"):
                continue
            labs = d["pdf"]["labels"]
            if name not in labs:
                continue
            po.append((e, d["charm_fraction_by_pdf"][labs.index(name)]))
        if len(po) >= 2:
            ax.plot([p[0] for p in po], [100*p[1] for p in po], "--",
                    color=col, lw=1.5, marker="s", ms=5.5, mfc="none",
                    label=f"{lab} (POWHEG-V2)")
        # ratio: each calculation against ITS OWN NNPDF4.0, so the panel shows
        # the PDF effect and not the offset between the two calculations
        base_y = yad.get(COMMON[0][0], {})
        for pts, style, mfc in ((ya, "-", col), (po, "--", "none")):
            r = [(e, v / base_y[e]) for e, v in pts if base_y.get(e)] \
                if pts is ya else None
            if pts is not ya:
                r = []
                for e, v in pts:
                    d = pw.get(e)
                    labs = d["pdf"]["labels"]
                    nom = d["charm_fraction_by_pdf"][
                        labs.index(COMMON[0][0])]
                    if nom:
                        r.append((e, v / nom))
            if r and len(r) >= 2:
                axr.plot([p[0] for p in r], [p[1] for p in r], style,
                         color=col, lw=1.5, marker="o" if style == "-" else "s",
                         ms=4.5, mfc=mfc)
    axr.axhline(1.0, color="#33383f", ls=":", lw=1.1)
    for a in (ax, axr):
        a.set_xscale("log")
        a.grid(alpha=0.25, lw=0.6)
        a.set_xlim(min(es)*0.8, max(es)*1.25)
    ax.set_ylabel(tex("charm fraction [%]"))
    ax.tick_params(labelbottom=False)
    axr.set_ylabel(tex("ratio to NNPDF4.0") + "\n"
                   + tex("(within each calculation)"),
                   fontsize=plotstyle.FS_YLABEL)
    axr.set_xlabel(tex(r"$E_{\nu}$ [GeV]"))
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])
    axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    _h, _l = ax.get_legend_handles_labels()
    axr.legend(_h, [tex(x) for x in _l], fontsize=plotstyle.FS_LEGEND - 3,
               frameon=True, ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, -0.32), handlelength=2.2)
    fig.suptitle(tex("PDF dependence of the CC charm fraction: "
                     "analytic against Monte Carlo") + "\n"
                 + tex("YADISM NLO in ZM-VFNS (solid) and POWHEG-V2 NLO by "
                       "LHE reweighting (dashed); parton-level charm"),
                 y=1.02, fontsize=plotstyle.FS_SUPTITLE)
    out = f"{RESULTS_NU}/cmp_pdf_crosscheck_E.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def mhou_figure():
    """The 7-point envelope on the charm fraction, both ways, per energy."""
    es = list(beams.ENERGIES)
    yad = yadism_mhou()
    pw = {e: powheg_file(e) for e in es}
    if not yad or not any(pw.values()):
        print("(no MHOU cross-check inputs)")
        return
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for label, colour, marker in (("YADISM NLO (ZM-VFNS)", "#2a78d6", "o"),
                                  ("POWHEG-V2 (reweighted)", "#eb6834", "s")):
        los, his, xs = [], [], []
        for e in es:
            if label.startswith("YADISM"):
                rows = yad.get(e) or []
                fr = [f for _, f in rows if f]
            else:
                d = pw.get(e)
                fr = d.get("charm_fraction_by_scale") if d else None
            if not fr:
                continue
            nom = fr[0]
            xs.append(e)
            los.append(100*(min(fr)/nom - 1))
            his.append(100*(max(fr)/nom - 1))
        if not xs:
            continue
        ax.plot(xs, his, "-", color=colour, marker=marker, ms=6, lw=1.7,
                label=f"{label}, upper")
        ax.plot(xs, los, "--", color=colour, marker=marker, ms=6, lw=1.7,
                mfc="none", label=f"{label}, lower")
        ax.fill_between(xs, los, his, color=colour, alpha=0.12, lw=0)
    ax.axhline(0.0, color="#33383f", ls=":", lw=1.1)
    ax.set_xscale("log")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_xlim(min(es)*0.8, max(es)*1.25)
    ax.set_xticks(es)
    ax.set_xticklabels([f"{e:g}" for e in es])
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(tex(r"$E_{\nu}$ [GeV]"))
    ax.set_ylabel(tex("7-point scale envelope on the") + "\n"
                  + tex("charm fraction [%]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=True, ncol=2,
              loc="best")
    ax.set_title(tex("Missing higher orders on the CC charm fraction: "
                     "analytic against Monte Carlo"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    out = f"{RESULTS_NU}/cmp_mhou_crosscheck_E.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    pdf_figure()
    mhou_figure()
