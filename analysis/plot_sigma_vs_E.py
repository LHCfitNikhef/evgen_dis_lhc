#!/usr/bin/env python3
"""Inclusive fiducial cross-section versus beam energy, both currents.

Two columns -- muon NC and neutrino CC -- each with a ratio panel underneath
normalising every generator to the YADISM reference at the SAME energy and
the SAME perturbative order.

WHY THE RATIO PANEL IS THE POINT.  sigma_fid rises by more than a factor of
two across the scan, so on the top panel every generator lies on top of every
other and a 2% disagreement is invisible.  The ratio is where the benchmark
actually lives: it shows whether a generator's offset from the analytic
reference is CONSTANT in energy (a normalisation or scheme effect) or GROWS
(a kinematic one).  Pythia's known massless-charm deficit, for instance, is
flat at about -2% on the muon side, which is what makes it a scheme problem
rather than a scan artefact.

Each generator is compared against its OWN order: the LO generators against
YADISM LO, POWHEG-RES and the NLO rows against YADISM NLO.  Comparing an NLO
sample to the LO reference would fold the K-factor into what is meant to be a
closure test.

Usage: plot_sigma_vs_E.py            -> results/cmp_sigma_vs_E.png
"""
import json
import os
import sys

import plotstyle          # noqa: E402 -- the house style (analysis/plotstyle.py)
plotstyle.apply()         # must precede the pyplot import
from plotstyle import tex  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams  # noqa: E402
import labels  # noqa: E402  -- the one display-name table
from labels import strip_order as _strip_order  # noqa: E402

SURFACE = "#ffffff"

# (label, result key, YADISM reference key, colour, marker)
# The reference is the same-order YADISM row, per the note above.
MU_ROWS = [
    ("Pythia 8.311 LO",      "pythia",     "yadism",     "#eb6834", "o"),
    ("Herwig 7.3.0 LO",      "herwig",     "yadism",     "#8e5bd0", "s"),
    ("Sherpa 3.0.5 LO",      "sherpa_lo",  "yadism",     "#3fa66a", "^"),
    ("POWHEG-RES + Pythia",  "powheg",     "yadism_nlo", "#2b6cb0", "D"),
    # the second, independent NLO matching -- the figure showed one NLO curve
    # where the benchmark has two, which is the comparison it exists to make
    ("Herwig NLO (POWHEG)",  "herwig_nlo_powheg_full", "yadism_nlo",
     "#7a3fb0", "v"),
    ("Sherpa 3.0.5 MC@NLO",  "sherpa",     "yadism_nlo", "#2a78d6", "*"),
    # the NNPDF4.0 variant of the FASER tune: it differs from the row below in
    # what GENIE does, not in which parton densities it uses, which is the
    # only way to separate the tune's model from its PDF
    (labels.GENIE_MU_NNPDF,  "genie_nnpdf", "yadism",    "#8c1f52", "d"),
    # GENIE's MUON entry is EMDIS with tune G18_02a -- Bodek-Yang on GRV98LO,
    # a LEADING-ORDER calculation -- so its reference is YADISM LO, the same
    # one the cross-section tables and the LO blurb use.  Dividing it by the
    # NLO reference here put it at 0.902 in this figure while every table on
    # the page said 0.831: the same sample, two different denominators, and
    # the figure flattered it.  (The NEUTRINO side is different and already
    # correct: HEDIS/BGR18 is NLO, the classic Bodek-Yang path is LO, and the
    # two rows below carry the matching references.)
    (labels.GENIE_MU,        "genie",      "yadism",     "#c2317b", "P"),
]
NU_ROWS = [
    ("Pythia 8.311 LO",      "pythia",     "yadism",     "#eb6834", "o"),
    ("Herwig 7.3.0 LO",      "herwig",     "yadism",     "#8e5bd0", "s"),
    ("Sherpa 3.0.5 LO",      "sherpa_lo",  "yadism",     "#3fa66a", "^"),
    ("POWHEG-V2 + Pythia",   "powheg_nu",  "yadism_nlo", "#2b6cb0", "D"),
    ("Herwig NLO (POWHEG)",  "herwig_nlo_full", "yadism_nlo", "#7a3fb0", "v"),
    ("Sherpa 3.0.5 MC@NLO",  "sherpa_nlo",  "yadism_nlo", "#2a78d6", "*"),
    (labels.GENIE_NU_HEDIS,  "genie",      "yadism_nlo", "#c2317b", "P"),
    (labels.GENIE_NU,        "genie_lo",  "yadism",     "#e88bb8", "X"),
    (labels.GENIE_NU_NNPDF,  "genie_nnpdf", "yadism",   "#8c1f52", "d"),
]

# Which perturbative order each row belongs to.  Kept as a lookup on the
# result key rather than a sixth tuple field, so DIMUON_ROWS -- which unpacks
# five -- keeps working and there is one place to answer the question.
# GENIE: the classic Bodek-Yang path (GRV98LO, and its NNPDF4.0 variant) is
# LO; HEDIS/BGR18 is NLO.  Same split the report tables use.
ROW_ORDER = {
    "pythia": "LO", "herwig": "LO", "sherpa_lo": "LO", "mg5_me": "LO",
    "genie_lo": "LO", "genie_nnpdf": "LO",
    "powheg": "NLO", "powheg_nu": "NLO", "sherpa": "NLO",
    "sherpa_nlo": "NLO", "herwig_nlo_powheg_full": "NLO",
    "herwig_nlo_full": "NLO",
}


def rows_at(rows, order, resdir):
    """The rows of `rows` belonging to `order`.

    GENIE's key is "genie" on BOTH sides but sits at a different order on
    each -- EMDIS/Bodek-Yang (LO) under results/, HEDIS/BGR18 (NLO) under
    results_nu/ -- so it cannot go in ROW_ORDER and is resolved here from the
    directory.  That ambiguity has already produced one wrong denominator in
    this figure (see the MU_ROWS note above), so it is spelled out rather
    than defaulted.
    """
    out = []
    for r in rows:
        key = r[1]
        o = ("NLO" if resdir.endswith("results_nu") else "LO") \
            if key == "genie" else ROW_ORDER.get(key)
        if o is None:
            raise KeyError(f"no perturbative order recorded for row {key!r} "
                           f"-- add it to ROW_ORDER")
        if o == order:
            out.append(r)
    return out
# The analytic curves.  (label, key, colour, style, which figure, ratio ref)
#
# BOTH SCHEMES ARE SHOWN (user, 2026-08-27).  ZM-VFNS is the benchmark
# convention -- massless charm and bottom -- and FONLL restores the mass
# effects in the heavy-quark part while leaving the light part alone.  Where a
# FONLL row carries a ratio reference it is the SAME-ORDER ZM row, so the
# ratio panel reads directly as the size of the scheme change; that is the
# quantity that has to vanish at high energy, where the masses stop mattering.
#
# FONLL HAS NO LO COUNTERPART: the massive contribution gamma* g -> c cbar
# first enters at O(alpha_s), so the LO figure carries ZM alone.
# The NNLO ladder is drawn on the NLO figure -- it is the next order up, and
# there is no NNLO generator for it to sit beside.
REF_ROWS = [
    ("YADISM LO, ZM-VFNS",   "yadism",     "#111111", "-",  "LO",  None),
    ("YADISM NLO, ZM-VFNS",  "yadism_nlo", "#111111", "--", "NLO", None),
    ("YADISM NNLO, ZM-VFNS", "yadism_nnlo","#111111", ":",  "NLO", None),
    ("YADISM NLO, FONLL",  "yadism_nlo_fonll_damp",
     "#c0392b", "--", "NLO", "yadism_nlo"),
    ("YADISM NNLO, FONLL", "yadism_nnlo_fonll_damp",
     "#c0392b", ":",  "NLO", "yadism_nnlo"),
]


def sigma(resdir, key, e_lab, sel="inclusive"):
    """(sigma_fid, err) in pb at one energy, or None.

    None covers two cases that must both keep a point off the figure: the
    sample was not run here, and the sample WAS run but the analysis flagged
    it `stat_insufficient` -- a handful of events, or a negative cross-section
    from incomplete negative-weight cancellation.  Plotting the second would
    put a meaningless number on an axis, which is exactly what the flag
    exists to prevent, so the flag is honoured here and not merely stored.

    `sel` names the event selection and is applied the same way everywhere
    else: the inclusive region keeps the bare key, a tier appends `_<tier>`,
    and beams.at_energy() then appends the energy tag.
    """
    name = key if sel == "inclusive" else f"{key}_{sel}"
    f = f"{resdir}/histos_{beams.at_energy(name, e_lab)}.json"
    try:
        d = json.load(open(f))
    except FileNotFoundError:
        return None
    if d.get("stat_insufficient"):
        return None
    s = d.get("sigma_fid_pb")
    return None if s is None else (s, d.get("sigma_fid_err_pb") or 0.0)


# The dimuon tier, NEUTRINO ONLY: FASER sees a neutrino beam, so an
# opposite-sign dimuon charm tag is a measurable signal only in the charged
# current.  The NC version exists in the results but describes a beam the
# experiment does not have.
DIMUON_ROWS = [r for r in NU_ROWS]

# The dimuon ratio panel is normalised to the tune FASER ACTUALLY RUNS, not to
# an analytic curve: the tier is defined on DECAY muons, so it lives past the
# hadronisation and decay stages and there is no YADISM counterpart to divide
# by.  G18_02a is therefore the only reference that answers the question a
# reader of this figure has -- how does each matching compare with what the
# experiment uses today.
DIMUON_REF_KEY = "genie_lo"
DIMUON_REF_COLOUR = "#e88bb8"


def dimuon_pdf_band(e):
    """(rel_plus, rel_minus) PDF band on the POWHEG-V2 dimuon rate, or None.

    From `powheg_selection_members.py`: the LHE member weights joined to the
    SHOWERED events on lhe_index, because the dimuon tier cuts on a decay muon
    that no Les Houches file contains.  NNPDF4.0's own band -- the benchmark's
    reference set -- not the spread between sets, which is a different and
    much larger number.

    THE BAND IS ~5x THE INCLUSIVE ONE (4.2-4.8% against 0.85%), which is the
    point of drawing it: the dimuon tier is a charm tag, so it weighs the
    strange density directly instead of the valence quarks that carry the
    inclusive rate.
    """
    tag = beams.Beams("nu", e).tag
    p = f"{BASE}/results_nu/powheg_selmembers_faser_dimuon_{tag}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    r = (d.get("sets") or {}).get("NNPDF40_nnlo_as_01180")
    if not r or r.get("rel_plus") is None:
        return None
    return r["rel_plus"], r["rel_minus"]


def dimuon_figure():
    """sigma, selection efficiency and ratio to GENIE versus beam energy.

    The companion figure to the two dimuon tables on the energy-scan tab.
    The efficiency panel is the one to read for the generators themselves:
    numerator and denominator come from the SAME sample, so every
    normalisation convention cancels and the generators are directly
    comparable.  The panel beneath it repeats the top panel against the FASER
    tune, which is the comparison that decides whether a difference matters
    for the experiment.
    """
    resdir = f"{BASE}/results_nu"
    es = list(beams.ENERGIES)
    fig, (ax_top, ax_eff, ax_rat) = plt.subplots(
        3, 1, figsize=(8.2, 9.6), sharex=True,
        gridspec_kw={"height_ratios": [1.55, 1.0, 1.0], "hspace": 0.0})

    # the denominator of the bottom panel, gathered before the loop so a
    # missing GENIE point drops that ENERGY from the ratio rather than
    # silently rescaling it against some other energy's number
    ref = {}
    for e in es:
        v = sigma(resdir, DIMUON_REF_KEY, e, "faser_dimuon")
        if v and v[0] > 0:
            ref[e] = v

    drew = False
    ratio_vals = []
    for label, key, _ref, colour, marker in DIMUON_ROWS:
        pts, eff, rat = [], [], []
        for e in es:
            v = sigma(resdir, key, e, "faser_dimuon")
            base = sigma(resdir, key, e, "inclusive")
            if not v:
                continue
            pts.append((e, v[0] / 1e-3, v[1] / 1e-3))     # pb -> fb
            if base and base[0] > 0:
                eff.append((e, 100.0 * v[0] / base[0],
                            100.0 * v[1] / base[0]))
            if e in ref:
                # the NUMERATOR error only; the reference's own error is the
                # band drawn around 1 below, and adding it to every point too
                # would double-count the one quantity common to them all
                rat.append((e, v[0] / ref[e][0], v[1] / ref[e][0]))
                ratio_vals.append(v[0] / ref[e][0])
        if not pts:
            continue
        drew = True
        ax_top.errorbar([p[0] for p in pts], [p[1] for p in pts],
                        yerr=[p[2] for p in pts], marker=marker, ms=6.0,
                        color=colour, lw=1.4, capsize=2.5, label=label)
        # THE PDF BAND, on the POWHEG-V2 row only -- it is the only sample
        # whose LHE weights have been joined to the showered events, which the
        # dimuon tier requires because it cuts on a decay muon.  Drawn as a
        # shaded envelope so it reads as a systematic rather than as the
        # Monte Carlo error bars beside it.
        if key == "powheg_nu":
            bp = [(e, v, dimuon_pdf_band(e)) for e, v, _err in pts]
            bp = [(e, v, b) for e, v, b in bp if b]
            if len(bp) >= 2:
                ax_top.fill_between(
                    [q[0] for q in bp],
                    [q[1] * (1 - q[2][1]) for q in bp],
                    [q[1] * (1 + q[2][0]) for q in bp],
                    color=colour, alpha=0.20, lw=0, zorder=1,
                    label=tex("POWHEG-V2 PDF band (NNPDF4.0, 68% CL)"))
                rb = [(e, r, dimuon_pdf_band(e))
                      for e, r, _er in rat]
                rb = [(e, r, b) for e, r, b in rb if b]
                if len(rb) >= 2:
                    ax_rat.fill_between(
                        [q[0] for q in rb],
                        [q[1] * (1 - q[2][1]) for q in rb],
                        [q[1] * (1 + q[2][0]) for q in rb],
                        color=colour, alpha=0.20, lw=0, zorder=1)
        if eff:
            ax_eff.errorbar([p[0] for p in eff], [p[1] for p in eff],
                            yerr=[p[2] for p in eff], marker=marker, ms=6.0,
                            color=colour, lw=1.4, capsize=2.5)
        if rat:
            ax_rat.errorbar([p[0] for p in rat], [p[1] for p in rat],
                            yerr=[p[2] for p in rat], marker=marker, ms=6.0,
                            color=colour, lw=1.4, capsize=2.5)
    if not drew:
        plt.close(fig)
        print("(no dimuon results to plot)")
        return

    # the reference's own statistical error, as a band around 1.  The dimuon
    # tier keeps a few hundred events even at 4 TeV, so this band is wide
    # enough that leaving it out would make agreement look like disagreement.
    if ref:
        xs = sorted(ref)
        rel = [ref[e][1] / ref[e][0] for e in xs]
        ax_rat.fill_between(xs, [1 - r for r in rel], [1 + r for r in rel],
                            color=DIMUON_REF_COLOUR, alpha=0.22, lw=0,
                            zorder=0)
    ax_rat.axhline(1.0, color=DIMUON_REF_COLOUR, linestyle="--", lw=1.3,
                   zorder=1)

    for ax in (ax_top, ax_eff, ax_rat):
        ax.set_xscale("log")
        # LINEAR y throughout the energy-scan figures (user, 2026-08-26).
        # NOTE for this one: the dimuon rate spans three orders of magnitude
        # across the scan (0.13 fb at 400 GeV to 183 fb at 4 TeV), so on a
        # linear axis the 400 GeV points sit on the baseline.  The efficiency
        # panel below carries the same information without that compression.
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_xlim(min(es) * 0.8, max(es) * 1.25)
    ax_top.set_ylabel(r"$\sigma_{\rm fid}$" + tex("(dimuon) [fb]"))
    ax_top.tick_params(labelbottom=False)
    ax_eff.tick_params(labelbottom=False)
    ax_eff.set_ylabel(tex("selection efficiency") + "\n"
                      + tex("(% of inclusive)"),
                      fontsize=plotstyle.FS_YLABEL)
    ax_rat.set_ylabel(tex("ratio to") + "\n" + tex("GENIE (GRV98LO)"),
                      fontsize=plotstyle.FS_YLABEL)
    # AUTO-RANGED, as in plot_charm_ratio.py: a fixed window either clips a
    # generator that is genuinely far off or squashes the rest onto one line.
    if ratio_vals:
        lo_r, hi_r = min(ratio_vals), max(ratio_vals)
        pad = max(0.06 * (hi_r - lo_r), 0.03)
        ax_rat.set_ylim(lo_r - pad, hi_r + pad)
    ax_rat.set_xlabel(tex("neutrino beam energy [GeV]"))
    ax_rat.set_xticks(es)
    ax_rat.set_xticklabels([f"{e:g}" for e in es])
    ax_rat.xaxis.set_minor_formatter(ticker.NullFormatter())
    # the legend sits BELOW the panels: with three stacked panels there is no
    # corner of the top one it can occupy without covering the 4 TeV points,
    # which are the ones the figure is about.
    _h, _l = ax_top.get_legend_handles_labels()
    ax_rat.legend(_h, [tex(x) for x in _l],
                  fontsize=plotstyle.FS_LEGEND - 0.8, frameon=False, ncol=3,
                  # centred on the FIGURE, not on the axes: the axes box is
                  # inset from the right, and the third column carries the
                  # longest label in the benchmark, which ran off the page.
                  loc="upper center", bbox_to_anchor=(0.455, -0.30),
                  handletextpad=0.5, columnspacing=1.0)
    fig.suptitle("Dimuon signal tier versus beam energy \u2014 "
                 "CC $\\nu_\\mu$ DIS\n"
                 "Tier S on both muons, opposite sign "
                 "($E' > 100$ GeV, $\\theta' < 25$ mrad)", fontsize=11.5)
    fig.text(0.5, 0.012,
             "Points a sample cannot measure are omitted, not drawn at zero: "
             "at 400 GeV the tier keeps only a handful of events.\n"
             "The efficiency rises two orders of magnitude across the scan "
             "because the tier requires both muons above 100 GeV.\n"
             "The band on the bottom panel is the statistical error of the "
             "GENIE reference itself.",
             ha="center", va="bottom", fontsize=8.2, style="italic",
             color="#5d6470")
    out = f"{resdir}/cmp_dimuon_vs_E.png"
    fig.subplots_adjust(top=0.905, bottom=0.235, left=0.13, right=0.97)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def panel(ax_top, ax_rat, resdir, rows, title, unit, div, order):
    es = list(beams.ENERGIES)
    # analytic references first, so the generators draw over them
    for label, key, colour, style, fig_order, rref in REF_ROWS:
        if fig_order != order:
            continue
        pts = [(e, sigma(resdir, key, e)) for e in es]
        pts = [(e, v[0] / div) for e, v in pts if v]
        if len(pts) < 2:
            continue
        ax_top.plot([p[0] for p in pts], [p[1] for p in pts], style,
                    color=colour, lw=1.3, label=_strip_order(label, order),
                    zorder=1)
        # A FONLL row also goes in the RATIO panel, against the same-order ZM
        # row, where it reads as the size of the scheme change itself.  Drawn
        # as a line with no markers so it stays visibly analytic beside the
        # generators' points.
        if rref:
            rr = [(e, sigma(resdir, key, e), sigma(resdir, rref, e))
                  for e in es]
            rr = [(e, v[0] / r[0]) for e, v, r in rr if v and r]
            if len(rr) >= 2:
                ax_rat.plot([p[0] for p in rr], [p[1] for p in rr], style,
                            color=colour, lw=1.3, zorder=2)
    for label, key, ref, colour, marker in rows:
        pts, rat = [], []
        for e in es:
            v, r = sigma(resdir, key, e), sigma(resdir, ref, e)
            if not v:
                continue
            pts.append((e, v[0] / div, v[1] / div))
            if r:
                rat.append((e, v[0] / r[0], v[1] / r[0]))
        if not pts:
            continue
        # the figure heading already says the order (user, 2026-08-28); a row
        # at a DIFFERENT order keeps its own, which is what strip_order does
        ax_top.errorbar([p[0] for p in pts], [p[1] for p in pts],
                        yerr=[p[2] for p in pts], marker=marker, ms=5.5,
                        color=colour, lw=1.4, capsize=2.5,
                        label=_strip_order(label, order),
                        zorder=3)
        if rat:
            ax_rat.errorbar([p[0] for p in rat], [p[1] for p in rat],
                            yerr=[p[2] for p in rat], marker=marker, ms=5.5,
                            color=colour, lw=1.4, capsize=2.5, zorder=3)
    ax_rat.axhline(1.0, color="#111111", lw=1.0, zorder=1)
    for ax in (ax_top, ax_rat):
        ax.set_xscale("log")
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_xlim(min(es) * 0.8, max(es) * 1.25)
    ax_top.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
    ax_top.set_ylabel(r"$\sigma_{\rm fid}$" + tex(f" [{unit}]"))
    ax_top.tick_params(labelbottom=False)
    ax_rat.set_ylabel(tex("ratio to YADISM") + "\n" + tex("(same order)"),
                      fontsize=plotstyle.FS_YLABEL)
    ax_rat.set_xlabel(tex("beam energy [GeV]"))
    ax_rat.set_xticks(es)
    ax_rat.set_xticklabels([f"{e:g}" for e in es])
    ax_rat.xaxis.set_minor_formatter(ticker.NullFormatter())


def main():
    # TWO FIGURES, ONE PER ORDER (user, 2026-08-27).  A single panel carried
    # every LO and NLO generator plus three analytic curves on both currents,
    # and at that density the 1-2% differences the benchmark is about were
    # unreadable.  Splitting by order also makes each ratio panel honest:
    # every point in a figure is divided by the same reference.
    for order, ladder in (("LO", ""), ("NLO", " (the NNLO ladder in grey)")):
        fig, axes = plt.subplots(2, 2, figsize=(10.6, 6.6), sharex="col",
                                 gridspec_kw={"height_ratios": [2.4, 1],
                                              "hspace": 0.0, "wspace": 0.26})
        mu_dir, nu_dir = f"{BASE}/results", f"{BASE}/results_nu"
        panel(axes[0][0], axes[1][0], mu_dir,
              rows_at(MU_ROWS, order, mu_dir),
              "Muon DIS (NC): $\\mu^-$ on p at rest", "nb", 1e3, order)
        panel(axes[0][1], axes[1][1], nu_dir,
              rows_at(NU_ROWS, order, nu_dir),
              "Neutrino DIS (CC): $\\nu_\\mu$ on p at rest", "pb", 1.0,
              order)
        # LEGENDS GO BELOW THEIR PANEL, NOT INSIDE IT.  The curves run
        # diagonally from lower-left to upper-right, so every in-axes corner
        # is either on the data or on the errorbars.  Anchored to the RATIO
        # panel, not the top one: anchoring to the top axes means the offset
        # has to clear the ratio panel beneath it, which depends on the
        # height ratio -- guessing it left a large void.
        for col in (0, 1):
            _h, _l = axes[0][col].get_legend_handles_labels()
            axes[1][col].legend(_h, [tex(x) for x in _l],
                                fontsize=plotstyle.FS_LEGEND - 2,
                                frameon=True, ncol=1, loc="upper center",
                                bbox_to_anchor=(0.5, -0.42),
                                borderaxespad=0.0, handlelength=1.6,
                                columnspacing=1.0, labelspacing=0.35)
        # TWO LINES: on one, the title plus the cut list ran past both edges
        # of the figure at every font size that kept it legible.
        fig.suptitle(tex(f"Inclusive fiducial cross-section versus beam "
                         f"energy \u2014 {order}{ladder}") + "\n"
                     + tex("$Q^2 > 4$ GeV$^2$, $0.2 < y < 0.9$"),
                     # lifted clear of the per-panel titles: the second line
                     # was printed straight through "Muon DIS (NC): ..."
                     y=1.07, fontsize=plotstyle.FS_SUPTITLE)
        if order == "NLO":
            # BELOW the legends, which hang off the ratio panel at -0.42 in
            # axes coordinates; at -0.02 this note was printed straight
            # through them.
            fig.text(0.5, -0.30,
                     tex("Red: FONLL, which restores charm-mass effects in "
                         "the heavy-quark part of the structure function. "
                         "In the ratio panel it is drawn against the "
                         "same-order ZM-VFNS curve, so it reads directly as "
                         "the size of the scheme change."),
                     ha="center", va="top", fontsize=8.2, style="italic",
                     color="#5d6470")
        out = f"{BASE}/results/cmp_sigma_vs_E_{order.lower()}.png"
        fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out}")
    dimuon_figure()


if __name__ == "__main__":
    main()
