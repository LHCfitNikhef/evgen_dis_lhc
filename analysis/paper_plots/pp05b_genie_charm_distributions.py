#!/usr/bin/env python3
"""Paper plot 5b: GENIE's DIFFERENTIAL charm production at 1 TeV, the
counterpart of plot 5 in the way plot 7 is the counterpart of plot 6
(user, 2026-10-04: "why don't we have Figure 5.5 at the differential level,
in analogy with Figure 5.4? ... if only for the neutrino differential
distributions").

The same rows and the same reference as plot 5 -- the three GENIE
configurations, POWHEG-V2 (massless) and POWHEG-V2mc (massive charm), YADISM
in FFNS n_f = 3, ratio to YADISM (FONLL) -- in the three observables of
plot 7 (x_Bj, E_h, Q2), with plot 7's binning and x windows on the neutrino
column so the two differential charm figures can be read side by side.

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184.
  * THE CHARM REGION, Q2 > 4 GeV2 and W > 5 GeV, no y cut (selection q4w5).
  * 1 TeV, inside G18_02a's declared validity, so all three GENIE rows are
    drawn.
  * THE BAND IS FONLL, the scheme of the reference it surrounds
    (analysis/mhou_diff.py --current nu --scheme fonll --flavour charm with
    BENCH_SELECTION=q4w5 BENCH_TARGET=W BENCH_TMC=3, computed 2026-10-04 for
    this figure), and the references carry target-mass corrections.

>>> NEUTRINO ONLY, AND THE REASON IS A MEASURED ZERO (rule 1b), as on plot
5. <<<  GENIE's muon neutral-current DIS has no charm matrix element and AGKY
makes no charm in fragmentation, so there is no muon column to draw.

Usage: analysis/paper_plots/pp05b_genie_charm_distributions.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import plotstyle                                            # noqa: E402
plotstyle.apply()
from plotstyle import tex                                   # noqa: E402
import matplotlib.pyplot as plt                             # noqa: E402
import matplotlib.ticker as ticker                          # noqa: E402
import beams                                                # noqa: E402

SLUG = "pp05b_genie_charm_distributions"
IN_PAPER = True
TITLE = "GENIE differential charm production distributions at 1 TeV"
FIG_TITLE = (r"GENIE differential charm production distributions at "
             r"$E_{\nu} = 1$ TeV")
OUTPUT = "pp_genie_charm_distributions.png"
RESULTS = "results_nu"
CAPTION = ("Differential charm-production distributions in Bjorken <i>x</i>, "
           "hadronic energy and momentum transfer in neutrino "
           "charged-current scattering, per nucleon on tungsten at 1 TeV "
           "for Q&sup2; &gt; 4 GeV&sup2; and W &gt; 5 GeV, for the GENIE "
           "configurations, with POWHEG-V2 (massless charm), POWHEG-V2mc "
           "(massive charm) and YADISM in FFNS n_f = 3. Ratio to YADISM "
           "(FONLL) beneath, with its seven-point scale band; the references "
           "include target-mass corrections. Neutral-current charm is not "
           "shown because GENIE's electromagnetic DIS has no charm channel.")

SURFACE = "#ffffff"
# THE CANVAS IS 15 in WIDE against pp07's 11.4, and both are set at the text
# width, so every font is scaled by the ratio to print at pp07's size (user,
# 2026-10-04: title, then axis labels and ticks, "consistent with all other
# plots in the paper").
FSCALE = 15.0 / 11.4
# THE CANVAS IS 15 in WIDE against pp07's 11.4, and both are set at the text
# width, so every font is scaled by the ratio to print at pp07's size (user,
# 2026-10-04: title, then axis labels and ticks, "consistent with all other
# plots in the paper").
FSCALE = 15.0 / 11.4
REGION = "q4w5"
TSUF = "_W"
TMC = "_tmc"
ENERGY = beams.ANCHOR_ENERGY
Q2MIN = 4.0
REF = "yadism_charm_nlo_fonll_damp"
REF_LABEL = "YADISM (FONLL)"
REF_COLOUR = "#111111"
FFNS3 = "yadism_charm_nlo_ffns3"
FFNS3_LABEL = r"YADISM (FFNS $n_f=3$)"
# grey, not V2mc's brown (user, 2026-10-04): the two were indistinguishable
FFNS3_COLOUR = "#7f7f7f"
MHOU_FILE = "mhou_diff_fonll_charm.json"
# The xiF = 1/2 points are KEPT, as on plot 7.
MHOU_KEEP_HALF = True

# Plot 5's rows, colours and styles, so a reader carries one legend from the
# integrated figure to this one.
ROWS = [
    ("GENIE (GRV98LO)",  "genie_lo_charmfinal",     "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", "genie_nnpdf_charmfinal",  "#56b4e9", "D", "--"),
    ("GENIE (HEDIS)",    "genie_charmfinal",        "#e69f00", "^", "-."),
    ("POWHEG-V2",        "powheg_nu_charmfinal",    "#0072b2", "s", ":"),
    ("POWHEG-V2mc",      "powheg_nu_mc_charmfinal", "#8c6d1f", "v",
     (0, (6.5, 1.6))),
]

# Plot 7's observables, windows and merges on its neutrino column.
# (key, x label, x scale, y scale, y quantity, unit, divisor): E_h in fb so
# its linear axis does not read 0.00100.
OBS = [
    ("xbj", r"$x_{\rm Bj}$", "log", "log", r"d\sigma_c/dx_{\rm Bj}", "pb", 1.0),
    ("nu",  r"$E_h$  [GeV]", "linear", "linear", r"d\sigma_c/dE_h",
     "fb/GeV", 1e-3),
    ("Q2",  r"$Q^2$  [GeV$^2$]", "log", "log", r"d\sigma_c/dQ^2",
     "pb/GeV$^2$", 1.0),
]
XLIM = {"xbj": (3e-3, 0.3), "Q2": (None, 400.0)}
REBIN = {"nu": 5, "Q2": 2}


def _suffix(key):
    tag = beams.Beams("mu", ENERGY).tag
    esuf = "" if ENERGY == beams.ANCHOR_ENERGY else f"_{tag}"
    return f"_{REGION}{TSUF}{TMC if key.startswith('yadism') else ''}{esuf}"


def _merge(edges, dsig, err, n):
    """Merge n adjacent bins, keeping dsig a density (errors weighted by
    bin width, since what adds is the integral over each bin)."""
    m = (len(dsig) // n) * n
    w = np.diff(edges)[:m].reshape(-1, n)
    y = (dsig[:m].reshape(-1, n) * w).sum(axis=1)
    e = np.sqrt(((err[:m].reshape(-1, n) * w) ** 2).sum(axis=1))
    wide = w.sum(axis=1)
    return edges[:m + 1:n], y / wide, e / wide


def hist(key, obs):
    """(edges, dsig, err) for one observable, or None if absent; the file's
    region and Q2 stamp are checked."""
    p = f"{BASE}/{RESULTS}/histos_{key}{_suffix(key)}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    sel = d.get("selection") or {}
    if sel.get("name") not in (None, REGION):
        raise SystemExit(f"{os.path.relpath(p, BASE)} is not a {REGION} result")
    q2 = sel.get("q2_min")
    if q2 is not None and abs(q2 - Q2MIN) > 1e-9:
        raise SystemExit(f"{os.path.relpath(p, BASE)} was computed at "
                         f"Q2 > {q2:g}, but this figure is Q2 > {Q2MIN:g}")
    if d.get("stat_insufficient") or obs not in d.get("hists", {}):
        return None
    h = d["hists"][obs]
    out = (np.asarray(h["edges"]), np.asarray(h["dsig"]),
           np.asarray(h["err"]))
    n = REBIN.get(obs)
    return _merge(*out, n) if n else out


def mhou(obs):
    """(edges, rel_hi, rel_lo) of the 7-point FONLL band, merged on the
    cross-sections before the ratio is taken, or None."""
    p = (f"{BASE}/{RESULTS}/"
         + MHOU_FILE.replace(".json", f"_{REGION}{TSUF}{TMC}.json"))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        h = json.load(f)["hists"].get(obs)
    if h is None:
        return None
    hi_key = "rel_hi" if MHOU_KEEP_HALF else "rel_hi_nohalf"
    lo_key = "rel_lo" if MHOU_KEEP_HALF else "rel_lo_nohalf"
    e0 = np.asarray(h["edges"])
    cen = np.asarray(h["dsig"])
    hi = cen * (1.0 + np.asarray(h[hi_key]))
    lo = cen * (1.0 + np.asarray(h[lo_key]))
    e = e0
    n = REBIN.get(obs)
    if n:
        z = np.zeros_like(cen)
        e, cen, _ = _merge(e0, cen, z, n)
        _, hi, _ = _merge(e0, hi, z, n)
        _, lo, _ = _merge(e0, lo, z, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (e, np.where(cen > 0, hi / cen - 1.0, 0.0),
                np.where(cen > 0, lo / cen - 1.0, 0.0))


def _snap(edges, v, upper):
    e = np.asarray(edges, dtype=float)
    if v is None:
        return float(e[-1] if upper else e[0])
    return float(e[int(np.argmin(np.abs(e - v)))])


def _x_kin_min():
    """The smallest x_Bj the region reaches, Q2min/(y_max 2k.P), as on plot 7."""
    return Q2MIN / (0.9 * beams.Beams("mu", ENERGY).two_kP)


def x_lo(obs, edges, y):
    v = XLIM.get(obs, (None, None))[0]
    if obs == "xbj":
        e = np.asarray(edges, dtype=float)
        want = max(v if v is not None else 0.0, _x_kin_min())
        above = e[e >= want - 1e-12]
        return float(above[0]) if above.size else float(e[0])
    if v is None:
        nz = np.nonzero(np.asarray(y) > 0)[0]
        if nz.size:
            return float(np.asarray(edges)[nz[0]])
    return _snap(edges, v, False)


def x_hi(obs, edges):
    return _snap(edges, XLIM.get(obs, (None, None))[1], True)


def _bulk(edges, y):
    """Bins carrying at least 0.1% of the distribution (drawn either way,
    but only these set the ratio range and enter the claims)."""
    w = y * np.diff(edges)
    return w > 0.001 * w.sum()


def _vis(obs):
    """Mask of the bins the figure shows -- one definition for the drawing
    and for the claims."""
    edges, ry, _ = hist(REF, obs)
    return (_bulk(edges, ry) & (ry > 0)
            & (edges[:-1] >= x_lo(obs, edges, ry))
            & (edges[1:] <= x_hi(obs, edges)))


def _log_ticks(a, lo, hi):
    """1-2-5 ticks on a log axis spanning less than two decades."""
    if hi / lo >= 100.0:
        return
    ticks = [m * 10.0 ** k for k in range(-6, 7) for m in (1, 2, 5)
             if lo <= m * 10.0 ** k <= hi]
    a.set_xticks(ticks)
    a.set_xticklabels([f"{t:g}" for t in ticks])
    a.xaxis.set_minor_formatter(ticker.NullFormatter())


def block(ax, axr, obs, xs, ys, ylab_q, unit, div):
    """One observable: absolute panel plus its ratio panel."""
    edges, ry, _ = hist(REF, obs)
    mids = 0.5 * (edges[:-1] + edges[1:])
    mask = _vis(obs)
    seen = []
    ax.stairs(ry / div, edges, color=REF_COLOUR, lw=1.7, ls=":", baseline=None,
              zorder=7, label=tex(REF_LABEL))
    mb = mhou(obs)
    if mb is not None and np.array_equal(mb[0], edges):
        axr.stairs(1.0 + mb[1], edges, baseline=1.0 + mb[2], fill=True,
                   color=REF_COLOUR, alpha=0.12, lw=0, zorder=0,
                   label=tex("NLO MHOU"))
        seen += [(1.0 + mb[1], mask), (1.0 + mb[2], mask)]
        # and in the absolute panel, around the reference (user,
        # 2026-10-04, as on pp07); the lower edge floored for the log axes
        ax.stairs(ry * (1.0 + mb[1]) / div, edges,
                  baseline=np.maximum(ry * (1.0 + mb[2]), 1e-12 * ry.max()) / div,
                  fill=True, color=REF_COLOUR, alpha=0.12, lw=0, zorder=0)
    ff = hist(FFNS3, obs)
    if ff is not None:
        ax.stairs(ff[1] / div, edges, color=FFNS3_COLOUR, lw=1.7, ls="-",
                  baseline=None, zorder=6, label=tex(FFNS3_LABEL))
        with np.errstate(divide="ignore", invalid="ignore"):
            fr = np.where(ry > 0, ff[1] / ry, np.nan)
        axr.stairs(fr, edges, color=FFNS3_COLOUR, lw=1.7, ls="-",
                   baseline=None, zorder=6)
        seen.append((fr, mask))
    for lab, key, colour, _marker, ls in ROWS:
        h = hist(key, obs)
        if h is None:
            continue
        _, y, e = h
        ax.stairs(y / div, edges, color=colour, lw=1.8, ls=ls, baseline=None,
                  zorder=5, label=tex(lab))
        ax.errorbar(mids, y / div, yerr=e / div, fmt="none", ecolor=colour,
                    elinewidth=1.0, alpha=0.75, zorder=5)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(ry > 0, y / ry, np.nan)
            re = np.where(ry > 0, e / ry, np.nan)
        axr.stairs(r, edges, color=colour, lw=1.8, ls=ls, baseline=None,
                   zorder=5)
        axr.errorbar(mids, r, yerr=re, fmt="none", ecolor=colour,
                     elinewidth=1.0, alpha=0.75, zorder=5)
        seen.append((r, mask))
    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-", zorder=1)
    x0, x1 = x_lo(obs, edges, ry), x_hi(obs, edges)
    for a in (ax, axr):
        a.set_xscale(xs)
        a.set_xlim(x0, x1)
        if xs == "log":
            _log_ticks(a, x0, x1)
        a.grid(alpha=0.22, lw=0.6)
    ax.set_yscale(ys)
    if ys == "log":
        pos = ry[mask & (ry > 0)] / div
        if pos.size:
            ax.set_ylim(0.3 * float(pos.min()), 3.0 * float(pos.max()))
    else:
        ax.set_ylim(0.0, None)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="lower"))
    axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="upper"))
    ax.set_ylabel(tex(rf"${ylab_q}$  [{unit}]"), fontsize=round((plotstyle.FS_YLABEL - 2) * FSCALE))
    ax.tick_params(labelbottom=False)
    return seen


def ratio_range(seen):
    """One ratio range for the whole row: the GRV98LO row sits at ~0.6 and
    the others near 1, and a shared axis lets the panels be compared."""
    lo, hi = 1.0, 1.0
    for r, m in seen:
        v = r[m & np.isfinite(r)]
        if v.size:
            lo, hi = min(lo, float(v.min())), max(hi, float(v.max()))
    pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


def main():
    fig = plt.figure(figsize=(15.0, 6.2))
    outer = fig.add_gridspec(1, len(OBS), wspace=0.32, top=0.765,
                             bottom=0.12, left=0.06, right=0.99)
    rat, seen_all = [], []
    for j, (obs, xlab, xs, ys, ylab_q, unit, div) in enumerate(OBS):
        inner = outer[0, j].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                        hspace=0.0)
        ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        for a_ in (ax, axr):
            a_.tick_params(which="both",
                           labelsize=round(plotstyle.FS_TICKS * FSCALE))
        seen_all += block(ax, axr, obs, xs, ys, ylab_q, unit, div)
        axr.set_xlabel(tex(xlab), fontsize=round((plotstyle.FS_XLABEL - 1) * FSCALE))
        if j == 0:
            axr.set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                           fontsize=round((plotstyle.FS_YLABEL - 4) * FSCALE))
        rat.append(axr)
    lo, hi = ratio_range(seen_all)
    for a in rat:
        a.set_ylim(lo, hi)
    h, lab = [], []
    for a in fig.axes[:2]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    # YADISM first, then POWHEG, then GENIE -- plot 5's order (user,
    # 2026-10-04).
    want = [tex(REF_LABEL), tex(FFNS3_LABEL), tex("NLO MHOU")] \
        + [tex(r[0]) for r in ROWS[3:]] + [tex(r[0]) for r in ROWS[:3]]
    order = [lab.index(w) for w in want if w in lab]
    order += [i for i in range(len(lab)) if i not in order]
    h, lab = [h[i] for i in order], [lab[i] for i in order]
    # the band shares the reference's entry (user, 2026-10-04)
    h, lab = plotstyle.merge_bands(h, lab,
                                   {tex("NLO MHOU"): [tex(REF_LABEL)]})
    # LEGEND FONT RAISED AGAIN (user, 2026-10-04: "quite a bit")
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.525, 0.995),
               ncol=4, frameon=True, handlelength=2.1, columnspacing=1.15,
               labelspacing=0.35, fontsize=plotstyle.FS_LEGEND + 8)
    # THE TITLE SCALES WITH THE CANVAS (user, 2026-10-04): this figure is
    # 15 in wide against ~11.4 for the others, all set at the text width, so
    # FS_SUPTITLE alone printed ~30% smaller than every other title.
    fig.suptitle(tex(FIG_TITLE), y=1.08,
                 fontsize=round(plotstyle.FS_SUPTITLE * FSCALE))
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------- claims --
def _ratio(key, obs):
    """key / FONLL over the visible bins."""
    m = _vis(obs)
    return hist(key, obs)[1][m] / hist(REF, obs)[1][m]


def _span(key, obs):
    r = _ratio(key, obs)
    return float(r.min()), float(r.max())


def _all(key):
    """(min, max) of key / FONLL over every visible bin of the three panels."""
    v = np.concatenate([_ratio(key, o) for o, *_ in OBS])
    return float(v.min()), float(v.max())


def _first_last(key, obs):
    r = _ratio(key, obs)
    return float(r[0]), float(r[-1])


def _rel(a, b, obs):
    """a / b bin by bin, (min, max), over the visible bins."""
    r = _ratio(a, obs) / _ratio(b, obs)
    return float(r.min()), float(r.max())


MESSAGE = (
    "<b>GENIE's charm deficit is a shape, not only a normalisation.</b> "
    "The default tune on GRV98 LO, 0.63 of the FONLL reference integrated "
    "(plot 5), is 0.43&ndash;0.72 of it bin by bin, rising from the lowest "
    "to the highest bin in all three variables: from 0.48 to 0.70 in "
    "x<sub>Bj</sub>, from 0.55 to 0.72 in E<sub>h</sub> and from 0.43 to "
    "0.72 in Q&sup2;. The deficit is worst at small Q&sup2; and small x, "
    "where the strange sea that charged-current charm is made of is "
    "smallest in GRV98 LO relative to NNPDF4.0."
    "\n\n"
    "<b>Substituting NNPDF4.0 removes most of it.</b> The same tune on "
    "NNPDF4.0 lies at 0.88&ndash;1.13 of the reference, the substitution "
    "multiplying the charm cross-section by 1.33&ndash;2.12 bin by bin in "
    "x<sub>Bj</sub> and Q&sup2; and by 1.53&ndash;1.62 in E<sub>h</sub>; a "
    "residual tilt in E<sub>h</sub>, from 0.88 to 1.13, is left by the "
    "tune's own charm model. <b>GENIE HEDIS</b> sits at 0.69&ndash;0.91, "
    "lowest again at small Q&sup2; and x."
    "\n\n"
    "<b>POWHEG-V2mc reproduces the massive calculation bin by bin</b>: "
    "within 3% of YADISM in FFNS n<sub>f</sub> = 3 in every visible bin of "
    "the three distributions, while massless POWHEG-V2 rises above FONLL "
    "towards small x, where the mass effect is largest."
)

CLAIMS = [
    {"what": "the reference is FONLL with target-mass corrections, in the "
             "charm region q4w5 on tungsten at 1 TeV, with a FONLL band",
     "check": lambda: REF.endswith("fonll_damp") and REGION == "q4w5"
     and TSUF == "_W" and TMC == "_tmc" and ENERGY == 1000.0
     and "fonll" in MHOU_FILE and all(mhou(o) is not None for o, *_ in OBS),
     "detail": lambda: f"{REF}_{REGION}{TSUF}{TMC}, {MHOU_FILE}"},
    {"what": "GENIE (GRV98LO) is 0.43-0.72 of FONLL bin by bin",
     "check": lambda: 0.40 < _all("genie_lo_charmfinal")[0] < 0.46
     and 0.715 < _all("genie_lo_charmfinal")[1] < 0.725,
     "detail": lambda: "%.3f-%.3f" % _all("genie_lo_charmfinal")},
    {"what": "GENIE (GRV98LO) rises from first to last bin: x 0.48->0.70, "
             "E_h 0.55->0.72, Q2 0.43->0.72",
     "check": lambda: all(abs(a - x) < 0.015 and abs(b - y) < 0.015
                          for (a, b), (x, y) in zip(
                              [_first_last("genie_lo_charmfinal", o)
                               for o, *_ in OBS],
                              [(0.48, 0.70), (0.55, 0.72), (0.43, 0.72)])),
     "detail": lambda: ", ".join("%s %.3f->%.3f" % ((o,) + _first_last(
         "genie_lo_charmfinal", o)) for o, *_ in OBS)},
    {"what": "GENIE (NNPDF4.0) is 0.88-1.13 of FONLL; the substitution "
             "multiplies charm by 1.33-2.12 in x and Q2 and 1.53-1.62 in E_h",
     "check": lambda: 0.87 < _all("genie_nnpdf_charmfinal")[0] < 0.895
     and 1.115 < _all("genie_nnpdf_charmfinal")[1] < 1.14
     and all(1.30 < _rel("genie_nnpdf_charmfinal", "genie_lo_charmfinal", o)[0]
             and _rel("genie_nnpdf_charmfinal", "genie_lo_charmfinal", o)[1]
             < 2.15 for o in ("xbj", "Q2"))
     and 1.52 < _rel("genie_nnpdf_charmfinal", "genie_lo_charmfinal", "nu")[0]
     and _rel("genie_nnpdf_charmfinal", "genie_lo_charmfinal", "nu")[1] < 1.63,
     "detail": lambda: "%.3f-%.3f; " % _all("genie_nnpdf_charmfinal")
     + ", ".join("%s x%.2f-%.2f" % ((o,) + _rel(
         "genie_nnpdf_charmfinal", "genie_lo_charmfinal", o)) for o, *_ in OBS)},
    {"what": "GENIE (NNPDF4.0) tilts from 0.88 to 1.13 in E_h",
     "check": lambda: abs(_first_last("genie_nnpdf_charmfinal", "nu")[0]
                          - 0.88) < 0.015
     and abs(_first_last("genie_nnpdf_charmfinal", "nu")[1] - 1.13) < 0.015,
     "detail": lambda: "%.3f->%.3f" % _first_last("genie_nnpdf_charmfinal",
                                                   "nu")},
    {"what": "GENIE (HEDIS) is 0.69-0.91 of FONLL",
     "check": lambda: 0.68 < _all("genie_charmfinal")[0] < 0.70
     and 0.90 < _all("genie_charmfinal")[1] < 0.92,
     "detail": lambda: "%.3f-%.3f" % _all("genie_charmfinal")},
    {"what": "POWHEG-V2mc is within 3% of YADISM FFNS n_f=3 in every visible "
             "bin; POWHEG-V2 rises above FONLL towards small x",
     "check": lambda: all(max(abs(_rel("powheg_nu_mc_charmfinal", FFNS3, o)[0]
                                  - 1), abs(_rel("powheg_nu_mc_charmfinal",
                                                 FFNS3, o)[1] - 1)) < 0.031
                          for o, *_ in OBS)
     and _first_last("powheg_nu_charmfinal", "xbj")[0]
     > _first_last("powheg_nu_charmfinal", "xbj")[1] > 1.0,
     "detail": lambda: ", ".join("%s %.3f-%.3f" % ((o,) + _rel(
         "powheg_nu_mc_charmfinal", FFNS3, o)) for o, *_ in OBS)},
]


if __name__ == "__main__":
    main()
