#!/usr/bin/env python3
"""Paper plot 4: GENIE's differential DIS distributions at 1 TeV, in the
FINAL benchmark region.

The redo of the earlier figure for the "Paper plots" tab (user, 2026-09-14:
"produce now [the new version] of pp_03 and pp_04,
if statistics are an issue, feel free to increase sample size (for 1 TeV of
course)").  It stands to pp03 as pp02 stands to pp01 : the structure
(four observables, binning, windows, fused panel pairs, file checks) is pp02
its, and WHICH CURVES are drawn is earlier pp04's.

  * Tungsten per nucleon, (74 sigma_p + 110 sigma_n)/184, every generator run
    on the proton and the neutron separately; Q2 > 4 GeV2 and W > 3 GeV with
    no y cut on BOTH currents; 1 TeV; references with target-mass corrections.
  * THE REFERENCE IS FONLL, as in the earlier production pp04 and pp03: GENIE carries charm masses.
    YADISM (ZM) is drawn beside it, so the charm-mass effect stays visible and
    the massless POWHEG entries have their own reference on the page.
  * THE ROWS (CONVENTIONS.md rule 1b: both GENIE rows count): G18_02a on GRV98 (the
    FASER tune), G18_02a with NNPDF4.0 (the PDF diagnostic), HEDIS GHE19_00a
    (neutrino only), and POWHEG-RES / POWHEG-V2 for scale.  The GENIE samples
    are genie/production/genie_job.sh: the earlier production's configurations with the target as the only
    change, 500k (G18_02a) and 400k (HEDIS) events per nucleon at 1 TeV.
  * NO NNLO CURVE, as in the earlier production pp04: the GENIE deviations are several times the
    NNLO correction.  The NLO band is kept, in FONLL, bin by bin, keeping the
    xiF = 1/2 points as pp02 does.

Usage: analysis/paper_plots/pp04_genie_dis_distributions.py
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
import target                                               # noqa: E402

SLUG = "pp04_genie_dis_distributions"
V1_NUMBER = 4
# IN THE PAPER in place of the earlier production's figure (user, 2026-09-14: every updated 
# plot goes into the paper).  Read by tools/check_paper_figures_used.py.
IN_PAPER = True
TITLE = "GENIE differential DIS distributions at 1 TeV"
# the title ON the figure (user, 2026-09-18); TITLE above names the report tab
FIG_TITLE = r"Differential DIS distributions at $E_{\ell} = 1$ TeV"
OUTPUT = "pp_genie_dis_distributions.png"
RESULTS = "results_nu"
CAPTION = ("Differential distributions per nucleon on tungsten at 1 TeV, for "
           "Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no y cut, in "
           "Bjorken <i>x</i>, scattered-lepton energy and angle, and momentum "
           "transfer, for the GENIE configurations with POWHEG-RES and "
           "POWHEG-V2 for scale; muon neutral current (left) and neutrino "
           "charged current (right). Ratio to YADISM (FONLL) with target-mass "
           "corrections beneath each, with YADISM (ZM) drawn as well. The "
           "shaded band is the reference's seven-point scale variation, bin "
           "by bin.")

SURFACE = "#ffffff"
ENERGY = 1000.0
REGION = "q4w3"
Q2MIN = 4.0
W_MIN = 3.0
TAG = f"_{REGION}_W"
TMC = "_tmc"

REF = f"yadism_nlo_fonll_damp{TAG}{TMC}"
REF_LABEL = "YADISM NLO (FONLL)"     # the order is stated (user, 2026-09-18)
ZM_REF = f"yadism_nlo{TAG}{TMC}"
ZM_LABEL = "YADISM NLO (ZM)"
REF_COLOUR = "#111111"
ZM_COLOUR = "#c0392b"

# The 7-point band on the reference, bin by bin, IN FONLL to match the curve
# it surrounds (BENCH_SELECTION=q4w3 BENCH_TARGET=W BENCH_TMC=3
# analysis/mhou_diff.py --current mu|nu --scheme fonll).
MHOU_FILE = f"mhou_diff_fonll{TAG}{TMC}.json"
MHOU_KEEP_HALF = True

# earlier pp04's rows and styles; the two POWHEGs and HEDIS name their current, as
# there is one legend above both columns.
POWHEG_RES_LS = (0, (6.5, 1.6))
POWHEG_V2_LS = (0, (6.5, 1.4, 1.2, 1.4))
MU_ROWS = [
    ("GENIE (GRV98LO)",          f"genie{TAG}",        "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", f"genie_nnpdf{TAG}",  "#56b4e9", "D", "--"),
    (r"POWHEG-RES ($\mu$)",      f"powheg{TAG}",       "#0072b2", "s",
     POWHEG_RES_LS),
]
NU_ROWS = [
    ("GENIE (GRV98LO)",          f"genie_lo{TAG}",     "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", f"genie_nnpdf{TAG}",  "#56b4e9", "D", "--"),
    ("GENIE (HEDIS)", f"genie{TAG}",  "#e69f00", "^", "-."),
    (r"POWHEG-V2 ($\nu$)",       f"powheg_nu{TAG}",    "#0072b2", "s",
     POWHEG_V2_LS),
]

# G18_02a is declared valid to 1000 GeV (user, 2026-09-01); this figure is AT
# 1 TeV, so nothing is cut, but the check is kept so moving ENERGY cannot draw
# a clamped spline.  `genie` on the neutrino side is HEDIS (valid to 1e12).
EMAX = {(f"genie{TAG}", "results"): 1000.0,
        (f"genie_nnpdf{TAG}", "results"): 1000.0,
        (f"genie_lo{TAG}", "results_nu"): 1000.0,
        (f"genie_nnpdf{TAG}", "results_nu"): 1000.0}

# pp02's windows (the user's earlier values), merges and log muon E_lep' panel,
# so the two figures overlay bin for bin.
XLIM = {
    "theta": {"results": (None, 3e-2), "results_nu": (None, 8e-2)},
    "Q2": {"results": (None, 100.0), "results_nu": (None, 400.0)},
    "xbj": {"results": (None, 0.6), "results_nu": (None, 0.6)},
}
REBIN = {("Emu", "results_nu"): 5, ("Emu", "results"): 5,
         ("Q2", "results_nu"): 2,
         ("xbj", "results"): 2, ("xbj", "results_nu"): 2}
YSCALE = {("Emu", "results"): "log"}

OBS = [
    ("xbj",   r"$x_{\rm Bj}$", "log", "log", r"d\sigma/dx_{\rm Bj}"),
    ("Emu",   r"$E_{\ell'}$  [GeV]", "linear", "linear",
     r"d\sigma/dE_{\ell'}"),
    ("theta", r"$\theta_{\ell'}$  [rad]", "log", "log",
     r"d\sigma/d\theta_{\ell'}"),
    ("Q2",    r"$Q^2$  [GeV$^2$]", "log", "log",
     r"d\sigma/dQ^2"),
]
# titles are the current alone (user, 2026-09-14)
COLUMNS = [("results", MU_ROWS, "Muon DIS", 1000.0, "nb"),
           ("results_nu", NU_ROWS, "Neutrino DIS", 1.0, "pb")]


def _snap(edges, v, upper):
    e = np.asarray(edges, dtype=float)
    if v is None:
        return float(e[-1] if upper else e[0])
    return float(e[int(np.argmin(np.abs(e - v)))])


def _merge(edges, dsig, err, n):
    """Merge n adjacent bins, keeping dsig a density (errors width-weighted)."""
    m = (len(dsig) // n) * n
    w = np.diff(edges)[:m].reshape(-1, n)
    y = (dsig[:m].reshape(-1, n) * w).sum(axis=1)
    e = np.sqrt(((err[:m].reshape(-1, n) * w) ** 2).sum(axis=1))
    wide = w.sum(axis=1)
    return edges[:m + 1:n], y / wide, e / wide


def _load(resdir, key):
    """The result JSON, CHECKED to be what the figure claims, or None.

    The q4w3 region at 1 TeV, a tungsten per-nucleon result, and for a
    reference the TMC mode its name says and the W184free set.  A earlier proton
    file, a TMC-off reference or another energy would each draw a plausible
    curve and say nothing.
    """
    p = f"{BASE}/{resdir}/histos_{key}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    rel = os.path.relpath(p, BASE)
    sel = d.get("selection") or {}
    if (sel.get("name") != REGION or abs(sel.get("q2_min", -1) - Q2MIN) > 1e-9
            or abs(sel.get("w2_min", -1) - W_MIN ** 2) > 1e-9
            or sel.get("y_min") != 0.0 or sel.get("y_max") != 1.0):
        raise SystemExit(f"{rel}: not the {REGION} region ({sel})")
    if d.get("target") != "W":
        raise SystemExit(f"{rel}: target {d.get('target')!r}, not tungsten")
    if abs(float(d.get("energy_gev", -1)) - ENERGY) > 1e-6:
        raise SystemExit(f"{rel}: energy {d.get('energy_gev')}, not {ENERGY:g}")
    if key.startswith("yadism"):
        want = target.TMC_CHOSEN if key.endswith(TMC) else 0
        if d.get("tmc") != want or d.get("pdf_set") != target.PDFSET_W:
            raise SystemExit(f"{rel}: tmc={d.get('tmc')} pdf="
                             f"{d.get('pdf_set')}, expected tmc={want} "
                             f"pdf={target.PDFSET_W}")
    return d


def hist(resdir, key, obs, rebin=True):
    """(edges, dsig, err) for one observable, or None if absent (or outside
    the tune's declared validity)."""
    lim = EMAX.get((key, resdir))
    if lim is not None and ENERGY > lim + 1e-9:
        return None
    d = _load(resdir, key)
    if d is None or d.get("stat_insufficient") or obs not in d.get("hists", {}):
        return None
    h = d["hists"][obs]
    out = (np.asarray(h["edges"], dtype=float), np.asarray(h["dsig"]),
           np.asarray(h.get("err", np.zeros(len(h["dsig"])))))
    n = REBIN.get((obs, resdir)) if rebin else None
    return _merge(*out, n) if n else out


def raw_hist(resdir, key, obs):
    return hist(resdir, key, obs, rebin=False)


def mhou(resdir, obs, fn=None):
    """(edges, rel_hi, rel_lo) of the 7-point band, or None.

    Merged on the cross-sections, the ratio taken afterwards.  The file's own
    stamps are checked like a histogram's.
    """
    fn = fn or MHOU_FILE
    p = f"{BASE}/{resdir}/{fn}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if (d.get("target") != "W" or d.get("tmc") != target.TMC_CHOSEN
            or (d.get("selection") or {}).get("name") != REGION
            or abs(float(d.get("energy_gev", -1)) - ENERGY) > 1e-6
            or d.get("pdf") != target.PDFSET_W):
        raise SystemExit(f"{resdir}/{fn}: not a 1 TeV tungsten TMC-on "
                         f"{REGION} band")
    h = d["hists"].get(obs)
    if h is None:
        return None
    hi_key = "rel_hi" if MHOU_KEEP_HALF else "rel_hi_nohalf"
    lo_key = "rel_lo" if MHOU_KEEP_HALF else "rel_lo_nohalf"
    e = np.asarray(h["edges"])
    cen = np.asarray(h["dsig"])
    hi = cen * (1.0 + np.asarray(h[hi_key]))
    lo = cen * (1.0 + np.asarray(h[lo_key]))
    n = REBIN.get((obs, resdir))
    if n:
        z = np.zeros_like(cen)
        e0 = np.asarray(h["edges"])
        e, cen, _ = _merge(e0, cen, z, n)
        _, hi, _ = _merge(e0, hi, z, n)
        _, lo, _ = _merge(e0, lo, z, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (e, np.where(cen > 0, hi / cen - 1.0, 0.0),
                np.where(cen > 0, lo / cen - 1.0, 0.0))


def _bulk(edges, y):
    """Bins carrying at least 0.1% of the distribution (they set the range)."""
    w = y * np.diff(edges)
    return w > 0.001 * w.sum()


def ratio_range(rows):
    dev = 0.0
    for r, m in rows:
        v = r[m & np.isfinite(r)] if m is not None else r[np.isfinite(r)]
        if v.size:
            dev = max(dev, float(np.max(np.abs(v - 1.0))))
    dev = max(dev, 0.01)
    return 1.0 - 1.25 * dev, 1.0 + 1.25 * dev


def _log_ticks(a, lo, hi):
    """1-2-5 ticks on a log axis spanning less than two decades."""
    if hi / lo >= 100.0:
        return
    ticks = [m * 10.0 ** k for k in range(-6, 7) for m in (1, 2, 5)
             if lo <= m * 10.0 ** k <= hi]
    a.set_xticks(ticks)
    a.set_xticklabels([f"{t:g}" for t in ticks])
    a.xaxis.set_minor_formatter(ticker.NullFormatter())


def _x_kin_min():
    """Smallest reachable x_Bj, Q2min/(y_max 2k.P), with y_max = 1 here."""
    return Q2MIN / beams.Beams("mu", ENERGY).two_kP


def x_lo(obs, resdir, edges, y=None):
    """The window's lower edge, defaulting to the first POPULATED bin."""
    v = XLIM.get(obs, {}).get(resdir, (None, None))[0]
    if obs == "xbj":
        e = np.asarray(edges, dtype=float)
        want = max(v if v is not None else 0.0, _x_kin_min())
        above = e[e >= want - 1e-12]
        return float(above[0]) if above.size else float(e[0])
    if v is None and y is not None:
        # the first bin carrying the distribution (_bulk), not merely a
        # non-zero one: the reference leaks a sliver below theta = 2 mrad, the
        # smallest angle Q2 > 4 allows, which opened the panel on empty axis
        nz = np.nonzero(_bulk(np.asarray(edges), np.asarray(y)))[0]
        if nz.size:
            return float(np.asarray(edges)[nz[0]])
    return _snap(edges, v, False)


def x_hi(obs, resdir, edges):
    return _snap(edges, XLIM.get(obs, {}).get(resdir, (None, None))[1], True)


def _vis(resdir, obs):
    """Mask of the bins this figure SHOWS, for the drawing and the claims."""
    edges, ry, _ = hist(resdir, REF, obs)
    return (_bulk(edges, ry) & (ry > 0)
            & (edges[:-1] >= x_lo(obs, resdir, edges, ry) - 1e-12)
            & (edges[1:] <= x_hi(obs, resdir, edges) + 1e-12))


def _same_edges(a, b):
    return a.shape == b.shape and np.allclose(a, b, rtol=1e-12, atol=0.0)


def block(ax, axr, resdir, rows, obs, div, unit, xs, ys, ylab_q):
    """One observable in one current: absolute panel plus its ratio panel."""
    ref = hist(resdir, REF, obs)
    if ref is None:
        return
    edges, ry, _ = ref
    mids = 0.5 * (edges[:-1] + edges[1:])
    mask = _vis(resdir, obs)
    seen = []

    def edges_of(h, what):
        if not _same_edges(h[0], edges):
            raise SystemExit(f"{resdir} {obs}: {what} is binned differently "
                             f"from {REF}")

    ax.stairs(ry / div, edges, color=REF_COLOUR, lw=1.7, ls=":",
              baseline=None, zorder=7, label=tex(REF_LABEL))
    mb = mhou(resdir, obs)
    if mb is not None:
        edges_of(mb, MHOU_FILE)
        axr.stairs(1.0 + mb[1], edges, baseline=1.0 + mb[2], fill=True,
                   color=REF_COLOUR, alpha=0.12, lw=0, zorder=0,
                   label=tex("NLO MHOU"))
    zm = hist(resdir, ZM_REF, obs)
    if zm is not None:
        edges_of(zm, ZM_REF)
        ax.stairs(zm[1] / div, edges, color=ZM_COLOUR, lw=1.7, ls="-",
                  baseline=None, zorder=6, label=tex(ZM_LABEL))
        with np.errstate(divide="ignore", invalid="ignore"):
            zr = np.where(ry > 0, zm[1] / ry, np.nan)
        axr.stairs(zr, edges, color=ZM_COLOUR, lw=1.7, ls="-",
                   baseline=None, zorder=6)
        seen.append((zr, mask))

    for lab, key, colour, _marker, ls in rows:
        h = hist(resdir, key, obs)
        if h is None:
            continue
        edges_of(h, key)
        _, y, e = h
        ax.stairs(y / div, edges, color=colour, lw=1.8, ls=ls,
                  baseline=None, zorder=5, label=tex(lab))
        ax.errorbar(mids, y / div, yerr=e / div, fmt="none", ecolor=colour,
                    elinewidth=1.0, alpha=0.75, zorder=5)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(ry > 0, y / ry, np.nan)
            re = np.where(ry > 0, e / ry, np.nan)
        axr.stairs(r, edges, color=colour, lw=1.8, ls=ls,
                   baseline=None, zorder=5)
        axr.errorbar(mids, r, yerr=re, fmt="none", ecolor=colour,
                     elinewidth=1.0, alpha=0.75, zorder=5)
        seen.append((r, mask))

    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-", zorder=1)
    x0, x1 = x_lo(obs, resdir, edges, ry), x_hi(obs, resdir, edges)
    for a in (ax, axr):
        a.set_xscale(xs)
        a.set_xlim(x0, x1)
        if xs == "log":
            _log_ticks(a, x0, x1)
        a.grid(alpha=0.22, lw=0.6)
    ys = YSCALE.get((obs, resdir), ys)
    ax.set_yscale(ys)
    if ys == "log":
        pos = ry[mask & (ry > 0)] / div
        if pos.size:
            ax.set_ylim(0.03 * float(pos.min()), 4.0 * float(pos.max()))
    else:
        # room at the bottom so the lowest tick label does not sit on the
        # spine shared with the ratio panel (earlier pp04)
        y0, y1 = ax.get_ylim()
        ax.set_ylim(y0 - 0.10 * (y1 - y0), y1)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="lower"))
    axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="upper"))
    ax.set_ylabel(tex(rf"${ylab_q}$  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL - 2)
    ax.tick_params(labelbottom=False)
    return seen


def main():
    fig = plt.figure(figsize=(11.4, 15.8))
    outer = fig.add_gridspec(4, 2, hspace=0.30, wspace=0.22,
                             top=0.895, bottom=0.045, left=0.085, right=0.985)
    for i, (obs, xlab, xs, ys, ylab_q) in enumerate(OBS):
        pair = []
        for j, (resdir, rows, title, div, unit) in enumerate(COLUMNS):
            inner = outer[i, j].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                            hspace=0.0)
            ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            seen = block(ax, axr, resdir, rows, obs, div, unit, xs, ys,
                         ylab_q)
            pair.append((axr, seen or []))
            if i == 0:
                ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
            axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 1)
            if j == 0:
                axr.set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                               fontsize=plotstyle.FS_YLABEL - 4)
        for a, seen_here in pair:
            a.set_ylim(*ratio_range(seen_here))

    h, lab = [], []
    for a in fig.axes[:4]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    want = [tex(REF_LABEL), tex(ZM_LABEL),
            tex(MU_ROWS[0][0]), tex(MU_ROWS[1][0]), tex(NU_ROWS[2][0]),
            tex(MU_ROWS[2][0]), tex(NU_ROWS[3][0]), tex("NLO MHOU")]
    order = [lab.index(w) for w in want if w in lab]
    order += [i for i in range(len(lab)) if i not in order]
    h, lab = [h[i] for i in order], [lab[i] for i in order]
    # each band shares its curve's entry (user, 2026-10-04)
    h, lab = plotstyle.merge_bands(h, lab, {tex("NLO MHOU"): [tex(REF_LABEL)]})
    if h:
        # pp02's legend (user, 2026-09-14: "pp04 same legend as pp02
        # "): full legend font, three rows (8 entries, so 3 columns), grid
        # lowered.
        fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.978),
                   ncol=3, frameon=True, handlelength=2.3,
                   columnspacing=1.4, labelspacing=0.4,
                   fontsize=plotstyle.FS_LEGEND)
    fig.suptitle(tex(FIG_TITLE), y=0.995, fontsize=plotstyle.FS_SUPTITLE)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------- claims --
_COLS = [("results", MU_ROWS), ("results_nu", NU_ROWS)]


def _ratio(resdir, key, obs):
    ref, h = hist(resdir, REF, obs), hist(resdir, key, obs)
    m = _vis(resdir, obs)
    return h[1][m] / ref[1][m]


def _dev(resdir, key, obs):
    return float(np.max(np.abs(_ratio(resdir, key, obs) - 1.0)))


def _median(resdir, key, obs):
    return float(np.median(_ratio(resdir, key, obs)))


def _band_mean(resdir, obs):
    b = mhou(resdir, obs)
    m = _vis(resdir, obs)
    return float(0.5 * (np.asarray(b[1])[m] - np.asarray(b[2])[m]).mean())


def _bins(resdir, key, obs):
    """[(lower edge, ratio)] over the visible bins, in order."""
    e, y, _ = hist(resdir, key, obs)
    r = hist(resdir, REF, obs)[1]
    return [(float(e[i]), float(y[i] / r[i]))
            for i in np.nonzero(_vis(resdir, obs))[0]]


def _ends(resdir, key, obs):
    b = _bins(resdir, key, obs)
    return b[0][1], b[-1][1]


def _monotone_up(resdir, key, obs, slack=0.012):
    v = [r for _e, r in _bins(resdir, key, obs)]
    return all(b >= a - slack for a, b in zip(v, v[1:]))


MU_GRV, MU_NN, MU_PWG = (k for _l, k, *_ in MU_ROWS)
NU_GRV, NU_NN, NU_HEDIS, NU_PWG = (k for _l, k, *_ in NU_ROWS)


def _fmtb(resdir, key, obs):
    return ", ".join("%g:%.3f" % (e, r) for e, r in _bins(resdir, key, obs))


MESSAGE = (
    "<b>GENIE differentially, at 1 TeV in the final region</b>: tungsten per "
    "nucleon, Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no y cut, "
    "500k events per nucleon for each G18_02a row and 400k for HEDIS. The "
    "reference is FONLL with target-mass corrections; YADISM (ZM) is drawn "
    "beside it, within 1% on the muon side and 3% on the neutrino side, "
    "and POWHEG-RES and POWHEG-V2 stay within about 1% of the reference in "
    "the median bin of every panel."
    "\n\n"
    "<b>The default tune's per-cent agreement on the charged-current rate "
    "hides a tilt.</b> Its x<sub>Bj</sub> spectrum rises from 0.85 of the "
    "reference in the lowest bin to 1.08 in the highest, and Q&sup2; from "
    "0.87 to 1.03, while the integrated rate agrees to 0.6%: a shape error "
    "of about 20% that cancels in the total. The lepton energy, which at "
    "fixed beam energy is y, is flat to 2.5%. <b>On the neutral current the "
    "tilt is larger</b>: x<sub>Bj</sub> from 0.84 to 1.22, Q&sup2; from "
    "0.88 to 1.09, and the scattered-lepton energy 16% low at the bottom of "
    "its range, coming back to 1.03 at E<sub>&#8467;'</sub> &gt; 900 GeV."
    "\n\n"
    "<b>NNPDF4.0 in the same machinery</b> raises the charged-current "
    "spectra by 3.5&ndash;6% and the neutral-current Q&sup2; tail by up to "
    "21%. <b>HEDIS is a normalisation, not a shape</b>: 5% low and flat to "
    "a per cent in the lepton energy, with a residual 9% slope in "
    "x<sub>Bj</sub> at the smallest x. The seven-point band is 3&ndash;8% "
    "per bin on the muon side and 1&ndash;4% on the neutrino side."
)

CLAIMS = [
    {"what": "every file drawn is a 1 TeV tungsten per-nucleon q4w3 result; "
             "the reference is FONLL with TMC",
     "check": lambda: (REF.startswith("yadism_nlo_fonll_damp") and REF.endswith("_tmc")
                       and all(hist(d, k, o) is not None for d, rows in _COLS
                               for _l, k, *_ in rows for o, *_ in OBS)),
     "detail": lambda: REF},
    {"what": "every drawn histogram integrates to its own sigma_fid (1e-4 "
             "in x, E_lep', Q2)",
     "check": lambda: all(
         abs(float(np.sum(np.asarray(_load(d, k)["hists"][o]["dsig"])
                          * np.diff(_load(d, k)["hists"][o]["edges"]))
                   / _load(d, k)["sigma_fid_pb"]) - 1.0) < 2e-4
         for d, rows in _COLS for k in [r[1] for r in rows] + [REF, ZM_REF]
         for o in ("xbj", "Emu", "Q2")),
     "detail": lambda: "all rows, both references"},
    {"what": "ZM within 1% of FONLL in every muon bin and 3% on the neutrino side",
     "check": lambda: (max(_dev("results", ZM_REF, o) for o, *_ in OBS) < 0.010
                       and max(_dev("results_nu", ZM_REF, o) for o, *_ in OBS) < 0.03),
     "detail": lambda: "mu %.4f, nu %.4f" % (
         max(_dev("results", ZM_REF, o) for o, *_ in OBS),
         max(_dev("results_nu", ZM_REF, o) for o, *_ in OBS))},
    {"what": "POWHEG-RES and POWHEG-V2 within about 1% (1.1%) of the reference "
             "in the median bin of every panel",
     "check": lambda: all(abs(_median(d, k, o) - 1) < 0.011
                          for d, k in (("results", MU_PWG), ("results_nu", NU_PWG))
                          for o, *_ in OBS),
     "detail": lambda: ", ".join("%.4f" % _median(d, k, o)
                                 for d, k in (("results", MU_PWG), ("results_nu", NU_PWG))
                                 for o, *_ in OBS)},
    {"what": "CC default tune tilted: x from 0.85 to 1.08 and Q2 from 0.87 to "
             "1.03, rising bin by bin, while its lepton energy is flat to 2.5%",
     "check": lambda: (
         abs(_ends("results_nu", NU_GRV, "xbj")[0] - 0.85) < 0.01
         and abs(_ends("results_nu", NU_GRV, "xbj")[1] - 1.08) < 0.01
         and abs(_ends("results_nu", NU_GRV, "Q2")[0] - 0.87) < 0.01
         and abs(_ends("results_nu", NU_GRV, "Q2")[1] - 1.03) < 0.01
         and _monotone_up("results_nu", NU_GRV, "xbj")
         and _monotone_up("results_nu", NU_GRV, "Q2")
         and _dev("results_nu", NU_GRV, "Emu") < 0.025),
     "detail": lambda: "x " + _fmtb("results_nu", NU_GRV, "xbj")
     + " | Q2 " + _fmtb("results_nu", NU_GRV, "Q2")},
    {"what": "NC default tune: x from 0.84 to 1.22, Q2 from 0.88 to 1.09, "
             "E_lep' 0.84 in the lowest bin back to 1.03 above 900 GeV",
     "check": lambda: (
         abs(_ends("results", MU_GRV, "xbj")[0] - 0.84) < 0.01
         and abs(_ends("results", MU_GRV, "xbj")[1] - 1.22) < 0.01
         and abs(_ends("results", MU_GRV, "Q2")[0] - 0.88) < 0.01
         and abs(_ends("results", MU_GRV, "Q2")[1] - 1.09) < 0.01
         and abs(_ends("results", MU_GRV, "Emu")[0] - 0.84) < 0.01
         and abs(_ends("results", MU_GRV, "Emu")[1] - 1.03) < 0.01
         and _monotone_up("results", MU_GRV, "Emu")),
     "detail": lambda: "x " + _fmtb("results", MU_GRV, "xbj")
     + " | Emu " + _fmtb("results", MU_GRV, "Emu")},
    {"what": "NNPDF4.0 raises the CC spectra by 3.5-6% (median, every panel) and "
             "the NC Q2 tail to 1.21",
     "check": lambda: (all(1.03 <= _median("results_nu", NU_NN, o) <= 1.065
                           for o, *_ in OBS)
                       and abs(_ends("results", MU_NN, "Q2")[1] - 1.207) < 0.01),
     "detail": lambda: ", ".join("%.3f" % _median("results_nu", NU_NN, o) for o, *_ in OBS)
     + " | NC Q2 last %.3f" % _ends("results", MU_NN, "Q2")[1]},
    {"what": "HEDIS 5% low and flat to a per cent in E_lep' (0.943-0.954), "
             "with x_Bj falling to 0.87 in the lowest bin",
     "check": lambda: (all(0.94 <= r <= 0.956 for _e, r in _bins("results_nu", NU_HEDIS, "Emu"))
                       and abs(_ends("results_nu", NU_HEDIS, "xbj")[0] - 0.871) < 0.01
                       and all(0.9 <= _median("results_nu", NU_HEDIS, o) <= 0.96
                               for o, *_ in OBS)),
     "detail": lambda: "Emu " + _fmtb("results_nu", NU_HEDIS, "Emu")},
    {"what": "the seven-point band is 3-8% per bin on the muon side and 1-4% "
             "on the neutrino side (panel means)",
     "check": lambda: (all(0.025 <= _band_mean("results", o) <= 0.08 for o, *_ in OBS)
                       and all(0.008 <= _band_mean("results_nu", o) <= 0.042
                               for o, *_ in OBS)),
     "detail": lambda: "mu " + ", ".join("%.3f" % _band_mean("results", o) for o, *_ in OBS)
     + " | nu " + ", ".join("%.3f" % _band_mean("results_nu", o) for o, *_ in OBS)},
]


if __name__ == "__main__":
    main()
