#!/usr/bin/env python3
"""Paper plot 2: NLO differential DIS distributions at 1 TeV, in the FINAL
benchmark region.

The redo of the earlier figure for the "Paper plots"
tab (user, 2026-09-14, after pp01).  The settings are pp01's
(PRODUCTION.md):

  * TUNGSTEN PER NUCLEON, (74 sigma_p + 110 sigma_n)/184, every generator run
    on a genuine proton and a genuine neutron beam; the references use the
    free-nucleon average set, NNPDF40_nnlo_as_01180_W184free.
  * Q2 > 4 GeV2 AND W > 3 GeV, NO y CUT, on BOTH currents -- the earlier production's muon floor
    of Q2 > 11 does not carry over, so the two columns share one region again.
  * ONE ENERGY, 1 TeV (user, 2026-09-13: "differential distributions only at
    1 TeV"), which is why the 1 TeV samples are the high-statistics ones.
  * The YADISM references carry TARGET-MASS CORRECTIONS (TMC on), and so do
    both scale bands.

Everything else is the earlier production's: the four leptonic observables the analytic
calculation can reach (x_Bj, E_lep', theta_lep', Q2), the three NLO matchings
per current in the earlier production's colours and styles, the ratio to YADISM (ZM) with YADISM
(FONLL) beside it, the NNLO curve in the scheme each current has at that order
(FONLL for the muon, ZM for the neutrino), and a seven-point band on each
order, bin by bin, keeping the xiF = 1/2 points (analysis/mhou_diff.py).

>>> THE LEPTON-ENERGY ROW NEEDED NEW BINS (found 2026-09-14). <<<  The y-like
axes (y, nu, E_lep', E_had) were literals for 0.2 < y < 0.9.  In this region
the generator E_lep' histograms held 35% (muon) and 70% (neutrino) of
sigma_fid, np.histogram dropping the rest, while the YADISM calculators clip
out-of-range entries into the edge bins -- its first y bin came out 24x too
high.  analyze.bins_for() now follows the selection's y window, and the 1 TeV
generator results, references and bands were all re-derived on it.  A claim
below checks that every histogram drawn integrates to its own sigma_fid, which
is the test that would have caught it.

GENIE IS DELIBERATELY ABSENT, as in the earlier production and for the earlier production's reason: this figure asks
whether the NLO matchings agree differentially, and GENIE, a leading-order
calculation with a 6-13% offset, would set the ratio range and hide the
per-cent effects the figure is about.  The paper-plots carve-out of CONVENTIONS.md
rule 1b applies, stated here as it requires.

Usage: analysis/paper_plots/pp02_dis_distributions.py
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

SLUG = "pp02_dis_distributions"
V1_NUMBER = 2
# IN THE PAPER, in place of the earlier production's figure (user, 2026-09-14: "update the paper
# ... when other plots are updated").  Read by
# tools/check_paper_figures_used.py.
IN_PAPER = True
TITLE = "NLO differential DIS distributions at 1 TeV"
FIG_TITLE = r"NLO differential DIS distributions at $E_{\ell} = 1$ TeV"
# NOT "paper_*": tools/sync_paper_figures.sh copies every paper_*.png into the
# tracked paper/figures/ by glob; a figure enters the paper by name.
OUTPUT = "pp_dis_distributions.png"
RESULTS = "results_nu"
CAPTION = ("Differential distributions per nucleon on tungsten at 1 TeV, for "
           "Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no y cut, in "
           "Bjorken <i>x</i>, scattered-lepton energy and angle, and momentum "
           "transfer, for the three NLO matchings; muon neutral current "
           "(left) and neutrino charged current (right). Ratio to YADISM "
           "(ZM) with target-mass corrections beneath each, with YADISM "
           "(FONLL) and the NNLO result drawn as well &mdash; the latter in "
           "FONLL on the muon side and in ZM on the neutrino side, where the "
           "massive coefficient functions, known at NNLO, are not implemented "
           "in the structure-function code used here. The shaded bands are the "
           "seven-point scale variation on each order, bin by bin.")

SURFACE = "#ffffff"
ENERGY = 1000.0
REGION = "q4w3"
Q2MIN = 4.0
W_MIN = 3.0
TAG = f"_{REGION}_W"            # generators: tungsten per nucleon
TMC = "_tmc"                    # references: TMC on (user, 2026-09-13)

REF = f"yadism_nlo{TAG}{TMC}"
REF_NOTMC = f"yadism_nlo{TAG}"
REF_LABEL = "YADISM (ZM)"
FONLL_REF = f"yadism_nlo_fonll_damp{TAG}{TMC}"
FONLL_LABEL = "YADISM (FONLL)"
REF_COLOUR = "#111111"
FONLL_COLOUR = "#c0392b"

NNLO = {"results": (f"yadism_nnlo_fonll_damp{TAG}{TMC}",
                    r"YADISM NNLO (FONLL, $\mu$)"),
        "results_nu": (f"yadism_nnlo{TAG}{TMC}", r"YADISM NNLO (ZM, $\nu$)")}
NNLO_COLOUR = "#6a3d9a"
NNLO_LS = {"results": "-.", "results_nu": (0, (5.0, 1.3, 1.0, 1.3, 1.0, 1.3))}

# The 7-point bands, bin by bin, each in the scheme of the curve it surrounds
# (analysis/mhou_diff.py with BENCH_SELECTION=q4w3 BENCH_TARGET=W BENCH_TMC=3).
MHOU_FILE = f"mhou_diff_zm{TAG}{TMC}.json"
NNLO_MHOU_FILE = {"results": f"mhou_diff_fonll_pto2{TAG}{TMC}.json",
                  "results_nu": f"mhou_diff_zm_pto2{TAG}{TMC}.json"}
MHOU_KEEP_HALF = True

# the earlier production's rows, the earlier production's styling -- one legend, so the two POWHEGs differ in style.
POWHEG_V2_LS = (0, (6.5, 1.6))
NU_ROWS = [
    (r"POWHEG-V2 ($\nu$)", f"powheg_nu{TAG}",             "#0072b2", "o",
     POWHEG_V2_LS),
    ("Herwig",     f"herwig_nlo_full{TAG}",               "#d55e00", "s", "--"),
    ("Sherpa",     f"sherpa_nlo{TAG}",                    "#009e73", "^", "-."),
]
MU_ROWS = [
    (r"POWHEG-RES ($\mu$)", f"powheg{TAG}",               "#0072b2", "o", "-"),
    ("Herwig",     f"herwig_nlo_powheg_full{TAG}",        "#d55e00", "s", "--"),
    ("Sherpa",     f"sherpa{TAG}",                        "#009e73", "^", "-."),
]

# X RANGES: the user's earlier windows (2026-08-30, revised 09-01), snapped to bin
# edges where used.  With no y cut the lepton-energy row now spans its full
# kinematic range, 0 to E.
XLIM = {
    "theta": {"results": (None, 3e-2), "results_nu": (None, 8e-2)},
    "Q2": {"results": (None, 100.0), "results_nu": (None, 400.0)},
    "xbj": {"results": (None, 0.6), "results_nu": (None, 0.6)},
}

# Merges applied to EVERY curve of a panel, reference included: the earlier production's.  E_lep'
# 20 GeV -> 100 GeV, x_Bj and neutrino Q2 in pairs; the muon Q2 row unmerged
# (14 bins over 4-100 GeV2; merged in pairs it had 7, too few for a fall of
# three decades).
REBIN = {("Emu", "results_nu"): 5, ("Emu", "results"): 5,
         ("Q2", "results_nu"): 2,
         ("xbj", "results"): 2, ("xbj", "results_nu"): 2}

# THE MUON LEPTON-ENERGY PANEL IS LOGARITHMIC.  With no y cut the neutral
# current piles up at y -> 0, so its last 100 GeV bin holds ~2.5x the next and
# a linear axis flattens every other bin onto zero.  The neutrino spectrum
# varies by 30% over the range and stays linear.
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

# Column titles are the current alone (user, 2026-09-14, on pp01): the
# target and the cuts are in the caption.
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
    """(edges, dsig, err) for one observable, or None if absent."""
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
        # every curve on a panel must be binned like the reference: a ratio
        # of differently binned histograms is not a ratio of anything
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

    nnlo_key, nnlo_label = NNLO[resdir]
    nh = hist(resdir, nnlo_key, obs)
    if nh is not None:
        edges_of(nh, nnlo_key)
        ax.stairs(nh[1] / div, edges, color=NNLO_COLOUR, lw=1.8,
                  ls=NNLO_LS[resdir], baseline=None, zorder=4,
                  label=tex(nnlo_label))
        with np.errstate(divide="ignore", invalid="ignore"):
            nr = np.where(ry > 0, nh[1] / ry, np.nan)
        axr.stairs(nr, edges, color=NNLO_COLOUR, lw=1.8,
                   ls=NNLO_LS[resdir], baseline=None, zorder=4)
        seen.append((nr, mask))
        # the NNLO band sits around the NNLO CURVE; a flat fill, not hatched
        # (user, 2026-09-18: drawn as in every other figure of the paper)
        nb = mhou(resdir, obs, NNLO_MHOU_FILE[resdir])
        if nb is not None:
            edges_of(nb, NNLO_MHOU_FILE[resdir])
            axr.stairs(nr * (1.0 + nb[1]), edges,
                       baseline=nr * (1.0 + nb[2]), fill=True,
                       color=NNLO_COLOUR, alpha=0.20, lw=0, zorder=2,
                       label=tex("NNLO MHOU"))
            seen.append((nr * (1.0 + nb[1]), mask))
            seen.append((nr * (1.0 + nb[2]), mask))

    fon = hist(resdir, FONLL_REF, obs)
    if fon is not None:
        edges_of(fon, FONLL_REF)
        ax.stairs(fon[1] / div, edges, color=FONLL_COLOUR, lw=1.7, ls="-",
                  baseline=None, zorder=6, label=tex(FONLL_LABEL))
        with np.errstate(divide="ignore", invalid="ignore"):
            fr = np.where(ry > 0, fon[1] / ry, np.nan)
        axr.stairs(fr, edges, color=FONLL_COLOUR, lw=1.7, ls="-",
                   baseline=None, zorder=6)
        seen.append((fr, mask))

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
    if ys != "log":
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=5, prune="lower"))
    axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="upper"))
    ax.set_ylabel(tex(rf"${ylab_q}$  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL - 2)
    ax.tick_params(labelbottom=False)
    return seen


def main():
    fig = plt.figure(figsize=(11.4, 15.8))
    outer = fig.add_gridspec(4, 2, hspace=0.30, wspace=0.22,
                             top=0.905, bottom=0.045, left=0.085, right=0.985)
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
        # every ratio panel sets its own range (user, 2026-09-01)
        for a, seen_here in pair:
            a.set_ylim(*ratio_range(seen_here))

    # ONE LEGEND above the figure, de-duplicated by label, in a fixed order
    h, lab = [], []
    for a in fig.axes[:4]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    want = [tex(REF_LABEL), tex(FONLL_LABEL),
            tex(NNLO["results"][1]), tex(NNLO["results_nu"][1]),
            tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0]),
            tex(MU_ROWS[1][0]), tex(MU_ROWS[2][0]),
            tex("NLO MHOU"), tex("NNLO MHOU")]
    order = [lab.index(w) for w in want if w in lab]
    order += [i for i in range(len(lab)) if i not in order]
    h, lab = [h[i] for i in order], [lab[i] for i in order]
    # each band shares its curve's entry (user, 2026-10-04)
    h, lab = plotstyle.merge_bands(h, lab, {
        tex("NLO MHOU"): [tex(REF_LABEL)],
        tex("NNLO MHOU"): [tex(NNLO["results"][1]), tex(NNLO["results_nu"][1])]})
    if h:
        # THREE ROWS AT A LARGER FONT (user, 2026-09-14: "increase the font
        # of the legends ... maybe three rows better than two"): 10 entries
        # in 4 columns, the main grid lowered to make room.
        fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.978),
                   ncol=4, frameon=True, handlelength=2.3,
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
    """Largest |ratio - 1| over the VISIBLE bins of one curve."""
    return float(np.max(np.abs(_ratio(resdir, key, obs) - 1.0)))


def _panel_median_of(resdir, obs):
    """Largest median |ratio - 1| among the generators of one panel."""
    return max(float(np.median(np.abs(_ratio(resdir, k, obs) - 1.0)))
               for _, k, *_r in dict(_COLS)[resdir])


def _worst_gen():
    """(deviation, current, generator, observable) of the worst visible bin."""
    return max((_dev(d, k, o), d, lab, o) for d, rows in _COLS
               for lab, k, *_r in rows for o, *_ in OBS)


def _worst_fonll(resdir=None):
    cols = [c for c in _COLS if resdir in (None, c[0])]
    return max(_dev(d, FONLL_REF, o) for d, _ in cols for o, *_ in OBS)


def _coverage(resdir, key):
    """min, max over the four observables of (histogram integral)/sigma_fid,
    on the RAW histograms and over the full range, not the window."""
    d = _load(resdir, key)
    out = []
    for o, *_ in OBS:
        h = d["hists"][o]
        out.append(float(np.sum(np.asarray(h["dsig"]) * np.diff(h["edges"]))
                         / d["sigma_fid_pb"]))
    return min(out), max(out)


def _cut_rate(resdir, obs):
    """Fraction of the reference rate above the panel's upper window edge."""
    edges, y, _ = hist(resdir, REF, obs)
    w = np.diff(edges)
    hi = edges[:-1] >= x_hi(obs, resdir, edges) - 1e-12
    return float((y * w)[hi].sum() / (y * w).sum())


def _band_mean(resdir, obs, fn=None):
    """Mean half-width of a 7-point band over the visible bins."""
    b = mhou(resdir, obs, fn)
    m = _vis(resdir, obs)
    return float(0.5 * (np.asarray(b[1])[m] - np.asarray(b[2])[m]).mean())


def _nnlo_band_narrowing(resdir):
    """NLO band width / NNLO band width per panel, over the visible bins."""
    return [_band_mean(resdir, o) / _band_mean(resdir, o, NNLO_MHOU_FILE[resdir])
            for o, *_ in OBS]


def _nnlo_range(resdir, obs):
    v = _ratio(resdir, NNLO[resdir][0], obs)
    return float(v.min()), float(v.max())


def _band_central_matches(resdir, fn, key):
    """The band's central histogram IS the curve it is drawn around."""
    with open(f"{BASE}/{resdir}/{fn}") as f:
        b = json.load(f)
    d = _load(resdir, key)
    worst = 0.0
    for o, *_ in OBS:
        c = np.asarray(b["hists"][o]["dsig"])
        r = np.asarray(d["hists"][o]["dsig"])
        m = r > 0
        worst = max(worst, float(np.max(np.abs(c[m] / r[m] - 1.0))))
    return worst


def _ebeam(resdir, key):
    """Beam energy a sample reconstructs, median(E_lep') / (1 - median(y))."""
    out = []
    for o in ("Emu", "y"):
        h = raw_hist(resdir, key, o)
        c = np.cumsum(h[1] * np.diff(h[0]))
        out.append(float(np.interp(0.5, c / c[-1], h[0][1:])))
    return out[0] / (1.0 - out[1])


def _bin_ratios(resdir, key, obs):
    """(lower edge, upper edge, ratio) for every visible bin of one curve."""
    e, y, _ = hist(resdir, key, obs)
    r = hist(resdir, REF, obs)[1]
    return [(e[i], e[i + 1], y[i] / r[i])
            for i in np.nonzero(_vis(resdir, obs))[0]]


def _pull(resdir, key, obs, i):
    """(ratio - 1) / relative error of the i-th visible bin of one curve."""
    e, y, er = hist(resdir, key, obs)
    r = hist(resdir, REF, obs)[1]
    b = np.nonzero(_vis(resdir, obs))[0][i]
    return (y[b] / r[b] - 1) / (er[b] / r[b])


def _tmc_max(resdir, obs):
    """Largest |TMC-on / TMC-off - 1| of the reference over visible bins."""
    m = _vis(resdir, obs)
    a, b = hist(resdir, REF, obs)[1], hist(resdir, REF_NOTMC, obs)[1]
    return float(np.max(np.abs(a[m] / b[m] - 1.0)))


def _stat(resdir, key):
    """Largest over the four panels of the median relative per-bin error."""
    out = 0.0
    for o, *_ in OBS:
        _, y, e = hist(resdir, key, o)
        m = _vis(resdir, o)
        out = max(out, float(np.median(e[m] / y[m])))
    return out


MESSAGE = (
    "<b>The final region at 1 TeV: tungsten per nucleon, Q&sup2; &gt; 4 "
    "GeV&sup2; and W &gt; 3 GeV with no y cut, on both currents.</b> The two "
    "columns share one region, "
    "and without a y window the lepton-energy row spans its full range, "
    "0 to 1 TeV. The reference and both bands include target-mass "
    "corrections, which move no visible bin by more than 1% (the muon's "
    "highest-x bin) and none by more than 0.2% on the neutrino side."
    "\n\n"
    "<b>The per-cent agreement survives differentially.</b> In the median "
    "bin of every panel each NLO matching is within 1.6% of YADISM (ZM) on "
    "the muon side and 0.8% on the neutrino side. Herwig is the closest on "
    "the muon side, within 0.8%, and on the neutrino side it and POWHEG-V2 "
    "are both within 0.5%. Herwig's x<sub>Bj</sub> spectrum follows the "
    "reference to 1.1% on the muon side, apart from 3% low in the "
    "highest-x bin, and to 2.1% on the neutrino side; it was tilted by up "
    "to 9% while Herwig replaced the events whose beam remnant it could not "
    "build by new ones, instead of re-showering them on the same hard "
    "process. <b>One shape difference is outside the statistics, and it "
    "comes from the shower.</b> POWHEG-V2 is 9% high in the lowest "
    "neutrino x bin and 3% high in the next, each more than three "
    "standard deviations, where Herwig and Sherpa agree with the "
    "reference. Its hard events follow the reference there to 2%: Pythia, "
    "giving the massless charm and bottom quarks of the hard event their "
    "mass, lowers the lepton's Q&sup2; in 7% of the events, every one of "
    "them with a heavy quark, and moves them to lower x. The samples are "
    "generated with POWHEG-V2's correction for events above its "
    "unweighting bound; without it the hard events were 6&ndash;8% short "
    "at x = 0.005&ndash;0.011, which had cancelled this excess. The largest "
    "single deviation elsewhere is "
    "statistical: the muon Q&sup2; and angle tails, where Sherpa's "
    "per-bin error reaches 3&ndash;4%."
    "\n\n"
    "<b>The charm-mass effect</b>, FONLL against ZM, stays below 1% in "
    "every muon bin and below 2.7% on the neutrino side, largest at the "
    "lowest x. <b>Both scale bands narrow from NLO to NNLO in all eight "
    "panels</b>, by 1.5&ndash;1.9 on the muon side and 1.2&ndash;1.9 on "
    "the neutrino side, where the lepton-energy band is about 1% wide at "
    "either order. The muon NNLO correction reaches &minus;12% at the "
    "lowest x and &minus;10% at the lowest lepton energy. Each band is in "
    "the scheme of the curve it surrounds and keeps the &xi;<sub>F</sub> = "
    "&frac12; points."
    "\n\n"
    "<b>What the windows hide</b>: under 1.5% of either rate above "
    "x = 0.6, 0.6% of the muon rate above Q&sup2; = 100 GeV&sup2;, and 12% "
    "of the neutrino rate above 400 GeV&sup2;."
)

_GEN = [(d, k) for d, rows in _COLS for _l, k, *_r in rows]

CLAIMS = [
    {"what": "every file drawn is a 1 TeV tungsten per-nucleon result in the "
             "q4w3 region, and the references carry TMC and the W184free set",
     "check": lambda: all(hist(d, k, o) is not None
                          for d, rows in _COLS
                          for k in [r[1] for r in rows]
                          + [REF, FONLL_REF, NNLO[d][0], REF_NOTMC]
                          for o, *_ in OBS),
     "detail": lambda: "2 currents x 7 curves x 4 observables"},
    # THE TEST THAT WOULD HAVE CAUGHT THE y-WINDOW BINS (2026-09-14)
    {"what": "every histogram drawn integrates to its own sigma_fid: to 1e-4 "
             "in x, E_lep' and Q2 and to 0.25% in theta, whose last edge "
             "sits below the largest angles",
     "check": lambda: all(
         _coverage(d, k)[0] > 0.9975 and _coverage(d, k)[1] < 1.0001
         and all(abs(float(np.sum(np.asarray(_load(d, k)["hists"][o]["dsig"])
                                  * np.diff(_load(d, k)["hists"][o]["edges"]))
                           / _load(d, k)["sigma_fid_pb"]) - 1.0) < 1e-4
                 for o in ("xbj", "Emu", "Q2"))
         for d, k in _GEN + [(d, x) for d, _r in _COLS
                             for x in (REF, FONLL_REF, NNLO[d][0])]),
     "detail": lambda: ", ".join("%s %.4f" % (k, _coverage(d, k)[0])
                                 for d, k in _GEN)},
    {"what": "each band's central histogram IS the curve it is drawn around",
     "check": lambda: all(_band_central_matches(d, MHOU_FILE, REF) < 1e-9
                          and _band_central_matches(d, NNLO_MHOU_FILE[d],
                                                    NNLO[d][0]) < 1e-9
                          for d, _r in _COLS),
     "detail": lambda: "NLO ZM and NNLO bands, both currents"},
    {"what": "every sample reconstructs its own 1 TeV beam to 0.2% (the "
             "median estimate is limited by the 20 GeV bins), so the "
             "lepton-energy row is binned against the right energy",
     "check": lambda: all(abs(_ebeam(d, k) / ENERGY - 1.0) < 2e-3
                          for d, k in _GEN),
     "detail": lambda: ", ".join("%.1f" % _ebeam(d, k) for d, k in _GEN)},
    {"what": "TMC moves no visible reference bin by more than 1% (muon) or "
             "0.2% (neutrino)",
     "check": lambda: (max(_tmc_max("results", o) for o, *_ in OBS) < 0.010
                       and max(_tmc_max("results_nu", o) for o, *_ in OBS)
                       < 0.002),
     "detail": lambda: " | ".join(d + " " + ", ".join(
         "%.2f%%" % (100 * _tmc_max(d, o)) for o, *_ in OBS)
         for d, _r in _COLS)},
    {"what": "median bin of every panel: each matching within 1.6% (muon) "
             "and 0.8% (neutrino); Herwig the closest on the muon side, "
             "within 0.8%, and it and POWHEG-V2 both within 0.5% on the "
             "neutrino side (until 2026-09-19 Herwig's x panels were the "
             "worst, 2.7% and 2.2%, and Herwig was the closest on both)",
     "check": lambda: (
         max(_panel_median_of("results", o) for o, *_ in OBS) < 0.016
         and max(_panel_median_of("results_nu", o) for o, *_ in OBS) < 0.008
         and max(float(np.median(np.abs(_ratio("results", MU_ROWS[1][1], o) - 1)))
                 for o, *_ in OBS) < 0.008
         and all(max(float(np.median(np.abs(_ratio("results_nu", NU_ROWS[i][1], o) - 1)))
                     for o, *_ in OBS) < 0.005 for i in (0, 1))
         and max(float(np.median(np.abs(_ratio("results", MU_ROWS[1][1], o) - 1)))
                 for o, *_ in OBS)
             <= min(max(float(np.median(np.abs(_ratio("results", k, o) - 1)))
                        for o, *_ in OBS)
                    for lab, k, *_r in MU_ROWS if k != MU_ROWS[1][1])),
     "detail": lambda: " | ".join(d + " " + ", ".join(
         "%.2f%%" % (100 * _panel_median_of(d, o)) for o, *_ in OBS)
         for d, _r in _COLS)},
    {"what": "Herwig's x spectrum follows the reference to 1.1% on the muon "
             "side apart from the highest-x bin (2-4% low) and to 2.2% on "
             "the neutrino side (until 2026-09-19: tilted, 2-3.5% high at "
             "low x and 4-9% low above x = 0.15)",
     "check": lambda: (
         all(0.989 <= r <= 1.012 for lo, hi, r in
             _bin_ratios("results", MU_ROWS[1][1], "xbj")[:-1])
         and 0.96 <= _bin_ratios("results", MU_ROWS[1][1], "xbj")[-1][2] <= 0.98
         and all(0.978 <= r <= 1.022 for lo, hi, r in
                 _bin_ratios("results_nu", NU_ROWS[1][1], "xbj"))),
     "detail": lambda: " | ".join(", ".join(
         "%.3g:%.3f" % (lo, r) for lo, hi, r in _bin_ratios(d, rows[1][1],
                                                           "xbj"))
         for d, rows in _COLS)},
    # THE DIP WAS GENERATION, THE EXCESS IS THE SHOWER (2026-09-19).  Until the
    # POWHEG-V2 samples were regenerated with ubexcess_correct 1 this claim read
    # "3.5-5% low at x = 0.005-0.008, 1.5-3% low in the next bin, 1-2% low
    # below Q2 = 10": events above the unweighting bound were under-produced.
    # With them restored the hard events follow YADISM (ZM) to 2% in these
    # bins (PRODUCTION.md), and the excess left is Pythia re-massing c/b,
    # which lowers the lepton's Q2 in 6.9% of events, all with a heavy quark.
    {"what": "POWHEG-V2 is 7-11% high in the lowest neutrino x bin and 2-4% "
             "high in the next, each more than 3 sigma, while Herwig and "
             "Sherpa are within 3% there; no POWHEG-V2 x bin is more than 1% "
             "low and its Q2 bins below 10 GeV2 are within 1.5% (the dip "
             "of the un-corrected samples is gone)",
     "check": lambda: (
         1.07 <= _bin_ratios("results_nu", NU_ROWS[0][1], "xbj")[0][2] <= 1.11
         and 1.02 <= _bin_ratios("results_nu", NU_ROWS[0][1], "xbj")[1][2] <= 1.04
         and all(_pull("results_nu", NU_ROWS[0][1], "xbj", i) > 3 for i in (0, 1))
         and all(abs(_bin_ratios("results_nu", k, "xbj")[i][2] - 1) < 0.03
                 for _l, k, *_r in NU_ROWS[1:] for i in (0, 1))
         and min(r for _a, _b, r in _bin_ratios("results_nu", NU_ROWS[0][1], "xbj")) > 0.99
         and all(abs(r - 1) < 0.015 for lo, hi, r in
                 _bin_ratios("results_nu", NU_ROWS[0][1], "Q2") if hi <= 10.1)),
     "detail": lambda: ", ".join("%.3g:%.3f (%.1f sigma)" % (
         lo, r, _pull("results_nu", NU_ROWS[0][1], "xbj", i)) for i, (lo, hi, r) in
         enumerate(_bin_ratios("results_nu", NU_ROWS[0][1], "xbj")[:3]))},
    {"what": "median per-bin statistical error: every matching under 1.4% "
             "(muon) and 0.9% (neutrino) in every panel",
     "check": lambda: (all(_stat("results", k) < 0.014 for _l, k, *_r in MU_ROWS)
                       and all(_stat("results_nu", k) < 0.009
                               for _l, k, *_r in NU_ROWS)),
     "detail": lambda: ", ".join("%.2f%%" % (100 * _stat(d, k))
                                 for d, k in _GEN)},
    {"what": "charm-mass effect below 1% in every muon bin and below 2.7% "
             "on the neutrino side, where it is largest at the lowest x",
     "check": lambda: (_worst_fonll("results") < 0.010
                       and _worst_fonll("results_nu") < 0.027
                       and _dev("results_nu", FONLL_REF, "xbj")
                       == _worst_fonll("results_nu")),
     "detail": lambda: "mu %.2f%%, nu %.2f%%" % (
         100 * _worst_fonll("results"), 100 * _worst_fonll("results_nu"))},
    {"what": "both bands narrow from NLO to NNLO in all eight panels, by "
             "1.5-1.9 (muon) and 1.2-1.95 (neutrino)",
     "check": lambda: (all(1.45 < f < 1.9
                           for f in _nnlo_band_narrowing("results"))
                       and all(1.15 < f < 1.95
                               for f in _nnlo_band_narrowing("results_nu"))),
     "detail": lambda: "mu " + ", ".join(
         "%.2f" % f for f in _nnlo_band_narrowing("results"))
         + "; nu " + ", ".join("%.2f" % f
                               for f in _nnlo_band_narrowing("results_nu"))},
    {"what": "the neutrino lepton-energy band is about 1% at both orders",
     "check": lambda: (0.007 < _band_mean("results_nu", "Emu") < 0.012
                       and 0.006 < _band_mean("results_nu", "Emu",
                                              NNLO_MHOU_FILE["results_nu"])
                       < 0.012),
     "detail": lambda: "%.2f%% / %.2f%%" % (
         100 * _band_mean("results_nu", "Emu"),
         100 * _band_mean("results_nu", "Emu", NNLO_MHOU_FILE["results_nu"]))},
    {"what": "muon NNLO/NLO reaches -12% at the lowest x and -10% at the "
             "lowest lepton energy",
     "check": lambda: (abs(_nnlo_range("results", "xbj")[0] - 0.88) < 0.01
                       and abs(_nnlo_range("results", "Emu")[0] - 0.90) < 0.01
                       and _ratio("results", NNLO["results"][0], "xbj")[0]
                       == _nnlo_range("results", "xbj")[0]
                       and _ratio("results", NNLO["results"][0], "Emu")[0]
                       == _nnlo_range("results", "Emu")[0]),
     "detail": lambda: "x %.3f, E_lep' %.3f" % (
         _nnlo_range("results", "xbj")[0], _nnlo_range("results", "Emu")[0])},
    {"what": "the windows hide under 1.5% of either rate above x = 0.6, 0.6% "
             "of the muon rate above Q2 = 100 and 12% of the neutrino rate "
             "above Q2 = 400",
     "check": lambda: (all(_cut_rate(d, "xbj") < 0.015 for d, _r in _COLS)
                       and _cut_rate("results", "Q2") < 0.007
                       and 0.11 < _cut_rate("results_nu", "Q2") < 0.13),
     "detail": lambda: "x %.2f%% / %.2f%%, Q2 %.2f%% / %.2f%%" % (
         100 * _cut_rate("results", "xbj"), 100 * _cut_rate("results_nu", "xbj"),
         100 * _cut_rate("results", "Q2"), 100 * _cut_rate("results_nu", "Q2"))},
]


if __name__ == "__main__":
    main()
