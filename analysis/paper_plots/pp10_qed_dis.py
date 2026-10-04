#!/usr/bin/env python3
"""Paper plot 10: QED corrections to the DIS observables, in the FINAL region.

The redo of the earlier figure for the "Paper plots" tab: the
same events re-showered with QED radiation switched on in stages
(powheg/qed/: fsr, fsrisr on the muon side only, full), each arm against the
sample with no QED.  The conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184, each arm run on p and n.
  * THE REGION is Q2 > 4 GeV2 and W > 3 GeV with no y cut (selection q4w3).
  * THE NEUTRINO ARMS RE-SHOWER POWHEG-V2 (user, 2026-09-14), the benchmark's
    own neutrino NLO row, where earlier used POWHEG-RES on both currents.  Each
    column is a ratio to its own no-QED sample, so the generator cancels.
  * THE DENOMINATOR IS THE SAME EVENTS: the arms re-shower a SUBSET of the 
    1 TeV production (first ten POWHEG-RES seeds, first two POWHEG-V2
    batches; powheg/production/qed_arms.sh), and the no-QED curve is the production
    restricted to those jobs ("powheg_sub", sample_layout.arm_files).  Pythia runs
    from a fixed seed, and the empty `off` overlay was checked to reproduce
    the production byte for byte on both currents, so the arms stay
    correlated with the denominator event by event until QED changes one.

Usage: analysis/paper_plots/pp10_qed_dis.py
"""
import json
import math
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
from matplotlib.legend_handler import HandlerBase           # noqa: E402
from matplotlib.lines import Line2D                         # noqa: E402
import beams                                                # noqa: E402
import mhou_band                                            # noqa: E402

SLUG = "pp10_qed_dis"
V1_NUMBER = 10
IN_PAPER = True
TITLE = "QED corrections to the DIS observables"
FIG_TITLE = (r"QED corrections to the DIS observables, "
             r"$E_{\ell} = 1$ TeV")
OUTPUT = "pp_qed_dis.png"
RESULTS = "results_nu"
CAPTION = ("Differential distributions in Bjorken <i>x</i>, the scattered "
           "lepton&rsquo;s energy and polar angle, and the momentum "
           "transfer, per nucleon on tungsten at 1 TeV in the region Q&sup2; &gt; 4 GeV&sup2;, W &gt; 3 GeV, for the same "
           "Les Houches events showered with QED radiation switched on in "
           "stages (POWHEG-RES for muons, POWHEG-V2 for neutrinos). Muon neutral current (left), neutrino charged current "
           "(right); the lower panels are the ratio to the sample with no "
           "QED radiation, with its statistical error as the band. "
           "Initial-state radiation off the lepton line exists only on the "
           "muon side, a charged incoming lepton being required for it, so "
           "the neutrino column carries two arms where the muon column "
           "carries three.")

SURFACE = "#ffffff"
ENERGY = beams.ANCHOR_ENERGY
SELECTION = "q4w3"
TSUF = "_W"
Q2MIN = 4.0

# (label, result key, colour, marker, linestyle).  Okabe-Ito blue / vermillion
# / bluish-green, separable under every common form of colour blindness, with
# the line style repeating the information so the figure survives greyscale.
# The two POWHEGs are the two columns' denominators and never share a panel,
# so they carry different line styles and both appear in the one legend --
# the arrangement paper plot 2 settled on.
POWHEG_V2_LS = (0, (6.5, 1.6))
MU_ROWS = [
    (r"no QED ($\mu$)",      "powheg_sub",          "#0072b2", "o", "-"),
    ("QED FSR",              "powheg_qedfsr",       "#d55e00", "s", "--"),
    ("QED FSR + ISR",        "powheg_qedfsrisr",    "#cc79a7", "v", ":"),
    ("QED full",             "powheg_qedfull",      "#009e73", "^", "-."),
]
NU_ROWS = [
    (r"no QED ($\nu$)",      "powheg_nu_sub",       "#0072b2", "o",
     POWHEG_V2_LS),
    ("QED FSR",              "powheg_nu_qedfsr",       "#d55e00", "s", "--"),
    ("QED full",             "powheg_nu_qedfull",      "#009e73", "^", "-."),
]
# The denominator of each column's ratio panel: its own POWHEG.
# THE DENOMINATOR IS THE NO-QED ARM OF EACH COLUMN, and on the neutrino side
# that is POWHEG-RES and not the benchmark's POWHEG-V2 row: the arms
# re-shower POWHEG-RES events, so anything else would be comparing two
# different matrix elements and calling the difference QED.
DENOM = {"results": "powheg_sub", "results_nu": "powheg_nu_sub"}
# ...and what to CALL it on the ratio axis.  Never a bare "POWHEG"
# (CONVENTIONS.md rule 3), and never the column's usual code where the figure
# uses another one.
DENOM_NAME = {"powheg_sub": "POWHEG-RES", "powheg_nu_sub": "POWHEG-V2",
              "powheg": "POWHEG-RES", "powheg_res_nu": "POWHEG-RES",
              "powheg_nu": "POWHEG-V2", "powheg_v2": "POWHEG-V2"}

# "MUON-NEUTRINO" (user, 2026-10-04): QED radiation off the final-state
# lepton depends on its mass, so the flavour matters on the QED figures.
COLUMNS = [("results", MU_ROWS, "Muon DIS", 1000.0, "nb"),
           ("results_nu", NU_ROWS, "Muon-neutrino DIS", 1.0, "pb")]

# (key, x label, x scale, y scale, dsigma label)
OBS = [
    ("xbj",   r"$x_{\rm Bj}$", "log", "log", r"d\sigma/dx_{\rm Bj}"),
    ("Emu",   r"$E_{\ell'}$  [GeV]", "linear", "linear",
     r"d\sigma/dE_{\ell'}"),
    ("theta", r"$\theta_{\ell'}$  [rad]", "log", "log",
     r"d\sigma/d\theta_{\ell'}"),
    ("Q2",    r"$Q^2$  [GeV$^2$]", "log", "log", r"d\sigma/dQ^2"),
]

# X WINDOWS.  Each is the region that carries the rate; anything unlisted
# keeps the full histogram.  The multiplicity axis stops where the samples
# run out of events rather than at the last bin edge, and the Delta phi axis
# starts at the cut, since the selection populates nothing below it.
# >>> FIG. 4.4's BINNING AND WINDOWS (user, 2026-10-04: "use the same
# binning and x-axis range as Fig. 4.4 for consistency"). <<<  Copied from
# pp04_genie_dis_distributions.py: per-current upper edges, merges, and the
# lower edge read from the data -- the kinematic x floor on x_Bj, the first
# bin carrying the rate elsewhere.
XLIM = {
    "theta": {"results": (None, 3e-2), "results_nu": (None, 8e-2)},
    "Q2": {"results": (None, 100.0), "results_nu": (None, 400.0)},
    "xbj": {"results": (None, 0.6), "results_nu": (None, 0.6)},
}
REBIN = {("Emu", "results_nu"): 5, ("Emu", "results"): 5,
         ("Q2", "results_nu"): 2,
         ("xbj", "results"): 2, ("xbj", "results_nu"): 2}


def _rb(obs, resdir):
    return REBIN.get((obs, resdir))


def _snap(edges, v, upper):
    e = np.asarray(edges, dtype=float)
    if v is None:
        return float(e[-1] if upper else e[0])
    return float(e[int(np.argmin(np.abs(e - v)))])


def _bulk(edges, y):
    """Bins carrying at least 0.1% of the distribution."""
    w = y * np.diff(edges)
    return w > 0.001 * w.sum()


def _x_kin_min():
    """Smallest reachable x_Bj, Q2min/(y_max 2k.P), with y_max = 1 here."""
    return Q2MIN / beams.Beams("mu", ENERGY).two_kP


def _xl(obs, resdir):
    """(lo, hi) of the drawn window, as pp04 computes it."""
    lo, hi = XLIM.get(obs, {}).get(resdir, (None, None))
    den = _raw_hist(resdir, DENOM[resdir], obs)
    if den is None:
        return lo, hi
    e, y = den[0], den[1]
    if obs == "xbj":
        want = max(lo if lo is not None else 0.0, _x_kin_min())
        above = e[e >= want - 1e-12]
        lo = float(above[0]) if above.size else float(e[0])
    elif lo is None:
        nz = np.nonzero(_bulk(e, y))[0]
        if nz.size:
            lo = float(e[nz[0]])
    if hi is not None:
        hi = _snap(e, hi, True)
    return lo, hi


def _suffix():
    """The selection and energy tags, as the result filenames carry them."""
    esuf = ("" if ENERGY == beams.ANCHOR_ENERGY
            else f"_{beams.Beams('mu', ENERGY).tag}")
    # THE INCLUSIVE REGION CARRIES NO SELECTION SUFFIX: histos_powheg.json,
    # not histos_powheg_inclusive.json.  Naming it anyway would look for a
    # file that has never existed.
    return f"_{SELECTION}{TSUF}{esuf}"


def _path(resdir, key):
    return f"{BASE}/{resdir}/histos_{key}{_suffix()}.json"


def _merge(edges, dsig, err, n):
    """Merge n adjacent bins, keeping dsig a density and err its error.

    Errors add in quadrature WEIGHTED BY BIN WIDTH, because what adds is the
    integral over each bin and not the density.  These axes happen to be
    uniform, where the two are the same; it is written correctly rather than
    correctly by accident, as paper plot 2's version is.
    """
    m = (len(dsig) // n) * n
    w = np.diff(edges)[:m].reshape(-1, n)
    y = (dsig[:m].reshape(-1, n) * w).sum(axis=1)
    e = np.sqrt(((err[:m].reshape(-1, n) * w) ** 2).sum(axis=1))
    wide = w.sum(axis=1)
    return edges[:m + 1:n], y / wide, e / wide


def _raw_hist(resdir, key, obs):
    return hist(resdir, key, obs, rebin=False)


def hist(resdir, key, obs, rebin=True):
    """(edges, dsig, err) for one observable, or None if absent.

    THE FILE'S OWN STAMP IS CHECKED, both the selection and the Q2 floor.  A
    result from the inclusive region or from a raised floor would draw a
    perfectly plausible curve in a panel titled "inclusive" -- the silent
    failure CONVENTIONS.md rule 2 is about -- so it is an error, not a warning.
    """
    p = _path(resdir, key)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    sel = d.get("selection") or {}
    if sel.get("name") != SELECTION:
        raise SystemExit(f"{p}: selection is {sel.get('name')!r}, not "
                         f"{SELECTION!r} -- this figure would be mislabelled")
    if sel.get("q2_min") not in (None, Q2MIN):
        raise SystemExit(f"{p}: Q2 floor is {sel.get('q2_min')}, not {Q2MIN}")
    h = (d.get("hists") or {}).get(obs)
    if not h:
        return None
    e = np.asarray(h["edges"], dtype=float)
    y = np.asarray(h["dsig"], dtype=float)
    r = np.asarray(h["err"], dtype=float)
    n = _rb(obs, resdir) if rebin else None
    if n:
        # >>> TRIM TO THE WINDOW BEFORE MERGING, NOT AFTER. <<<  The lepton
        # grid starts at 100 GeV in 20 GeV steps and the selection begins at 200, so
        # merging three-to-one from the grid's own start puts the boundary
        # inside a bin: the 160-220 bin holds only its top fifth and draws as
        # a spurious dip at the threshold -- the sharp edge that is the point
        # of that row rendered as a feature.  Dropping the bins below the
        # window first makes the merged edges start ON it.
        lo = XLIM.get(obs, {}).get(resdir, (None, None))[0]
        if lo is not None:
            k = int(np.searchsorted(e, lo - 1e-9))
            e, y, r = e[k:], y[k:], r[k:]
        e, y, r = _merge(e, y, r, n)
    return e, y, r


def _vis(edges, obs, resdir):
    """Boolean mask of the bins inside this observable's window."""
    lo, hi = _xl(obs, resdir)
    left, right = edges[:-1], edges[1:]
    keep = np.ones(len(left), dtype=bool)
    if lo is not None:
        keep &= left >= lo - 1e-9
    if hi is not None:
        keep &= right <= hi + 1e-9
    return keep


def block(ax, axr, resdir, rows, obs, div, unit, xs, ys, ylab_q):
    """One observable in one column: absolute panel above, ratio below.

    >>> BOTH PANELS ARE HISTOGRAMS (user, 2026-09-08: "in the bottom panels
    use histograms ... we should be consistent in our plotting style through
    the paper"). <<<  Every binned figure in this set draws a step and hangs
    the error bars off the bin centres -- `stairs` plus `errorbar(fmt="none")`
    -- and the ratio panel is drawn the same way as the absolute one above it.
    This figure opened with markers in the ratio panel, which read as points
    of a measurement rather than as the binned distribution they are.  Paper
    plot 2 is the one figure that keeps markers, and it has to: it is a
    cross-section against beam energy at three points, not a histogram.
    """
    den = hist(resdir, DENOM[resdir], obs)
    for label, key, colour, marker, ls in rows:
        h = hist(resdir, key, obs)
        if h is None:
            continue
        e, y, r = h
        mid = 0.5 * (e[:-1] + e[1:])
        ax.stairs(y / div, e, color=colour, lw=1.8, ls=ls, baseline=None,
                  zorder=5, label=tex(label))
        ax.errorbar(mid, y / div, yerr=r / div, fmt="none", ecolor=colour,
                    elinewidth=1.0, alpha=0.75, zorder=5)
        if den is None or key == DENOM[resdir]:
            continue
        de, dy, dr = den
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(dy > 0, y / dy, np.nan)
            rerr = np.where(dy > 0, r / np.where(dy > 0, dy, 1.0), np.nan)
        axr.stairs(ratio, e, color=colour, lw=1.8, ls=ls, baseline=None,
                   zorder=5)
        axr.errorbar(mid, ratio, yerr=rerr, fmt="none", ecolor=colour,
                     elinewidth=1.0, alpha=0.75, zorder=5)
    if den is not None:
        de, dy, dr = den
        with np.errstate(divide="ignore", invalid="ignore"):
            band = np.where(dy > 0, dr / np.where(dy > 0, dy, 1.0), np.nan)
        axr.stairs(1.0 + band, de, baseline=1.0 - band, fill=True,
                   color="#0072b2", alpha=0.18, lw=0, zorder=1,
                   label=tex("stat. error"))
        axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-", zorder=1)
        mhou(ax, axr, resdir, obs, de, dy, div)
    ax.set_xscale(xs)
    ax.set_yscale(ys)
    axr.set_xscale(xs)
    lo, hi = _xl(obs, resdir)
    for a in (ax, axr):
        if lo is not None:
            a.set_xlim(left=lo)
        if hi is not None:
            a.set_xlim(right=hi)
    # THE LOG AXIS IS SET BY THE CURVES, not by the scale band's fill,
    # whose lower edge can reach zero in an unresolved tail bin and dragged
    # the muon Q2 panel down to 1e-18.
    if ys == "log" and den is not None:
        de, dy, _dr = den
        keep = (dy > 0) & _vis(de, obs, resdir)
        if keep.any():
            ax.set_ylim(float(dy[keep].min()) / div * 0.3,
                        float(dy[keep].max()) / div * 3.0)
    if xs == "log":
        # no minor-tick labels: on a short log range they overprint
        axr.xaxis.set_minor_formatter(ticker.NullFormatter())
        # 1-2-5 ticks on a range under two decades, as on pp04
        x0, x1 = axr.get_xlim()
        if x1 / x0 < 100.0:
            t = [m * 10.0 ** k for k in range(-6, 7) for m in (1, 2, 5)
                 if x0 <= m * 10.0 ** k <= x1]
            axr.set_xticks(t)
            axr.set_xticklabels([f"{v:g}" for v in t])
            axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    if ys == "log":
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_ylabel(tex(rf"${ylab_q}$  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL - 3)
    # PRUNE THE ABSOLUTE PANEL'S LOWEST TICK LABEL.  The two panels of a block
    # are fused, so that label sits on the boundary and collides with the
    # ratio panel's y-axis label -- which is exactly where the eye goes.
    # Only on a linear axis: on a log one the locator is not this one and
    # pruning would drop a decade.
    if ys == "linear":
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune="lower"))
    plotstyle.ticks(ax, labelbottom=False)
    plotstyle.ticks(axr)
    return None


def mhou(ax, axr, resdir, obs, de, dy, div):
    """The 7-point scale band of the column's own no-QED reference.

    >>> A QED CORRECTION IS ONLY INTERESTING BESIDE THE UNCERTAINTY OF THE
    CALCULATION IT CORRECTS (user, 2026-09-08: "can you also please show
    MHOUs in the POWHEG-prediction?  After all, an important aspect of NLO
    calculations is a proper estimate of MHOUs, so we should not forget about
    this"). <<<  The whole message of this figure is a per-cent-level shift,
    and whether a per cent matters is a question the missing higher orders
    answer.

    Computed by analysis/mhou_hadron.py, which joins POWHEG's seven scale
    weights to the SHOWERED events by `lhe_index`.

    BOTH COLUMNS' BANDS ARE MEASURED ON POWHEG-V2 AND THE LEGEND SAYS SO.
    The curves here are POWHEG-RES on both currents, and POWHEG-RES runs with
    `manyseeds` over a hundred-odd Les Houches files, which cannot be
    reweighted in one pass.  The band is the scale uncertainty of the same
    next-to-leading-order calculation of the same process, measured in the
    code that can be reweighted.
    """
    b = mhou_band.rel(resdir, obs, SELECTION, esuffix=TSUF,
                      trim_lo=(XLIM.get(obs, {}).get(resdir, (None, None))[0]
                               if _rb(obs, resdir) else None),
                      rebin=_rb(obs, resdir))
    if b is None:
        return
    e, lo, hi = b
    if not np.all([np.any(np.isclose(x, de)) for x in e]):
        raise SystemExit(f"{obs}: the scale band's edges are not a subset of "
                         f"the curve's -- refusing to draw it")
    # >>> AN OUTLINE, NOT A HATCH, ON THIS FIGURE. <<<  Paper plot 8 hatches
    # the band because there it is small beside the curves.  Here it is
    # SEVERAL TIMES the effect the figure exists to show, so a filled or
    # hatched band covers the whole panel and buries the QED arms under the
    # uncertainty they are being compared with -- and a sparse hatch is worse
    # still, since its long diagonals read as curves.  Two dotted envelope
    # lines say the same thing and compete with nothing.
    lab = tex(band_label(resdir))
    axr.stairs(1.0 + hi, e, baseline=None, color="#7a7a7a", lw=1.1,
               ls=(0, (2.0, 1.6)), zorder=2, label=lab)
    axr.stairs(1.0 + lo, e, baseline=None, color="#7a7a7a", lw=1.1,
               ls=(0, (2.0, 1.6)), zorder=2)
    # AND ONLY IN THE RATIO PANEL.  On paper plot 8 the band is also drawn
    # around the reference above, where it is a thin ribbon.  Here it is wider
    # than the spacing between the four arms, so on the absolute panel it
    # would cross every one of them and read as two more curves; the ratio
    # panel is where a per-cent band can be seen at all.


class _RefBandHandler(HandlerBase):
    """Legend key for (reference line, band envelope): the reference line in
    the middle, one dotted envelope line above it and one below, as the band
    is drawn in the ratio panels."""

    def create_artists(self, legend, orig, xd, yd, w, h, fs, trans):
        ref, env = orig
        out = []
        for y, src in ((0.5 * h - yd, ref), (0.95 * h - yd, env),
                       (0.05 * h - yd, env)):
            ln = Line2D([-xd, w - xd], [y, y], color=src.get_edgecolor()
                        if hasattr(src, "get_edgecolor") else src.get_color())
            ln.set_linewidth(src.get_linewidth())
            ln.set_linestyle(src.get_linestyle())
            ln.set_transform(trans)
            out.append(ln)
        return out


def band_label(resdir):
    """The band is POWHEG's own on each current (no cross-variant)."""
    # ONE LEGEND ENTRY FOR BOTH COLUMNS (user, 2026-09-18, as on paper
    # plots 8 and 9): each band sits in its own column, so the column already
    # says which current it is; the label is the paper's for the same band.
    return "NLO MHOU"


def ratio_range(resdir, obs):
    """A window that contains every drawn ratio point, padded, and clipped.

    Read from the data rather than typed, so a panel cannot silently crop the
    curve it exists to show.  The clip stops one empty tail bin, where a
    ratio can be anything at all, from setting the range for the rest.
    """
    den = hist(resdir, DENOM[resdir], obs)
    if den is None:
        return 0.5, 1.5
    de, dy, _ = den
    keep = _carries_rate(resdir, obs)
    if keep is None or not keep.any():
        keep = _vis(de, obs, resdir)
    vals = []
    rows = MU_ROWS if resdir == "results" else NU_ROWS
    for _l, key, _c, _m, _ls in rows:
        if key == DENOM[resdir]:
            continue
        h = hist(resdir, key, obs)
        if h is None:
            continue
        _e, y, _r = h
        with np.errstate(divide="ignore", invalid="ignore"):
            rat = np.where(dy > 0, y / dy, np.nan)
        vals += [v for v in rat[keep] if np.isfinite(v)]
    # AND THE SCALE BAND, which is part of what the panel shows.  Left out,
    # the range is set by the curves alone and the band is cropped by the
    # axis -- a theory uncertainty half drawn reads as a smaller one.
    b = mhou_band.rel(resdir, obs, SELECTION, esuffix=TSUF,
                      trim_lo=(XLIM.get(obs, {}).get(resdir, (None, None))[0]
                               if _rb(obs, resdir) else None),
                      rebin=_rb(obs, resdir))
    # >>> THE BAND MAY OPEN THE RANGE, BUT NOT SET IT. <<<  On this figure the
    # scale band is several times the QED shift, so a range large enough to
    # contain all of it flattens the arms into one line -- and the arms are
    # what the figure is for.  The window is therefore grown to hold the band
    # up to twice the span the curves need, and no further; the band is drawn
    # either way and its true size is quoted in the text.
    if b is not None and vals:
        lo_c, hi_c = min(vals), max(vals)
        span = max(hi_c - lo_c, 0.02)
        be, blo, bhi = b
        inside = _vis(be, obs, resdir)
        bl = [1.0 + v for v in blo[inside] if np.isfinite(v)]
        bh = [1.0 + v for v in bhi[inside] if np.isfinite(v)]
        if bl and bh:
            vals += [max(min(bl), lo_c - 0.5 * span),
                     min(max(bh), hi_c + 0.5 * span)]
    if not vals:
        return 0.5, 1.5
    lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
    lo, hi = max(lo, 0.0), min(hi, 4.0)
    pad = 0.12 * max(hi - lo, 0.05)
    return max(0.0, lo - pad), hi + pad


def _carries_rate(resdir, obs, frac=1e-2):
    """Bins holding more than `frac` of the peak density in the denominator.

    THE RANGE IS SET BY THE BINS THAT CARRY THE DISTRIBUTION, not by its last
    populated one.  A tail bin with a handful of events has a ratio that can
    be anything, and on the multiplicity row it was setting a range three
    times too wide for the bulk -- so the 10-20% differences the figure is
    about were a few pixels.  The threshold is on the DENOMINATOR alone, so
    the same bins are used for every curve and no generator can widen the
    window by having a longer tail than the others.
    """
    den = hist(resdir, DENOM[resdir], obs)
    if den is None:
        return None
    de, dy, _ = den
    peak = float(np.nanmax(dy)) if len(dy) else 0.0
    return _vis(de, obs, resdir) & (dy > frac * peak)


# PANEL TITLES ENLARGED to FS_PANEL_TITLE + 4 (user, 2026-10-04), on all
# four Sect. 6 figures alike.
def column_title(name):
    # THE PANEL TITLE NAMES THE CUT AND NOT THE SELECTION: the selection is
    # in the figure title, where it is said once (user, 2026-09-08 -- "replace
    # FASER Tier E by FASER$\nu$ Selection, which is more accurate").  Spelling
    # it out in both column titles as well would take two lines each.
    # the region is NOT repeated in the panel title (user, 2026-09-18:
    # "keep things symmetric with Figure 4.2"); the caption states it
    return name


def main():
    fig = plt.figure(figsize=(11.4, 15.8))
    # Four observable blocks, each an absolute panel FUSED to its ratio panel
    # (hspace = 0 inside a block, plotstyle's house rule), with real space
    # between blocks so the rows do not read as one eight-panel grid.  Same
    # geometry as paper plot 2, by user instruction ("use the same 4*2 layout
    # of pp_03"), so the two figures read as a pair.
    outer = fig.add_gridspec(4, 2, hspace=0.30, wspace=0.22,
                             top=0.878, bottom=0.045, left=0.085, right=0.985)
    for i, (obs, xlab, xs, ys, ylab_q) in enumerate(OBS):
        for j, (resdir, rows, title, div, unit) in enumerate(COLUMNS):
            inner = outer[i, j].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                            hspace=0.0)
            ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            block(ax, axr, resdir, rows, obs, div, unit, xs, ys, ylab_q)
            if i == 0:
                ax.set_title(tex(column_title(title)),
                             fontsize=plotstyle.FS_PANEL_TITLE + 4)
            axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 2)
            # THE RATIO LABEL NAMES THE COLUMN'S OWN DENOMINATOR, on BOTH
            # columns.  A single "ratio to POWHEG" on the left would stand for
            # two different codes -- POWHEG-RES here, POWHEG-V2 there -- which
            # is exactly the ambiguity CONVENTIONS.md rule 3 forbids in prose and
            # which is no better on an axis.
            # >>> AND IT READS THE DENOMINATOR, NOT THE COLUMN. <<<  Both
            # columns of this figure re-shower POWHEG-RES events -- that is
            # the whole design, so the QED arms differ from their reference
            # in QED alone -- so the neutrino panel's denominator is
            # POWHEG-RES too.  It said POWHEG-V2, which is the benchmark's
            # neutrino row but not this figure's reference.
            axr.set_ylabel(tex("ratio to") + "\n"
                           + tex(DENOM_NAME[DENOM[resdir]]),
                           fontsize=plotstyle.FS_YLABEL - 4)
            # EACH PANEL SETS ITS OWN RATIO RANGE.  The two columns are
            # different currents and the four rows are different observables;
            # one shared range would be set by whichever panel spreads most
            # and would flatten the rest to a line.
            axr.set_ylim(*ratio_range(resdir, obs))
    # ONE LEGEND, ABOVE THE FIGURE, assembled from the top row and
    # de-duplicated by label -- the arrangement paper plot 2 settled on after
    # the two-legend version was found to duplicate five of seven entries.
    h, lab = [], []
    for a in fig.axes[:4]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    # References, then the QED variants in order, then the two bands: the
    # three columns read (no QED mu, no QED nu, FSR) (FSR + ISR, full)
    # (stat. error, NLO MHOU).  "QED full" used to fall off the end of this
    # list and land between the bands.
    want = [tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0]),
            tex(MU_ROWS[1][0]), tex(MU_ROWS[2][0]), tex(MU_ROWS[3][0]),
            tex("stat. error"), tex(band_label("results"))]
    want = list(dict.fromkeys(want))
    order = [lab.index(w) for w in want if w in lab]
    order += [k for k in range(len(lab)) if k not in order]
    h, lab = [h[k] for k in order], [lab[k] for k in order]
    # THE BAND AND ITS CENTRAL VALUE ARE ONE LEGEND ENTRY (user, 2026-09-30:
    # "integrate the band MHOU and the central value ... for POWHEG in the
    # legend").  The band is the scale uncertainty of each column's no-QED
    # reference, so each reference entry is drawn as its line between the
    # two dotted envelope lines, and the separate band entry goes.
    mlab = tex(band_label("results"))
    if mlab in lab:
        k = lab.index(mlab)
        env = h.pop(k)
        lab.pop(k)
        for raw in (MU_ROWS[0][0], NU_ROWS[0][0]):
            if tex(raw) in lab:
                j = lab.index(tex(raw))
                h[j] = (h[j], env)
                lab[j] = tex(raw + plotstyle.MHOU_SUFFIX)
    if h:
        # THE LEGEND NAMES THE GENERATORS AND IS READ FIRST, so it is set at
        # the house size rather than shrunk to fit (user, 2026-09-08:
        # "increase the font of the legend indicating the various MC
        # generators").  Three columns rather than five: at this size five
        # would run past the figure, and two rows of three is easier to scan
        # than one long strip.
        fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.958),
                   handler_map={tuple: _RefBandHandler()},
                   ncol=3, frameon=True, handlelength=2.4,
                   columnspacing=1.6, labelspacing=0.45,
                   fontsize=plotstyle.FS_LEGEND + 2)
    fig.suptitle(tex(FIG_TITLE), y=0.99, fontsize=plotstyle.FS_SUPTITLE)
    # NO PROSE INSIDE THE FIGURE (user rule, 2026-08-29): the explanation
    # lives in MESSAGE and in the caption.
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ------------------------------------------------------------------ claims --
# Executable, and they re-derive themselves from the same JSONs the figure
# plots.  A claim comparing a hand-typed number to itself checks nothing.
def _sel(resdir):
    p = _path(resdir, DENOM[resdir])
    with open(p) as f:
        return (json.load(f).get("selection") or {})


def _mean(resdir, key, obs):
    h = hist(resdir, key, obs)
    if h is None:
        return None
    e, y, _r = h
    keep = _vis(e, obs, resdir)
    w = np.diff(e)[keep]
    mid = (0.5 * (e[:-1] + e[1:]))[keep]
    n = (y[keep] * w).sum()
    return float((y[keep] * w * mid).sum() / n) if n else None


def _integral(resdir, key, obs, lo=None, hi=None, raw=False):
    """The rate in a range.  `raw` ignores the DRAWN window.

    A claim usually has to stay inside what the figure shows.  The
    leading-hadron panel is the exception: its window starts at 10 GeV by
    request, and the softest bins -- where the models differ most -- are
    measured but not drawn.  `raw=True` says so at the call site instead of
    silently returning zero, which is what the visible mask would do.
    """
    h = hist(resdir, key, obs)
    if h is None:
        return None
    e, y, _r = h
    keep = np.ones(len(y), dtype=bool) if raw else _vis(e, obs, resdir)
    if lo is not None:
        keep &= e[:-1] >= lo - 1e-9
    if hi is not None:
        keep &= e[1:] <= hi + 1e-9
    return float((y[keep] * np.diff(e)[keep]).sum())


def _tot(resdir, key):
    """Integrated fiducial cross-section of an arm, from its own result."""
    with open(_path(resdir, key)) as f:
        return json.load(f)["sigma_fid_pb"]


def _eff(resdir, key):
    return _tot(resdir, key) / _tot(resdir, DENOM[resdir]) - 1.0


def _edge_ratio(resdir, key, obs, lo=None, hi=None):
    a, b = _integral(resdir, key, obs, lo, hi), _integral(resdir, DENOM[resdir], obs, lo, hi)
    return a / b if b else None


def _qed_vs_mhou(r, obs):
    """(max |full QED - 1|, bins outside the NLO MHOU, bins) over the drawn bins."""
    rows = MU_ROWS if r == "results" else NU_ROWS
    full = [k for _l, k, *_ in rows if "full" in k][0]
    de, dy, _ = hist(r, DENOM[r], obs)
    e, y, _ = hist(r, full, obs)
    vis = _vis(de, obs, r) & (dy > 0)
    q = y[vis] / dy[vis] - 1.0
    be, lo, hi = mhou_band.rel(
        r, obs, SELECTION, esuffix=TSUF,
        trim_lo=(XLIM.get(obs, {}).get(r, (None, None))[0]
                 if _rb(obs, r) else None),
        rebin=_rb(obs, r))
    mid = 0.5 * (de[:-1] + de[1:])
    i = np.clip(np.searchsorted(be, mid) - 1, 0, len(lo) - 1)
    L, H = lo[i][vis], hi[i][vis]
    out = ((q > 0) & (q > H)) | ((q < 0) & (q < L))
    return float(np.abs(q).max()), int(out.sum()), len(q)


CLAIMS = [
    # THE PAPER'S QED-vs-MHOU SENTENCE (2026-10-04): inside the band for every
    # muon distribution and for nu theta and Q2; outside it in most nu x_Bj
    # and E_l' bins, ~10% at the smallest x and ~4% at the largest E_l'.
    {"what": "full QED stays inside the NLO MHOU for every muon distribution "
             "and for nu theta, Q2, but exceeds it in most nu x_Bj and E_l' "
             "bins, reaching ~10% (smallest x) and ~4% (largest E_l')",
     "check": lambda: all(_qed_vs_mhou("results", o)[1] == 0 for o, *_ in OBS)
     and all(_qed_vs_mhou("results_nu", o)[1] == 0 for o in ("theta", "Q2"))
     and all(_qed_vs_mhou("results_nu", o)[1] > _qed_vs_mhou("results_nu", o)[2] / 2
             for o in ("xbj", "Emu"))
     and 0.09 < _qed_vs_mhou("results_nu", "xbj")[0] < 0.12
     and 0.035 < _qed_vs_mhou("results_nu", "Emu")[0] < 0.05,
     "detail": lambda: "; ".join("%s %s max %.1f%% out %d/%d" % (
         (r, o, 100 * _qed_vs_mhou(r, o)[0]) + _qed_vs_mhou(r, o)[1:])
         for r in ("results", "results_nu") for o, *_ in OBS)},
    {"what": "the 7-point scale band of the matched calculation is drawn on "
             "both columns",
     "check": lambda: all(
         (mhou_band.rel(r, "Emu", SELECTION, esuffix=TSUF) is not None
          and float(np.max(np.abs(mhou_band.rel(r, "Emu", SELECTION, esuffix=TSUF)[2])))
          > 0.0) for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s %s" % (r, "absent" if mhou_band.rel(r, "Emu", SELECTION, esuffix=TSUF) is None
                    else "present") for r in ("results", "results_nu"))},
    {"what": "the arms are the SAME events -- each column re-showers one set "
             "of Les Houches files, so no arm moves the integrated rate by "
             "more than 0.2%",
     "check": lambda: all(abs(_eff(r, k)) < 0.002
                          for r in ("results", "results_nu")
                          for _l, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS)),
     "detail": lambda: "; ".join(
         "%s %s %+.3f%%" % (r, lab.split(" (")[0], 100 * _eff(r, k))
         for r in ("results", "results_nu")
         for lab, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS))},
    {"what": "integrated over Q2 > 4, W > 3 GeV the full QED effect is "
             "-0.07% on the charged current and 0.01% on the neutral one: "
             "without a y window nothing is lost at a y edge",
     "check": lambda: (-0.0015 < _eff("results_nu", "powheg_nu_qedfull") < 0.0
                       and abs(_eff("results", "powheg_qedfull")) < 0.0005
                       and abs(_eff("results_nu", "powheg_nu_qedfull"))
                       > abs(_eff("results", "powheg_qedfull"))),
     "detail": lambda: "mu %+.4f%%, nu %+.4f%%" % (
         100 * _eff("results", "powheg_qedfull"),
         100 * _eff("results_nu", "powheg_nu_qedfull"))},
    {"what": "differentially the charged-current lepton is SOFTENED: the "
             "lowest scattered-lepton energies gain 2-4% and the top 100 GeV "
             "loses 3-6%, while the neutral-current spectrum moves by under 1%",
     "check": lambda: (1.015 < _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", hi=100.0) < 1.05
                       and 0.93 < _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", lo=900.0) < 0.975
                       and abs(_edge_ratio("results", "powheg_qedfull", "Emu", lo=100.0, hi=900.0) - 1) < 0.01),
     "detail": lambda: "nu E'<100 %.4f, E'>900 %.4f; mu 100-900 %.4f" % (
         _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", hi=100.0),
         _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", lo=900.0),
         _edge_ratio("results", "powheg_qedfull", "Emu", lo=100.0, hi=900.0))},
    {"what": "initial-state radiation off the lepton line is a muon-only arm, "
             "and the neutrino side has no such result at all",
     "check": lambda: (len(MU_ROWS) == 4 and len(NU_ROWS) == 3
                       and not os.path.exists(
                           _path("results_nu", "powheg_nu_qedfsrisr"))),
     "detail": lambda: "muon arms %d, neutrino arms %d"
                       % (len(MU_ROWS), len(NU_ROWS))},
]

MESSAGE = (
    "<b>Every curve in a column is the same event sample.</b> The arms "
    "re-shower the same Les Houches files -- POWHEG-RES on the muon side, "
    "POWHEG-V2 on the neutrino side -- with QED radiation switched on in "
    "stages, so the matrix element, the matching, the parton distribution and "
    "the scale are identical by construction and what is drawn is the "
    "radiation alone. The no-QED re-shower reproduces the production sample "
    "byte for byte on both currents."
    "\n\n"
    "<b>On the integrated rate QED is now negligible.</b> In the region "
    "Q&sup2; &gt; 4 GeV&sup2;, W &gt; 3 GeV the full effect is &minus;0.07% on "
    "the charged current and +0.01% on the neutral one. With a y window the "
    "same radiation moved events across the y edge and the charged-current "
    "rate changed by half a per cent; without one the reconstructed "
    "kinematics migrate inside the region rather than out of it."
    "\n\n"
    "<b>Differentially the charged current still moves, and the neutral one "
    "does not.</b> Final-state radiation softens the scattered muon of the "
    "charged current: the lowest lepton energies gain a few per cent, the "
    "top of the spectrum loses 4&ndash;5%, and the reconstructed x and "
    "Q&sup2; shift with it. The neutral current has a charged lepton on both "
    "legs, the dipole between them largely cancels, and its spectra move by "
    "under a per cent."
    "\n\n"
    "<b>Initial-state radiation off the lepton line is a muon-only arm.</b> "
    "It needs a charged incoming lepton and a neutrino beam has none, so the "
    "neutrino column carries two arms where the muon column carries three."
)


if __name__ == "__main__":
    main()
