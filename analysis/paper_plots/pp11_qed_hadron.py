#!/usr/bin/env python3
"""Paper plot 11: QED corrections on the hadron-level observables, in the FINAL region.

The redo of the earlier figure for the "Paper plots" tab: the
same events re-showered with QED radiation switched on in stages
(powheg/qed/: fsr, fsrisr on the muon side only, full), each arm against the
sample with no QED.  The conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184, each arm run on p and n.
  * THE REGION is q4w3_faser_e: Q2 > 4 GeV2, W > 3 GeV, no y window, plus Tier E.
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

Usage: analysis/paper_plots/pp11_qed_hadron.py
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

SLUG = "pp11_qed_hadron"
V1_NUMBER = 11
IN_PAPER = True
TITLE = "QED corrections to the hadron-level observables"
FIG_TITLE = (r"QED corrections, hadron level, FASER$\nu$ selection, "
             r"$E_{\ell} = 1$ TeV")
OUTPUT = "pp_qed_hadron.png"
RESULTS = "results_nu"
CAPTION = ("Charged-hadron multiplicity inside the emulsion track "
           "acceptance, leading charged-hadron energy, the azimuth between "
           "the lepton and the summed hadrons, and the lepton momentum, per nucleon on tungsten at 1 TeV under the FASER&nu; selection in the region Q&sup2; &gt; 4 GeV&sup2;, W &gt; 3 GeV, for the same "
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
SELECTION = "q4w3_faser_e"
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
    ("nch05", r"$N_{\rm ch}$  ($E > 1$ GeV, $\tan\theta < 0.5$)",
     "linear", "linear", r"d\sigma/dN_{\rm ch}"),
    ("Elead", r"$E_{\rm lead}$  [GeV]", "log", "log",
     r"d\sigma/dE_{\rm lead}"),
    ("dphix", r"$\Delta\phi(\ell', \sum h^\pm)$  [rad]",
     "linear", "linear", r"d\sigma/d\Delta\phi"),
    ("Emu",   r"$p_{\ell'}$  [GeV]", "linear", "linear",
     r"d\sigma/dp_{\ell'}"),
]
# X WINDOWS.  Each is the region that carries the rate; anything unlisted
# keeps the full histogram.  The multiplicity axis stops where the samples
# run out of events rather than at the last bin edge, and the Delta phi axis
# starts at the cut, since the selection populates nothing below it.
# the window the quoted <N_ch> are averaged over (99.9% of the rate)
NCH_MEAN_WINDOW = (4.5, 19.5)

XLIM = {
    # FROM N_ch = 5 (user, 2026-10-04): the selection leaves the lower
    # bins empty by construction; 4.5 is the lower edge of the N = 5 bin.
    # DRAWN TO N_ch = 14 (user, 2026-10-04); the mean multiplicities the text
    # quotes are still taken over 5-19 (NCH_MEAN_WINDOW), so the cut is visual.
    "nch05": (4.5, 14.5),
    # FROM DELTA PHI = 2.5 (user, 2026-10-04: below it the MC errors are
    # too large for the prediction to be informative); 38 pi / 48 = 2.487 is
    # the stored edge nearest 2.5, and the 2:1 merge then gives five bins.
    "dphix": (38.0 * math.pi / 48.0, math.pi),
    "Elead": (10.0, 900.0),
    "Emu": (200.0, 1000.0),
}

REBIN = {"dphix": 2, "Emu": 3}


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


def hist(resdir, key, obs):
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
    n = REBIN.get(obs)
    if n:
        # >>> TRIM TO THE WINDOW BEFORE MERGING, NOT AFTER. <<<  The lepton
        # grid starts at 100 GeV in 20 GeV steps and the selection begins at 200, so
        # merging three-to-one from the grid's own start puts the boundary
        # inside a bin: the 160-220 bin holds only its top fifth and draws as
        # a spurious dip at the threshold -- the sharp edge that is the point
        # of that row rendered as a feature.  Dropping the bins below the
        # window first makes the merged edges start ON it.
        lo, _hi = XLIM.get(obs, (None, None))
        if lo is not None:
            k = int(np.searchsorted(e, lo - 1e-9))
            e, y, r = e[k:], y[k:], r[k:]
        e, y, r = _merge(e, y, r, n)
    return e, y, r


def _vis(edges, obs):
    """Boolean mask of the bins inside this observable's window."""
    lo, hi = XLIM.get(obs, (None, None))
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
    lo, hi = XLIM.get(obs, (None, None))
    for a in (ax, axr):
        if lo is not None:
            a.set_xlim(left=lo)
        if hi is not None:
            a.set_xlim(right=hi)
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
                      trim_lo=(XLIM.get(obs, (None, None))[0]
                               if REBIN.get(obs) else None),
                      rebin=REBIN.get(obs))
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
    # A GREY FILLED BAND, as on pp08/pp09 (user, 2026-10-04: the dotted
    # envelopes looked inconsistent with Figs. 6.1/6.2), in the ratio panel
    # and hugging the reference in the absolute one, beneath every curve.
    axr.stairs(1.0 + hi, e, baseline=1.0 + lo, fill=True, color="#111111",
               alpha=0.14, lw=0, zorder=0, label=tex(band_label(resdir)))
    with np.errstate(invalid="ignore"):
        i = np.searchsorted(e, 0.5 * (de[:-1] + de[1:])) - 1
        ok = (i >= 0) & (i < len(lo))
        f_lo = np.where(ok, lo[np.clip(i, 0, len(lo) - 1)], 0.0)
        f_hi = np.where(ok, hi[np.clip(i, 0, len(hi) - 1)], 0.0)
        ax.stairs(dy * (1.0 + f_hi) / div, de,
                  baseline=dy * (1.0 + f_lo) / div,
                  fill=True, color="#111111", alpha=0.14, lw=0, zorder=0)
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
        keep = _vis(de, obs)
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
                      trim_lo=(XLIM.get(obs, (None, None))[0]
                               if REBIN.get(obs) else None),
                      rebin=REBIN.get(obs))
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
        inside = _vis(be, obs)
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
    return _vis(de, obs) & (dy > frac * peak)


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
    # THE BAND AND ITS CENTRAL VALUE ARE ONE ENTRY (user, 2026-10-04), drawn
    # as on pp10: the reference line between the two dotted envelope lines.
    h, lab = plotstyle.merge_bands(h, lab, {
        tex(band_label("results")): [tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0])]})
    if h:
        # THE LEGEND NAMES THE GENERATORS AND IS READ FIRST, so it is set at
        # the house size rather than shrunk to fit (user, 2026-09-08:
        # "increase the font of the legend indicating the various MC
        # generators").  Three columns rather than five: at this size five
        # would run past the figure, and two rows of three is easier to scan
        # than one long strip.
        fig.legend(h, lab,                    loc="upper center", bbox_to_anchor=(0.535, 0.958),
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
    keep = _vis(e, obs)
    if obs == "nch05":
        keep = ((e[:-1] >= NCH_MEAN_WINDOW[0] - 1e-9)
                & (e[1:] <= NCH_MEAN_WINDOW[1] + 1e-9))
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
    keep = np.ones(len(y), dtype=bool) if raw else _vis(e, obs)
    if lo is not None:
        keep &= e[:-1] >= lo - 1e-9
    if hi is not None:
        keep &= e[1:] <= hi + 1e-9
    return float((y[keep] * np.diff(e)[keep]).sum())


def _tot(resdir, key):
    with open(_path(resdir, key)) as f:
        d = json.load(f)
    return d["sigma_fid_pb"], d.get("sigma_fid_err_pb") or 0.0


def _eff(resdir, key):
    return _tot(resdir, key)[0] / _tot(resdir, DENOM[resdir])[0] - 1.0


def _edge_ratio(resdir, key, obs, lo=None, hi=None):
    a, b = _integral(resdir, key, obs, lo, hi), _integral(resdir, DENOM[resdir], obs, lo, hi)
    return a / b if b else None


CLAIMS = [
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
             "of Les Houches files, so no arm moves the selected rate by more "
             "than 1%",
     "check": lambda: all(abs(_eff(r, k)) < 0.01
                          for r in ("results", "results_nu")
                          for _l, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS)),
     "detail": lambda: "; ".join(
         "%s %s %+.3f%%" % (r, lab.split(" (")[0], 100 * _eff(r, k))
         for r in ("results", "results_nu")
         for lab, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS))},
    {"what": "under the FASER&nu; selection QED lowers the charged-current "
             "rate by 0.5-1.0%, resolved against the statistical error, "
             "while on the neutral current no arm is resolved",
     "check": lambda: (-0.010 < _eff("results_nu", "powheg_nu_qedfull") < -0.005
                       and abs(_eff("results_nu", "powheg_nu_qedfull"))
                       > 3 * _tot("results_nu", DENOM["results_nu"])[1]
                       / _tot("results_nu", DENOM["results_nu"])[0]
                       and all(abs(_eff("results", k))
                               < 1.5 * _tot("results", DENOM["results"])[1]
                               / _tot("results", DENOM["results"])[0]
                               for _l, k, _c, _m, _ls in MU_ROWS[1:])),
     "detail": lambda: "nu full %+.3f%% (stat %.3f%%); mu %s (stat %.3f%%)" % (
         100 * _eff("results_nu", "powheg_nu_qedfull"),
         100 * _tot("results_nu", DENOM["results_nu"])[1] / _tot("results_nu", DENOM["results_nu"])[0],
         ", ".join("%+.3f%%" % (100 * _eff("results", k)) for _l, k, _c, _m, _ls in MU_ROWS[1:]),
         100 * _tot("results", DENOM["results"])[1] / _tot("results", DENOM["results"])[0])},
    {"what": "the charged-current lepton is softened: above 800 GeV the "
             "selected rate falls by 3-8%",
     "check": lambda: 0.92 < _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", lo=800.0) < 0.97,
     "detail": lambda: "nu p_lep > 800: %.4f" % _edge_ratio("results_nu", "powheg_nu_qedfull", "Emu", lo=800.0)},
    {"what": "initial-state radiation off the lepton line is a muon-only arm, "
             "and the neutrino side has no such result at all",
     "check": lambda: (len(MU_ROWS) == 4 and len(NU_ROWS) == 3
                       and not os.path.exists(
                           _path("results_nu", "powheg_nu_qedfsrisr"))),
     "detail": lambda: "muon arms %d, neutrino arms %d"
                       % (len(MU_ROWS), len(NU_ROWS))},
]

MESSAGE = (
    "<b>The same question as the previous figure, asked where FASER "
    "measures.</b> The arms re-shower the same Les Houches files -- POWHEG-RES "
    "on the muon side, POWHEG-V2 on the neutrino side -- with QED radiation "
    "switched on in stages, and the curves are the hadron-level observables "
    "under the FASER&nu; selection."
    "\n\n"
    "<b>On the charged current QED costs the selection just under a per "
    "cent.</b> The selection cuts on the lepton momentum, and final-state "
    "radiation softens it: the full arm lowers the selected rate by "
    "0.7%, several times its statistical error, and removes 4&ndash;6% of "
    "the events above 800 GeV. The hadronic observables do not move beyond "
    "that normalisation."
    "\n\n"
    "<b>On the neutral current no arm is resolved.</b> With a charged lepton "
    "on both legs the dipole largely cancels, and the shifts sit inside the "
    "statistical error of the sample, bin by bin as well as integrated."
    "\n\n"
    "<b>Initial-state radiation off the lepton line is a muon-only arm</b>, "
    "for the reason the previous figure gives."
)


if __name__ == "__main__":
    main()
