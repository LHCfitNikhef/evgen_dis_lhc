#!/usr/bin/env python3
"""Paper plot 8: hadron-level distributions under the FASER$\nu$
selection, in the FINAL region.

The redo of the earlier figure for the "Paper plots"
tab.  Everything is the earlier production's figure -- the same four observables in the
form the emulsion measures them (arXiv:2403.12520), the same generators, each
column normalised to its own matched calculation -- with the conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184, every generator run on the
    proton and the neutron separately (PRODUCTION.md).
  * THE SELECTION IS NESTED IN THE  REGION (user ruling, paper Sec. 3):
    Q2 > 4 GeV2 and W > 3 GeV with no y window, plus Tier E's cuts
    (selection q4w3_faser_e).  With no y < 0.9 the lepton reaches the beam
    energy, so the lepton row runs to 1 TeV.
  * THE SCALE BAND IS EACH COLUMN'S OWN.  earlier measured the muon band on the
    POWHEG-V2 cross-variant because POWHEG-RES could not be reweighted; 
    reweights POWHEG-RES seed by seed (powheg/production/reweight_scale.sh) and the
    band is joined to the showered events per job, per nucleon
    (analysis/mhou_hadron.py), so both bands sit on the curve they
    were measured on.

GENIE IS SHOWN, in the default tune FASER runs (rule 1b).

Usage: analysis/paper_plots/pp08_hadron_level_faser_e.py
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
import beams                                                # noqa: E402
import mhou_band                                            # noqa: E402

SLUG = "pp08_hadron_level_faser_e"
V1_NUMBER = 8
# IN THE PAPER in place of the earlier production's figure (user, 2026-09-14/15: every updated 
# plot goes into the paper).  Read by tools/check_paper_figures_used.py.
IN_PAPER = True
TITLE = "Hadron-level distributions under the FASER$\nu$ selection"
FIG_TITLE = (r"Hadron-level distributions, FASER$\nu$ selection, "
             r"$E_{\ell} = 1$ TeV")
OUTPUT = "pp_hadron_level_faser_e.png"
RESULTS = "results_nu"
CAPTION = ("Charged-hadron multiplicity above 1 GeV inside the emulsion "
           "track acceptance, leading charged-hadron energy, azimuthal "
           "separation between the lepton and the summed charged-hadron "
           "system, and outgoing-lepton momentum, per nucleon on tungsten "
           "at 1 TeV under the FASER&nu; selection in the region Q&sup2; "
           "&gt; 4 GeV&sup2;, W &gt; 3 GeV; muon neutral current (left) and "
           "neutrino charged "
           "current (right). POWHEG-RES on the muon side and POWHEG-V2 on "
           "the neutrino side, against GENIE in the default tune FASER runs "
           "Sherpa MC@NLO and Herwig. Beneath each panel the ratio to that column&rsquo;s "
           "own matched calculation &mdash; POWHEG-RES on the left, "
           "POWHEG-V2 on the right &mdash; with its statistical error and its "
           "seven-point scale band. The selection is that of arXiv:2403.12520: at least "
           "five charged hadrons "
           "with tan&thinsp;&theta; &lt; 0.5, at least four with "
           "tan&thinsp;&theta; &lt; 0.1, a lepton above 200 GeV at "
           "&theta; &gt; 5 mrad, and &Delta;&phi; &gt; &pi;/2, so the third "
           "row is populated above &pi;/2 only.")

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
    (r"POWHEG-RES ($\mu$)", "powheg",     "#0072b2", "o", "-"),
    ("GENIE (GRV98LO)",     "genie",      "#d55e00", "s", "--"),
    ("Sherpa",              "sherpa",     "#009e73", "^", "-."),
    # HERWIG ADDED (user, 2026-10-04: "surely we want to include them as
    # well, for completeness") -- it had never been run through Tier E.
    ("Herwig",              "herwig_nlo_powheg_full", "#7b3294", "D", ":"),
]
NU_ROWS = [
    (r"POWHEG-V2 ($\nu$)",  "powheg_nu",  "#0072b2", "o", POWHEG_V2_LS),
    ("GENIE (GRV98LO)",     "genie_lo",   "#d55e00", "s", "--"),
    ("Sherpa",              "sherpa_nlo", "#009e73", "^", "-."),
    ("Herwig",              "herwig_nlo_full", "#7b3294", "D", ":"),
]
# The denominator of each column's ratio panel: its own POWHEG.
DENOM = {"results": "powheg", "results_nu": "powheg_nu"}

COLUMNS = [("results", MU_ROWS, "Muon DIS", 1000.0, "nb"),
           ("results_nu", NU_ROWS, "Neutrino DIS", 1.0, "pb")]

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
    # 99.9% of the rate is below 20 tracks on both currents; the bins beyond
    # it hold a handful of events and would set the ratio range for the rest.
    # >>> THE MULTIPLICITY PANEL STARTS AT THREE (user, 2026-09-08). <<<
    # The selection needs five tracks including the lepton, so four charged
    # hadrons is the lowest multiplicity that can appear at all and the bins
    # below it hold only migration; drawing them puts three empty columns
    # where the eye starts reading.
    # FROM N_ch = 5 (user, 2026-10-04): the selection leaves the lower
    # bins empty by construction; 4.5 is the lower edge of the N = 5 bin.
    # DRAWN TO N_ch = 14 (user, 2026-10-04); the mean multiplicities the text
    # quotes are still taken over 5-19 (NCH_MEAN_WINDOW), so the cut is visual.
    "nch05": (4.5, 14.5),
    "dphix": (math.pi / 2.0, math.pi),
    # >>> THE LEADING-HADRON PANEL STARTS AT 10 GeV (user, 2026-09-08). <<<
    # Below it the samples have a handful of events per bin and the ratio
    # panel is all error bar; the region is still MEASURED -- the claims
    # below read it straight from the histogram with `raw=True` -- but it is
    # not drawn, so the panel shows the part that carries the rate.
    "Elead": (10.0, 900.0),
    # THE LEPTON PANEL STARTS AT THE CUT.  Tier E requires 200 GeV and the
    # histogram grid starts at 100, so without this the panel opens with a
    # decade of empty axis and the sharp edge at the threshold -- which is
    # the point of the row -- reads as a feature somewhere in the middle.
    "Emu": (200.0, 1000.0),
}

# REBINNING, applied to EVERY curve on a panel including the denominator, so
# no ratio is ever taken between differently binned things.  Only the two
# fine axes need it: dphix is 48 bins over the full range and Tier E fills 24
# of them, which is more than the eye can read in a panel this size.
REBIN = {"dphix": 2,
         # THE LEPTON ROW IS REBINNED THREE TO ONE (user, 2026-09-08: "the
         # p_lep distributions could do with some rebinning to tame MC
         # integration errors").  The stored grid is 20 GeV wide from 100 to
         # 800, which leaves 30 bins inside the Tier E window and a visible
         # bin-to-bin jitter on a row whose whole content is that the ratio is
         # FLAT -- scatter is exactly what obscures flatness.  60 GeV bins cut
         # the error by 1.7 and leave ten across the window.
         "Emu": 3}


def _suffix():
    """The selection and energy tags, as the result filenames carry them."""
    esuf = ("" if ENERGY == beams.ANCHOR_ENERGY
            else f"_{beams.Beams('mu', ENERGY).tag}")
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
    perfectly plausible curve in a panel titled Tier E -- the silent failure
    CONVENTIONS.md rule 2 is about -- so it is an error rather than a warning.
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
    # THE MULTIPLICITY GRID IS INTEGER, AND A FILE ON THE OLD ONE IS REFUSED.
    # nch05 was first histogrammed on the width-2 grid the other multiplicity
    # observables use; a result written before 2026-09-08 therefore has half
    # as many bins, and drawing it beside a re-parsed one would put two
    # different binnings on the same panel and take a ratio between them.
    if obs == "nch05" and abs(float(np.diff(e)[0]) - 1.0) > 1e-9:
        raise SystemExit(f"{p}: nch05 is binned in steps of "
                         f"{np.diff(e)[0]:g}, not 1 -- re-run "
                         f"tools/tier_e_hadron_passes.sh")
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
    # A LOG PANEL IS SCALED ON THE DENOMINATOR, not on whichever curve has
    # the emptiest last bin: with Herwig added (2026-10-04) one sparse muon
    # E_lead tail bin pulled the axis down to 1e-19.
    if ys == "log" and len(dy):
        pos = dy[dy > 0] / div
        if pos.size:
            ax.set_ylim(0.3 * float(pos.min()), 5.0 * float(pos.max()))
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
    """The 7-point scale band of the column's own matched calculation.

    >>> A MATCHED NLO PREDICTION WITHOUT ITS SCALE BAND IS HALF A PREDICTION
    (user, 2026-09-08: "an important aspect of NLO calculations is a proper
    estimate of MHOUs, so we should not forget about this"). <<<  The band is
    computed by analysis/mhou_hadron.py, which joins POWHEG's seven scale
    weights to the SHOWERED events through `lhe_index` -- the only route that
    survives a selection cutting on the charged-hadron multiplicity.

    It is drawn twice: hugging the POWHEG curve in the absolute panel, and
    around unity in the ratio panel, where it is what the other generators
    are to be read against.

    THE MUON COLUMN'S BAND IS MEASURED ON POWHEG-V2 AND SAYS SO IN THE
    LEGEND.  POWHEG-RES runs with `manyseeds` over a hundred-odd Les Houches
    files and cannot be reweighted in one pass, so the muon band comes from
    the same code as the neutrino one, run on the neutral current.  Drawing
    it unlabelled would put a POWHEG-V2 number on a POWHEG-RES curve without
    saying so.
    """
    b = mhou_band.rel(resdir, obs, SELECTION, esuffix=TSUF,
                      trim_lo=(XLIM.get(obs, (None, None))[0]
                               if REBIN.get(obs) else None),
                      rebin=REBIN.get(obs))
    if b is None:
        return
    e, lo, hi = b
    # THE BAND'S EDGES MUST BE A SUBSET OF THE CURVE'S.  mhou_band merges
    # adjacent bins where the envelope is not resolved above its own Monte
    # Carlo error, so the band is coarser than the curve -- but never on a
    # different grid, and a band on a different grid would misplace the
    # uncertainty of one bin onto part of another.
    if not np.all([np.any(np.isclose(x, de)) for x in e]):
        raise SystemExit(f"{obs}: the scale band's edges are not a subset of "
                         f"the curve's -- refusing to draw it")
    # A GREY, UNHATCHED BAND (user, 2026-09-18), as the MHOU bands are drawn
    # in every other figure of the paper.
    axr.stairs(1.0 + hi, e, baseline=1.0 + lo, fill=True, color="#111111",
               alpha=0.14, lw=0, zorder=2, label=tex(band_label(resdir)))
    with np.errstate(invalid="ignore"):
        # and on the absolute panel, hugging the curve.  The band is coarse
        # and the curve is not, so the relative band is spread back onto the
        # curve's own bins before it is drawn.
        i = np.searchsorted(e, 0.5 * (de[:-1] + de[1:])) - 1
        ok = (i >= 0) & (i < len(lo))
        f_lo = np.where(ok, lo[np.clip(i, 0, len(lo) - 1)], 0.0)
        f_hi = np.where(ok, hi[np.clip(i, 0, len(hi) - 1)], 0.0)
        ax.stairs(dy * (1.0 + f_hi) / div, de,
                  baseline=dy * (1.0 + f_lo) / div,
                  fill=True, color="#111111", alpha=0.14, lw=0, zorder=2)


def band_label(resdir):
    """Each band is measured on its own column's matched calculation, so
    the current alone names it (no cross-variant caveat, unlike the earlier production)."""
    # ONE LEGEND ENTRY FOR BOTH COLUMNS (user, 2026-09-18: "one grey band
    # for MHOU is sufficient").  The two bands are drawn alike and each sits
    # in its own column, so the column already says which current it is; the
    # label is the one every other paper figure uses for the same band.
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
    if b is not None:
        be, blo, bhi = b
        inside = _vis(be, obs)
        vals += [1.0 + v for v in blo[inside] if np.isfinite(v)]
        vals += [1.0 + v for v in bhi[inside] if np.isfinite(v)]
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
            axr.set_ylabel(tex("ratio to") + "\n"
                           + tex("POWHEG-RES" if resdir == "results"
                                 else "POWHEG-V2"),
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
    want = [tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0]),
            tex(MU_ROWS[1][0]), tex(MU_ROWS[2][0]), tex(MU_ROWS[3][0]),
            tex("stat. error"),
            tex(band_label("results_nu")),
            tex(band_label("results"))]
    want = list(dict.fromkeys(want))
    order = [lab.index(w) for w in want if w in lab]
    order += [k for k in range(len(lab)) if k not in order]
    h, lab = [h[k] for k in order], [lab[k] for k in order]
    # the band shares POWHEG's entry, on each current (user, 2026-10-04)
    h, lab = plotstyle.merge_bands(h, lab, {
        tex(band_label("results")): [tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0])]})
    if h:
        # THE LEGEND NAMES THE GENERATORS AND IS READ FIRST, so it is set at
        # the house size rather than shrunk to fit (user, 2026-09-08:
        # "increase the font of the legend indicating the various MC
        # generators").  Three columns rather than five: at this size five
        # would run past the figure, and two rows of three is easier to scan
        # than one long strip.
        fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.958),
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


def _dphix_below_cut(resdir):
    """Rate below pi/2, which the selection forbids."""
    h = hist(resdir, DENOM[resdir], "dphix")
    if h is None:
        return None
    e, y, _r = h
    below = e[1:] <= math.pi / 2.0 + 1e-9
    return float((y[below] * np.diff(e)[below]).sum())


def _mult_spread(resdir):
    ms = [_mean(resdir, k, "nch05")
          for _l, k, _c, _m, _ls in (MU_ROWS if resdir == "results"
                                     else NU_ROWS)]
    ms = [m for m in ms if m]
    return (max(ms) / min(ms)) if len(ms) > 1 else None



def _mhou(resdir):
    """(hi, lo) of the integrated 7-point band, as fractions."""
    with open(mhou_band.path(resdir, SELECTION, TSUF)) as f:
        d = json.load(f)
    return d["sigma_rel_hi"], d["sigma_rel_lo"]


CLAIMS = [
    # AND THE SIZE OF THE BAND, which the message quotes: the muon neutral
    # current's scale uncertainty is five times the neutrino charged
    # current's, and that ratio is the message, not the two numbers alone.
    {"what": "integrated over the selection the 7-point band is "
             "+4.3/-3.3% on the muon neutral current (POWHEG-RES) and "
             "+0.9/-0.7% on the neutrino charged current (POWHEG-V2), a "
             "factor of four to five",
     "check": lambda: (abs(100 * _mhou("results")[0] - 4.29) < 0.2
                       and abs(100 * _mhou("results")[1] + 3.25) < 0.2
                       and abs(100 * _mhou("results_nu")[0] - 0.89) < 0.15
                       and abs(100 * _mhou("results_nu")[1] + 0.73) < 0.15
                       and 3.5 < (_mhou("results")[0] - _mhou("results")[1])
                       / (_mhou("results_nu")[0]
                          - _mhou("results_nu")[1]) < 5.5),
     "detail": lambda: "mu +%.2f/%.2f%%, nu +%.2f/%.2f%%, ratio %.1f" % (
         100 * _mhou("results")[0], 100 * _mhou("results")[1],
         100 * _mhou("results_nu")[0], 100 * _mhou("results_nu")[1],
         (_mhou("results")[0] - _mhou("results")[1])
         / (_mhou("results_nu")[0] - _mhou("results_nu")[1]))},
    # >>> THE SCALE BAND IS ON THE FIGURE, AND ITS ABSENCE IS AN ERROR. <<<
    # Paper plot 6 drew for days without its band because the file it read
    # had never been produced and nothing checked; this claim is that
    # oversight written down.  It requires the band file for BOTH columns and
    # that the band be non-trivial on the observable the figure is about.
    {"what": "the 7-point scale band of the matched calculation is drawn on "
             "both columns",
     "check": lambda: all(
         (mhou_band.rel(r, "nch05", SELECTION, esuffix=TSUF) is not None
          and float(np.max(np.abs(mhou_band.rel(r, "nch05", SELECTION, esuffix=TSUF)[2])))
          > 0.0) for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s %s" % (r, "absent" if mhou_band.rel(r, "nch05", SELECTION, esuffix=TSUF) is None
                    else "+%.2f/-%.2f%% on N_ch"
                    % (100 * float(np.max(mhou_band.rel(r, "nch05", SELECTION, esuffix=TSUF)[2])),
                       100 * abs(float(np.min(mhou_band.rel(r, "nch05", SELECTION, esuffix=TSUF)[1]))))
                    ) for r in ("results", "results_nu"))},
    # THE SELECTION IS FASER'S, and it is checked against the numbers in
    # arXiv:2403.12520 rather than asserted in prose.
    {"what": "the FASER&nu; selection is that of arXiv:2403.12520: a lepton "
             "above 200 GeV at theta > 5 mrad, at least five charged hadrons "
             "inside tan(theta) < 0.5 of which four inside 0.1, and "
             "Delta phi > pi/2",
     "check": lambda: all(
         (_sel(r).get("e_lep_min") == 200.0
          and _sel(r).get("w2_min") == 9.0 and _sel(r).get("y_min") == 0.0
          and _sel(r).get("y_max") == 1.0
          and _sel(r).get("theta_min") == 0.005
          and _sel(r).get("n05_min") == 5
          and _sel(r).get("n01_min") == 4
          and abs(_sel(r).get("dphi_x_min", 0) - math.pi / 2) < 1e-9)
         for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s E'>%g, theta>%g, n05>=%d, n01>=%d, dphi>%.4f"
         % (r, _sel(r).get("e_lep_min"), _sel(r).get("theta_min"),
            _sel(r).get("n05_min"), _sel(r).get("n01_min"),
            _sel(r).get("dphi_x_min")) for r in ("results", "results_nu"))},
    # THE FIGURE'S OWN CONSISTENCY: the Delta phi row is the variable the
    # selection cuts on, so nothing may appear below the cut.  If this fails,
    # the figure is drawing the OTHER dphi.
    {"what": "the Delta phi row is the variable the selection cuts on, so no "
             "rate at all appears below pi/2 on either current",
     "check": lambda: all((_dphix_below_cut(r) or 0.0) == 0.0
                          for r in ("results", "results_nu")),
     "detail": lambda: ", ".join(
         "%s %.3g" % (r, _dphix_below_cut(r))
         for r in ("results", "results_nu"))},
    {"what": "every event on the figure has a lepton above 200 GeV, the "
             "emulsion's momentum threshold",
     "check": lambda: all(
         (_integral(r, DENOM[r], "Emu", hi=200.0) or 0.0) == 0.0
         for r in ("results", "results_nu")),
     "detail": lambda: ", ".join(
         "%s %.3g below 200 GeV" % (r, _integral(r, DENOM[r], "Emu",
                                                 hi=200.0))
         for r in ("results", "results_nu"))},
    {"what": "the mean charged-hadron multiplicity inside the emulsion "
             "acceptance spreads by no more than 26% across the four "
             "generators, on both currents",
     "check": lambda: all(1.0 <= (_mult_spread(r) or 0) <= 1.26
                          for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s: " % r + ", ".join(
             "%s %.2f" % (lab.split(" (")[0], _mean(r, k, "nch05"))
             for lab, k, _c, _m, _ls in (MU_ROWS if r == "results"
                                         else NU_ROWS)
             if _mean(r, k, "nch05"))
         for r in ("results", "results_nu"))},
    {"what": "GENIE's mean multiplicity is BELOW both matched generators on "
             "each current -- the hadronic activity is where the leading-order "
             "tune differs most",
     "check": lambda: all(
         _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[1][1], "nch05")
         < min(_mean(r, (MU_ROWS if r == "results" else NU_ROWS)[i][1],
                     "nch05") for i in (0, 2))
         for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s GENIE %.2f vs POWHEG %.2f, Sherpa %.2f"
         % (r, _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[1][1],
                     "nch05"),
            _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[0][1], "nch05"),
            _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[2][1], "nch05"))
         for r in ("results", "results_nu"))},
    {"what": "Herwig has the highest multiplicity on each current",
     "check": lambda: all(
         _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[3][1], "nch05")
         > max(_mean(r, (MU_ROWS if r == "results" else NU_ROWS)[i][1],
                     "nch05") for i in (0, 1, 2))
         for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s Herwig %.2f" % (r, _mean(r, (MU_ROWS if r == "results"
                                          else NU_ROWS)[3][1], "nch05"))
         for r in ("results", "results_nu"))},
    {"what": "below 10 GeV Herwig's leading charged hadron is 2.1 times the "
             "matched rate on both currents",
     "check": lambda: all(
         2.0 < _integral(r, (MU_ROWS if r == "results" else NU_ROWS)[3][1],
                         "Elead", hi=10.0, raw=True)
         / _integral(r, DENOM[r], "Elead", hi=10.0, raw=True) < 2.2
         for r in ("results", "results_nu")),
     "detail": lambda: ", ".join("%s %.2f" % (r, _integral(
         r, (MU_ROWS if r == "results" else NU_ROWS)[3][1], "Elead",
         hi=10.0, raw=True) / _integral(r, DENOM[r], "Elead", hi=10.0,
                                        raw=True))
         for r in ("results", "results_nu"))},
    {"what": "Sherpa is ABOVE POWHEG and GENIE on each current, so the "
             "generators are ordered the same way on the muon and the "
             "neutrino side",
     "check": lambda: all(
         _mean(r, (MU_ROWS if r == "results" else NU_ROWS)[2][1], "nch05")
         > max(_mean(r, (MU_ROWS if r == "results" else NU_ROWS)[i][1],
                     "nch05") for i in (0, 1))
         for r in ("results", "results_nu")),
     "detail": lambda: "; ".join(
         "%s Sherpa %.2f" % (r, _mean(r, (MU_ROWS if r == "results"
                                          else NU_ROWS)[2][1], "nch05"))
         for r in ("results", "results_nu"))},
    # THE SOFT LEADING HADRON IS WHERE GENIE PARTS COMPANY, and the hard end
    # is where it does not -- a claim on both halves, since only the pair of
    # them says the effect is in the fragmentation rather than in the rate.
    # Until the 1 GeV track threshold entered Tier E (2026-09-19) the muon
    # side read GENIE 1.0 and Sherpa 1.3: the threshold removes soft-track
    # events, most in POWHEG-RES, so the ratios converge on the muon side.
    {"what": "below 10 GeV GENIE's leading charged hadron is 1.2 times "
             "POWHEG's rate on the muon side and 1.4 times it on the neutrino "
             "side, Sherpa at 1.2 and 1.8; above 100 GeV every generator "
             "agrees with POWHEG to better than 10%",
     "check": lambda: (
         1.15 < _integral("results", "genie", "Elead", hi=10.0, raw=True)
         / _integral("results", DENOM["results"], "Elead", hi=10.0, raw=True) < 1.3
         and 1.15 < _integral("results", MU_ROWS[2][1], "Elead", hi=10.0, raw=True)
         / _integral("results", DENOM["results"], "Elead", hi=10.0, raw=True) < 1.3
         and 1.7 < _integral("results_nu", NU_ROWS[2][1], "Elead", hi=10.0, raw=True)
         / _integral("results_nu", DENOM["results_nu"], "Elead", hi=10.0, raw=True) < 1.9
         and 1.35 < _integral("results_nu", "genie_lo", "Elead", hi=10.0, raw=True)
         / _integral("results_nu", DENOM["results_nu"], "Elead", hi=10.0, raw=True) < 1.45
         and all(abs(_integral(r, k, "Elead", lo=100.0)
                     / _integral(r, DENOM[r], "Elead", lo=100.0) - 1.0) < 0.10
                 for r in ("results", "results_nu")
                 for _l, k, _c, _m, _ls in (MU_ROWS if r == "results"
                                            else NU_ROWS))),
     "detail": lambda: "; ".join(
         "%s soft GENIE %.2f, Sherpa %.2f | hard GENIE %.3f, Sherpa %.3f"
         % (r,
            _integral(r, (MU_ROWS if r == "results" else NU_ROWS)[1][1],
                      "Elead", hi=10.0, raw=True) / _integral(r, DENOM[r], "Elead",
                                                    hi=10.0, raw=True),
            _integral(r, (MU_ROWS if r == "results" else NU_ROWS)[2][1],
                      "Elead", hi=10.0, raw=True) / _integral(r, DENOM[r], "Elead",
                                                    hi=10.0, raw=True),
            _integral(r, (MU_ROWS if r == "results" else NU_ROWS)[1][1],
                      "Elead", lo=100.0) / _integral(r, DENOM[r], "Elead",
                                                     lo=100.0),
            _integral(r, (MU_ROWS if r == "results" else NU_ROWS)[2][1],
                      "Elead", lo=100.0) / _integral(r, DENOM[r], "Elead",
                                                     lo=100.0))
         for r in ("results", "results_nu"))},
    # THE LEPTON ROW IS FLAT: the ratio in the first half of the momentum
    # range equals the ratio in the second to a few per cent, so what the
    # panel shows is a normalisation and not a shape.
    {"what": "the outgoing-lepton ratio is flat -- its low and high halves "
             "agree to better than 10% -- so the generators differ in how "
             "many events pass the selection, not in the lepton spectrum",
     "check": lambda: all(
         abs((_integral(r, k, "Emu", lo=500.0)
              / _integral(r, DENOM[r], "Emu", lo=500.0))
             / (_integral(r, k, "Emu", hi=500.0)
                / _integral(r, DENOM[r], "Emu", hi=500.0)) - 1.0) < 0.10
         for r in ("results", "results_nu")
         for _l, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS)),
     "detail": lambda: "; ".join(
         "%s %s low %.3f high %.3f"
         % (r, lab.split(" (")[0],
            _integral(r, k, "Emu", hi=500.0)
            / _integral(r, DENOM[r], "Emu", hi=500.0),
            _integral(r, k, "Emu", lo=500.0)
            / _integral(r, DENOM[r], "Emu", lo=500.0))
         for r in ("results", "results_nu")
         for lab, k, _c, _m, _ls in (MU_ROWS if r == "results" else NU_ROWS)
         if k != DENOM[r])},
]


MESSAGE = (
    "<b>This is the half of the comparison an analytic calculation cannot "
    "reach.</b> Paper plot 2 draws the observables structure functions "
    "predict, where the matched generators agree with YADISM (ZM-VFNS) to a "
    "per cent or two. Everything on this page needs a shower and a "
    "hadronisation model instead, so there is no reference curve; each "
    "column is normalised to its own matched calculation, POWHEG-RES on the "
    "muon side and POWHEG-V2 on the neutrino side, per nucleon on tungsten."
    "\n\n"
    "<b>The selection is FASER&rsquo;s own</b>, nested in the benchmark region "
    "(Q&sup2; &gt; 4 GeV&sup2;, W &gt; 3 GeV): at least five charged hadrons "
    "inside tan&thinsp;&theta; &lt; 0.5 of which four inside 0.1, a lepton "
    "above 200 GeV at &theta; &gt; 5 mrad, and &Delta;&phi; &gt; &pi;/2 between "
    "the lepton and the summed charged-hadron system, after arXiv:2403.12520. "
    "The third row is populated above &pi;/2 only."
    "\n\n"
    "<b>The four generators order themselves the same way on both currents, "
    "and the ordering is the multiplicity.</b> The mean charged-hadron "
    "multiplicity inside the acceptance is 7.74 (POWHEG-RES), 7.05 (GENIE), "
    "8.23 (Sherpa) and 8.85 (Herwig) on the muon side and 7.80, 7.26, 8.08, "
    "8.24 on the neutrino side: GENIE below the matched generators, Sherpa "
    "above POWHEG-RES and POWHEG-V2 and Herwig highest, a spread of up to "
    "26%. The "
    "vertex selection is a cut on exactly this distribution, so its "
    "efficiency inherits the spread."
    "\n\n"
    "<b>The soft leading hadron separates the generators on the neutrino "
    "side.</b> Below 10 GeV &mdash; measured, though below the range the panel "
    "draws &mdash; GENIE puts 1.4 times POWHEG-V2&rsquo;s rate and Sherpa 1.8 "
    "times on the neutrino side, while on the muon side both are 1.2 times "
    "POWHEG-RES; Herwig puts 2.1 times the matched rate there on both "
    "currents. Above 100 GeV all agree to within "
    "7%."
    "\n\n"
    "<b>The outgoing-lepton row is nearly flat.</b> Its ratio panels sit at "
    "0.88&ndash;0.95 (GENIE), 1.07&ndash;1.23 (Sherpa) and 0.98&ndash;1.13 "
    "(Herwig): the generators "
    "agree on the shape of the lepton spectrum and disagree on how many "
    "events pass the selection, which enters through the hadronic "
    "requirements."
    "\n\n"
    "<b>The matched calculation carries its own scale band.</b> The seven "
    "scale points are joined to the showered events through the Les Houches "
    "index, the only route that survives a selection cutting on the "
    "multiplicity, and each band is measured on the curve it is drawn on. "
    "Integrated over the selection it is +4.3/&minus;3.3% on the muon neutral "
    "current and +0.9/&minus;0.7% on the neutrino charged current, a factor "
    "of four to five; differentially every generator difference above is "
    "larger than it."
)


if __name__ == "__main__":
    main()
