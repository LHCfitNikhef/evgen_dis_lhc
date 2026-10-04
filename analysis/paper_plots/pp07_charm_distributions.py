#!/usr/bin/env python3
"""Paper plot 7: NLO differential charm distributions at 1 TeV, in the
CHARM region of the final paper.

The redo of the earlier figure for the "Paper plots"
tab: charm distributions in x, E_h and Q2 for the NLO matchings
(plus POWHEG-V2mc on the charged current) against YADISM charm in ZM-VFNS,
with FONLL beside it and NNLO with its band.  The conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184.
  * THE CHARM REGION, Q2 > 4 GeV2 and W > 5 GeV, no y cut, BOTH CURRENTS
    (user, 2026-09-14; selection q4w5): the earlier muon floor Q2 > 11 does not
    carry over.
  * THE REFERENCES AND THEIR PER-BIN BANDS carry target-mass corrections
    (analysis/mhou_diff.py --flavour charm with BENCH_SELECTION=q4w5
    BENCH_TARGET=W BENCH_TMC=3).
  * With no y window the hadronic energy E_h = y E runs over the whole beam
    energy.

Usage: analysis/paper_plots/pp07_charm_distributions.py
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
import selection                                            # noqa: E402

SLUG = "pp07_charm_distributions"
V1_NUMBER = 7
IN_PAPER = True
# TITLE and CAPTION are rendered as HTML on the report page, so they carry no
# LaTeX; the figure's own suptitle is FIG_TITLE, which goes through tex().
TITLE = "NLO differential charm production distributions at 1 TeV"
FIG_TITLE = (r"NLO differential charm production distributions at "
             r"$E_{\ell} = 1$ TeV")
OUTPUT = "pp_charm_distributions.png"
RESULTS = "results_nu"
CAPTION = ("Differential charm distributions in Bjorken <i>x</i>, hadronic "
           "energy and momentum transfer, per nucleon on tungsten at 1 TeV "
           "for Q&sup2; &gt; 4 GeV&sup2; and W &gt; 5 GeV, for the three NLO "
           "matchings; muon neutral current (left) and neutrino charged "
           "current (right). Ratio to YADISM charm (ZM) beneath each, with "
           "YADISM charm (FONLL) drawn as well &mdash; the gap between them "
           "is the charm-mass effect, which every generator here is missing "
           "by construction. The NNLO result is shown where it exists. "
           "POWHEG-V2mc, the one massive-charm entry, appears on the "
           "charged current only. The references and their scale bands "
           "include target-mass corrections. Herwig&rsquo;s muon charm falls "
           "10&ndash;40% short above x &asymp; 0.25 (about 2% of the charm "
           "rate): events whose beam remnant it can build only by displacing "
           "the scattered lepton, which is not allowed here.")

SURFACE = "#ffffff"
ENERGY = beams.ANCHOR_ENERGY

# THE Q2 FLOOR IS PER CURRENT (user, 2026-09-01), exactly as on paper plot 2,
# and the reason bites hardest here: the neutral-current charm series is the
# observable most sensitive to the bottom of the Q2 range, because the charm
# threshold sits inside it.  At Q2 > 4 the muon charm K-factor is 0.53 at
# 400 GeV and the NLO band 95% wide; at Q2 > 11 they are 0.76 and 49%.  The
# charged current is well behaved at Q2 > 4 and keeps it.
#
# Going back is one edit -- set "results" to 4.0 -- and every result computed
# above the default carries the floor in its filename, so both regions sit on
# disk at once.  Re-derive with:  tools/make_q2min_muon.sh 11
Q2MIN = {"results": 4.0, "results_nu": 4.0}
REGION = "q4w5"
TSUF = "_W"
TMC = "_tmc"


def _suffix(resdir, key=""):
    """_q4w5_W[_tmc][_tag]: the references carry the TMC tag, generators none."""
    tag = beams.Beams("mu", ENERGY).tag
    esuf = "" if ENERGY == beams.ANCHOR_ENERGY else f"_{tag}"
    return f"_{REGION}{TSUF}{TMC if key.startswith('yadism') else ''}{esuf}"

# THE REFERENCE IS THE MASSLESS CHARM CALCULATION, for the reason in the
# docstring: every generator on this figure computes massless charm, so a
# FONLL denominator would charge them for a mass they never put in.
REF = "yadism_charm_nlo"
REF_LABEL = "YADISM (ZM)"
FONLL_REF = "yadism_charm_nlo_fonll_damp"
FONLL_LABEL = "YADISM (FONLL)"
REF_COLOUR = "#111111"
FONLL_COLOUR = "#c0392b"
# POWHEG-V2mc's own scheme, neutrino column only -- as on plot 6, in plot 6's
# colour and label (user, 2026-10-04: "inconsistent with Fig. 5.3").
FFNS3_REF = "yadism_charm_nlo_ffns3"
FFNS3_LABEL = r"YADISM (FFNS $n_f\!=\!3$)"
FFNS3_COLOUR = "#cc79a7"

# NNLO, in the scheme each current HAS at that order -- FONLL on the muon
# side, ZM on the neutrino side, where the massive charged-current
# coefficient functions do not exist.  Same arrangement as plots 1 and 3.
NNLO = {"results": ("yadism_charm_nnlo_fonll_damp",
                    r"YADISM NNLO (FONLL, $\mu$)"),
        "results_nu": ("yadism_charm_nnlo", r"YADISM NNLO (ZM, $\nu$)")}
NNLO_COLOUR = "#6a3d9a"
NNLO_LS = {"results": "-.", "results_nu": (0, (5.0, 1.3, 1.0, 1.3, 1.0, 1.3))}

# THE SEVEN-POINT BANDS, at both orders, exactly as on plot 3 (user,
# 2026-08-31).  They needed a calculation that did not exist: mhou_diff.py
# integrated the TOTAL structure functions into bins, and these are the CHARM
# ones -- `--flavour charm`, which is the same code with the charm observable
# names, since both charm calc modules integrate with the same make_sigma_red
# and integrate_* as the inclusive ones.
#
# EACH BAND IS IN THE SCHEME OF THE CURVE IT SURROUNDS, the rule this whole
# set of figures follows: the NLO band is ZM, matching the massless reference
# the ratio is taken to, and the NNLO band is FONLL on the muon side and ZM on
# the neutrino side, matching the NNLO curve drawn there.
MHOU_FILE = "mhou_diff_zm_charm.json"
NNLO_MHOU_FILE = {"results": "mhou_diff_fonll_charm_pto2.json",
                  "results_nu": "mhou_diff_zm_charm_pto2.json"}
# The xiF = 1/2 points are KEPT, as on plot 3: a per-bin band is read as a
# shape rather than quoted as a number, and dropping them would narrow it
# precisely in the bins where the instability lives.
MHOU_KEEP_HALF = True

# (label, result key, colour, marker, linestyle) -- the keys are pp02's, so
# the two charm figures name the same samples, and the colours are pp01's and
# pp03's, so a reader carries one legend across the whole set.
# ONE LEGEND, so the two POWHEGs carry different line styles -- see the note
# in pp03; the same user instruction (2026-09-01) and the same fix.
POWHEG_V2_LS = (0, (6.5, 1.6))
NU_ROWS = [
    (r"POWHEG-V2 ($\nu$)", "powheg_nu_charmfinal",     "#0072b2", "o",
     POWHEG_V2_LS),
    ("Herwig",      "herwig_nlo_full_charmfinal",        "#d55e00", "s", "--"),
    ("Sherpa",      "sherpa_nlo_charmfinal",             "#009e73", "^", "-."),
    # THE MASSIVE-CHARM ENTRY, charged current only.  Drawn in the same purple
    # family as the NNLO curve would be misleading, so it gets its own colour
    # and a dotted style; it is a GENERATOR, not a reference.
    ("POWHEG-V2mc", "powheg_nu_mc_charmfinal",           "#8c6d1f", "D", ":"),
]
MU_ROWS = [
    (r"POWHEG-RES ($\mu$)", "powheg_charmfinal",       "#0072b2", "o", "-"),
    ("Herwig",      "herwig_nlo_powheg_full_charmfinal", "#d55e00", "s", "--"),
    ("Sherpa",      "sherpa_charmfinal",                 "#009e73", "^", "-."),
]

# X RANGES, the user's values (2026-09-01), written as the PHYSICAL numbers
# and snapped to a bin edge where they are used (_snap below) rather than
# hand-computed into the table -- the hand-computed edges went stale the
# moment the muon Q2 floor moved, which it just did.
#
#   x_Bj    0.003 < x < 0.6 (muon), 0.003 < x < 0.3 (neutrino)
#   Q2      < 100 GeV2 (muon), < 400 GeV2 (neutrino)
#
# THE x > 0.003 FLOOR IS NEW AND IS THE USER'S ("for pp04 and the x_bj
# distributions, please impose x_bj > 0.003").  It removes the lowest bins,
# where the charm samples -- a tenth of the inclusive rate -- carry error bars
# several times the effect the panel is about.  On the muon column it changes
# nothing that the Q2 > 11 floor has not already removed: the smallest x
# reachable there is 11/(0.9 * 2k.P) = 0.0065.
# {observable: {resdir: (lo, hi)}}; anything unlisted keeps the full range.
# THE NEUTRINO x WINDOW ENDS AT 0.3 (user, 2026-09-01), not 0.6.  On the
# charged current the charm rate above 0.3 is a tail of wide error bars --
# the samples are a tenth of the inclusive rate and the last bins hold a
# handful of weighted events -- and it was that tail, not the bulk, that set
# the panel's ratio range.  The muon column keeps 0.6: the Q2 > 11 floor has
# already ended it at 0.24 through the bulk cut, so the window does not bite
# there and stating it differently would only invite the reader to think the
# two columns were treated by different rules.
XLIM = {
    "xbj": {"results": (3e-3, 0.6), "results_nu": (3e-3, 0.3)},
    "Q2": {"results": (None, 100.0), "results_nu": (None, 400.0)},
}


def _snap(edges, v, upper):
    """v moved to the nearest edge of the (possibly merged) grid."""
    e = np.asarray(edges, dtype=float)
    if v is None:
        return float(e[-1] if upper else e[0])
    return float(e[int(np.argmin(np.abs(e - v)))])


# REBINNING, on the same principle as plot 3: the merge is applied to EVERY
# curve on a panel, the reference included, so no ratio is ever taken between
# differently binned quantities.  E_h is written in 35 bins of 20 GeV and Q2
# in 26 geometric ones, and on the charm samples -- which are a tenth of the
# inclusive rate -- the bin-to-bin scatter is correspondingly worse.  Both
# rows are merged on both currents, so the two columns of a row always share
# one binning.
# THE MUON Q2 ROW IS NO LONGER MERGED (2026-09-01), for the reason given in
# pp03: with the floor at 11 and the window at 100 the merge left six bins
# across a distribution that falls two decades.  Twelve unmerged bins carry
# the shape, and on charm the scatter is still readable because the window
# now excludes the low-Q2 corner where the charm sample is thinnest.
REBIN = {("nu", "results"): 5, ("nu", "results_nu"): 5,
         ("Q2", "results_nu"): 2}

# (histogram key, x label, x scale, y scale, y quantity).  THREE rows: no
# theta_l' (a leptonic variable, already on plot 3) and the hadronic energy
# in place of the lepton energy -- see the docstring.
OBS = [
    ("xbj", r"$x_{\rm Bj}$", "log", "log", r"d\sigma_c/dx_{\rm Bj}"),
    ("nu",  r"$E_h$  [GeV]", "linear", "linear", r"d\sigma_c/dE_h"),
    ("Q2",  r"$Q^2$  [GeV$^2$]", "log", "log", r"d\sigma_c/dQ^2"),
]

# Per-bin units, as on plot 3: the neutrino charm cross-section is sub-pb and
# the muon one is nb, so a shared unit would put one column in scientific
# notation throughout.
# THE PANEL TITLE NAMES THE CUT AND SAYS "charm production" (user,
# 2026-09-01), built from Q2MIN so it cannot go stale against the data.
COLUMNS = [("results", MU_ROWS, "Muon DIS", 1000.0, "nb"),
           ("results_nu", NU_ROWS, "Neutrino DIS", 1.0, "pb")]


def column_title(name, resdir):
    # NO REGION IN THE TITLE (user, 2026-09-21): the cuts are in the caption.
    return rf"{name}, charm production"

def _merge(edges, dsig, err, n):
    """Merge n adjacent bins, keeping dsig a density and err its error.

    Errors add in quadrature WEIGHTED BY BIN WIDTH, because what adds is the
    integral over each bin, not the density.  With uniform bins the two are
    the same and the distinction would never show; the neutrino
    lepton-energy bins are uniform, so this is written correctly rather than
    correctly by accident.
    """
    m = (len(dsig) // n) * n
    w = np.diff(edges)[:m].reshape(-1, n)
    y = (dsig[:m].reshape(-1, n) * w).sum(axis=1)
    e = np.sqrt(((err[:m].reshape(-1, n) * w) ** 2).sum(axis=1))
    wide = w.sum(axis=1)
    return edges[:m + 1:n], y / wide, e / wide


def hist(resdir, key, obs):
    """(edges, dsig, err) for one observable, or None if absent.

    THE FILE'S OWN STAMP IS CHECKED against Q2MIN, as on paper plot 2: a
    Q2 > 4 file read into a Q2 > 11 panel would draw a plausible curve and
    say nothing.
    """
    p = f"{BASE}/{resdir}/histos_{key}{_suffix(resdir, key)}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    q2 = (d.get("selection") or {}).get("q2_min")
    if q2 is not None and abs(q2 - Q2MIN[resdir]) > 1e-9:
        raise SystemExit(f"{os.path.relpath(p, BASE)} was computed at "
                         f"Q2 > {q2:g}, but this figure's {resdir} column is "
                         f"Q2 > {Q2MIN[resdir]:g}")
    if d.get("stat_insufficient") or obs not in d.get("hists", {}):
        return None
    h = d["hists"][obs]
    out = (np.asarray(h["edges"]), np.asarray(h["dsig"]),
           np.asarray(h["err"]))
    n = REBIN.get((obs, resdir))
    return _merge(*out, n) if n else out


def raw_hist(resdir, key, obs):
    """The histogram BEFORE any rebinning.

    The mirror check compares the lepton-energy and energy-transfer rows,
    and one of them is merged for the figure while the other is not.  Merging
    is a presentation choice; the identity being checked is not, so it is
    checked on the numbers the analysis wrote.
    """
    saved = dict(REBIN)
    REBIN.clear()
    try:
        return hist(resdir, key, obs)
    finally:
        REBIN.update(saved)


def mhou(resdir, obs, fn=None):
    """(edges, rel_hi, rel_lo) of the 7-point band, or None.

    The band is a RATIO per bin, so merging it means averaging the edges
    weighted by the bin's own cross-section -- which is what re-deriving it
    from the varied and central histograms does.  Both are stored, so the
    merge is done on the cross-sections and the ratio taken afterwards,
    never on the ratios themselves.
    """
    p = (f"{BASE}/{resdir}/"
         + (fn or MHOU_FILE).replace(".json", f"_{REGION}{TSUF}{TMC}.json"))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        h = json.load(f)["hists"].get(obs)
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
        e, cen, _ = _merge(e, cen, z, n)
        _, hi, _ = _merge(np.asarray(h["edges"]), hi, z, n)
        _, lo, _ = _merge(np.asarray(h["edges"]), lo, z, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (e, np.where(cen > 0, hi / cen - 1.0, 0.0),
                np.where(cen > 0, lo / cen - 1.0, 0.0))


def _bulk(edges, y):
    """Mask of bins carrying at least 0.1% of the distribution.

    The ratio in a bin holding a thousandth of the cross-section is a
    statistical statement about a handful of events, not a physics one; the
    tails of x_Bj and of theta are made of such bins.  They are still DRAWN --
    hiding them would be a different figure -- but they do not set the ratio
    range, which otherwise collapses every per-cent effect the panel exists to
    show into a flat line.
    """
    w = y * np.diff(edges)
    return w > 0.001 * w.sum()


def ratio_range(rows):
    """Symmetric ratio range from the bulk bins of one observable row."""
    dev = 0.0
    for r, m in rows:
        v = r[m & np.isfinite(r)] if m is not None else r[np.isfinite(r)]
        if v.size:
            dev = max(dev, float(np.max(np.abs(v - 1.0))))
    dev = max(dev, 0.01)
    return 1.0 - 1.25 * dev, 1.0 + 1.25 * dev


def _log_ticks(a, lo, hi):
    """1-2-5 ticks on a log axis that spans less than two decades.

    Matplotlib's default puts a label on every minor tick there, and on the
    muon theta panel -- 0.002 to 0.04, barely more than one decade -- they
    overlap into an unreadable smear.  Placing 1, 2 and 5 per decade and
    labelling them plainly is what a publication figure does.
    """
    if hi / lo >= 100.0:
        return
    ticks = [m * 10.0 ** k for k in range(-6, 7) for m in (1, 2, 5)
             if lo <= m * 10.0 ** k <= hi]
    a.set_xticks(ticks)
    a.set_xticklabels([f"{t:g}" for t in ticks])
    a.xaxis.set_minor_formatter(ticker.NullFormatter())


def _x_kin_min(resdir):
    """The smallest x_Bj the fiducial region can reach, x = Q2min/(y_max 2k.P).

    A histogram bin that STRADDLES this boundary is only partly inside the
    region, so its density is right but it holds a fraction of the events the
    bin next to it does -- which is why it arrives with an error bar several
    times everyone else's and, being included, sets the ratio range for the
    whole panel.  It is not a measurement worth that.  Excluding it is a
    kinematic statement, not a cosmetic one, and it moves with the Q2 floor
    rather than being a literal that goes stale.
    """
    return Q2MIN[resdir] / (0.9 * beams.Beams("mu", ENERGY).two_kP)


def x_lo(obs, resdir, edges, y=None):
    """The window's lower edge, defaulting to the first POPULATED bin."""
    v = XLIM.get(obs, {}).get(resdir, (None, None))[0]
    if obs == "xbj":
        # the requested floor OR the kinematic one, whichever is higher, and
        # snapped UP so no straddling bin survives
        e = np.asarray(edges, dtype=float)
        want = max(v if v is not None else 0.0, _x_kin_min(resdir))
        above = e[e >= want - 1e-12]
        return float(above[0]) if above.size else float(e[0])
    if v is None and y is not None:
        nz = np.nonzero(np.asarray(y) > 0)[0]
        if nz.size:
            return float(np.asarray(edges)[nz[0]])
    return _snap(edges, v, False)


def x_hi(obs, resdir, edges):
    return _snap(edges, XLIM.get(obs, {}).get(resdir, (None, None))[1], True)


def block(ax, axr, resdir, rows, obs, div, unit, xs, ys, ylab_q):
    """One observable in one current: absolute panel plus its ratio panel."""
    ref = hist(resdir, REF, obs)
    if ref is None:
        return
    edges, ry, _ = ref
    mids = 0.5 * (edges[:-1] + edges[1:])
    mask = _vis(resdir, obs)
    seen = []

    # the analytic reference: dotted black, a style none of the generators
    # uses, so it reads as a different KIND of object rather than a fourth
    # generator
    ax.stairs(ry / div, edges, color=REF_COLOUR, lw=1.7, ls=":",
              baseline=None, zorder=7, label=tex(REF_LABEL))

    # the 7-point band on the reference, drawn FIRST so everything sits on
    # top of it, and only in the ratio panel -- in the absolute panel it
    # would be a shaded strip under curves that already coincide
    mb = mhou(resdir, obs)
    if mb is not None and np.array_equal(mb[0], edges):
        axr.stairs(1.0 + mb[1], edges, baseline=1.0 + mb[2], fill=True,
                   color=REF_COLOUR, alpha=0.12, lw=0, zorder=0,
                   label=tex("NLO MHOU"))
        # AND IN THE ABSOLUTE PANEL (user, 2026-10-04), around the reference
        # it belongs to; unlabelled, the ratio panel's entry names it.  The
        # lower edge is floored above zero for the log axes.
        ax.stairs(ry * (1.0 + mb[1]) / div, edges,
                  baseline=np.maximum(ry * (1.0 + mb[2]), 1e-12 * ry.max()) / div,
                  fill=True, color=REF_COLOUR, alpha=0.12, lw=0, zorder=0)

    nnlo_key, nnlo_label = NNLO[resdir]
    nh = hist(resdir, nnlo_key, obs)
    if nh is not None:
        ax.stairs(nh[1] / div, edges, color=NNLO_COLOUR, lw=1.8,
                  ls=NNLO_LS[resdir], baseline=None, zorder=4,
                  label=tex(nnlo_label))
        with np.errstate(divide="ignore", invalid="ignore"):
            nr = np.where(ry > 0, nh[1] / ry, np.nan)
        axr.stairs(nr, edges, color=NNLO_COLOUR, lw=1.8,
                   ls=NNLO_LS[resdir], baseline=None, zorder=4)
        seen.append((nr, mask))
        # The NNLO band goes around the NNLO CURVE, not around unity: this
        # panel is normalised to the NLO reference, so the band's centre sits
        # at nr.  Multiplying the relative envelope by nr is what puts it
        # there; drawn about 1.0 it would be an uncertainty on the wrong
        # calculation and would look perfectly reasonable.
        # A FLAT FILL, not hatched (user, 2026-09-18: the bands are drawn
        # the same way in every figure of the paper, as in pp01 and pp06).
        nb = mhou(resdir, obs, NNLO_MHOU_FILE.get(resdir))
        if nb is not None and np.array_equal(nb[0], edges):
            axr.stairs(nr * (1.0 + nb[1]), edges,
                       baseline=nr * (1.0 + nb[2]), fill=True,
                       color=NNLO_COLOUR, alpha=0.20, lw=0, zorder=2,
                       label=tex("NNLO MHOU"))
            seen.append((nr * (1.0 + nb[1]), mask))
            seen.append((nr * (1.0 + nb[2]), mask))
            # and around the NNLO curve in the absolute panel (user,
            # 2026-10-04)
            ax.stairs(nh[1] * (1.0 + nb[1]) / div, edges,
                      baseline=np.maximum(nh[1] * (1.0 + nb[2]),
                                          1e-12 * nh[1].max()) / div,
                      fill=True, color=NNLO_COLOUR, alpha=0.20, lw=0,
                      zorder=1)

    fon = hist(resdir, FONLL_REF, obs)
    if fon is not None:
        ax.stairs(fon[1] / div, edges, color=FONLL_COLOUR, lw=1.7, ls="-",
                  baseline=None, zorder=6, label=tex(FONLL_LABEL))
        with np.errstate(divide="ignore", invalid="ignore"):
            fr = np.where(ry > 0, fon[1] / ry, np.nan)
        axr.stairs(fr, edges, color=FONLL_COLOUR, lw=1.7, ls="-",
                   baseline=None, zorder=6)
        seen.append((fr, mask))

    ff = hist(resdir, FFNS3_REF, obs) if resdir == "results_nu" else None
    if ff is not None:
        ax.stairs(ff[1] / div, edges, color=FFNS3_COLOUR, lw=1.5, ls="-",
                  baseline=None, zorder=6, label=tex(FFNS3_LABEL))
        with np.errstate(divide="ignore", invalid="ignore"):
            f3 = np.where(ry > 0, ff[1] / ry, np.nan)
        axr.stairs(f3, edges, color=FFNS3_COLOUR, lw=1.5, ls="-",
                   baseline=None, zorder=6)
        seen.append((f3, mask))

    for lab, key, colour, _marker, ls in rows:
        h = hist(resdir, key, obs)
        if h is None:
            continue
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
    # THROUGH x_lo/x_hi, so the axis and _vis()'s mask cannot disagree.
    x0, x1 = x_lo(obs, resdir, edges, ry), x_hi(obs, resdir, edges)
    for a in (ax, axr):
        a.set_xscale(xs)
        a.set_xlim(x0, x1)
        if xs == "log":
            _log_ticks(a, x0, x1)
        a.grid(alpha=0.22, lw=0.6)
    ax.set_yscale(ys)
    if ys == "log":
        # ONE tail bin must not set the axis.  The muon theta distribution
        # has a FONLL bin nine decades below the peak, which on a shared log
        # axis compresses the whole shape into the top decade and makes the
        # figure unreadable.  The floor is taken from the bins that carry the
        # cross-section; anything below is drawn and runs off the bottom,
        # which is the honest way to show a bin that small.
        pos = ry[mask & (ry > 0)] / div
        if pos.size:
            ax.set_ylim(0.03 * float(pos.min()), 4.0 * float(pos.max()))
    # THE FUSED PAIR SHARES AN EDGE, so the absolute panel's lowest tick label
    # and the ratio panel's highest sit on top of one another -- visibly, on
    # the linear-axis rows.  Pruning one from each is the standard fix and
    # costs nothing: the pruned value is still on the axis, only unlabelled.
    if ys != "log":
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=5, prune="lower"))
    axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="upper"))
    ax.set_ylabel(tex(rf"${ylab_q}$  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL - 2)
    ax.tick_params(labelbottom=False)
    return seen


def main():
    fig = plt.figure(figsize=(11.4, 12.4))
    # Four observable blocks, each an absolute panel FUSED to its ratio panel
    # (hspace = 0 inside a block, plotstyle's house rule), with real space
    # between blocks so the rows do not read as one eight-panel grid.
    outer = fig.add_gridspec(len(OBS), 2, hspace=0.30, wspace=0.22,
                             top=0.845, bottom=0.055, left=0.085, right=0.985)
    for i, (obs, xlab, xs, ys, ylab_q) in enumerate(OBS):
        row_seen = []
        pair = []
        for j, (resdir, rows, title, div, unit) in enumerate(COLUMNS):
            # RATIO PANEL RAISED (user, 2026-08-31): was 2.3 : 1.  The ratio
            # panel is where every number this figure quotes is read, and at
            # 2.3 : 1 a 5% band and a 10% deviation were a few pixels apart.
            # The absolute panel only has to show the shape, which it still
            # does; the ratio panel has to be measurable by eye.
            inner = outer[i, j].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                            hspace=0.0)
            ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            seen = block(ax, axr, resdir, rows, obs, div, unit, xs, ys,
                         ylab_q)
            row_seen += seen or []
            pair.append((axr, seen or []))
            if i == 0:
                ax.set_title(tex(column_title(title, resdir)),
                             fontsize=plotstyle.FS_PANEL_TITLE - 1)
            axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 1)
            if j == 0:
                axr.set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                               fontsize=plotstyle.FS_YLABEL - 4)
        # EVERY RATIO PANEL SETS ITS OWN RANGE (user, 2026-09-01).  It matters
        # more here than on plot 3: the muon charm NNLO curve and its band run
        # to 0.4 of the reference while the neutrino column lives inside a few
        # per cent, so one shared range collapses the charged current to a
        # flat line -- the panel whose message is that its series DOES
        # converge.
        for a, seen_here in pair:
            lo, hi = ratio_range(seen_here)
            a.set_ylim(lo, hi)

    # ONE LEGEND, ABOVE THE FIGURE (user, 2026-09-01) -- see the note in
    # pp03.  The column-specific curves carry their own line style and their
    # own current in the label, so a single legend can name them all without
    # mislabelling either column.  Assembled from the top row's four axes
    # (the bands live in the ratio panel), de-duplicated by label, then put
    # into a fixed reading order: references, generators, bands.
    h, lab = [], []
    for a in fig.axes[:4]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    want = [tex(REF_LABEL), tex(FONLL_LABEL), tex(FFNS3_LABEL),
            tex(NNLO["results"][1]), tex(NNLO["results_nu"][1]),
            tex(MU_ROWS[0][0]), tex(NU_ROWS[0][0]),
            tex(MU_ROWS[1][0]), tex(MU_ROWS[2][0]), tex(NU_ROWS[3][0]),
            tex("NLO MHOU"), tex("NNLO MHOU")]
    order = [lab.index(w) for w in want if w in lab]
    order += [i for i in range(len(lab)) if i not in order]
    h, lab = [h[i] for i in order], [lab[i] for i in order]
    # each band shares its curve's entry (user, 2026-10-04)
    h, lab = plotstyle.merge_bands(h, lab, {
        tex("NLO MHOU"): [tex(REF_LABEL)],
        tex("NNLO MHOU"): [tex(NNLO["results"][1]), tex(NNLO["results_nu"][1])]})
    if h:
        # LEGEND FONT RAISED (user, 2026-10-04: "way too small", then "a
        # bit more"), from FS_LEGEND - 3; three columns so the entries fit.
        fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.962),
                   ncol=3, frameon=True, handlelength=2.1,
                   columnspacing=1.15, labelspacing=0.35,
                   fontsize=plotstyle.FS_LEGEND + 2)
    fig.suptitle(tex(FIG_TITLE), y=0.993, fontsize=plotstyle.FS_SUPTITLE)
    # NO PROSE INSIDE THE FIGURE (user rule, 2026-08-29): the explanation
    # lives in MESSAGE and in the caption.
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------- claims --

# ----------------------------------------------------------------- claims --
# Executable, and they re-derive themselves from the same JSONs the figure
# plots.  A claim that compared a hand-typed number to itself would check
# nothing (see tools/check_paper_plots.py).
_COLS = [("results", MU_ROWS), ("results_nu", NU_ROWS)]


def _vis(resdir, obs):
    """Mask of the bins this figure actually SHOWS for one panel.

    The bulk cut and the x range are applied in the same place for the
    drawing and for the claims -- on plot 3 they were not, briefly, and the
    prose then described bins outside the plotted window.
    """
    ref = hist(resdir, REF, obs)
    edges, ry, _ = ref
    return (_bulk(edges, ry) & (ry > 0)
            & (edges[:-1] >= x_lo(obs, resdir, edges, ry))
            & (edges[1:] <= x_hi(obs, resdir, edges)))


def _dev(resdir, key, obs):
    """Largest |ratio - 1| over the VISIBLE bins of one curve."""
    ref, h = hist(resdir, REF, obs), hist(resdir, key, obs)
    if ref is None or h is None:
        return None
    m = _vis(resdir, obs)
    return float(np.max(np.abs(h[1][m] / ref[1][m] - 1.0)))


def _panel_median_of(resdir, obs):
    """Largest median |ratio - 1| among the generators of one panel."""
    ref = hist(resdir, REF, obs)
    m = _vis(resdir, obs)
    return max(float(np.median(np.abs(hist(resdir, k, obs)[1][m]
                                      / ref[1][m] - 1.0)))
               for _, k, _c, _mk, _ls in dict(_COLS)[resdir]
               if hist(resdir, k, obs) is not None)


def _panel_median():
    """Largest median |ratio - 1| over the six panels, worst generator."""
    return max(_panel_median_of(d, o) for d, _r in _COLS for o, *_ in OBS)


def _herwig_highx(lo=0.238):
    """(Herwig/ZM in the visible muon x bins above `lo`, their share of the
    ZM charm rate)."""
    ref = hist("results", REF, "xbj")
    h = hist("results", "herwig_nlo_powheg_full_charmfinal", "xbj")
    e, vis = ref[0], _vis("results", "xbj")
    w = np.diff(e)
    sel = vis & (e[:-1] >= lo)
    return (list((h[1] / ref[1])[sel]),
            float((ref[1] * w)[sel].sum() / (ref[1] * w)[vis].sum()))


def _worst_gen():
    return max(_dev(d, k, o) for d, rows in _COLS
               for _, k, _c, _mk, _ls in rows for o, *_ in OBS
               if _dev(d, k, o) is not None)


def _mass_effect(resdir, obs):
    """|FONLL/ZM - 1| over the visible bins: the charm-mass effect itself."""
    ref, f = hist(resdir, REF, obs), hist(resdir, FONLL_REF, obs)
    m = _vis(resdir, obs)
    r = f[1][m] / ref[1][m] - 1.0
    return float(np.max(np.abs(r))), float(np.median(np.abs(r)))


def _nnlo_range(resdir, obs):
    """(min, max) of the NNLO/NLO ratio over the visible bins."""
    ref = hist(resdir, REF, obs)
    nh = hist(resdir, NNLO[resdir][0], obs)
    if nh is None:
        return None
    m = _vis(resdir, obs)
    r = nh[1][m] / ref[1][m]
    return float(r.min()), float(r.max())


def _incl_nnlo_ratio(resdir, obs):
    """The same ratio on the INCLUSIVE distribution, for the contrast."""
    ref = hist(resdir, "yadism_nlo", obs)
    key = ("yadism_nnlo_fonll_damp" if resdir == "results" else "yadism_nnlo")
    nh = hist(resdir, key, obs)
    if ref is None or nh is None:
        return None
    m = _vis(resdir, obs)
    r = nh[1][m] / ref[1][m]
    return float(r.min()), float(r.max())


def _med_ratio(resdir, key, ref, obs):
    """Median ratio of one curve to a chosen reference, visible bins only."""
    m = _vis(resdir, obs)
    a, b = hist(resdir, key, obs), hist(resdir, ref, obs)
    return float(np.median(a[1][m] / b[1][m]))


def _sigma(resdir, key):
    """Integrated sigma_fid, from the SAME fiducial region the column plots.

    Through _suffix(): a bare filename here quoted the Q2 > 4 integrated
    charm correction beside a Q2 > 11 figure, and the two differ by a factor
    of two -- which is the whole point of the raised floor.
    """
    with open(f"{BASE}/{resdir}/histos_{key}{_suffix(resdir, key)}.json") as f:
        return json.load(f)["sigma_fid_pb"]


def _band_mean(resdir, obs, fn=None):
    """Mean half-width of a band over the visible bins."""
    b = mhou(resdir, obs, fn)
    m = _vis(resdir, obs)
    return float(0.5 * (np.asarray(b[1])[m] - np.asarray(b[2])[m]).mean())


def _band_narrowing(resdir):
    """NLO band width / NNLO band width, per panel."""
    return [_band_mean(resdir, o) / _band_mean(resdir, o,
                                               NNLO_MHOU_FILE[resdir])
            for o, *_ in OBS]


def _nnlo_inside(resdir, obs):
    """(bins where the NNLO curve lies inside the NLO band, bins shown)."""
    m = _vis(resdir, obs)
    ref = hist(resdir, REF, obs)
    nr = hist(resdir, NNLO[resdir][0], obs)[1][m] / ref[1][m]
    b = mhou(resdir, obs)
    hi, lo = 1.0 + np.asarray(b[1])[m], 1.0 + np.asarray(b[2])[m]
    return int(((nr <= hi) & (nr >= lo)).sum()), int(len(nr))


def _panel_median_of_gen(resdir, key):
    """Median |ratio - 1| of ONE generator, worst over that column's panels."""
    out = 0.0
    for o, *_ in OBS:
        ref, h = hist(resdir, REF, o), hist(resdir, key, o)
        if h is None:
            continue
        mk = _vis(resdir, o)
        out = max(out, float(np.median(np.abs(h[1][mk] / ref[1][mk] - 1.0))))
    return out


MESSAGE = (
    "<b>Both columns are now cut the same way</b>: Q&sup2; &gt; 4 GeV&sup2; "
    "and W &gt; 5 GeV on tungsten, the charm region of the paper. The "
    "low-Q&sup2; corner, where the neutral-current charm series converges "
    "worst, is shown, "
    "and it is what the NNLO column of the muon side is about."
    "\n\n"
    "<b>Everything that was a per-cent effect on the inclusive figures is a "
    "ten-per-cent effect here.</b> The worst panel median spread of the "
    "matchings against the massless charm reference is 12%, and the worst "
    "single bin is 42% out, at high x on the muon side. POWHEG-RES is the "
    "furthest of the three on the muon column, as on the integrated figure. "
    "Herwig&rsquo;s events whose beam remnant it cannot build are re-showered "
    "on the same hard process and its scattered lepton is never moved to "
    "make room; the few that cannot be built at all sit at high x, which is "
    "where its muon charm falls below the others. The charm mass, the gap "
    "between the massless and FONLL curves, is about 7% in the median muon "
    "bin (up to 13%) and 2&ndash;4% on the neutrino side; POWHEG-V2mc, the "
    "one massive-charm entry, sits below its massless twin and closer to "
    "FONLL in every row."
    "\n\n"
    "<b>The neutral-current charm series does not converge at the bottom of "
    "the region</b>: the muon NNLO charm cross-section is 0.33&ndash;0.97 of "
    "the NLO one across the bins, against 0.88&ndash;1.02 for the inclusive "
    "rate, and &minus;36% integrated in FONLL at both orders. The NLO band "
    "on muon charm is 26&ndash;44% wide on average and narrows only by 1.2 at "
    "NNLO. The charged-current charm series converges: NNLO is "
    "0.88&ndash;0.99 of NLO bin by bin, tracking its own inclusive rate, and "
    "&minus;3% integrated, with an NLO band of 2.5&ndash;6% narrowing by "
    "1.5&ndash;1.7."
    "\n\n"
    "<b>On charm the NLO band does not contain the NNLO result</b> in most "
    "bins: at most six in ten on the muon side, a third on the neutrino "
    "side, and none in the neutrino hadronic-energy row. A seven-point scale "
    "variation is not a faithful uncertainty on charm production here: "
    "scale-varied coefficient functions are combinations of the ones already "
    "present, and cannot anticipate a channel that first appears at the next "
    "order."
)


CLAIMS = [
    {"what": "POWHEG-V2mc follows YADISM in its own scheme (FFNS n_f = 3) "
             "to 3% in every visible neutrino bin",
     "check": lambda: all(float(np.max(np.abs(
         hist("results_nu", "powheg_nu_mc_charmfinal", o)[1][_vis("results_nu", o)]
         / hist("results_nu", FFNS3_REF, o)[1][_vis("results_nu", o)] - 1)))
         < 0.031 for o, *_ in OBS),
     "detail": lambda: ", ".join("%s %.3f" % (o, float(np.max(np.abs(
         hist("results_nu", "powheg_nu_mc_charmfinal", o)[1][_vis("results_nu", o)]
         / hist("results_nu", FFNS3_REF, o)[1][_vis("results_nu", o)] - 1))))
         for o, *_ in OBS)},
    {"what": "both columns are the charm region Q2 > 4, W > 5 GeV on "
             "tungsten, and every file read is the one computed there",
     "check": lambda: (Q2MIN["results"] == Q2MIN["results_nu"] == 4.0
                       and REGION == "q4w5" and TSUF == "_W"
                       and all(hist(d, REF, o) is not None
                               for d, _r in _COLS for o, *_ in OBS)),
     "detail": lambda: f"{REGION}{TSUF}, Q2 > 4 on both"},
    {"what": "the worst panel median spread against the massless charm "
             "reference is 10-20%",
     "check": lambda: 0.10 < _panel_median() < 0.20,
     "detail": lambda: "worst panel median %.2f%%; per panel " % (
         100 * _panel_median()) + ", ".join(
         "%s %s %.1f%%" % (d[:2], o, 100 * _panel_median_of(d, o))
         for d, _r in _COLS for o, *_ in OBS)},
    {"what": "the worst single bin is about 40% out",
     "check": lambda: 0.30 < _worst_gen() < 0.55,
     "detail": lambda: "worst bin %.1f%%" % (100 * _worst_gen())},
    {"what": "the charm-mass effect is about 7% in the median muon bin (up "
             "to 13%) and 2-4% on the neutrino side",
     "check": lambda: (
         all(0.06 < _mass_effect("results", o)[1] < 0.085 for o, *_ in OBS)
         and max(_mass_effect("results", o)[0] for o, *_ in OBS) < 0.15
         and all(0.02 < _mass_effect("results_nu", o)[1] < 0.04
                 for o, *_ in OBS)),
     "detail": lambda: ", ".join(
         "%s %s med %.1f%% max %.1f%%" % (d[:2], o,
                                          100 * _mass_effect(d, o)[1],
                                          100 * _mass_effect(d, o)[0])
         for d, _r in _COLS for o, *_ in OBS)},
    {"what": "POWHEG-V2mc, the one massive-charm entry, sits BELOW its own "
             "massless twin against the massless reference and closer to "
             "FONLL -- the mass shift is in the right direction",
     "check": lambda: all(
         _med_ratio("results_nu", "powheg_nu_mc_charmfinal", REF, o)
         < _med_ratio("results_nu", "powheg_nu_charmfinal", REF, o)
         and (_med_ratio("results_nu", "powheg_nu_mc_charmfinal",
                         FONLL_REF, o)
              > _med_ratio("results_nu", "powheg_nu_mc_charmfinal", REF, o))
         for o, *_ in OBS),
     "detail": lambda: ", ".join(
         "%s V2mc/ZM %.3f vs V2/ZM %.3f, V2mc/FONLL %.3f" % (
             o,
             _med_ratio("results_nu", "powheg_nu_mc_charmfinal", REF, o),
             _med_ratio("results_nu", "powheg_nu_charmfinal", REF, o),
             _med_ratio("results_nu", "powheg_nu_mc_charmfinal",
                        FONLL_REF, o))
         for o, *_ in OBS)},
    {"what": "the MUON charm series does not converge at the bottom of the "
             "region: NNLO is 0.33-0.97 of NLO across the bins, against "
             "0.88-1.02 for the inclusive rate",
     "check": lambda: (
         0.30 < min(_nnlo_range("results", o)[0] for o, *_ in OBS) < 0.40
         and max(_nnlo_range("results", o)[1] for o, *_ in OBS) < 0.99
         and min(_incl_nnlo_ratio("results", o)[0] for o, *_ in OBS) > 0.86),
     "detail": lambda: ", ".join(
         "%s charm %.3f-%.3f incl %.3f-%.3f" % (
             o, *_nnlo_range("results", o), *_incl_nnlo_ratio("results", o))
         for o, *_ in OBS)},
    {"what": "the NEUTRINO charm series converges: NNLO is 0.88-0.99 of NLO "
             "and tracks its own inclusive rate to 6%",
     "check": lambda: (
         min(_nnlo_range("results_nu", o)[0] for o, *_ in OBS) > 0.85
         and all(abs(_nnlo_range("results_nu", o)[0]
                     - _incl_nnlo_ratio("results_nu", o)[0]) < 0.06
                 for o, *_ in OBS)),
     "detail": lambda: ", ".join(
         "%s charm %.3f-%.3f incl %.3f-%.3f" % (
             o, *_nnlo_range("results_nu", o),
             *_incl_nnlo_ratio("results_nu", o))
         for o, *_ in OBS)},
    {"what": "integrated, the muon charm correction in FONLL at both orders "
             "is -36%, and the neutrino one in ZM is -3%",
     "check": lambda: (abs(_sigma("results", "yadism_charm_nnlo_fonll_damp")
                           / _sigma("results", "yadism_charm_nlo_fonll_damp")
                           - 0.638) < 0.02
                       and abs(_sigma("results_nu", "yadism_charm_nnlo")
                               / _sigma("results_nu", "yadism_charm_nlo")
                               - 0.969) < 0.01),
     "detail": lambda: "mu FONLL %.4f, nu ZM %.4f" % (
         _sigma("results", "yadism_charm_nnlo_fonll_damp")
         / _sigma("results", "yadism_charm_nlo_fonll_damp"),
         _sigma("results_nu", "yadism_charm_nnlo")
         / _sigma("results_nu", "yadism_charm_nlo"))},
    {"what": "the NNLO curve is FONLL on the muon side and ZM on the "
             "neutrino side, where the massive charged-current coefficient "
             "functions at NNLO are known but not implemented in YADISM",
     "check": lambda: ("fonll" in NNLO["results"][0]
                       and "fonll" not in NNLO["results_nu"][0]
                       and all(hist(d, NNLO[d][0], "Q2") is not None
                               for d, _r in _COLS)),
     "detail": lambda: ", ".join(f"{d}: {NNLO[d][0]}" for d, _r in _COLS)},
    {"what": "the NLO scale band on muon charm is 26-44% wide on average, "
             "against 2.5-6% on the neutrino side",
     "check": lambda: (all(0.22 < _band_mean("results", o) < 0.48
                           for o, *_ in OBS)
                       and all(0.02 < _band_mean("results_nu", o) < 0.07
                               for o, *_ in OBS)),
     "detail": lambda: ", ".join(
         "%s %s %.1f%%" % (d[:2], o, 100 * _band_mean(d, o))
         for d, _r in _COLS for o, *_ in OBS)},
    {"what": "the charm bands narrow from NLO to NNLO by only about 1.2 on "
             "the muon side and 1.5-1.7 on the neutrino side",
     "check": lambda: (all(1.1 < f < 1.35 for f in _band_narrowing("results"))
                       and all(1.4 < f < 1.8 for f in _band_narrowing("results_nu"))),
     "detail": lambda: "; ".join(
         d[:2] + " " + ", ".join("%.2f" % f for f in _band_narrowing(d))
         for d, _r in _COLS)},
    {"what": "ON CHARM THE NLO BAND DOES NOT CONTAIN THE NNLO RESULT in most "
             "bins -- at most six in ten on the muon side, a third on the neutrino "
             "side, none in the neutrino hadronic-energy row",
     "check": lambda: (all(_nnlo_inside("results", o)[0] <= 0.62 * _nnlo_inside("results", o)[1]
                           for o, *_ in OBS)
                       and all(_nnlo_inside("results_nu", o)[0] <= 0.4 * _nnlo_inside("results_nu", o)[1]
                               for o, *_ in OBS)
                       and _nnlo_inside("results_nu", "nu")[0] == 0),
     "detail": lambda: ", ".join(
         "%s %s %d/%d" % (d[:2], o, *_nnlo_inside(d, o))
         for d, _r in _COLS for o, *_ in OBS)},
    # THE HERWIG HIGH-x RESIDUAL (TODO item 6, 2026-09-19), stated in the
    # caption: with remnant retries and DISRemnantOption NoLepton the events
    # Herwig can only build by moving the lepton are dropped after 100 tries;
    # they sit in the top x bins.  Pinned so the caption cannot drift.
    {"what": "Herwig's muon charm falls 10-40% short of the ZM reference in "
             "every bin above x = 0.24, bins holding about 2% of the charm rate",
     "check": lambda: (lambda r, sh: all(0.55 < v < 0.9 for v in r)
                       and 0.015 < sh < 0.03)(*_herwig_highx()),
     "detail": lambda: "ratios %s, rate share %.4f" % (
         " ".join("%.3f" % v for v in _herwig_highx()[0]), _herwig_highx()[1])},
    {"what": "POWHEG-RES is the furthest of the three from the muon charm "
             "reference (Herwig was, until its remnant failures were "
             "re-showered instead of replaced, 2026-09-19)",
     "check": lambda: (
         _panel_median_of_gen("results", "powheg_charmfinal")
         > max(_panel_median_of_gen("results", k)
               for _l, k, _c, _m, _s in MU_ROWS
               if k != "powheg_charmfinal")),
     "detail": lambda: ", ".join(
         "%s %.1f%%" % (lab.split(" (")[0],
                        100 * _panel_median_of_gen("results", k))
         for lab, k, _c, _m, _s in MU_ROWS)},
    {"what": "each band is in the scheme of the curve it surrounds: NLO in "
             "ZM matching the reference, NNLO in FONLL on the muon side and "
             "ZM on the neutrino side matching the NNLO curve",
     "check": lambda: ("zm" in MHOU_FILE
                       and "fonll" in NNLO_MHOU_FILE["results"]
                       and "fonll" not in NNLO_MHOU_FILE["results_nu"]
                       and all(mhou(d, "Q2", NNLO_MHOU_FILE[d]) is not None
                               for d, _r in _COLS)),
     "detail": lambda: f"NLO {MHOU_FILE}; NNLO " + ", ".join(
         f"{d}: {NNLO_MHOU_FILE[d]}" for d, _r in _COLS)},
    {"what": "every generator on this figure is a MASSLESS-charm "
             "calculation, which is why the ratio denominator is the "
             "massless reference",
     "check": lambda: REF == "yadism_charm_nlo" and "fonll" in FONLL_REF,
     "detail": lambda: f"denominator {REF}, mass curve {FONLL_REF}"},
]


if __name__ == "__main__":
    main()
