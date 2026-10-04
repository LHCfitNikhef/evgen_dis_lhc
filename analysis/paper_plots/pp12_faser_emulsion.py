#!/usr/bin/env python3
"""Paper plot 12: our predictions against the FASER$\nu$ emulsion
measurement, with the samples.

The redo of the earlier figure for the "Paper plots"
tab.  pp12 STAYS AT FASER'S OWN CUTS (user, 2026-09-14) -- no Q2 window
beyond the ladders' generation floor, no y window, the lepton counted among
the tracks, flux-folded over neutrinos and antineutrinos on tungsten and
normalised to the 33 and 7 observed candidates -- so its region is not the 
benchmark region.  What changed is the samples:

  * POWHEG-V2 is re-showered with LesHouches:matchInOut = off, the charm
    fix (Pythia had dropped ~2% of the events, nearly all charm).
  * THE NEUTRON LADDER POINTS ARE GENUINE NEUTRONS (user: "Regenerate n
    points").  earlier ran POWHEG-V2 and Sherpa on a proton beam carrying the
    isospin-mirrored PDF -- right cross-section, proton remnant, one unit of
    charge off in the hadronic system that the track multiplicity sees.  
    regenerates them on a 2112 beam (tools/faser_emulsion_ladder.sh);
    GENIE always used the neutron target.  The multiplicity shapes move by up
    to 9% per bin.
  * THE FIXED-ENERGY SCALE BAND is imported from the  1 TeV Tier E join
    (mhou_hadron_q4w3_faser_e_W.json); the energy-remix piece is unchanged.

Inputs: results_nu/faser_emulsion_shapes.json
(analysis/faser_emulsion_shapes.py combine).

Usage: analysis/paper_plots/pp12_faser_emulsion.py
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SLUG = "pp12_faser_emulsion"
V1_NUMBER = 12
IN_PAPER = True
TITLE = "Comparison with the FASER$\nu$ emulsion measurement"
OUTPUT = "pp_faser_emulsion.png"
RESULTS = "results_nu"

SRC = f"{BASE}/results_nu/faser_emulsion_shapes.json"
OUT = f"{BASE}/results_nu/{OUTPUT}"
SURFACE = "#ffffff"

# The exposure behind the measured points: 9.5 fb^-1 of 2022 running on the
# 681.1 kg subset of the target, CERN-FASER-CONF-2026-002.  It replaces the
# note's "Preliminary" tag on this figure at the user's request (2026-09-09).
LUMI = r"FASER, 9.5 fb$^{-1}$"

ROWS = [("ntracks", r"$N_{\rm tr}$", 1.0),
        ("tanlep", r"$\tan\theta_{\ell'}$", 1.0),
        ("plep", r"$p_{\ell'}$  [GeV]", 1.0),
        ("dphi", r"$\Delta\phi$  [deg]", 180.0 / math.pi)]
# nu_e ON THE LEFT, nu_mu ON THE RIGHT (user, 2026-09-30).
COLS = [("nue", r"$\nu_e + \bar\nu_e$ CC, 7 candidates"),
        ("numu", r"$\nu_\mu + \bar\nu_\mu$ CC, 33 candidates")]
# THE SAME THREE GENERATORS, AND THE SAME COLOURS, AS PAPER PLOT 8 -- one
# reader looking at both pages should not have to relearn which curve is
# which.  Sherpa is here by user instruction (2026-09-08: "let's also add
# SHERPA NLO, I understand this requires regeneration but this is an important
# plot"); its samples had to be regenerated because the production card cuts
# 0.2 < y < 0.9 AT GENERATION and FASER imposes no y window at all.
GENS = [("powheg_v2", "POWHEG-V2", "#0072b2", "-"),
        ("genie", "GENIE (GRV98LO)", "#d55e00", "--"),
        ("sherpa_nlo", "Sherpa", "#009e73", "-."),
        # Herwig (user, 2026-10-04), its FASER ladder (tools/herwig_faser_ladder.sh)
        ("herwig_nlo", "Herwig", "#7b3294", ":")]


CAPTION = (
    "The four distributions FASER measures in its emulsion detector, against "
    "the predictions of this benchmark, for muon (left) and electron (right) "
    "neutrino charged-current candidates. Points: the 33 and 7 candidates of "
    "CERN-FASER-CONF-2026-002, digitised from its Figs. 9 and 10, with "
    "Poisson errors. Curves: POWHEG-V2, GENIE in the default tune, and Sherpa "
    "and Herwig at next-to-leading order, folded over the neutrino and antineutrino "
    "fluxes on tungsten inside FASER&rsquo;s own selection &mdash; no "
    "Q&sup2; cut, no y window, the lepton counted among the tracks &mdash; "
    "and normalised to the observed count, so what is compared is the shape. "
    "Neutron targets are genuine neutrons in every generator, and POWHEG-V2 "
    "is showered without Pythia&rsquo;s charm-dropping LHE matching. "
    "The band is POWHEG-V2&rsquo;s seven-point scale envelope on the shape. "
    "Lower panels: GENIE, Sherpa and Herwig as a ratio to POWHEG-V2.")

MESSAGE = """The four predictions separate exactly where the hadronic final state enters, and nowhere else. The lepton momentum and the lepton angle are the same in all of them: bin by bin, GENIE, Sherpa and Herwig sit within 4.0%, 2.3% and 3.8% of POWHEG-V2 on the momentum, and within 4.8%, 7.1% and 5.5% on the angle. That is the inclusive agreement of the earlier sections surviving the fold over the flux, and it is what one should expect of an observable a structure-function calculation could have predicted.

The two hadronic observables are a different matter. On the track multiplicity GENIE places more of its rate at low multiplicity and the three matched calculations more in the tail; they cross GENIE between nine and ten tracks, and above twelve POWHEG-V2 is 2.7 times it, Sherpa 3.7 and Herwig 4.0. The azimuth between the lepton and the summed hadrons &mdash; hadronic despite being an angle &mdash; has GENIE 11% above POWHEG-V2 well below the back-to-back peak and a third below it just short of the peak.

Missing higher orders are not what limits this, and the band says so directly: over the bins that carry the rate the seven-point envelope is at most 1.1% on the multiplicity and 2.6% on the lepton momentum, against a generator-to-generator difference on the multiplicity that reaches a factor of three and a Poisson error of twenty per cent at best on 33 candidates. On the hadronic final state it is the modelling that dominates, not the perturbative order.

With 33 and 7 candidates the data cannot arbitrate: the region where the predictions differ by more than a factor of two holds one event. What the comparison shows is where a larger sample would bite. The multiplicity is not an incidental observable for this measurement &mdash; the vertex selection is a cut on it, and the hadronisation uncertainty on the selection efficiency is among the largest entries in FASER's systematic budget &mdash; so a measurement of the multiplicity itself would constrain the modelling that uncertainty stands for."""


_D = {}


def _d():
    if not _D:
        with open(SRC) as f:
            _D.update(json.load(f))
    return _D


def _shape(gen, fl, obs):
    return np.asarray(_d()["generators"][gen][fl]["shapes"][obs], dtype=float)


def _edges(fl, obs):
    return np.asarray(_d()["generators"]["powheg_v2"][fl]["bins"][obs],
                      dtype=float)


def _spread(fl, obs, gen):
    """max |gen/POWHEG-V2 - 1| over the bins that carry the rate.

    A BIN HOLDING A PER MILLE OF THE SAMPLE gives a ratio that is mostly its
    own Monte Carlo error, so the threshold is one per cent of the rate --
    the same one the band claims use.
    """
    p, g = _shape("powheg_v2", fl, obs), _shape(gen, fl, obs)
    m = (p > 0.01 * p.sum()) & (g > 0)
    return float(np.max(np.abs(g[m] / p[m] - 1.0)))


def _above(fl, x, gen):
    """The multiplicity ratio to GENIE, integrated above x tracks."""
    e = _edges(fl, "ntracks")
    m = e[:-1] >= x
    return float(_shape(gen, fl, "ntracks")[m].sum()
                 / _shape("genie", fl, "ntracks")[m].sum())


def _cross(fl):
    """The multiplicity bin in which POWHEG-V2 crosses GENIE, as (lo, hi)."""
    e = _edges(fl, "ntracks")
    r = _shape("powheg_v2", fl, "ntracks") / _shape("genie", fl, "ntracks")
    i = int(np.argmax(r > 1.0))
    return float(e[i]), float(e[i + 1])


def _band(fl, obs):
    """The largest scale band over the bins carrying the rate, in per cent."""
    rel = _d()["generators"]["powheg_v2"][fl]["mhou_rel"][obs]
    p = _shape("powheg_v2", fl, obs)
    m = p > 0.01 * p.sum()
    return 100.0 * float(max(np.max(np.abs(np.asarray(rel[0])[m])),
                             np.max(np.abs(np.asarray(rel[1])[m]))))


def _dphi_ratio(fl, lo_deg, hi_deg):
    """GENIE over POWHEG-V2, averaged over an azimuthal range."""
    e = _edges(fl, "dphi") * 180.0 / math.pi
    m = (e[:-1] >= lo_deg - 1e-6) & (e[1:] <= hi_deg + 1e-6)
    return float(_shape("genie", fl, "dphi")[m].sum()
                 / _shape("powheg_v2", fl, "dphi")[m].sum())


def private_available():
    """The digitised FASER data (CERN-FASER-CONF-2026-002) are not part of the
    public release; claims marked private are skipped without them."""
    return "data" in _d()


CLAIMS = [
    {"what": "the digitised measurement closes: each of the eight panels "
             "sums to the published candidate count, 33 and 7",
     "private": True,
     "check": lambda: all(
         abs(sum(_d()["data"][fl][obs]["counts"]) - _d()["observed"][fl]) < 1e-9
         for fl in ("numu", "nue") for obs in ("ntracks", "tanlep", "plep",
                                               "dphi")),
     "detail": lambda: ", ".join(
         f"{fl} {obs} {sum(_d()['data'][fl][obs]['counts']):.0f}"
         for fl in ("numu", "nue") for obs in ("ntracks", "tanlep", "plep",
                                               "dphi"))},
    {"what": "it is FASER's selection and not the benchmark's: no Q2 cut and "
             "no y window, the lepton counted among the tracks",
     "check": lambda: ("q2_min" not in _d()["cuts"]
                       and "y_min" not in _d()["cuts"]
                       and _d()["cuts"]["n_tracks_min"] == 5
                       and _d()["cuts"]["lep_p_min"] == 200.0),
     "detail": lambda: ", ".join(f"{k}={v}" for k, v in
                                 sorted(_d()["cuts"].items()))},
    {"what": "the leptonic observables agree: GENIE and Sherpa within 4.0% "
             "and 2.3% of POWHEG-V2 on the lepton momentum, and within 4.8% "
             "and 7.1% on the lepton angle",
     "check": lambda: (abs(_spread("numu", "plep", "genie") - 0.040) < 0.002
                       and abs(_spread("numu", "plep", "sherpa_nlo") - 0.023) < 0.002
                       and abs(_spread("numu", "tanlep", "genie") - 0.048) < 0.002
                       and abs(_spread("numu", "tanlep", "sherpa_nlo") - 0.071) < 0.002),
     "detail": lambda: ", ".join(
         f"{obs} {gen} {100 * _spread('numu', obs, gen):.1f}%"
         for obs in ("plep", "tanlep") for gen in ("genie", "sherpa_nlo"))},
    {"what": "Herwig's leptonic observables agree as well: within 3.8% of "
             "POWHEG-V2 on the lepton momentum and 5.5% on the lepton angle",
     "check": lambda: (abs(_spread("numu", "plep", "herwig_nlo") - 0.038) < 0.002
                       and abs(_spread("numu", "tanlep", "herwig_nlo") - 0.055) < 0.002),
     "detail": lambda: (f"plep {100 * _spread('numu', 'plep', 'herwig_nlo'):.1f}%, "
                        f"tanlep {100 * _spread('numu', 'tanlep', 'herwig_nlo'):.1f}%")},
    {"what": "above twelve tracks Herwig is 4.0 times GENIE, the most of "
             "the three matched calculations",
     "check": lambda: (abs(_above("numu", 12.5, "herwig_nlo") - 3.97) < 0.05
                       and _above("numu", 12.5, "herwig_nlo")
                       > max(_above("numu", 12.5, "powheg_v2"),
                             _above("numu", 12.5, "sherpa_nlo"))),
     "detail": lambda: f"Herwig {_above('numu', 12.5, 'herwig_nlo'):.2f}"},
    {"what": "on the multiplicity the matched calculations cross GENIE "
             "between nine and ten tracks, and above twelve POWHEG-V2 is 2.7 "
             "times it and Sherpa 3.7",
     "check": lambda: (_cross("numu") == (8.5, 10.5)
                       and abs(_above("numu", 12.5, "powheg_v2") - 2.74) < 0.05
                       and abs(_above("numu", 12.5, "sherpa_nlo") - 3.67) < 0.05),
     "detail": lambda: (f"crosses in {_cross('numu')}, above 12.5 "
                        f"POWHEG-V2 {_above('numu', 12.5, 'powheg_v2'):.2f}, "
                        f"Sherpa {_above('numu', 12.5, 'sherpa_nlo'):.2f}")},
    {"what": "on the azimuth GENIE is 11% above POWHEG-V2 well below the "
             "back-to-back peak and a third below it just short of the peak",
     "check": lambda: (abs(_dphi_ratio("numu", 90.0, 120.0) - 1.11) < 0.02
                       and abs(1.0 / _dphi_ratio("numu", 160.0, 170.0) - 1.39) < 0.03),
     "detail": lambda: (f"90-120 deg {_dphi_ratio('numu', 90.0, 120.0):.3f}, "
                        f"160-170 deg {_dphi_ratio('numu', 160.0, 170.0):.3f}")},
    {"what": "the scale band is at most 1.1% on the multiplicity and 2.6% on "
             "the lepton momentum, an order of magnitude under the spread "
             "between the generators on the multiplicity",
     "check": lambda: (abs(_band("numu", "ntracks") - 1.14) < 0.15
                       and abs(_band("numu", "plep") - 2.62) < 0.15
                       and _spread("numu", "ntracks", "genie")
                       > 10 * _band("numu", "ntracks") / 100.0),
     "detail": lambda: (f"band N {_band('numu', 'ntracks'):.2f}%, "
                        f"p_lep {_band('numu', 'plep'):.2f}%, "
                        f"GENIE spread on N "
                        f"{100 * _spread('numu', 'ntracks', 'genie'):.0f}%")},
    {"what": "only POWHEG-V2 carries a band: GENIE has no scale weights and "
             "the Sherpa and Herwig ladders were not reweighted",
     "check": lambda: (["mhou_rel" in _d()["generators"][k]["numu"]
                        for k in ("powheg_v2", "genie", "sherpa_nlo", "herwig_nlo")]
                       == [True, False, False, False]),
     "detail": lambda: ", ".join(
         f"{k} {'band' if 'mhou_rel' in _d()['generators'][k]['numu'] else 'none'}"
         for k in ("powheg_v2", "genie", "sherpa_nlo", "herwig_nlo"))},
]


def _drawband(a, edges, band, colour, label=None):
    """The scale band as a step region: filled, and edged so it is legible.

    A few per cent is thinner than the curve that runs through it, so the
    fill alone disappears; the two thin edges are what make it readable.
    """
    lo = np.append(band[0], band[0][-1])
    hi = np.append(band[1], band[1][-1])
    a.fill_between(edges, lo, hi, step="post", color=colour, alpha=0.32,
                   lw=0, zorder=2, label=label)
    for e in (lo, hi):
        a.step(edges, e, where="post", color=colour, lw=0.8, alpha=0.85,
               zorder=3)


def main():
    # IMPORTED HERE, NOT AT MODULE LEVEL: make_report.py imports every script
    # in this directory to read its TITLE, MESSAGE and CLAIMS, and applying
    # the matplotlib style as a side effect of that import is how a report
    # build starts depending on the order its figures were loaded in.
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    if not os.path.exists(SRC):
        sys.exit(f"{SRC} not found -- run "
                 f"analysis/faser_emulsion_shapes.py combine")
    with open(SRC) as f:
        d = json.load(f)
    # >>> A RATIO PANEL UNDER EVERY DISTRIBUTION, TO POWHEG-V2 (user,
    # 2026-09-08 for the panel, 2026-09-09 for the denominator: "let's
    # normalise to POWHEG instead, it is more natural here"). <<<  It was
    # GENIE, on the argument that GENIE is what FASER's own simulation uses;
    # POWHEG-V2 is the better reference because it is the one prediction here
    # that carries an uncertainty, so its band sits AROUND UNITY where it can
    # be read against the spread of the other two rather than riding on a
    # curve that is itself displaced.
    # Same fused geometry as the paper plots: no space between a panel and its
    # ratio, real space between rows.
    fig = plt.figure(figsize=(10.6, 14.6))
    # THE LEGEND SITS ABOVE THE PANELS AND THE COLUMN TITLES SIT ON THEM, so
    # the grid has to start low enough for both.  With three generators the
    # legend grew and began overprinting the titles.
    outer = fig.add_gridspec(4, 2, hspace=0.34, wspace=0.24,
                             top=0.900, bottom=0.045, left=0.085, right=0.985)
    for i, (obs, xlab, xscale) in enumerate(ROWS):
        for j, (fl, coltitle) in enumerate(COLS):
            inner = outer[i, j].subgridspec(2, 1, height_ratios=[1.9, 1.0],
                                            hspace=0.0)
            ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            edges = None
            den = None
            gen_curves = {}
            band = None                     # (lo, hi) absolute, POWHEG-V2
            for key, lab, colour, ls in GENS:
                g = d["generators"].get(key, {}).get(fl)
                if not g:
                    continue
                edges = np.array(g["bins"][obs]) * xscale
                y = np.array(g["shapes"][obs])
                gen_curves[key] = y
                # >>> THE MHOU BAND, ON POWHEG-V2 ONLY (user, 2026-09-08: "I
                # don't see MHOUs in the POWHEG predictions?"). <<<  It is the
                # 7-point envelope, built in two pieces by
                # faser_emulsion_shapes._mhou -- the energy remix exactly, the
                # shape at fixed energy imported from the benchmark's own
                # emulsion-like tier at 1 TeV.  GENIE has no scale weights and
                # Sherpa's ladder was not reweighted, so only one curve
                # carries a band and the caption says which.
                rel = g.get("mhou_rel", {}).get(obs)
                if rel is not None:
                    band = (y * (1.0 + np.array(rel[0])),
                            y * (1.0 + np.array(rel[1])))
                    # IT IS A FEW PER CENT WIDE, so at alpha 0.20 and no edge
                    # it was invisible beside a 1.9 pt curve and read as
                    # absent (user, 2026-09-09: "the MHOU band is missing").
                    # Drawn darker and with its own edges it can be seen for
                    # what it is without being drawn any wider than it is.
                    _drawband(ax, edges, band, colour)
                ax.stairs(y, edges, color=colour, lw=1.9, ls=ls,
                          baseline=None, label=tex(lab), zorder=4)
            den = gen_curves.get("powheg_v2")
            # THE DIGITISED DATA ARE NOT PUBLIC (CONVENTIONS.md: private
            # inputs): without them the predictions are drawn alone
            if "data" in d:
                blk = d["data"][fl][obs]
                de = np.array(blk["edges"], dtype=float)
                cnt = np.array(blk["counts"], dtype=float)
                mid = 0.5 * (de[:-1] + de[1:])
                drawn = cnt > 0
                ax.errorbar(mid[drawn], cnt[drawn], yerr=np.sqrt(cnt[drawn]),
                            xerr=0.5 * np.diff(de)[drawn], fmt="o", ms=4.5,
                            color="#111111", elinewidth=1.3, capsize=0,
                            zorder=6)
            else:
                de, cnt = edges, np.zeros(0)
                drawn = np.zeros(0, dtype=bool)
            # ---------------------------------------------------- the ratio
            if den is not None:
                ok = den > 0
                ratios = []
                rband = None
                if band is not None:
                    b0 = np.where(ok, band[0] / np.where(ok, den, 1.0), np.nan)
                    b1 = np.where(ok, band[1] / np.where(ok, den, 1.0), np.nan)
                    rband = (b0, b1)
                    _drawband(axr, edges, rband, GENS[0][2])
                for key, _lab, colour, ls in GENS:
                    if key == "powheg_v2" or key not in gen_curves:
                        continue
                    with np.errstate(divide="ignore", invalid="ignore"):
                        r = np.where(ok, gen_curves[key] / den, np.nan)
                    ratios.append(r)
                    axr.stairs(r, edges, color=colour, lw=1.9, ls=ls,
                               baseline=None, zorder=4)
                axr.axhline(1.0, color=GENS[0][2], lw=1.4, ls=GENS[0][3],
                            zorder=3)
                # >>> NO DATA HERE (user, 2026-09-09). <<<  The ratio panel
                # is generator against generator, with the scale band; the
                # measurement stays in the panel above it.
                # THE RANGE IS SET BY THE BINS THAT CARRY THE DISTRIBUTION.
                # A bin where POWHEG-V2 has a per cent of its peak gives a
                # ratio that can be anything, and letting it set the scale
                # flattens the bins the comparison is actually about.  The
                # bins themselves are still DRAWN; only the range ignores
                # them.
                carry = ok & (den > 0.03 * float(np.nanmax(den)))
                extra = ([rband[0], rband[1]] if rband is not None else [])
                vals = [v for v in
                        np.concatenate([r[carry & np.isfinite(r)]
                                        for r in ratios]
                                       + [b[carry & np.isfinite(b)]
                                          for b in extra])]
                lo, hi = (min(vals), max(vals)) if vals else (0.0, 2.0)
                # ... BUT A CURVE RUNNING ALONG THE PANEL EDGE READS AS A
                # BROKEN PLOT, and with the data gone the range is no longer
                # stretched by their error bars, so the sparse bins now fall
                # outside far more often.  The range is therefore widened to
                # take them in, by at most 60% of its own span on each side:
                # enough for Delta phi, where the three low bins hold 0.9% of
                # the rate and sit just under, and not enough for the
                # multiplicity tail, where the ratio reaches thirty.
                span = max(hi - lo, 0.1)
                allv = np.concatenate(
                    [r[np.isfinite(r)] for r in ratios]
                    + [b[np.isfinite(b)] for b in extra]) if ratios else []
                if len(allv):
                    lo = max(float(np.min(allv)), lo - 0.6 * span)
                    hi = min(float(np.max(allv)), hi + 0.6 * span)
                pad = 0.10 * max(hi - lo, 0.1)
                axr.set_ylim(max(0.0, lo - pad), hi + pad)
                axr.yaxis.set_major_locator(
                    ticker.MaxNLocator(nbins=4, prune="upper"))
            axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 3)
            ax.set_ylabel(tex("events per bin"),
                          fontsize=plotstyle.FS_YLABEL - 4)
            if j == 0:
                # ON TWO LINES: "Ratio to POWHEG-V2" on one runs the whole
                # height of the ratio panel and collides with the "events
                # per bin" of the panel above it.  (It cannot be shortened
                # to a bare "POWHEG" -- that is the one name this repository
                # never writes, and tools/check_powheg_naming.py enforces it.)
                axr.set_ylabel(tex("Ratio to") + "\n" + tex("POWHEG-V2"),
                               fontsize=plotstyle.FS_YLABEL - 5,
                               linespacing=1.2)
            if edges is None:
                edges = de
            for a_ in (ax, axr):
                a_.set_xlim(edges[0], edges[-1])
            # THE TOP PANEL MUST HOLD THE DATA'S ERROR BARS.  With 33 and 7
            # candidates the Poisson errors are the largest thing on the
            # page, and a range set by the curves alone clips them -- which
            # reads as agreement the statistics do not support.
            top = max(float(np.nanmax(cnt[drawn] + np.sqrt(cnt[drawn])))
                      if drawn.any() else 0.0,
                      max((float(np.nanmax(y)) for y in gen_curves.values()),
                          default=0.0),
                      float(np.nanmax(band[1])) if band is not None else 0.0)
            ax.set_ylim(0.0, 1.08 * top)
            ax.yaxis.set_major_locator(
                ticker.MaxNLocator(nbins=6, prune="lower"))
            plotstyle.ticks(ax, labelbottom=False)
            plotstyle.ticks(axr)
            if i == 0:
                ax.set_title(tex(coltitle),
                             fontsize=plotstyle.FS_PANEL_TITLE - 1)
            # The note's "Preliminary" tag used to sit in this corner; the
            # exposure now rides on the data's own legend entry instead
            # (user, 2026-09-09), which puts it beside the points it
            # qualifies and on every panel rather than on one.
    # >>> FOUR ENTRIES: THE BAND RIDES ON POWHEG-V2's OWN (user, 2026-09-09:
    # "integrate POWHEG MHOU into the POWHEG legend"). <<<  A fifth entry for
    # a band that belongs to one curve, drawn in that curve's colour, told
    # the reader nothing the swatch does not.  The handles are built here
    # rather than harvested from the axes, because a composite swatch --
    # the line drawn over its band -- has no single artist to harvest.
    pw, gn, sh, hw = GENS
    _data_handle = next((c for c in fig.axes[0].containers
                         if c.__class__.__name__ == "ErrorbarContainer"),
                        Line2D([], [], color="#111111", marker="o", ms=4.5,
                               ls="none"))
    hband = Patch(facecolor=pw[2], alpha=0.32, edgecolor=pw[2], lw=0.8)
    handles = [(hband, Line2D([], [], color=pw[2], lw=1.9, ls=pw[3])),
               Line2D([], [], color=gn[2], lw=1.9, ls=gn[3]),
               Line2D([], [], color=sh[2], lw=1.9, ls=sh[3]),
               Line2D([], [], color=hw[2], lw=1.9, ls=hw[3]),
               # the data key carries its error bar, as on pp12b (user,
               # 2026-10-04): harvested from the first panel's errorbar
               _data_handle]
    labels = [tex(pw[1] + plotstyle.MHOU_SUFFIX), tex(gn[1]), tex(sh[1]),
              tex(hw[1]), tex(LUMI)]
    if "data" not in _d():
        handles, labels = handles[:-1], labels[:-1]
    fig.legend(handles, labels, loc="upper center",
               bbox_to_anchor=(0.53, 0.980), ncol=5, frameon=True,
               fontsize=plotstyle.FS_LEGEND + 2,
               handler_map={tuple: HandlerTuple(ndivide=None)})
    fig.suptitle(tex(r"Comparison with FASER$\nu$ data "
                     r"(normalised to observed yields)"),
                 y=0.992, fontsize=plotstyle.FS_SUPTITLE - 2)
    fig.savefig(OUT, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
