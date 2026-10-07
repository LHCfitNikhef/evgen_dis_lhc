#!/usr/bin/env python3
"""Appendix plot 2b: nuclear-PDF effects on the hadron-level distributions
under the FASERnu selection, NEUTRINO DIS ONLY.

>>> ONE CURRENT, BY USER DECISION (2026-10-07). <<<  The figure is an
illustration of how nPDF uncertainties propagate to hadron-level observables,
and the user restricted it to nu_mu CC on tungsten: a deliberate carve-out
from CONVENTIONS.md rule 2b, like the single-inclusive pion study's.  The muon
current's nPDF impact on the leptonic distributions is in Fig. C.1 (ppA2).

The question (F. Kling's reading of App. C, 2026-10-07): Fig. C.1 shows the
nPDF impact on the leptonic DIS variables only; how does it propagate to the
hadronic observables FASERnu actually selects on -- the charged-track
multiplicity, the leading-hadron energy, the lepton-hadron azimuthal
separation?

  * THE OBSERVABLES AND SELECTION are those of Fig. 6.1 (pp08): N_ch, E_lead,
    Delta phi and p_l' at 1 TeV, per nucleon on tungsten, in Q2 > 4 GeV2,
    W > 3 GeV plus the FASERnu emulsion cuts (selection q4w3_faser_e), with
    pp08's windows and rebinning.
  * THE CALCULATION: the NLO-matched POWHEG-V2 + Pythia 8 neutrino sample of
    Fig. 6.1, reweighted event by event with every member of nNNPDF3.0, EPPS21 and
    nCTEQ15HQ and joined to the showered events by lhe_index
    (analysis/npdf_hadron.py, powheg/reweight/run_reweight_npdf.sh).  The
    proton sample is reweighted with the AVERAGE-nucleon tungsten grids, which
    makes it a sample of the average nucleon exactly at parton level (see
    npdf_hadron.py on why the neutron sample cannot be used).
  * THE LOWER PANELS are each set over its OWN free-nucleon baseline, with its
    68% CL band, exactly as Fig. C.1, so the two figures read alike.  The
    upper panels show the free-nucleon NNPDF4.0 prediction (normalised to the
    published tungsten sigma_fid of the selection) and the three nPDF
    predictions with their bands.
  * MONTE CARLO ERROR: every ratio is two weightings of the SAME events, so the
    statistical error of the ratio is far below that of either curve; it is
    computed per bin from the per-event weight pairs and checked in CLAIMS.

GENIE IS ABSENT, under rule 1b's paper-plot carve-out: GENIE cannot be
reweighted per PDF member, and the figure is about the PDF, not a generator.
Its hadron-level comparison is Fig. 6.1.

Inputs: results{,_nu}/npdf_hadron_q4w3_faser_e.json (analysis/npdf_hadron.py).
Usage:  analysis/paper_plots/ppA2b_npdf_hadron.py
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

PAPER_SECTION = "appendix"
SLUG = "ppA2b_npdf_hadron"
V1_NUMBER = "A2b"
# IN THE PAPER as Fig. C.2 (user approved, 2026-10-07).  Neutrino DIS only.
IN_PAPER = True
TITLE = "Impact of nuclear PDFs on hadron-level distributions"
OUTPUT = "pp_npdf_hadron.png"
RESULTS = "results_nu"
SURFACE = "#ffffff"
FILE = f"{BASE}/results_nu/npdf_hadron_q4w3_faser_e.json"

NUCLEAR = ["nNNPDF30_nlo_as_0118_A184_Z74", "EPPS21nlo_CT18Anlo_W184",
           "nCTEQ15HQ_FullNuc_184_74"]
LABEL = {"base": "NNPDF4.0, free nucleons",
         "nNNPDF30_nlo_as_0118_A184_Z74": "nNNPDF3.0",
         "EPPS21nlo_CT18Anlo_W184": "EPPS21",
         "nCTEQ15HQ_FullNuc_184_74": "nCTEQ15HQ"}
# the colours and line styles of Fig. C.1 (ppA2), so the two figures read alike
COLOUR = {"base": "#111111", "nNNPDF30_nlo_as_0118_A184_Z74": "#0072b2",
          "EPPS21nlo_CT18Anlo_W184": "#d55e00",
          "nCTEQ15HQ_FullNuc_184_74": "#009e73"}
LS = {"base": "-", "nNNPDF30_nlo_as_0118_A184_Z74": "--",
      "EPPS21nlo_CT18Anlo_W184": "-.",
      "nCTEQ15HQ_FullNuc_184_74": (0, (5.0, 1.3, 1.0, 1.3, 1.0, 1.3))}

OBS = [("nch05", r"$N_{\rm ch}$  ($E > 1$ GeV, $\tan\theta < 0.5$)",
        "linear", "linear", r"d\sigma/dN_{\rm ch}"),
       ("Elead", r"$E_{\rm lead}$  [GeV]", "log", "log", r"d\sigma/dE_{\rm lead}"),
       ("dphix", r"$\Delta\phi(\ell', \sum h^\pm)$  [rad]", "linear", "linear",
        r"d\sigma/d\Delta\phi"),
       ("Emu", r"$p_{\ell'}$  [GeV]", "linear", "linear", r"d\sigma/dp_{\ell'}")]
# the windows and rebinning of Fig. 6.1 (pp08), copied rather than imported so
# that a change there is a decision here too
XLIM = {"nch05": (4.5, 14.5), "dphix": (math.pi / 2.0, math.pi),
        "Elead": (10.0, 900.0), "Emu": (200.0, 1000.0)}
# Delta phi is merged FOUR to one here, not pp08's two: below 2.6 rad each
# of pp08's bins holds under 1% of the rate and the ratio's MC error there
# (1-3%) would draw as wiggles comparable to the nuclear effect.
REBIN = {"dphix": 4, "Emu": 3}
# THE DISPLAYED WINDOWS (user, 2026-10-07): p_l' < 900 GeV, Delta phi > 2.4,
# 20 < E_lead < 500 GeV.  XLIM above still anchors the rebinning, so the bins
# are those of Fig. 6.1; a bin is drawn (and enters the checked numbers) if it
# overlaps the window, and the axis clips it there.
AXLIM = {"nch05": (4.5, 14.5), "dphix": (2.4, math.pi),
         "Elead": (20.0, 500.0), "Emu": (200.0, 900.0)}
# BINS HOLDING LESS THAN 0.6% OF THE RATE ARE NOT DRAWN (Fig. C.1, ppA2,
# uses 0.3% on an analytic calculation): below it the ratio's MC error
# reaches 1-5% and the replica spread of nNNPDF3.0 there doubles from
# sampling alone (8% half-width in the 13 GeV E_lead bin, 0.5% of the rate).
MIN_SHARE = 6e-3

CAPTION = (
    "Illustration of the propagation of nuclear-PDF uncertainties to "
    "hadron-level observables: the distributions of Fig. 6.1 for "
    "muon-neutrino charged-current scattering at 1 TeV, per nucleon of "
    "tungsten under the FASER&nu; selection -- charged-hadron multiplicity, "
    "leading charged-hadron energy, lepton-hadron azimuthal separation and "
    "outgoing-lepton momentum. POWHEG-V2 matched to Pythia 8, reweighted "
    "event by event to NNPDF4.0 on free nucleons and to nNNPDF3.0, EPPS21 "
    "and nCTEQ15HQ for tungsten. Beneath each distribution, the ratio of "
    "each set to its own free-nucleon baseline with its 68% CL band, as in "
    "Fig. C.1.")

MESSAGE = """<b>Illustration: nuclear PDFs barely change the shape of the hadronic final state.</b> Under the FASER&nu; selection at 1 TeV the nuclear modification of the neutrino charged-current rate, each set against its own free-nucleon baseline, is -0.0% (nNNPDF3.0), -1.2% (EPPS21) and +1.1% (nCTEQ15HQ), with 68% CL bands of 4.6%, +3.4/-2.7% and 2.3%: the same picture as the leptonic distributions of Fig. C.1, small shifts well inside bands of a few per cent. Differentially the ratios stay within 3% of unity in N_ch, E_lead and p_l' and within 4.4% in Delta phi, and are nearly flat, so the nPDF uncertainty is close to a normalisation; it is small compared with the spread between generators of Fig. 6.1, up to 13% in the mean multiplicity, which is set by the shower and hadronisation."""


def _load():
    with open(FILE) as f:
        return json.load(f)


def _merge_sum(edges, s, n, lo):
    """Merge n adjacent bins of a bin CONTENT array (not a density), after
    trimming to the window start, as pp08 does."""
    e, s = np.asarray(edges, float), np.asarray(s, float)
    if lo is not None:
        k = int(np.searchsorted(e, lo - 1e-9))
        e, s = e[k:], s[..., k:]
    m = (s.shape[-1] // n) * n
    s = s[..., :m].reshape(*s.shape[:-1], -1, n).sum(axis=-1)
    return e[:m + 1:n], s


def curves(obs):
    """{name: (edges, central, lo, hi)} of bin contents, and the ratios
    {name: (ratio, ratio_lo, ratio_hi, mc_err)} to each set's own baseline,
    rebinned BEFORE any ratio is taken."""
    d = _load()["obs"][obs]
    e0 = d["edges"]
    n = REBIN.get(obs, 1)
    lo = XLIM.get(obs, (None, None))[0] if n > 1 else None
    out, rat = {}, {}
    e, b = _merge_sum(e0, d["baseline"], n, lo)
    out["base"] = (e, b, b, b)
    for name in NUCLEAR:
        x = d[name]
        c = np.asarray(x["central"]); p = np.asarray(x["errplus"]); m = np.asarray(x["errminus"])
        ob = np.asarray(x["own_baseline"])
        mc = np.asarray(x["ratio_mc_error"]) * c        # absolute MC error of c
        _e, c2 = _merge_sum(e0, c, n, lo)
        # band edges add linearly when bins merge: members move bins together
        _e, p2 = _merge_sum(e0, p, n, lo)
        _e, m2 = _merge_sum(e0, m, n, lo)
        _e, ob2 = _merge_sum(e0, ob, n, lo)
        _e, mc2 = _merge_sum(e0, np.nan_to_num(mc) ** 2, n, lo)
        out[name] = (e, c2, c2 - m2, c2 + p2)
        with np.errstate(divide="ignore", invalid="ignore"):
            rat[name] = (c2 / ob2, (c2 - m2) / ob2, (c2 + p2) / ob2,
                         np.sqrt(mc2) / ob2)
    return out, rat


def _vis(edges, obs, content=None):
    """Bins inside the window and, given the baseline contents, holding at
    least MIN_SHARE of the rate inside it."""
    lo, hi = XLIM.get(obs, (None, None))
    keep = np.ones(len(edges) - 1, dtype=bool)
    if lo is not None:
        keep &= edges[:-1] >= lo - 1e-9
    if hi is not None:
        keep &= edges[1:] <= hi + 1e-9
    alo, ahi = AXLIM[obs]
    keep &= (edges[1:] > alo + 1e-9) & (edges[:-1] < ahi - 1e-9)
    if content is not None:
        share = np.asarray(content) / np.asarray(content)[keep].sum()
        keep &= share >= MIN_SHARE
    return keep


def _integ(name):
    return _load()["integrated"][name]


def _ratio_range(obs):
    """(min, max) of every drawn ratio curve's central value."""
    _o, rat = curves(obs)
    e = _o["base"][0]
    v = _vis(e, obs, _o["base"][1])
    vals = np.concatenate([rat[n][0][v] for n in NUCLEAR])
    return float(np.nanmin(vals)), float(np.nanmax(vals))


def _worst_mc():
    """Largest MC error of a drawn ratio relative to its 68% half-band."""
    w = 0.0
    for obs, *_x in OBS:
        _o, rat = curves(obs)
        v = _vis(_o["base"][0], obs, _o["base"][1])
        for n in NUCLEAR:
            r, lo, hi, mc = rat[n]
            half = 0.5 * (hi - lo)
            w = max(w, float(np.nanmax((mc / half)[v])))
    return w


def _worst_mc_abs():
    w = 0.0
    for obs, *_x in OBS:
        _o, rat = curves(obs)
        v = _vis(_o["base"][0], obs, _o["base"][1])
        for n in NUCLEAR:
            w = max(w, float(np.nanmax(rat[n][3][v])))
    return w


def _own(name):
    return 100 * (_integ(name)["to_own_baseline"] - 1.0)


def _half(name):
    v = _integ(name)
    return 100 * v["errplus"], 100 * v["errminus"]


def _maxdev(obs):
    """max |ratio - 1| over the drawn bins and the three sets, in per cent."""
    cv, rat = curves(obs)
    v = _vis(cv["base"][0], obs, cv["base"][1])
    return 100 * max(float(np.nanmax(np.abs(rat[n][0][v] - 1.0))) for n in NUCLEAR)


CLAIMS = [
    {"what": "the join of weights to showered events closes (closure weight "
             "= showered weight to 1e-4)",
     "check": lambda: _load()["join_closure"] < 1e-4,
     "detail": lambda: f"{_load()['join_closure']:.1e}"},
    {"what": "the input is the FASERnu selection q4w3_faser_e at 1 TeV, "
             "neutrino current",
     "check": lambda: (_load()["selection"] == "q4w3_faser_e"
                       and _load()["energy_gev"] == 1000.0
                       and _load()["current"] == "nu"),
     "detail": lambda: f"{_load()['selection']}, {_load()['n_events']} events"},
    {"what": "integrated nuclear modification against each set's own "
             "baseline: -0.0% (nNNPDF3.0), -1.2% (EPPS21), +1.1% (nCTEQ15HQ)",
     "check": lambda: all(abs(_own(n) - q) < 0.06 for n, q in
                          zip(NUCLEAR, (-0.04, -1.2, 1.1))),
     "detail": lambda: ", ".join(f"{_own(n):+.2f}%" for n in NUCLEAR)},
    {"what": "68% CL bands 4.6%, +3.4/-2.7%, 2.3%",
     "check": lambda: all(abs(_half(n)[0] - a) < 0.06 and abs(_half(n)[1] - b) < 0.06
                          for n, a, b in zip(NUCLEAR, (4.6, 3.4, 2.3), (4.6, 2.7, 2.3))),
     "detail": lambda: ", ".join("+%.2f/-%.2f%%" % _half(n) for n in NUCLEAR)},
    {"what": "every drawn ratio within 3% of unity in N_ch, E_lead and p_l', "
             "within 4.4% in Delta phi",
     "check": lambda: (max(_maxdev(o) for o in ("nch05", "Elead", "Emu")) < 3.0
                       and _maxdev("dphix") < 4.45),
     "detail": lambda: ", ".join(f"{o} {_maxdev(o):.2f}%" for o, *_x in OBS)},
    {"what": "the ratios' Monte Carlo error is below 1% and below 30% of the "
             "68% half-band in every drawn bin",
     "check": lambda: _worst_mc_abs() < 0.01 and _worst_mc() < 0.30,
     "detail": lambda: f"abs {_worst_mc_abs():.4f}, relative to band {_worst_mc():.2f}"},
    {"what": "the generator spread of Fig. 6.1 in the mean neutrino N_ch is "
             "13% (GENIE 7.26 to Herwig 8.24)",
     "check": lambda: abs(8.24 / 7.26 - 1.13) < 0.01,
     "detail": lambda: "quoted from pp08's claims, not re-derived here"},
]


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    div, unit = 1.0, "pb"
    fig = plt.figure(figsize=(11.5, 10.6))
    outer = fig.add_gridspec(2, 2, wspace=0.28, hspace=0.30, left=0.09,
                             right=0.99, top=0.885, bottom=0.07)
    for c, (obs, xlab, xs, ys, ylab) in enumerate(OBS):
        inner = outer[c // 2, c % 2].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                                 hspace=0.0)
        ax = fig.add_subplot(inner[0])
        ar = fig.add_subplot(inner[1], sharex=ax)
        cv, rat = curves(obs)
        e = cv["base"][0]
        w = np.diff(e)
        vis = _vis(e, obs, cv["base"][1])
        for name in ["base"] + NUCLEAR:
            _e, cen, lo, hi = cv[name]
            col = COLOUR[name]
            ax.stairs(np.where(vis, cen / w / div, np.nan), e, color=col,
                      lw=1.6, ls=LS[name], baseline=None, zorder=5)
            if name == "base":
                continue
            ax.stairs(np.where(vis, hi / w / div, np.nan), e,
                      baseline=np.where(vis, lo / w / div, np.nan),
                      fill=True, color=col, alpha=0.16, lw=0, zorder=2)
            rr, rlo, rhi, _mc = rat[name]
            ar.stairs(np.where(vis, rr, np.nan), e, color=col, lw=1.6,
                      ls=LS[name], baseline=None, zorder=5)
            ar.stairs(np.where(vis, rhi, np.nan), e,
                      baseline=np.where(vis, rlo, np.nan), fill=True,
                      color=col, alpha=0.16, lw=0, zorder=2)
        ar.axhline(1.0, color="#666666", lw=0.8, ls=":")
        ax.set_xscale(xs); ax.set_yscale(ys); ar.set_xscale(xs)
        ax.set_xlim(*AXLIM[obs])
        if obs == "Elead":
            # plain labels on the short log axis: matplotlib's minor labels
            # (2 x 10^1, 3 x 10^1, ...) collide at this panel width
            ar.xaxis.set_major_locator(ticker.FixedLocator([20, 50, 100, 200, 500]))
            ar.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
            ar.xaxis.set_minor_formatter(ticker.NullFormatter())
        b = (cv["base"][1] / w / div)[vis]
        if ys == "log":
            pos = b[b > 0]
            ax.set_ylim(0.3 * pos.min(), 5.0 * pos.max())
        else:
            ax.set_ylim(0.0, 1.25 * b.max())
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune="lower"))
        vals = np.concatenate([np.concatenate([rat[n][1][vis], rat[n][2][vis]])
                               for n in NUCLEAR])
        lo_r, hi_r = np.nanmin(vals), np.nanmax(vals)
        pad = 0.10 * (hi_r - lo_r)
        ar.set_ylim(min(lo_r - pad, 0.985), max(hi_r + pad, 1.015))
        ax.set_ylabel(tex(rf"${ylab}$  [{unit}]"), fontsize=plotstyle.FS_YLABEL - 3)
        if c % 2 == 0:
            ar.set_ylabel(tex("ratio to own") + "\n" + tex("free nucleons"),
                          fontsize=plotstyle.FS_YLABEL - 5)
        ar.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 1)
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(ar)
    hs = [Line2D([], [], color=COLOUR["base"], lw=1.6, ls=LS["base"])]
    hs += [(Patch(color=COLOUR[n], alpha=0.16, lw=0),
            Line2D([], [], color=COLOUR[n], lw=1.6, ls=LS[n])) for n in NUCLEAR]
    labs = [tex(LABEL[n]) for n in ["base"] + NUCLEAR]
    fig.legend(hs, labs, loc="upper center", bbox_to_anchor=(0.54, 0.955),
               ncol=4, frameon=False, handlelength=3.0,
               fontsize=plotstyle.FS_LEGEND + 1)
    fig.suptitle(tex(r"Nuclear PDF corrections, neutrino DIS, FASER$\nu$ "
                     r"selection, $E_{\nu} = 1$ TeV"), y=0.995,
                 fontsize=plotstyle.FS_SUPTITLE)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
