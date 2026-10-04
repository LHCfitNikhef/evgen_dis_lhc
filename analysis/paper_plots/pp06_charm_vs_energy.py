#!/usr/bin/env python3
"""Paper plot 6: fiducial CHARM cross-section against beam energy, in the
CHARM region of the final paper.

The redo of the earlier figure for the "Paper plots"
tab: the NLO matchings (plus POWHEG-V2mc on the charged current) against
YADISM charm in ZM-VFNS, with FONLL beside it and NNLO with its band, both
currents.  The conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184, every generator on p and n.
  * THE CHARM REGION, Q2 > 4 GeV2 and W > 5 GeV with no y cut, ON BOTH
    CURRENTS (user, 2026-09-14; selection q4w5).  The earlier muon floor Q2 > 11
    does not carry over.  W > 5 removes the heavy-quark remnant threshold
    below which Pythia and Herwig cannot hadronise NC charm-initiated events
    (POWHEG-RES charm-tagged loss -3.3% at 1 TeV and -7.3% at 400 GeV at
    W > 3, below half a per cent at W > 5; PRODUCTION.md).
  * THE REFERENCES carry target-mass corrections, as every reference does.
  * FIVE ENERGIES, 400-4000 GeV.
  * HERWIG WITH REMNANT RETRIES (2026-09-19).  Until then Herwig read
    0.74-0.98 (mu) / 0.81-0.95 (nu) of the reference, rising with energy:
    HwRemDecayer cannot put the remnant on shell after a SEA quark's forced
    splitting at high x, ThePEG replaced such events by NEW ones, and the
    normalisation handed their (charm-rich) rate to valence d/u.  The samples
    are now generated with EventHandler:MaxEventErrorRetries 100
    (patches/thepeg-2.3.0-event-error-retry.diff) and RemnantDecayer:
    DISRemnantOption NoLepton (never move the lepton), and read 0.97-1.03
    (mu) / 1.00 (nu).

GENIE is absent, as in the earlier production, under the paper-plot carve-out of rule 1b: it has
its own charm figure (paper plot 5).

Usage: analysis/paper_plots/pp06_charm_vs_energy.py
"""
import json
import os
import sys

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

SLUG = "pp06_charm_vs_energy"
V1_NUMBER = 6
IN_PAPER = True
TITLE = "Fiducial NLO charm production cross-section vs beam energy"
OUTPUT = "pp_charm_vs_energy.png"
RESULTS = "results_nu"
CAPTION = ("Fiducial charm-production cross-section per nucleon on tungsten "
           "against beam energy, for Q&sup2; &gt; 4 GeV&sup2; and W &gt; 5 GeV, "
           "muon neutral current (left) and neutrino charged current "
           "(right), for the three NLO matchings, plus POWHEG-V2mc on the "
           "charged current, beside YADISM in the fixed-flavour scheme with "
           "three light flavours and no charm in the proton (FFNS n_f = 3), "
           "the scheme POWHEG-V2mc implements. The reference is YADISM in the massless "
           "(ZM) scheme, because every generator shown computes massless "
           "charm; the FONLL curve is drawn beside it, and the gap between "
           "the two is the effect of the charm mass. Every curve in the "
           "lower panels is divided by the massless reference. Both panels "
           "also carry the NNLO result with its own scale band, as in the "
           "inclusive figure. The references include target-mass "
           "corrections.")

SURFACE = "#ffffff"
# CONSISTENT WITH PAPER PLOT 1 (user, 2026-08-29): the reference is the
# MASSLESS one, because every generator drawn here computes massless charm.
# FONLL goes in the ratio panel, where its distance from unity is the
# charm-mass effect -- 7% on the muon side, 1-4% on the neutrino side, so an
# order of magnitude larger than on the inclusive rate.
REF = "yadism_charm_nlo"
FONLL_REF = "yadism_charm_nlo_fonll_damp"
# The ORDER is dropped from the label (user): the title already
# says NLO and there is no other order on the figure.  The SCHEME
# stays -- CONVENTIONS.md rule 3, and tools/check_yadism_scheme.py
# enforces it: "YADISM" alone does not identify a calculation.
REF_LABEL = "YADISM (ZM)"
# POWHEG-V2mc's own scheme: FFNS with three light flavours and no charm PDF.
FFNS3_REF = "yadism_charm_nlo_ffns3"
FFNS3_LABEL = r"YADISM (FFNS $n_f\!=\!3$)"
# Okabe-Ito reddish purple: separable from POWHEG-V2mc's brown, from FONLL's
# red and from the NNLO purple (#6a3d9a), which is much darker.
FFNS3_COLOUR = "#cc79a7"

# (label, result key, colour, marker, linestyle).
#
# MAXIMALLY DISTINCT COLOURS AND LINE STYLES (user, 2026-08-28).  The three
# were previously three blues/purples, which is unreadable in a printed
# figure and invisible to the ~8% of male readers with a red-green
# deficiency.  These are blue / vermillion / bluish-green from Okabe-Ito, a
# palette built to stay separable under every common form of colour blindness
# -- and the LINE STYLE carries the same information independently, so the
# figure survives being printed in greyscale.
#
# The POWHEG entries keep their variant: POWHEG-RES and POWHEG-V2 are
# different codes and the difference changes results (CONVENTIONS.md rule 3).
MU_ROWS = [
    ("POWHEG-RES", "powheg_charmfinal",                 "#0072b2", "o", "-"),
    ("Herwig",     "herwig_nlo_powheg_full_charmfinal", "#d55e00", "s", "--"),
    ("Sherpa",     "sherpa_charmfinal",                 "#009e73", "^", "-."),
]
NU_ROWS = [
    ("POWHEG-V2",  "powheg_nu_charmfinal",              "#0072b2", "o", "-"),
    ("Herwig",     "herwig_nlo_full_charmfinal",        "#d55e00", "s", "--"),
    ("Sherpa",     "sherpa_nlo_charmfinal",             "#009e73", "^", "-."),
]

# POWHEG-V2mc ON THE NEUTRINO PANEL ONLY (user, 2026-08-29).  It is the one
# entry here that HAS the charm mass -- the only massive-charm calculation
# this benchmark has.  It exists on the CHARGED CURRENT only: the muon NC card
# forbids a non-zero qmass outright.
#
# >>> IT IS NORMALISED TO ZM LIKE EVERYTHING ELSE (user, 2026-08-30). <<<
# It was briefly drawn against FONLL, on the argument that FONLL is the one
# reference it can be compared with like-for-like, and its legend entry said
# so.  That is worse: a ratio panel with two different denominators cannot be
# read by eye, since a reader comparing two curves in it would be comparing
# them to different things.  ONE denominator, named once on the axis.  The
# like-for-like comparison with FONLL is still made -- as a number, in
# MESSAGE, and as a claim -- which is where a statement that cannot be drawn
# without breaking the panel belongs.
# BROWN, not the purple it used to be: the NNLO curve added to this panel is
# purple, and two purples in one ratio panel -- one solid, one dotted -- read
# as the same object at a glance.  Brown is separable from all five other
# colours here and stays separable in greyscale.
NU_MC_ROW = ("POWHEG-V2mc", "powheg_nu_mc_charmfinal",
             "#8c564b", "D", ":")

# The MASSLESS-charm reference, drawn alongside.  See the docstring: every
# generator here computes massless charm, so this is the like-for-like curve
# and the gap between the two analytic lines is the charm-mass effect the
# generators do not have.



# THE Q2 FLOOR IS PER CURRENT (user, 2026-09-01) -- see pp01's note for the
# reason and pp04's for why charm feels it hardest.  Going back is one edit:
# set "results" to 4.0.  Re-derive with tools/make_q2min_muon.sh 11.
Q2MIN = {"results": 4.0, "results_nu": 4.0}
REGION = "q4w5"
TSUF = "_W"
TMC = "_tmc"


def _esuf(e):
    tag = beams.Beams("mu", e).tag
    return "" if e == beams.ANCHOR_ENERGY else f"_{tag}"


def _mhou_name(resdir, fn):
    """An MHOU filename moved to this column's Q2 floor (see pp01)."""
    tag = f"_{REGION}{TSUF}{TMC}"
    if fn.startswith("mhou_sigma_fonll_mu"):
        return fn.replace("mhou_sigma_fonll_mu", "mhou_sigma_fonll_mu" + tag, 1)
    return fn.replace("mhou_sigma_fonll", "mhou_sigma_fonll" + tag, 1)


def sigma(resdir, key, e):
    """(sigma_fid, err) in pb, or None if absent or statistically unusable.

    THE FILE'S OWN STAMP IS CHECKED against Q2MIN, as on every figure here.
    """
    tmc = TMC if key.startswith("yadism") else ""
    p = f"{BASE}/{resdir}/histos_{key}_{REGION}{TSUF}{tmc}{_esuf(e)}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    q2 = (d.get("selection") or {}).get("q2_min")
    if q2 is not None and abs(q2 - Q2MIN[resdir]) > 1e-9:
        raise SystemExit(f"{os.path.relpath(p, BASE)} was computed at "
                         f"Q2 > {q2:g}, but this figure's {resdir} column is "
                         f"Q2 > {Q2MIN[resdir]:g}")
    if d.get("stat_insufficient"):
        return None
    return d.get("sigma_fid_pb"), d.get("sigma_fid_err_pb") or 0.0


def _r(resdir, key, e):
    v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
    return None if not (v and r) else v[0] / r[0]


# ----------------------------------------------------------------- claims --
# Executable, and they re-derive themselves from the same JSONs the figure
# plots.  A claim that compared a hand-typed number to itself would check
# nothing (see tools/check_paper_plots.py).
_ES = list(beams.BENCH_ENERGIES)
_NLO = [("results", ["powheg_charmfinal",
                     "herwig_nlo_powheg_full_charmfinal",
                     "sherpa_charmfinal"]),
        ("results_nu", ["powheg_nu_charmfinal",
                        "herwig_nlo_full_charmfinal",
                        "sherpa_nlo_charmfinal"])]


def _rz(d, k, e):
    """Ratio to the massless reference, which is now REF itself."""
    return _r(d, k, e)


def _nlo_dev():
    """Largest |ratio - 1| over the NLO matchings, both currents."""
    return max(abs(_r(d, k, e) - 1.0)
               for d, ks in _NLO for k in ks for e in _ES)


def _nlo_spread(only=None):
    """Largest max/min among the NLO matchings at a single point
    (optionally on one current only: "results" or "results_nu")."""
    out = 0.0
    for d, ks in _NLO:
        if only and d != only:
            continue
        for e in _ES:
            v = [_r(d, k, e) for k in ks]
            out = max(out, max(v) / min(v) - 1.0)
    return out


MESSAGE = (
    "<b>On charm the matched calculations agree on the neutrino side and "
    "not on the muon side.</b> In the charm region (Q&sup2; &gt; 4 GeV&sup2;, "
    "W &gt; 5 GeV) on tungsten the three neutrino NLO matchings agree to 3% "
    "at every energy: POWHEG-V2 at 1.01&ndash;1.03 of the massless "
    "reference, Herwig at 1.00, Sherpa at 0.98&ndash;1.00. On the muon side "
    "they span up to 17% at a single point, against 2.1% on the inclusive "
    "cross-section of the same region: Herwig sits at 0.97&ndash;1.03 and "
    "Sherpa at 0.93&ndash;1.05, while POWHEG-RES is at 0.90&ndash;0.92."
    "\n\n"
    "<b>Herwig's remnant failures are regenerated on the same hard "
    "process.</b> At these energies Herwig cannot build the beam remnant for "
    "a few per cent of events, almost all of them initiated by a sea quark at "
    "high x, which makes them largely charm events. Replacing such an event "
    "by a new one depleted Herwig's charm by up to 26% at 400 GeV and made "
    "it rise with beam energy. Every Herwig sample here instead re-showers "
    "the failed event on the same hard kinematics until it succeeds, and "
    "never moves the scattered lepton to make room for the remnant."
    "\n\n"
    "<b>The charm mass is a smaller effect than the muon-side spread.</b> "
    "The gap between the massless and the FONLL curves is 3&ndash;7% on the "
    "muon side and 1&ndash;4% on the neutrino side; POWHEG-V2mc, the one "
    "entry with a massive charm, sits 5&ndash;6% below the FONLL curve and "
    "within 1% of YADISM in its own scheme (three light flavours, no charm "
    "in the proton): the gap to FONLL is the charm-initiated scattering "
    "FONLL takes from NNPDF4.0's fitted charm, not a defect of either."
    "\n\n"
    "<b>The neutral-current charm series does not converge at Q&sup2; = 4 "
    "GeV&sup2;.</b> In FONLL at both orders the muon NNLO charm cross-section "
    "is 0.51 of NLO at 400 GeV and 0.74 at 4 TeV, and the NLO seven-point "
    "band is 80&ndash;100% wide, narrowing to 41&ndash;94% at NNLO. The "
    "neutrino charm series behaves: NNLO is 3&ndash;3.5% below NLO at every "
    "energy and the band narrows from 5% to 3&ndash;4%. The low-Q&sup2; corner "
    "is what drives the neutral-current behaviour; the region keeps the "
    "common floor of 4 GeV&sup2; on both currents and shows it."
)


def _nnlo_over_nlo(resdir, e):
    v = nnlo_charm(resdir, e)
    r = sigma(resdir, REF, e)
    return None if not (v and r) else v / r[0]


def _charm_band_width(resdir, fn, e):
    with open(f"{BASE}/{resdir}/{_mhou_name(resdir, fn)}") as f:
        pt = json.load(f)["points"][f"{e:g}"]
    return pt["charm_rel_hi"] - pt["charm_rel_lo"]


def _same_scheme(e):
    """Muon NNLO/NLO charm with BOTH orders in FONLL.

    The figure divides by the massless reference, so its ratio folds a small
    scheme change into a large order change.  For a statement about the
    perturbative series that mixing is not acceptable, so the message quotes
    this number and the plotted one side by side.
    """
    return nnlo_charm("results", e) / sigma("results", FONLL_REF, e)[0]


def _charm_edge(fn, e):
    """Lower edge of the charm envelope, as a fraction."""
    with open(f"{BASE}/results/{_mhou_name('results', fn)}") as f:
        return json.load(f)["points"][f"{e:g}"]["charm_rel_lo"]


def _sherpa_ratio(resdir, key, e):
    v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
    return None if not (v and r) else v[0] / r[0]


def _incl_spread():
    """Largest max/min among the NLO matchings on the INCLUSIVE rate of the
    same region, the comparison the charm spread is quoted against."""
    rows = [("results", ["powheg", "herwig_nlo_powheg_full", "sherpa"]),
            ("results_nu", ["powheg_nu", "herwig_nlo_full", "sherpa_nlo"])]
    out = 0.0
    for d, ks in rows:
        for e in _ES:
            v = [sigma(d, k, e)[0] for k in ks]
            out = max(out, max(v) / min(v) - 1.0)
    return out


CLAIMS = [
    {"what": "the inclusive rate of the same region spreads by under 3.5% "
             "across the matchings, under a fifth of the charm spread",
     "check": lambda: _incl_spread() < 0.035 and _nlo_spread() > 5 * _incl_spread(),
     "detail": lambda: "%.2f%%" % (100 * _incl_spread())},
    {"what": "the charm region is q4w5 on tungsten, Q2 > 4 on BOTH currents, "
             "references with target-mass corrections",
     "check": lambda: REGION == "q4w5" and TSUF == "_W" and TMC == "_tmc"
     and Q2MIN["results"] == Q2MIN["results_nu"] == 4.0,
     "detail": lambda: f"{REGION}{TSUF}{TMC}, Q2MIN {Q2MIN}"},
    {"what": "the NLO matchings span up to 17% on muon charm and agree to 3% "
             "on neutrino charm",
     "check": lambda: 0.15 < _nlo_spread("results") < 0.20
     and _nlo_spread("results_nu") < 0.035,
     "detail": lambda: "muon %.1f%%, neutrino %.1f%%" % (
         100 * _nlo_spread("results"), 100 * _nlo_spread("results_nu"))},
    {"what": "Herwig at 0.97-1.03 of the reference on the muon side and 1.00 "
             "on the neutrino side (it read 0.74-0.98 when its remnant "
             "failures were replaced by new events); POWHEG-RES 0.90-0.92, "
             "POWHEG and Sherpa all within 12% of the reference and of "
             "themselves",
     "check": lambda: all(0.96 < _r("results", "herwig_nlo_powheg_full_charmfinal", e) < 1.04
                          for e in _ES)
     and all(0.995 < _r("results_nu", "herwig_nlo_full_charmfinal", e) < 1.005
             for e in _ES)
     and all(0.0 < _r(d, k, e)
                          for d, k in (("results", "herwig_nlo_powheg_full_charmfinal"),
                                       ("results_nu", "herwig_nlo_full_charmfinal"))
                          for e in _ES)
     and all(0.895 < _r("results", "powheg_charmfinal", e) < 0.925 for e in _ES)
     and all(abs(_r(d, k, e) - 1.0) < 0.12 for d, ks in _NLO for k in ks
             if "herwig" not in k for e in _ES)
     and all(max(_r(d, k, e) for e in _ES) - min(_r(d, k, e) for e in _ES) < 0.12
             for d, ks in _NLO for k in ks if "herwig" not in k),
     "detail": lambda: "; ".join("%s %s" % (k.split("_charm")[0], " ".join(
         "%.3f" % _r(d, k, e) for e in _ES)) for d, ks in _NLO for k in ks)},
    {"what": "the charm-mass effect ZM/FONLL is 3-7% on the muon side and "
             "1-4% on the neutrino side",
     "check": lambda: all(1.02 < sigma("results", REF, e)[0] / sigma("results", FONLL_REF, e)[0] < 1.08
                          and 1.005 < sigma("results_nu", REF, e)[0] / sigma("results_nu", FONLL_REF, e)[0] < 1.045
                          for e in _ES),
     "detail": lambda: "mu " + ", ".join("%.3f" % (sigma("results", REF, e)[0] / sigma("results", FONLL_REF, e)[0]) for e in _ES)
     + "; nu " + ", ".join("%.3f" % (sigma("results_nu", REF, e)[0] / sigma("results_nu", FONLL_REF, e)[0]) for e in _ES)},
    {"what": "POWHEG-V2mc, the only entry with the charm mass, sits 5-6% "
             "below the FONLL reference and within 1% of YADISM FFNS n_f = 3 "
             "with the charm PDF set to zero, its own process (2026-09-19)",
     "check": lambda: all(0.935 < sigma("results_nu", "powheg_nu_mc_charmfinal", e)[0]
                          / sigma("results_nu", FONLL_REF, e)[0] < 0.955
                          and abs(sigma("results_nu", "powheg_nu_mc_charmfinal", e)[0]
                                  / sigma("results_nu", "yadism_charm_nlo_ffns3", e)[0] - 1) < 0.01
                          for e in _ES),
     "detail": lambda: ", ".join("%.4f" % (sigma("results_nu", "powheg_nu_mc_charmfinal", e)[0]
                                           / sigma("results_nu", FONLL_REF, e)[0]) for e in _ES)},
    {"what": "in FONLL at both orders the muon NNLO/NLO charm ratio is "
             "0.51 at 400 GeV rising to 0.74 at 4 TeV",
     "check": lambda: 0.48 < _same_scheme(400.0) < 0.54 and 0.71 < _same_scheme(4000.0) < 0.77
     and all(_same_scheme(a) < _same_scheme(b) for a, b in zip(_ES, _ES[1:])),
     "detail": lambda: ", ".join("%.3f" % _same_scheme(e) for e in _ES)},
    {"what": "the muon NLO charm band is 80-100% wide and narrows to 41-94% "
             "at NNLO at every energy",
     "check": lambda: all(0.78 < _charm_band_width("results", "mhou_sigma_fonll_mu_zm.json", e) < 1.02
                          and 0.39 < _charm_band_width("results", "mhou_sigma_fonll_mu_pto2.json", e) < 0.96
                          and _charm_band_width("results", "mhou_sigma_fonll_mu_pto2.json", e)
                          < _charm_band_width("results", "mhou_sigma_fonll_mu_zm.json", e) for e in _ES),
     "detail": lambda: ", ".join("%.0f%%->%.0f%%" % (
         100 * _charm_band_width("results", "mhou_sigma_fonll_mu_zm.json", e),
         100 * _charm_band_width("results", "mhou_sigma_fonll_mu_pto2.json", e)) for e in _ES)},
    {"what": "the neutrino NNLO charm correction is -3 to -3.5% and its band "
             "narrows from 5% to 3-4%",
     "check": lambda: all(0.96 < _nnlo_over_nlo("results_nu", e) < 0.975
                          and 0.045 < _charm_band_width("results_nu", "mhou_sigma_fonll_zm.json", e) < 0.055
                          and 0.028 < _charm_band_width("results_nu", "mhou_sigma_fonll_pto2_zm.json", e) < 0.042
                          for e in _ES),
     "detail": lambda: ", ".join("%.4f (%.1f%%->%.1f%%)" % (
         _nnlo_over_nlo("results_nu", e),
         100 * _charm_band_width("results_nu", "mhou_sigma_fonll_zm.json", e),
         100 * _charm_band_width("results_nu", "mhou_sigma_fonll_pto2_zm.json", e)) for e in _ES)},
    {"what": "the NNLO curve and its band come from ONE file, and where a "
             "standalone NNLO charm result also exists it agrees exactly",
     "check": lambda: all(
         abs(nnlo_charm("results", e)
             / sigma("results", "yadism_charm_nnlo_fonll_damp", e)[0] - 1.0) < 1e-3
         for e in _ES if sigma("results", "yadism_charm_nnlo_fonll_damp", e)),
     "detail": lambda: "muon NNLO charm: MHOU file == histos file"},
]


MHOU = {"results": "mhou_sigma_fonll_mu_zm.json",
        "results_nu": "mhou_sigma_fonll_zm.json"}


def _symmetrise(resdir, hi, lo):
    """Symmetrise the MUON band about the central value (user, 2026-08-29).

    On the muon neutral current the band is badly lopsided once the
    xiF = 1/2 points are removed -- +13.0/-1.4% at 4 TeV -- because the
    direction that was removed is the downward one.  Quoting that asymmetry
    would read as a physical statement about the sign of the missing higher
    orders, which it is not: it is an artefact of which points survived a
    reliability cut.  The larger side is therefore taken on both.

    The NEUTRINO band is left alone: it is nearly symmetric already, nothing
    was removed from it in practice, and symmetrising it would discard real
    information.
    """
    if resdir != "results":
        return hi, lo
    d = max(abs(hi), abs(lo))
    return d, -d


def _drop_half(resdir, d):
    """Must the xiF = 1/2 points be excluded?  Same rule as paper plot 1.

    They are dropped ONLY when they leave the parton distribution's grid,
    which is a statement about the Q2 floor: muF^2 = Q2min/4 is 1.0 at
    Q2min = 4, below NNPDF4.0's Q0^2 = 2.7225, and 2.75 at Q2min = 11, which
    is on it.  So the exclusion lifts itself on the muon column here rather
    than being carried over out of habit.
    """
    q0 = d.get("q2min_pdf")
    return q0 is None or Q2MIN[resdir] / 4.0 < q0 - 1e-9


def charm_band(resdir, e):
    """(rel_hi, rel_lo) of the 7-point envelope on sigma_CHARM.

    The xiF = 1/2 points are kept or dropped by _drop_half, exactly as on
    paper plot 1, so the two figures cannot end up quoting bands built
    differently for the same calculation.
    """
    p = f"{BASE}/{resdir}/{_mhou_name(resdir, MHOU[resdir])}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    pt = d["points"].get(f"{e:g}")
    if not pt or "charm_central_pb" not in pt:
        return None
    cen = pt["charm_central_pb"]
    drop = _drop_half(resdir, d)
    v = [r["sigma_charm_pb"] for r in pt["rows"]
         if not (drop and r["xiF"] == 0.5)]
    hi, lo = max(v) / cen - 1.0, min(v) / cen - 1.0
    return _symmetrise(resdir, hi, lo) if drop else (hi, lo)


# ALIGNED WITH PAPER PLOT 1 (user, 2026-08-30): the same information, in the
# same places, with POWHEG-V2mc as the one addition.  So this figure too
# carries the NNLO curve and its scale band on both panels.
#
# THE CURVE AND ITS BAND COME FROM THE SAME FILE.  The MHOU output stores
# `charm_central_pb` beside the envelope, so taking both from there makes it
# impossible for a band to end up drawn around a curve from another run.  The
# muon side also has a standalone NNLO charm result, and the two agree to the
# last digit -- checked, not assumed, in the claims below.  The NEUTRINO side
# has no standalone file at all, which is the other reason to read the MHOU
# one.
NNLO = {
    "results":    ("mhou_sigma_fonll_mu_pto2.json", "YADISM NNLO (FONLL)"),
    "results_nu": ("mhou_sigma_fonll_pto2_zm.json", "YADISM NNLO (ZM)"),
}
NNLO_COLOUR = "#6a3d9a"


def _nnlo_pt(resdir, e):
    with open(f"{BASE}/{resdir}/{_mhou_name(resdir, NNLO[resdir][0])}") as f:
        return json.load(f)["points"].get(f"{e:g}")


def nnlo_charm(resdir, e):
    """NNLO charm cross-section in pb, from the MHOU file's central row."""
    pt = _nnlo_pt(resdir, e)
    return None if not pt else pt["charm_central_pb"]


def nnlo_charm_band(resdir, e):
    """(rel_hi, rel_lo) of the FULL 7-point NNLO envelope on charm.

    Nothing excluded, as in paper plot 1 and for the same reason: the
    xiF = 1/2 instability that forces the exclusion at NLO is a feature of
    that order.
    """
    pt = _nnlo_pt(resdir, e)
    return None if not pt else (pt["charm_rel_hi"], pt["charm_rel_lo"])


def panel(ax, axr, resdir, rows, title, xlab, div, unit):
    es = list(beams.BENCH_ENERGIES)
    pts = [(e, sigma(resdir, REF, e)) for e in es]
    pts = [(e, v[0] / div) for e, v in pts if v]
    if len(pts) >= 2:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ":",
                color="#111111", lw=1.6, label=tex(REF_LABEL), zorder=1)
    # FONLL against the massless reference: the gap to unity is the
    # charm-mass effect, and on CHARM it is an order of magnitude larger than
    # on the inclusive rate.
    # FONLL as an absolute curve in the TOP panel -- it is a cross-section
    # like everything else there -- and its ratio in the lower panel with no
    # legend entry of its own, since the colour already identifies it.
    fabs = [(e, sigma(resdir, FONLL_REF, e)) for e in es]
    fabs = [(e, v[0] / div) for e, v in fabs if v]
    if len(fabs) >= 2:
        ax.plot([q[0] for q in fabs], [q[1] for q in fabs], "-",
                color="#c0392b", lw=1.6, zorder=2,
                label=tex("YADISM (FONLL)"))
    fp = [(e, sigma(resdir, FONLL_REF, e), sigma(resdir, REF, e)) for e in es]
    fp = [(e, v[0] / r[0]) for e, v, r in fp if v and r]
    if len(fp) >= 2:
        axr.plot([q[0] for q in fp], [q[1] for q in fp], "-",
                 color="#c0392b", lw=1.8, zorder=4)
    # THE FIXED-FLAVOUR CURVE, neutrino panel only (user, 2026-09-19: "yes"
    # to drawing it).  YADISM FFNS n_f = 3 with the charm PDF set to zero is
    # POWHEG-V2mc's own process (analysis/yadism_cc_charm_calc.py nlo ffns3),
    # solid, in both panels: V2mc sits on it to <1%, and its gap to FONLL (6%)
    # is FONLL's charm-initiated channels.
    # ITS OWN COLOUR, not V2mc's brown (user, 2026-09-21): drawn in the same
    # colour the analytic curve and the generator read as one object.
    if resdir == "results_nu":
        mcol = FFNS3_COLOUR
        fa = [(e, sigma(resdir, FFNS3_REF, e)) for e in es]
        fa = [(e, v[0] / div) for e, v in fa if v]
        if len(fa) >= 2:
            ax.plot([q[0] for q in fa], [q[1] for q in fa], "-",
                    color=mcol, lw=1.5, zorder=2, label=tex(FFNS3_LABEL))
        fr = [(e, sigma(resdir, FFNS3_REF, e), sigma(resdir, REF, e)) for e in es]
        fr = [(e, v[0] / r[0]) for e, v, r in fr if v and r]
        if len(fr) >= 2:
            axr.plot([q[0] for q in fr], [q[1] for q in fr], "-",
                     color=mcol, lw=1.7, zorder=4)
    # the 7-point scale band on CHARM, in the same scheme as the reference.
    # DRAWN IN BOTH PANELS (user, 2026-09-01): on charm the band is 15-25%
    # wide, wide enough to read off the absolute panel, and it is the only
    # thing there that says how much of the generator spread the reference
    # itself cannot resolve.
    bands = [(e, charm_band(resdir, e)) for e in es]
    bands = [(e, b) for e, b in bands if b]
    if len(bands) >= 2:
        axr.fill_between([b[0] for b in bands],
                         [1.0 + b[1][1] for b in bands],
                         [1.0 + b[1][0] for b in bands],
                         color="#111111", alpha=0.12, lw=0, zorder=0,
                         label=tex("NLO MHOU"))
        ref_abs = dict((e, sigma(resdir, REF, e)) for e in es)
        ab = [(e, b, ref_abs[e][0] / div) for e, b in bands if ref_abs.get(e)]
        if len(ab) >= 2:
            ax.fill_between([q[0] for q in ab],
                            [c * (1.0 + b[1]) for _e, b, c in ab],
                            [c * (1.0 + b[0]) for _e, b, c in ab],
                            color="#111111", alpha=0.12, lw=0, zorder=0,
                            label=tex("NLO MHOU"))

    # NNLO, with its own band, exactly as in paper plot 1
    nkey, nlabel = NNLO[resdir]
    nabs = [(e, nnlo_charm(resdir, e)) for e in es]
    nabs = [(e, v / div) for e, v in nabs if v]
    if len(nabs) >= 2:
        ax.plot([q[0] for q in nabs], [q[1] for q in nabs], "-.",
                color=NNLO_COLOUR, lw=1.8, zorder=2, label=tex(nlabel))
        nr = [(e, nnlo_charm(resdir, e), sigma(resdir, REF, e)) for e in es]
        nr = [(e, v / r[0]) for e, v, r in nr if v and r]
        axr.plot([q[0] for q in nr], [q[1] for q in nr], "-.",
                 color=NNLO_COLOUR, lw=1.9, zorder=4)
        nb = [(e, nnlo_charm_band(resdir, e), dict(nr).get(e)) for e in es]
        nb = [(e, b, c) for e, b, c in nb if b and c]
        if len(nb) >= 2:
            axr.fill_between([q[0] for q in nb],
                             [c * (1.0 + b[1]) for _e, b, c in nb],
                             [c * (1.0 + b[0]) for _e, b, c in nb],
                             color=NNLO_COLOUR, alpha=0.16, lw=0, zorder=0,
                             label=tex("NNLO MHOU"))

    # the massive-charm row.  Its ratio uses the SAME denominator as every
    # other curve in the panel (user, 2026-08-30) -- see NU_MC_ROW.
    if resdir == "results_nu":
        lab, key, colour, marker, ls = NU_MC_ROW
        mp, mr = [], []
        for e in es:
            v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
            if not v:
                continue
            mp.append((e, v[0] / div, v[1] / div))
            if r:
                mr.append((e, v[0] / r[0], v[1] / r[0]))
        if mp:
            ax.errorbar([p[0] for p in mp], [p[1] for p in mp],
                        yerr=[p[2] for p in mp], marker=marker, ms=6.5,
                        color=colour, lw=1.7, ls=ls, capsize=2.5,
                        label=tex(lab), zorder=3)
            if mr:
                axr.errorbar([p[0] for p in mr], [p[1] for p in mr],
                             yerr=[p[2] for p in mr], marker=marker, ms=6.5,
                             color=colour, lw=1.7, ls=ls, capsize=2.5,
                             zorder=3)
    for lab, key, colour, marker, ls in rows:
        pts, rat = [], []
        for e in es:
            v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
            if not v:
                continue
            pts.append((e, v[0] / div, v[1] / div))
            if r:
                rat.append((e, v[0] / r[0], v[1] / r[0]))
        if not pts:
            continue
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts],
                    yerr=[p[2] for p in pts], marker=marker, ms=6.5,
                    color=colour, lw=1.7, ls=ls, capsize=2.5, label=tex(lab),
                    zorder=3)
        if rat:
            axr.errorbar([p[0] for p in rat], [p[1] for p in rat],
                         yerr=[p[2] for p in rat], marker=marker, ms=6.5,
                         color=colour, lw=1.7, ls=ls, capsize=2.5, zorder=3)
    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-")
    for a in (ax, axr):
        a.set_xscale("log")
        a.grid(alpha=0.22, lw=0.6)
        a.set_xlim(min(es) * 0.82, max(es) * 1.22)
    # LINEAR y (user).  The scan spans about an order of magnitude, which a
    # linear axis holds without compressing the low-energy end into the
    # frame edge, and it makes the per-cent agreement the ratio panel
    # quantifies visible as a single band rather than as coincident curves.
    ax.set_ylabel(tex(rf"$\sigma_{{\rm charm}}$  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
    ax.tick_params(labelbottom=False)
    axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])
    axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    # ONE LEGEND, in the top panel, carrying the ratio panel's entries too.
    # A legend inside the ratio panel sat on top of the very curves it was
    # naming (user, 2026-08-29); the top panel has empty space and the ratio
    # panel does not.
    # top legend: absolute cross-sections only.  The ratio panel carries the
    # one entry that exists only as a ratio -- the band -- and nothing else,
    # so neither legend mixes two kinds of quantity.
    _h, _l = ax.get_legend_handles_labels()
    # EACH BAND SHARES ITS CURVE'S ENTRY (user, 2026-10-04): NLO MHOU with
    # the reference, NNLO MHOU (ratio panel only) with the NNLO curve, so the
    # ratio-panel legend below is left with nothing.
    for hh, ss in zip(*axr.get_legend_handles_labels()):
        if ss not in _l:
            _h.append(hh)
            _l.append(ss)
    _h, _l = plotstyle.merge_bands(_h, _l, {
        tex("NLO MHOU"): [tex(REF_LABEL)],
        tex("NNLO MHOU"): [tex(NNLO[resdir][1])]})
    leg = ax.legend(_h, _l, fontsize=plotstyle.FS_LEGEND - 2, frameon=True,
                    loc="upper left", handlelength=1.9, labelspacing=0.30)
    # BOTH LEGENDS IN THE TOP PANEL, as in paper plot 1: the ratio panel is
    # now filled by two bands and their curves, and a legend there lands on
    # the band it is naming.
    _rh, _rl = [], []
    # The NLO band is drawn in BOTH panels since 2026-09-01, so its entry is
    # already in the top legend; keep the lower legend to what exists only as
    # a ratio, rather than naming the same band twice.
    if _rh:
        _keep = [(h, s) for h, s in zip(_rh, _rl) if s not in _l]
        _rh = [h for h, _s in _keep]
        _rl = [s for _h, s in _keep]
    # NLO BEFORE NNLO in the legend (user, 2026-08-30), whatever order the
    # bands happened to be drawn in -- the reading order should be the order
    # of the perturbative expansion, not an artefact of z-order.
    if _rh:
        _pairs = sorted(zip(_rl, _rh), key=lambda t: ("NNLO" in t[0], t[0]))
        _rl, _rh = [t[0] for t in _pairs], [t[1] for t in _pairs]
    if _rh:
        ax.add_artist(leg)
        ax.legend(_rh, _rl, fontsize=plotstyle.FS_LEGEND - 3, frameon=True,
                  loc="lower right", handlelength=1.6, labelspacing=0.3)


def main():
    fig, axes = plt.subplots(
        2, 2, figsize=(11.0, 6.6), sharex="col",
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.0,
                     "wspace": 0.24})
    panel(axes[0][0], axes[1][0], "results", MU_ROWS,
          "Muon DIS",
          r"$E_{\mu}$  [GeV]", 1000.0, "nb")
    panel(axes[0][1], axes[1][1], "results_nu", NU_ROWS,
          "Neutrino DIS",
          r"$E_{\nu}$  [GeV]", 1.0, "pb")
    for a in (axes[1][0], axes[1][1]):
        a.set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                     fontsize=plotstyle.FS_YLABEL - 3)
    # THE TWO RATIO PANELS NO LONGER SHARE A RANGE (user, 2026-08-30), for
    # the reason paper plot 1 gives: with the NNLO band on the figure the two
    # panels now hold quantities of very different size, and a shared range
    # would either clip the muon one or flatten the neutrino one.  Each is
    # scaled to its own content and labelled separately.
    for a in (axes[1][0], axes[1][1]):
        lo, hi = a.dataLim.intervaly
        half = max(abs(hi - 1.0), abs(1.0 - lo)) * 1.14
        a.set_ylim(1.0 - half, 1.0 + half)
    fig.suptitle(tex(TITLE), y=0.98, fontsize=plotstyle.FS_SUPTITLE)
    # NO PROSE INSIDE THE FIGURE (user rule, 2026-08-29):
    # the explanation lives in MESSAGE and in the caption.
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
