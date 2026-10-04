#!/usr/bin/env python3
"""Paper plot 5: GENIE's charm production against the NLO reference, in
the CHARM region of the final paper.

The redo of the earlier figure for the "Paper plots"
tab.  the earlier production's figure -- the GENIE configurations, POWHEG-V2 (massless)
and POWHEG-V2mc (massive charm), ratio to YADISM (FONLL) -- with the 
conventions:

  * TUNGSTEN PER NUCLEON, (74 p + 110 n)/184, every generator on p and n.
  * THE CHARM REGION: Q2 > 4 GeV2 and W > 5 GeV, no y cut (user, 2026-09-14:
    the W floor is raised for charm so that no heavy-quark-initiated event
    falls below the threshold the showered generators cannot hadronise;
    paper Sec. 3).  Selection q4w5, applied to generators and references.
  * THE REFERENCE carries target-mass corrections, as every reference does.
  * FIVE ENERGIES, 400-4000 GeV; G18_02a stops at 1 TeV, its declared validity.

>>> NEUTRINO ONLY, AND THE REASON IS A MEASURED ZERO (rule 1b). <<<  GENIE's
muon neutral-current DIS has no charm matrix element and AGKY makes no charm
in fragmentation, so its NC charm cross-section is exactly zero; the claims
below check that it still is, on tungsten, at every G18_02a energy.

Usage: analysis/paper_plots/pp05_genie_charm_vs_energy.py
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
import numpy as np                                          # noqa: E402
import beams                                                # noqa: E402
import selection                                            # noqa: E402

SLUG = "pp05_genie_charm_vs_energy"
V1_NUMBER = 5
IN_PAPER = True
# The user's wording (2026-09-18).  The region is NOT in the figure: it is
# in the caption, as on every other paper plot.
TITLE = "Inclusive charm production vs beam energy"
OUTPUT = "pp_genie_charm_vs_energy.png"
RESULTS = "results_nu"
CAPTION = ("Fiducial charm-production cross-section per nucleon on tungsten "
           "against beam energy in neutrino charged-current scattering, for "
           "Q&sup2; &gt; 4 GeV&sup2; and W &gt; 5 GeV, for the GENIE "
           "configurations, with POWHEG-V2 (massless charm) and POWHEG-V2mc "
           "(massive charm), and YADISM in the fixed-flavour scheme with "
           "three light flavours and no charm in the proton (FFNS n_f = 3), "
           "the scheme POWHEG-V2mc implements. Ratio to YADISM (FONLL) "
           "beneath, with the "
           "reference's own seven-point scale band; the reference includes "
           "target-mass corrections. The G18_02a curves stop at 1 TeV, the "
           "tune's declared validity. Neutral-current charm is "
           "not shown because GENIE's electromagnetic DIS has no charm "
           "channel at all.")

SURFACE = "#ffffff"
REGION = "q4w5"
TSUF = "_W"
TMC = "_tmc"
REF = "yadism_charm_nlo_fonll_damp"
REF_LABEL = "YADISM (FONLL)"
# >>> THE BAND FILE, AND WHY IT WAS THE WRONG ONE. <<<  This read
# `mhou_diff_fonll_charm.json` until 2026-09-08, a file that has never
# existed: `mhou_diff.py` writes the per-BIN band on a differential
# distribution, and the neutrino side of it has only the ZM-VFNS charm
# variants anyway.  What a cross-section-against-energy figure needs is the
# band on the INTEGRATED sigma_fid, one number per energy, which is what
# `mhou_sigma_fonll.py --current nu` writes -- and its points already carry
# `charm_rel_hi`/`charm_rel_lo` beside the inclusive pair, put there for
# paper plot 6.  The loader below returns None when the file is absent, so
# the figure DREW WITHOUT ITS BAND and said nothing: exactly the silent
# failure CONVENTIONS.md rule 2 is about, and the reason there is now a claim
# testing that the band was found.
MHOU = f"mhou_sigma_fonll_{REGION}{TSUF}{TMC}.json"
Q2MIN = 4.0

# >>> A TUNE IS NOT DRAWN OUTSIDE ITS VALIDITY RANGE (user, 2026-09-01). <<<
# Same rule and same reason as paper plot 3: G18_02a is Bodek-Yang, declared
# valid to 1000 GeV, and GENIE clamps the spline and buries the warning
# rather than refusing the request -- so a 4 TeV G18_02a point is produced
# and looks healthy.  HEDIS GHE19_00a declares 1e12 GeV and keeps its point.
# The two GENIE charm curves that stop here are the two that must.
EMAX = {"genie_lo_charmfinal": 1000.0, "genie_nnpdf_charmfinal": 1000.0}


def _in_range(key, e):
    lim = EMAX.get(key)
    return lim is None or e <= lim + 1e-9

ROWS = [
    ("GENIE (GRV98LO)",          "genie_lo_charmfinal",    "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", "genie_nnpdf_charmfinal", "#56b4e9", "D", "--"),
    ("GENIE (HEDIS)",            "genie_charmfinal",       "#e69f00", "^", "-."),
    ("POWHEG-V2",                "powheg_nu_charmfinal",   "#0072b2", "s", ":"),
    ("POWHEG-V2mc",              "powheg_nu_mc_charmfinal", "#8c6d1f", "v",
     (0, (6.5, 1.6))),
]


def _esuf(e):
    tag = beams.Beams("mu", e).tag
    return "" if e == beams.ANCHOR_ENERGY else f"_{tag}"


def _path(key, e, resdir=RESULTS):
    """histos_<key>_q4w5_W[_tmc][_tag].json: the references carry the TMC
    tag after the target, the generators none."""
    tmc = TMC if key.startswith("yadism") else ""
    return f"{BASE}/{resdir}/histos_{key}_{REGION}{TSUF}{tmc}{_esuf(e)}.json"


def sigma(key, e):
    """(sigma_fid, err) in pb, or None; the file's Q2 stamp is checked.

    None outside the tune's validity range, so the point is absent from the
    figure AND from every claim rather than drawn with a caveat.
    """
    if not _in_range(key, e):
        return None
    p = _path(key, e)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if (d.get("selection") or {}).get("name") != REGION:
        raise SystemExit(f"{os.path.relpath(p, BASE)} is not a {REGION} result")
    q2 = (d.get("selection") or {}).get("q2_min")
    if q2 is not None and abs(q2 - Q2MIN) > 1e-9:
        raise SystemExit(f"{os.path.relpath(p, BASE)} was computed at "
                         f"Q2 > {q2:g}, but this figure is Q2 > {Q2MIN:g}")
    if d.get("stat_insufficient"):
        return None
    return d.get("sigma_fid_pb"), d.get("sigma_fid_err_pb") or 0.0


def _r(key, e):
    v, r = sigma(key, e), sigma(REF, e)
    return None if not (v and r) else v[0] / r[0]


def charm_band(e):
    """(rel_hi, rel_lo) of the 7-point NLO envelope on the charm reference."""
    p = f"{BASE}/{RESULTS}/{MHOU}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    pt = (d.get("points") or {}).get(f"{e:g}")
    if not pt:
        return None
    return pt.get("charm_rel_hi"), pt.get("charm_rel_lo")


# ----------------------------------------------------------------- claims --
_ES = list(beams.BENCH_ENERGIES)
_GENIE = ["genie_lo_charmfinal", "genie_nnpdf_charmfinal", "genie_charmfinal"]


def _es(key):
    """The energies this row is DRAWN at -- validity range applied."""
    return [e for e in _ES if _in_range(key, e)]


def _nc_charm_is_zero():
    """GENIE's muon NC charm, which is a measured zero rather than a gap."""
    out = []
    for e in _ES:
        if e > 1000.0:
            continue          # G18_02a's validity: the EM tune is not run above
        p = _path("genie_charmfinal", e, "results")
        if os.path.exists(p):
            with open(p) as f:
                out.append(json.load(f).get("sigma_fid_pb"))
    return out


def _incl_pdf_ratio(e):
    """GENIE NNPDF4.0 / GRV98 on the INCLUSIVE rate, same region."""
    a, b = sigma("genie_nnpdf", e), sigma("genie_lo", e)
    return None if not (a and b) else a[0] / b[0]


def _zm_over_fonll(e):
    """The analytic charm-mass effect: YADISM ZM-VFNS / FONLL on charm."""
    a, b = sigma("yadism_charm_nlo", e), sigma(REF, e)
    return None if not (a and b) else a[0] / b[0]


def _zm_over_ffns3(e):
    """The like-for-like analytic mass effect (2026-09-19): YADISM ZM-VFNS over
    FFNS n_f = 3 with the charm PDF set to zero -- POWHEG-V2mc's own process
    (analysis/yadism_cc_charm_calc.py nlo ffns3).  FONLL is NOT the massive
    twin of V2mc: with NNPDF4.0's fitted charm it also carries charm-initiated
    channels, 6% of the charm rate here, which FFNS n_f = 3 excludes."""
    a, b = sigma("yadism_charm_nlo", e), sigma("yadism_charm_nlo_ffns3", e)
    return None if not (a and b) else a[0] / b[0]


def _v2mc_over_ffns3(e):
    a, b = sigma("powheg_nu_mc_charmfinal", e), sigma("yadism_charm_nlo_ffns3", e)
    return None if not (a and b) else a[0] / b[0]


def _v2_over_v2mc(e):
    a, b = sigma("powheg_nu_charmfinal", e), sigma("powheg_nu_mc_charmfinal", e)
    return None if not (a and b) else a[0] / b[0]


MESSAGE = (
    "<b>GENIE's charm deficit is its parton distribution.</b> On tungsten, "
    "in the charm region (Q&sup2; &gt; 4 GeV&sup2;, W &gt; 5 GeV), the default "
    "tune reads 0.63 of the FONLL reference at 400 GeV, 700 GeV "
    "and 1 TeV &mdash; more than a third low on the quantity a dimuon measurement is "
    "made of. Substituting NNPDF4.0 for GRV98 LO inside the same tune "
    "multiplies the charm cross-section by <b>1.58&ndash;1.60</b> and takes it "
    "onto the reference (0.998&ndash;1.003), within 5% of POWHEG-V2, while the "
    "inclusive rate moves by only 4% under the same substitution: "
    "charged-current charm is CKM-weighted s &rarr; c off the strange sea, "
    "and GRV98 LO's strange density is about half NNPDF4.0's in the x range "
    "this region probes. A generator tuned to the total rate can still be "
    "wrong by sixty per cent on its charm."
    "\n\n"
    "<b>GENIE HEDIS sits between the two</b>, at 0.85&ndash;0.89 from 400 GeV "
    "to 4 TeV; it is a different calculation throughout and not the tune "
    "FASER runs. The two G18_02a rows stop at 1 TeV, the validity range that "
    "tune declares. The G18_02a rows are flat in energy (they rose by 4% with the "
    "uncorrected charm propagator, see genie/README.md) and HEDIS rises "
    "towards the reference; the POWHEG-V2 and POWHEG-V2mc ratios are flat to within two per cent over the whole "
    "range."
    "\n\n"
    "<b>The generator reproduces the analytic charm-mass effect.</b> "
    "POWHEG-V2, massless as run, lies 3&ndash;5% above FONLL and within 3% "
    "of the massless YADISM curve; POWHEG-V2mc, with massive charm, lies "
    "5&ndash;6% below FONLL but within 1% of YADISM in its own scheme, "
    "fixed-flavour with three light flavours and no charm in the proton. "
    "FONLL is not that scheme: with NNPDF4.0's fitted charm it also counts "
    "charm-initiated scattering, 6% of the charm rate here. Compared like "
    "for like, the massless/massive ratio is 1.09&ndash;1.11 between "
    "POWHEG-V2 and POWHEG-V2mc and 1.07&ndash;1.11 between ZM-VFNS and the fixed-flavour "
    "YADISM calculation, the same to 0.6% up to 1 TeV."
    "\n\n"
    "<b>There is no muon panel because there is nothing to draw</b>: GENIE's "
    "electromagnetic DIS has no charm matrix element and its AGKY "
    "hadronisation produces no charm in fragmentation, so its "
    "neutral-current charm cross-section is a measured zero at every energy "
    "the tune is run at, rather than a small number or a missing result."
)

CLAIMS = [
    {"what": "the reference is the MASSIVE one (FONLL, with target-mass "
             "corrections), in the charm region q4w5 on tungsten",
     "check": lambda: REF.endswith("fonll_damp") and REGION == "q4w5"
     and TSUF == "_W" and TMC == "_tmc",
     "detail": lambda: f"{REF}_{REGION}{TSUF}{TMC}"},
    {"what": "GENIE's default tune is a third below the charm reference "
             "(0.63) at every energy it is run at",
     "check": lambda: all(0.60 < _r("genie_lo_charmfinal", e) < 0.75
                          for e in _es("genie_lo_charmfinal")),
     "detail": lambda: ", ".join("%.4f" % _r("genie_lo_charmfinal", e)
                                 for e in _es("genie_lo_charmfinal"))},
    {"what": "NNPDF4.0 in the same tune multiplies the charm cross-section "
             "by 1.58-1.60, to 0.998-1.003 of the reference, while the "
             "inclusive rate moves by 4%",
     "check": lambda: all(1.5 < (_r("genie_nnpdf_charmfinal", e)
                                 / _r("genie_lo_charmfinal", e)) < 1.65
                          and 0.99 < _r("genie_nnpdf_charmfinal", e) < 1.01
                          and 1.035 < _incl_pdf_ratio(e) < 1.05
                          for e in _es("genie_lo_charmfinal")),
     "detail": lambda: ", ".join(
         "charm x%.3f (%.4f), incl x%.3f" % (
             _r("genie_nnpdf_charmfinal", e) / _r("genie_lo_charmfinal", e),
             _r("genie_nnpdf_charmfinal", e), _incl_pdf_ratio(e))
         for e in _es("genie_lo_charmfinal"))},
    {"what": "with NNPDF4.0 GENIE's charm is within 5% of POWHEG-V2",
     "check": lambda: all(abs(_r("genie_nnpdf_charmfinal", e)
                              / _r("powheg_nu_charmfinal", e) - 1.0) < 0.05
                          for e in _es("genie_nnpdf_charmfinal")),
     "detail": lambda: ", ".join(
         "%.4f" % (_r("genie_nnpdf_charmfinal", e) / _r("powheg_nu_charmfinal", e))
         for e in _es("genie_nnpdf_charmfinal"))},
    {"what": "GENIE HEDIS does better on charm than the default tune, at "
             "0.85-0.89",
     "check": lambda: all(0.83 < _r("genie_charmfinal", e) < 0.91 for e in _ES)
     and all(abs(_r("genie_charmfinal", e) - 1.0)
             < abs(_r("genie_lo_charmfinal", e) - 1.0)
             for e in _es("genie_lo_charmfinal")),
     "detail": lambda: "HEDIS " + ", ".join(
         "%.4f" % _r("genie_charmfinal", e) for e in _ES)},
    {"what": "the G18_02a rows are flat in energy to 1% (they rose by 4% before the "
             "2026-10-01 charm-propagator fix) and HEDIS rises towards the "
             "reference, while both POWHEG ratios are flat to 2% across 400-4000 GeV",
     "check": lambda: all(abs(_r(k, max(_es(k))) / _r(k, 400.0) - 1) < 0.01
                          for k in _GENIE if k != "genie_charmfinal")
     and _r("genie_charmfinal", max(_es("genie_charmfinal"))) > _r("genie_charmfinal", 400.0)
     and all(max(_r(k, e) for e in _ES) - min(_r(k, e) for e in _ES) < 0.02
             for k in ("powheg_nu_charmfinal", "powheg_nu_mc_charmfinal")),
     "detail": lambda: ", ".join(
         "%s %.3f->%.3f" % (k.split("_charm")[0], _r(k, 400.0), _r(k, max(_es(k))))
         for k in _GENIE + ["powheg_nu_charmfinal", "powheg_nu_mc_charmfinal"])},
    {"what": "POWHEG-V2 (massless) is 3-5% above FONLL and within 3% of "
             "ZM-VFNS; POWHEG-V2mc (massive) is 5-6% below FONLL",
     "check": lambda: all(1.025 < _r("powheg_nu_charmfinal", e) < 1.05
                          and abs(_r("powheg_nu_charmfinal", e) / _zm_over_fonll(e) - 1) < 0.0305
                          and 0.935 < _r("powheg_nu_mc_charmfinal", e) < 0.955
                          for e in _ES),
     "detail": lambda: ", ".join(
         "V2 %.4f (V2/ZM %.4f), V2mc %.4f" % (
             _r("powheg_nu_charmfinal", e),
             _r("powheg_nu_charmfinal", e) / _zm_over_fonll(e),
             _r("powheg_nu_mc_charmfinal", e)) for e in _ES)},
    {"what": "POWHEG-V2mc is within 1% of YADISM in its own scheme (FFNS "
             "n_f = 3, no charm PDF) at every energy; like for like the "
             "massless/massive ratio is 1.09-1.11 in POWHEG and 1.07-1.11 in "
             "YADISM (ZM-VFNS/FFNS), equal to 0.6% up to 1 TeV -- not the "
             "1.01-1.04 of ZM-VFNS/FONLL, which is a different comparison",
     "check": lambda: all(abs(_v2mc_over_ffns3(e) - 1) < 0.01
                          and 1.085 < _v2_over_v2mc(e) < 1.115
                          and 1.065 < _zm_over_ffns3(e) < 1.12
                          and 1.005 < _zm_over_fonll(e) < 1.045 for e in _ES)
     and all(abs(_v2_over_v2mc(e) / _zm_over_ffns3(e) - 1) < 0.006
             for e in _ES if e <= 1000.0),
     "detail": lambda: ", ".join("%g GeV POWHEG %.3f, ZM/FFNS3 %.3f, ZM/FONLL %.3f, V2mc/FFNS3 %.3f" % (
         e, _v2_over_v2mc(e), _zm_over_ffns3(e), _zm_over_fonll(e), _v2mc_over_ffns3(e)) for e in _ES)},
    {"what": "the two G18_02a charm rows stop at 1 TeV, the tune's declared "
             "validity range, while HEDIS keeps its 4 TeV point",
     "check": lambda: (sigma("genie_lo_charmfinal", 4000.0) is None
                       and sigma("genie_nnpdf_charmfinal", 4000.0) is None
                       and sigma("genie_charmfinal", 4000.0) is not None),
     "detail": lambda: "G18_02a to %g GeV, HEDIS to %g"
                       % (max(_es("genie_lo_charmfinal")),
                          max(_es("genie_charmfinal")))},
    {"what": "there is no muon panel because GENIE's NC charm is a MEASURED "
             "zero at every G18_02a energy, not a missing result",
     "check": lambda: (len(_nc_charm_is_zero()) == 3
                       and all(v == 0.0 for v in _nc_charm_is_zero())),
     "detail": lambda: "muon NC charm = " + ", ".join(
         "%g" % v for v in _nc_charm_is_zero())},
    {"what": "the NLO seven-point band on the FONLL charm reference is "
             "actually present at every energy",
     "check": lambda: all(charm_band(e) is not None
                          and 0.005 < charm_band(e)[0] < 0.08
                          and -0.08 < charm_band(e)[1] < -0.005
                          for e in _ES),
     "detail": lambda: ", ".join(
         "%g GeV +%.2f/%.2f%%" % (e, 100 * charm_band(e)[0],
                                  100 * charm_band(e)[1]) for e in _ES)},
]


def _fig53_ylim():
    """The linear y range of Fig. 5.3's neutrino panel, taken from pp06 itself.

    Same axis as Fig. 5.3 so the two compare like with like (user,
    2026-09-30).  The panel is built off-screen by pp06's own code rather
    than copied as numbers, so a change to pp06 carries over here.
    """
    import pp06_charm_vs_energy as pp06
    f, axes = plt.subplots(2, 1, gridspec_kw={"height_ratios": [2.2, 1.0]})
    pp06.panel(axes[0], axes[1], "results_nu", pp06.NU_ROWS,
               "Neutrino DIS", r"$E_{\nu}$  [GeV]", 1.0, "pb")
    lim = axes[0].get_ylim()
    plt.close(f)
    return lim


def _check_frame(fig, ax, leg):
    """Refuse to write if a point leaves the fixed range or hits the legend.

    The range is Fig. 5.3's and is not adjusted here, so the check is all
    that stands between a new sample and a clipped or hidden point.
    Tested on every drawn vertex and segment, in display space, after a draw.
    """
    lo, hi = ax.get_ylim()
    for ln in ax.get_lines():
        y = ln.get_ydata()
        if len(y) and (min(y) < lo or max(y) > hi):
            raise SystemExit(f"{ln.get_label()} leaves the Fig. 5.3 range "
                             f"{lo:.3g}-{hi:.3g} pb -- refusing to write")
    fig.canvas.draw()
    box = leg.get_window_extent().expanded(1.03, 1.08)
    for ln in ax.get_lines():
        xy = ax.transData.transform(ln.get_xydata())
        # the segments too, not only the vertices: matplotlib draws them
        # straight in display space, so sampling them there is exact
        t = np.linspace(0.0, 1.0, 60)[:, None]
        pts = [xy[k] + t * (xy[k + 1] - xy[k]) for k in range(len(xy) - 1)]
        pts = np.vstack([xy] + pts) if pts else xy
        if any(box.contains(x, y) for x, y in pts):
            raise SystemExit("legend overlaps the curves -- refusing to write")


def main():
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.6, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.0})
    es = list(beams.BENCH_ENERGIES)
    pts = [(e, sigma(REF, e)) for e in es]
    pts = [(e, v[0]) for e, v in pts if v]
    if len(pts) >= 2:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ":",
                color="#111111", lw=1.6, label=tex(REF_LABEL), zorder=1)
    # THE FIXED-FLAVOUR CURVE (user, 2026-09-19): YADISM FFNS n_f = 3 with the
    # charm PDF set to zero, POWHEG-V2mc's own process, in V2mc's colour and
    # solid; V2mc sits on it to <1% while FONLL lies 6% above it (FONLL also
    # counts charm-initiated scattering).  See _zm_over_ffns3.
    # ITS OWN COLOUR, NOT V2mc's (user, 2026-10-04: in the same brown the two
    # could not be told apart): grey, with the other YADISM curve's black.
    mcol = "#7f7f7f"
    fa = [(e, sigma("yadism_charm_nlo_ffns3", e), sigma(REF, e)) for e in es]
    fa = [(e, v[0], r[0]) for e, v, r in fa if v and r]
    if len(fa) >= 2:
        ax.plot([q[0] for q in fa], [q[1] for q in fa], "-", color=mcol,
                lw=1.5, zorder=2, label=tex(r"YADISM (FFNS $n_f\!=\!3$)"))
        axr.plot([q[0] for q in fa], [q[1] / q[2] for q in fa], "-",
                 color=mcol, lw=1.7, zorder=4)
    bands = [(e, charm_band(e)) for e in es]
    bands = [(e, b) for e, b in bands if b and None not in b]
    if len(bands) >= 2:
        axr.fill_between([b[0] for b in bands],
                         [1.0 + b[1][1] for b in bands],
                         [1.0 + b[1][0] for b in bands],
                         color="#111111", alpha=0.12, lw=0, zorder=0,
                         label=tex("NLO MHOU"))
        ref_abs = dict((e, sigma(REF, e)) for e in es)
        ab = [(e, b, ref_abs[e][0]) for e, b in bands if ref_abs.get(e)]
        ax.fill_between([q[0] for q in ab],
                        [c * (1.0 + b[1]) for _e, b, c in ab],
                        [c * (1.0 + b[0]) for _e, b, c in ab],
                        color="#111111", alpha=0.12, lw=0, zorder=0,
                        label=tex("NLO MHOU"))
    for lab, key, colour, marker, ls in ROWS:
        pp, rat = [], []
        for e in es:
            v, r = sigma(key, e), sigma(REF, e)
            if not v:
                continue
            pp.append((e, v[0], v[1]))
            if r:
                rat.append((e, v[0] / r[0], v[1] / r[0]))
        if not pp:
            continue
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp],
                    yerr=[p[2] for p in pp], marker=marker, ms=6.0,
                    color=colour, lw=1.7, ls=ls, capsize=2.5, label=tex(lab),
                    zorder=3)
        if rat:
            axr.errorbar([p[0] for p in rat], [p[1] for p in rat],
                         yerr=[p[2] for p in rat], marker=marker, ms=6.0,
                         color=colour, lw=1.7, ls=ls, capsize=2.5, zorder=3)
    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-")
    for a in (ax, axr):
        a.set_xscale("log")
        a.grid(alpha=0.22, lw=0.6)
        a.set_xlim(min(es) * 0.82, max(es) * 1.22)
    # LINEAR y, SAME RANGE AS FIG. 5.3 (user, 2026-09-30), for a like-with-
    # like comparison with that figure's neutrino panel.
    ax.set_ylim(*_fig53_ylim())
    ax.set_ylabel(tex(r"$\sigma_{\rm charm}$  [pb]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex("Neutrino DIS"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.tick_params(labelbottom=False)
    axr.set_xlabel(tex(r"$E_{\nu}$  [GeV]"), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])
    axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    axr.set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                   fontsize=plotstyle.FS_YLABEL - 1)
    axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune="upper"))
    _h, _l = ax.get_legend_handles_labels()
    _seen, _hh, _ll = set(), [], []
    for h, s in zip(_h, _l):
        if s not in _seen:
            _seen.add(s)
            _hh.append(h)
            _ll.append(s)
    # ORDER (user, 2026-10-04): YADISM first, then POWHEG, then GENIE.
    _rank = lambda s: (0 if "YADISM" in s or "MHOU" in s else
                       1 if "POWHEG" in s else 2)
    _o = sorted(range(len(_ll)), key=lambda i: (_rank(_ll[i]), i))
    _hh, _ll = [_hh[i] for i in _o], [_ll[i] for i in _o]
    # the band shares the reference's entry (user, 2026-10-04)
    _hh, _ll = plotstyle.merge_bands(_hh, _ll,
                                     {tex("NLO MHOU"): [tex(REF_LABEL)]})
    # THE LEGEND SITS BELOW THE CURVES, IN THE CORNER THEY LEAVE EMPTY (user,
    # 2026-09-18: no overlap between curves and legend).  Every curve rises
    # with energy, so on the linear axis of Fig. 5.3 the upper left is the
    # free corner, as it is there; the range is fixed, so this is checked
    # after drawing rather than made room for.
    leg = ax.legend(_hh, _ll, fontsize=plotstyle.FS_LEGEND - 1, frameon=True,
                    loc="upper left", handlelength=1.9, labelspacing=0.32)
    _check_frame(fig, ax, leg)
    fig.suptitle(tex(TITLE), fontsize=plotstyle.FS_TITLE, y=0.98)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
