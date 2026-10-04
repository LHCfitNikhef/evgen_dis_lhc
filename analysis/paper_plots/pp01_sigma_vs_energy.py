#!/usr/bin/env python3
"""Paper plot 1: inclusive fiducial cross-section against beam energy,
in the FINAL benchmark region.

The redo of the earlier figure for the "Paper plots"
tab (user, 2026-09-13): "repeat 'Paper Plots' but with the very final
settings: W target (with NNPDF4.0 NNLO as PDF, so no nuclear effects),
Q2 > 4 GeV2 and W > 3 GeV at generation level.  Everything else unchanged."
Settings fixed with the user the same day (PRODUCTION.md):

  * TUNGSTEN PER NUCLEON, (74 sigma_p + 110 sigma_n)/184.  Every generator
    ran the proton and the neutron SEPARATELY, each on its own beam particle
    (analysis/combine_target.py combines them); the analytic references
    use the free-nucleon average set NNPDF40_nnlo_as_01180_W184free, which is
    the same thing exactly because structure functions are linear in the PDF.
  * Q2 > 4 GeV2 AND W > 3 GeV, NO y CUT, on BOTH currents -- the earlier muon
    floor of Q2 > 11 does not carry over.
  * 400, 700, 1000, 2000 and 4000 GeV (beams.BENCH_ENERGIES).
  * The YADISM references carry TARGET-MASS CORRECTIONS (Georgi-Politzer,
    yadism TMC mode 3), and the TMC-off results are kept on disk so the
    effect can be quoted: it is tiny here, which the claims check.

Everything else is the earlier production's: the three NLO matchings per current, the massless
(ZM) reference with FONLL in the ratio panel, the NNLO curve with its own
seven-point band, the NLO band with its xiF = 1/2 points dropped where they
leave the PDF grid (they do at Q2min/4 = 1 GeV2, on both currents now) and
symmetrised on the muon side.

GENIE IS DELIBERATELY ABSENT, exactly as in the earlier production and for the earlier production's reason: this
figure benchmarks whether the NLO matchings agree with each other and with
the reference, and GENIE's offset would set the axis range.  The paper-plots
carve-out of CONVENTIONS.md rule 1b applies (stated here, as it requires).

Usage: analysis/paper_plots/pp01_sigma_vs_energy.py
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
import target                                               # noqa: E402

SLUG = "pp01_sigma_vs_energy"
V1_NUMBER = 1
# IN THE PAPER, in place of the earlier production's figure (user, 2026-09-14).  Read by
# tools/check_paper_figures_used.py, which then requires THIS output in the
# document and releases the earlier script of the same SLUG from that requirement.
IN_PAPER = True
TITLE = "Fiducial NLO total cross-section vs beam energy"
# NOT "paper_*": tools/sync_paper_figures.sh copies every paper_*.png into the
# tracked paper/figures/ by glob; a figure enters the paper BY NAME
# (its EXTRA list), which is how the user put this one in on 2026-09-14.
OUTPUT = "pp_sigma_vs_energy.png"
RESULTS = "results_nu"
CAPTION = ("Fiducial cross-section per nucleon on tungsten against beam "
           "energy, for Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no y "
           "cut: muon neutral current (left) and neutrino charged current "
           "(right), for the three NLO matchings. Ratio to YADISM (ZM) with "
           "target-mass corrections beneath, with YADISM (FONLL) drawn there "
           "too. Each panel also carries the NNLO result with its own "
           "seven-point scale band &mdash; in FONLL on the muon side and in "
           "ZM on the neutrino side.")

SURFACE = "#ffffff"
REGION = "q4w3"
ES = list(beams.BENCH_ENERGIES)
Q2MIN = 4.0
W_MIN = 3.0

# The reference: massless, per nucleon on tungsten, TMC ON (user, 2026-09-13).
TMC = "_tmc"
REF = f"yadism_nlo_{REGION}_W{TMC}"
REF_NOTMC = f"yadism_nlo_{REGION}_W"
REF_LABEL = "YADISM (ZM)"
FONLL_REF = f"yadism_nlo_fonll_damp_{REGION}_W{TMC}"

# (label, result key, colour, marker, linestyle) -- the earlier production's rows, the earlier production's styling.
MU_ROWS = [
    ("POWHEG-RES", f"powheg_{REGION}_W",                 "#0072b2", "o", "-"),
    ("Herwig",     f"herwig_nlo_powheg_full_{REGION}_W", "#d55e00", "s", "--"),
    ("Sherpa",     f"sherpa_{REGION}_W",                 "#009e73", "^", "-."),
]
NU_ROWS = [
    ("POWHEG-V2",  f"powheg_nu_{REGION}_W",              "#0072b2", "o", "-"),
    ("Herwig",     f"herwig_nlo_full_{REGION}_W",        "#d55e00", "s", "--"),
    ("Sherpa",     f"sherpa_nlo_{REGION}_W",             "#009e73", "^", "-."),
]
ROWS = {"results": MU_ROWS, "results_nu": NU_ROWS}

# Bands: NLO in the reference's own scheme (ZM), NNLO in the scheme of the
# NNLO curve.  The region/target/TMC tag sits right after the current, as
# mhou_sigma_fonll.py writes it.
_TAG = f"_{REGION}_W{TMC}"
MHOU = {"results": f"mhou_sigma_fonll_mu{_TAG}_zm.json",
        "results_nu": f"mhou_sigma_fonll{_TAG}_zm.json"}
NNLO = {
    "results":    (f"yadism_nnlo_fonll_damp_{REGION}_W{TMC}",
                   "YADISM NNLO (FONLL)", f"mhou_sigma_fonll_mu{_TAG}_pto2.json"),
    "results_nu": (f"yadism_nnlo_{REGION}_W{TMC}", "YADISM NNLO (ZM)",
                   f"mhou_sigma_fonll{_TAG}_pto2_zm.json"),
}
NNLO_COLOUR = "#6a3d9a"


def _load(resdir, key, e):
    p = f"{BASE}/{resdir}/histos_{beams.at_energy(key, e)}.json"
    if not os.path.exists(p):
        return None, p
    with open(p) as f:
        return json.load(f), p


def sigma(resdir, key, e):
    """(sigma_fid, err) in pb per nucleon, or None if absent or unusable.

    EVERY FILE IS CHECKED to be what the figure claims before it is drawn:
    the q4w3 region with Q2 > 4 and W2 > 9 and no y window, and a TUNGSTEN
    per-nucleon result (target W).  A earlier proton file, a single-nucleon file
    or a TMC-off reference under the TMC-on name would each draw a plausible
    point and say nothing.
    """
    d, p = _load(resdir, key, e)
    if d is None:
        return None
    sel = d.get("selection") or {}
    rel = os.path.relpath(p, BASE)
    if (sel.get("name") != REGION or abs(sel.get("q2_min", -1) - Q2MIN) > 1e-9
            or abs(sel.get("w2_min", -1) - W_MIN ** 2) > 1e-9
            or sel.get("y_min") != 0.0 or sel.get("y_max") != 1.0):
        raise SystemExit(f"{rel}: not the {REGION} region ({sel})")
    if d.get("target") != "W":
        raise SystemExit(f"{rel}: target {d.get('target')!r}, not tungsten")
    if key.startswith("yadism"):
        want = target.TMC_CHOSEN if key.endswith(TMC) else 0
        if d.get("tmc") != want or d.get("pdf_set") != target.PDFSET_W:
            raise SystemExit(f"{rel}: tmc={d.get('tmc')} pdf="
                             f"{d.get('pdf_set')}, expected tmc={want} "
                             f"pdf={target.PDFSET_W}")
    if d.get("stat_insufficient"):
        return None
    return d.get("sigma_fid_pb"), d.get("sigma_fid_err_pb") or 0.0


def _r(resdir, key, e, ref=REF):
    v, r = sigma(resdir, key, e), sigma(resdir, ref, e)
    return None if not (v and r) else v[0] / r[0]


def _mhou(resdir, fn):
    p = f"{BASE}/{resdir}/{fn}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if (d.get("target") != "W" or d.get("tmc") != target.TMC_CHOSEN
            or (d.get("selection") or {}).get("name") != REGION):
        raise SystemExit(f"{resdir}/{fn}: not a tungsten TMC-on {REGION} band")
    return d


def _symmetrise(resdir, hi, lo):
    """Symmetrise the MUON band, as the earlier production did, when points were dropped."""
    if resdir != "results":
        return hi, lo
    dd = max(abs(hi), abs(lo))
    return dd, -dd


def _drop_half(d):
    """xiF = 1/2 leaves the PDF grid when Q2min/4 < Q0^2 -- true here (1 GeV2
    against 2.72), on both currents, exactly as in the earlier production before the muon floor
    was raised to 11."""
    q0 = d.get("q2min_pdf")
    return q0 is None or Q2MIN / 4.0 < q0 - 1e-9


def scale_band(resdir, e):
    """(rel_hi, rel_lo) of the NLO envelope, the earlier production's prescription."""
    d = _mhou(resdir, MHOU[resdir])
    if d is None or f"{e:g}" not in d["points"]:
        return None
    pt = d["points"][f"{e:g}"]
    drop = _drop_half(d)
    v = [r["sigma_fid_pb"] for r in pt["rows"]
         if not (drop and r["xiF"] == 0.5)]
    hi, lo = max(v) / pt["central_pb"] - 1.0, min(v) / pt["central_pb"] - 1.0
    return _symmetrise(resdir, hi, lo) if drop else (hi, lo)


def nnlo_band(resdir, e):
    """(rel_hi, rel_lo) of the FULL 7-point NNLO envelope, as before."""
    d = _mhou(resdir, NNLO[resdir][2])
    if d is None or f"{e:g}" not in d["points"]:
        return None
    pt = d["points"][f"{e:g}"]
    v = [r["sigma_fid_pb"] for r in pt["rows"]]
    return max(v) / pt["central_pb"] - 1.0, min(v) / pt["central_pb"] - 1.0


def _band_central_matches(resdir, fn, key, e):
    d = _mhou(resdir, fn)
    return abs(d["points"][f"{e:g}"]["central_pb"]
               / sigma(resdir, key, e)[0] - 1.0) < 1e-9


# ----------------------------------------------------------------- claims --
_CUR = [("results", [k for _l, k, *_ in MU_ROWS]),
        ("results_nu", [k for _l, k, *_ in NU_ROWS])]


def _dev(resdir, keys=None):
    keys = keys or dict(_CUR)[resdir]
    return max(abs(_r(resdir, k, e) - 1.0) for k in keys for e in ES)


def _spread(resdir):
    ks = dict(_CUR)[resdir]
    return max(max(_r(resdir, k, e) for k in ks)
               - min(_r(resdir, k, e) for k in ks) for e in ES)


def _tmc(resdir, e):
    return sigma(resdir, REF, e)[0] / sigma(resdir, REF_NOTMC, e)[0] - 1.0


def _nnlo_ratio(resdir, e):
    return _r(resdir, NNLO[resdir][0], e)


def _fonll(resdir, e):
    return _r(resdir, FONLL_REF, e)


def _pct(x):
    return f"{100 * x:+.2f}%"


def _in_nlo_band(resdir, e):
    lo, hi = scale_band(resdir, e)[1], scale_band(resdir, e)[0]
    return lo <= _nnlo_ratio(resdir, e) - 1.0 <= hi


MESSAGE = (
    "<b>The final region: tungsten per nucleon, Q&sup2; &gt; 4 GeV&sup2; and "
    "W &gt; 3 GeV with no y cut, on both currents.</b> Every generator ran "
    "the proton and the neutron separately, each on a genuine neutron or "
    "proton beam, combined as (74p + 110n)/184; the reference uses the "
    "free-nucleon average of NNPDF4.0 and includes target-mass corrections."
    "\n\n"
    "<b>On the neutrino charged current the three NLO matchings agree with "
    "the reference to half a per cent</b> &mdash; within 0.55% of YADISM "
    "(ZM) at every energy and within 0.45% of each other. <b>On the muon "
    "neutral current they agree to 2%:</b> Sherpa and Herwig within 0.35%, "
    "POWHEG-RES 0.8&ndash;1.9% below, the largest offset at 400 GeV."
    "\n\n"
    "<b>Target-mass corrections are negligible for the inclusive rate "
    "here</b>, at most 0.15% (muon, 400 GeV) and below 0.05% on the "
    "neutrino side: the W &gt; 3 GeV cut removes the region where they "
    "would matter. The charm-mass effect, the distance of the FONLL curve "
    "from unity, stays below 0.6%."
    "\n\n"
    "<b>The NNLO correction</b> is &minus;4.1% at 400 GeV growing to "
    "&minus;6.8% at 4 TeV on the muon side, with an NNLO band 1.7 to 2 "
    "times narrower than the NLO one; it sits just outside the NLO band at "
    "400 and 700 GeV (by under 0.2 percentage points) and inside it from "
    "1 TeV up. On the neutrino side the correction is &minus;1.2 to "
    "&minus;1.5%, about half a percentage point below the NLO band at every "
    "energy, whose seven-point width is only about &plusmn;1% and does not "
    "shrink at NNLO. The two NNLO curves are in different schemes &mdash; "
    "FONLL for the muon, ZM for the neutrino, where the massive "
    "charged-current coefficient functions, known at NNLO, are not "
    "implemented in the structure-function code used here."
)

CLAIMS = [
    {"what": "every file drawn is a tungsten per-nucleon result in the q4w3 "
             "region, and the references carry TMC and the W184free set",
     "check": lambda: all(sigma(d, k, e) is not None
                          for d, ks in _CUR for k in ks + [REF, FONLL_REF,
                                                           NNLO[d][0]]
                          for e in ES),
     "detail": lambda: f"{len(ES)} energies x 2 currents x 6 curves"},
    {"what": "each band's central value IS the reference curve it is drawn "
             "around (same calculation, to 1e-9)",
     "check": lambda: all(_band_central_matches(d, MHOU[d], REF, e)
                          and _band_central_matches(d, NNLO[d][2],
                                                    NNLO[d][0], e)
                          for d, _k in _CUR for e in ES),
     "detail": lambda: "NLO ZM and NNLO bands, both currents"},
    {"what": "neutrino CC: every NLO matching within 0.55% of the reference, "
             "and within 0.45% of each other",
     "check": lambda: _dev("results_nu") < 0.0055
     and _spread("results_nu") < 0.0045,
     "detail": lambda: "max |r-1| %.4f, max spread %.4f" % (
         _dev("results_nu"), _spread("results_nu"))},
    {"what": "muon NC: all three within 2% of the reference",
     "check": lambda: _dev("results") < 0.02,
     "detail": lambda: "max |r-1| %.4f" % _dev("results")},
    {"what": "muon NC: Sherpa and Herwig within 0.35% (Herwig was 0.5-0.8% "
             "above before its remnant failures were re-showered, "
             "2026-09-19), POWHEG-RES 0.8-1.9% below with its largest "
             "offset at 400 GeV",
     "check": lambda: (
         _dev("results", [MU_ROWS[2][1]]) < 0.0035
         and _dev("results", [MU_ROWS[1][1]]) < 0.0035
         and all(0.0075 < 1 - _r("results", MU_ROWS[0][1], e) < 0.019
                 for e in ES)
         and min(ES, key=lambda e: _r("results", MU_ROWS[0][1], e)) == 400.0),
     "detail": lambda: " | ".join(
         lab + " " + ", ".join("%.4f" % _r("results", k, e) for e in ES)
         for lab, k, *_ in MU_ROWS)},
    {"what": "TMC moves the inclusive reference by at most 0.15% (muon) and "
             "under 0.05% (neutrino)",
     "check": lambda: (max(abs(_tmc("results", e)) for e in ES) < 0.0015
                       and max(abs(_tmc("results_nu", e)) for e in ES)
                       < 0.0005),
     "detail": lambda: "mu " + ", ".join(_pct(_tmc("results", e)) for e in ES)
     + " | nu " + ", ".join(_pct(_tmc("results_nu", e)) for e in ES)},
    {"what": "the charm-mass effect (FONLL/ZM) stays below 0.6%",
     "check": lambda: all(abs(_fonll(d, e) - 1) < 0.006
                          for d, _k in _CUR for e in ES),
     "detail": lambda: ", ".join("%.4f" % _fonll(d, e)
                                 for d, _k in _CUR for e in ES)},
    {"what": "muon NNLO correction -4.1% at 400 GeV growing to -6.8% at "
             "4 TeV, monotonically",
     "check": lambda: (abs(_nnlo_ratio("results", 400.0) - 1 + 0.041) < 0.002
                       and abs(_nnlo_ratio("results", 4000.0) - 1 + 0.068)
                       < 0.002
                       and all(_nnlo_ratio("results", a)
                               > _nnlo_ratio("results", b)
                               for a, b in zip(ES, ES[1:]))),
     "detail": lambda: ", ".join(_pct(_nnlo_ratio("results", e) - 1)
                                 for e in ES)},
    {"what": "muon NNLO band 1.7-2.0x narrower than the NLO band",
     "check": lambda: all(
         1.65 < (scale_band("results", e)[0] - scale_band("results", e)[1])
         / (nnlo_band("results", e)[0] - nnlo_band("results", e)[1]) < 2.05
         for e in ES),
     "detail": lambda: ", ".join(
         "%.2fx" % ((scale_band("results", e)[0] - scale_band("results", e)[1])
                    / (nnlo_band("results", e)[0] - nnlo_band("results", e)[1]))
         for e in ES)},
    {"what": "muon NNLO just outside the NLO band at 400 and 700 GeV (by "
             "under 0.2 points), inside from 1 TeV up",
     "check": lambda: (
         all(not _in_nlo_band("results", e)
             and 0 < scale_band("results", e)[1]
             - (_nnlo_ratio("results", e) - 1) < 0.002 for e in (400.0, 700.0))
         and all(_in_nlo_band("results", e) for e in (1000.0, 2000.0, 4000.0))),
     "detail": lambda: ", ".join(
         "%g: %s vs low %s" % (e, _pct(_nnlo_ratio("results", e) - 1),
                              _pct(scale_band("results", e)[1])) for e in ES)},
    {"what": "neutrino NNLO correction -1.2 to -1.5%, 0.3-0.6 points below "
             "the NLO band at every energy",
     "check": lambda: all(
         -0.015 < _nnlo_ratio("results_nu", e) - 1 < -0.012
         and 0.003 < scale_band("results_nu", e)[1]
         - (_nnlo_ratio("results_nu", e) - 1) < 0.006 for e in ES),
     "detail": lambda: ", ".join(
         "%g: %s vs low %s" % (e, _pct(_nnlo_ratio("results_nu", e) - 1),
                              _pct(scale_band("results_nu", e)[1])) for e in ES)},
    {"what": "the neutrino band does not shrink from NLO to NNLO (width "
             "ratio 0.95-1.05) and is about +-1% wide",
     "check": lambda: all(
         0.95 < (scale_band("results_nu", e)[0] - scale_band("results_nu", e)[1])
         / (nnlo_band("results_nu", e)[0] - nnlo_band("results_nu", e)[1])
         < 1.05
         and scale_band("results_nu", e)[0] - scale_band("results_nu", e)[1]
         < 0.02 for e in ES),
     "detail": lambda: ", ".join(
         "%.2fx" % ((scale_band("results_nu", e)[0]
                     - scale_band("results_nu", e)[1])
                    / (nnlo_band("results_nu", e)[0]
                       - nnlo_band("results_nu", e)[1])) for e in ES)},
]


def panel(ax, axr, resdir, rows, title, xlab, div, unit):
    pts = [(e, sigma(resdir, REF, e)) for e in ES]
    pts = [(e, v[0] / div) for e, v in pts if v]
    if len(pts) >= 2:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ":",
                color="#111111", lw=1.6, label=tex(REF_LABEL), zorder=1)
    fabs = [(e, sigma(resdir, FONLL_REF, e)) for e in ES]
    fabs = [(e, v[0] / div) for e, v in fabs if v]
    if len(fabs) >= 2:
        ax.plot([q[0] for q in fabs], [q[1] for q in fabs], "-",
                color="#c0392b", lw=1.6, zorder=2,
                label=tex("YADISM (FONLL)"))
    fp = [(e, _r(resdir, FONLL_REF, e)) for e in ES]
    fp = [(e, v) for e, v in fp if v]
    if len(fp) >= 2:
        axr.plot([q[0] for q in fp], [q[1] for q in fp], "-",
                 color="#c0392b", lw=1.8, zorder=4)

    nnlo_key, nnlo_label, _ = NNLO[resdir]
    nabs = [(e, sigma(resdir, nnlo_key, e)) for e in ES]
    nabs = [(e, v[0] / div) for e, v in nabs if v]
    if len(nabs) >= 2:
        ax.plot([q[0] for q in nabs], [q[1] for q in nabs], "-.",
                color=NNLO_COLOUR, lw=1.8, zorder=2, label=tex(nnlo_label))
        np_ = [(e, _r(resdir, nnlo_key, e)) for e in ES]
        np_ = [(e, v) for e, v in np_ if v]
        axr.plot([q[0] for q in np_], [q[1] for q in np_], "-.",
                 color=NNLO_COLOUR, lw=1.9, zorder=4)
        nb = [(e, nnlo_band(resdir, e), _r(resdir, nnlo_key, e)) for e in ES]
        nb = [(e, b, c) for e, b, c in nb if b and c]
        if len(nb) >= 2:
            axr.fill_between([q[0] for q in nb],
                             [c * (1.0 + b[1]) for _e, b, c in nb],
                             [c * (1.0 + b[0]) for _e, b, c in nb],
                             color=NNLO_COLOUR, alpha=0.16, lw=0, zorder=0,
                             label=tex("NNLO MHOU"))

    for lab, key, colour, marker, ls in rows:
        pts, rat = [], []
        for e in ES:
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
    bands = [(e, scale_band(resdir, e)) for e in ES]
    bands = [(e, b) for e, b in bands if b]
    if len(bands) >= 2:
        axr.fill_between([b[0] for b in bands],
                         [1.0 + b[1][1] for b in bands],
                         [1.0 + b[1][0] for b in bands],
                         color="#111111", alpha=0.12, lw=0, zorder=0,
                         label=tex("NLO MHOU"))
        ref_abs = dict((e, sigma(resdir, REF, e)) for e in ES)
        ab = [(e, b, ref_abs[e][0] / div) for e, b in bands if ref_abs.get(e)]
        if len(ab) >= 2:
            ax.fill_between([q[0] for q in ab],
                            [c * (1.0 + b[1]) for _e, b, c in ab],
                            [c * (1.0 + b[0]) for _e, b, c in ab],
                            color="#111111", alpha=0.12, lw=0, zorder=0,
                            label=tex("NLO MHOU"))
    axr.axhline(1.0, color="#9aa1a9", lw=0.9, ls="-")
    for a in (ax, axr):
        a.set_xscale("log")
        a.grid(alpha=0.22, lw=0.6)
        a.set_xlim(min(ES) * 0.82, max(ES) * 1.22)
    ax.set_ylabel(tex(rf"$\sigma_{{\rm fid}}$ / nucleon  [{unit}]"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
    ax.tick_params(labelbottom=False)
    axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(ES)
    axr.set_xticklabels([f"{e:g}" for e in ES])
    axr.xaxis.set_minor_formatter(ticker.NullFormatter())
    _h, _l = ax.get_legend_handles_labels()
    _rh, _rl = axr.get_legend_handles_labels()
    # EACH BAND SHARES ITS CURVE'S ENTRY (user, 2026-10-04): NLO MHOU goes
    # with the YADISM reference, NNLO MHOU (drawn in the ratio panel only)
    # with the NNLO curve, so the ratio-panel legend has nothing left to say.
    for hh, ss in zip(_rh, _rl):
        if ss not in _l:
            _h.append(hh)
            _l.append(ss)
    _h, _l = plotstyle.merge_bands(_h, _l, {
        tex("NLO MHOU"): [tex(REF_LABEL)],
        tex("NNLO MHOU"): [tex(NNLO[resdir][1])]})
    # TWO COLUMNS AND HEADROOM (user, 2026-10-04: the one-column legend ran
    # into the muon panel's band and curves): the legend takes the top strip
    # and the y range is stretched so nothing drawn reaches it.
    leg = ax.legend(_h, _l, fontsize=plotstyle.FS_LEGEND - 3, frameon=True,
                    loc="upper left", handlelength=1.4, handletextpad=0.5,
                    labelspacing=0.3, ncol=2, columnspacing=0.7,
                    borderaxespad=0.3)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, lo + (hi - lo) * 1.30)
    _rh, _rl = [], []
    # Entries only the ratio panel draws (the NNLO band) are legended IN the
    # ratio panel, upper right (user, 2026-09-14); main() leaves headroom.
    if _rh:
        axr.legend(_rh, _rl, fontsize=plotstyle.FS_LEGEND - 3, frameon=True,
                   loc="upper right", handlelength=1.6, labelspacing=0.3)


def main():
    fig, axes = plt.subplots(
        2, 2, figsize=(11.0, 6.6), sharex="col",
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.0,
                     "wspace": 0.24})
    panel(axes[0][0], axes[1][0], "results", MU_ROWS,
          "Muon DIS", r"$E_{\mu}$  [GeV]", 1000.0, "nb")
    panel(axes[0][1], axes[1][1], "results_nu", NU_ROWS,
          "Neutrino DIS", r"$E_{\nu}$  [GeV]", 1.0, "pb")
    axes[1][0].set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                          fontsize=plotstyle.FS_YLABEL - 1)
    for a in (axes[1][0], axes[1][1]):
        lo, hi = a.dataLim.intervaly
        half = max(abs(hi - 1.0), abs(1.0 - lo)) * 1.14
        # extra room on top: the NNLO MHOU legend sits upper right
        a.set_ylim(1.0 - half, 1.0 + 1.45 * half)
    fig.suptitle(tex(TITLE), y=0.98, fontsize=plotstyle.FS_SUPTITLE)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
