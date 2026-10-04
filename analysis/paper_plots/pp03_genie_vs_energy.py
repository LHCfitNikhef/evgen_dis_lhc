#!/usr/bin/env python3
"""Paper plot 3: GENIE against the analytic reference versus beam energy,
in the FINAL benchmark region.

The redo of the earlier figure for the "Paper plots"
tab (user, 2026-09-14: "produce now of pp_03 and pp_04").  The region
and conventions are pp01's (PRODUCTION.md): tungsten per nucleon,
(74 sigma_p + 110 sigma_n)/184 with every generator run on the proton and the
neutron separately; Q2 > 4 GeV2 and W > 3 GeV with no y cut on BOTH currents;
400, 700, 1000, 2000 and 4000 GeV; references with target-mass corrections.

Everything else is the earlier production's figure: the same GENIE rows, the same POWHEG entries
for scale, and the same reference.

  * THE REFERENCE IS FONLL, as in the earlier production: GENIE carries charm masses, so the
    general-mass calculation is the like-for-like reference for it.
  * THE GENIE ROWS (CONVENTIONS.md rule 1b: both GENIE rows count):
      - G18_02a on GRV98, the classic tune and the one FASER runs;
      - G18_02a with NNPDF4.0 substituted into Bodek-Yang, the diagnostic
        that isolates the parton distribution;
      - HEDIS GHE19_00a on BGR18, neutrino only.
    the earlier production's configurations with the target as the only change
    (genie/production/genie_job.sh), each nucleon with its own cross-section spline.
  * G18_02a IS NOT DRAWN ABOVE 1 TeV (user, 2026-09-01), on either current:
    the tune declares validity to 1000 GeV.  HEDIS keeps every energy.
  * THE NLO BAND is the reference's seven-point variation with pp01's
    prescription: at Q2min = 4 the xiF = 1/2 points leave the PDF grid and
    are dropped, and the muon band is symmetrised.

POWHEG-RES and POWHEG-V2 are drawn, for scale, as in the earlier production; Herwig and Sherpa are
the subject of paper plot 1 and are not repeated here.

Usage: analysis/paper_plots/pp03_genie_vs_energy.py
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

SLUG = "pp03_genie_vs_energy"
V1_NUMBER = 3
# IN THE PAPER in place of the earlier production's figure (user, 2026-09-14: every updated 
# plot goes into the paper).  Read by tools/check_paper_figures_used.py.
IN_PAPER = True
TITLE = "GENIE against the NLO matchings: total cross-section vs beam energy"
# the title ON the figure is pp01's (user, 2026-09-18: "change the title to
# the same as pp_01"); TITLE above names the report tab
FIG_TITLE = "Fiducial NLO total cross-section vs beam energy"
OUTPUT = "pp_genie_vs_energy.png"
RESULTS = "results_nu"
CAPTION = ("Fiducial cross-section per nucleon on tungsten against beam "
           "energy, for Q&sup2; &gt; 4 GeV&sup2; and W &gt; 3 GeV with no y "
           "cut, for the GENIE configurations: muon neutral current (left) "
           "and neutrino charged current (right), with POWHEG-RES and "
           "POWHEG-V2 for scale. Ratio to YADISM (FONLL) with target-mass "
           "corrections beneath, with the reference's seven-point scale "
           "band. The G18_02a curves stop at 1 TeV, the tune's declared "
           "validity. The dashed orange curve on the neutrino panel is the "
           "same reference evaluated with the parton set GENIE HEDIS itself "
           "uses, NNPDF3.1sx rather than NNPDF4.0; the HEDIS points lie on "
           "it, so their offset from the benchmark reference is that choice "
           "of parton distribution and not a disagreement between two NLO "
           "calculations.")

SURFACE = "#ffffff"
REGION = "q4w3"
ES = list(beams.BENCH_ENERGIES)
Q2MIN = 4.0
W_MIN = 3.0
TMC = "_tmc"
REF = f"yadism_nlo_fonll_damp_{REGION}_W{TMC}"
REF_LABEL = "YADISM (FONLL)"
MHOU = {"results": f"mhou_sigma_fonll_mu_{REGION}_W{TMC}.json",
        "results_nu": f"mhou_sigma_fonll_{REGION}_W{TMC}.json"}

# >>> THE SAME REFERENCE, WITH GENIE HEDIS'S OWN PARTON DENSITIES (user,
# 2026-09-21: "show this agreement explicitely ... GENIE HEDIS can reproduce
# the analytical reference, once settings are made uniform"). <<<
# GHE19_00a evaluates its structure functions with NNPDF3.1sx, not NNPDF4.0,
# and that is the whole of its 5% offset.  Drawing the reference a SECOND
# time, in GENIE's set and in the same FONLL scheme with the same target-mass
# corrections, puts the HEDIS points on a curve instead of below one --
# which is the difference between "GENIE disagrees with NLO QCD" and "GENIE
# uses a different PDF", and only the second is true.
# It is written by tools/genie_hedis_pdf.sh; the decomposition behind it is
# analysis/genie_hedis_pdf.py.  NEUTRINO PANEL ONLY, because HEDIS is a
# charged-current module and has no muon counterpart to draw.
HEDIS_PDF = "NNPDF31sx_nlo_as_0118_LHCb_nf_6"
HEDIS_REF_FILE = (f"mhou_sigma_fonll_{REGION}{TMC}_"
                  "nnpdf31sxnloas0118lhcbnf6w184free.json")
HEDIS_REF_LABEL = "YADISM (FONLL), HEDIS PDF"
HEDIS_REF_COLOUR = "#444444"      # grey, the YADISM family (user, 2026-10-04: in GENIE (HEDIS)'s orange the two could not be told apart)


def hedis_ref(resdir, e):
    """The reference recomputed with GENIE HEDIS's PDF, pb per nucleon."""
    if resdir != "results_nu":
        return None
    p = f"{BASE}/{resdir}/{HEDIS_REF_FILE}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    r = d["points"].get(f"{e:g}")
    return r["central_pb"] if r else None

# (label, result key, colour, marker, linestyle) -- the earlier production's rows and styling
MU_ROWS = [
    ("GENIE (GRV98LO)",          f"genie_{REGION}_W",        "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", f"genie_nnpdf_{REGION}_W",  "#56b4e9", "D", "--"),
    ("POWHEG-RES",               f"powheg_{REGION}_W",       "#0072b2", "s", ":"),
]
NU_ROWS = [
    ("GENIE (GRV98LO)",          f"genie_lo_{REGION}_W",     "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)", f"genie_nnpdf_{REGION}_W",  "#56b4e9", "D", "--"),
    ("GENIE (HEDIS)",            f"genie_{REGION}_W",        "#e69f00", "^", "-."),
    ("POWHEG-V2",                f"powheg_nu_{REGION}_W",    "#0072b2", "s", ":"),
]

# G18_02a is declared valid to 1000 GeV -- a property of the TUNE, so it binds
# both currents; `genie` on the neutrino side is HEDIS and is not limited.
EMAX = {(f"genie_{REGION}_W", "results"): 1000.0,
        (f"genie_nnpdf_{REGION}_W", "results"): 1000.0,
        (f"genie_lo_{REGION}_W", "results_nu"): 1000.0,
        (f"genie_nnpdf_{REGION}_W", "results_nu"): 1000.0}


def _in_range(resdir, key, e):
    lim = EMAX.get((key, resdir))
    return lim is None or e <= lim + 1e-9


def sigma(resdir, key, e):
    """(sigma_fid, err) in pb per nucleon, or None.

    Every file is CHECKED to be the q4w3 region on tungsten, and a reference
    to carry TMC and the W184free set; outside a tune's validity nothing is
    read at all.
    """
    if not _in_range(resdir, key, e):
        return None
    p = f"{BASE}/{resdir}/histos_{beams.at_energy(key, e)}.json"
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
    if key.startswith("yadism"):
        if d.get("tmc") != target.TMC_CHOSEN or d.get("pdf_set") != target.PDFSET_W:
            raise SystemExit(f"{rel}: tmc={d.get('tmc')} pdf={d.get('pdf_set')}")
    if d.get("stat_insufficient"):
        return None
    return d.get("sigma_fid_pb"), d.get("sigma_fid_err_pb") or 0.0


def _r(resdir, key, e):
    v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
    return None if not (v and r) else v[0] / r[0]


def _band_file(resdir):
    p = f"{BASE}/{resdir}/{MHOU[resdir]}"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if (d.get("target") != "W" or d.get("tmc") != target.TMC_CHOSEN
            or (d.get("selection") or {}).get("name") != REGION
            or not d.get("scheme", "").startswith("FONLL")):
        raise SystemExit(f"{resdir}/{MHOU[resdir]}: not a FONLL tungsten "
                         f"TMC-on {REGION} band")
    return d


def scale_band(resdir, e):
    """(rel_hi, rel_lo) of the NLO envelope, pp01's prescription."""
    d = _band_file(resdir)
    if d is None or f"{e:g}" not in d["points"]:
        return None
    pt = d["points"][f"{e:g}"]
    q0 = d.get("q2min_pdf")
    drop = q0 is None or Q2MIN / 4.0 < q0 - 1e-9
    v = [r["sigma_fid_pb"] for r in pt["rows"] if not (drop and r["xiF"] == 0.5)]
    hi, lo = max(v) / pt["central_pb"] - 1.0, min(v) / pt["central_pb"] - 1.0
    if drop and resdir == "results":
        dd = max(abs(hi), abs(lo))
        hi, lo = dd, -dd
    return hi, lo


# ----------------------------------------------------------------- claims --
_GENIE = [("results", [k for _l, k, *_ in MU_ROWS[:2]]),
          ("results_nu", [k for _l, k, *_ in NU_ROWS[:3]])]


def _es(resdir, key):
    return [e for e in ES if _in_range(resdir, key, e)]


def _band_edge(resdir, key, e):
    """How far a point sits outside the NLO band: 0 inside, else the excess."""
    hi, lo = scale_band(resdir, e)
    r = _r(resdir, key, e) - 1.0
    return max(0.0, r - hi, lo - r)


MU_GRV, MU_NN, MU_PWG = (k for _l, k, *_ in MU_ROWS)
NU_GRV, NU_NN, NU_HEDIS, NU_PWG = (k for _l, k, *_ in NU_ROWS)


def _rs(resdir, key):
    return [_r(resdir, key, e) for e in _es(resdir, key)]


def _fmt(resdir, key):
    return ", ".join("%g: %.4f" % (e, _r(resdir, key, e)) for e in _es(resdir, key))


MESSAGE = (
    "<b>This is the figure a FASER reader comes for, now in the final "
    "region</b>: tungsten per nucleon, Q&sup2; &gt; 4 GeV&sup2; and W &gt; "
    "3 GeV with no y cut, every GENIE configuration run on the proton and "
    "the neutron separately. POWHEG-RES and POWHEG-V2 are drawn for scale "
    "and sit within 2% of the reference at every energy on both currents; "
    "the reference is FONLL, with target-mass corrections, because GENIE "
    "carries charm masses. The G18_02a rows stop at 1 TeV, the tune's "
    "declared validity; HEDIS keeps every energy."
    "\n\n"
    "<b>On the charged current GENIE's default tune agrees with the NLO "
    "reference to 0.6% at every energy where it is valid</b>, inside the "
    "reference's own scale band, <b>and that agreement is a cancellation</b>: "
    "substituting NNPDF4.0 into the same Bodek-Yang "
    "machinery moves it 3.5&ndash;3.7% high, outside the band, so the parton "
    "distribution alone is worth four per cent. <b>GENIE HEDIS, the NLO "
    "configuration, is a flat 5% low</b> (0.946 to 0.951 from 400 GeV to "
    "4 TeV), well below the band."
    "\n\n"
    "<b>The dashed orange curve is that same NLO reference evaluated with "
    "HEDIS's own parton set, and the HEDIS points lie on it.</b> "
    "GHE19_00a evaluates its structure functions with NNPDF3.1sx (the NLO "
    "LHCb nf&thinsp;=&thinsp;6 member) where this benchmark uses NNPDF4.0 at "
    "NNLO, and that one difference is the whole of the 5%: the curve sits at "
    "0.947 to 0.951 from 400 GeV to 4 TeV against the 0.946 to 0.951 the "
    "points show, reproducing even the slow rise with energy. With the "
    "parton set matched, GENIE HEDIS and the NLO calculation agree to 0.2% "
    "at every energy and to 0.04% from 1 TeV up, inside the GENIE "
    "sample&rsquo;s own statistics. What little is left grows towards low "
    "energy &mdash; 0.16% at 400 GeV &mdash; and that is the last "
    "unmatched setting: the reference&rsquo;s FONLL is threshold-damped and "
    "APFEL&rsquo;s FONLL-B is not, which is a low-Q&sup2; difference. <b>It "
    "is a choice of parton distribution, not a disagreement between two NLO "
    "calculations.</b>"
    "\n\n"
    "It is none of the other candidates, each measured rather than assumed: "
    "not the perturbative order, which is NLO on both sides; not alpha_s, "
    "which HEDIS takes from its own PDF at the same 0.1180; not the "
    "heavy-quark masses, worth 0.4% here; not the target-mass corrections, "
    "worth 0.01%; and not an artefact of the structure-function grid, which "
    "the agreement leaves no room for. On the tungsten nucleon the HEDIS set "
    "carries 4&ndash;8% fewer d&thinsp;+&thinsp;s&thinsp;+&thinsp;b quarks "
    "than NNPDF4.0 over the 0.02&thinsp;&lt;&thinsp;x&thinsp;&lt;&thinsp;0.5 "
    "that holds the rate, crossing above only past "
    "x&thinsp;&asymp;&thinsp;0.6. The decomposition in "
    "<code>results_nu/genie_hedis_pdf.json</code> puts 0.6% on the order of "
    "the fit, 1.8% on NNPDF3.1 against NNPDF4.0 and 0.1% on small-x "
    "resummation, so what is left &mdash; and it is most of the effect "
    "&mdash; is the sx set itself."
    "\n\n"
    "<b>On the neutral current the default tune is not flat in energy</b>: "
    "0.998 at 400 GeV, 0.985 at 700 GeV and 0.973 at 1 TeV, all inside a "
    "band that is 4&ndash;5% wide here. With NNPDF4.0 it is 6.3% high at "
    "400 GeV, falling to 3.2% at 1 TeV, above the band at 400 GeV and inside "
    "it at 1 TeV."
)

def _hedis_pdf():
    """The GENIE-HEDIS-PDF cross-check (analysis/genie_hedis_pdf.py)."""
    import json as _json
    p = f"{BASE}/results_nu/genie_hedis_pdf.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return _json.load(f)


CLAIMS = [
    {"what": "every file drawn is a tungsten per-nucleon q4w3 result, and the "
             "reference is FONLL with TMC",
     "check": lambda: (REF.startswith("yadism_nlo_fonll_damp") and REF.endswith("_tmc")
                       and all(sigma(d, k, e) is not None
                               for d, rows in (("results", MU_ROWS),
                                               ("results_nu", NU_ROWS))
                               for _l, k, *_ in rows for e in _es(d, k))),
     "detail": lambda: REF},
    {"what": "G18_02a rows stop at 1 TeV on both currents; HEDIS and POWHEG "
             "keep all five energies",
     "check": lambda: (all(max(_es(d, k)) == 1000.0 for d, k in
                           (("results", MU_GRV), ("results", MU_NN),
                            ("results_nu", NU_GRV), ("results_nu", NU_NN)))
                       and _es("results_nu", NU_HEDIS) == ES
                       and all(sigma(d, k, 4000.0) is None for d, k in
                               (("results", MU_GRV), ("results_nu", NU_GRV)))),
     "detail": lambda: "HEDIS energies %s" % _es("results_nu", NU_HEDIS)},
    {"what": "both POWHEG entries within 2% of the reference at every energy",
     "check": lambda: (max(abs(x - 1) for x in _rs("results", MU_PWG)) < 0.02
                       and max(abs(x - 1) for x in _rs("results_nu", NU_PWG)) < 0.02),
     "detail": lambda: "mu " + _fmt("results", MU_PWG) + " | nu " + _fmt("results_nu", NU_PWG)},
    {"what": "CC default tune within 0.6% of the reference at 400/700/1000 GeV "
             "(0.4% before the 2026-10-01 charm-propagator fix), "
             "inside the scale band",
     "check": lambda: (max(abs(x - 1) for x in _rs("results_nu", NU_GRV)) < 0.006
                       and all(_band_edge("results_nu", NU_GRV, e) == 0
                               for e in _es("results_nu", NU_GRV))),
     "detail": lambda: _fmt("results_nu", NU_GRV)},
    {"what": "and with NNPDF4.0 the same tune is 3.5-3.7% high, outside the band",
     "check": lambda: (all(0.034 <= x - 1 <= 0.038 for x in _rs("results_nu", NU_NN))
                       and all(_band_edge("results_nu", NU_NN, e) > 0
                               for e in _es("results_nu", NU_NN))),
     "detail": lambda: _fmt("results_nu", NU_NN)},
    {"what": "the dashed curve is drawn at every energy and the HEDIS points "
             "sit on it: the reference in the SAME scheme with GENIE's own "
             "PDF is within 0.2% of the delivered HEDIS cross-section at all "
             "five energies, and within 0.05% from 1 TeV up",
     "check": lambda: (_hedis_pdf() is not None
                       and all(hedis_ref("results_nu", e) is not None for e in ES)
                       and all(abs(r["genie_over_fonll_same_pdf"] - 1.0) < 0.002
                               for r in _hedis_pdf()["rows"])
                       and all(abs(r["genie_over_fonll_same_pdf"] - 1.0) < 0.0005
                               for r in _hedis_pdf()["rows"]
                               if r["energy_gev"] >= 1000)),
     "detail": lambda: ", ".join(
         f"{r['energy_gev']:g} {r['genie_over_fonll_same_pdf']:.4f}"
         for r in _hedis_pdf()["rows"]) if _hedis_pdf() else "missing"},
    {"what": "the HEDIS offset IS its PDF: the same YADISM (ZM-VFNS) "
             "calculation with NNPDF3.1sx reproduces the published ratio to "
             "0.2% at every energy, and with the flavour scheme matched too "
             "GENIE HEDIS agrees with the NLO reference to 0.1%",
     "check": lambda: (_hedis_pdf() is not None
                       and all(abs(r["pdf_only"] - r["ratio_published"]) < 0.002
                               for r in _hedis_pdf()["rows"])
                       and all(abs(r["genie_over_like_for_like"] - 1.0) < 0.001
                               for r in _hedis_pdf()["rows"])
                       and _hedis_pdf()["hedis_pdf"].startswith("NNPDF31sx")
                       and all(0.92 <= v <= 0.97
                               for k, v in _hedis_pdf()["quark_ratio_x"].items()
                               if 0.02 <= float(k) <= 0.5)
                       and _hedis_pdf()["quark_ratio_x"]["0.7"] > 1.0
                       and abs(_hedis_pdf()["fonll_factorisation_check"]
                               ["genie_over_direct_fonll"] - 1.0) < 0.0006),
     "detail": lambda: (", ".join(
         f"{r['energy_gev']:g} {r['pdf_only']:.4f} vs {r['ratio_published']:.4f}"
         f" -> {r['genie_over_like_for_like']:.4f}"
         for r in _hedis_pdf()["rows"]) if _hedis_pdf() else
         "results_nu/genie_hedis_pdf.json missing "
         "(run tools/genie_hedis_pdf.sh)")},
    {"what": "HEDIS a flat 5% low, 0.946-0.951 from 400 GeV to 4 TeV, below the band",
     "check": lambda: (all(0.945 <= x <= 0.952 for x in _rs("results_nu", NU_HEDIS))
                       and all(_band_edge("results_nu", NU_HEDIS, e) > 0 for e in ES)),
     "detail": lambda: _fmt("results_nu", NU_HEDIS)},
    {"what": "NC default tune 0.998, 0.985, 0.973 -- falling with energy, "
             "inside the band",
     "check": lambda: (all(abs(_r("results", MU_GRV, e) - v) < 0.0015
                           for e, v in ((400.0, 0.998), (700.0, 0.985), (1000.0, 0.973)))
                       and all(_band_edge("results", MU_GRV, e) == 0
                               for e in _es("results", MU_GRV))),
     "detail": lambda: _fmt("results", MU_GRV)},
    {"what": "NC with NNPDF4.0: 6.3% high at 400 GeV (above the band) falling "
             "to 3.2% at 1 TeV (inside it)",
     "check": lambda: (abs(_r("results", MU_NN, 400.0) - 1.063) < 0.002
                       and abs(_r("results", MU_NN, 1000.0) - 1.032) < 0.002
                       and _band_edge("results", MU_NN, 400.0) > 0
                       and _band_edge("results", MU_NN, 1000.0) == 0),
     "detail": lambda: _fmt("results", MU_NN)},
]


def panel(ax, axr, resdir, rows, title, xlab, div, unit):
    pts = [(e, sigma(resdir, REF, e)) for e in ES]
    pts = [(e, v[0] / div) for e, v in pts if v]
    if len(pts) >= 2:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ":",
                color="#111111", lw=1.6, label=tex(REF_LABEL), zorder=1)
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
        ax.fill_between([q[0] for q in ab],
                        [c * (1.0 + b[1]) for _e, b, c in ab],
                        [c * (1.0 + b[0]) for _e, b, c in ab],
                        color="#111111", alpha=0.12, lw=0, zorder=0,
                        label=tex("NLO MHOU"))
    # the same calculation with GENIE HEDIS's set, drawn before the points so
    # the HEDIS markers sit ON it rather than under it
    hp = [(e, hedis_ref(resdir, e), sigma(resdir, REF, e)) for e in ES]
    hp = [(e, v / div, v / r[0]) for e, v, r in hp if v and r]
    if len(hp) >= 2:
        # A REGULAR DASHED LINE, dark grey and ON TOP of the GENIE (HEDIS)
        # curve it agrees with to 0.2% (user, 2026-10-04: the wide halo
        # "looks funny")
        ax.plot([q[0] for q in hp], [q[1] for q in hp], "--",
                color=HEDIS_REF_COLOUR, lw=1.7,
                label=tex(HEDIS_REF_LABEL), zorder=6)
        axr.plot([q[0] for q in hp], [q[2] for q in hp], "--",
                 color=HEDIS_REF_COLOUR, lw=1.7, zorder=6)
    for lab, key, colour, marker, ls in rows:
        pp, rat = [], []
        for e in ES:
            v, r = sigma(resdir, key, e), sigma(resdir, REF, e)
            if not v:
                continue
            pp.append((e, v[0] / div, v[1] / div))
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
    _seen, _hh, _ll = set(), [], []
    for h, s in zip(_h, _l):
        if s not in _seen:
            _seen.add(s)
            _hh.append(h)
            _ll.append(s)
    # the band shares the reference's entry (user, 2026-10-04)
    _hh, _ll = plotstyle.merge_bands(_hh, _ll,
                                     {tex("NLO MHOU"): [tex(REF_LABEL)]})
    ax.legend(_hh, _ll, fontsize=plotstyle.FS_LEGEND - 1, frameon=True,
              loc="upper left", handlelength=1.9, labelspacing=0.32)


def main():
    fig, axes = plt.subplots(
        2, 2, figsize=(11.0, 6.6), sharex="col",
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.0,
                     "wspace": 0.24})
    # titles are the current alone (user, 2026-09-14, on pp01)
    panel(axes[0][0], axes[1][0], "results", MU_ROWS, "Muon DIS",
          r"$E_{\mu}$  [GeV]", 1000.0, "nb")
    panel(axes[0][1], axes[1][1], "results_nu", NU_ROWS, "Neutrino DIS",
          r"$E_{\nu}$  [GeV]", 1.0, "pb")
    for c in (0, 1):
        axes[1][c].set_ylabel(tex("ratio to") + "\n" + tex(REF_LABEL),
                              fontsize=plotstyle.FS_YLABEL - 1)
        axes[0][c].yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=6, prune="lower"))
        axes[1][c].margins(y=0.12)   # no marker on the fused spine
        axes[1][c].yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=5, prune="upper"))
    fig.suptitle(tex(FIG_TITLE), fontsize=plotstyle.FS_TITLE, y=0.98)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
