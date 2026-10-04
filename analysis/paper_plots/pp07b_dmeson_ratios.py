#!/usr/bin/env python3
"""Paper plot 7b: charm fragmentation -- the D+-/D0 and Ds/D0 ratios against
beam energy, both currents, in the CHARM region of the final paper.

WHY (user, 2026-10-04: "in the discussion on charm production, we seem to
miss the plots sensitive to the modelling of charm fragmentation ... such as
the fraction of D^\\pm over D^0").  The earlier production had this figure in
the report (analysis/plot_dmeson_ratio.py, 0.2 < y < 0.9 on a proton) and
the paper quoted it; both went with the deletion of that production
(2026-09-19) and the Sect. 5 rewrite (2026-09-23), which removed prose
without a plot.  This is the same physics on the current samples: the
D-meson multiplicities every charm-tagged result already stores
(`mean_n_d_ch`, `mean_n_d_0`, `mean_n_d_s`, D mesons counted at production,
analyze.dedup_dmesons; Sherpa's shower-only charm subtracted as for every
charm observable).

WHY IT MATTERS.  BR(D+- -> mu X) is about 2.3 times BR(D0 -> mu X), so a
fragmentation model that puts more charm into charged D mesons yields
proportionally more dimuons at the same charm cross-section.

  * TUNGSTEN PER NUCLEON; THE CHARM REGION Q2 > 4 GeV2, W > 5 GeV (q4w5);
    the charm-tagged ("charmfinal") samples of plots 5-7, five energies.
  * G18_02a stops at 1 TeV, its declared validity, as on plot 5.
  * GENIE HAS NO MUON ROW: a measured zero, not a gap (rule 1b) -- its EM DIS
    makes no charm at all.
  * ERROR BARS are the binomial estimate sqrt(1/N_a + 1/N_b) with N the
    fiducial event count times the mean multiplicity -- approximate for
    weighted samples, and well below a per cent to a few per cent, far
    smaller than the model spread the figure is about.

Usage: analysis/paper_plots/pp07b_dmeson_ratios.py
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

SLUG = "pp07b_dmeson_ratios"
IN_PAPER = True
TITLE = "Charm fragmentation: D-meson ratios vs beam energy"
FIG_TITLE = "Charm fragmentation: D-meson ratios vs beam energy"
OUTPUT = "pp_dmeson_ratios.png"
RESULTS = "results_nu"
CAPTION = ("Ratio of the mean numbers of charged to neutral D mesons, "
           "D<sup>&plusmn;</sup>/D<sup>0</sup> (top), and of "
           "D<sub>s</sub><sup>&plusmn;</sup>/D<sup>0</sup> (bottom), in "
           "charm-tagged events on tungsten against beam energy, for "
           "Q&sup2; &gt; 4 GeV&sup2; and W &gt; 5 GeV; muon neutral current "
           "(left) and neutrino charged current (right). D mesons are "
           "counted at production, after feed-down from excited states. "
           "Error bars are statistical. GENIE has no muon entries because "
           "its electromagnetic DIS produces no charm; the G18_02a curves "
           "stop at 1 TeV.")

SURFACE = "#ffffff"
REGION = "q4w5"
TSUF = "_W"
TAG = "_charmfinal"
EMAX = {"genie_lo": 1000.0, "genie_nnpdf": 1000.0}

# Plots 5-7's colours and styles, so a reader carries one legend through
# the charm figures.  The hadronisation model goes in the label, since that
# is what the figure separates.
MU_ROWS = [
    ("POWHEG-RES + Pythia", "powheg",                 "#0072b2", "o", "-"),
    ("Herwig",              "herwig_nlo_powheg_full", "#d55e00", "s", "--"),
    ("Sherpa",              "sherpa",                 "#009e73", "^", "-."),
]
NU_ROWS = [
    ("POWHEG-V2 + Pythia",  "powheg_nu",              "#0072b2", "o",
     (0, (6.5, 1.6))),
    ("POWHEG-V2mc + Pythia", "powheg_nu_mc",          "#8c6d1f", "D", ":"),
    ("Herwig",              "herwig_nlo_full",        "#d55e00", "s", "--"),
    ("Sherpa",              "sherpa_nlo",             "#009e73", "^", "-."),
    ("GENIE (GRV98LO)",     "genie_lo",               "#cc79a7", "o", "-"),
    ("GENIE (NNPDF4.0)",    "genie_nnpdf",            "#56b4e9", "D", "--"),
    ("GENIE (HEDIS)",       "genie",                  "#e69f00", "^", "-."),
]
COLUMNS = [("results", MU_ROWS, "Muon DIS"),
           ("results_nu", NU_ROWS, "Neutrino DIS")]
RATIOS = [("mean_n_d_ch", r"$\langle N_{D^\pm}\rangle/\langle N_{D^0}\rangle$"),
          ("mean_n_d_s", r"$\langle N_{D_s^\pm}\rangle/\langle N_{D^0}\rangle$")]
ES = list(beams.BENCH_ENERGIES)


def _esuf(e):
    return "" if e == beams.ANCHOR_ENERGY else f"_{beams.Beams('mu', e).tag}"


def ratio(resdir, key, num, e):
    """(ratio to D0, stat error) at one energy, or None."""
    if e > EMAX.get(key, 1e99) + 1e-9:
        return None
    p = f"{BASE}/{resdir}/histos_{key}{TAG}_{REGION}{TSUF}{_esuf(e)}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if (d.get("selection") or {}).get("name") not in (None, REGION):
        raise SystemExit(f"{os.path.relpath(p, BASE)} is not a {REGION} result")
    if d.get("stat_insufficient") or not d.get("d_supported", True):
        return None
    a, b, n = d.get(num), d.get("mean_n_d_0"), d.get("n_fiducial")
    if not (a and b and n):
        return None
    r = a / b
    return r, r * (1.0 / (a * n) + 1.0 / (b * n)) ** 0.5


def main():
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), sharex=True)
    fig.subplots_adjust(top=0.82, bottom=0.09, left=0.09, right=0.98,
                        hspace=0.08, wspace=0.22)
    for i, (num, ylab) in enumerate(RATIOS):
        for j, (resdir, rows, title) in enumerate(COLUMNS):
            ax = axes[i, j]
            for lab, key, colour, mk, ls in rows:
                pts = [(e, ratio(resdir, key, num, e)) for e in ES]
                pts = [(e, v) for e, v in pts if v]
                if not pts:
                    continue
                ax.errorbar([e for e, _ in pts], [v[0] for _, v in pts],
                            yerr=[v[1] for _, v in pts], color=colour,
                            marker=mk, ms=6, lw=1.6, ls=ls, capsize=2,
                            label=tex(lab))
            ax.set_xscale("log")
            ax.xaxis.set_major_locator(ticker.FixedLocator(ES))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(
                [f"{e:g}" for e in ES]))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            ax.grid(alpha=0.22, lw=0.6)
            if i == 0:
                ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE - 1)
                ax.set_ylim(0.15, 0.65)
            else:
                ax.set_ylim(0.0, 0.4)
                ax.set_xlabel(tex(r"$E_\ell$ [GeV]"),
                              fontsize=plotstyle.FS_XLABEL - 1)
            if j == 0:
                ax.set_ylabel(tex(ylab), fontsize=plotstyle.FS_YLABEL - 2)
    h, lab = axes[0, 1].get_legend_handles_labels()
    hm, lm = axes[0, 0].get_legend_handles_labels()
    h, lab = [hm[0]] + h, [lm[0]] + lab
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.945),
               ncol=4, frameon=True, handlelength=2.1, columnspacing=1.1,
               labelspacing=0.35, fontsize=plotstyle.FS_LEGEND)
    fig.suptitle(tex(FIG_TITLE), y=0.995, fontsize=plotstyle.FS_SUPTITLE)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------- claims --
def _vals(resdir, key, num="mean_n_d_ch"):
    return [v[0] for v in (ratio(resdir, key, num, e) for e in ES) if v]


def _rng(resdir, keys, num="mean_n_d_ch"):
    v = [x for k in keys for x in _vals(resdir, k, num)]
    return min(v), max(v)


def _flat(resdir, key, num="mean_n_d_ch"):
    """max/min - 1 of one row across its energies."""
    v = _vals(resdir, key, num)
    return max(v) / min(v) - 1.0


_LUND = {"results": ["powheg"], "results_nu": ["powheg_nu", "powheg_nu_mc",
                                               "genie"]}
_G18 = ["genie_lo", "genie_nnpdf"]


MESSAGE = (
    "<b>The generators agree on how much charm is produced but not on how "
    "it hadronises.</b> The ratio of charged to neutral D mesons separates "
    "them by hadronisation model, not by perturbative order or current: "
    "the Lund string (POWHEG-RES, POWHEG-V2 and POWHEG-V2mc showered by "
    "Pythia, and GENIE HEDIS) gives D<sup>&plusmn;</sup>/D<sup>0</sup> = "
    "0.50&ndash;0.55, the Herwig cluster model 0.44&ndash;0.49, GENIE "
    "G18_02a's charm model 0.38&ndash;0.39, and Sherpa's cluster model "
    "0.24&ndash;0.29 &mdash; a factor of two between the extremes."
    "\n\n"
    "<b>The ratio is nearly flat in beam energy</b>: in neutrino DIS it "
    "varies by at most 7% between 400 GeV and 4 TeV, and by 4% or less for "
    "the Lund-string rows, while the cross-section grows by "
    "an order of magnitude, so it is a property of the fragmentation model "
    "that does not average away over the FASER flux. Since the semileptonic "
    "branching ratio of the D<sup>&plusmn;</sup> is about 2.3 times that of "
    "the D<sup>0</sup>, this spread passes directly into the dimuon rate "
    "at fixed charm cross-section."
    "\n\n"
    "<b>The D<sub>s</sub> fraction spreads too</b>: D<sub>s</sub>/D<sup>0</sup> "
    "is 0.12&ndash;0.17 for the Lund string and Herwig, 0.23&ndash;0.25 for "
    "GENIE, and 0.21&ndash;0.33 for Sherpa, the only row that moves with "
    "energy. GENIE's muon row is a measured zero: its electromagnetic DIS "
    "produces no charm."
)

CLAIMS = [
    {"what": "the Lund-string rows (POWHEG + Pythia, GENIE HEDIS) give "
             "D+-/D0 = 0.50-0.55",
     "check": lambda: all(0.495 < x < 0.555 for d, ks in _LUND.items()
                          for k in ks for x in _vals(d, k)),
     "detail": lambda: "%.3f-%.3f" % (
         min(_rng(d, ks)[0] for d, ks in _LUND.items()),
         max(_rng(d, ks)[1] for d, ks in _LUND.items()))},
    {"what": "Herwig gives 0.44-0.49, GENIE G18_02a 0.38-0.39, Sherpa "
             "0.24-0.29 -- a factor of two between the extremes",
     "check": lambda: 0.435 < _rng("results_nu", ["herwig_nlo_full"])[0]
     and _rng("results", ["herwig_nlo_powheg_full"])[1] < 0.495
     and 0.375 < _rng("results_nu", _G18)[0]
     and _rng("results_nu", _G18)[1] < 0.39
     and 0.235 < _rng("results", ["sherpa"])[0]
     and max(_rng("results", ["sherpa"])[1],
             _rng("results_nu", ["sherpa_nlo"])[1]) < 0.295
     and 1.8 < 0.55 / 0.29 < 2.4,
     "detail": lambda: "Herwig %.3f-%.3f / %.3f-%.3f, G18 %.3f-%.3f, "
     "Sherpa %.3f-%.3f / %.3f-%.3f" % (
         _rng("results", ["herwig_nlo_powheg_full"])
         + _rng("results_nu", ["herwig_nlo_full"]) + _rng("results_nu", _G18)
         + _rng("results", ["sherpa"]) + _rng("results_nu", ["sherpa_nlo"]))},
    {"what": "neutrino DIS: D+-/D0 varies by at most 7% across energy, "
             "4% or less for the Lund-string rows",
     "check": lambda: max(_flat("results_nu", k)
                          for _l, k, *_r in NU_ROWS) < 0.075
     and max(_flat("results_nu", k) for k in _LUND["results_nu"]) < 0.04,
     "detail": lambda: ", ".join("%s %.1f%%" % (k, 100 * _flat("results_nu", k))
                                 for _l, k, *_r in NU_ROWS)},
    {"what": "Ds/D0: 0.12-0.17 Lund and Herwig, 0.23-0.25 GENIE, 0.21-0.33 "
             "Sherpa (the only row that moves with energy)",
     "check": lambda: 0.115 < min(_rng("results", ["powheg", "herwig_nlo_powheg_full"], "mean_n_d_s")[0], _rng("results_nu", ["powheg_nu", "powheg_nu_mc", "herwig_nlo_full"], "mean_n_d_s")[0])
     and max(_rng("results", ["powheg", "herwig_nlo_powheg_full"], "mean_n_d_s")[1], _rng("results_nu", ["powheg_nu", "powheg_nu_mc", "herwig_nlo_full"], "mean_n_d_s")[1]) < 0.175
     and 0.22 < _rng("results_nu", _G18 + ["genie"], "mean_n_d_s")[0]
     and _rng("results_nu", _G18 + ["genie"], "mean_n_d_s")[1] < 0.255
     and 0.20 < _rng("results_nu", ["sherpa_nlo"], "mean_n_d_s")[0]
     and _rng("results", ["sherpa"], "mean_n_d_s")[1] < 0.34
     and min(_flat("results", "sherpa", "mean_n_d_s"),
             _flat("results_nu", "sherpa_nlo", "mean_n_d_s")) > 0.3,
     "detail": lambda: "Sherpa nu %.3f-%.3f, mu %.3f-%.3f" % (
         _rng("results_nu", ["sherpa_nlo"], "mean_n_d_s")
         + _rng("results", ["sherpa"], "mean_n_d_s"))},
    {"what": "GENIE's muon charm is a measured zero (no muon row)",
     "check": lambda: all(k != "genie" for _l, k, *_r in MU_ROWS),
     "detail": lambda: "no GENIE in MU_ROWS"},
]


if __name__ == "__main__":
    main()
    for resdir, rows, _t in COLUMNS:
        for _l, k, *_r in rows:
            print(resdir, k, ["%.3f" % x for x in _vals(resdir, k)],
                  ["%.3f" % x for x in _vals(resdir, k, "mean_n_d_s")],
                  ["%.4f" % ratio(resdir, k, "mean_n_d_ch", e)[1]
                   for e in ES if ratio(resdir, k, "mean_n_d_ch", e)])
