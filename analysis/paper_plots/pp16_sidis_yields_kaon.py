#!/usr/bin/env python3
"""Paper plot 16: the expected charged-kaon yields in z at FASERnu,
300 fb^-1, in the region.

The kaon companion of pp15 (user, 2026-09-21: "the same but for identified
kaons in the final state").  Identified charged kaons -- |PDG| = 321, the
species analysis/faser_pions.py calls "K" -- inside the emulsion track
acceptance in events passing FASER Tier E nested in the benchmark region
(Q2 > 4 GeV2, W > 3 GeV, no y window, no cut on x), counted per bin of
z = E_h/nu at 300 fb^-1 on tungsten, for the nu_e (left) and nu_mu (right)
charged currents, and the five generators.

Every count carries the 60% hadron selection efficiency of
faser_pions.HADRON_SELECTION_EFF (user, 2026-09-21); the K/pi ratio and both
charge ratios are untouched by it.

Same layout as pp15, deliberately: charge-summed yields on top, the ratio
to GENIE's default tune in the middle, the K-/K+ charge ratio at the
bottom, the z axis stopping at 0.8.  The one figure this carries that the
pion one cannot is the strangeness of the hadronisation model, which is
what the K/pi ratio and the K-/K+ ratio measure.

IN THE PAPER, NEXT TO THE PION FIGURE (user, 2026-09-21): "pp16 goes to
paper, next to the pion figure."  It sits immediately after
Fig. sidis-yields in bench-sec-rates.tex and its caption is written against
that one.  tools/check_paper_figures_used.py enforces IN_PAPER in both
directions, so the flag, the \\includegraphics and the PNG's entry in
tools/sync_paper_figures.sh move together.

Neutrino current only, the 2026-09-07 carve-out from rule 2b that the whole
single-inclusive study sits under; both GENIE tunes drawn.

Inputs: results_nu/faser_sidis.json (analysis/faser_sidis.py).
Usage: analysis/paper_plots/pp16_sidis_yields_kaon.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SLUG = "pp16_sidis_yields_kaon"
V1_NUMBER = None
# NOT IN THE PAPER since 2026-09-30: merged into pp15's inclusive-hadron
# figure (user: "FASERnu cannot tell apart pions from kaons"); kept for the
# report.
IN_PAPER = False
TITLE = "Charged-kaon yields in z at FASER$\\nu$, 300 fb$^{-1}$"
OUTPUT = "pp_sidis_yields_kaon.png"
RESULTS = "results_nu"

SRC = f"{BASE}/results_nu/faser_sidis.json"
REF = "powheg_nu"
RATIO_REF = "genie_lo"
Z_MAX = 0.8
# nu + nubar since 2026-10-04, as pp15
COLUMNS = [("nu_e", r"$\nu_e+\bar\nu_e$ CC"), ("nu_mu", r"$\nu_\mu+\bar\nu_\mu$ CC")]
GEN_ORDER = [("powheg_nu", "POWHEG-V2", "#1f5fa9", "-"),
             ("herwig_nlo_full", "Herwig", "#3fa66a", "-"),
             ("sherpa_nlo", "Sherpa", "#8e5bd0", "-"),
             ("genie_lo", "GENIE (GRV98LO)", "#c0392b", "-"),
             ("genie", "GENIE (HEDIS)", "#c0392b", "--")]

CAPTION = (
    "Identified charged kaons inside the emulsion track acceptance, counted "
    "per bin of z = E<sub>h</sub>/&nu; in events passing the FASER Tier E "
    "selection nested in the benchmark region (Q<sup>2</sup> &gt; 4 "
    "GeV<sup>2</sup>, W &gt; 3 GeV, no cut on y or x), at 300 "
    "fb<sup>&minus;1</sup> on tungsten, for the &nu;<sub>e</sub> + &nu;&#772;<sub>e</sub> (left) "
    "and &nu;<sub>&mu;</sub> + &nu;&#772;<sub>&mu;</sub> (right) charged currents. The two charges are "
    "summed and a 60% hadron selection efficiency is applied. The shaded "
    "band is the z &lt; 0.1 that is cut. Middle panels: "
    "the ratio to GENIE&rsquo;s default tune. Lower panels: the "
    "K<sup>&minus;</sup>/K<sup>+</sup> charge ratio.")

MESSAGE = """Every yield here is neutrino plus antineutrino charged-current scattering (since 2026-10-04), as on the pion figure. At 300 fb^-1 on the tungsten target the 2813 nu_mu + nubar_mu charged-current events POWHEG-V2 expects to pass Tier E in the region give, at a 60% selection efficiency, 538 identified charged kaons above z = 0.1 inside the emulsion track acceptance -- 288 K+ and 250 K- -- with a statistical uncertainty of 4.4%; for nu_e + nubar_e it is 218 kaons at 7.0%. The generators span 467 to 575 kaons. Integrated, the kaon yield is a fifth of the pion one: 19.7 to 20.2% for the three NLO matchings against 21.3% and 22.8% for GENIE's two tunes, the most direct handle on the strangeness suppression of the fragmentation.

The z shape repeats the pion pattern. Summed over 0.1 < z < 0.8 the NLO matchings give 1.16 (POWHEG-V2), 1.22 (Herwig) and 1.24 (Sherpa) times GENIE's default kaon yield, but they cross it near z = 0.4 and fall to 0.80, 0.71 and 0.72 by z = 0.5: GENIE fragments harder, in kaons as in pions. GENIE's second tune stays within 23% of the default across the drawn range.

The kaon charge ratio is where the models disagree most on this figure. Above z = 0.1 the K+/K- ratio is 1.15 for POWHEG-V2, 1.40 for Herwig, 1.05 for Sherpa and 1.42 and 1.20 for GENIE's tunes -- a 35% spread, against the 5% spread the same generators show on the pion charge ratio. In z, K-/K+ falls from about 1 at the cut to below 0.25 in the last two bins for POWHEG-V2 and the GENIE default, while Sherpa's stays above 1 from z = 0.3 to 0.7. What limits the measurement is the statistics: the intra-event correlation inflates the error by only 1.05, but there are five times fewer kaons than pions -- in bins of 0.05 in z the relative error stays below 10% only up to z = 0.20 for nu_mu, and for nu_e no bin reaches it, the best being 12.9%. A kaon measurement binned in z runs out of events before the models separate, whereas the integrated K/pi ratio does not."""


_D = {}


def _d():
    if not _D:
        with open(SRC) as f:
            _D.update(json.load(f))
    return _D


def _r(fl, key, field):
    return _d()["currents"]["nu"]["flavours"][fl]["generators"][key][
        "regions"][_d()["yield_region"]][field]


def _bins():
    """The drawn window: the first bin above the z cut, and the first bin at
    or above Z_MAX, which is where the figure stops."""
    edges = np.array(_d()["z_edges"])
    lo = int(np.argmax(edges[:-1] >= _d()["z_min"] - 1e-12))
    hi = int(np.argmax(edges[:-1] >= Z_MAX - 1e-12))
    return edges, lo, hi


def _bin_at(z):
    """The index of the bin containing z."""
    edges, _lo, _hi = _bins()
    return int(np.searchsorted(edges, z + 1e-12) - 1)


def _gratio(fl, key, which=None):
    """A generator's kaon spectrum over GENIE's default tune, summed or per
    bin, inside the drawn z window."""
    _edges, lo, hi = _bins()
    y = np.array(_r(fl, key, "K_z"))[lo:hi]
    ref = np.array(_r(fl, RATIO_REF, "K_z"))[lo:hi]
    if which == "sum":
        return y.sum() / ref.sum()
    return y / ref


def _gratio_at(fl, key, z):
    """The same ratio in the bin containing z -- quoted instead of the last
    drawn bin, which holds a handful of kaons and jumps."""
    i = _bin_at(z)
    return _r(fl, key, "K_z")[i] / _r(fl, RATIO_REF, "K_z")[i]


def _chratio(fl, key):
    """K-/K+ per bin inside the drawn z window."""
    _edges, lo, hi = _bins()
    p = np.array(_r(fl, key, "Kp_z"))[lo:hi]
    m = np.array(_r(fl, key, "Km_z"))[lo:hi]
    return m / p


def _kpi(fl, key):
    """Charged kaons as a percentage of the charged pions above the z cut."""
    return 100 * _r(fl, key, "K_total_zcut") / _r(fl, key, "pi_total_zcut")


def _cluster():
    """The intra-event correlation's inflation of the statistical error in
    the first kaon bin above the z cut."""
    _edges, lo, _hi = _bins()
    return _r("nu_mu", REF, "K_z_err")[lo] / _r("nu_mu", REF, "K_z_err_poisson")[lo]


def _reach(fl):
    """The upper edge of the last z bin whose statistical error is under 10%,
    or nan when no bin reaches it -- which is the nu_e case here, and the
    reason this returns nan rather than raising."""
    y = np.array(_r(fl, REF, "K_z"))
    e = np.array(_r(fl, REF, "K_z_err"))
    edges = np.array(_d()["z_edges"])
    ok = [i for i in range(len(y))
          if edges[i] >= _d()["z_min"] - 1e-12 and y[i] > 0 and 100 * e[i] / y[i] < 10.0]
    return float(edges[max(ok) + 1]) if ok else float("nan")


def _best_err(fl):
    """The smallest relative statistical error of any drawn kaon bin."""
    _edges, lo, hi = _bins()
    y = np.array(_r(fl, REF, "K_z"))[lo:hi]
    e = np.array(_r(fl, REF, "K_z_err"))[lo:hi]
    return float(np.min(100 * e[y > 0] / y[y > 0]))


CLAIMS = [
    {"what": "POWHEG-V2, nu_mu + nubar_mu: 2813 Tier E events, 288 K+ and "
             "250 K- selected above z = 0.1, 538 kaons at 4.4%",
     "check": lambda: (abs(_r("nu_mu", REF, "events") - 2813) < 1
                       and abs(_r("nu_mu", REF, "Kp_total_zcut") - 288) < 1
                       and abs(_r("nu_mu", REF, "Km_total_zcut") - 250) < 1
                       and abs(_r("nu_mu", REF, "K_total_zcut") - 538) < 1
                       and abs(100 * _r("nu_mu", REF, "K_total_zcut_err")
                               / _r("nu_mu", REF, "K_total_zcut") - 4.4) < 0.05),
     "detail": lambda: (f"{_r('nu_mu', REF, 'events'):.0f} events, "
                        f"{_r('nu_mu', REF, 'Kp_total_zcut'):.0f} / "
                        f"{_r('nu_mu', REF, 'Km_total_zcut'):.0f}, "
                        f"{100 * _r('nu_mu', REF, 'K_total_zcut_err') / _r('nu_mu', REF, 'K_total_zcut'):.2f}%")},
    {"what": "nu_e + nubar_e: 218 kaons at 7.0%",
     "check": lambda: (abs(_r("nu_e", REF, "K_total_zcut") - 218) < 1
                       and abs(100 * _r("nu_e", REF, "K_total_zcut_err")
                               / _r("nu_e", REF, "K_total_zcut") - 6.95) < 0.05),
     "detail": lambda: (f"{_r('nu_e', REF, 'K_total_zcut'):.0f} at "
                        f"{100 * _r('nu_e', REF, 'K_total_zcut_err') / _r('nu_e', REF, 'K_total_zcut'):.2f}%")},
    {"what": "the generators span 467 to 575 kaons (nu_mu + nubar_mu)",
     "check": lambda: (abs(min(_r("nu_mu", k, "K_total_zcut") for k, _b, _c, _l in GEN_ORDER) - 467) < 1
                       and abs(max(_r("nu_mu", k, "K_total_zcut") for k, _b, _c, _l in GEN_ORDER) - 575) < 1),
     "detail": lambda: ", ".join(f"{k} {_r('nu_mu', k, 'K_total_zcut'):.0f}"
                                 for k, _b, _c, _l in GEN_ORDER)},
    {"what": "kaons are 19.7-20.2% of the pions for the three NLO matchings "
             "and 21.3% and 22.8% for GENIE's two tunes",
     "check": lambda: (all(19.7 - 0.05 <= _kpi("nu_mu", k) <= 20.2 + 0.05
                           for k in ("powheg_nu", "herwig_nlo_full", "sherpa_nlo"))
                       and abs(_kpi("nu_mu", "genie_lo") - 21.3) < 0.05
                       and abs(_kpi("nu_mu", "genie") - 22.8) < 0.05),
     "detail": lambda: ", ".join(f"{k} {_kpi('nu_mu', k):.1f}%"
                                 for k, _b, _c, _l in GEN_ORDER)},
    {"what": "over GENIE's default tune, summed over 0.1 < z < 0.8: "
             "POWHEG-V2 1.16, Herwig 1.22, Sherpa 1.24; they cross it near "
             "z = 0.4 and are at 0.80 / 0.71 / 0.72 by z = 0.5",
     "check": lambda: (abs(_gratio("nu_mu", "powheg_nu", "sum") - 1.16) < 0.01
                       and abs(_gratio("nu_mu", "herwig_nlo_full", "sum") - 1.22) < 0.01
                       and abs(_gratio("nu_mu", "sherpa_nlo", "sum") - 1.24) < 0.01
                       and all(_gratio_at("nu_mu", k, 0.32) > 1.0
                               > _gratio_at("nu_mu", k, 0.47)
                               for k in ("powheg_nu", "herwig_nlo_full", "sherpa_nlo"))
                       and abs(_gratio_at("nu_mu", "powheg_nu", 0.52) - 0.80) < 0.01
                       and abs(_gratio_at("nu_mu", "herwig_nlo_full", 0.52) - 0.71) < 0.01
                       and abs(_gratio_at("nu_mu", "sherpa_nlo", 0.52) - 0.72) < 0.01),
     "detail": lambda: ", ".join(
         f"{k} {_gratio('nu_mu', k, 'sum'):.3f} "
         f"(z=0.3 {_gratio_at('nu_mu', k, 0.32):.2f}, "
         f"z=0.4 {_gratio_at('nu_mu', k, 0.42):.2f}, "
         f"z=0.5 {_gratio_at('nu_mu', k, 0.52):.2f})"
         for k, _b, _c, _l in GEN_ORDER if k != RATIO_REF)},
    {"what": "GENIE HEDIS stays within 23% of the default tune across the "
             "drawn range",
     "check": lambda: float(np.max(np.abs(_gratio("nu_mu", "genie") - 1.0))) < 0.23,
     "detail": lambda: ("max |HEDIS/default - 1| = "
                        f"{float(np.max(np.abs(_gratio('nu_mu', 'genie') - 1.0))):.3f}")},
    {"what": "K+/K- above z = 0.1 is 1.15 (POWHEG-V2), 1.40 (Herwig), 1.05 "
             "(Sherpa), 1.42 and 1.20 (GENIE) -- a 35% spread against 5% on "
             "the pion charge ratio",
     "check": lambda: (abs(_r("nu_mu", "powheg_nu", "Kp_total_zcut")
                           / _r("nu_mu", "powheg_nu", "Km_total_zcut") - 1.15) < 0.01
                       and abs(_r("nu_mu", "herwig_nlo_full", "Kp_total_zcut")
                               / _r("nu_mu", "herwig_nlo_full", "Km_total_zcut") - 1.40) < 0.01
                       and abs(_r("nu_mu", "sherpa_nlo", "Kp_total_zcut")
                               / _r("nu_mu", "sherpa_nlo", "Km_total_zcut") - 1.05) < 0.01
                       and abs(_r("nu_mu", "genie_lo", "Kp_total_zcut")
                               / _r("nu_mu", "genie_lo", "Km_total_zcut") - 1.42) < 0.01
                       and abs(_r("nu_mu", "genie", "Kp_total_zcut")
                               / _r("nu_mu", "genie", "Km_total_zcut") - 1.20) < 0.01
                       and abs(_spread("Kp", "Km") - 35) < 1
                       and abs(_spread("pip", "pim") - 5) < 1),
     "detail": lambda: (", ".join(
         f"{k} {_r('nu_mu', k, 'Kp_total_zcut') / _r('nu_mu', k, 'Km_total_zcut'):.3f}"
         for k, _b, _c, _l in GEN_ORDER)
         + f"; spread K {_spread('Kp', 'Km'):.0f}% vs pi {_spread('pip', 'pim'):.0f}%")},
    {"what": "in z, K-/K+ falls below 0.25 in the last two drawn bins for "
             "POWHEG-V2 and the GENIE default, while Sherpa's stays above 1 "
             "from z = 0.3 to 0.7",
     "check": lambda: (max(_chratio("nu_mu", "powheg_nu")[-2:]) < 0.25
                       and max(_chratio("nu_mu", "genie_lo")[-2:]) < 0.25
                       and _sherpa_above_one()),
     "detail": lambda: ", ".join(
         f"{k} ({_chratio('nu_mu', k)[0]:.2f} to {_chratio('nu_mu', k)[-1]:.2f})"
         for k, _b, _c, _l in GEN_ORDER)},
    {"what": "clustering factor 1.05 on the statistical error at a 60% "
             "selection efficiency; under 10% per bin only up to z = 0.20 "
             "(nu_mu) and in no nu_e bin, the best being 12.9%",
     "check": lambda: (abs(_cluster() - 1.05) < 0.005
                       and abs(_reach("nu_mu") - 0.20) < 1e-9
                       and np.isnan(_reach("nu_e"))
                       and abs(_best_err("nu_e") - 12.9) < 0.05),
     "detail": lambda: (f"{_cluster():.3f}; nu_mu {_reach('nu_mu'):.2f}, "
                        f"nu_e {_reach('nu_e')} (best bin "
                        f"{_best_err('nu_e'):.1f}%)")},
]


def _spread(plus, minus):
    """The percentage spread of a charge ratio across the five generators."""
    r = [_r("nu_mu", k, f"{plus}_total_zcut") / _r("nu_mu", k, f"{minus}_total_zcut")
         for k, _b, _c, _l in GEN_ORDER]
    return 100 * (max(r) / min(r) - 1)


def _sherpa_above_one():
    """Sherpa's K-/K+ exceeds 1 in every bin from z = 0.3 to z = 0.7."""
    edges, lo, hi = _bins()
    r = _chratio("nu_mu", "sherpa_nlo")
    sel = (edges[lo:hi] >= 0.3 - 1e-12) & (edges[lo:hi] <= 0.65 + 1e-12)
    return bool(r[sel].min() > 1.0)


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    d = _d()
    reg = d["yield_region"]
    edges, _lo, _hi = _bins()
    zmin = d["z_min"]
    fig, axes = plt.subplots(3, len(COLUMNS), figsize=(12.6, 9.2), sharex="col",
                             gridspec_kw={"height_ratios": (2.4, 1.0, 1.0),
                                          "hspace": 0.07, "wspace": 0.24})
    for col, (fname, beam) in enumerate(COLUMNS):
        ax, axg, axr = axes[0, col], axes[1, col], axes[2, col]
        gens = d["currents"]["nu"]["flavours"][fname]["generators"]
        ref = np.array(gens[RATIO_REF]["regions"][reg]["K_z"])
        for key, label, colour, ls in GEN_ORDER:
            g = gens.get(key)
            r = g["regions"].get(reg) if g else None
            if not r:
                continue
            y = np.array(r["K_z"])
            ax.step(edges, np.append(y, y[-1]), where="post",
                    color=colour, lw=1.7, ls=ls, label=tex(label))
            with np.errstate(divide="ignore", invalid="ignore"):
                gr = np.where(ref > 0, y / ref, np.nan)
            axg.step(edges, np.append(gr, gr[-1]), where="post",
                     color=colour, lw=1.4, ls=ls)
            kp = np.array(r["Kp_z"])
            km = np.array(r["Km_z"])
            with np.errstate(divide="ignore", invalid="ignore"):
                rat = np.where(kp > 0, km / kp, np.nan)
            axr.step(edges, np.append(rat, rat[-1]), where="post",
                     color=colour, lw=1.4, ls=ls)
        for a in (ax, axg, axr):
            a.axvspan(0.0, zmin, color="#000000", alpha=0.07, lw=0)
            a.set_xlim(0.0, Z_MAX)
        ax.set_yscale("log")
        # The log axis is set from the drawn window, not autoscaled: the
        # histogram runs to z = 1 while the figure stops at Z_MAX, and
        # matplotlib reads the bins past the x limit.
        shown = [np.array(g["regions"][reg]["K_z"])[:_hi]
                 for k, _b, _c, _l in GEN_ORDER
                 if (g := gens.get(k)) is not None]
        above = np.concatenate([y[_lo:] for y in shown])
        ax.set_ylim(10 ** np.floor(np.log10(above[above > 0].min())),
                    10 ** (np.ceil(np.log10(max(y.max() for y in shown))) + 0.4))
        # the hadrons counted are named in the title (user, 2026-09-21)
        ax.set_title(tex(r"%s, FASER$\nu$: $K^+ + K^-$" % beam),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.set_ylabel(tex(r"Events per bin, 300 fb$^{-1}$"),
                      fontsize=plotstyle.FS_YLABEL - 1)
        axg.set_ylabel(tex(r"ratio to GENIE"), fontsize=plotstyle.FS_YLABEL - 2)
        axg.set_ylim(0.0, 1.8)
        axg.axhline(1.0, color="#111111", lw=1.0, ls="--")
        axr.set_ylabel(tex(r"$K^-/K^+$"), fontsize=plotstyle.FS_YLABEL - 2)
        axr.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        axr.set_ylim(0.0, 1.6)
        axr.axhline(1.0, color="#111111", lw=1.0, ls="--")
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(axg, labelbottom=False)
        plotstyle.ticks(axr)
        ax.legend(fontsize=plotstyle.FS_LEGEND - 4, frameon=False,
                  loc="upper right")
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
