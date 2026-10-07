#!/usr/bin/env python3
"""Paper plot 15: the expected charged-pion yields in z at FASERnu,
300 fb^-1, in the region.

The yields of the SIDIS study moved to the region (user, 2026-09-18):
charged pions inside the emulsion track acceptance in events passing FASER
Tier E nested in the benchmark region (Q2 > 4 GeV2, W > 3 GeV, no y window,
no cut on x), counted per bin of z = E_h/nu at 300 fb^-1 on tungsten, for
the nu_e (left) and nu_mu (right) charged currents -- the same
cross-section against two fluxes -- and the five generators.

LAYOUT SET BY THE USER, 2026-09-21: nu_e left and nu_mu right; the yields
are CHARGE SUMMED (pi+ + pi-) so the top panel carries one curve per
generator; an intermediate panel holds the ratio to GENIE's default tune,
which is the generator FASER itself runs; the bottom panel keeps the
pi-/pi+ charge ratio; the z axis stops at 0.8, where the yields are spent;
the panel titles name the detector rather than repeating the fiducial cuts
(they are in the caption); and the annotation that used to explain the
solid/dotted charge split is gone with the split itself.

Every count carries the 80% hadron selection efficiency of
faser_pions.HADRON_SELECTION_EFF (user, 2026-09-21); the ratios -- to GENIE
and between the charges -- are untouched by it, and so is every statement
this figure makes about the models.

This is analysis/plot_faser_sidis.py's yields figure built to the
paper-plot contract.  Neutrino current only (user, 2026-09-07); both GENIE
tunes drawn.  The companion figure for identified charged kaons is
pp16_sidis_yields_kaon.py.

Inputs: results_nu/faser_sidis.json (analysis/faser_sidis.py).
Usage: analysis/paper_plots/pp15_sidis_yields.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SLUG = "pp15_sidis_yields"
V1_NUMBER = None
IN_PAPER = True
TITLE = "Inclusive charged-hadron production in z at FASER$\\nu$, 300 fb$^{-1}$"
OUTPUT = "pp_sidis_yields.png"
RESULTS = "results_nu"

SRC = f"{BASE}/results_nu/faser_sidis.json"
REF = "powheg_nu"
# THE RATIO PANEL'S DENOMINATOR IS GENIE'S DEFAULT TUNE (user, 2026-09-21) --
# G18_02a, the generator FASER actually runs, so the panel reads as "what a
# matched NLO calculation would change".  Not POWHEG-V2, which is the
# reference everywhere else in this study.
RATIO_REF = "genie_lo"
# >>> THE MERGED FIGURE (user, 2026-09-30). <<<  "FASERnu cannot tell apart
# pions from kaons", nor the charge: the yield is EVERY STABLE CHARGED HADRON
# ("h", faser_pions.SPECIES), the pi-/pi+ panel is gone, the ratio is to
# POWHEG-V2 with its 7-point band (analysis/mhou_hadron_z.py).  RATIO_REF
# above stays for the CLAIMS that still quote GENIE.
SPECIES_KEY = "h_z"
MHOU = f"{BASE}/results_nu/mhou_hadron_z_q4w3_faser_e_W.json"
Z_MAX = 0.8
# nu_e LEFT, nu_mu RIGHT (user, 2026-09-21).
# nu + nubar since 2026-10-04 (user: "add the nubar contribution")
COLUMNS = [("nu_e", r"$\nu_e+\bar\nu_e$ CC"), ("nu_mu", r"$\nu_\mu+\bar\nu_\mu$ CC")]
# (key, legend label, colour, linestyle).  THE LABELS ARE THE PAPER PLOTS'
# OWN, not the result files' (user, 2026-09-21): "Herwig" and "Sherpa" as in
# pp01/pp02, rather than "Herwig 7 (POWHEG)" and "Sherpa MC@NLO".
# THE MHOU BAND IS GREY, as in every other paper figure (user, 2026-10-06):
# a blue band reads as the statistical error, which is blue elsewhere.
MHOU_GREY, MHOU_ALPHA = "#111111", 0.14
GEN_ORDER = [("powheg_nu", "POWHEG-V2", "#1f5fa9", "-"),
             ("herwig_nlo_full", "Herwig", "#3fa66a", "-"),
             ("sherpa_nlo", "Sherpa", "#8e5bd0", "-"),
             ("genie_lo", "GENIE (GRV98LO)", "#c0392b", "-"),
             ("genie", "GENIE (HEDIS)", "#c0392b", "--")]

CAPTION = (
    "Charged pions inside the emulsion track acceptance, counted per bin of "
    "z = E<sub>h</sub>/&nu; in events passing the FASER Tier E selection "
    "nested in the benchmark region (Q<sup>2</sup> &gt; 4 GeV<sup>2</sup>, "
    "W &gt; 3 GeV, no cut on y or x), at 300 fb<sup>&minus;1</sup> on "
    "tungsten, for the &nu;<sub>e</sub> + &nu;&#772;<sub>e</sub> (left) and "
    "&nu;<sub>&mu;</sub> + &nu;&#772;<sub>&mu;</sub> (right) charged currents. The two charges are summed and an 80% hadron "
    "selection efficiency is applied. The shaded band "
    "is the z &lt; 0.1 that is cut. Middle panels: the ratio to GENIE&rsquo;s "
    "default tune. Lower panels: the &pi;<sup>&minus;</sup>/"
    "&pi;<sup>+</sup> charge ratio.")

MESSAGE = """<b>The figure shows every stable charged hadron</b> (pions, kaons, protons and the charged hyperons), charge summed, since the emulsion identifies neither the species nor the charge. Each yield is the sum of neutrino and antineutrino charged-current scattering, every generator run with its own antineutrino samples against the antineutrino flux of the same flavour. POWHEG-V2 expects 3689 of them above z = 0.1 for nu_mu + nubar_mu at 300 fb^-1 (1.9% statistical, scale band +0.9/-0.7%) and 885 for nu_e + nubar_e, with the EPOS-LHC light + POWHEG charm flux of arXiv:2402.13318; over 0.1 < z < 0.8 Herwig gives 1.03 of that, Sherpa 1.11 and GENIE 0.84 (default) and 0.79 (HEDIS). The paragraphs below describe the charged pions within that sum.

At 300 fb^-1 on the tungsten target POWHEG-V2 expects 2173 nu_mu + nubar_mu charged-current events passing Tier E in the region, 21% of them from the antineutrino, and at an 80% hadron selection efficiency 2811 selected charged pions above z = 0.1 inside the emulsion track acceptance -- 1569 pi+ and 1242 pi- -- with a statistical uncertainty of 2.1%. Charged-current DIS is the same calculation for an electron and for a muon in the final state, so the nu_e prediction differs only through the flux, which is ten times smaller but considerably harder -- a mean energy of 446 GeV against 258 GeV, because much of it is fed by charm and kaon decays rather than by pions -- so Tier E keeps 0.32 of it against 0.26 (neutrinos): 522 events, 26% from the antineutrino, and 675 pions. The other generators span 1766 to 2474 events and 2184 to 3123 pions, the Tier E rate and the multiplicity both entering.

Against GENIE's default tune, the generator FASER runs, the three NLO matchings deliver 1.22 (POWHEG-V2), 1.29 (Herwig) and 1.35 (Sherpa) times as many pions summed over 0.1 < z < 0.8 -- but the ratio is nowhere near flat: it falls from 1.37, 1.46 and 1.66 in the first bin above the cut to 0.47, 0.38 and 0.25 at z = 0.75, a factor of three across the range. GENIE's second tune stays within 8% of the default in every bin, so what the middle panel measures is the hadronisation model and not the DIS model: AGKY fragments harder than Pythia and than either cluster model.

The charge of the pion carries the flavour information: nu p -> mu- pi+ X proceeds through the valence d -> u transition already at Born level while a pi- requires a sea combination, and the antineutrino does the opposite (nubar u -> mu+ d favours pi-), so the antineutrino share dilutes the asymmetry: above z = 0.1 the charge ratio pi+/pi- is 1.26 for POWHEG-V2 and 1.22 to 1.27 across the generators. In z the pi-/pi+ ratio falls from 0.85-0.87 just above the cut to 0.44-0.50 in the last drawn bin for POWHEG-V2, Herwig and both GENIE tunes, while Sherpa stays between 0.70 and 0.86 above z = 0.4 -- the one place where the charge ratio separates the models as strongly as the spectrum does. How far such a measurement reaches is set by the statistics rather than by the acceptance: pions in a bin are not independent counts, several come from the same event, so the bin total has variance N_ev <n^2> and the relative error exceeds the naive 1/sqrt(N). The selection efficiency works the other way: keeping each pion with probability 0.8 decorrelates the pions within an event, and the inflation is only 1.20, so the error grows by less than the naive 1/sqrt(0.8). In bins of 0.05 in z it stays below 10% up to z = 0.40 for nu_mu and z = 0.20 for nu_e, still inside the region where the generators differ by tens of per cent."""


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


def _gratio(fl, key, field="pi_z", which=None):
    """A generator over GENIE's default tune, summed or per bin, inside the
    drawn z window."""
    _edges, lo, hi = _bins()
    y = np.array(_r(fl, key, field))[lo:hi]
    ref = np.array(_r(fl, RATIO_REF, field))[lo:hi]
    if which == "sum":
        return y.sum() / ref.sum()
    r = y / ref
    return r if which is None else (r[0] if which == "first" else r[-1])


def _chratio(fl, key):
    """pi-/pi+ per bin inside the drawn z window."""
    _edges, lo, hi = _bins()
    p = np.array(_r(fl, key, "pip_z"))[lo:hi]
    m = np.array(_r(fl, key, "pim_z"))[lo:hi]
    return m / p


def _flat_sherpa():
    """Sherpa's pi-/pi+ stays between 0.70 and 0.86 above z = 0.4, where
    every other generator is still falling."""
    edges, lo, hi = _bins()
    r = _chratio("nu_mu", "sherpa_nlo")
    sel = edges[lo:hi] >= 0.4 - 1e-12
    return bool(r[sel].min() >= 0.70 - 0.005 and r[sel].max() <= 0.86 + 0.005)


def _cluster():
    """The intra-event correlation's inflation of the statistical error in
    the first bin above the z cut."""
    r = _r("nu_mu", REF, "pi_z_err")
    p = _r("nu_mu", REF, "pi_z_err_poisson")
    _edges, lo, _hi = _bins()
    return r[lo] / p[lo]


def _reach(fl):
    """The upper edge of the last z bin whose statistical error is under 10%."""
    y = np.array(_r(fl, REF, "pi_z"))
    e = np.array(_r(fl, REF, "pi_z_err"))
    edges = np.array(_d()["z_edges"])
    ok = [i for i in range(len(y))
          if edges[i] >= _d()["z_min"] - 1e-12 and y[i] > 0 and 100 * e[i] / y[i] < 10.0]
    return float(edges[max(ok) + 1]) if ok else float("nan")


def _flux_mean(pid):
    """Mean energy of one FASER neutrino flux, GeV."""
    import numpy as np
    import faser_rates as fr
    e, phi, _n = fr.flux_weights("nu", pid)
    e, phi = np.asarray(e), np.asarray(phi)
    return float((e * phi).sum() / phi.sum()), float(phi.sum())


def _tier_eff(fl):
    """Tier E event efficiency inside the region, from
    analysis/faser_sidis_efficiency.py."""
    p = f"{BASE}/results_nu/faser_sidis_efficiency.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    return d["flux_averaged"][fl][REF]["tier_e_event_efficiency"]


def _hsum(fl, key):
    """All charged hadrons over POWHEG-V2, summed over 0.1 < z < 0.8."""
    _edges, lo, hi = _bins()
    y = np.array(_r(fl, key, "h_z"))[lo:hi]
    ref = np.array(_r(fl, REF, "h_z"))[lo:hi]
    return y.sum() / ref.sum()


def _hband():
    with open(MHOU) as f:
        b = json.load(f)
    return 100 * b["zcut_rel_hi"], 100 * b["zcut_rel_lo"]


# Every yield below moved with the flux switch of 2026-10-07: the flux is now
# EPOS-LHC light + POWHEG charm (data/faser_flux_2025), whose nu_e is about
# half the 2021 average's (DPMJET dominated that), so the nu_e yields halved.
CLAIMS = [
    {"what": "all charged hadrons, nu_mu + nubar_mu, POWHEG-V2: 3689 above "
             "z = 0.1 at 1.9%; nu_e + nubar_e 885; scale band +0.9/-0.7%",
     "check": lambda: (abs(_r("nu_mu", REF, "h_total_zcut") - 3689) < 1
                       and abs(100 * _r("nu_mu", REF, "h_total_zcut_err")
                               / _r("nu_mu", REF, "h_total_zcut") - 1.87) < 0.05
                       and abs(_r("nu_e", REF, "h_total_zcut") - 885) < 1
                       and abs(_hband()[0] - 0.9) < 0.05 and abs(_hband()[1] + 0.7) < 0.05),
     "detail": lambda: (f"{_r('nu_mu', REF, 'h_total_zcut'):.0f}, "
                        f"{_r('nu_e', REF, 'h_total_zcut'):.0f}; band "
                        f"+{_hband()[0]:.2f}/{_hband()[1]:.2f}%")},
    {"what": "over POWHEG-V2 in 0.1 < z < 0.8, nu_mu: Herwig 1.03, Sherpa "
             "1.11, GENIE default 0.84, HEDIS 0.79",
     "check": lambda: (abs(_hsum("nu_mu", "herwig_nlo_full") - 1.034) < 0.005
                       and abs(_hsum("nu_mu", "sherpa_nlo") - 1.114) < 0.005
                       and abs(_hsum("nu_mu", "genie_lo") - 0.843) < 0.005
                       and abs(_hsum("nu_mu", "genie") - 0.794) < 0.005),
     "detail": lambda: ", ".join(f"{k} {_hsum('nu_mu', k):.3f}"
                                 for k, _b, _c, _l in GEN_ORDER if k != REF)},
    {"what": "POWHEG-V2, nu_mu + nubar_mu: 2173 Tier E events (21% from "
             "nubar), 1569 pi+ and 1242 pi- selected above z = 0.1, 2811 "
             "pions at 2.1%",
     "check": lambda: (abs(_r("nu_mu", REF, "events") - 2173) < 1
                       and abs(_r("nu_mu", REF, "events_nubar")
                               / _r("nu_mu", REF, "events") - 0.206) < 0.005
                       and abs(_r("nu_mu", REF, "pip_total_zcut") - 1569) < 1
                       and abs(_r("nu_mu", REF, "pim_total_zcut") - 1242) < 1
                       and abs(_r("nu_mu", REF, "pi_total_zcut") - 2811) < 1
                       and abs(100 * _r("nu_mu", REF, "pi_total_zcut_err")
                               / _r("nu_mu", REF, "pi_total_zcut") - 2.10) < 0.05),
     "detail": lambda: (f"{_r('nu_mu', REF, 'events'):.0f} events, "
                        f"{_r('nu_mu', REF, 'pip_total_zcut'):.0f} / "
                        f"{_r('nu_mu', REF, 'pim_total_zcut'):.0f}, "
                        f"{100 * _r('nu_mu', REF, 'pi_total_zcut_err') / _r('nu_mu', REF, 'pi_total_zcut'):.2f}%")},
    {"what": "the nu_e flux is ten times smaller but harder -- 446 GeV "
             "against 258 GeV -- so Tier E keeps 0.32 of it against 0.26",
     "check": lambda: (abs(_flux_mean("12")[0] - 446) < 1
                       and abs(_flux_mean("14")[0] - 258) < 1
                       and abs(_flux_mean("14")[1] / _flux_mean("12")[1] - 9.6) < 0.5
                       and _tier_eff("nu_e") is not None
                       and abs(_tier_eff("nu_e") - 0.315) < 0.005
                       and abs(_tier_eff("nu_mu") - 0.259) < 0.005),
     "detail": lambda: (f"<E> nu_e {_flux_mean('12')[0]:.0f} GeV, nu_mu "
                        f"{_flux_mean('14')[0]:.0f} GeV, ratio of fluxes "
                        f"{_flux_mean('14')[1] / _flux_mean('12')[1]:.1f}; "
                        f"Tier E {_tier_eff('nu_e'):.4f} / "
                        f"{_tier_eff('nu_mu'):.4f}")},
    {"what": "nu_e + nubar_e: 522 events (26% from nubar) and 675 pions",
     "check": lambda: (abs(_r("nu_e", REF, "events") - 522) < 1
                       and abs(_r("nu_e", REF, "events_nubar")
                               / _r("nu_e", REF, "events") - 0.263) < 0.005
                       and abs(_r("nu_e", REF, "pi_total_zcut") - 675) < 1),
     "detail": lambda: f"{_r('nu_e', REF, 'events'):.0f} / {_r('nu_e', REF, 'pi_total_zcut'):.0f}"},
    {"what": "the generators span 1766 to 2474 events and 2184 to 3123 pions",
     "check": lambda: (abs(min(_r("nu_mu", k, "events") for k, _b, _c, _l in GEN_ORDER) - 1766) < 1
                       and abs(max(_r("nu_mu", k, "events") for k, _b, _c, _l in GEN_ORDER) - 2474) < 1
                       and abs(min(_r("nu_mu", k, "pi_total_zcut") for k, _b, _c, _l in GEN_ORDER) - 2184) < 1
                       and abs(max(_r("nu_mu", k, "pi_total_zcut") for k, _b, _c, _l in GEN_ORDER) - 3123) < 1),
     "detail": lambda: ", ".join(f"{k} {_r('nu_mu', k, 'events'):.0f} / "
                                 f"{_r('nu_mu', k, 'pi_total_zcut'):.0f}"
                                 for k, _b, _c, _l in GEN_ORDER)},
    {"what": "over GENIE's default tune, summed over 0.1 < z < 0.8: "
             "POWHEG-V2 1.22, Herwig 1.29, Sherpa 1.35; per bin they fall "
             "from 1.37 / 1.46 / 1.66 to 0.47 / 0.38 / 0.25",
     "check": lambda: (abs(_gratio("nu_mu", "powheg_nu", which="sum") - 1.22) < 0.01
                       and abs(_gratio("nu_mu", "herwig_nlo_full", which="sum") - 1.29) < 0.01
                       and abs(_gratio("nu_mu", "sherpa_nlo", which="sum") - 1.35) < 0.01
                       and abs(_gratio("nu_mu", "powheg_nu", which="first") - 1.37) < 0.01
                       and abs(_gratio("nu_mu", "herwig_nlo_full", which="first") - 1.46) < 0.01
                       and abs(_gratio("nu_mu", "sherpa_nlo", which="first") - 1.66) < 0.01
                       and abs(_gratio("nu_mu", "powheg_nu", which="last") - 0.47) < 0.01
                       and abs(_gratio("nu_mu", "herwig_nlo_full", which="last") - 0.38) < 0.01
                       and abs(_gratio("nu_mu", "sherpa_nlo", which="last") - 0.25) < 0.01),
     "detail": lambda: ", ".join(
         f"{k} {_gratio('nu_mu', k, which='sum'):.3f} "
         f"({_gratio('nu_mu', k, which='first'):.2f} to "
         f"{_gratio('nu_mu', k, which='last'):.2f})"
         for k, _b, _c, _l in GEN_ORDER if k != RATIO_REF)},
    {"what": "GENIE HEDIS stays within 8% of the default tune in every drawn "
             "bin",
     "check": lambda: float(np.max(np.abs(_gratio("nu_mu", "genie") - 1.0))) < 0.08,
     "detail": lambda: ("max |HEDIS/default - 1| = "
                        f"{float(np.max(np.abs(_gratio('nu_mu', 'genie') - 1.0))):.3f}")},
    {"what": "charge ratio pi+/pi- 1.26 for POWHEG-V2, 1.22 to 1.27 across "
             "the generators (nu + nubar); in z pi-/pi+ falls from 0.85-0.87 "
             "just above the cut to 0.44-0.50 in the last drawn bin for every "
             "generator but Sherpa, which stays between 0.70 and 0.86 above "
             "z = 0.4",
     "check": lambda: (abs(_r("nu_mu", REF, "pip_total_zcut")
                           / _r("nu_mu", REF, "pim_total_zcut") - 1.263) < 0.005
                       and all(1.22 - 0.005 <= _r("nu_mu", k, "pip_total_zcut")
                               / _r("nu_mu", k, "pim_total_zcut") <= 1.27 + 0.005
                               for k, _b, _c, _l in GEN_ORDER)
                       and all(0.85 - 0.005 <= _chratio("nu_mu", k)[0] <= 0.87 + 0.005
                               for k, _b, _c, _l in GEN_ORDER)
                       and all(0.44 - 0.005 <= _chratio("nu_mu", k)[-1] <= 0.50 + 0.005
                               for k, _b, _c, _l in GEN_ORDER if k != "sherpa_nlo")
                       and _flat_sherpa()),
     "detail": lambda: ", ".join(
         f"{k} {_r('nu_mu', k, 'pip_total_zcut') / _r('nu_mu', k, 'pim_total_zcut'):.3f}"
         f" ({_chratio('nu_mu', k)[0]:.2f} to {_chratio('nu_mu', k)[-1]:.2f})"
         for k, _b, _c, _l in GEN_ORDER)},
    {"what": "clustering factor 1.20 on the statistical error at an 80% "
             "selection efficiency; under 10% per bin up to z = 0.40 (nu_mu) "
             "and 0.20 (nu_e)",
     "check": lambda: (abs(_cluster() - 1.196) < 0.005
                       and abs(_reach("nu_mu") - 0.40) < 1e-9
                       and abs(_reach("nu_e") - 0.20) < 1e-9),
     "detail": lambda: f"{_cluster():.3f}; {_reach('nu_mu'):.2f} / {_reach('nu_e'):.2f}"},
]


def _band():
    """POWHEG-V2's 7-point scale band on the charged-hadron z spectrum,
    relative, per z bin (analysis/mhou_hadron_z.py, 1 TeV, tungsten)."""
    with open(MHOU) as f:
        b = json.load(f)
    if b["z_edges"] != list(_d()["z_edges"]):
        raise SystemExit("the MHOU band's z bins are not the yields' -- refusing")
    return np.array(b["rel_lo"]), np.array(b["rel_hi"])


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    d = _d()
    reg = d["yield_region"]
    edges, _lo, _hi = _bins()
    zmin = d["z_min"]
    blo, bhi = _band()
    fig, axes = plt.subplots(2, len(COLUMNS), figsize=(12.6, 7.4), sharex="col",
                             gridspec_kw={"height_ratios": (2.4, 1.1),
                                          "hspace": 0.07, "wspace": 0.24})
    for col, (fname, beam) in enumerate(COLUMNS):
        ax, axg = axes[0, col], axes[1, col]
        gens = d["currents"]["nu"]["flavours"][fname]["generators"]
        ref = np.array(gens[REF]["regions"][reg][SPECIES_KEY])
        # POWHEG-V2's scale band, in both panels (user, 2026-09-30)
        ax.stairs(ref * (1 + bhi), edges, baseline=ref * (1 + blo),
                  fill=True, color=MHOU_GREY, alpha=MHOU_ALPHA, lw=0,
                  zorder=1)
        axg.stairs(1 + bhi, edges, baseline=1 + blo, fill=True,
                   color=MHOU_GREY, alpha=MHOU_ALPHA, lw=0, zorder=1)
        for key, label, colour, ls in GEN_ORDER:
            g = gens.get(key)
            r = g["regions"].get(reg) if g else None
            if not r:
                continue
            y = np.array(r[SPECIES_KEY])
            ax.step(edges, np.append(y, y[-1]), where="post",
                    color=colour, lw=1.7, ls=ls, label=tex(label))
            with np.errstate(divide="ignore", invalid="ignore"):
                gr = np.where(ref > 0, y / ref, np.nan)
            axg.step(edges, np.append(gr, gr[-1]), where="post",
                     color=colour, lw=1.4, ls=ls)
        # z < zmin shaded AMBER, not grey: grey is the MHOU band (user,
        # 2026-10-06), and amber is none of the generators' colours
        for a in (ax, axg):
            a.axvspan(0.0, zmin, color="#e69f00", alpha=0.15, lw=0)
            a.set_xlim(0.0, Z_MAX)
        ax.set_yscale("log")
        # THE LOG AXIS IS SET FROM THE DRAWN WINDOW, not autoscaled: the
        # histogram runs to z = 1 while the figure stops at Z_MAX, and
        # matplotlib's autoscale reads the bins past the x limit.
        shown = [np.array(g["regions"][reg][SPECIES_KEY])[:_hi]
                 for k, _b, _c, _l in GEN_ORDER
                 if (g := gens.get(k)) is not None]
        above = np.concatenate([y[_lo:] for y in shown])
        ax.set_ylim(10 ** np.floor(np.log10(above[above > 0].min())),
                    10 ** (np.ceil(np.log10(max(y.max() for y in shown))) + 0.4))
        # panel title fixed by the user, 2026-10-05
        ax.set_title(tex(r"%s SIDIS (charged-hadron production)" % beam),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.set_ylabel(tex(r"Events per bin, 300 fb$^{-1}$"),
                      fontsize=plotstyle.FS_YLABEL - 1)
        axg.set_ylabel(tex("ratio to") + "\n" + tex("POWHEG-V2"),
                       fontsize=plotstyle.FS_YLABEL - 3)
        axg.axhline(1.0, color="#111111", lw=1.0, ls="--")
        # POWHEG-V2's expected statistical error, sqrt(N) on the selected
        # count of each bin, as error bars on its ratio line (Felix Kling via
        # the user, 2026-10-06: "only add them to the POWHEG calculation")
        mid = 0.5 * (edges[1:] + edges[:-1])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(ref > 0, 1.0 / np.sqrt(ref), np.nan)
        win = mid < Z_MAX
        stat = axg.errorbar(mid[win], np.ones(win.sum()), yerr=rel[win],
                            fmt="none", ecolor="#111111",
                            elinewidth=1.3, capsize=2.5, zorder=5)
        vis = [np.array(g["regions"][reg][SPECIES_KEY])[_lo:_hi] / ref[_lo:_hi]
               for k, _b, _c, _l in GEN_ORDER if (g := gens.get(k)) is not None]
        vis += [1.0 - rel[_lo:_hi], 1.0 + rel[_lo:_hi]]
        vis = np.concatenate(vis)
        vis = vis[np.isfinite(vis)]
        axg.set_ylim(max(0.0, vis.min() - 0.1), vis.max() + 0.1)
        axg.set_xlabel(tex(r"$z = E_{h^\pm}/\nu$"), fontsize=plotstyle.FS_XLABEL)
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(axg)
        h_, l_ = ax.get_legend_handles_labels()
        # the band rides on POWHEG-V2's own entry, as on the other figures
        band_patch = ax.patches[0] if ax.patches else None
        if band_patch is not None and l_ and l_[0] == tex(GEN_ORDER[0][1]):
            h_[0], l_[0] = (band_patch, h_[0]), tex(GEN_ORDER[0][1] + plotstyle.MHOU_SUFFIX)
        h_.append(stat)
        l_.append(tex("Stat. errors"))  # black, user 2026-10-06
        ax.legend(h_, l_, fontsize=plotstyle.FS_LEGEND + 1, frameon=False,
                  loc="upper right")
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
