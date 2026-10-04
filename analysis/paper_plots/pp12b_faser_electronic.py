#!/usr/bin/env python3
"""Paper plot 12b: our predictions against FASER's ELECTRONIC-detector
muon-neutrino measurements.

Promoted from the report's "Comparison with data -> FASER electronic
detector" sub-tab (analysis/plot_faser_electronic.py) at the user's request,
2026-09-21, to sit beside the emulsion comparison of paper plot 12.  Two
panels on FASER's own axis, -L/E_nu, whose SIGN is the neutrino charge (L is
the lepton number) and whose magnitude is 1/E:

  left    the flux-weighted charged-current cross-section per nucleon in each
          bin, GENIE and Sherpa as a ratio to POWHEG-V2.  FASER's per-bin number
          (Table III of arXiv:2412.03186) is its SIMULATION's, so it is not
          drawn (user, 2026-09-30: "we only want FASER data"); it stays in
          the CLAIMS as a check of the fold.
  right   the number of interactions in the fiducial volume at 186 fb^-1,
          ours against the unfolded data (Fig. 5) of
          CERN-FASER-CONF-2026-005, as a ratio to POWHEG-V2.  FASER's
          simulation (Table 7) is likewise not drawn.

What changed against the report figure:

  * THE POWHEG-V2 LADDER IS THE CURRENT ONE: the ten-energy ladder of
    tools/faser_emulsion_ladder_ubexcess.sh (ubexcess_correct 1, genuine
    2112 neutron beams, 100000 events a point), folded by
    analysis/faser_powheg_rates.py and analysis/faser_electronic.py.
    Against the 2026-09-08 ladder the report figure used, the
    cross-sections move by 0.8% at most.  The seven-point band is imported
    from the old ladder's reweighting as a relative shift, since the 
    Les Houches files carry no scale weights.
  * NO PROSE IN THE FIGURE: the nu / nubar labels on the canvas moved into
    the x-axis label, and the region caveats into CAPTION and MESSAGE.

THE REGION IS FASER'S, NOT THE BENCHMARK'S, as for paper plot 12: the
cross-section is FASER's total charged-current one, so there is no W cut and
no y window.  POWHEG-V2 is generated above Q2 > 4 GeV2 and is low by the rate
below that floor, by construction; GENIE's splines carry the whole
cross-section and are the like-for-like row.  Sherpa is absent: no Sherpa
ladder exists over this energy range in FASER's region.

The two panels carry different amounts of detector.  The yields need the
fiducial cylinder (r = 100 mm) as a share of the flux files' 25 x 25 cm
aperture, which assumes a uniform flux; the flux peaks on axis, and the size
of the miss is measured against FASER's own simulated flux (Table III).

Inputs: results_nu/faser_electronic.json, results_nu/faser_powheg_rates.json,
data/faser_electronic/faser_electronic.json.

Usage: analysis/paper_plots/pp12b_faser_electronic.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SLUG = "pp12b_faser_electronic"
V1_NUMBER = "12b"
IN_PAPER = True
TITLE = "Comparison with the FASER electronic-detector measurements"
OUTPUT = "pp_faser_electronic.png"
RESULTS = "results_nu"

SRC = f"{BASE}/results_nu/faser_electronic.json"
RATES = f"{BASE}/results_nu/faser_powheg_rates.json"
DATA = f"{BASE}/data/faser_electronic/faser_electronic.json"
SURFACE = "#ffffff"

# the colours and line styles of paper plots 8 and 12
POWHEG_C, GENIE_C = "#0072b2", "#d55e00"
# Sherpa NLO (user, 2026-09-30), in the colour and style of the emulsion figure
SHERPA_C, SHERPA_LS = "#009e73", "-."
# Herwig (user, 2026-10-04), in its Sect. 6 colour
HERWIG_C, HERWIG_LS = "#7b3294", ":"

# >>> POWHEG-V2 AND SHERPA EXTENDED BELOW Q2 = 4 GeV2 WITH GENIE (user,
# 2026-10-04: "compare like with like"). <<<  Both are generated above the
# floor; FASER's numbers are total CC.  The *_lowq2 columns of
# faser_electronic.json add GENIE's Q2 < 4 piece (genie_lowq2_sigma there).
CAPTION = (
    "Our predictions against FASER&rsquo;s electronic-detector measurements "
    "of muon-neutrino charged-current scattering on tungsten, in bins of "
    "&minus;L/E<sub>&nu;</sub>, whose sign is the neutrino charge and whose "
    "magnitude is the inverse energy; the central bin, above 1 TeV, holds "
    "both charges. Left: the flux-weighted cross-section per nucleon, "
    "against the value FASER quotes from its simulation in arXiv:2412.03186. "
    "Right: the number of interactions in the fiducial volume at 186 "
    "fb<sup>&minus;1</sup>, against FASER&rsquo;s simulation and its "
    "unfolded data (CERN-FASER-CONF-2026-005). POWHEG-V2, Sherpa and Herwig are "
    "folded over the flux from ladders of beam energies generated above "
    "Q&sup2; = 4 GeV&sup2;, extended below that floor with GENIE&rsquo;s "
    "cross-section there; GENIE is the default tune, with the whole "
    "charged-current cross-section. The band on POWHEG-V2 is its "
    "seven-point scale envelope, and the band in the lower right panel "
    "FASER&rsquo;s total uncertainty on its simulation. The yields assume the "
    "flux uniform across the flux files&rsquo; aperture.")

MESSAGE = """FASER publishes, bin by bin in &minus;L/E<sub>&nu;</sub>, the flux-weighted charged-current cross-section per nucleon, a number that can be compared with a prediction without any model of the detector. POWHEG-V2, Sherpa and Herwig are generated above Q2 = 4 GeV2 while FASER's numbers are the total charged-current ones, so all three are extended below that floor with GENIE's Q2 < 4 GeV2 cross-section (its weighted event fraction times its total), which carries 1-23% of the neutrino and 1-37% of the antineutrino rate between 6.8 TeV and 100 GeV. Compared like with like, POWHEG-V2 agrees with GENIE to within 3% in every bin, Sherpa with POWHEG-V2 to within 0.6% and Herwig to within 0.3%; GENIE in its default tune lands 3.7-8.2% above FASER's own value, and POWHEG-V2 4-10% above. FASER's simulation is GENIE too, so this is a check of our flux fold rather than of the physics.

The fiducial yields carry one approximation: the fiducial cylinder is taken as a uniform share of the flux files' aperture, while the flux peaks on axis, and our folded flux is 0.79-0.93 of FASER's in five bins of six. Summed over the six bins POWHEG-V2 predicts 3074 interactions at 186 fb-1 and GENIE 3031, against 3052 unfolded and 2596 in FASER's simulation. Bin by bin the unfolded data differ from POWHEG-V2 by -53% to +74%, against a POWHEG-V2 scale band of at most 1.4%: at this exposure the measurement is limited by the flux, not by the cross-section."""


_J = {}


def _j(p):
    if p not in _J:
        with open(p) as f:
            _J[p] = json.load(f)
    return _J[p]


def _rows():
    return _j(SRC)["rows"]


def _col(key):
    return np.array([r[key] for r in _rows()], float)


def _prl(key):
    return np.array(_j(DATA)["prl"][key], float)


def _conf(key):
    return np.array(_j(DATA)["conf"][key], float)


def _p_over_f():
    return _col("sigma_1e38_cm2_per_nucleon_lowq2") / _prl("sigma_1e38_cm2_per_nucleon")


def _g_over_f():
    return _col("genie_sigma_1e38_cm2_per_nucleon") / _prl("sigma_1e38_cm2_per_nucleon")


def _flux_ratio():
    return _col("flux_1e6_fb_cm2") / _prl("flux_1e6_fb_cm2")


def _s_over_p():
    return (_col("sherpa_sigma_1e38_cm2_per_nucleon_lowq2")
            / _col("sigma_1e38_cm2_per_nucleon_lowq2"))


def _h_over_p():
    return (_col("herwig_sigma_1e38_cm2_per_nucleon_lowq2")
            / _col("sigma_1e38_cm2_per_nucleon_lowq2"))


def _band():
    n = _col("n_fid_conf_lowq2")
    return 0.5 * (_col("n_fid_conf_lowq2_scale_hi") - _col("n_fid_conf_lowq2_scale_lo")) / n


def private_available():
    """The CONF-note block of the electronic-detector data
    (CERN-FASER-CONF-2026-005) is not part of the public release; claims
    marked private are skipped without it."""
    return "conf" in _j(DATA)


def _data_over_sim():
    return _conf("fig5") / _conf("sim_nominal")


# bin order (data/faser_electronic): nu 100-300, nu 300-600, nu 600-1000,
# both > 1000, nubar 300-1000, nubar 100-300
CLAIMS = [
    {"what": "Sherpa NLO, folded on its ladder extended to 100 GeV and "
             "6.8 TeV and like POWHEG-V2 extended below Q2 = 4 with GENIE, "
             "is within 0.7% of POWHEG-V2 in every bin",
     "check": lambda: all(abs(v - 1.0) < 0.007 for v in _s_over_p()),
     "detail": lambda: " ".join(f"{v:.4f}" for v in _s_over_p())},
    {"what": "Herwig, folded on its own FASER ladder (100 GeV to 6.8 TeV, "
             "positive minus negative weights) and extended below Q2 = 4 with "
             "GENIE in the same way, is within 0.3% of POWHEG-V2 in every bin",
     "check": lambda: all(abs(v - 1.0) < 0.003 for v in _h_over_p()),
     "detail": lambda: " ".join(f"{v:.4f}" for v in _h_over_p())},
    {"what": "the POWHEG-V2 rows are the current ladder (ubexcess_correct 1, "
             "genuine neutrons), not the 2026-09-08 one",
     "check": lambda: (_j(SRC)["ladder"] == os.path.basename(RATES)
                       and _j(RATES)["ladder_dir"] == "ladder-faser-v2"),
     "detail": lambda: f"{_j(SRC)['ladder']} <- {_j(RATES)['ladder_dir']}"},
    {"what": "GENIE lands 3.7-8.2% above FASER's cross-section in every bin",
     "check": lambda: all(1.036 < v < 1.083 for v in _g_over_f()),
     "detail": lambda: " ".join(f"{v:.3f}" for v in _g_over_f())},
    {"what": "extended below Q2 = 4 GeV2 with GENIE, POWHEG-V2 agrees with "
             "GENIE to within 3% in every bin and is 4-10% above FASER's "
             "cross-section",
     "check": lambda: (all(abs(v - 1) < 0.03 for v in _p_over_f() / _g_over_f())
                       and all(1.035 < v < 1.105 for v in _p_over_f())),
     "detail": lambda: ("P/G " + " ".join(f"{v:.3f}" for v in _p_over_f() / _g_over_f())
                        + "; P/F " + " ".join(f"{v:.3f}" for v in _p_over_f()))},
    {"what": "GENIE's Q2 < 4 GeV2 share runs 1-23% (nu) and 1-37% (nubar) "
             "from 6.8 TeV down to 100 GeV",
     "check": lambda: (0.0 < min(_j(SRC)["genie_lowq2_fraction"]["nu_n"]["6800.0"],
                                 _j(SRC)["genie_lowq2_fraction"]["nu_p"]["6800.0"]) < 0.02
                       and 0.22 < _j(SRC)["genie_lowq2_fraction"]["nu_p"]["100.0"] < 0.24
                       and 0.36 < _j(SRC)["genie_lowq2_fraction"]["nubar_n"]["100.0"] < 0.38),
     "detail": lambda: str({k: (round(v["100.0"], 3), round(v["6800.0"], 3))
                            for k, v in _j(SRC)["genie_lowq2_fraction"].items()})},
    {"what": "our folded flux is 0.79-0.93 of FASER's in five bins of six",
     "check": lambda: sum(0.785 < v < 0.935 for v in _flux_ratio()) == 5,
     "detail": lambda: " ".join(f"{v:.3f}" for v in _flux_ratio())},
    {"what": "summed yields at 186 fb-1: POWHEG-V2 3074, GENIE 3031, unfolded "
             "3052, FASER simulation 2596",
     "private": True,
     "check": lambda: (round(_col("n_fid_conf_lowq2").sum()) == 3074
                       and round(_col("genie_n_fid_conf").sum()) == 3031
                       and round(_conf("fig5").sum()) == 3052
                       and round(_conf("sim_nominal").sum()) == 2596),
     "detail": lambda: (f"{_col('n_fid_conf_lowq2').sum():.1f} "
                        f"{_col('genie_n_fid_conf').sum():.1f} "
                        f"{_conf('fig5').sum():.1f} "
                        f"{_conf('sim_nominal').sum():.1f}")},
    {"what": "unfolded data over POWHEG-V2 spans -53% to +74% bin by "
             "bin, while the POWHEG-V2 scale band is at most 1.4%",
     "private": True,
     "check": lambda: (0.465 < (_conf("fig5") / _col("n_fid_conf_lowq2")).min() < 0.48
                       and 1.735 < (_conf("fig5") / _col("n_fid_conf_lowq2")).max() < 1.745
                       and _band().max() < 0.014),
     "detail": lambda: (" ".join(f"{v:.3f}" for v in _conf("fig5") / _col("n_fid_conf_lowq2"))
                        + f"; band max {100 * _band().max():.2f}%")},
]


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt

    dat = _j(DATA)
    edges = np.array(dat["edges_minus_L_over_E"], float)
    mid = 0.5 * (edges[:-1] + edges[1:])
    half = 0.5 * np.diff(edges)
    # the CONF-note block (the unfolded data) is not part of the public
    # release: without it the right panel shows the predictions alone
    prl, conf = dat["prl"], dat.get("conf")

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.6),
                             gridspec_kw={"height_ratios": [1.9, 1.0],
                                          "hspace": 0.0, "wspace": 0.26})
    (a1, a2), (r1, r2) = axes
    fig.subplots_adjust(top=0.80, bottom=0.105, left=0.085, right=0.985)

    # ---------------------------------------------------- the cross-section
    ours_s = _col("sigma_1e38_cm2_per_nucleon_lowq2")
    gen_s = _col("genie_sigma_1e38_cm2_per_nucleon")
    sh_s = _col("sherpa_sigma_1e38_cm2_per_nucleon_lowq2")
    # NO FASER ENTRY IN THIS PANEL (user, 2026-09-30: "we only want FASER
    # data, remove their simulations"; option 2 of the three offered).  The
    # cross-section FASER quotes per bin is its simulation's, not a
    # measurement, so the panel compares the generators alone, as a ratio to
    # POWHEG-V2.  FASER's number still enters the CLAIMS, as a check.
    # POWHEG-V2's scale band, in both panels as on the right.  Within one bin
    # <sigma> and the yield are the same flux-weighted sum over different
    # constants (the flux in the bin; lumi x column x geometry), so the
    # yield's relative band IS <sigma>'s -- exactly, not approximately.
    _n = _col("n_fid_conf_lowq2")
    s_hi = ours_s * _col("n_fid_conf_lowq2_scale_hi") / _n
    s_lo = ours_s * _col("n_fid_conf_lowq2_scale_lo") / _n
    a1.stairs(ours_s, edges, color=POWHEG_C, lw=2.0, baseline=None,
              label=tex("POWHEG-V2"))
    a1.stairs(s_hi, edges, baseline=s_lo, fill=True, color=POWHEG_C,
              alpha=0.22, lw=0)
    a1.stairs(gen_s, edges, color=GENIE_C, lw=1.8, ls="--", baseline=None,
              label=tex("GENIE (GRV98LO)"))
    a1.set_yscale("log")
    # the range from the curves and the band: FASER's points used to set it,
    # and without them the lowest POWHEG-V2 bin sat on the frame
    a1.stairs(sh_s, edges, color=SHERPA_C, lw=1.8, ls=SHERPA_LS,
              baseline=None, label=tex("Sherpa"))
    hw_s = _col("herwig_sigma_1e38_cm2_per_nucleon_lowq2")
    a1.stairs(hw_s, edges, color=HERWIG_C, lw=1.8, ls=HERWIG_LS,
              baseline=None, label=tex("Herwig"))
    _v = np.concatenate([ours_s, gen_s, sh_s, hw_s, s_hi, s_lo])
    a1.set_ylim(_v.min() / 1.35, _v.max() * 1.35)
    a1.set_ylabel(tex(r"$\langle\sigma\rangle$  "
                      r"[$10^{-38}$ cm$^2$/nucleon]"),
                  fontsize=plotstyle.FS_YLABEL - 3)
    a1.set_title(tex("Flux-weighted cross-section"),
                 fontsize=plotstyle.FS_PANEL_TITLE + 3)
    r1.stairs(np.ones_like(ours_s), edges, color=POWHEG_C, lw=2.0,
              baseline=None)
    r1.stairs(s_hi / ours_s, edges, baseline=s_lo / ours_s, fill=True,
              color=POWHEG_C, alpha=0.22, lw=0)
    r1.stairs(gen_s / ours_s, edges, color=GENIE_C, lw=1.8, ls="--",
              baseline=None)
    r1.axhline(1.0, color="#9aa1a9", lw=0.9, zorder=1)
    r1.stairs(sh_s / ours_s, edges, color=SHERPA_C, lw=1.8, ls=SHERPA_LS,
              baseline=None)
    r1.stairs(hw_s / ours_s, edges, color=HERWIG_C, lw=1.8, ls=HERWIG_LS,
              baseline=None)
    rv = np.concatenate([gen_s / ours_s, sh_s / ours_s, hw_s / ours_s, s_hi / ours_s,
                         s_lo / ours_s])
    r1.set_ylim(min(0.95, rv.min()) - 0.04, max(1.05, rv.max()) + 0.04)
    r1.set_ylabel(tex("ratio to") + "\n" + tex("POWHEG-V2"),
                  fontsize=plotstyle.FS_YLABEL - 4)

    # ----------------------------------------------------------- the yields
    n_ours = _col("n_fid_conf_lowq2")
    hi, lo = _col("n_fid_conf_lowq2_scale_hi"), _col("n_fid_conf_lowq2_scale_lo")
    gen_n = _col("genie_n_fid_conf")
    sh_n = _col("sherpa_n_fid_conf_lowq2")
    if conf:
        d5 = np.array(conf["fig5"], float)
        d5e = np.array([e if e else 0.0 for e in conf["fig5_err_hi"]], float)
    # FASER DATA ONLY IN THIS PANEL (user, 2026-09-30: "we only want FASER
    # data, remove their simulations which are not relevant here"): the
    # simulation curve and its uncertainty band are gone, and the ratio panel
    # is now to the data.  The simulation still enters the CLAIMS, as numbers.
    a2.stairs(n_ours, edges, color=POWHEG_C, lw=2.0, baseline=None)
    a2.stairs(hi, edges, baseline=lo, fill=True, color=POWHEG_C,
              alpha=0.22, lw=0)
    a2.stairs(gen_n, edges, color=GENIE_C, lw=1.8, ls="--", baseline=None)
    a2.stairs(sh_n, edges, color=SHERPA_C, lw=1.8, ls=SHERPA_LS,
              baseline=None)
    hw_n = _col("herwig_n_fid_conf_lowq2")
    a2.stairs(hw_n, edges, color=HERWIG_C, lw=1.8, ls=HERWIG_LS,
              baseline=None)
    if conf:
        a2.errorbar(mid, d5, yerr=d5e, xerr=half, fmt="o", ms=5,
                    color="#111111", elinewidth=1.3, capsize=0,
                    label=tex("FASER data"))
    a2.set_ylabel(tex(r"interactions in the fiducial volume"),
                  fontsize=plotstyle.FS_YLABEL - 4)
    a2.set_title(tex(r"Fiducial yield, $\mathcal{L}=186$ fb$^{-1}$"),
                 fontsize=plotstyle.FS_PANEL_TITLE + 3)
    # THE RATIO IS TO POWHEG-V2's CENTRAL VALUE, as in the left panel (user,
    # 2026-09-30: "for consistency with the left panel"); the data are drawn
    # as data / POWHEG-V2, with their error scaled by the same number.
    r2.stairs(np.ones_like(n_ours), edges, color=POWHEG_C, lw=2.0,
              baseline=None)
    r2.stairs(hi / n_ours, edges, baseline=lo / n_ours, fill=True,
              color=POWHEG_C, alpha=0.22, lw=0)
    r2.stairs(gen_n / n_ours, edges, color=GENIE_C, lw=1.8, ls="--",
              baseline=None)
    r2.stairs(sh_n / n_ours, edges, color=SHERPA_C, lw=1.8, ls=SHERPA_LS,
              baseline=None)
    r2.stairs(hw_n / n_ours, edges, color=HERWIG_C, lw=1.8, ls=HERWIG_LS,
              baseline=None)
    if conf:
        r2.errorbar(mid, d5 / n_ours, yerr=d5e / n_ours, xerr=half, fmt="o",
                    ms=4, color="#111111", elinewidth=1.2, capsize=0)
    r2.axhline(1.0, color="#9aa1a9", lw=0.9, zorder=1)
    rv = np.concatenate([gen_n / n_ours, sh_n / n_ours]
                        + ([(d5 + d5e) / n_ours, (d5 - d5e) / n_ours] if conf else []))
    r2.set_ylim(min(0.9, rv.min()) - 0.08, max(1.1, rv.max()) + 0.08)
    r2.set_ylabel(tex("ratio to") + "\n" + tex("POWHEG-V2"),
                  fontsize=plotstyle.FS_YLABEL - 4)

    for ax, rx in ((a1, r1), (a2, r2)):
        for a in (ax, rx):
            a.set_xlim(edges[0], edges[-1])
            a.axvline(0.0, color="#9aa1a9", lw=0.8, ls=":", zorder=1)
        # the charge is the sign of the axis: said in the axis label, since
        # no prose goes on the canvas
        rx.set_xlabel(tex(r"$-L/E_\nu$  [GeV$^{-1}$]   "
                          r"($\nu_\mu$ left, $\bar\nu_\mu$ right)"),
                      fontsize=plotstyle.FS_XLABEL - 3)
        # the seven edges do not fit on a linear axis: the three widely
        # spaced ones are labelled, the four central ones ticked only
        rx.set_xticks(list(edges))
        rx.set_xticklabels([tex(t) for t in
                            [r"$-\frac{1}{100}$", r"$-\frac{1}{300}$",
                             "", "", "", r"$\frac{1}{300}$",
                             r"$\frac{1}{100}$"]], fontsize=9)
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(rx)
    # ONE legend above both columns: the panels are full edge to edge
    h, lab = [], []
    for a in (a1, a2, r2):
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    # THE BAND RIDES ON POWHEG-V2's OWN ENTRY (user, 2026-09-30), as on the
    # emulsion figure: one swatch, the band with the central line through it.
    k = lab.index(tex("POWHEG-V2"))
    h[k], lab[k] = (a1.patches[1], h[k]), tex("POWHEG-V2" + plotstyle.MHOU_SUFFIX)
    # ONE ROW, larger (user, 2026-09-30: "more legible"); five entries since
    # Herwig joined (2026-10-04), so a little tighter to stay on one row
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.535, 0.93),
               ncol=len(h), frameon=True, fontsize=plotstyle.FS_LEGEND + 2,
               handlelength=2.0, columnspacing=1.1)
    fig.suptitle(tex(TITLE), y=0.985, fontsize=plotstyle.FS_SUPTITLE)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
