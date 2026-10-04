#!/usr/bin/env python3
"""Appendix plot 1: where GENIE's non-DIS channels sit, on TUNGSTEN.

The redo of the earlier figure for the "Paper plots"
tab.  User, 2026-09-14: "regenerate GENIE's non-DIS sample on p+n (so
corresponding to a tungsten nucleus)".  Everything else is the earlier production's figure.

  * THE TARGET: GENIE G18_02a run separately on a free proton and a free
    neutron (genie/run_genie_nondis.sh <cur> <E> ... <p|n>), each nucleon
    normalised to its own spline cross-sections, combined per nucleon as
    (74 p + 110 n)/184 (analysis/genie_nondis_diff.py).
  * THE REGION: every panel is INSIDE the benchmark region, Q2 > 4 GeV2
    AND W > 3 GeV applied together (user, 2026-09-21: "Both cuts should be
    applied at once.  And also in the x_Bj plot ... I want to quantify
    impact of non-DIS processes for the fiducial cuts of the benchmark.  If
    these are small, so much the better").  The "<name>_region" histograms
    of genie_nondis_diff.py; Q2 and W are booked from their cut
    upwards.  Inside the region only neutrino diffractive scattering
    survives, at 0.004% of the rate, and NO non-DIS event of the muon
    samples falls in it, so the muon lower panels are empty: a measured
    zero, stated in the caption.
  * MUON ROW ON TOP (user, 2026-09-21), as in every other paper figure.
    The lower panels are labelled "non-DIS [%]", without "share".
  * THE MUON ROW IS Q2 > 4 GeV2 ONLY, as in the earlier production: this GENIE build floors every
    electromagnetic channel at Q2 = 4 GeV2.  The neutrino row is the total
    cross-section.  Both are stated in the caption.

BOTH CURRENTS ARE SHOWN, so rule 1b's carve-out is not used here.

Inputs: results{,_nu}/genie_nondis_diff_W[_400GeV].json.

Usage: analysis/paper_plots/ppA1_genie_nondis.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

PAPER_SECTION = "appendix"
SLUG = "ppA1_genie_nondis"
V1_NUMBER = "A1"
# IN THE PAPER in place of the earlier production's figure (user, 2026-09-14: every updated 
# plot goes into the paper).  Read by tools/check_paper_figures_used.py.
IN_PAPER = True
TITLE = "Impact of non-DIS processes in GENIE"
OUTPUT = "pp_genie_nondis.png"
RESULTS = "results_nu"

FILES = {"nu": f"{BASE}/results_nu/genie_nondis_diff_W.json",
         "mu": f"{BASE}/results/genie_nondis_diff_W.json",
         "nu400": f"{BASE}/results_nu/genie_nondis_diff_W_400GeV.json",
         "mu400": f"{BASE}/results/genie_nondis_diff_W_400GeV.json"}

GROUPS = ["DIS", "DIS charm", "RES", "QEL", "DFR"]
NONDIS = ["RES", "QEL", "DFR"]
COLOURS = {"DIS": "#2b6cb0", "DIS charm": "#c2317b", "RES": "#3fa66a",
           "QEL": "#eb6834", "DFR": "#8e5bd0"}
TOTAL = "#1f1f1f"
# DRAWN GROUPS: GENIE's "DIS" and "DIS charm" are one physical channel and are
# drawn as one (user, 2026-09-30: "add together DIS and DIS-charm, which are
# actually the same thing").  The split stays in the JSON and the CLAIMS.
DRAW = [("DIS", ["DIS", "DIS charm"]), ("RES", ["RES"]), ("QEL", ["QEL"]),
        ("DFR", ["DFR"])]
# x_Bj, not y, in the third column; all three inside the region, both cuts
# at once (user, 2026-09-21).
OBS = ["Q2_region", "W_region", "x_region"]
# larger than the house style (user, 2026-09-21: "now it is too small"): six
# panels across a 14-inch canvas are shrunk to the text width in the paper
FS_TICK, FS_LAB, FS_TITLE, FS_LEG = 16, 18, 20, 18

CAPTION = (
    "GENIE G18_02a on tungsten (per nucleon, free protons and neutrons) at "
    "1 TeV, fully inclusive, channel by channel: momentum transfer, hadronic "
    "invariant mass and Bjorken x for muon electromagnetic scattering (top) "
    "and &nu;<sub>&mu;</sub> "
    "charged current (bottom), inside the benchmark region, Q<sup>2</sup> "
    "&gt; 4 GeV<sup>2</sup> and W &gt; 3 GeV, applied to all three; x is "
    "shown above 0.001, the range relevant for FASER. Beneath each distribution, "
    "the percentage of the rate the non-DIS channels carry per bin: only "
    "neutrino diffractive scattering enters the region, and no non-DIS event "
    "of the muon samples does.")

MESSAGE = """On tungsten, as on the proton, the explicit non-DIS channels of GENIE's default tune live where the benchmark does not look. Resonance production sits entirely below W = 2 GeV, at inelasticities of a few per mille, and quasi-elastic scattering at W equal to the mass of the recoiling baryon -- the nucleon, or a charmed baryon near 2.3-2.5 GeV in the neutrino case -- so the W > 3 GeV floor removes both outright. What survives is a sliver of neutrino diffractive scattering: the non-DIS share of the fiducial rate is 0.005% at 400 GeV and 0.004% at 1 TeV for neutrinos, and zero for muons, against 0.3-0.7% of the total cross-section.

What GENIE labels DIS is not all DIS in the benchmark's sense, and the piece that is not lies in Q2, not in W: 3-6% of the neutrino DIS channel sits below Q2 = 4 GeV2, while more than 99.5% of it is above W = 2 GeV at both energies. The benchmark's Q2 floor is what separates the two definitions.

A delivery trap on the way: GENIE's resonance kinematics generator gives up on a growing fraction of events and silently redraws them from the whole list, so the delivered resonance rate is 10-35% below the tune's own cross-section, more at higher energy and on both nucleons. Immaterial here, where the channel does not enter the region, but not for anyone reading GENIE's channel mix at TeV energies."""


def _file(cur, e):
    """The tungsten result at beam energy e (GeV) for one current."""
    import beams
    return (f"{BASE}/{'results_nu' if cur == 'nu' else 'results'}/"
            f"{beams.at_energy('genie_nondis_diff_W', e)}.json")


def _load(k):
    with open(FILES[k]) as f:
        return json.load(f)


def _share_fid(k):
    return _load(k)["regions"]["share_of_fiducial"]["non-DIS"]


def _res_delivered(k):
    d = _load(k)["delivery"]
    return min(d["p"]["res_delivered_over_spline"],
               d["n"]["res_delivered_over_spline"])


CLAIMS = [
    {"what": "the inputs are the tungsten combination, (74 p + 110 n)/184, "
             "in the benchmark region",
     "check": lambda: all(_load(k)["target"] == "W" and _load(k)["w_min"] == 3.0
                          and _load(k)["Z"] == 74 and _load(k)["A"] == 184
                          for k in FILES),
     "detail": lambda: ", ".join(f"{k} {_load(k)['region']}" for k in FILES)},
    {"what": "the non-DIS share of the fiducial rate is below 0.01% on both "
             "currents at both energies (0.004-0.005% neutrino, zero muon)",
     "check": lambda: all(_share_fid(k) < 1e-4 for k in FILES)
     and _share_fid("mu") == 0.0 and _share_fid("mu400") == 0.0,
     "detail": lambda: ", ".join(f"{k} {100*_share_fid(k):.4f}%" for k in FILES)},
    {"what": "the non-DIS share of the TOTAL cross-section is 0.3-0.7% "
             "(so the fiducial region removes it, not the energy)",
     "check": lambda: all(0.0025 < 1 - (_load(k)["sigma_pb"]["DIS"]
                                        + _load(k)["sigma_pb"]["DIS charm"])
                          / _load(k)["sigma_pb"]["total"] < 0.007 for k in FILES),
     "detail": lambda: ", ".join(
         f"{k} {100*(1-(_load(k)['sigma_pb']['DIS']+_load(k)['sigma_pb']['DIS charm'])/_load(k)['sigma_pb']['total']):.3f}%"
         for k in FILES)},
    {"what": "resonances sit entirely below W = 2 GeV, at inelasticities "
             "of a few per mille",
     "check": lambda: all(_load(k)["regions"]["RES"]["w_above_2"] < 1e-3
                          and _load(k)["regions"]["RES"]["mean_y"] < 0.02
                          for k in FILES),
     "detail": lambda: ", ".join(
         f"{k} W>2 {100*_load(k)['regions']['RES']['w_above_2']:.3f}%, <y> "
         f"{_load(k)['regions']['RES']['mean_y']:.4f}" for k in FILES)},
    {"what": "3-6% of the neutrino DIS channel lies below the Q2 floor, "
             "more than 99.5% of it above W = 2 GeV",
     "check": lambda: all(0.025 < _load(k)["dis_channel"]["fraction_below_q2_floor"] < 0.07
                          and _load(k)["dis_channel"]["fraction_below_w2"] < 0.005
                          for k in ("nu", "nu400")),
     "detail": lambda: ", ".join(
         f"{k} Q2<4 {100*_load(k)['dis_channel']['fraction_below_q2_floor']:.2f}%, "
         f"W<2 {100*_load(k)['dis_channel']['fraction_below_w2']:.2f}%"
         for k in ("nu", "nu400"))},
    {"what": "the delivered resonance rate is 10-35% below the spline one on "
             "both nucleons, and the deficit grows with energy on both currents",
     "check": lambda: all(0.6 < _res_delivered(k) < 0.92 for k in FILES)
     and _res_delivered("nu") < _res_delivered("nu400")
     and _res_delivered("mu") < _res_delivered("mu400"),
     "detail": lambda: ", ".join(f"{k} {_res_delivered(k):.3f}" for k in FILES)},
]


def main(energy=1000.0, output=OUTPUT):
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.ticker import FormatStrFormatter, NullFormatter

    data = {}
    for cur in ("nu", "mu"):
        with open(_file(cur, energy)) as f:
            data[cur] = json.load(f)
        if abs(data[cur]["energy_gev"] - energy) > 1e-6:
            raise SystemExit(f"{_file(cur, energy)} is at {data[cur]['energy_gev']} GeV")
    fig = plt.figure(figsize=(14.0, 9.6))
    outer = gridspec.GridSpec(2, 3, figure=fig, hspace=0.34, wspace=0.36,
                              left=0.08, right=0.975, top=0.89, bottom=0.08)
    # the beam named in each row title (user, 2026-09-30)
    num, unit = ((f"{energy / 1000:g}", "TeV") if energy >= 1000
                 else (f"{energy:g}", "GeV"))
    rowlab = {"nu": r"$\nu_\mu$ CC, $E_\nu = %s$ %s" % (num, unit),
              "mu": r"$\mu^-$ NC, $E_\mu = %s$ %s" % (num, unit)}
    handles = {}
    for r, cur in enumerate(("mu", "nu")):
        d = data[cur]
        for c, name in enumerate(OBS):
            inner = gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=outer[r, c], height_ratios=[2.0, 1.0], hspace=0.0)
            ax = fig.add_subplot(inner[0])
            ar = fig.add_subplot(inner[1], sharex=ax)
            h = d["hist"][name]
            ed = np.array(h["edges"])
            tot = sum(np.array(h[g]["dsig"]) for g in GROUPS)
            ax.step(ed, np.r_[tot, tot[-1]], where="post", color=TOTAL, lw=1.8,
                    label=tex("total"))
            drawn = [tot[tot > 0]]
            for lab, members in DRAW:
                v = sum(np.array(h[g]["dsig"]) for g in members)
                if not v.any():
                    continue
                drawn.append(v[v > 0])
                ax.step(ed, np.r_[v, v[-1]], where="post", color=COLOURS[lab],
                        lw=1.4, label=tex(lab))
            nd = sum(np.array(h[g]["dsig"]) for g in NONDIS)
            # fewer than five non-DIS events in a bin is noise, not a share
            nnd = sum(np.array(h[g]["n"]) for g in NONDIS)
            if h["log"]:
                ax.set_xscale("log")
            ax.set_yscale("log")
            if nnd.sum() == 0:
                # >>> A MEASURED ZERO, DRAWN AS ONE (user, 2026-09-30: "are
                # non-DIS effects exactly zero?"). <<<  Yes: GENIE's muon
                # non-DIS channels are RES (W < 2 GeV) and QEL (W = M_N), and
                # not one of their events enters W > 3 GeV.  A log axis cannot
                # show a zero, so this panel is linear, with the share at 0.
                zero = np.where(tot > 0, 0.0, np.nan)
                ar.step(ed, np.r_[zero, zero[-1]], where="post", color=TOTAL, lw=1.8)
                ar.set_ylim(-1.0, 1.0)
                ar.set_yticks([0.0])
                ar.yaxis.set_major_formatter(FormatStrFormatter("%g"))
            else:
                sh = np.where((tot > 0) & (nnd >= 5),
                              100.0 * nd / np.where(tot > 0, tot, 1), np.nan)
                ar.step(ed, np.r_[sh, sh[-1]], where="post", color=TOTAL, lw=1.8)
                shares = [sh]
                # each non-DIS channel that has events, DASHED over the total
                # so that a channel equal to it (DFR here) stays visible
                for g in NONDIS:
                    n_g = np.array(h[g]["n"])
                    if not n_g.any():
                        continue
                    sg = np.where((tot > 0) & (n_g >= 5),
                                  100.0 * np.array(h[g]["dsig"]) / np.where(tot > 0, tot, 1),
                                  np.nan)
                    ar.step(ed, np.r_[sg, sg[-1]], where="post", color=COLOURS[g],
                            lw=1.6, ls=(0, (3.0, 2.0)))
                    shares.append(sg)
                ar.set_yscale("log")
                # THE RANGE HOLDS EVERY DRAWN SHARE (user, 2026-09-30: "display
                # all contributions, even if they are tiny")
                allv = np.concatenate(shares)
                allv = allv[np.isfinite(allv) & (allv > 0)]
                lo_e = np.floor(np.log10(allv.min())) - 0.3
                hi_e = np.ceil(np.log10(allv.max())) + 0.3
                ar.set_ylim(10 ** lo_e, 10 ** hi_e)
                ar.set_yticks([10.0 ** e for e in range(int(np.ceil(lo_e)),
                                                        int(np.floor(hi_e)) + 1, 2)])
            # x > 0.001, the region relevant for FASER (user, 2026-09-21)
            ax.set_xlim(max(ed[0], 1e-3) if name == "x_region" else ed[0], ed[-1])
            # the upper range holds EVERY drawn channel, DFR included
            pos = np.concatenate(drawn)
            ax.set_ylim(pos.min() * 0.3, pos.max() * 3.0)
            ax.set_ylabel(tex(r"d$\sigma$/d" + h["label"].split(" [")[0] + " [pb/unit]"),
                          fontsize=FS_LAB)
            ar.set_ylabel(tex("non-DIS [%]"), fontsize=FS_LAB)
            ar.set_xlabel(tex(h["label"]), fontsize=FS_LAB + 1)
            if c == 1:
                ax.set_title(tex(rowlab[cur]), fontsize=FS_TITLE)
            plotstyle.ticks(ax, labelbottom=False)
            plotstyle.ticks(ar)
            for a in (ax, ar):
                a.tick_params(which="both", labelsize=FS_TICK)
            # the legend is gathered from every panel: the muon row, now on
            # top, has no DIS charm and no DFR
            for hh, ll in zip(*ax.get_legend_handles_labels()):
                handles.setdefault(ll, hh)
            if name == "W_region":
                # a log axis spanning barely a decade labels its minor ticks
                # and they collide; name a few round values instead
                # only the round values inside the range: a tick set outside
                # it stretches the axis (W stops near 14 GeV at 100 GeV)
                ar.set_xticks([v for v in (3, 5, 10, 20, 40) if ed[0] <= v <= ed[-1]])
                ar.xaxis.set_major_formatter(FormatStrFormatter("%g"))
                ar.xaxis.set_minor_formatter(NullFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
    order = [tex(g) for g in ["total"] + [lab for lab, _m in DRAW]]
    ls = [lab for lab in order if lab in handles]
    hs = [handles[lab] for lab in ls]
    fig.legend(hs, ls, loc="upper center", ncol=len(hs), frameon=False,
               fontsize=FS_LEG, bbox_to_anchor=(0.5, 0.985))
    out = f"{BASE}/{RESULTS}/{output}"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
