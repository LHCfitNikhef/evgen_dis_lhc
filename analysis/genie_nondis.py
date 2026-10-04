#!/usr/bin/env python3
"""Impact of non-DIS processes: GENIE's default tune split process by process.

THE QUESTION.  The benchmark treats every lepton-nucleon interaction as deep
inelastic scattering, and every generator in it computes only that.  GENIE
does not: its default tune G18_02a sums quasi-elastic, resonance, diffractive
and deep-inelastic channels (and, for the charged current, their charm
counterparts), and it is the generator FASER runs.  So GENIE is the one
place in this benchmark where "how much of the inclusive rate is NOT DIS"
can be asked, and this script asks it as a function of the beam energy on
BOTH currents (CONVENTIONS.md rule 2b).

WHAT IS READ.  The cross-section splines gmkspl writes for the FULL default
event-generator list, process by process:

  * neutrino CC: genie/splines/faser/{14,-14}_{p,n}_cc_e8000.xml, built by
    tools/genie_faser_splines.sh for the FASER rate reproduction, which SUMS
    them.  Here they are kept apart.  The list is CC: QEL, RES, DIS, COH,
    DIS-CHARM, QEL-CHARM, MEC, DFR, QEL-LAMBDA -- of which COH and MEC need a
    nucleus and produce no spline on a free nucleon.
  * muon EM: genie/splines/faser/13_{p,n}_em_e8000.xml, built by
    tools/genie_nondis_splines.sh with the list EM: QEL-EM, RES-EM, DIS-EM,
    MEC-EM (the last again nuclear only).

>>> THE TWO CURRENTS ARE NOT ASKED THE SAME QUESTION, AND IT IS THIS BUILD'S
    DOING. <<<  The neutrino splines are TOTAL cross-sections, every Q2 down
    to zero.  The muon ones are not: this GENIE build raises
    `electromagnetic::kMinQ2Limit` from the upstream 0.02 to 4 GeV2
    (patches/genie-em-q2-floor-and-pythia8-teardown.diff, genie/README.md),
    for every EM channel, because the 1/Q^4 photon propagator makes
    full-phase-space EM generation ~0.5% efficient inside the fiducial
    region.  So the muon-side fractions are fractions of the rate ABOVE the
    benchmark's own Q2 floor, where the elastic and resonance form factors
    have already done most of their suppressing, and NOT of the total EM
    rate -- which is dominated by Q2 -> 0 and is not a DIS quantity at all.
    The figure and the report say so in place.

WHAT "DIS" MEANS TO GENIE.  Its DIS channel is the whole inelastic continuum
above the resonance region: every Q2 (down to the floor above) and, below
its W cut of 1.7 GeV, the non-resonant background that the RES/DIS join
scheme assigns to DIS.  The benchmark's DIS region is Q2 > 4 GeV2.  The
process split here therefore answers "how much of the rate do the EXPLICIT
non-DIS channels carry"; the shallow-inelastic corner INSIDE GENIE's DIS
channel (low Q2, low W) is a different and larger number, and needs the
differential -- a later layer of this study.

WHAT IS WRITTEN.
  results/genie_nondis.json      muon side      (target p, n, tungsten/nucleon)
  results_nu/genie_nondis.json   neutrino side  (nu_mu and nubar_mu, same targets)
  results_nu/cmp_genie_nondis_E.png   the figure, three columns
Each JSON carries sigma per process group on a fine log grid in E, the
fractions at the three benchmark energies, and the fraction of the FASERnu
EVENT RATE the non-DIS channels carry, from the vendored fluxes convoluted on
tungsten (74 p + 110 n; the column density cancels in a fraction).

VALIDITY.  G18_02a is declared valid to 1 TeV; the splines extend it to 8 TeV
through config/hienergy exactly as the FASER rate prediction does (memory:
genie-energy-validity).  The figure shades the extrapolated region.

Usage: analysis/genie_nondis.py
"""
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import beams  # noqa: E402
import target  # noqa: E402
import faser_rates as fr  # noqa: E402

SPLDIR = f"{BASE}/genie/splines/faser"
GEV2_TO_PB = 0.389379e9          # 1 GeV^-2 = 0.389379 mb = 3.89379e8 pb
EMAX = 8000
E_LO, E_HI = 10.0, 5000.0        # the plotted range; the flux reaches 5.7 TeV
E_GRID = np.logspace(np.log10(E_LO), np.log10(E_HI), 240)
# + 100 GeV (user, 2026-09-30): the non-DIS appendix also shows E = 100 GeV,
# where the non-DIS channels could be larger than at the benchmark energies.
E_BENCH = sorted({float(e) for e in beams.ENERGIES} | {100.0})
GVLD_EMAX = 1000.0               # G18_02a's declared validity

# Process groups, in the FIXED order they are drawn and tabulated.  The
# grouping key is GENIE's own `proc:` tag (interaction type) plus its charm
# flag; the algorithm name is kept in the JSON for the record.
GROUPS = ["DIS", "DIS charm", "RES", "QEL", "DFR", "MEC", "COH"]
NONDIS = [g for g in GROUPS if not g.startswith("DIS")]
COLOURS = {"DIS": "#2b6cb0", "DIS charm": "#c2317b", "RES": "#3fa66a",
           "QEL": "#eb6834", "DFR": "#8e5bd0", "MEC": "#6b6b6b",
           "COH": "#a08a00"}
TOTAL_COLOUR = "#1f1f1f"

CURRENTS = {
    # key: (pid, list, spline stem, display)
    "mu":    (13,  "em", "13_{t}_em_e{emax}.xml",  r"$\mu^-$ EM, $Q^2 > 4$ GeV$^2$ only"),
    "nu":    (14,  "cc", "14_{t}_cc_e{emax}.xml",  r"$\nu_\mu$ CC"),
    "nubar": (-14, "cc", "-14_{t}_cc_e{emax}.xml", r"$\bar\nu_\mu$ CC"),
}


def _group(name):
    """Process group of one spline from its GENIE name string."""
    m = re.search(r"proc:([^;]*)", name)
    if not m:
        raise ValueError(f"no proc tag in spline name {name!r}")
    itype = m.group(1).split(",")[-1].strip()
    charm = "charm:" in name
    if itype == "DIS":
        return "DIS charm" if charm else "DIS"
    if itype == "QES":
        return "QEL"           # charm and Lambda quasi-elastic folded in
    if itype in ("RES", "DFR", "MEC", "COH"):
        return itype
    raise ValueError(f"unclassified interaction type {itype!r} in {name!r}")


# A knot-to-knot DROP by more than this factor above 1 TeV is not physics:
# every channel here is flat or rising in E by then.  It is GENIE's
# quasi-elastic Q2 integral losing the narrow low-Q2 peak once Q2max ~ 2ME
# reaches thousands of GeV2 -- the Llewellyn-Smith spline falls from 2.28e-11
# to 1.28e-16 GeV^-2 between its 1779 and 2205 GeV knots, on both nucleons.
# The spline is HELD at its last sound knot above it, and the hold is
# recorded and drawn on the figure: the channel is flat in E to 0.3% from
# 300 GeV to the collapse (2.277e-11 -> 2.283e-11), so a constant is the
# physical continuation, whereas a zero would drop up to 40% of the
# antineutrino non-DIS sum above 1.8 TeV and a cliff would draw as physics.
COLLAPSE_FACTOR = 100.0
COLLAPSE_ABOVE_GEV = 1000.0


def _truncate_collapse(e, x):
    """(e, x, e_cut): the knots up to a numerical collapse, or all of them."""
    for i in range(1, len(e)):
        if e[i] > COLLAPSE_ABOVE_GEV and x[i - 1] > 0 and \
                x[i] < x[i - 1] / COLLAPSE_FACTOR:
            return e[:i], x[:i], float(e[i - 1])
    return e, x, None


def _interp(e, x, grid):
    """Log-log interpolation of a spline onto `grid`, zero below its first
    positive knot, held at the last knot above it -- which, after
    _truncate_collapse, is the last knot BEFORE a numerical collapse."""
    e, x, e_cut = _truncate_collapse(e, x)
    ok = x > 0
    if ok.sum() < 2:
        return np.zeros_like(grid)
    le, lx = np.log(e[ok]), np.log(x[ok])
    out = np.exp(np.interp(np.log(grid), le, lx, left=-np.inf, right=lx[-1]))
    out[grid < e[ok][0]] = 0.0
    return out


def read_spline(fn, grid=E_GRID):
    """{group: sigma [pb] on grid}, plus the per-spline record."""
    if not os.path.exists(fn):
        raise SystemExit(f"missing {fn} -- run tools/genie_faser_splines.sh "
                         f"(CC) or tools/genie_nondis_splines.sh (EM)")
    txt = open(fn).read()
    groups = {g: np.zeros_like(grid) for g in GROUPS}
    record = []
    last_knot = None
    for m in re.finditer(r'<spline name="([^"]+)" nknots="(\d+)">(.*?)</spline>',
                         txt, flags=re.S):
        name, body = m.group(1), m.group(3)
        e = np.array([float(v) for v in re.findall(r"<E>\s*([0-9.eE+-]+)", body)])
        x = np.array([float(v) for v in re.findall(r"<xsec>\s*([0-9.eE+-]+)", body)])
        g = _group(name)
        groups[g] += _interp(e, x, grid) * GEV2_TO_PB
        _e, _x, e_cut = _truncate_collapse(e, x)
        record.append({"name": name, "group": g,
                       "alg": name.split("::")[1].split("/")[0],
                       "nknots": len(e), "max_xsec_gev2": float(x.max()),
                       "collapse_held_above_gev": e_cut})
        if e_cut is not None:
            print(f"  NB {os.path.basename(fn)}: {g} spline "
                  f"({name.split('::')[1].split('/')[0]}) collapses "
                  f"numerically above {e_cut:.0f} GeV; held at its value "
                  f"there ({_x[-1] * GEV2_TO_PB:.3g} pb)")
        last_knot = e[-1] if last_knot is None else max(last_knot, e[-1])
    # a clamped spline (memory: genie-energy-validity) must not draw as if it
    # reached the plotted range
    if last_knot is None or last_knot < E_HI:
        raise SystemExit(f"{fn}: last knot {last_knot} GeV is below the "
                         f"plotted range ({E_HI} GeV) -- the spline was "
                         f"clamped by GVLD-Emax; rebuild with config/hienergy")
    return groups, record


def current_table(cur):
    """Per-target sigma by group on the grid, and the fractions."""
    pid, _lst, stem, _disp = CURRENTS[cur]
    per = {}
    for t in ("p", "n"):
        groups, rec = read_spline(f"{SPLDIR}/{stem.format(t=t, emax=EMAX)}")
        per[t] = groups
        per[f"{t}_record"] = rec
    per["W"] = {g: target.per_nucleon(per["p"][g], per["n"][g]) for g in GROUPS}
    out = {}
    for t in ("p", "n", "W"):
        tot = sum(per[t][g] for g in GROUPS)
        nondis = sum(per[t][g] for g in NONDIS)
        out[t] = {"sigma_pb": {g: per[t][g] for g in GROUPS},
                  "total_pb": tot, "nondis_pb": nondis,
                  "frac": {g: np.where(tot > 0, per[t][g] / np.where(tot > 0, tot, 1), 0)
                           for g in GROUPS},
                  "nondis_frac": np.where(tot > 0, nondis / np.where(tot > 0, tot, 1), 0)}
    out["records"] = {"p": per["p_record"], "n": per["n_record"]}
    return out


def _at(grid_vals, e):
    """Log-log interpolation at one energy; exactly zero for an absent
    channel rather than the 1e-300 floor the log needs."""
    v = float(np.exp(np.interp(np.log(e), np.log(E_GRID),
                               np.log(np.maximum(grid_vals, 1e-300)))))
    return 0.0 if v < 1e-200 else v


def flux_weighted(cur, tab):
    """The non-DIS share of the FASERnu EVENT RATE, on tungsten per nucleon.

    The rate is sum_i phi_i sigma(E_i) x column density, and the column
    density cancels in a ratio, so only the flux shape enters.  The neutrino
    flux is a histogram of counts (data/faser_flux, 150 fb^-1); the muon flux
    is an LHAPDF grid of x f(x) with f = n_T L_T dN/dx, so the integral over
    dx is the integral of x f over d ln x (analysis/faser_rates.muon_flux).
    """
    W = tab["W"]
    if cur == "mu":
        e, xf = fr.muon_flux()
        lo = e >= E_LO
        e, xf = e[lo], xf[lo]
        w = xf
        lne = np.log(e)
        def integral(s):
            return float(np.trapezoid(w * s, lne))
    else:
        pid = CURRENTS[cur][0]
        e, phi = fr.flux(str(pid))
        keep = (e >= E_LO)
        e, w = e[keep], phi[keep]
        def integral(s):
            return float(np.sum(w * s))
    tot = np.array([_at(W["total_pb"], x) for x in e])
    res = {"energy_min_gev": E_LO,
           "flux_source": ("data/faser_muon_flux (25x30 cm)" if cur == "mu"
                           else f"data/faser_flux/FASER_{CURRENTS[cur][0]}.txt")}
    den = integral(tot)
    for g in GROUPS + ["nondis"]:
        s = W["nondis_pb"] if g == "nondis" else W["sigma_pb"][g]
        num = integral(np.array([_at(s, x) for x in e]))
        res[g if g == "nondis" else f"frac_{g}"] = num / den if den > 0 else 0.0
    # the mean beam energy of the interaction rate, for the reader
    res["mean_energy_gev"] = integral(tot * e) / den if den > 0 else None
    return res


def to_json(cur, tab):
    d = {"current": cur, "pid": CURRENTS[cur][0],
         "generator_list": CURRENTS[cur][1].upper(),
         "tune": "G18_02a_00_000",
         "gvld_emax_declared_gev": GVLD_EMAX,
         "spline_emax_gev": EMAX,
         "em_q2_floor_gev2": 4.0 if cur == "mu" else None,
         "note": ("EM splines integrate only Q2 > 4 GeV2: this build raises "
                  "electromagnetic::kMinQ2Limit from 0.02 to 4 for every EM "
                  "channel (patches/genie-em-q2-floor-and-pythia8-teardown.diff)"
                  if cur == "mu" else
                  "CC splines are total cross-sections, every Q2 and every y"),
         "groups": GROUPS, "nondis_groups": NONDIS,
         "energies_gev": [float(x) for x in E_GRID],
         "targets": {}, "benchmark": {}, "records": tab["records"]}
    for t in ("p", "n", "W"):
        T = tab[t]
        d["targets"][t] = {
            "sigma_pb": {g: [float(v) for v in T["sigma_pb"][g]] for g in GROUPS},
            "total_pb": [float(v) for v in T["total_pb"]],
            "nondis_frac": [float(v) for v in T["nondis_frac"]]}
        d["benchmark"][t] = {}
        for e in E_BENCH:
            tot = _at(T["total_pb"], e)
            d["benchmark"][t][f"{e:g}"] = {
                "total_pb": tot,
                "sigma_pb": {g: _at(T["sigma_pb"][g], e) for g in GROUPS},
                "frac": {g: _at(T["sigma_pb"][g], e) / tot for g in GROUPS},
                "nondis_frac": _at(T["nondis_pb"], e) / tot}
    d["flux_weighted_W"] = flux_weighted(cur, tab)
    return d


def figure(tabs, out_png):
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    cols = ["mu", "nu", "nubar"]
    fig = plt.figure(figsize=(15.0, 7.6))
    gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.35, 1.0],
                           hspace=0.0, wspace=0.22)
    for j, cur in enumerate(cols):
        tab = tabs[cur]
        ax = fig.add_subplot(gs[0, j])
        ar = fig.add_subplot(gs[1, j], sharex=ax)
        P = tab["p"]
        ax.plot(E_GRID, P["total_pb"], color=TOTAL_COLOUR, lw=1.8,
                label=tex("total, proton"))
        drawn = []
        for g in GROUPS:
            s = P["sigma_pb"][g]
            if not np.any(s > 0):
                continue
            drawn.append(g)
            ax.plot(E_GRID, s, color=COLOURS[g], lw=1.5, label=tex(g))
            ar.plot(E_GRID, 100.0 * P["frac"][g], color=COLOURS[g], lw=1.5)
        # the non-DIS SUM on all three targets, so the nucleon dependence is
        # visible without a second figure
        ar.plot(E_GRID, 100.0 * P["nondis_frac"], color=TOTAL_COLOUR, lw=2.0,
                label=tex("non-DIS, proton"))
        ar.plot(E_GRID, 100.0 * tab["n"]["nondis_frac"], color=TOTAL_COLOUR,
                lw=1.4, ls="--", label=tex("non-DIS, neutron"))
        ar.plot(E_GRID, 100.0 * tab["W"]["nondis_frac"], color=TOTAL_COLOUR,
                lw=1.4, ls=":", label=tex("non-DIS, tungsten/nucleon"))
        for a in (ax, ar):
            a.set_xscale("log")
            a.set_yscale("log")
            a.set_xlim(E_LO, E_HI)
            a.axvspan(GVLD_EMAX, E_HI, color="#000000", alpha=0.045, lw=0)
            for e in E_BENCH:
                a.axvline(e, color="#888888", lw=0.7, ls="-.", zorder=0)
        ax.set_title(tex(CURRENTS[cur][3]), fontsize=plotstyle.FS_PANEL_TITLE)
        ax.set_ylabel(tex(r"$\sigma$ per nucleon [pb]"))
        ar.set_ylabel(tex("fraction of total [%]"))
        ar.set_xlabel(tex("beam energy [GeV]"))
        ar.set_ylim(1e-3, 150.0)
        ar.axhline(100.0, color="#bbbbbb", lw=0.6)
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(ar)
        ax.legend(loc="lower right", fontsize=8.5, frameon=False, ncol=1)
        if j == 0:
            ar.legend(loc="lower left", fontsize=8.5, frameon=False)
        note = {0: "shaded: beyond the tune's declared validity (1 TeV)",
                1: "QEL held at its 1.8 TeV value above: its spline collapses numerically there",
                2: ""}[j]
        ax.text(0.03, 0.95, tex(note), transform=ax.transAxes, fontsize=8,
                va="top", color="#555555")
    fig.suptitle(tex("GENIE (GRV98LO) process by process, "
                     "free nucleons"), y=0.995)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.09)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {out_png}")


def main():
    tabs = {cur: current_table(cur) for cur in CURRENTS}
    outs = {"mu": f"{BASE}/results/genie_nondis.json",
            "nu": f"{BASE}/results_nu/genie_nondis.json"}
    payload = {"mu": to_json("mu", tabs["mu"])}
    nu = to_json("nu", tabs["nu"])
    nu["antineutrino"] = to_json("nubar", tabs["nubar"])
    payload["nu"] = nu
    for cur, fn in outs.items():
        with open(fn, "w") as f:
            json.dump(payload[cur], f, indent=1)
        print(f"wrote {fn}")
    # the console summary: the fractions at the benchmark energies
    for cur in CURRENTS:
        d = payload["nu"]["antineutrino"] if cur == "nubar" else payload[cur]
        print(f"\n{cur}: {CURRENTS[cur][1].upper()} list, per-process share of "
              f"the total [%] on the proton / neutron / tungsten per nucleon")
        for e in E_BENCH:
            row = []
            for t in ("p", "n", "W"):
                b = d["benchmark"][t][f"{e:g}"]
                row.append(" ".join(f"{g}={100*b['frac'][g]:.3f}" for g in GROUPS
                                    if b["sigma_pb"][g] > 0)
                           + f" | non-DIS={100*b['nondis_frac']:.3f}")
            print(f"  E={e:>5g} GeV  p: {row[0]}\n{'':15s}n: {row[1]}\n{'':15s}W: {row[2]}")
        fw = d["flux_weighted_W"]
        print(f"  flux-weighted (tungsten, E > {fw['energy_min_gev']:g} GeV): "
              f"non-DIS = {100*fw['nondis']:.3f}% of the rate, "
              f"<E> = {fw['mean_energy_gev']:.0f} GeV")
    figure(tabs, f"{BASE}/results_nu/cmp_genie_nondis_E.png")


if __name__ == "__main__":
    main()
