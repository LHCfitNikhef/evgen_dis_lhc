#!/usr/bin/env python3
"""Appendix plot 3: GENIE beside the dedicated neutrino generators NuWro and
GiBUU, on a free proton and on tungsten, at 200 GeV.

User, 2026-10-01: "a small appendix comparing GENIE with ... the other
tailored neutrino event generators that we implemented ... GIBUU and NuWro?
Keep it short, this is just meant to communicate that these two generators
are also part of the project pipeline."  Then: "compare the three MCs GENIE
vs NuWro and GiBUU both on a proton target and on a tungsten target ... No
need to go above 1 TeV"; "200 GeV for both targets"; "reference is FASER, I
would not show POWHEG or any other MC here, this is really only the three
neutrino generators."

  * THE PROCESS: nu_mu CC DIS at E_nu = 200 GeV, Q2 > 4 GeV2 and W > 3 GeV
    with no y cut (the benchmark region), DIS channel only in every
    generator.  Top row a free proton; bottom row the tungsten NUCLEUS with
    each generator's own nuclear model (GENIE G18_02a: local Fermi gas + hA
    intranuclear cascade; NuWro: local Fermi gas, Pauli blocking, cascade
    on), per nucleon.  The invariants are built with the STRUCK nucleon
    (genie/gtohepmc3.cc, nuwro2hepmc), so x and W are per nucleon.
  * THE REFERENCE is GENIE G18_02a with GRV98, FASER's default tune; both
    GENIE rows are drawn (CONVENTIONS.md rule 1b).
  * >>> GiBUU IS ON THE FREE-PROTON ROW ONLY (user decision "c",
    2026-10-01). <<<  Its tungsten transport does not complete at 200 GeV:
    the perturbative-particle buffer overflows (at 100 GeV with the default
    6000 per ensemble, at 200 GeV even with 30000), and its collision
    criterion reports probabilities above one hundreds of times, i.e. the
    time step is too coarse for these multiplicities.  On tungsten it
    completes at 50 GeV, the top of its documented neutrino range
    (arXiv:1205.1061).  Said in MESSAGE and the caption, never silent.
  * POWHEG-V2 and the analytic reference are deliberately absent (user): the
    figure is about the three neutrino generators; GENIE against NLO is
    pp03/pp04.

Inputs: results_nu/histos_nugen_genie{_lo,}_{p,W}_q4w3_200GeV.json,
        results_nu/histos_nuwro{,_W}_q4w3_200GeV.json,
        results_nu/histos_gibuu_q4w3_200GeV.json
        (BENCH_SELECTION=q4w3 BENCH_ENERGY=200 analysis/analyze_nu.py <key>)

Usage: analysis/paper_plots/ppA3_nu_generators.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

PAPER_SECTION = "appendix"
SLUG = "ppA3_nu_generators"
V1_NUMBER = "A3"
IN_PAPER = True
TITLE = "GENIE, NuWro and GiBUU on a proton and on tungsten"
FIG_TITLE = r"Neutrino event generators, $\nu_\mu$ CC DIS at $E_\nu = 200$ GeV"
OUTPUT = "pp_nu_generators.png"
RESULTS = "results_nu"

SURFACE = "#ffffff"
ENERGY = 200.0
TAG = "_q4w3_200GeV"
REGION = "q4w3"
Q2MIN = 4.0
W_MIN = 3.0

# (label, key stem, colour, linestyle); the first is the reference
GENS = [
    ("GENIE (GRV98LO)", {"p": "nugen_genie_lo_p", "W": "nugen_genie_lo_W"},
     "#cc79a7", "-"),
    ("GENIE (HEDIS)", {"p": "nugen_genie_p", "W": "nugen_genie_W"},
     "#e69f00", "-."),
    ("NuWro", {"p": "nuwro", "W": "nuwro_W"}, "#009e73", "--"),
    ("GiBUU", {"p": "gibuu"}, "#8e5bd0", (0, (1.2, 1.2))),
]
REF = GENS[0]
TARGETS = [("p", "Free proton"), ("W", "Tungsten, per nucleon")]
# what each key's label must say about its target (analyze_nu.py writes it)
TARGET_WORD = {"p": "free proton", "W": "tungsten"}

# (obs, x label, x scale, d sigma label, x window, rebin)
OBS = [
    ("xbj", r"$x_{\rm Bj}$", "log", r"d\sigma/dx_{\rm Bj}", (None, 0.7), 1),
    ("y", r"$y$", "linear", r"d\sigma/dy", (0.0, 1.0), 4),
    ("nch", r"$N_{\rm ch}$", "linear", r"d\sigma/dN_{\rm ch}", (-0.5, 15.5), 1),
    ("nprot", r"$N_p$", "linear", r"d\sigma/dN_p", (-0.5, 5.5), 1),
]

CAPTION = (
    "Neutrino charged-current DIS at 200 GeV for Q&sup2; &gt; 4 GeV&sup2; "
    "and W &gt; 3 GeV with no y cut, on a free proton (top) and on tungsten "
    "per nucleon (bottom): Bjorken <i>x</i>, inelasticity, charged-hadron "
    "multiplicity and proton multiplicity from GENIE (default G18_02a tune "
    "with GRV98, and HEDIS), NuWro and GiBUU, with the DIS channel alone. "
    "On tungsten each generator applies its own nuclear model, and the "
    "kinematic variables refer to the struck nucleon. GiBUU is shown on the "
    "proton only: its tungsten transport does not complete at this energy. "
    "Beneath each, the ratio to GENIE G18_02a; error bars are statistical.")


def key(stem):
    return f"{stem}{TAG}"


def _load(stem, t):
    """The result JSON, CHECKED to be the 200 GeV q4w3 result on target t."""
    p = f"{BASE}/{RESULTS}/histos_{key(stem)}.json"
    with open(p) as f:
        d = json.load(f)
    sel = d.get("selection") or {}
    if (sel.get("name") != REGION or abs(sel.get("q2_min", -1) - Q2MIN) > 1e-9
            or abs(sel.get("w2_min", -1) - W_MIN ** 2) > 1e-9
            or sel.get("y_min") != 0.0 or sel.get("y_max") != 1.0):
        raise SystemExit(f"{p}: not the {REGION} region ({sel})")
    if abs(float(d.get("energy_gev", -1)) - ENERGY) > 1e-6:
        raise SystemExit(f"{p}: energy {d.get('energy_gev')}, not {ENERGY:g}")
    if d.get("generator") != stem:
        raise SystemExit(f"{p}: generator {d.get('generator')!r} != {stem!r}")
    if TARGET_WORD[t] not in d.get("label", "").lower():
        raise SystemExit(f"{p}: label {d.get('label')!r} does not say "
                         f"{TARGET_WORD[t]}")
    return d


def stem_of(gen_index, t):
    return GENS[gen_index][1].get(t)


def sigma(stem, t):
    return float(_load(stem, t)["sigma_fid_pb"])


def rel(stem, t):
    """sigma_fid relative to GENIE G18_02a on the same target, minus one."""
    return sigma(stem, t) / sigma(REF[1][t], t) - 1.0


def mean(stem, t, obs):
    d = _load(stem, t)
    if obs in d.get("means", {}):
        return float(d["means"][obs])
    h = d["hists"][obs]
    e = np.asarray(h["edges"])
    y = np.asarray(h["dsig"]) * np.diff(e)
    return float((y * 0.5 * (e[1:] + e[:-1])).sum() / y.sum())


def _merge(edges, y, e, n):
    if n <= 1:
        return edges, y, e
    m = (len(y) // n) * n
    w = np.diff(edges)[:m].reshape(-1, n)
    s = (y[:m].reshape(-1, n) * w).sum(axis=1)
    ee = np.sqrt(((e[:m].reshape(-1, n) * w) ** 2).sum(axis=1))
    wide = w.sum(axis=1)
    return edges[:m + 1:n], s / wide, ee / wide


def hist(stem, t, obs, rebin=1):
    h = _load(stem, t)["hists"][obs]
    return _merge(np.asarray(h["edges"], float), np.asarray(h["dsig"], float),
                  np.asarray(h["err"], float), rebin)


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    fig = plt.figure(figsize=(16.0, 9.6))
    outer = fig.add_gridspec(2, 4, hspace=0.36, wspace=0.30,
                             top=0.86, bottom=0.07, left=0.06, right=0.99)
    xmin_kin = Q2MIN / (2.0 * ENERGY * 0.938272)
    for r, (t, ttitle) in enumerate(TARGETS):
        for c, (obs, xlab, xs, ylab, (x0, x1), rb) in enumerate(OBS):
            inner = outer[r, c].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                            hspace=0.0)
            ax, axr = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            de, dy, _ = hist(REF[1][t], t, obs, rb)
            lo = x0 if x0 is not None else xmin_kin
            lo = float(de[de >= lo - 1e-12][0])
            share = dy * np.diff(de) / float((dy * np.diff(de)).sum())
            vis = ((de[:-1] >= lo - 1e-12) & (de[1:] <= x1 + 1e-12)
                   & (dy > 0) & (share > 5e-3))
            dev = 0.0
            for lab, stems, colour, ls in GENS:
                if t not in stems:
                    continue
                e, y, err = hist(stems[t], t, obs, rb)
                mid = 0.5 * (e[:-1] + e[1:])
                ax.stairs(y, e, color=colour, lw=1.8, ls=ls, baseline=None,
                          zorder=5, label=tex(lab))
                ax.errorbar(mid, y, yerr=err, fmt="none", ecolor=colour,
                            elinewidth=1.0, alpha=0.75, zorder=5)
                with np.errstate(divide="ignore", invalid="ignore"):
                    rr = np.where(dy > 0, y / dy, np.nan)
                    re = np.where(dy > 0, err / dy, np.nan)
                axr.stairs(rr, e, color=colour, lw=1.8, ls=ls, baseline=None,
                           zorder=5)
                if stems is not REF[1]:
                    axr.errorbar(mid, rr, yerr=re, fmt="none", ecolor=colour,
                                 elinewidth=1.0, alpha=0.75, zorder=5)
                v = rr[vis & np.isfinite(rr)]
                if v.size:
                    dev = max(dev, float(np.max(np.abs(v - 1.0))))
            axr.axhline(1.0, color="#9aa1a9", lw=0.9, zorder=1)
            for a in (ax, axr):
                a.set_xscale(xs)
                a.set_xlim(lo, x1)
                a.grid(alpha=0.22, lw=0.6)
                if xs == "log":
                    tk = [v for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
                          if lo <= v <= x1]
                    a.set_xticks(tk)
                    a.set_xticklabels([f"{v:g}" for v in tk])
                    a.xaxis.set_minor_formatter(ticker.NullFormatter())
            if obs == "xbj":
                ax.set_yscale("log")
                pos = dy[vis]
                ax.set_ylim(0.2 * float(pos.min()), 3.0 * float(pos.max()))
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            else:
                ax.set_ylim(0.0, 1.10 * ax.get_ylim()[1])
                ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5,
                                                              prune="lower"))
            dev = min(max(dev, 0.05), 0.9)
            axr.set_ylim(1.0 - 1.15 * dev, 1.0 + 1.15 * dev)
            axr.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4,
                                                           prune="upper"))
            ax.set_ylabel(tex(rf"${ylab}$  [pb]"),
                          fontsize=plotstyle.FS_YLABEL - 3)
            if c == 0:
                axr.set_ylabel(tex("ratio to") + "\n" + tex("G18_02a"),
                               fontsize=plotstyle.FS_YLABEL - 4)
            axr.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 1)
            ax.set_title(tex(ttitle), fontsize=plotstyle.FS_PANEL_TITLE - 1)
            plotstyle.ticks(ax, labelbottom=False)
            plotstyle.ticks(axr)
    h, lab = fig.axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.525, 0.955),
               ncol=4, frameon=True, handlelength=2.6, columnspacing=1.8,
               fontsize=plotstyle.FS_LEGEND + 3)
    fig.suptitle(tex(FIG_TITLE), y=0.995, fontsize=plotstyle.FS_SUPTITLE + 3)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


MESSAGE = """<b>NuWro and GiBUU run through the same pipeline as every other generator in this benchmark.</b> On a free proton at 200 GeV the three neutrino generators agree on the shapes: relative to GENIE's default G18_02a tune, the fiducial cross-section is 7.4% lower in GENIE HEDIS, 9.0% in NuWro and 9.9% in GiBUU, a deficit almost flat in y, and the mean charged-hadron multiplicity is 5.85 (G18_02a), 5.70 (HEDIS), 6.00 (NuWro) and 5.87 (GiBUU).

<b>On tungsten the nuclear models separate them.</b> The rate per nucleon is 5.9% (HEDIS) and 11.8% (NuWro) below G18_02a. The intranuclear cascade raises the mean charged multiplicity from 5.85 to 12.5 in G18_02a and from 6.00 to 10.6 in NuWro, and knocks out on average 5.2 and 4.0 protons respectively; GENIE HEDIS has no cascade, so its multiplicity (5.23) and proton count (0.65) stay at the free-nucleon level.

<b>GiBUU appears on the free proton only.</b> Its tungsten transport does not complete at 200 GeV: the buffer for the particles it propagates overflows, and its collision criterion reports probabilities above one, i.e. the time step is too coarse for the hadron multiplicities at this energy. On tungsten it completes at 50 GeV, the top of the range its neutrino module is documented for."""


def _flat_y(stem, t):
    """Max deviation of the y ratio from the integrated one, 0.1 < y < 0.9."""
    e, y, _ = hist(stem, t, "y", 4)
    _, yd, _ = hist(REF[1][t], t, "y", 4)
    m = (e[:-1] >= 0.1 - 1e-9) & (e[1:] <= 0.9 + 1e-9)
    return float(np.max(np.abs((y[m] / yd[m]) / (1 + rel(stem, t)) - 1.0)))


_P = [("nugen_genie_p", -7.4), ("nuwro", -9.0), ("gibuu", -9.9)]
_W = [("nugen_genie_W", -5.9), ("nuwro_W", -11.8)]

CLAIMS = [
    {"what": "every input is the q4w3 region at 200 GeV, on the target its "
             "key names, and the generator its key says",
     "check": lambda: all(_load(s, t) is not None for _l, st, _c, _s in GENS
                          for t, s in st.items()),
     "detail": lambda: ", ".join(s for _l, st, _c, _s in GENS
                                 for s in st.values())},
    {"what": "free proton: HEDIS -7.4%, NuWro -9.0%, GiBUU -9.9% against "
             "G18_02a; tungsten: HEDIS -5.9%, NuWro -11.8%",
     "check": lambda: (all(abs(100 * rel(s, "p") - v) < 0.06 for s, v in _P)
                       and all(abs(100 * rel(s, "W") - v) < 0.06
                               for s, v in _W)),
     "detail": lambda: ", ".join(f"{s} {100 * rel(s, t):.2f}%"
                                 for t, l in (("p", _P), ("W", _W))
                                 for s, _v in l)},
    {"what": "on the free proton the deficit is almost flat in y: every "
             "y bin in 0.1-0.9 within 6% of the integrated ratio",
     "check": lambda: all(_flat_y(s, "p") < 0.06 for s, _v in _P),
     "detail": lambda: ", ".join(f"{s} {_flat_y(s, 'p'):.3f}"
                                 for s, _v in _P)},
    {"what": "free-proton mean N_ch 5.85 / 5.70 / 6.00 / 5.87",
     "check": lambda: all(abs(mean(s, "p", "nch") - v) < 0.005 for s, v in (
         ("nugen_genie_lo_p", 5.85), ("nugen_genie_p", 5.70),
         ("nuwro", 6.00), ("gibuu", 5.87))),
     "detail": lambda: ", ".join(f"{mean(g[1]['p'], 'p', 'nch'):.3f}"
                                 for g in GENS)},
    {"what": "tungsten: mean N_ch 12.5 (G18_02a), 10.6 (NuWro), 5.23 "
             "(HEDIS); mean N_p 5.2, 4.0, 0.65",
     "check": lambda: (abs(mean("nugen_genie_lo_W", "W", "nch") - 12.5) < 0.05
                       and abs(mean("nuwro_W", "W", "nch") - 10.6) < 0.05
                       and abs(mean("nugen_genie_W", "W", "nch") - 5.23) < 0.005
                       and abs(mean("nugen_genie_lo_W", "W", "nprot") - 5.2) < 0.05
                       and abs(mean("nuwro_W", "W", "nprot") - 4.0) < 0.05
                       and abs(mean("nugen_genie_W", "W", "nprot") - 0.65) < 0.005),
     "detail": lambda: "; ".join(
         f"{s} nch {mean(s, 'W', 'nch'):.3f} np {mean(s, 'W', 'nprot'):.3f}"
         for s in ("nugen_genie_lo_W", "nuwro_W", "nugen_genie_W"))},
    {"what": "HEDIS has no cascade: its tungsten multiplicities stay within "
             "15% of its free-proton ones, while G18_02a's and NuWro's rise "
             "by more than 70%",
     "check": lambda: (abs(mean("nugen_genie_W", "W", "nch")
                           / mean("nugen_genie_p", "p", "nch") - 1) < 0.15
                       and mean("nugen_genie_lo_W", "W", "nch")
                       / mean("nugen_genie_lo_p", "p", "nch") > 1.7
                       and mean("nuwro_W", "W", "nch")
                       / mean("nuwro", "p", "nch") > 1.7),
     "detail": lambda: ", ".join(f"{a} {mean(a, 'W', 'nch') / mean(b, 'p', 'nch'):.2f}"
                                 for a, b in (("nugen_genie_W", "nugen_genie_p"),
                                              ("nugen_genie_lo_W", "nugen_genie_lo_p"),
                                              ("nuwro_W", "nuwro")))},
]


if __name__ == "__main__":
    main()
