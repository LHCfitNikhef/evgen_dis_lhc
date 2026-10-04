#!/usr/bin/env python3
"""Charged-pion yields in z at FASERnu, against the NNLO SIDIS calculation of
arXiv:2504.05376 (Bonino, Gehrmann, Loechner, Schoenwald, Stagnitto).

Three figures, all from results_nu/faser_sidis.json:

  cmp_faser_sidis_yields  the EVENT YIELDS: charged pions per bin of z at
                          300 fb^-1 in events passing FASER Tier E with
                          W > 3 GeV and NO y window and NO x cut, one column
                          per neutrino flavour.  pi+ solid, pi- dotted;
                          z < 0.1 shaded because it is cut.  Lower panels:
                          pi-/pi+, the charge asymmetry that carries the
                          flavour message.

  cmp_faser_sidis_stat    the same yields for the reference generator with
                          the STATISTICAL ERROR a measurement would have,
                          and the relative error bin by bin -- which is what
                          says where in z the measurement runs out.

  cmp_faser_sidis_mult    the MULTIPLICITY dM/dz per DIS event in the paper's
                          own region, at THEIR beam energy, under their
                          published LO, NLO and NNLO curves.  Right panel:
                          what their x > 0.1 costs a multiplicity.

>>> NEUTRINO ONLY, BY USER DECISION (2026-09-07). <<<  "to simplify the
discussion (I know is a rule violation, but accepted) show only neutrino DIS
here, remove the muon DIS in this specific study".  The muon-side numbers are
still COMPUTED and still in results_nu/faser_sidis.json; they are simply not
drawn.  This is the one documented exception to CONVENTIONS.md rule 2b.

>>> THE VERTICAL AXIS IS A COUNT IN THE BIN, NOT A DENSITY <<< (user,
2026-09-07: "replace 'pions per unit z' by 'Events per bin' ... this is
easier to interpret").  The one exception is the multiplicity figure, which
is compared against a published dM/dz and has to stay a density.

Usage: analysis/plot_faser_sidis.py
Writes results_nu/cmp_faser_sidis_{yields,stat,mult}.png
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plotstyle  # noqa: E402

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

tex = plotstyle.tex
SURFACE = "white"
RES = f"{BASE}/results_nu"

COLOURS = {"powheg_nu": "#1f5fa9", "powheg": "#1f5fa9",
           "herwig_nlo_full": "#3fa66a", "herwig_nlo_powheg_full": "#3fa66a",
           "sherpa_nlo": "#8e5bd0", "sherpa": "#8e5bd0",
           "genie_lo": "#c0392b", "genie": "#c0392b",
           "genie_nnpdf": "#eb6834"}
BEAM = {"nu_mu": r"$\nu_\mu$ CC", "nu_e": r"$\nu_e$ CC", "mu": r"$\mu^-$ NC"}
BEAM_COLOUR = {"nu_mu": "#1f5fa9", "nu_e": "#3fa66a", "mu": "#c0392b"}
# NEUTRINO ONLY (user, 2026-09-07) -- see the module docstring.
COLUMNS = [("nu", "nu_mu"), ("nu", "nu_e")]


def load():
    p = f"{RES}/faser_sidis.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run analysis/faser_sidis.py")
    with open(p) as f:
        return json.load(f)


def dashed(key, cur):
    """GENIE HEDIS is the non-default tune, as everywhere else on the tab."""
    return key == "genie" and cur == "nu"


# ----------------------------------------------------------------- yields
def yields_figure(d):
    reg = d["yield_region"]
    edges = np.array(d["z_edges"])
    zmin = d["z_min"]
    fig, axes = plt.subplots(2, len(COLUMNS), figsize=(12.6, 6.8), sharex="col",
                             gridspec_kw={"height_ratios": (2.4, 1.0),
                                          "hspace": 0.06, "wspace": 0.24})
    for col, (cur, fname) in enumerate(COLUMNS):
        ax, axr = axes[0, col], axes[1, col]
        fl = d["currents"][cur]["flavours"][fname]
        for key, g in fl["generators"].items():
            r = g["regions"].get(reg)
            if not r:
                continue
            c = COLOURS.get(key, "#444444")
            base = "--" if dashed(key, cur) else "-"
            for sp, ls, lw in (("pip", base, 1.7), ("pim", ":", 1.2)):
                y = np.array(r[f"{sp}_z"])
                ax.step(edges, np.append(y, y[-1]), where="post",
                        color=c, lw=lw, ls=ls,
                        label=tex(g["label"]) if sp == "pip" else None)
            pp = np.array(r["pip_z"])
            pm = np.array(r["pim_z"])
            with np.errstate(divide="ignore", invalid="ignore"):
                rat = np.where(pp > 0, pm / pp, np.nan)
            axr.step(edges, np.append(rat, rat[-1]), where="post",
                     color=c, lw=1.4, ls=base)
        for a in (ax, axr):
            a.axvspan(0.0, zmin, color="#000000", alpha=0.07, lw=0)
            a.set_xlim(0.0, 1.0)
        ax.set_yscale("log")
        ax.set_title(tex(r"%s, FASER Tier E $+\ W > 3$ GeV" % BEAM[fname]),
                     fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
        ax.set_ylabel(tex(r"Events per bin, 300 fb$^{-1}$"),
                      fontsize=plotstyle.FS_YLABEL - 1)
        axr.set_ylabel(tex(r"$\pi^-/\pi^+$"), fontsize=plotstyle.FS_YLABEL - 2)
        axr.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        axr.set_ylim(0.0, 1.2)
        axr.axhline(1.0, color="#111111", lw=1.0, ls="--")
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(axr)
        ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False,
                  loc="upper right")
    axes[0, 0].text(0.03, 0.06,
                    tex(r"solid $\pi^+$, dotted $\pi^-$;"
                        r" shaded: $z < 0.1$, cut."
                        "\n"
                        r"No $y$ window and no cut on $x$."),
                    transform=axes[0, 0].transAxes,
                    fontsize=plotstyle.FS_LEGEND - 3, color="#333333")
    out = f"{RES}/cmp_faser_sidis_yields.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# ------------------------------------------------------- statistical errors
def stat_figure(d):
    """The reference generator's yields with the error a measurement would
    have, and that error as a fraction, bin by bin."""
    reg = d["yield_region"]
    edges = np.array(d["z_edges"])
    zc = d["z_min"]
    keep = edges[:-1] >= zc - 1e-12
    xc = 0.5 * (edges[:-1] + edges[1:])
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(13.4, 5.6),
                                  gridspec_kw={"wspace": 0.22})
    for cur, fname in COLUMNS:
        key = d["reference"][cur]
        g = d["currents"][cur]["flavours"][fname]["generators"].get(key)
        if not g or reg not in g["regions"]:
            continue
        r = g["regions"][reg]
        y = np.array(r["pi_z"])
        e = np.array(r["pi_z_err"])
        ep = np.array(r["pi_z_err_poisson"])
        c = BEAM_COLOUR[fname]
        lab = f"{BEAM[fname]}, {g['label']}"
        ax.errorbar(xc[keep], y[keep], yerr=e[keep], fmt="o", ms=4.0, lw=1.3,
                    capsize=2.5, color=c, label=tex(lab))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(y > 0, 100.0 * e / y, np.nan)
            relp = np.where(y > 0, 100.0 * ep / y, np.nan)
        axr.step(edges, np.append(rel, rel[-1]), where="post", color=c, lw=1.8,
                 label=tex(lab))
        axr.step(edges, np.append(relp, relp[-1]), where="post", color=c,
                 lw=1.0, ls=":")
    for a in (ax, axr):
        a.set_xlim(zc, 1.0)
        a.set_yscale("log")
        a.grid(True, which="major", axis="both", alpha=0.25, lw=0.6)
        a.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
        plotstyle.ticks(a)
    ax.set_ylabel(tex(r"Events per bin, 300 fb$^{-1}$"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"expected yield and statistical error, "
                     r"Tier E $+\ W>3$ GeV"),
                 fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    ax.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=False, loc="upper right")
    axr.axhline(10.0, color="#888888", lw=1.0, ls="--")
    axr.set_ylim(1.0, 300.0)
    # NB the per-cent sign goes in UNESCAPED: plotstyle.tex() escapes it, and
    # a pre-escaped one comes out as "\\%", which LaTeX reads as a line break
    # followed by a comment -- it ate the second half of this label once.
    axr.set_ylabel(tex("statistical error [%]"), fontsize=plotstyle.FS_YLABEL)
    axr.set_title(tex(r"relative statistical error per bin"),
                  fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    axr.legend(fontsize=plotstyle.FS_LEGEND - 2, frameon=False, loc="upper left")
    axr.text(0.03, 0.72,
             tex("solid: pions from one event are correlated,\n"
                 r"$\sigma = \sqrt{N_{\rm ev}\langle n^2\rangle}$."
                 "\n"
                 r"dotted: the naive $\sqrt{N_\pi}$."
                 "\n"
                 r"dashed line: $10\%$."),
             transform=axr.transAxes, fontsize=plotstyle.FS_LEGEND - 3,
             color="#333333")
    out = f"{RES}/cmp_faser_sidis_stat.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# ------------------------------------------------------------ multiplicity
def _ref_step(ax, ref, colour, label, band=True):
    edges = np.append(ref["z_lo"], ref["z_hi"][-1])
    y = np.array(ref["dmdz"])
    ax.step(edges, np.append(y, y[-1]), where="post", color=colour, lw=1.9,
            label=tex(label), zorder=5)
    if band and "dmdz_min" in ref:
        ax.fill_between(edges, np.append(ref["dmdz_min"], ref["dmdz_min"][-1]),
                        np.append(ref["dmdz_max"], ref["dmdz_max"][-1]),
                        step="post", color=colour, alpha=0.20, lw=0, zorder=1)


def mult_figure(d):
    refs = d.get("reference_curves") or {}
    cmp_ = d.get("comparison") or {}
    if not cmp_:
        print("no comparison block; skipping the multiplicity figure")
        return
    e_lab = cmp_["energy_gev"]
    zc = d["z_min"]
    redges = np.append(cmp_["z_lo"], cmp_["z_hi"][-1])
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(13.8, 5.8),
                                  gridspec_kw={"wspace": 0.22})

    for lab, colour in (("LO", "#999999"), ("NLO", "#d99a2b"),
                        ("NNLO", "#111111")):
        if lab in refs:
            _ref_step(ax, refs[lab], colour, "Bonino et al. " + lab)
    for key, g in cmp_["generators"].items():
        y = np.array(g["dmdz"])
        ax.step(redges, np.append(y, y[-1]), where="post",
                color=COLOURS.get(key, "#444444"), lw=1.5,
                ls="--" if dashed(key, "nu") else "-",
                label=tex(g["label"]))
    ax.set_yscale("log")
    ax.set_xlim(0.05, 1.0)
    ax.set_ylim(5e-3, 40.0)
    ax.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
    ax.set_ylabel(tex(r"${\rm d}M^{\pi^+}/{\rm d}z$ per DIS event"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"$\nu_\mu p \to \mu^- \pi^+ X$: $x > 0.1$, $W > 3$ GeV, "
                     r"$E_\nu = %g$ GeV" % e_lab),
                 fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False, loc="upper right")
    plotstyle.ticks(ax)

    # THE INTEGRATED MULTIPLICITY, in the panel (user, 2026-09-07): the
    # headline of the comparison is a single number per generator, and it
    # belongs beside the curves rather than only in a table.
    lines = [r"$\int_{%g}^{1}\!{\rm d}z\ {\rm d}M^{\pi^+}/{\rm d}z$"
             r" per DIS event:" % zc,
             r"  Bonino et al.\ NNLO   $%.3f$" % cmp_["reference_total_zcut"]]
    for key, g in cmp_["generators"].items():
        lines.append(r"  %s   $%.3f$  $(%.2f)$"
                     % (g["label"], g["total_zcut"], g["ratio_total_zcut"]))
    ax.text(0.02, 0.03, tex("\n".join(lines)), transform=ax.transAxes,
            fontsize=plotstyle.FS_LEGEND - 4, color="#222222", va="bottom",
            bbox=dict(boxstyle="round,pad=0.35", fc="#f7f7f7", ec="#cccccc",
                      lw=0.8))

    nnlo = refs.get("NNLO", {})
    if "dmdz_min" in nnlo:
        c = np.array(nnlo["dmdz"])
        axr.fill_between(redges, np.append(np.array(nnlo["dmdz_min"]) / c, 1.0),
                         np.append(np.array(nnlo["dmdz_max"]) / c, 1.0),
                         step="post", color="#111111", alpha=0.18, lw=0)
    axr.axhline(1.0, color="#111111", lw=1.2, ls="--")
    for key, g in cmp_["generators"].items():
        rat = np.array(g["ratio"])
        axr.step(redges, np.append(rat, rat[-1]), where="post",
                 color=COLOURS.get(key, "#444444"), lw=1.5,
                 ls="--" if dashed(key, "nu") else "-",
                 label=tex(g["label"]))
    axr.set_xlim(0.05, 1.0)
    axr.set_ylim(0.0, 3.0)
    axr.set_xlabel(tex(r"$z = E_h/\nu$"), fontsize=plotstyle.FS_XLABEL)
    axr.set_ylabel(tex("generator / NNLO"), fontsize=plotstyle.FS_YLABEL)
    axr.set_title(tex(r"ratio to the published NNLO, $E_\nu = %g$ GeV" % e_lab),
                  fontsize=plotstyle.FS_PANEL_TITLE, loc="left")
    axr.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False, loc="upper left")
    plotstyle.ticks(axr)

    out = f"{RES}/cmp_faser_sidis_mult.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    d = load()
    yields_figure(d)
    stat_figure(d)
    mult_figure(d)


if __name__ == "__main__":
    main()
