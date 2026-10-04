#!/usr/bin/env python3
"""The arXiv:2402.13318 reproduction, in one figure.

WHAT IT IS FOR (user, 2026-09-04): "produce a nice plot showing that our
pipeline reproduces the benchmark numbers, and maybe add also the POWHEG
results there to show that our pipeline is flexible enough to add other NLO MC
event generators."

THREE PANELS, and they are three separate claims rather than three views of
one:

  LEFT -- THE CROSS-SECTION IS OURS AND IT AGREES.  Our GENIE G18_02a total
  charged-current cross-section on tungsten, divided by the flux authors' own,
  for the four beams.  Two independent versions of the same comparison are
  drawn: the solid lines are file against file, against the GENIE tungsten
  cross-sections arXiv:2309.12793 ships beside its neutrino events; the
  markers are the indirect closure, our convolution over the interaction
  counts that come with the arXiv:2105.08270 flux histograms.  Both sit inside
  +-2% over three decades in energy.  This is the half of the calculation the
  benchmark contributes, and it is the half that had to be shown to be right
  before any yield could be argued about.

  CENTRE -- THE YIELDS.  Charged-current interactions in FASERnu at
  250 fb^-1, per flavour, against Table I and its uncertainty band.  Two
  entries per flavour (user, 2026-09-04): our GENIE total CC on the flux
  WITHOUT DPMJET -- our cross-section on the vendored four-generator flux,
  carried to the no-DPMJET flux by the flux authors' own count ratio, since
  only that count is published -- and POWHEG-V2 at NLO with no y window,
  whose bar is the 7-point scale envelope from reweighting every point of
  its energy ladder.  The POWHEG entry is the point of the panel beyond the
  reproduction: an NLO matched generator dropped into the same convolution
  as a cross-section table, with its own missing-higher-order uncertainty.
  The vendored-flux GENIE entry and the benchmark-region POWHEG entry were
  dropped from the figure on 2026-09-04; both remain in the report's table.

  RIGHT -- THE CHARM COMPONENT, WHICH IS WHERE nu_e LIVES.  Table I's charm
  rows against ours, computed from the POWHEG neutrino events of
  arXiv:2309.12793 -- the very calculation that table uses -- with our own
  cross-section, on the two target masses the two papers state.  The scale
  envelope is drawn as the error bar, and its RATIO to the central value is
  what actually tests the calculation: 0.55-2.75 against Table I's 0.56-2.67.

WHY THE ELECTRON FLAVOUR IS THE INTERESTING ONE.  The vendored flux averages
four hadronic models for the charm component, one of them DPMJET, whose nu_e
prediction is 3457 against SIBYLL's 901.  That is the whole of our nu_e excess
against Table I, and the right panel is the demonstration: with the charm flux
Table I actually uses, the charm rate comes out right.

Usage: analysis/plot_faser_reproduction.py
Writes results_nu/cmp_faser_reproduction.png
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

BEAMS = [("12", r"$\nu_e$", "#1f5fa9", "-"),
         ("-12", r"$\bar\nu_e$", "#1f5fa9", "--"),
         ("14", r"$\nu_\mu$", "#c0392b", "-"),
         ("-14", r"$\bar\nu_\mu$", "#c0392b", "--")]


def load(name):
    p = f"{RES}/{name}.json"
    if not os.path.exists(p):
        raise SystemExit(f"missing {p}")
    with open(p) as f:
        return json.load(f)


def panel_xsec(ax, g, fc):
    """Our GENIE over theirs, both ways of asking."""
    for pid, lab, col, ls in BEAMS:
        c = fc.get("cross_section_closure", {}).get(pid)
        if c:
            # theirs/ours is what the file comparison stores; invert so that
            # every curve on this panel is OURS over THEIRS
            ax.plot(c["E_GeV"], 1.0 / np.array(c["theirs_over_ours"]),
                    color=col, ls=ls, lw=1.8,
                    label=tex(f"{lab}, cross-section files"))
        f = g["flavours"][pid]
        ax.plot(f["closure_e"], f["closure_ratio"], color=col, ls="none",
                marker="o" if ls == "-" else "s", ms=3.0, alpha=0.55,
                mew=0)
    ax.axhline(1.0, color="#111111", lw=1.0)
    ax.axhspan(0.98, 1.02, color="#111111", alpha=0.07, lw=0)
    ax.set_xscale("log")
    ax.set_xlim(20, 6000)
    ax.set_ylim(0.945, 1.055)
    ax.set_xlabel(tex(r"$E_\nu$ [GeV]"), fontsize=plotstyle.FS_XLABEL)
    ax.set_ylabel(tex(r"our $\sigma_{\rm CC}^{W}$ / the flux authors'"),
                  fontsize=plotstyle.FS_YLABEL)
    ax.set_title(tex(r"the cross-section, ours against theirs"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False, ncol=2,
              loc="upper left")
    ax.text(0.97, 0.05, tex(r"markers: interaction counts"
                            "\n" r"lines: cross-section tables"),
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=plotstyle.FS_LEGEND - 3, color="#444444")


def panel_yields(ax, g, w):
    """Total CC interactions against Table I, plus the POWHEG-V2 entries."""
    # (label, value, (lo, hi) or None, colour, marker)
    rows = [
        ("GENIE total CC",
         lambda k: g[k]["events_on_no_dpmjet_flux"], lambda k: None,
         "#c0392b", "o"),
        (r"POWHEG-V2 NLO, $Q^2>4$, all $y$",
         lambda k: w[k]["events_fid_ally"],
         lambda k: ((w[k]["events_fid_ally_scale_lo"],
                     w[k]["events_fid_ally_scale_hi"])
                    if "events_fid_ally_scale_lo" in w[k] else None),
         "#1f5fa9", "s"),
    ]
    flav = [("nue", r"$\nu_e + \bar\nu_e$"),
            ("numu", r"$\nu_\mu + \bar\nu_\mu$")]
    for i, (k, lab) in enumerate(flav):
        c, (up, dn) = g[k]["paper"], g[k]["paper_err"]
        base = i * (len(rows) + 1.4)
        ys = base + np.arange(len(rows))
        ax.axhspan(base - 0.7, base + len(rows) - 0.3, color="#111111",
                   alpha=0.035, lw=0)
        # Table I and its band
        ax.axvline(c, color="#111111", lw=1.3, ls="--",
                   ymin=0, ymax=1, alpha=0.0)   # placeholder, drawn per band
        ax.plot([c, c], [base - 0.7, base + len(rows) - 0.3], color="#111111",
                lw=1.4, ls="--", zorder=1)
        ax.fill_betweenx([base - 0.7, base + len(rows) - 0.3], c - dn, c + up,
                         color="#111111", alpha=0.10, lw=0, zorder=0)
        for y, (name, fn, band, col, mk) in zip(ys, rows):
            v = fn(k)
            b = band(k)
            ax.plot([v], [y], marker=mk, ms=7, color=col, ls="none", zorder=3,
                    label=tex(name) if i == 0 else None)
            txt = f"{v:.0f}"
            if b is not None:
                # the 7-point band is about 1% of the yield, a few pixels on
                # this axis, so the bar is drawn OVER the marker with tall
                # caps and the band is written into the label as well
                lo, hi = b
                ax.errorbar([v], [y], xerr=[[v - lo], [hi - v]], color="#111111",
                            lw=1.0, capsize=7, capthick=1.0, ls="none",
                            zorder=4)
                txt = f"{v:.0f}$^{{+{hi - v:.0f}}}_{{-{v - lo:.0f}}}$"
            # the label goes on the side AWAY from the reference line, so it
            # never lands on top of it
            off, ha = ((-12, "right") if v <= c else (12, "left"))
            ax.annotate(tex(txt), (v, y), textcoords="offset points",
                        xytext=(off, -3), ha=ha,
                        fontsize=plotstyle.FS_LEGEND - 3, color=col)
        ax.text(0.015, (base + len(rows) / 2.0 - 0.5), tex(lab),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=plotstyle.FS_LEGEND, color="#111111")
    ax.set_yticks([])
    ax.set_ylim(-1.2, 2 * len(rows) + 1.6)
    ax.set_xlim(0, 11000)
    ax.set_xlabel(tex(r"CC interactions in FASER$\nu$ at 250 fb$^{-1}$"),
                  fontsize=plotstyle.FS_XLABEL)
    ax.set_title(tex(r"the yields, against arXiv:2402.13318 Table I"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False,
              loc="lower right")
    ax.text(0.985, 0.965, tex(r"dashed: Table I, band: its uncertainty; "
                              r"bar: 7-point scale envelope"),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=plotstyle.FS_LEGEND - 3, color="#444444")


def panel_charm(ax, fc):
    """The charm component, from the POWHEG neutrino events themselves."""
    ref = fc["reference"]["charm"]
    keys = [k for k in ("nue", "numu", "nutau")
            if k in fc["species"]["charm"]["run3_1.1t/nominal"]]
    labs = {"nue": r"$\nu_e + \bar\nu_e$", "numu": r"$\nu_\mu + \bar\nu_\mu$",
            "nutau": r"$\nu_\tau + \bar\nu_\tau$"}
    for i, k in enumerate(keys):
        r = ref[k]
        a = fc["species"]["charm"]["run3_1.1t/nominal"][k]
        b = fc["species"]["charm"]["source_1.2t/nominal"][k]
        base = i * 3.0
        # Table I charm: central plus the POWHEG Max/Min envelope
        ax.plot([r["central"]], [base + 1.0], marker="D", ms=8,
                color="#111111", ls="none", zorder=3,
                label=tex(r"arXiv:2402.13318 Table I, charm") if i == 0 else None)
        ax.plot([r["min"], r["max"]], [base + 1.0] * 2, color="#111111",
                lw=1.6, zorder=2)
        for y, v, col, mk, name in (
                (base + 0.2, a, "#c0392b", "o",
                 r"ours, 1.1 t (arXiv:2402.13318)"),
                (base - 0.55, b, "#1f5fa9", "s",
                 r"ours, 1.2 t (arXiv:2309.12793)")):
            lo, hi = v["envelope"]
            ax.plot([v["interacting"]], [y], marker=mk, ms=8, color=col,
                    ls="none", zorder=3, label=tex(name) if i == 0 else None)
            ax.plot([lo, hi], [y] * 2, color=col, lw=1.6, zorder=2)
            ax.annotate(tex(f"{v['interacting']:.0f}"), (v["interacting"], y),
                        textcoords="offset points", xytext=(0, -15),
                        ha="center", fontsize=plotstyle.FS_LEGEND - 3,
                        color=col)
        ax.text(0.015, base + 0.2, tex(labs[k]),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=plotstyle.FS_LEGEND, color="#111111")
    ax.set_yticks([])
    ax.set_xscale("log")
    ax.set_xlim(8, 4000)
    ax.set_ylim(-4.4, 3.0 * len(keys) - 0.6)
    ax.set_xlabel(tex(r"charm-induced CC interactions at 250 fb$^{-1}$"),
                  fontsize=plotstyle.FS_XLABEL)
    ax.set_title(tex(r"the charm flux, from the POWHEG events themselves"),
                 fontsize=plotstyle.FS_PANEL_TITLE)
    ax.legend(fontsize=plotstyle.FS_LEGEND - 3, frameon=False,
              loc="lower left", ncol=1)
    ax.text(0.985, 0.965, tex(r"bars: 7-point scale envelope"),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=plotstyle.FS_LEGEND - 3, color="#444444")


def main():
    g = load("faser_genie_rates")
    w = load("faser_powheg_rates")
    fc = load("forward_charm_flux")
    fig, axes = plt.subplots(1, 3, figsize=(18.6, 5.0),
                             gridspec_kw={"wspace": 0.26})
    panel_xsec(axes[0], g, fc)
    panel_yields(axes[1], g, w)
    panel_charm(axes[2], fc)
    for ax in axes:
        plotstyle.ticks(ax)
    out = f"{RES}/cmp_faser_reproduction.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
