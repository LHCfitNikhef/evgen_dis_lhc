#!/usr/bin/env python3
"""Missing-higher-order uncertainty: the 7-point scale variation, with YADISM.

The renormalisation and factorisation scales are varied around the benchmark's
choice muR = muF = Q by the standard seven points

    (xiR, xiF) in {(1,1), (2,2), (1/2,1/2), (2,1), (1,2), (1/2,1), (1,1/2)}

-- the two opposite extremes (2, 1/2) and (1/2, 2) are excluded as usual,
since they put a spuriously large logarithm in the ratio of the two scales.
The envelope of the seven is the missing-higher-order estimate.

>>> WHY THIS IS ALMOST FREE, AND IT IS THE SAME REASON THE PDF SCAN IS. <<<
A YADISM run splits into an operator build (~150 s, expensive) and a
convolution (~20 s).  The operator stores the coefficients of each power of
log(muF^2/Q^2) and log(muR^2/Q^2) SEPARATELY -- in yadism/esf/result.py the
application reads

    muF2 = Q2 * xiF**2
    a_s  = alpha_s(sqrt(Q2) * xiR)
    lnF  = log(1/xiF**2) ** o[3]
    lnR  = log(1/xiR**2) ** o[2]

so the whole scale dependence is reconstructed at convolution time from one
operator.  Seven scale points therefore cost seven convolutions, not seven
runs.  NOTE the callable is handed the ALREADY-SCALED muR, so it must be
`pdf.alphasQ(muR)` and not `pdf.alphasQ(muR * xiR)`, which would apply xiR
twice.

WHAT IS AND IS NOT COVERED.  This is the scale uncertainty of the FIXED-ORDER
NLO structure functions.  It says nothing about the parton shower, the
hadronisation model or the matching scheme, which is why the generator spread
in the benchmark is a separate and complementary number.  It also uses one
PDF (the benchmark's own): the PDF uncertainty is the companion study in
pdf_dependence.py, and the two are meant to be read side by side, because on
this observable they are not the same size.

Usage:
  mhou.py                 # all energies  (~5 min per energy)
  mhou.py --plot-only
Writes results_nu/mhou_charm_fraction.json
   and results_nu/cmp_mhou_charmfrac_E.png
"""
import json
import math
import os
import sys
import time

import numpy as np
import plotstyle          # noqa: E402 -- the house style (analysis/plotstyle.py)
plotstyle.apply()         # must precede the pyplot import
from plotstyle import tex  # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.ticker as ticker       # noqa: E402
from scipy.interpolate import RectBivariateSpline  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")
import beams  # noqa: E402

RESULTS_NU = f"{BASE}/results_nu"
OUT_JSON = f"{RESULTS_NU}/mhou_charm_fraction.json"
OUT_PNG = f"{RESULTS_NU}/cmp_mhou_charmfrac_E.png"
PDFSET = "NNPDF40_nnlo_as_01180"

# the standard 7 points; (2, 1/2) and (1/2, 2) are deliberately absent
SCALES = [(1.0, 1.0), (2.0, 2.0), (0.5, 0.5),
          (2.0, 1.0), (1.0, 2.0), (0.5, 1.0), (1.0, 0.5)]


def compute():
    import lhapdf
    lhapdf.setVerbosity(0)
    import yadism
    records = {}
    for e_gev in beams.ENERGIES:
        os.environ["BENCH_ENERGY"] = f"{e_gev:g}"
        for m in ("yadism_calc", "yadism_cc_calc", "yadism_cc_charm_calc",
                  "analyze", "beams", "runmeta", "selection"):
            sys.modules.pop(m, None)
        from yadism_calc import sf_grid_nodes
        from yadism_cc_calc import make_cards, ALPHA, make_sigma_red
        from yadism_cc_charm_calc import run

        x_nodes, q2_nodes = sf_grid_nodes()
        theory, obs = make_cards(1, x_nodes, q2_nodes)          # NLO
        kins = obs["observables"]["F2_total"]
        names = ["F2_charm", "FL_charm", "F3_charm",
                 "F2_total", "FL_total", "F3_total"]
        obs["observables"] = {n: kins for n in names}
        t0 = time.time()
        print(f"[mhou] {e_gev:g} GeV: building the YADISM operator ...",
              flush=True)
        out = yadism.run_yadism(theory, obs)
        print(f"[mhou] {e_gev:g} GeV: operator built in {time.time()-t0:.0f} s",
              flush=True)
        pdf = lhapdf.mkPDF(PDFSET, 0)
        lx, lq = np.log(x_nodes), np.log(q2_nodes)

        pts = []
        for xir, xif in SCALES:
            res = out.apply_pdf_alphas_alphaqed_xir_xif(
                pdf,
                # yadism passes the ALREADY-SCALED muR in; do not rescale here
                lambda muR: pdf.alphasQ(muR),
                lambda _muR: ALPHA, xir, xif)

            def spl(which):
                d = {}
                for n in names:
                    if not n.endswith(which):
                        continue
                    v = np.array([p["result"] for p in res[n]])
                    g = v.reshape(len(q2_nodes), len(x_nodes)).T
                    d[n[:2]] = RectBivariateSpline(lx, lq, g, kx=3, ky=3)
                return d

            s_ch = run(*make_sigma_red(spl("charm")))[0]
            s_in = run(*make_sigma_red(spl("total")))[0]
            pts.append({"xiR": xir, "xiF": xif, "sigma_charm_pb": s_ch,
                        "sigma_incl_pb": s_in, "fraction": s_ch / s_in})
            print(f"[mhou] {e_gev:g} GeV  xiR={xir:<4g} xiF={xif:<4g}  "
                  f"sigma_incl={s_in:.4f}  f_charm={100*s_ch/s_in:.3f}%",
                  flush=True)
        records[f"{e_gev:g}"] = pts

    meta = {
        "what": "7-point muR/muF variation, YADISM NLO, CC nu_mu DIS",
        "pdf": PDFSET,
        "scales": [{"xiR": a, "xiF": b} for a, b in SCALES],
        "note": ("envelope of the seven points around muR = muF = Q; the two "
                 "opposite extremes (2, 1/2) and (1/2, 2) are excluded as "
                 "usual. Fixed-order scale uncertainty only -- it does not "
                 "cover the shower, hadronisation or matching."),
        "energies_gev": list(beams.ENERGIES),
        "points": records,
    }
    os.makedirs(RESULTS_NU, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(meta, f, indent=1)
    print(f"wrote {OUT_JSON}")
    return meta


def envelope(pts, key):
    """(central, max, min) over the seven points for one quantity."""
    c = next(p[key] for p in pts if p["xiR"] == 1.0 and p["xiF"] == 1.0)
    vals = [p[key] for p in pts]
    return c, max(vals), min(vals)


def plot(meta):
    es = [float(e) for e in meta["energies_gev"]]
    frac = [envelope(meta["points"][f"{e:g}"], "fraction") for e in es]
    incl = [envelope(meta["points"][f"{e:g}"], "sigma_incl_pb") for e in es]

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.8), sharex="col",
                             gridspec_kw={"height_ratios": [2.0, 1.0],
                                          "hspace": 0.0, "wspace": 0.26})
    for col, (data, ylab, scale, title) in enumerate((
            (incl, r"$\sigma_{\rm fid}$ [pb]", 1.0,
             "Inclusive CC cross-section"),
            (frac, r"$\sigma_{\rm charm}/\sigma_{\rm inclusive}$ [%]", 100.0,
             "CC charm fraction"))):
        c = np.array([d[0] for d in data]) * scale
        hi = np.array([d[1] for d in data]) * scale
        lo = np.array([d[2] for d in data]) * scale
        ax, axr = axes[0][col], axes[1][col]
        ax.plot(es, c, "-o", color="#2a78d6", ms=6, lw=1.6,
                label=r"$\mu_R=\mu_F=Q$")
        ax.fill_between(es, lo, hi, color="#2a78d6", alpha=0.22, lw=0,
                        label="7-point envelope")
        axr.plot(es, np.ones_like(c), "-", color="#2a78d6", lw=1.2)
        axr.fill_between(es, lo / c, hi / c, color="#2a78d6", alpha=0.22, lw=0)
        axr.axhline(1.0, color="#111111", lw=0.9, ls="--")
        ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
        ax.set_ylabel(tex(ylab))
        ax.tick_params(labelbottom=False)
        _h, _l = ax.get_legend_handles_labels()
        ax.legend(_h, [tex(x) for x in _l],
                  fontsize=plotstyle.FS_LEGEND, frameon=True,
                  loc="upper left")
        axr.set_ylabel(tex("ratio to central"),
                       fontsize=plotstyle.FS_YLABEL)
        axr.set_xlabel(tex(r"$E_{\nu}$ [GeV]"),
                       fontsize=plotstyle.FS_XLABEL)
        for a in (ax, axr):
            a.set_xscale("log")
            a.set_xlim(min(es) * 0.8, max(es) * 1.25)
            a.grid(alpha=0.25, which="both")
            a.xaxis.set_minor_formatter(ticker.NullFormatter())
        axr.set_xticks(es)
        axr.set_xticklabels([f"{e:g}" for e in es])
    fig.suptitle(tex("Missing-higher-order uncertainty: "
                     "7-point scale variation") + "\n"
                 + tex("YADISM NLO, CC $\\nu_\\mu$ DIS, NNPDF4.0, "
                       "$Q^2>4$ GeV$^2$, $0.2<y<0.9$"),
                 fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, 0.015,
             "Envelope of $(\\xi_R,\\xi_F)$ over the standard seven points. "
             "Fixed-order scale uncertainty only: it does not cover the "
             "shower,\nhadronisation or matching, which is what the generator "
             "spread measures.",
             ha="center", va="bottom", fontsize=8.4, style="italic",
             color="#5d6470")
    fig.subplots_adjust(top=0.87, bottom=0.135, left=0.085, right=0.975)
    fig.savefig(OUT_PNG, dpi=150, facecolor="#ffffff")
    plt.close(fig)
    print(f"wrote {OUT_PNG}")


def main():
    if "--plot-only" in sys.argv:
        with open(OUT_JSON) as f:
            plot(json.load(f))
        return
    plot(compute())


if __name__ == "__main__":
    main()
