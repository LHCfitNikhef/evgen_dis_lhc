#!/usr/bin/env python3
"""PDF dependence of the charm fraction, with PDF uncertainties.

THE SCHEME IS FONLL (user, 2026-08-28), which is the benchmark's convention
everywhere charm appears, so this study no longer stands apart from the charm
tables and the charm-fraction figures.  It used to be ZM-VFNS, labelled as an
exception; `--scheme zm` still reproduces that and writes to its own file.

WHAT MADE THE MOVE AFFORDABLE.  The earlier note here estimated three times
the cost, on the reasoning that FONLL needs three runs per member instead of
one AND that the denominator has to move with the numerator.  Both halves are
true and the conclusion was still wrong, for two reasons:

  * PineAPPL.  The three schemes differ only in the theory card, so each is
    one grid build (~2 min, PDF-INDEPENDENT) and the per-member cost is three
    convolutions at ~5 ms instead of one.  What used to be three yadism runs
    per member is now three table lookups.

  * The denominator is FREE.  In FONLL only the HEAVY-quark part changes
    scheme -- the light contributions are the ZM ones, untouched -- so with
    d = f_thr * (FFNS_charm - FFN0_charm) the mass correction is ONE array
    and it enters both:

        F_charm^FONLL = F_charm^ZM + d
        F_incl^FONLL  = F_incl^ZM  + d

    the same identity `yadism_fonll_inclusive.py` uses for the reference
    curves.  So moving the denominator costs no extra grid at all, and the
    numerator and denominator stay in ONE scheme, which is the whole point --
    a fraction formed across two schemes means nothing.

Measured: 12 grids per energy instead of 6, and the run went from about two
hours to about seven minutes.

ONLY CHARM IS MOVED TO FONLL.  No bottom FONLL reference exists in this
benchmark, so bottom stays massless inside the ZM terms above.  At FASER
kinematics that is a per-mille effect on the inclusive rate, but it IS an
approximation and is recorded in the output as `heavy_in_fonll: ["charm"]`
rather than left implicit -- same convention as the inclusive reference.

FONLL HERE ALWAYS MEANS DAMPED (CONVENTIONS.md rule 3): the threshold damping
multiplies the mass correction d, and no undamped variant is produced.


sigma_charm/sigma_inclusive in CC neutrino DIS at NLO, computed with YADISM
for several PDF sets across the beam-energy scan, each with its own PDF
uncertainty band, and everything shown as a ratio to the benchmark's own
NNPDF4.0 reference.

WHY THIS IS AFFORDABLE.  A YADISM run splits into two very unequal halves:

    out = yadism.run_yadism(theory, observables)     # ~80 s, PDF-INDEPENDENT
    res = out.apply_pdf_alphas_alphaqed_xir_xif(...) # ~9 s,  per PDF member

The operator depends only on the theory card and the kinematic grid, so it is
built ONCE PER ENERGY and then convoluted with every member of every set.
Naively re-running YADISM per member would be ~280 x 90 s x 3 energies, i.e.
a day; this way it is one build plus the convolutions.

AND THE CHARM AND INCLUSIVE OBSERVABLES ARE REQUESTED IN THE SAME RUN, which
is not merely an optimisation.  The fraction's uncertainty is NOT obtainable
from the two uncertainties separately: numerator and denominator move together
under a PDF variation, and strongly, since they share the same quark
densities.  The ratio has to be formed MEMBER BY MEMBER and the uncertainty
taken of the ratio.  Doing it the other way round would roughly double the
band.

UNCERTAINTIES ARE TAKEN FROM LHAPDF ITSELF (`PDFSet.uncertainty`), because the
sets do not agree on how to compute them: NNPDF4.0 is a Monte Carlo replica
set, CT18/MSHT20/ATLASpdf21 are Hessian, and ATLASpdf21 carries extra
parameter variations appended to its error type.  Everything is rescaled to a
COMMON 68% confidence level, which matters: CT18NNLO is published at 90%, so
quoting its raw band beside a 68% one would make it look 1.64x more uncertain
than it is.

GRV98 carries NO uncertainty at all, and that is not a gap in this study.  It
is a 1998 set, from before PDF uncertainties existed as a deliverable: it was
published as a single best-fit parametrisation with no error members, so there
is nothing to compute rather than something omitted.  It is included because
it is what GENIE's default tune uses, so it shows how far the FASER baseline
sits from a modern set -- but it is also a LO set being used with an NLO
calculation, which is a second mismatch worth remembering when reading it.

ONE THEORY CARD THROUGHOUT.  Only the PDF is varied: the flavour thresholds,
alpha_em and the scale choice stay at the benchmark's values, while alpha_s
comes from each set (LHAPDF's own alphasQ).  That is the point of a PDF
dependence study -- varying the scheme at the same time would confound it.

Usage:
  pdf_dependence.py                 # all energies, all sets  (~2 h)
  pdf_dependence.py --quick         # central members only, no bands (~5 min)
  pdf_dependence.py --plot-only     # re-plot from the stored JSON
Writes results_nu/pdf_dependence_charm_fraction.json
   and results_nu/cmp_pdf_dependence_charmfrac_E.png
"""
import importlib
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

# BOTH CURRENTS (CONVENTIONS.md rule 2b).  The CC study came first and keeps its
# original output names so nothing published moves; the muon NC one is tagged.
# Only the module names differ -- yadism_cc_calc/yadism_cc_charm_calc against
# yadism_calc/yadism_charm_calc -- so the current is a parameter here rather
# than a second copy of the script, which is how the two halves of this
# repository have drifted before.
# BENCH_CURRENT lets another module import this one for its helpers without
# having to forge sys.argv; the command line still wins when it is given.
CURRENT = (sys.argv[sys.argv.index("--current") + 1]
           if "--current" in sys.argv
           else os.environ.get("BENCH_CURRENT", "nu"))
if CURRENT not in ("mu", "nu"):
    sys.exit("--current takes mu or nu")
_CC = CURRENT == "nu"
CALC_MOD = "yadism_cc_calc" if _CC else "yadism_calc"
CHARM_MOD = "yadism_cc_charm_calc" if _CC else "yadism_charm_calc"

RESULTS_NU = f"{BASE}/results_nu" if _CC else f"{BASE}/results"
_SUF = "" if _CC else "_mu"
OUT_JSON = f"{RESULTS_NU}/pdf_dependence_charm_fraction{_SUF}.json"
OUT_PNG = f"{RESULTS_NU}/cmp_pdf_dependence_charmfrac{_SUF}_E.png"
BEAM_LABEL = r"$E_{\nu}$ [GeV]" if _CC else r"$E_{\mu}$ [GeV]"
TITLE_CUR = "CC" if _CC else "muon NC"

# FONLL BY DEFAULT (user, 2026-08-28).  `--scheme zm` writes to its own file
# so the two can coexist and neither silently overwrites the other.
SCHEME = (sys.argv[sys.argv.index("--scheme") + 1]
          if "--scheme" in sys.argv else "fonll")
if SCHEME not in ("fonll", "zm"):
    sys.exit("--scheme must be fonll or zm")
_FONLL = SCHEME == "fonll"
# The scheme is named in every user-visible string this module writes --
# CONVENTIONS.md rule 3.
SCHEME_LABEL = "FONLL" if _FONLL else "ZM-VFN"
if not _FONLL:
    OUT_JSON = OUT_JSON.replace(".json", "_zm.json")
    OUT_PNG = OUT_PNG.replace(".png", "_zm.png")

# The three theory cards FONLL is assembled from: F = ZM + f_thr*(FFNS - FFN0).
# Only the charm flavour needs the two massive runs -- see the docstring on why
# the same correction serves the inclusive denominator.
FONLL_PARTS = [("zm", "ZM-VFNS", 5), ("ffns", "FONLL-FFNS", 3),
               ("ffn0", "FONLL-FFN0", 3)]
CHARM_NAMES = ["F2_charm", "FL_charm", "F3_charm"]

# GRV98 IS DROPPED ON THE NEUTRAL CURRENT ONLY, and the asymmetry is physical
# rather than arbitrary.  Measured at 400 GeV: the NC charm fraction comes out
# at -3.469%, and the 4 TeV x-spectrum study finds a negative charm structure
# function in every bin.  The CC value is a perfectly ordinary 11.797%.
#
# The reason is which channel carries the charm.  CC charm is quark-initiated,
# s -> c, and stays positive with any sensible density.  NC charm at NLO leans
# on gamma* g -> c cbar together with its subtraction term, and a 1998
# LEADING-ORDER set whose heavy-flavour content and alpha_s do not match an
# nf = 5 NLO configuration cannot make those cancel.
#
# It is KEPT on the CC side, where it is meaningful and shows how far GENIE's
# default tune sits from a modern set -- which is why it was included at all.
SKIP_ON_NC = {"GRV98lo"}

# The benchmark's own set is the reference every ratio is taken to.
REFERENCE = "NNPDF40_nnlo_as_01180"

# (LHAPDF name, display label, colour)
PDF_SETS = [
    ("NNPDF40_nnlo_as_01180", "NNPDF4.0 NNLO", "#2a78d6"),
    ("CT18NNLO",              "CT18 NNLO",     "#eb6834"),
    ("MSHT20nnlo_as118",      "MSHT20 NNLO",   "#3fa66a"),
    ("ATLASpdf21_T1",         "ATLASpdf21",    "#8e5bd0"),
    # ABMP16 (user, 2026-08-27).  The 5-flavour NNLO member, so it matches the
    # others in both order and flavour scheme; ABMP16_3_* would be a
    # 3-flavour fit and would not be comparable here.  It matters for this
    # study in particular because ABMP determines the strange density from
    # dimuon data, which is exactly the channel the CC charm fraction probes.
    ("ABMP16als118_5_nnlo",   "ABMP16 NNLO",   "#d4a017"),   # alpha_s = 0.118 (user, 2026-10-05)
    ("GRV98lo",               "GRV98 LO",      "#c2317b"),
]

# One common confidence level.  CT18 is published at 90% and everything else
# at 68%; without this the CT18 band would read 1.64x too wide beside them.
CL68 = 100.0 * math.erf(1.0 / math.sqrt(2.0))     # 68.268...


def observables_card(pto, x_nodes, q2_nodes):
    """Theory + observable cards carrying BOTH the charm and total SFs.

    Requested together so a single convolution yields the numerator and the
    denominator of the fraction for the same PDF member -- see the module
    docstring on why the ratio must be formed member by member.
    """
    make_cards = importlib.import_module(CALC_MOD).make_cards
    theory, obs = make_cards(pto, x_nodes, q2_nodes)
    kins = obs["observables"]["F2_total"]
    names = ["F2_charm", "FL_charm", "F3_charm",
             "F2_total", "FL_total", "F3_total"]
    obs["observables"] = {n: kins for n in names}
    return theory, obs, names


def fonll_combine(res_zm, res_ffns, res_ffn0, names, q2_nodes, nx):
    """ZM results plus the damped charm mass correction, in `res` shape.

        d             = f_thr * (FFNS_charm - FFN0_charm)
        F_charm^FONLL = F_charm^ZM + d
        F_total^FONLL = F_total^ZM + d

    THE SAME d ENTERS BOTH, because in FONLL only the heavy-quark part of the
    structure function changes scheme; the light part is the ZM one and is
    untouched.  That identity is what makes the denominator free, and it is
    the one `yadism_fonll_inclusive.py` already uses for the reference curves.

    THE DAMPING IS A FUNCTION OF Q2 ONLY, and the flat result arrays run
    Q2-major (index = iq*nx + ix -- the ordering `splines_from` inverts with
    `.reshape(len(q2), len(x)).T`), so f_thr is repeated nx times per node
    rather than tiled.  Getting that backwards would damp along the wrong
    axis and still produce a smooth, plausible surface.
    """
    from yadism_charm_calc import damping_factor
    f_thr = np.repeat(damping_factor(np.asarray(q2_nodes, dtype=float)), nx)
    out = {}
    for name in names:
        zm = np.array([p["result"] for p in res_zm[name]])
        cname = name.replace("_total", "_charm")
        d = (np.array([p["result"] for p in res_ffns[cname]])
             - np.array([p["result"] for p in res_ffn0[cname]])) * f_thr
        out[name] = [{"result": float(v)} for v in (zm + d)]
    return out


def splines_from(res, names, x_nodes, q2_nodes, which):
    """RectBivariateSplines for one flavour group ('charm' or 'total')."""
    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    out = {}
    for name in names:
        if not name.endswith(which):
            continue
        vals = np.array([p["result"] for p in res[name]])
        grid = vals.reshape(len(q2_nodes), len(x_nodes)).T
        out[name[:2]] = RectBivariateSpline(lx, lq, grid, kx=3, ky=3)
    return out


def fraction_for(res, names, x_nodes, q2_nodes):
    """(sigma_charm, sigma_inclusive, fraction) in pb for one PDF member."""
    make_sigma_red = importlib.import_module(CALC_MOD).make_sigma_red
    run = importlib.import_module(CHARM_MOD).run
    # both currents' run() return sigma first, so [0] is the fiducial
    # cross-section in each; the NC one returns more elements after it
    s_ch = run(*make_sigma_red(splines_from(res, names, x_nodes,
                                            q2_nodes, "charm")))[0]
    s_in = run(*make_sigma_red(splines_from(res, names, x_nodes,
                                            q2_nodes, "total")))[0]
    return s_ch, s_in, (s_ch / s_in if s_in else float("nan"))


def band(pset, values):
    """(central, err_plus, err_minus) at a COMMON 68% CL, via LHAPDF.

    LHAPDF knows each set's error type -- Monte Carlo replicas for NNPDF4.0,
    Hessian for CT18/MSHT20, Hessian plus appended parameter variations for
    ATLASpdf21 -- and the confidence level it was published at.  Asking it
    rather than reimplementing the formulas is what keeps CT18's 90% band
    comparable with everyone else's 68%.
    """
    # A set with one member has no error band to compute.  GRV98 is the case
    # here: it predates PDF uncertainties, so it ships a single best-fit
    # parametrisation and there is nothing to evaluate.
    if len(values) < 2:
        return values[0], 0.0, 0.0
    try:
        u = pset.uncertainty(list(values), CL68)
        return u.central, u.errplus, u.errminus
    except Exception as exc:                                  # noqa: BLE001
        # Never silently fall back to "no error" -- a missing band is
        # indistinguishable from a zero one on a plot.
        print(f"    !! LHAPDF uncertainty() failed ({exc}); "
              f"using the symmetric Hessian formula instead")
        c = values[0]
        d = math.sqrt(sum((values[i] - c) ** 2
                          for i in range(1, len(values)))) 
        return c, d, d


def compute(quick=False):
    import lhapdf
    lhapdf.setVerbosity(0)
    from yadism_calc import sf_grid_nodes
    ALPHA = importlib.import_module(CALC_MOD).ALPHA
    import yadism

    records = {}
    for e_gev in beams.ENERGIES:
        os.environ["BENCH_ENERGY"] = f"{e_gev:g}"
        # the kinematics and the operator depend on the beam, so both are
        # rebuilt per energy; the modules read BENCH_ENERGY at import, so
        # they are re-imported here rather than cached across energies
        for m in ("yadism_calc", "yadism_cc_calc", "yadism_cc_charm_calc",
                  "yadism_charm_calc",
                  "analyze", "beams", "runmeta", "selection"):
            sys.modules.pop(m, None)
        from yadism_calc import sf_grid_nodes as _nodes
        x_nodes, q2_nodes = _nodes()
        theory, obs, names = observables_card(1, x_nodes, q2_nodes)   # NLO
        t0 = time.time()
        # PINEAPPL BY DEFAULT (user, 2026-08-28).  The grid stores the same
        # PDF-independent information and convolutes it ~46x faster, matching
        # yadism to 5e-12 -- verified per observable by
        # `pineappl_grids.py check`, not assumed.  --no-grids falls back to
        # calling yadism per member, which is what this did before and is kept
        # so the two can be compared at any time.
        use_grids = "--no-grids" not in sys.argv
        # one entry per scheme part; ZM alone when --scheme zm
        parts = FONLL_PARTS if _FONLL else FONLL_PARTS[:1]
        gkey, outs = {}, {}
        if use_grids:
            import pineappl_grids as pg
            print(f"[pdf] {e_gev:g} GeV: building/loading PineAPPL grids "
                  f"({len(parts)} scheme(s)) ...", flush=True)
            for tag, fns, nfff in parts:
                # the massive runs are only ever asked for the CHARM
                # observables -- the inclusive correction is the charm one
                want = names if tag == "zm" else CHARM_NAMES
                gkey[tag], _xn, _qn = pg.build(CURRENT, e_gev, 1, want,
                                               x_nodes=x_nodes,
                                               q2_nodes=q2_nodes,
                                               fns=fns, nfff=nfff)
        else:
            print(f"[pdf] {e_gev:g} GeV: building the YADISM operator(s) "
                  f"(PDF-independent) ...", flush=True)
            for tag, fns, nfff in parts:
                th, ob = observables_card(1, x_nodes, q2_nodes)[:2]
                th["FNS"], th["NfFF"] = fns, nfff
                if tag != "zm":
                    ob["observables"] = {n: ob["observables"][n]
                                         for n in CHARM_NAMES}
                outs[tag] = yadism.run_yadism(th, ob)
        print(f"[pdf] {e_gev:g} GeV: ready in "
              f"{time.time()-t0:.0f} s", flush=True)

        for setname, label, colour in PDF_SETS:
            if not _CC and setname in SKIP_ON_NC:
                continue
            pset = lhapdf.getPDFSet(setname)
            nmem = 1 if quick else pset.size
            vals, sch, sin = [], [], []
            t1 = time.time()
            for imem in range(nmem):
                pdf = pset.mkPDF(imem)
                per = {}
                for tag, _fns, _nfff in parts:
                    want = names if tag == "zm" else CHARM_NAMES
                    if use_grids:
                        per[tag] = pg.res_like(gkey[tag], want, pdf)
                    else:
                        per[tag] = outs[tag].apply_pdf_alphas_alphaqed_xir_xif(
                            pdf, lambda muR, _p=pdf: _p.alphasQ(muR),
                            lambda _muR: ALPHA, 1.0, 1.0)
                res = (fonll_combine(per["zm"], per["ffns"], per["ffn0"],
                                     names, q2_nodes, len(x_nodes))
                       if _FONLL else per["zm"])
                a, b, f = fraction_for(res, names, x_nodes, q2_nodes)
                sch.append(a); sin.append(b); vals.append(f)
            c, ep, em = band(pset, vals)
            c_ch, _, _ = band(pset, sch)
            c_in, _, _ = band(pset, sin)
            records.setdefault(setname, {"label": label, "colour": colour,
                                         "members": nmem, "points": {}})
            records[setname]["points"][f"{e_gev:g}"] = {
                "fraction": c, "err_plus": ep, "err_minus": em,
                "sigma_charm_pb": c_ch, "sigma_incl_pb": c_in,
                "members_used": nmem,
            }
            print(f"[pdf] {e_gev:g} GeV  {label:16s} "
                  f"f = {100*c:.3f} +{100*ep:.3f} -{100*em:.3f} %   "
                  f"({nmem} members, {time.time()-t1:.0f} s)", flush=True)

    meta = {
        "what": (f"sigma_charm/sigma_inclusive, {TITLE_CUR} DIS, "
                 f"YADISM NLO in "
                 f"{'FONLL' if _FONLL else 'ZM-VFNS'}"),
        "scheme": ("FONLL general-mass, damped; charm only, bottom stays "
                   "massless" if _FONLL else "ZM-VFNS"),
        "heavy_in_fonll": ["charm"] if _FONLL else [],
        "current": CURRENT,
        "reference": REFERENCE,
        "confidence_level_percent": CL68,
        "note": ("one theory card throughout; only the PDF is varied. "
                 "alpha_s is taken from each set. The fraction is formed "
                 "member by member before the uncertainty is evaluated, "
                 "because numerator and denominator are correlated."),
        "energies_gev": list(beams.ENERGIES),
        "quick": quick,
        "sets": records,
    }
    os.makedirs(RESULTS_NU, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(meta, f, indent=1)
    print(f"wrote {OUT_JSON}")
    return meta


def plot(meta):
    es = [float(e) for e in meta["energies_gev"]]
    ref = meta["sets"].get(REFERENCE)
    if ref is None:
        sys.exit("the reference set is missing from the stored results")

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8.4, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.0})

    def series(rec):
        x, y, ep, em = [], [], [], []
        for e in es:
            p = rec["points"].get(f"{e:g}")
            if not p:
                continue
            x.append(e); y.append(p["fraction"])
            ep.append(p["err_plus"]); em.append(p["err_minus"])
        return (np.array(x), np.array(y), np.array(ep), np.array(em))

    rx, ry, rep, rem = series(ref)
    for setname, rec in meta["sets"].items():
        x, y, ep, em = series(rec)
        if not len(x):
            continue
        c = rec["colour"]
        ax.plot(x, 100 * y, "-o", color=c, ms=6, lw=1.5, label=rec["label"])
        if np.any(ep + em):
            ax.fill_between(x, 100 * (y - em), 100 * (y + ep),
                            color=c, alpha=0.20, lw=0)
        # ratio to the reference, with THIS set's own band carried over.  The
        # reference's band is drawn once, around 1, as a grey envelope.
        r = y / ry
        axr.plot(x, r, "-o", color=c, ms=6, lw=1.5)
        if np.any(ep + em):
            axr.fill_between(x, (y - em) / ry, (y + ep) / ry,
                             color=c, alpha=0.20, lw=0)
    axr.fill_between(rx, 1 - rem / ry, 1 + rep / ry, color="#888888",
                     alpha=0.25, lw=0, zorder=0)
    axr.axhline(1.0, color="#111111", lw=1.0, ls="--")

    for a in (ax, axr):
        a.set_xscale("log")
        a.set_xlim(min(es) * 0.8, max(es) * 1.25)
        a.grid(alpha=0.25, which="both")
        a.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_ylabel(tex(r"$\sigma_{\rm charm}/\sigma_{\rm inclusive}$  [%]"),
                  fontsize=plotstyle.FS_YLABEL)
    _h, _l = ax.get_legend_handles_labels()
    ax.legend(_h, [tex(x) for x in _l], fontsize=plotstyle.FS_LEGEND,
              frameon=True, loc="upper left")
    ax.tick_params(labelbottom=False)
    axr.set_ylabel(tex("ratio to") + "\n" + tex("NNPDF4.0 NNLO"),
                   fontsize=plotstyle.FS_YLABEL)
    axr.set_xlabel(tex(BEAM_LABEL), fontsize=plotstyle.FS_XLABEL)
    axr.set_xticks(es)
    axr.set_xticklabels([f"{e:g}" for e in es])
    # THE SCHEME IS NAMED IN THE TITLE (user, 2026-08-27; CONVENTIONS.md rule 3).
    # It is read off the stored result rather than off this process, so a
    # --plot-only rerun cannot label a ZM file as FONLL or the reverse.
    sch = "FONLL" if "FONLL" in meta.get("what", "") else "ZM-VFN"
    fig.suptitle(tex(f"PDF dependence of the {TITLE_CUR} charm fraction")
                 + "\n"
                 + tex(f"YADISM NLO ({sch}), "
                       "$Q^2>4$ GeV$^2$, $0.2<y<0.9$; "
                       "68% CL PDF errors"),
                 fontsize=plotstyle.FS_SUPTITLE)
    fig.text(0.5, 0.015,
             "One theory card throughout: only the PDF is varied, with "
             "$\\alpha_s$ taken from each set. The ratio is formed member by "
             "member,\nso the correlation between numerator and denominator "
             "is kept. GRV98 predates PDF uncertainties and has none to show."
             + ("\nCharm is massive: numerator and denominator are both in "
                "FONLL, so the fraction is formed inside one scheme. Bottom "
                "stays massless." if sch == "FONLL" else
                "\nMassless charm, unlike the charm tables elsewhere: the "
                "absolute fractions are a few per cent high, the ratios "
                "between sets are not affected."),
             ha="center", va="bottom", fontsize=8.4, style="italic",
             color="#5d6470")
    fig.subplots_adjust(top=0.89, bottom=0.135, left=0.115, right=0.97)
    fig.savefig(OUT_PNG, dpi=150, facecolor="#ffffff")
    plt.close(fig)
    print(f"wrote {OUT_PNG}")


def main():
    if "--plot-only" in sys.argv:
        with open(OUT_JSON) as f:
            plot(json.load(f))
        return
    plot(compute(quick="--quick" in sys.argv))


if __name__ == "__main__":
    main()
