#!/usr/bin/env python3
"""The 7-point scale band on the DIFFERENTIAL distributions, bin by bin.

Paper plot 2 shows differential distributions against the analytic reference
but had no theory uncertainty on them (user, 2026-08-30: "can we add MHOUs of
the YADISM calculation, in full analogy to what we do in paper plot 1").
`mhou_sigma_fonll.py` gives the band on the INTEGRATED cross-section only --
one number per energy -- so this is its differential counterpart.

WHAT IT DOES.  Convolves the same stored PineAPPL operator at each of the
seven scale points, integrates each into the analysis's own bins, and records
the envelope per bin.  Nothing new is computed: the operator is
scale-independent by construction, which is exactly why the seven points cost
a convolution each rather than a run each.

THE SCHEME IS A CHOICE, AND IT FOLLOWS THE CURVE THE BAND SURROUNDS.  A band
in one scheme drawn around a curve in another is an uncertainty on a different
calculation, so `--scheme` exists and the output records what it did:

  * `--scheme zm` (default) for the NLO reference, which paper plot 2
    normalises to.  That reference is massless because every generator on the
    figure computes massless charm.
  * `--scheme fonll` for the muon NNLO curve, which IS in FONLL -- paper plot
    3 draws NNLO in the scheme each current has at that order.  Assembled
    from the same three cards as everywhere else in this benchmark
    (ZM-VFNS/5, FONLL-FFNS/3, FONLL-FFN0/3) via pdf_dependence.fonll_combine,
    so the seven scale points are seven combinations, not seven runs.

NO NEUTRINO DIS AT NNLO IN FONLL (user, 2026-08-29): the massive
charged-current coefficient functions do not exist at that order, so
`--scheme fonll --current nu --pto 2` is REFUSED rather than quietly
returning a plausible band for a calculation nobody performed.  The neutrino
NNLO band is in ZM-VFN, which is what there is, and paper plot 2 says so in
its legend.

THE SCALES ARE FLOORED, muF and muR both, at the parton distribution's own
Q0^2.  This is the standing rule (see mhou_sigma_fonll.py): without it the
downward variations leave the grid, and the band becomes an extrapolation
artefact rather than a scale variation -- which once faked a 45% collapse of
F2 at the bottom of a fiducial region.  Nodes needing the same effective
(xiR, xiF) are grouped, so this costs one convolution per DISTINCT pair.

>>> THE xiF = 1/2 POINTS ARE KEPT HERE, and the output says so. <<<  On the
INTEGRATED muon cross-section they are dropped, because there the lower edge
moves by a factor of three with the parton distribution and the prescription
and cannot be quoted.  A per-bin band is a different object: it is read as a
shape uncertainty rather than a number to quote, and dropping the points would
make the band narrower in exactly the bins where the instability lives, which
is the worst place to hide it.  Both envelopes are written -- `rel_hi`/`rel_lo`
with all seven, `rel_hi_nohalf`/`rel_lo_nohalf` without -- so the figure can
choose and the choice is visible.

Usage:  analysis/mhou_diff.py [--current mu|nu] [--pto 1] [--scheme zm|fonll]
                             [--flavour total|charm]
        -> results{,_nu}/mhou_diff_{zm,fonll}[_charm][_pto2]
                             [_<region>_<target>[_tmc]][_<tag>][_q2minN].json

REGION, TARGET AND TMC come from the environment, as for every calculator
($BENCH_SELECTION, $BENCH_TARGET, $BENCH_TMC).  The paper-plots band is
    BENCH_SELECTION=q4w3 BENCH_TARGET=W BENCH_TMC=3 analysis/mhou_diff.py ...
and the PDF is target.pdfset() (the free-nucleon tungsten average there).
None of the three was in the filename or the PDF before 2026-09-14, so a
q4w3 tungsten run would have written a PROTON-PDF band over the earlier file.
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

CURRENT = (sys.argv[sys.argv.index("--current") + 1]
           if "--current" in sys.argv else "mu")
if CURRENT not in ("mu", "nu"):
    sys.exit("--current takes mu or nu")
PTO = int(sys.argv[sys.argv.index("--pto") + 1]) if "--pto" in sys.argv else 1
SCHEME = (sys.argv[sys.argv.index("--scheme") + 1]
          if "--scheme" in sys.argv else "zm")
if SCHEME not in ("zm", "fonll"):
    sys.exit("--scheme takes zm or fonll")
# The refusal above, in code.  It fires before any grid is built, because the
# failure this guards against is not an error but a plausible-looking number.
if SCHEME == "fonll" and CURRENT == "nu" and PTO >= 2:
    sys.exit("no neutrino DIS at NNLO in FONLL: the massive charged-current\n"
             "coefficient functions do not exist at that order.  The neutrino\n"
             "NNLO band is ZM-VFN -- run without --scheme fonll.")
SCHEME_LABEL = "FONLL" if SCHEME == "fonll" else "ZM-VFN"
# --flavour charm gives the band on the CHARM cross-section, for paper plot 7
# (user, 2026-08-31: "we need to add the corresponding scale bands").  It is
# the same calculation with the charm structure functions in place of the
# total ones: both charm calc modules integrate with the SAME make_sigma_red
# and integrate_* as the inclusive ones, so only the observable names change.
# The CKM row in yadism_cc_charm_calc belongs to its direct LHAPDF
# cross-check, NOT to the yadism path, so nothing extra is owed on the CC.
FLAVOUR = (sys.argv[sys.argv.index("--flavour") + 1]
           if "--flavour" in sys.argv else "total")
if FLAVOUR not in ("total", "charm"):
    sys.exit("--flavour takes total or charm")

SCALES = [(1.0, 1.0), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5),
          (2.0, 2.0), (0.5, 0.5)]
OBS = ("Q2", "xbj", "y", "Emu", "theta", "nu")


def main():
    import importlib
    import lhapdf
    lhapdf.setVerbosity(0)
    import beams
    import pineappl_grids as pg
    import target

    e_gev = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
    os.environ["BENCH_ENERGY"] = f"{e_gev:g}"
    for m in ("yadism_calc", "yadism_cc_calc", "analyze", "beams",
              "runmeta", "selection"):
        sys.modules.pop(m, None)
    calc = importlib.import_module(
        "yadism_cc_calc" if CURRENT == "nu" else "yadism_calc")
    from yadism_calc import sf_grid_nodes
    from analyze import BINS

    x_nodes, q2_nodes = sf_grid_nodes()
    names = (["F2_charm", "FL_charm", "F3_charm"] if FLAVOUR == "charm"
             else ["F2_total", "FL_total", "F3_total"])
    # FONLL wants the SAME observables on the SAME nodes in three schemes; ZM
    # wants one.  pdf_dependence owns the part list and the combination, so
    # the assembly is not written a second time here (it has diverged between
    # copies before -- see tools/check_parser_drift.py for the same lesson).
    pdd = importlib.import_module("pdf_dependence")
    parts = pdd.FONLL_PARTS if SCHEME == "fonll" else pdd.FONLL_PARTS[:1]
    keys = {}
    for tag, fns, nfff in parts:
        want = names if tag == "zm" else pdd.CHARM_NAMES
        keys[tag], _x, _q = pg.build(CURRENT, e_gev, PTO, want,
                                     x_nodes=x_nodes, q2_nodes=q2_nodes,
                                     fns=fns, nfff=nfff)

    # THE TARGET PICKS THE SET (NNPDF40_nnlo_as_01180 for the default proton)
    pdfset = target.pdfset()
    pdf = lhapdf.mkPDF(pdfset, 0)
    q2a = np.asarray(q2_nodes, dtype=float)
    nx = len(x_nodes)
    floor = np.sqrt(pdf.q2Min / q2a)

    def res_floored(gkey, want, xir, xif):
        eff_r = np.round(np.maximum(xir, floor), 12)
        eff_f = np.round(np.maximum(xif, floor), 12)
        out = {n: [None] * (len(q2a) * nx) for n in want}
        for pr, pf in sorted(set(zip(eff_r.tolist(), eff_f.tolist()))):
            r = pg.res_like(gkey, want, pdf, float(pr), float(pf))
            sel = np.where((eff_r == pr) & (eff_f == pf))[0]
            for n in want:
                for iq in sel:
                    for ix in range(nx):
                        out[n][iq * nx + ix] = r[n][iq * nx + ix]
        return out

    def scheme_res(xir, xif):
        """Structure functions at one scale point, in the chosen scheme."""
        per = {tag: res_floored(keys[tag],
                                names if tag == "zm" else pdd.CHARM_NAMES,
                                xir, xif)
               for tag, _f, _n in parts}
        if SCHEME != "fonll":
            return per["zm"]
        # fonll_combine forms cname = name.replace("_total", "_charm") to
        # find the massive pieces.  With charm names that replace is a no-op
        # and the identity it applies, F^FONLL = F^ZM + f_thr*(FFNS - FFN0),
        # is the same one -- which is why one call serves both flavours.
        return pdd.fonll_combine(per["zm"], per["ffns"], per["ffn0"],
                                 names, q2_nodes, nx)

    def splines(res):
        lx, lq = np.log(x_nodes), np.log(q2_nodes)
        return {n[:2]: __import__(
            "scipy.interpolate", fromlist=["RectBivariateSpline"]
        ).RectBivariateSpline(
            lx, lq,
            np.array([p["result"] for p in res[n]]).reshape(
                len(q2_nodes), nx).T, kx=3, ky=3) for n in names}

    per_scale, sigmas = [], []
    for xir, xif in SCALES:
        sp = splines(scheme_res(xir, xif))
        dsig_dxdq2, dsig_dq2dy = calc.make_sigma_red(sp)
        hists = {}
        sig, _means = calc.integrate_q2y(dsig_dq2dy, hists)
        calc.integrate_xq2(dsig_dxdq2, hists)
        calc.integrate_ytheta(dsig_dq2dy, hists)
        per_scale.append({k: np.asarray(hists[k], dtype=float) for k in OBS})
        sigmas.append(sig)
        print(f"[mhou_diff] {CURRENT} {e_gev:g} GeV  ({xir:g}, {xif:g})  "
              f"sigma = {sig:.5g}", flush=True)

    out = {
        "what": (f"7-point muR/muF band on the DIFFERENTIAL "
                 f"{'CHARM ' if FLAVOUR == 'charm' else ''}distributions, "
                 f"{CURRENT} DIS, YADISM PTO={PTO} in {SCHEME_LABEL}"),
        "flavour": FLAVOUR,
        "scheme": ("FONLL general-mass, damped; charm only"
                   if SCHEME == "fonll" else "ZM-VFN"),
        "heavy_in_fonll": ["charm"] if SCHEME == "fonll" else [],
        "current": CURRENT, "energy_gev": e_gev, "pto": PTO,
        "pdf": pdfset, "q2min_pdf": pdf.q2Min,
        "scale_floor": "muF^2 = max(xiF^2 Q^2, Q0^2), muR likewise",
        "xif_half": ("KEPT in rel_hi/rel_lo; also given without them in "
                     "rel_hi_nohalf/rel_lo_nohalf. On the integrated muon "
                     "cross-section they are dropped as unquotable, but a "
                     "per-bin band is a shape statement and dropping them "
                     "would narrow it exactly where the instability is"),
        "scales": [{"xiR": a, "xiF": b} for a, b in SCALES],
        "sigma_fid_pb": sigmas[0],
        "hists": {},
    }
    # Region, target and TMC stamped whenever any is away from its default
    # (absent otherwise, so a earlier re-run is byte-identical).
    import selection
    sel = selection.get()
    if (target.current() != target.DEFAULT or target.tmc() != 0
            or sel.name != selection.DEFAULT):
        out["selection"] = sel.as_dict()
        out["pdf_set"] = pdfset
        out["target"] = target.current()
        out["tmc"] = target.tmc()
    keep = [i for i, (_r, f) in enumerate(SCALES) if f != 0.5]
    for k in OBS:
        w = np.diff(np.asarray(BINS[k], dtype=float))
        cen = per_scale[0][k]
        stack = np.array([p[k] for p in per_scale])
        with np.errstate(divide="ignore", invalid="ignore"):
            hi = np.where(cen > 0, stack.max(axis=0) / cen - 1.0, 0.0)
            lo = np.where(cen > 0, stack.min(axis=0) / cen - 1.0, 0.0)
            hin = np.where(cen > 0, stack[keep].max(axis=0) / cen - 1.0, 0.0)
            lon = np.where(cen > 0, stack[keep].min(axis=0) / cen - 1.0, 0.0)
        out["hists"][k] = {
            "edges": np.asarray(BINS[k], dtype=float).tolist(),
            "dsig": (cen / w).tolist(),
            "rel_hi": hi.tolist(), "rel_lo": lo.tolist(),
            "rel_hi_nohalf": hin.tolist(), "rel_lo_nohalf": lon.tolist(),
        }

    d = f"{BASE}/results_nu" if CURRENT == "nu" else f"{BASE}/results"
    suf = ("_charm" if FLAVOUR == "charm" else "") + (
        "" if PTO == 1 else f"_pto{PTO}")
    tag = beams.Beams("mu", e_gev).tag
    esuf = "" if e_gev == beams.ANCHOR_ENERGY else f"_{tag}"
    # region, target, TMC: all "" at the defaults, so earlier files keep their names
    suf += (("" if sel.name == selection.DEFAULT else f"_{sel.name}")
            + target.suffix() + target.tmc_suffix())
    fn = f"{d}/mhou_diff_{SCHEME}{suf}{esuf}{selection.q2_suffix()}.json"
    with open(fn, "w") as f:
        json.dump(out, f)
    print(f"wrote {fn}")


if __name__ == "__main__":
    main()
