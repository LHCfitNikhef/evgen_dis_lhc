#!/usr/bin/env python3
"""The 7-point scale band on the inclusive fiducial cross-section, in FONLL.

WHY THIS EXISTS.  The paper figure comparing the NLO matchings needs a
yardstick: the matchings differ from the analytic reference by about 1%, and
whether that is a discrepancy or simply the size of the missing higher orders
cannot be read off the figure without the scale band drawn on it.

The existing `mhou_charm_fraction.json` does not serve: it is ZM-VFNS and it
varies the CHARM FRACTION.  What the figure needs is sigma_fid itself, in the
same FONLL scheme the figure's reference curve uses -- a band computed in a
different scheme from the curve it is drawn around would be meaningless.

THE SEVEN POINTS are the standard ones: each of xiR, xiF in {1, 2, 1/2} with
the two opposite extremes (2, 1/2) and (1/2, 2) excluded.

WHY IT IS CHEAP NOW.  A PineAPPL grid takes the scale factors directly
(`Grid.convolve(..., xi=[(xiR, xiF)])`), so all seven points come from the
same stored grids rather than from seven YADISM runs.  FONLL needs three grids
per point -- ZM, FFNS, FFN0 -- and every one of them is already built and
cached for the PDF study.

Usage:
  analysis/mhou_sigma_fonll.py [--current mu|nu] [--zm] [--pto 2]
                               [--pdf SET] [--energies 400,700,1000,2000,4000]
Writes results{,_nu}/mhou_sigma_fonll{,_mu}[_<sel>][_<target>][_tmc]
       [_q2min<N>][_pto<N>][_zm][_<pdf>].json

The region, the target and the TMC mode come from $BENCH_SELECTION,
$BENCH_TARGET and $BENCH_TMC, like every calculator, and each is in the name
when it is not the default -- e.g. the final-region tungsten TMC-on NNLO ZM
band on the neutrino side is results_nu/mhou_sigma_fonll_q4w3_W_tmc_pto2_zm.json.
"""
import importlib
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CURRENT = (sys.argv[sys.argv.index("--current") + 1]
           if "--current" in sys.argv else "nu")
# --energies: the ladder has five points (beams.BENCH_ENERGIES), the earlier scan
# three; the default stays the earlier scan so every existing band is re-derived
# exactly as it was.
ENERGIES = ([float(e) for e in
             sys.argv[sys.argv.index("--energies") + 1].split(",")]
            if "--energies" in sys.argv else None)
# --pdf SETTLES WHETHER THE xiF GROWTH IS PHYSICS OR THE GRID EDGE.
# NNPDF4.0 stops at Q2 = 2.722, so at xiF = 1/2 the whole muon fiducial region
# is being extrapolated.  MSHT20 reaches Q2 = 1, which keeps muF^2 = Q^2/4
# ON-GRID for every Q2 >= 4 -- the entire region -- so the floor never
# engages and the band is a scale variation and nothing else.  If the growth
# with energy survives that, it is physical; if it collapses, it was the
# boundary.
PDF_SET = (sys.argv[sys.argv.index("--pdf") + 1]
           if "--pdf" in sys.argv else None)
# --zm computes the band in ZM-VFNS instead of FONLL.  Needed because a band
# must be drawn around a reference computed in the SAME scheme: paper plot 1
# compares massless-charm generators, so its reference is ZM-VFNS and its
# band has to be too (user, 2026-08-29).
ZM_ONLY = "--zm" in sys.argv
# --pto SELECTS THE PERTURBATIVE ORDER (1 = NLO, the default; 2 = NNLO).
# The point of running both is that a scale band MUST shrink as the order
# rises -- if it does not, the variation is not estimating what it claims to.
PTO = int(sys.argv[sys.argv.index("--pto") + 1]) if "--pto" in sys.argv else 1

# NO NEUTRINO DIS AT NNLO IN FONLL (user, 2026-08-29).  The general-mass
# construction needs massive charged-current coefficient functions at second
# order, which this benchmark does not have -- the neutrino side has an NNLO
# result in ZM-VFNS alone, and `results_nu/histos_yadism_nnlo_fonll_damp.json`
# has never existed while the muon one has.
#
# THE REASON THIS IS A HARD REFUSAL RATHER THAN A NOTE: asked for it anyway,
# the chain does not fail.  It builds grids, convolves them and prints a
# perfectly ordinary band -- +1.26/-0.98% at 400 GeV on the run that prompted
# this guard -- which is a number for a calculation that was never performed.
# A silent plausible answer is the worst outcome available here.
if PTO >= 2 and CURRENT == "nu" and not ZM_ONLY:
    sys.exit("no neutrino DIS at NNLO in FONLL: the massive charged-current\n"
             "  coefficient functions are not available at this order.\n"
             "  Use --zm for the massless NNLO band, or --current mu.")
if CURRENT not in ("mu", "nu"):
    sys.exit("--current takes mu or nu")
os.environ["BENCH_CURRENT"] = CURRENT

# the standard seven, with the two opposite extremes dropped
SCALES = [(1.0, 1.0), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5),
          (2.0, 2.0), (0.5, 0.5)]

# THE FACTORISATION SCALE IS FLOORED AT THE PDF'S OWN GRID MINIMUM, and
# without this the band is not a scale variation at all.
#
# The fiducial region starts at Q2 = 4 GeV^2, so xiF = 1/2 asks for
# muF^2 = 1 GeV^2 -- far below NNPDF4.0's grid, which stops at 2.72.  LHAPDF
# does not refuse: it extrapolates, and the extrapolation collapses.  Measured
# on F2 at x = 0.01, the xiF = 1/2 shift runs
#
#     Q2 = 100 (muF2 = 25, on-grid)    -4.45%
#     Q2 =  30 (muF2 = 7.5, on-grid)   -6.12%
#     Q2 =  10 (muF2 = 2.5, OFF-grid) -14.96%
#     Q2 =   4 (muF2 = 1.0, OFF-grid) -44.96%
#
# -- a clean jump at the boundary, not a physical trend.  Integrated over a
# fiducial region the photon propagator pins at low Q2, that gave a muon NC
# "band" of -21% that GREW with beam energy, which is backwards: more energy
# reaches higher Q2 and the missing higher orders shrink.  The growth was the
# low-x extrapolation getting worse, not physics.
#
# The fix is the standard frozen-scale prescription, applied PER Q2 NODE so
# only the nodes that would leave the grid are affected:
#
#     muF^2 = max(xiF^2 Q^2, Q2min)      and      muR^2 = max(xiR^2 Q^2, Q2min)
#
# BOTH scales, not only muF (user, 2026-08-28: "never evaluate the DIS
# structure functions below Q0^2, also when varying the scales").  Flooring
# muF alone still left the muon band GROWING with beam energy, because
# alpha_s(muR) at muR = Q/2 = 1 GeV is evaluated below the coupling's own
# reference scale and towards the Landau pole -- the same class of error one
# step along.  It is a practical modification of the 7-point prescription and
# is recorded in the output rather than left implicit.
FLOOR_SCALES = True


# Historic --pdf filename tags, kept so the file on disk keeps its name.
_LEGACY_PDF_TAGS = {"MSHT20nnlo_as118": "msht20nnloas"}


def pdf_tag(setname):
    """Filename tag for --pdf: the WHOLE compacted set name.

    It was the first twelve characters, which cannot tell
    NNPDF40_nnlo_as_01180, its isospin mirror _n and the tungsten average
    _W184free apart -- all three are "nnpdf40nnloa" -- so a band for one
    would have overwritten a band for another under a name that fits both.
    """
    return _LEGACY_PDF_TAGS.get(setname,
                                setname.replace("_", "").lower())


def main():
    import beams
    import lhapdf
    lhapdf.setVerbosity(0)
    import pineappl_grids as pg
    import selection
    import target
    pdd = importlib.import_module("pdf_dependence")

    # THE TARGET PICKS THE SET, as in every calculator (target.pdfset() is
    # NNPDF40_nnlo_as_01180 for the default proton, i.e. pdd.REFERENCE).
    # --pdf with a non-proton target is refused: the name would carry both
    # and the calculation only one of them.
    if PDF_SET and target.current() != target.DEFAULT:
        sys.exit("--pdf and BENCH_TARGET are exclusive: the target already "
                 f"selects the set ({target.pdfset()})")
    setname = PDF_SET or target.pdfset()
    pdf = lhapdf.mkPDF(setname, 0)
    print(f"[mhou] PDF {setname}, grid Q2min = {pdf.q2Min:.4g} GeV^2; "
          f"at xiF=1/2 that is on-grid for Q2 >= {4*pdf.q2Min:.3g}",
          flush=True)
    out = {"what": (f"7-point muR/muF band on sigma_fid, {pdd.TITLE_CUR} DIS, "
                    f"YADISM NLO in FONLL"),
           "q2_min": float(os.environ.get("BENCH_Q2MIN", 4.0)),
           "scheme": ("ZM-VFNS" if ZM_ONLY else
                      "FONLL general-mass, damped; charm only"),
           "scale_floor": ("muF^2 = max(xiF^2 Q^2, Q2min_PDF) and "
                           "muR^2 = max(xiR^2 Q^2, Q2min_PDF); without it the "
                           "downward points leave the PDF grid and the alphas "
                           "reference scale, and the band is an extrapolation "
                           "artefact rather than a scale variation"),
           "pdf": setname, "q2min_pdf": pdf.q2Min, "current": CURRENT,
           "scales": [{"xiR": a, "xiF": b} for a, b in SCALES],
           "points": {}}
    # Region, target and TMC, stamped whenever any is away from its default
    # (absent otherwise, so a earlier re-run is byte-identical).
    sel = selection.get()
    # (the the earlier production --pdf runs already name their set in "pdf" and stay as they were)
    if (target.current() != target.DEFAULT or target.tmc() != 0
            or sel.name != selection.DEFAULT):
        out["selection"] = sel.as_dict()
        out["pdf_set"] = setname
        out["target"] = target.current()
        out["tmc"] = target.tmc()
    for e_gev in (ENERGIES or beams.ENERGIES):
        os.environ["BENCH_ENERGY"] = f"{e_gev:g}"
        for m in ("yadism_calc", "yadism_cc_calc", "yadism_cc_charm_calc",
                  "yadism_charm_calc", "analyze", "beams", "runmeta",
                  "selection"):
            sys.modules.pop(m, None)
        from yadism_calc import sf_grid_nodes
        x_nodes, q2_nodes = sf_grid_nodes()
        theory, obs, names = pdd.observables_card(PTO, x_nodes, q2_nodes)
        keys = {}
        parts = pdd.FONLL_PARTS[:1] if ZM_ONLY else pdd.FONLL_PARTS
        for tag, fns, nfff in parts:
            want = names if tag == "zm" else pdd.CHARM_NAMES
            keys[tag], _x, _q = pg.build(CURRENT, e_gev, PTO, want,
                                         x_nodes=x_nodes, q2_nodes=q2_nodes,
                                         fns=fns, nfff=nfff)
        import numpy as np
        nx, q2a = len(x_nodes), np.asarray(q2_nodes, dtype=float)
        q2min = pdf.q2Min

        def res_floored(key, want, xir, xif):
            """res_like with BOTH scales floored at the PDF's grid minimum.

            Nodes needing the same effective (xiR, xiF) pair are grouped, so
            this costs one convolution per DISTINCT pair rather than one per
            node -- typically two or three.
            """
            if not FLOOR_SCALES:
                return pg.res_like(key, want, pdf, xir, xif)
            floor = np.sqrt(q2min / q2a)
            eff_r = np.round(np.maximum(xir, floor), 12)
            eff_f = np.round(np.maximum(xif, floor), 12)
            out = {n: [None] * (len(q2a) * nx) for n in want}
            for pair in sorted(set(zip(eff_r.tolist(), eff_f.tolist()))):
                r = pg.res_like(key, want, pdf, float(pair[0]),
                                float(pair[1]))
                sel = np.where((eff_r == pair[0]) & (eff_f == pair[1]))[0]
                for n in want:
                    for iq in sel:
                        lo = int(iq) * nx
                        out[n][lo:lo + nx] = r[n][lo:lo + nx]
            return out

        rows = []
        for xir, xif in SCALES:
            per = {}
            for tag, _f, _n in parts:
                want = names if tag == "zm" else pdd.CHARM_NAMES
                per[tag] = res_floored(keys[tag], want, xir, xif)
            res = (per["zm"] if ZM_ONLY else
                   pdd.fonll_combine(per["zm"], per["ffns"], per["ffn0"],
                                     names, q2_nodes, len(x_nodes)))
            s_ch, s_in, _f = pdd.fraction_for(res, names, x_nodes, q2_nodes)
            # CHARM TOO: paper plot 6 needs a band on the charm cross-section,
            # and it must come from the same run as the inclusive one so the
            # two are the same scale points on the same grids.
            rows.append({"xiR": xir, "xiF": xif, "sigma_fid_pb": s_in,
                         "sigma_charm_pb": s_ch})
            print(f"[mhou] {e_gev:g} GeV  ({xir:g}, {xif:g})  "
                  f"sigma = {s_in:.6g} pb", flush=True)
        cen = [r["sigma_fid_pb"] for r in rows
               if r["xiR"] == 1.0 and r["xiF"] == 1.0][0]
        vals = [r["sigma_fid_pb"] for r in rows]
        cch = [r["sigma_charm_pb"] for r in rows
               if r["xiR"] == 1.0 and r["xiF"] == 1.0][0]
        vch = [r["sigma_charm_pb"] for r in rows]
        out["points"][f"{e_gev:g}"] = {
            "rows": rows, "central_pb": cen,
            "rel_hi": max(vals) / cen - 1.0,
            "rel_lo": min(vals) / cen - 1.0,
            "charm_central_pb": cch,
            "charm_rel_hi": max(vch) / cch - 1.0,
            "charm_rel_lo": min(vch) / cch - 1.0}
        print(f"[mhou] {e_gev:g} GeV  band "
              f"+{100*(max(vals)/cen-1):.2f}% / "
              f"{100*(min(vals)/cen-1):.2f}%", flush=True)
    d = f"{BASE}/results_nu" if CURRENT == "nu" else f"{BASE}/results"
    suf = "" if CURRENT == "nu" else "_mu"
    # >>> THE REGION, THE TARGET AND THE TMC MODE GO RIGHT AFTER THE CURRENT
    # (2026-09-13). <<<  None of them was in the name, so a q4w3 tungsten run
    # would have written its band straight over mhou_sigma_fonll_mu_zm.json,
    # the band pp01 draws around the inclusive proton reference.  All three
    # are "" at the defaults, so every existing file keeps its name.
    suf += (("" if sel.name == selection.DEFAULT else f"_{sel.name}")
            + target.suffix() + target.tmc_suffix())
    q2min = os.environ.get("BENCH_Q2MIN")
    if q2min and float(q2min) != 4.0:
        suf += f"_q2min{float(q2min):g}"
    if PTO != 1:
        suf += f"_pto{PTO}"
    if ZM_ONLY:
        suf += "_zm"
    if PDF_SET:
        suf += "_" + pdf_tag(PDF_SET)
    fn = f"{d}/mhou_sigma_fonll{suf}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {fn}")


if __name__ == "__main__":
    main()
