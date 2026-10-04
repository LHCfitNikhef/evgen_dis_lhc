#!/usr/bin/env python3
"""Why GENIE HEDIS sits 5% below the NLO reference: it is the PDF.

THE QUESTION (user, 2026-09-21): "why HEDIS GENIE is 5% lower than the NLO
YADISM calculation.  Is this because of the use of different PDF, heavy quark
mass effects, some artifact, or some other reason?  The GENIE HEDIS
calculation is based on APFEL NLO and should match the YADISM NLO calculation
if the same PDF and theory settings are used."

THE ANSWER: THE SAME THEORY SETTINGS ARE NOT USED, AND THE ONE THAT MATTERS IS
THE PARTON DISTRIBUTION.  GENIE's GHE19_00a tune builds its structure
functions with APFEL from

    NNPDF31sx_nlo_as_0118_LHCb_nf_6

(read off the tune's CommonParam.xml, param_set HEDIS-SF), while this
benchmark's reference uses NNPDF40_nnlo_as_01180 throughout.  That set is a
small-x / high-energy-neutrino member of the NNPDF3.1sx family, and in the
FASER region -- where the mean x is about 0.21 -- its light-quark densities
are 6 to 12% below NNPDF4.0's.  Run the SAME YADISM calculation with GENIE's
own set and the reference falls onto GENIE.

The rest of the tune's settings are NOT the discrepancy, and this script
records that too: HEDISStrucFunc.cxx sets APFEL to perturbative order 1 (NLO,
as the user expected), takes alpha_s from the PDF set itself (0.1180, ours
too), and maps Scheme = "BGR" to APFEL's FONLL-B with the PDF's own pole
masses.  FONLL-B against our ZM-VFNS is worth 0.4% here, the target-mass
correction 0.01%, and the W mass difference (80.385 against 80.379) 0.03%.

WHAT THIS SCRIPT DOES.  It collects, into one checked JSON:

  * the delivered GENIE HEDIS cross-sections and the tracked YADISM NLO
    references, both in the final region on tungsten per nucleon;
  * the same YADISM calculation convolved with GENIE's own PDF, from
    `mhou_sigma_fonll.py --pdf <set>` (see tools/genie_hedis_pdf.sh, which is
    how those files are produced);
  * the proton-level decomposition of the PDF difference -- NNPDF4.0 at NNLO
    and NLO, NNPDF3.1 at NNLO, NNPDF3.1sx at NLO and NLO+NLLx -- which shows
    that the shift is the sx SET rather than the 3.1-to-4.0 vintage, the
    perturbative order of the fit, or the small-x resummation;
  * the light-quark ratio in x from LHAPDF directly, which is what the
    x-differential ratio of the two calculations follows.

Inputs: results_nu/histos_{genie,yadism_nlo*}_q4w3_W*.json (tracked) and
        results_nu/mhou_sigma_fonll_q4w3_tmc[_zm]_<pdf>.json
        (tools/genie_hedis_pdf.sh).
Output: results_nu/genie_hedis_pdf.json
Usage:  analysis/genie_hedis_pdf.py
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

RES = f"{BASE}/results_nu"
REGION = "q4w3"

# The tune's own set, read off $GENIE_DIR/config/GHE19_00a/CommonParam.xml.
HEDIS_PDF = "NNPDF31sx_nlo_as_0118_LHCb_nf_6"
BENCH_PDF = "NNPDF40_nnlo_as_01180"
# ... and their tungsten free-nucleon averages, (74 p + 110 n)/184, which is
# what the reference is computed on (tools/make_isoscalar_pdf.py).
W_SUFFIX = "_W184free"

# the energy ladder the figure uses
ENERGIES = [400, 700, 1000, 2000, 4000]
TAGS = {400: "_400GeV", 700: "_700GeV", 1000: "", 2000: "_2TeV", 4000: "_4TeV"}

# The proton-level decomposition: each set against the benchmark's own.
DECOMP = [
    (BENCH_PDF, "the benchmark's set"),
    ("NNPDF40_nlo_as_01180", "same fit, NLO instead of NNLO"),
    ("NNPDF31_nnlo_as_0118", "the previous NNPDF global fit"),
    (HEDIS_PDF, "GENIE HEDIS's set"),
    ("NNPDF31sx_nlonllx_as_0118_LHCb_nf_6", "the same with small-x resummation"),
]

# x values at which the quark ratio is recorded; the region's mean x is 0.21.
X_POINTS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
Q2_REF = 100.0


def _load(path):
    with open(path) as f:
        return json.load(f)


def _sigma(name, tag):
    p = f"{RES}/histos_{name}_{REGION}_W{tag}.json"
    if not os.path.exists(p):
        return None
    d = _load(p)
    return d["sigma_fid_pb"], d.get("sigma_fid_err_pb")


def _mhou(pdf, zm=True):
    """sigma_fid from a `mhou_sigma_fonll.py --pdf` run, per energy."""
    import importlib
    tag = importlib.import_module("mhou_sigma_fonll").pdf_tag(pdf)
    p = f"{RES}/mhou_sigma_fonll_{REGION}_tmc{'_zm' if zm else ''}_{tag}.json"
    if not os.path.exists(p):
        return {}
    d = _load(p)
    return {float(e): r["central_pb"] for e, r in d["points"].items()}


def quark_ratio():
    """(d + s + b) of GENIE's set over the benchmark's, on the tungsten free
    nucleon: the combination a neutrino charged current sees at Born level."""
    import target                                            # noqa: F401
    import lhapdf
    lhapdf.setVerbosity(0)
    out = {}
    pdfs = {k: lhapdf.mkPDF(k + W_SUFFIX, 0) for k in (BENCH_PDF, HEDIS_PDF)}

    def q(p, x):
        return (p.xfxQ2(1, x, Q2_REF) + p.xfxQ2(3, x, Q2_REF)
                + p.xfxQ2(5, x, Q2_REF))

    for x in X_POINTS:
        out[f"{x:g}"] = q(pdfs[HEDIS_PDF], x) / q(pdfs[BENCH_PDF], x)
    return out


def main():
    zm_hedis = _mhou(HEDIS_PDF + W_SUFFIX)
    fo_hedis = _mhou(HEDIS_PDF + W_SUFFIX, zm=False)
    fo_bench = _mhou(BENCH_PDF + W_SUFFIX, zm=False)
    zm_bench = _mhou(BENCH_PDF + W_SUFFIX)
    if not zm_hedis:
        sys.exit("run tools/genie_hedis_pdf.sh first: the YADISM reference "
                 "with GENIE's own PDF is not on disk")

    rows = []
    for e in ENERGIES:
        g = _sigma("genie", TAGS[e])
        y_zm = _sigma("yadism_nlo", "_tmc" + TAGS[e])
        y_fo = _sigma("yadism_nlo_fonll_damp", "_tmc" + TAGS[e])
        sx = zm_hedis.get(float(e))
        row = {"energy_gev": e,
               "genie_hedis_pb": g[0] if g else None,
               "genie_hedis_err_pb": g[1] if g else None,
               "yadism_nlo_zm_pb": y_zm[0] if y_zm else None,
               "yadism_nlo_fonll_pb": y_fo[0] if y_fo else None,
               "yadism_nlo_zm_hedis_pdf_pb": sx,
               # the curve pp03 draws: the SAME reference (FONLL, TMC) with
               # GENIE HEDIS's own parton densities
               "yadism_nlo_fonll_hedis_pdf_pb": fo_hedis.get(float(e))}
        # THE THREE NUMBERS THE ANSWER IS MADE OF.
        #  * what the figure shows: GENIE over the benchmark's own reference;
        #  * the PDF alone: the same YADISM calculation with GENIE's set;
        #  * what is left once BOTH the PDF and the flavour scheme are
        #    GENIE's.  The scheme factor is FONLL/ZM measured on the
        #    benchmark's own set -- a ratio of schemes is far less
        #    PDF-sensitive than either term, and the direct FONLL run with
        #    GENIE's set at the anchor energy is recorded below to show by
        #    how little (fonll_factorisation_check).
        row["ratio_published"] = (row["genie_hedis_pb"] / row["yadism_nlo_fonll_pb"]
                                  if row["genie_hedis_pb"] and row["yadism_nlo_fonll_pb"]
                                  else None)
        row["pdf_only"] = (sx / row["yadism_nlo_zm_pb"]
                           if sx and row["yadism_nlo_zm_pb"] else None)
        if row["yadism_nlo_fonll_pb"] and row["yadism_nlo_zm_pb"]:
            row["scheme_factor_fonll_over_zm"] = (row["yadism_nlo_fonll_pb"]
                                                  / row["yadism_nlo_zm_pb"])
        if row["yadism_nlo_fonll_hedis_pdf_pb"] and row["yadism_nlo_fonll_pb"]:
            row["pdf_only_fonll"] = (row["yadism_nlo_fonll_hedis_pdf_pb"]
                                     / row["yadism_nlo_fonll_pb"])
            row["genie_over_fonll_same_pdf"] = (
                row["genie_hedis_pb"] / row["yadism_nlo_fonll_hedis_pdf_pb"]
                if row["genie_hedis_pb"] else None)
        if sx and row.get("scheme_factor_fonll_over_zm"):
            lfl = sx * row["scheme_factor_fonll_over_zm"]
            row["yadism_like_for_like_pb"] = lfl
            row["genie_over_like_for_like"] = (row["genie_hedis_pb"] / lfl
                                               if row["genie_hedis_pb"] else None)
        rows.append(row)

    check = None
    if fo_hedis.get(1000.0) and zm_hedis.get(1000.0):
        direct = fo_hedis[1000.0]
        row = [r for r in rows if r["energy_gev"] == 1000][0]
        check = {"energy_gev": 1000,
                 "fonll_with_hedis_pdf_pb": direct,
                 "factorised_pb": row.get("yadism_like_for_like_pb"),
                 "ratio": (row["yadism_like_for_like_pb"] / direct
                           if row.get("yadism_like_for_like_pb") else None),
                 "genie_over_direct_fonll": (row["genie_hedis_pb"] / direct
                                             if row["genie_hedis_pb"] else None)}

    # Closure: the same driver on the benchmark's own set must reproduce the
    # tracked reference, or the comparison above is between two pipelines
    # rather than between two PDFs.
    closure = None
    if zm_bench.get(1000.0):
        t = _sigma("yadism_nlo", "_tmc")
        closure = {"mhou_pb": zm_bench[1000.0], "tracked_pb": t[0] if t else None,
                   "ratio": (zm_bench[1000.0] / t[0]) if t else None}

    decomp = []
    for setname, why in DECOMP:
        pts = _mhou(setname)
        decomp.append({"pdf": setname, "what": why,
                       "sigma_1tev_pb": pts.get(1000.0)})
    ref = decomp[0]["sigma_1tev_pb"]
    for d in decomp:
        d["ratio_to_benchmark"] = (d["sigma_1tev_pb"] / ref
                                   if d["sigma_1tev_pb"] and ref else None)

    out = {
        "what": ("GENIE HEDIS (GHE19_00a) against the YADISM NLO reference in "
                 "the final region, and the PDF that accounts for the "
                 "difference"),
        "region": REGION, "target": "tungsten per nucleon (74 p + 110 n)",
        "hedis_tune": "GHE19_00a_00_000",
        "hedis_pdf": HEDIS_PDF, "benchmark_pdf": BENCH_PDF,
        "hedis_settings": {
            "perturbative_order": "NLO (APFEL::SetPerturbativeOrder(1))",
            "scheme": "BGR -> APFEL FONLL-B, pole masses from the PDF set",
            "alphas": "taken from the PDF set at MZ (0.1180, as ours)",
            "MW_gev": 80.385, "benchmark_MW_gev": 80.379,
            "sf_grid": "400 x 400 nodes, x >= 1e-9, Q2 in [2.6896, 1e10]",
        },
        "rows": rows,
        "fonll_factorisation_check": check,
        "pipeline_closure": closure,
        "proton_decomposition_1tev": decomp,
        "quark_ratio_x": quark_ratio(),
        "quark_ratio_note": (f"(d + s + b) of {HEDIS_PDF} over {BENCH_PDF}, "
                             f"tungsten free nucleon, Q2 = {Q2_REF:g} GeV2"),
    }
    p = f"{RES}/genie_hedis_pdf.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {p}\n")
    print(f"{'E [GeV]':>8} {'GENIE':>9} {'YADISM':>9} {'shown':>8} "
          f"{'PDF only':>9} {'FONLL PDF':>10} {'GENIE/FONLL same PDF':>21}")
    for r in rows:
        print(f"{r['energy_gev']:8g} {r['genie_hedis_pb'] or 0:9.4f} "
              f"{r['yadism_nlo_fonll_pb'] or 0:9.4f} "
              f"{r.get('ratio_published') or 0:8.4f} "
              f"{r.get('pdf_only') or 0:9.4f} "
              f"{r.get('pdf_only_fonll') or 0:10.4f} "
              f"{r.get('genie_over_fonll_same_pdf') or 0:21.4f}")
    if check:
        print(f"\n  FONLL factorisation at 1 TeV: direct "
              f"{check['fonll_factorisation_check'] if False else check['fonll_with_hedis_pdf_pb']:.4f} pb, "
              f"factorised {check['factorised_pb']:.4f} pb "
              f"({check['ratio']:.5f}); GENIE/direct "
              f"{check['genie_over_direct_fonll']:.5f}")
    if closure:
        print(f"  pipeline closure on the benchmark set: "
              f"{closure['ratio']:.6f}")
    print()
    for d in decomp:
        print(f"  {d['pdf']:38s} {d['sigma_1tev_pb'] or 0:8.4f} pb  "
              f"{d['ratio_to_benchmark'] or 0:.4f}   {d['what']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
