#!/usr/bin/env python3
"""YADISM (FONLL) against the NOMAD dimuon measurement.

The benchmark's first comparison against DATA rather than against another
calculation (user, 2026-08-30).  Everything else on the page compares
generators with an analytic reference or with one another; this asks whether
the analytic reference reproduces a measurement, in the one place where a
charged-current charm measurement exists.

WHAT IS MEASURED.  NOMAD reports

    R(E_nu) = sigma_mumu / sigma_CC ,

the fraction of charged-current interactions producing a second, opposite-sign
muon, binned in reconstructed neutrino energy over 6-300 GeV -- see
data/nomad/README.md.  A second muon in charged-current DIS comes from the
semileptonic decay of a charmed hadron, so R is a charm-production fraction
times a branching ratio.

WHAT IS COMPUTED HERE.

    R_th(E) = B(E) * sigma_charm(E) / sigma_incl(E)

with both cross-sections from YADISM in the FONLL general-mass scheme at NLO,
integrated over the same region the measurement's own theory used:

    Q2 > Q2_min ,   x from Q2_min/(2 M E) to 1 ,   y <= 1  (no y cut),

and B(E) = a/(1 + b/E), a = 0.097, b = 6.7, the effective branching ratio
fitted by NOMAD.

>>> B(E) IS NOT THE SEMILEPTONIC BRANCHING RATIO. <<<  It is an EFFECTIVE
branching ratio, and its energy dependence is where the second-muon
ACCEPTANCE lives -- the threshold on the second muon's momentum, the
reconstruction efficiency, the charm hadron mix.  The consequence matters for
anything beyond this script: a Monte Carlo that decays its own charm and
applies a muon momentum cut must NOT also multiply by B(E), or the acceptance
is counted twice.  Here there is no cascade decay to speak of, so B(E) carries
all of it and is applied once.

TWO APPROXIMATIONS, EACH DELIBERATE AND EACH STATED IN THE OUTPUT.
The target is NOT one of them -- see TARGET below: it is iron, as the
measurement's own theory used, built from free-proton parton distributions.
"No nuclear corrections" in that analysis means no nuclear PDFs, not a free
proton, and the difference between the two readings is a factor of two.

1. Q2 > 2.7225 GeV2, THE PDF'S OWN GRID MINIMUM, where NOMAD reconstructs
   from Q2 > 1 GeV2.  The benchmark's standing rule is never to evaluate a
   structure function below Q0^2, so the lower value is not available to us
   with this parton distribution.  It cancels to a large extent in the RATIO
   -- the cut removes the same low-Q2 events from numerator and denominator --
   and that cancellation is MEASURED here rather than assumed: the whole
   comparison is repeated at Q2_min = 4 and 6 GeV2 and the spread recorded.

2. THE BIN CENTRE, not an integral over the bin.  The measurement's own
   theory is evaluated at the bin centre of column 3, so doing the same keeps
   this like-for-like with the published curve.  It is an approximation the
   published calculation also makes, and it is worst in the first bin
   (6-22 GeV, represented by 15.91) where R is steepest.

NOT DONE, AND NOT PRETENDED: no PDF uncertainty band.  The structure-function
operator is PDF-independent, so a band is affordable, but it is a member loop
and this is the first pass.  The scale band IS computed, because the same
operator gives it for the cost of a convolution.

WHY THE STRUCTURE FUNCTIONS ARE COMPUTED ONCE FOR ALL NINETEEN ENERGIES.
F2, FL and xF3 do not depend on the beam energy -- only the region integrated
over does.  So one grid spanning the union of the nineteen regions serves them
all, which turns nineteen YADISM runs into one.  Getting this wrong is not
subtle, it is just slow.

Usage:  analysis/nomad_dimuon.py
        -> results_nu/nomad_dimuon_fonll.json
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.interpolate import RectBivariateSpline

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

DATA = f"{BASE}/data/nomad/nomad_table6_enu.txt"
OUT = f"{BASE}/results_nu/nomad_dimuon_fonll.json"

# The NOMAD effective branching ratio, from the experiment's paper by way of
# the reference implementation (modules/branching.f).  See the docstring: the
# energy dependence is the acceptance.
BR_A, BR_B = 0.097, 6.7

M_P = 0.938272
PDFSET = "NNPDF40_nnlo_as_01180"
PTO = 1                                  # NLO

# >>> THE TARGET IS IRON, BUILT FROM FREE-PROTON PARTON DISTRIBUTIONS. <<<
# This is what the measurement's own theory did, and what "no nuclear
# corrections" means in that analysis: no nuclear PDFs, but the right
# nucleon.  It is not a detail.  For a NEUTRINO charged current the inclusive
# rate on a proton is driven by the valence d, while on a near-isoscalar
# nucleus it is driven by the average of u and d -- and u_v is about twice
# d_v, so the proton denominator is far too small.  Charm production, being
# mostly s -> c off the sea, barely moves.  The charm FRACTION on a free
# proton therefore comes out roughly twice the iron one: the first run of
# this script, on a proton, sat a factor 1.8 to 3.0 above the data with the
# right shape, and it was the target rather than the physics.
TARGET = "iron"

# Q2 floors to compare.  The first is the nominal -- the parton
# distribution's own grid minimum, which is as low as the standing rule
# allows -- and the other two exist to measure how much of the answer the
# choice is responsible for.
Q2MINS = [2.7225, 4.0, 6.0]

# The 7-point envelope, on the same operator, for the cost of a convolution.
SCALES = [(1.0, 1.0), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5),
          (2.0, 2.0), (0.5, 0.5)]

XMAX = 0.999


def load_data():
    """The 19 NOMAD points: (E_lo, E_hi, E_centre, R, stat, sys), R in 1e-3."""
    rows = []
    with open(DATA) as f:
        for line in f:
            c = line.split()
            if len(c) < 6:
                continue
            rows.append({"e_lo": float(c[0]), "e_hi": float(c[1]),
                         "e": float(c[2]), "r": float(c[3]) * 1e-3,
                         "stat": float(c[4]) * 1e-3,
                         "sys": float(c[5]) * 1e-3})
    return rows


def br(e):
    """NOMAD's effective dimuon branching ratio at neutrino energy e."""
    return BR_A / (1.0 + BR_B / e)


# ----------------------------------------------------------- the SF grid --
def grid_nodes(e_max, q2_min_lowest):
    """(x, Q2) nodes spanning the union of all nineteen fiducial regions.

    Q2 reaches 2 M E x <= 2 M E_max at x = 1; x reaches down to
    q2_min/(2 M E_max).  A margin on each side keeps the cubic spline off its
    own boundary, where a RectBivariateSpline extrapolates rather than
    refusing.
    """
    two_kp_max = 2.0 * M_P * e_max
    x_lo = 0.5 * q2_min_lowest / two_kp_max
    q2_hi = 1.3 * two_kp_max
    return (np.geomspace(x_lo, XMAX, 60).tolist(),
            np.geomspace(0.9 * q2_min_lowest, q2_hi, 40).tolist())


def splines_from(res, names, x_nodes, q2_nodes):
    """Cubic splines in (ln x, ln Q2), one per structure function."""
    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    out = {}
    for n in names:
        v = np.array([p["result"] for p in res[n]])
        out[n] = RectBivariateSpline(
            lx, lq, v.reshape(len(q2_nodes), len(x_nodes)).T, kx=3, ky=3)
    return out


# ------------------------------------------------------------ integration --
def _gl(a, b, order, npanel, log):
    """Gauss-Legendre nodes and weights on [a, b], npanel panels."""
    t, w = np.polynomial.legendre.leggauss(order)
    edges = (np.geomspace(a, b, npanel + 1) if log
             else np.linspace(a, b, npanel + 1))
    xs, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if log:
            la, lb = math.log(lo), math.log(hi)
            n = np.exp(0.5 * (la + lb) + 0.5 * (lb - la) * t)
            xs.append(n)
            ws.append(0.5 * (lb - la) * w * n)
        else:
            xs.append(0.5 * (lo + hi) + 0.5 * (hi - lo) * t)
            ws.append(0.5 * (hi - lo) * w)
    return np.concatenate(xs), np.concatenate(ws)


GF = 1.1663787e-05
MW = 80.398
MW2 = MW * MW
GEV2_PB = 3.893793e8


def sigma(sp, tag, e_nu, q2_min):
    """Integrated cross-section in pb over the NOMAD region at one energy.

    tag is "total" or "charm"; the two differ only in which structure
    functions are read, so numerator and denominator cannot drift apart.

        d2sigma/dx dQ2 = GF^2/(4 pi x) (MW^2/(MW^2+Q2))^2
                         [ Ypc F2 - y^2 FL + Ym xF3 ]

    with Ypc = 1 + (1-y)^2 - 2 (M x y)^2 / Q2, the exact coefficient rather
    than the massless Y_+.  It reaches 3.6% at the bottom of the region here,
    where the benchmark's own kinematics never goes, and it very nearly
    cancels in the ratio -- but the reference implementation carries it and
    dropping it would be a difference between the two for no reason.
    """
    two_kp = 2.0 * M_P * e_nu
    x_lo = q2_min / two_kp
    if x_lo >= XMAX:
        return 0.0
    xs, xw = _gl(x_lo, XMAX, 12, 8, log=True)
    total = 0.0
    for x, wx in zip(xs, xw):
        q2_hi = two_kp * x                      # y <= 1
        if q2_hi <= q2_min:
            continue
        q2s, q2w = _gl(q2_min, q2_hi, 10, 4, log=True)
        y = q2s / (two_kp * x)
        ypc = 1.0 + (1.0 - y) ** 2 - 2.0 * (M_P * x * y) ** 2 / q2s
        ym = 1.0 - (1.0 - y) ** 2
        lx, lq = np.log(x), np.log(q2s)
        f2 = sp[f"F2_{tag}"].ev(lx, lq)
        fl = sp[f"FL_{tag}"].ev(lx, lq)
        f3 = sp[f"F3_{tag}"].ev(lx, lq)
        comb = ypc * f2 - y * y * fl + ym * f3
        prop = (MW2 / (MW2 + q2s)) ** 2
        integ = GF * GF / (4.0 * math.pi * x) * prop * comb * GEV2_PB
        total += wx * float(np.sum(q2w * integ))
    return total


def r_theory(sp, e_nu, q2_min):
    """B(E) sigma_charm / sigma_incl, and the two cross-sections."""
    sc = sigma(sp, "charm", e_nu, q2_min)
    si = sigma(sp, "total", e_nu, q2_min)
    return (br(e_nu) * sc / si if si > 0 else float("nan")), sc, si


def chi2(rows, key):
    """chi2 with stat and sys added in quadrature -- no covariance exists."""
    return sum((r[key] - r["r"]) ** 2 / (r["stat"] ** 2 + r["sys"] ** 2)
               for r in rows if np.isfinite(r[key]))


# ------------------------------------------------------------------ main --
def main():
    import lhapdf
    lhapdf.setVerbosity(0)
    import pdf_dependence as pdd
    # the EW coupling the theory card carries; it lives in the calculator, not
    # in pdf_dependence, which only re-exports what it happens to use
    from yadism_cc_calc import ALPHA

    # NO PINEAPPL GRIDS HERE, and the reason is a cache collision waiting to
    # happen: pineappl_grids.grid_key() hashes the current, the energy, the
    # order, the nodes and the flavour scheme -- but NOT the target.  A grid
    # built for iron would be handed back to a proton request, or the other
    # way round, with no error and a plausible number.  The direct YADISM
    # path costs a few minutes more and cannot do that.
    use_grids = False
    rows = load_data()
    e_max = max(r["e"] for r in rows)
    x_nodes, q2_nodes = grid_nodes(e_max, min(Q2MINS))

    pdf = lhapdf.mkPDF(PDFSET, 0)
    if min(Q2MINS) < pdf.q2Min - 1e-9:
        sys.exit(f"Q2_min {min(Q2MINS)} is below the PDF grid minimum "
                 f"{pdf.q2Min:.4g}: the standing rule is never to evaluate a "
                 f"structure function below Q0^2.")
    print(f"[nomad] {len(rows)} points, E = {rows[0]['e']:.2f} - "
          f"{e_max:.1f} GeV; grid x >= {x_nodes[0]:.2e}, "
          f"Q2 <= {q2_nodes[-1]:.0f} GeV2", flush=True)

    names = ["F2_charm", "FL_charm", "F3_charm",
             "F2_total", "FL_total", "F3_total"]
    t0 = time.time()
    keys, outs = {}, {}
    if use_grids:
        import pineappl_grids as pg
        for tag, fns, nfff in pdd.FONLL_PARTS:
            want = names if tag == "zm" else pdd.CHARM_NAMES
            keys[tag], _x, _q = pg.build("nu", 1000.0, PTO, want,
                                         x_nodes=x_nodes, q2_nodes=q2_nodes,
                                         fns=fns, nfff=nfff,
                                         force="--force" in sys.argv)
    else:
        import yadism
        for tag, fns, nfff in pdd.FONLL_PARTS:
            th, ob = pdd.observables_card(PTO, x_nodes, q2_nodes)[:2]
            th["FNS"], th["NfFF"] = fns, nfff
            ob["TargetDIS"] = TARGET
            if tag != "zm":
                ob["observables"] = {n: ob["observables"][n]
                                     for n in pdd.CHARM_NAMES}
            outs[tag] = yadism.run_yadism(th, ob)
    print(f"[nomad] structure functions ready in {time.time()-t0:.0f} s",
          flush=True)

    q2a = np.asarray(q2_nodes, dtype=float)
    nx = len(x_nodes)
    q2min_pdf = pdf.q2Min

    def convolved(xir, xif):
        """FONLL structure functions at one scale choice, scales FLOORED.

        The floor is the standing rule (see analysis/mhou_sigma_fonll.py):
        muF^2 = max(xiF^2 Q^2, Q0^2) and the same for muR.  Without it the
        downward variations leave the parton distribution's grid and the
        band becomes an extrapolation artefact -- which once faked a -45%
        collapse of F2 at the bottom of a fiducial region.
        """
        floor = np.sqrt(q2min_pdf / q2a)
        eff_r = np.round(np.maximum(xir, floor), 12)
        eff_f = np.round(np.maximum(xif, floor), 12)
        per = {}
        for tag, _fns, _nfff in pdd.FONLL_PARTS:
            want = names if tag == "zm" else pdd.CHARM_NAMES
            acc = {n: [None] * (len(q2a) * nx) for n in want}
            for pr, pf in sorted(set(zip(eff_r.tolist(), eff_f.tolist()))):
                if use_grids:
                    r = pg.res_like(keys[tag], want, pdf, float(pr), float(pf))
                else:
                    r = outs[tag].apply_pdf_alphas_alphaqed_xir_xif(
                        pdf, lambda muR, _p=pdf: _p.alphasQ(muR),
                        lambda _muR: ALPHA, float(pr), float(pf))
                sel = np.where((eff_r == pr) & (eff_f == pf))[0]
                for n in want:
                    for iq in sel:
                        for ix in range(nx):
                            acc[n][iq * nx + ix] = r[n][iq * nx + ix]
            per[tag] = acc
        res = pdd.fonll_combine(per["zm"], per["ffns"], per["ffn0"],
                                names, q2_nodes, nx)
        sp_fonll = splines_from(res, names, x_nodes, q2_nodes)
        # >>> THE SAME DENOMINATOR, THREE DIFFERENT NUMERATORS. <<<  The
        # question this comparison kept failing on is which charm structure
        # function belongs at NOMAD's kinematics, so the three candidates are
        # built side by side and the ratio computed with each.  Only the
        # NUMERATOR changes: the inclusive cross-section is light-quark
        # dominated and its scheme moves it by per cent, while the charm
        # piece moves by a factor -- so holding the denominator fixed
        # isolates exactly the thing in question.
        #   fonll  ZM + damped (FFNS - FFN0), the benchmark's convention
        #   zm     massless charm, no threshold suppression at all
        #   ffns   nf = 3 MASSIVE charm, full threshold suppression
        alt = {}
        for tag in ("zm", "ffns"):
            r2 = dict(res)
            for n in pdd.CHARM_NAMES:
                r2[n] = per[tag][n]
            alt[tag] = splines_from(r2, names, x_nodes, q2_nodes)
        return sp_fonll, alt["zm"], alt["ffns"]

    if use_grids:
        import pineappl_grids as pg                          # noqa: F811

    # --- central, at each Q2 floor ---------------------------------------
    sp, sp_zm, sp_ffns = convolved(1.0, 1.0)
    # THE SCHEME SCAN, which is the point of this pass: the same ratio with a
    # massless charm numerator and with a massive one.
    for r in rows:
        r["r_zm"] = r_theory(sp_zm, r["e"], Q2MINS[0])[0]
        r["r_ffns"] = r_theory(sp_ffns, r["e"], Q2MINS[0])[0]
    for q2m in Q2MINS:
        k = f"r_q2min_{q2m:g}"
        for r in rows:
            r[k], sc, si = r_theory(sp, r["e"], q2m)
            if q2m == Q2MINS[0]:
                r["sigma_charm_pb"], r["sigma_incl_pb"] = sc, si
    for r in rows:
        r["r_th"] = r[f"r_q2min_{Q2MINS[0]:g}"]

    # --- the 7-point scale band, at the nominal floor ---------------------
    print("[nomad] 7-point scale envelope ...", flush=True)
    env = {r["e"]: [] for r in rows}
    for xir, xif in SCALES:
        spv = sp if (xir, xif) == (1.0, 1.0) else convolved(xir, xif)[0]
        for r in rows:
            env[r["e"]].append(r_theory(spv, r["e"], Q2MINS[0])[0])
    for r in rows:
        v = [x for x in env[r["e"]] if np.isfinite(x)]
        r["r_th_hi"], r["r_th_lo"] = max(v), min(v)

    out = {
        "what": ("YADISM NLO in FONLL against the NOMAD dimuon ratio "
                 "sigma_mumu/sigma_CC, binned in neutrino energy"),
        "scheme": "FONLL general-mass, damped",
        "pdf": PDFSET, "pto": PTO,
        "target": (f"{TARGET} nucleon built from free-proton parton "
                   f"distributions, as the measurement's own theory did; no "
                   f"nuclear PDFs, which that analysis argued cancel in the "
                   f"ratio. Using a free PROTON instead inflates the charm "
                   f"fraction by roughly a factor of two, because the "
                   f"inclusive neutrino CC rate on a proton is driven by the "
                   f"valence d while charm production is not"),
        "q2_min_nominal": Q2MINS[0],
        "q2_min_variants": Q2MINS,
        "q2_min_note": ("the PDF's own grid minimum; NOMAD reconstructs from "
                        "Q2 > 1 GeV2, which the standing rule against "
                        "evaluating below Q0^2 puts out of reach. The cut "
                        "largely cancels in the ratio, and the variants "
                        "measure how much"),
        "branching": {"form": "a / (1 + b/E)", "a": BR_A, "b": BR_B,
                      "note": ("NOMAD's EFFECTIVE branching ratio: its energy "
                               "dependence carries the second-muon "
                               "acceptance, so a Monte Carlo that applies its "
                               "own muon cut must not also apply this")},
        "covariance": ("none published; stat and sys added in quadrature and "
                       "the points treated as uncorrelated, which is an "
                       "approximation"),
        "bin_treatment": ("theory evaluated at the bin centre, as the "
                          "measurement's own theory does"),
        "scales": [{"xiR": a, "xiF": b} for a, b in SCALES],
        "scale_floor": ("muF^2 = max(xiF^2 Q^2, Q0^2), muR likewise"),
        "points": rows,
    }
    out["chi2"] = chi2(rows, "r_th")
    out["ndat"] = len(rows)
    out["chi2_per_point"] = out["chi2"] / len(rows)
    for q2m in Q2MINS:
        out[f"chi2_q2min_{q2m:g}"] = chi2(rows, f"r_q2min_{q2m:g}")
    for k in ("r_zm", "r_ffns"):
        out[f"chi2_{k}"] = chi2(rows, k)
    out["scheme_scan"] = (
        "the same ratio with three charm numerators and ONE denominator: "
        "r_th is FONLL (ZM + damped massive correction), r_zm is the massless "
        "charm alone, r_ffns is the nf=3 MASSIVE charm alone.  The FONLL "
        "damping factor is 0.03 at Q2 = 2.72 GeV2 and 0.19 at 4, so over "
        "NOMAD's kinematics FONLL is very nearly the massless result -- which "
        "is why the three are worth separating here and nowhere else in this "
        "benchmark")
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[nomad] chi2/N = {out['chi2_per_point']:.2f} over "
          f"{len(rows)} points")
    for q2m in Q2MINS:
        print(f"           Q2min = {q2m:g}: chi2/N = "
              f"{out[f'chi2_q2min_{q2m:g}']/len(rows):.2f}")
    for k, lab in (("r_zm", "massless charm"), ("r_ffns", "massive (nf=3)")):
        print(f"           {lab:16s}: chi2/N = "
              f"{out[f'chi2_{k}']/len(rows):.2f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
