#!/usr/bin/env python3
"""YADISM analytic LO CHARM-PRODUCTION reference for the muon-DIS benchmark.

*** What "charm production" means here ***

In the ZM-VFNS at LO, charm production in NC DIS is **charm-INITIATED
scattering**: gamma*/Z + c -> c off the charm PDF, i.e. a charm quark in BOTH
the initial and the final state of the hard process.  It is the charm-quark
piece of the structure functions,

    F2^c = A_c(Q2) * x (c + cbar),   xF3^c = B_c(Q2) * x (c - cbar),  FL^c = 0

with A_c/B_c the usual gamma + gamma-Z + Z coupling combinations for a charge
2/3, T3 = +1/2 quark.  It is **NOT** photon-gluon fusion gamma* g -> c cbar,
which is O(alpha_s) and belongs to charm-PAIR production (a different strand).

Everything else -- kinematics, fiducial region, EW parameters, PDF, scales,
binning, lab-frame maps, quadrature -- is imported unchanged from
yadism_calc.py, so this is literally the inclusive calculation with the
structure functions restricted to the charm flavour.

Two INDEPENDENT calculations are performed and compared:
  (a) yadism  F2_charm / FL_charm / F3_charm at PTO = 0, ZM-VFNS;
  (b) a direct LHAPDF parton-model integral over the charm PDF only.
Both are pushed through the SAME three Gauss-Legendre integrators as the
inclusive reference, so the check is bin-by-bin and not only on sigma_fid.

The same machinery serves the bottom-initiated channel (`bottom` argument),
which is needed to decompose the inclusive Pythia/YADISM deficit by channel.

PERTURBATIVE ORDER.  The LO number is charm-INITIATED scattering only.  From
O(alpha_s) on, yadism's `F2_charm` observable additionally carries the charm
share of the GLUON/singlet terms, i.e. photon-gluon fusion gamma* g -> c cbar,
in which charm appears in the FINAL STATE ONLY.  That is the correct reference
for an NLO/NNLO generator sample tagged on hard-process charm, because such a
sample contains both channels -- but it means **the NLO charm number is not
the LO one plus a correction to the same process; it is a wider process**.
Never compare an NLO charm sample against the LO charm reference.

The direct LHAPDF cross-check (b) is a LO parton-model integral by
construction, so it is only run -- and only meaningful -- at LO.  At NLO/NNLO
the validation falls back to the agreement of the three independent
quadratures.

Usage:  yadism_charm_calc.py [charm|bottom] [lo|nlo|nnlo]   (default charm lo)
        -> results/histos_yadism_charm.json      (LO, name kept for
           back-compatibility), _charm_nlo.json, _charm_nnlo.json
"""
import json
import math
import hashlib
import os
import sys

import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams as beams_mod  # noqa: E402
import selection  # noqa: E402  -- the one Q2 floor
from analyze import BINS, sel_suffix  # noqa: E402
from yadism_calc import (sf_grid_nodes,  # noqa: E402
                         ALPHA, BASE, EBEAM, GEV2_PB, MZ, PDFSET, S2TW,  # noqa: E402
                         TWO_KP, XMAX, integrate_q2y, integrate_xq2,
                         integrate_ytheta, make_cards, make_sigma_red)
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp

FLAVOURS = {"charm": 4, "bottom": 5}
ORDERS = {"lo": 0, "nlo": 1, "nnlo": 2}




# --------------------------------------------------------------- (a) yadism
def _flavour_grid(pto, pdf, flavour, x_nodes, q2_nodes, fns, nfff):
    """One yadism run: F2/FL/F3 for a single flavour, as raw [x, Q2] grids."""
    import yadism

    theory, observables = make_cards(pto, x_nodes, q2_nodes)
    theory["FNS"] = fns
    theory["NfFF"] = nfff
    kins = observables["observables"]["F2_total"]
    names = [f"{k}_{flavour}" for k in ("F2", "FL", "F3")]
    observables["observables"] = {n: kins for n in names}
    out = yadism.run_yadism(theory, observables)
    res = out.apply_pdf_alphas_alphaqed_xir_xif(
        pdf, lambda muR: pdf.alphasQ(muR), lambda _muR: ALPHA, 1.0, 1.0)
    grids = {}
    for name in names:
        vals = np.array([p["result"] for p in res[name]])
        grids[name[:2]] = vals.reshape(len(q2_nodes), len(x_nodes)).T
    return grids


MC_GEV = 1.51          # matches make_cards()'s "mc"; m_c^2 = 2.2801 GeV^2
DAMP_POWER = 2         # conventional; the docs call the exponent arbitrary
CACHE = f"{BASE}/results/fonll_cache"


def _fonll_components(pto, pdf, flavour, x_nodes, q2_nodes):
    """The three FONLL ingredients, cached to disk.

    The yadism runs are the expensive part and the damped and undamped
    variants share them exactly, so they are computed once per (flavour, pto)
    and reused.  Delete results/fonll_cache to force a recomputation.
    """
    os.makedirs(CACHE, exist_ok=True)
    # >>> THE PDF AND THE TMC MODE ARE PART OF THE KEY TOO (2026-09-13). <<<
    # These are PDF-CONVOLVED structure functions, and the key named neither
    # the set nor the target-mass correction: a tungsten (W184free) or a
    # TMC-on run would have loaded the proton TMC-off arrays, printed
    # "[cache] reusing ..." and written a proton charm reference under a
    # tungsten name.  The pre-existing proton / TMC-off files keep their
    # names (no tag), so every earlier re-run still hits its own cache.
    # >>> THE GRID IS PART OF THE KEY. <<<
    # These arrays are structure functions evaluated ON x_nodes/q2_nodes, and
    # that grid is the FIDUCIAL grid of one beam energy.  The key was
    # (flavour, pto) only, so the 400 GeV and 4 TeV runs of 2026-08-25 loaded
    # the 1 TeV cache written on 2026-08-20 and combined 1 TeV structure
    # functions with their own kinematics.  It printed "[cache] reusing ..."
    # and produced FONLL/ZM = 1.76 at 400 GeV and 0.47 at 4 TeV -- backwards,
    # since mass effects must VANISH at high Q2, not grow.
    #
    # Hashing the grid itself rather than the beam energy is deliberate: the
    # grid is what the arrays actually depend on, so any future change to the
    # binning invalidates the cache too, without anyone having to remember.
    grid_id = hashlib.sha1(
        np.concatenate([np.asarray(x_nodes, dtype=float),
                        np.asarray(q2_nodes, dtype=float)]).tobytes()
    ).hexdigest()[:12]
    fn = f"{CACHE}/{flavour}_pto{pto}_{grid_id}{cache_tag(pdf)}.npz"
    if os.path.exists(fn):
        z = np.load(fn)
        print(f"      [cache] reusing {os.path.relpath(fn, BASE)}", flush=True)
        return ({k: z[f"zm_{k}"] for k in ("F2", "FL", "F3")},
                {k: z[f"ffns_{k}"] for k in ("F2", "FL", "F3")},
                {k: z[f"ffn0_{k}"] for k in ("F2", "FL", "F3")})
    zm, ffns, ffn0 = grid_components("mu", pto, pdf, flavour, x_nodes,
                                     q2_nodes)
    np.savez_compressed(
        fn, **{f"zm_{k}": v for k, v in zm.items()},
        **{f"ffns_{k}": v for k, v in ffns.items()},
        **{f"ffn0_{k}": v for k, v in ffn0.items()})
    return zm, ffns, ffn0


def cache_tag(pdf):
    """"" for the proton set with TMC off (the pre-cache names), else
    "_<set>_tmc<N>" -- see the note in _fonll_components."""
    import target
    name = pdf.set().name
    if name == target.PDFSET_P and target.tmc() == 0:
        return ""
    return f"_{name}_tmc{target.tmc()}"


def grid_components(current, pto, pdf, flavour, x_nodes, q2_nodes):
    """(zm, ffns, ffn0) single-flavour [x, Q2] arrays from PineAPPL grids.

    ON A CACHE MISS THE MASSIVE RUNS GO THROUGH THE GRIDS (2026-09-13), not a
    second direct yadism run.  The grid is the same card on the same nodes --
    verified to 1e-14 against yadism by `pineappl_grids.py check` -- and it is
    PDF-INDEPENDENT, so the scale-band script and this one share ONE build per
    (current, energy, order, scheme, TMC) instead of each paying for the
    massive coefficient functions, which dominate the cost (~12 min a scheme
    per energy, 1.5x that with TMC on).  A pre-existing .npz still wins, so
    every published re-run is unchanged.
    """
    import pineappl_grids as pg
    names = [f"{k}_{flavour}" for k in ("F2", "FL", "F3")]
    nq, nx = len(q2_nodes), len(x_nodes)
    out = []
    for fns, nfff in (("ZM-VFNS", 5), ("FONLL-FFNS", 3), ("FONLL-FFN0", 3)):
        key, _x, _q = pg.build(current, EBEAM, pto, names, x_nodes=x_nodes,
                               q2_nodes=q2_nodes, fns=fns, nfff=nfff)
        print(f"      [grid] {key} x {pdf.set().name}", flush=True)
        out.append({n[:2]: pg.convolve(key, n, pdf).reshape(nq, nx).T
                    for n in names})
    return tuple(out)


def damping_factor(q2_nodes, mass=MC_GEV, power=DAMP_POWER):
    """FONLL threshold damping f_thr(Q2), applied to F^(d).

        f_thr = theta(Q2 - m^2) * (1 - m^2/Q2)^p

    NOTE the yadism docs render this as (1 - Q2/m^2)^2, which GROWS without
    bound above threshold instead of damping; that is a typo, confirmed by
    the user (an author of the FONLL paper).  The form above is the one in
    the paper: it vanishes at Q2 = m^2 and tends to 1 as Q2 -> infinity.

    Why it matters here more than usual: our fiducial region starts at
    Q2 = 4 GeV2 against m_c^2 = 2.28, so f_thr is only
    (1 - 2.28/3.9)^2 = 0.17 at the bottom of the grid and 0.69 at the charm
    <Q2> of 13.25.  Damping is therefore a LARGE effect over most of the
    region, not a threshold detail.
    """
    q2 = np.asarray(q2_nodes, dtype=float)
    m2 = mass * mass
    return np.where(q2 > m2, np.clip(1.0 - m2 / q2, 0.0, None) ** power, 0.0)


def fonll_structure_function_splines(pto, pdf, flavour="charm", damp=False):
    """FONLL general-mass charm structure functions.

    This yadism build has no single "FONLL-A/B/C" card: valid FNS values are
    ZM-VFNS, FONLL-FFNS, FONLL-FFN0, FFNS and FFN0, and FONLL is assembled by
    COMBINING runs.  Per the yadism docs (theory/fns.html),

        F_FONLL = F^(nl+1) + F^(d),    F^(d) = F^(nl) - F^(nl,0)

    with F^(nf) the FFNS result carrying the mass corrections, F^(nf+1) the
    ZM-VFNS result that resums the collinear logs, and F^(nf,0) the FFN0
    asymptotic (massless) limit of the FFNS calculation, subtracted so the
    overlap is not double counted.  Implemented here as

        F_FONLL = F_ZM-VFNS + F_FONLL-FFNS - F_FONLL-FFN0

    The massive pieces run with NfFF = 3 (u, d, s light and charm massive);
    the ZM piece keeps NfFF = 5, which is harmless because ZM-VFNS runs
    nf = nf(Q2) from the thresholds rather than from NfFF.

    TWO DELIBERATE SIMPLIFICATIONS, both documented rather than hidden:
    * **No damping.**  The card sets DAMP = 0, so F^(d) is not suppressed
      near threshold.  The docs give the damped variant
      f_thr = theta(Q2 - m2) (1 - m2/Q2)^p with p a user choice; adding it
      would change the low-Q2 end, which is exactly where our fiducial region
      sits, so this is the main caveat on the numbers.
    * **No switch above the bottom threshold.**  The docs apply the matched
      formula for Q2 < Q2_thr,nf+2 (= mb^2 = 24.2 GeV2 here) and pure ZM-VFNS
      above it; this code applies the matched formula at all Q2.  The effect
      is confined to the tail, where FONLL converges to ZM anyway.

    Sanity behaviour, checked on a coarse grid: FONLL sits ~16% below ZM at
    Q2 = 5 GeV2, where the mass effects bite, and converges to within ~4% of
    it by Q2 = 100 -- the expected shape.

    EXPLORATORY.  The benchmark's stated convention is massless charm, and
    this is provided for completeness alongside the ZM numbers, not as a
    replacement reference.
    """
    # one definition, in yadism_calc.sf_grid_nodes()
    x_nodes, q2_nodes = sf_grid_nodes()
    zm, ffns, ffn0 = _fonll_components(pto, pdf, flavour, x_nodes, q2_nodes)

    # f_thr depends on Q2 only, so it multiplies whole columns of the [x, Q2]
    # grid.  Undamped is the same expression with f_thr == 1.
    f_thr = damping_factor(q2_nodes) if damp else np.ones(len(q2_nodes))
    if damp:
        print(f"      damping f_thr: {f_thr[0]:.4f} at Q2={q2_nodes[0]:.2f}"
              f" ... {f_thr[-1]:.4f} at Q2={q2_nodes[-1]:.0f}", flush=True)

    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    splines = {}
    for k in ("F2", "FL", "F3"):
        # F_FONLL = F^(nl+1) + f_thr * F^(d),  F^(d) = F^(nl) - F^(nl,0)
        #         =   ZM     + f_thr * (FFNS  -  FFN0)
        #
        # WHICH TERM GETS DAMPED MATTERS.  Undamped the grouping is
        # irrelevant -- the sum FFNS + ZM - FFN0 is the same either way -- but
        # f_thr multiplies F^(d) only, so the two groupings diverge as soon as
        # damping is on.  F^(d) is the MASS CORRECTION, massive minus its own
        # massless limit (FONLL paper, arXiv:1001.2312); it is the piece that
        # is unreliable near threshold and therefore the piece to suppress.
        # Damping then interpolates FONLL -> ZM at threshold and -> full FONLL
        # at large Q2, POINTWISE IN Q2.
        #
        # It does NOT follow that the integrated sigma lies between the two.
        # That holds only while F^(d) keeps a fixed sign: at NLO it is
        # negative at all 42 Q2 nodes and sigma_damped duly sits between ZM
        # and undamped FONLL, but at NNLO F^(d) CHANGES SIGN -- it is +206%
        # of ZM at Q2 = 3.9 and negative above Q2 ~ 8 -- so damping suppresses
        # the positive spike ~6x while leaving the negative tail, and
        # sigma_damped (1.8415) falls BELOW both ZM (1.9155) and undamped
        # FONLL (2.0206).  That is arithmetic, not a bug.
        #
        # The yadism docs, as rendered, assign F^(d) = F^(nf+1) - F^(nf,0)
        # (i.e. ZM - FFN0) and add F^(nf) to it.  Damping THAT sends the
        # result to the pure massive FFNS at threshold, and the NLO number
        # came out at 3.544 nb -- OUTSIDE the range spanned by ZM (2.947) and
        # undamped FONLL (2.505), which an interpolation between the two
        # schemes cannot be.  The docs' role descriptions are also internally
        # inconsistent there (F^(nf+1) is described as the FFN0 asymptotic
        # limit), and their damping formula is separately a confirmed typo, so
        # the paper's grouping is used here.
        grid = zm[k] + f_thr[None, :] * (ffns[k] - ffn0[k])
        splines[k] = RectBivariateSpline(lx, lq, grid, kx=3, ky=3)
    return splines


def flavour_structure_function_splines(pto, pdf, flavour="charm"):
    """Run yadism for ONE quark flavour only; return splines in (ln x, ln Q2).

    In ZM-VFNS the `charm` observable is assembled by yadism's Combiner from
    `generate_single_flavor_light(esf, nf, ihq=4)` -- the massless (light)
    coefficient functions carrying the charm couplings and, from O(alpha_s) on,
    the charm share of the gluon/singlet terms.  At PTO = 0 only the
    non-singlet kernel survives, with partons {4: w, -4: w}: exactly
    F2^c = A_c x (c + cbar).
    """
    import yadism

    # one definition, in yadism_calc.sf_grid_nodes()
    x_nodes, q2_nodes = sf_grid_nodes()
    theory, observables = make_cards(pto, x_nodes, q2_nodes)
    kins = observables["observables"]["F2_total"]
    names = [f"{k}_{flavour}" for k in ("F2", "FL", "F3")]
    observables["observables"] = {n: kins for n in names}
    out = yadism.run_yadism(theory, observables)
    res = out.apply_pdf_alphas_alphaqed_xir_xif(
        pdf, lambda muR: pdf.alphasQ(muR), lambda _muR: ALPHA, 1.0, 1.0)

    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    splines = {}
    for name in names:
        vals = np.array([p["result"] for p in res[name]])
        grid = vals.reshape(len(q2_nodes), len(x_nodes)).T   # [x, Q2]
        splines[name[:2]] = RectBivariateSpline(lx, lq, grid, kx=3, ky=3)
    return splines


# ------------------------------------------------------- (b) direct LHAPDF
def direct_flavour_dsig(pdf, flavour="charm"):
    """LO parton-model dsigma from the charm PDF alone, straight out of LHAPDF.

    Same coupling algebra as yadism_calc.manual_lo_sigma() but with the flavour
    sum reduced to pid = +-4, and returned as the two differential functions
    the shared integrators expect, so the closure is bin-by-bin.
    """
    pid = FLAVOURS[flavour]
    el, vl, al = -1.0, -0.5 + 2.0 * S2TW, -0.5
    eq = 2.0 / 3.0 if pid % 2 == 0 else -1.0 / 3.0
    t3 = 0.5 if pid % 2 == 0 else -0.5
    vq, aq = t3 - 2.0 * eq * S2TW, t3
    c1 = 2.0 * el * eq * vl * vq
    c2 = (vl**2 + al**2) * (vq**2 + aq**2)
    d1 = 2.0 * el * al * eq * aq
    d2 = 2.0 * vl * al * 2.0 * vq * aq

    def xf(pid, x, q2):
        x = np.asarray(x, float)
        q2 = np.asarray(q2, float)
        flat = np.array([pdf.xfxQ2(pid, float(xx), float(qq)) if xx < 1.0 else 0.0
                         for xx, qq in zip(x.ravel(), q2.ravel())])
        return flat.reshape(x.shape)

    def dsig_dxdq2(x, q2, y):
        x, q2, y = np.broadcast_arrays(np.asarray(x, float),
                                       np.asarray(q2, float),
                                       np.asarray(y, float))
        kap = q2 / (q2 + MZ * MZ) / (4.0 * S2TW * (1.0 - S2TW))
        a2 = eq * eq + c1 * kap + c2 * kap * kap
        a3 = d1 * kap + d2 * kap * kap
        c = xf(pid, x, q2)
        cb = xf(-pid, x, q2)
        f2 = a2 * (c + cb)
        xf3 = a3 * (c - cb)
        yp = 1.0 + (1.0 - y) ** 2
        ym = 1.0 - (1.0 - y) ** 2
        out = (2.0 * math.pi * ALPHA**2 * yp / (x * q2 * q2)
               * (f2 + ym / yp * xf3) * GEV2_PB)
        return np.where(x > XMAX, 0.0, out)

    def dsig_dq2dy(q2, y):
        q2, y = np.broadcast_arrays(np.asarray(q2, float), np.asarray(y, float))
        x = q2 / (y * TWO_KP)
        return dsig_dxdq2(x, q2, y) * x / y

    return dsig_dxdq2, dsig_dq2dy


# ------------------------------------------------------------------- driver
def run(dsig_dxdq2, dsig_dq2dy):
    hists = {}
    sigma1, means = integrate_q2y(dsig_dq2dy, hists)
    sigma2 = integrate_xq2(dsig_dxdq2, hists)
    sigma3 = integrate_ytheta(dsig_dq2dy, hists)
    return sigma1, sigma2, sigma3, means, hists


def main():
    flavour = (sys.argv[1] if len(sys.argv) > 1 else "charm").lower()
    if flavour not in FLAVOURS:
        sys.exit(f"unknown flavour '{flavour}' (use {'|'.join(FLAVOURS)})")
    order = (sys.argv[2] if len(sys.argv) > 2 else "lo").lower()
    scheme = (sys.argv[3] if len(sys.argv) > 3 else "zm").lower()
    if scheme not in ("zm", "fonll", "fonll_damp"):
        sys.exit("scheme must be zm, fonll or fonll_damp")
    if order not in ORDERS:
        sys.exit(f"unknown order '{order}' (use {'|'.join(ORDERS)})")
    pto = ORDERS[order]
    tag = f"[{flavour} {order} {scheme}]"

    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDFSET, 0)

    print(f"{tag} (a) yadism {order.upper()} (PTO={pto}), "
          f"F2/FL/F3_{flavour}, ZM-VFNS ...", flush=True)
    if pto > 0:
        print(f"{tag}     NOTE: at PTO > 0 this includes the charm share of "
              f"the gluon/singlet terms (gamma* g -> c cbar), so it is a "
              f"WIDER process than the LO charm-initiated reference.")
    if scheme.startswith("fonll"):
        splines = fonll_structure_function_splines(
            pto, pdf, flavour, damp=(scheme == "fonll_damp"))
    else:
        splines = flavour_structure_function_splines(pto, pdf, flavour)
    dxq_a, dqy_a = make_sigma_red(splines)
    sa1, sa2, sa3, means_a, hists_a = run(dxq_a, dqy_a)
    print(f"{tag} (a) yadism   sigma_fid = {sa1/1e3:.5f} nb "
          f"(cross-checks: x-Q2 {sa2/1e3:.5f}, y-theta {sa3/1e3:.5f})")

    sb1 = None
    if pto == 0 and scheme == "zm":
        print(f"{tag} (b) direct LHAPDF {flavour}-PDF integral ...", flush=True)
        dxq_b, dqy_b = direct_flavour_dsig(pdf, flavour)
        sb1, sb2, sb3, means_b, hists_b = run(dxq_b, dqy_b)
        print(f"{tag} (b) direct   sigma_fid = {sb1/1e3:.5f} nb "
              f"(cross-checks: x-Q2 {sb2/1e3:.5f}, y-theta {sb3/1e3:.5f})")
        print(f"{tag} ratio (a)/(b) = {sa1/sb1:.6f}   "
              f"[{1e3*(sa1/sb1 - 1):+.3f} permille]")

        # bin-by-bin closure of the two methods
        worst = []
        for key in ("Q2", "xbj", "y", "Emu", "theta", "nu"):
            a, b = hists_a[key], hists_b[key]
            m = b > 0
            rel = np.abs(a[m] / b[m] - 1.0)
            worst.append((key, float(rel.max()), int(np.argmax(rel))))
        print(f"{tag} worst bin-by-bin (a)/(b) deviation per observable:")
        for key, r, i in worst:
            print(f"          {key:6s} {1e3*r:7.3f} permille (bin {i})")
    else:
        print(f"{tag} (b) SKIPPED: the direct LHAPDF integral is a LO "
              f"parton-model calculation and is not a valid check above LO.")

    means = means_a
    print(f"{tag} means: <Q2> = {means['Q2']:.3f} GeV2, "
          f"<x> = {means['xbj']:.5f}, <y> = {means['y']:.4f}")

    terms = [abs(sa1 - sa2), abs(sa1 - sa3), 2e-4 * sa1]
    if sb1 is not None:
        terms.append(abs(sa1 - sb1))
    err = max(terms)
    suffix = "" if order == "lo" else f"_{order}"
    if scheme == "fonll":
        suffix += "_fonll"
    elif scheme == "fonll_damp":
        suffix += "_fonll_damp"
    out = {
        "generator": f"yadism_{flavour}{suffix}",
        "label": (f"YADISM {order.upper()} ({flavour}"
                  + {"fonll": ", FONLL", "fonll_damp": ", FONLL"}.get(
                      scheme, "") + ")"),
        "scheme": {"fonll": "FONLL general-mass (undamped)",
                   "fonll_damp": (f"FONLL general-mass, damped "
                                  f"(1-mc^2/Q2)^{DAMP_POWER}, mc={MC_GEV}"),
                   }.get(scheme, "ZM-VFNS"),
        "n_parsed": None,
        "n_fiducial": None,
        "sigma_fid_pb": sa1,
        "sigma_fid_err_pb": err,
        "sigma_fid_pb_direct_lhapdf": sb1,
        "pto": pto,
        "includes_gluon_channel": pto > 0,
        "means": {k: float(v) for k, v in means.items()},
        "hists": {},
    }
    for key in ("Q2", "xbj", "y", "Emu", "theta", "nu"):
        edges = np.asarray(BINS[key], dtype=float)
        widths = np.diff(edges)
        out["hists"][key] = {
            "edges": edges.tolist(),
            "dsig": (hists_a[key] / widths).tolist(),
            "err": np.zeros(len(widths)).tolist(),
        }

    os.makedirs(f"{BASE}/results", exist_ok=True)
    # per-energy result key; untagged at the 1 TeV anchor
    # THE SELECTION, THE TARGET AND THE TMC MODE ARE IN THE NAME (2026-09-13),
    # in the order yadism_calc uses, histos_<gen>[_<sel>][_<target>][_tmc]
    # [_<energy>].  They were not: a charm reference computed in any other
    # region or on any other target was written straight over the inclusive
    # proton one, which yadism_fonll_inclusive.py then combined with a
    # correctly named inclusive result.  All three are "" at the defaults,
    # so nothing published moves.
    import target
    ofn = (f"{BASE}/results/"
           f"{beams_mod.at_energy(f'histos_yadism_{flavour}{suffix}' + sel_suffix() + target.suffix() + target.tmc_suffix(), EBEAM)}"
           f"{selection.q2_suffix()}.json")
    with open(ofn, "w") as f:
        stamp(out, current="mu")
        prov = target.provenance(PDFSET)
        if prov:              # absent at proton / TMC off: earlier stays byte-identical
            out.update(prov)
        json.dump(out, f)
    print(f"{tag} wrote {ofn}")


if __name__ == "__main__":
    main()
