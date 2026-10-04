#!/usr/bin/env python3
"""YADISM inclusive (analytic) CC neutrino-DIS reference for the generator benchmark.

Charged-current counterpart of analysis/yadism_calc.py: 1 TeV nu_mu on a proton
at rest, W+ exchange, outgoing lepton mu-.  Computes d^2sigma/dQ2dy from the CC
structure functions F2/FL/xF3 convolved with NNPDF40_nnlo_as_01180, integrates
over the fiducial region Q2 > 4 GeV2, 0.2 < y < 0.9, and writes
results_nu/histos_yadism[_<order>].json in the same schema and with the SAME
bin edges as analysis/analyze.py (BINS imported, not copied), so the neutrino
and muon studies are directly comparable.

Conventions (matched to the Sherpa/Pythia/Herwig CC setup):
  * Gmu scheme: G_F = 1.1663787e-5 GeV^-2, M_W = 80.379 GeV, M_Z = 91.1876,
    sin2thetaW = 0.22291.  For CC the whole EW normalisation is G_F^2 M_W^4;
    the equivalent form pi alpha^2 / (8 x sin^4tw (Q2+M_W^2)^2) with
    alpha = 1/132.119 differs by ~0.2% (see NEUTRINO_NOTES.md) -- G_F is the
    quantity the Gmu scheme fixes exactly, so it is the one used here.
  * W propagator: real, 1/(Q2 + M_W^2)  (generators use the spacelike
    t-channel propagator without the width; including M_W*Gamma_W would shift
    sigma by -0.07%).
  * muF = muR = Q (XIF = XIR = 1)
  * ZM-VFNS: massless charm and bottom coefficient functions, NNPDF4.0 NNLO PDF
  * CKM from the standard theory card, top row masked out by yadism's nf = 5
    ZM treatment -> the b channel is suppressed to |Vub|^2 + |Vcb|^2 = 0.0017
    and no top is produced, exactly as in the generators.
  * alpha_s(muR) from the LHAPDF set itself at every order
  * TMC from $BENCH_TMC (target.tmc(); default 0 = none)

Cross-section formula (neutrino, i.e. +xF3):

  d2sigma/dx dQ2 = G_F^2 / (4 pi x) * [M_W^2/(M_W^2+Q2)]^2
                   * (Y+ F2 - y^2 FL + Y- xF3)
  Y+- = 1 +- (1-y)^2

verified to 1e-5 against yadism's own XSFPFCC_total observable (which is
literally this expression in pb/GeV^2), and its LO F2/xF3 verified against a
direct LHAPDF parton-model evaluation
  F2  = 2x [ sum_d CKM_d (d)      + sum_ubar CKM_u (ubar) ]
  xF3 = 2x [ sum_d CKM_d (d)      - sum_ubar CKM_u (ubar) ]
to 1e-4.

x = Q2 / (y * 2 k.P), 2 k.P = 2 E_nu M_P = 1876.544 GeV^2 (target at rest,
massless neutrino -> exact).

Lab-frame observables (target at rest, exact):
  E_mu' = (1 - y) * E_beam            (from P.k' = (1-y) P.k)
  nu    = y * E_beam
  theta: Q2 = 2 (E E' - |k||k'| cos th) - m_mu^2,  |k| = E (massless nu),
         |k'| = sqrt(E'^2 - m_mu^2)

Usage:  yadism_cc_calc.py [lo|nlo|nnlo]     (default lo)
"""
import json
import math
import os
import sys

import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import (BINS, M_MU, M_P, Q2_MIN, Y_MAX,  # noqa: E402
                     Y_MIN, sel_suffix)
# ONE SOURCE FOR THE W^2 BOUND, imported rather than re-derived: the muon and
# neutrino integrators used to carry their own copies of the kinematic limits
# and that is exactly the duplication CONVENTIONS.md rule 2 asks to avoid.  W2_MIN
# is None for every published selection on either current, so these are
# no-ops here today; they exist so that a charged-current comparison against a
# (Q, W) region costs a selection and not a code change (rule 2b).
from yadism_calc import (q2_max_from_w2, q2_min_from_w2,  # noqa: E402
                         y_min_from_w2, y_lattice)
import beams as beams_mod  # noqa: E402
import selection  # noqa: E402  -- the one Q2 floor
import target  # noqa: E402  -- the one target nucleon
from beams import Beams  # noqa: E402
# the (x, Q2) interpolation grid has ONE definition, and it scales with
# the beam energy -- see yadism_calc.sf_grid_nodes()
from yadism_calc import install_nan_audit, sf_grid_nodes  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = f"{BASE}/results_nu"

# ------------------------------------------------------------------ constants
# GeV, lab (nu_mu).  From the environment so the CC reference follows the
# scan exactly as the NC one does -- CONVENTIONS.md rule 2b, the two currents stay
# in step.  TWO_KP and everything derived from it follow.
EBEAM = float(os.environ.get("BENCH_ENERGY") or beams_mod.ANCHOR_ENERGY)
ALPHA = 1.0 / 132.119                # only used for the theory card's alphaqed
S2TW = 0.22291
MZ = 91.1876
MW = 80.379                          # Sherpa/Pythia default W mass
MW2 = MW * MW
GF = 1.1663787e-05                   # GeV^-2
# 2 k.P comes from analysis/beams.py, the one place that turns
# "E GeV lepton on a proton at rest" into each code's beam inputs.
# Numerically identical to the literal it replaces at 1 TeV; the point
# is that an energy scan changes ONE line rather than fifteen cards.
TWO_KP = Beams("nu", EBEAM).two_kP   # 2 k.P = 1876.544 GeV^2 (exact)
GEV2_PB = 3.893793e8                 # (1 GeV^-2) in pb
XMAX = 0.999
# THE KINEMATIC y CEILING, and the neutrino half of a fix the muon side
# already had (CONVENTIONS.md rule 2b).  The scattered lepton is a MUON on this
# current too, so it cannot carry less than its own mass: y <= 1 - m_mu/E.
# integrate_ytheta evaluates sqrt(E'^2 - m_mu^2), which dies with "math domain
# error" once the region reaches the endpoint -- exactly what happened the
# first time the all-y region was run below 100 GeV, and exactly what
# yadism_calc.py had already been taught for dis2506 while this file was left
# with Y_MAX.  At 1 TeV the ceiling is 0.99989 and at 10 GeV 0.98943, so it is
# a no-op wherever Y_MAX is 0.9 and every existing result is untouched.
Y_HI = min(Y_MAX, 1.0 - M_MU / EBEAM)
Q2_KIN_MAX = Y_HI * TWO_KP * XMAX
# THE TARGET NUCLEON, and with it the PDF set: BENCH_TARGET=n swaps to the
# isospin-mirrored NNPDF4.0 and stamps "_n" on the filename.  Default is the
# proton, so every existing result is bit-for-bit what it was.  See
# analysis/target.py for why tungsten needs the neutron at all.
PDFSET = target.pdfset()
CKM_STR = ("0.97428 0.22530 0.003470 0.22520 0.97345 0.041000 "
           "0.00862 0.04030 0.999152")

ORDERS = {"lo": 0, "nlo": 1, "nnlo": 2}


# ------------------------------------------------------------------- yadism
def make_cards(pto, x_nodes, q2_nodes):
    install_nan_audit()        # see yadism_calc.install_nan_audit
    theory = {
        "ID": 0, "PTO": pto, "PTODIS": None, "FNS": "ZM-VFNS", "NfFF": 5,
        "DAMP": 0, "IC": 0, "IB": 0,
        "ModEv": "EXA", "ModSV": None, "XIR": 1.0, "XIF": 1.0,
        "fact_to_ren_scale_ratio": 1.0, "RenScaleVar": True, "FactScaleVar": True,
        "MaxNfAs": 5, "MaxNfPdf": 5,
        "Q0": 1.65, "nf0": 4, "alphas": 0.1180, "Qref": MZ, "nfref": 5,
        "QED": 0, "alphaqed": ALPHA, "Qedref": 1.777,
        "SIN2TW": S2TW, "MZ": MZ, "MW": MW, "GF": GF,
        "CKM": CKM_STR,
        "mc": 1.51, "Qmc": 1.51, "kcThr": 1.0,
        "mb": 4.92, "Qmb": 4.92, "kbThr": 1.0,
        "mt": 172.5, "Qmt": 172.5, "ktThr": 1.0,
        # TMC from the one knob -- see yadism_calc.make_cards
        "MP": M_P, "TMC": target.tmc(), "n3lo_cf_variation": 0,
        "EScaleVar": 1,
    }
    kins = [{"x": float(x), "Q2": float(q2), "y": 0.5}
            for q2 in q2_nodes for x in x_nodes]
    xgrid = np.concatenate([np.geomspace(1e-4, 0.1, 40, endpoint=False),
                            np.linspace(0.1, 1.0, 41)]).tolist()
    observables = {
        "observables": {"F2_total": kins, "FL_total": kins, "F3_total": kins},
        "prDIS": "CC",
        "ProjectileDIS": "neutrino",   # nu_mu: identical CC couplings to nu_e
        "TargetDIS": "proton",
        "PolarizationDIS": 0.0,
        "PropagatorCorrection": 0.0,
        "NCPositivityCharge": None,
        "interpolation_xgrid": xgrid,
        "interpolation_polynomial_degree": 4,
        "interpolation_is_log": True,
    }
    return theory, observables


def structure_function_splines(pto, pdf):
    """Run yadism on an (x, Q2) grid, return cubic splines in (ln x, ln Q2)."""
    import yadism

    # one definition, in yadism_calc.sf_grid_nodes()
    x_nodes, q2_nodes = sf_grid_nodes()
    theory, observables = make_cards(pto, x_nodes, q2_nodes)
    out = yadism.run_yadism(theory, observables)
    res = out.apply_pdf_alphas_alphaqed_xir_xif(
        pdf, lambda muR: pdf.alphasQ(muR), lambda _muR: ALPHA, 1.0, 1.0)

    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    splines = {}
    for name in ("F2_total", "FL_total", "F3_total"):
        vals = np.array([p["result"] for p in res[name]])
        grid = vals.reshape(len(q2_nodes), len(x_nodes)).T   # [x, Q2]
        splines[name[:2]] = RectBivariateSpline(lx, lq, grid, kx=3, ky=3)
    return splines


def make_sigma_red(splines):
    """d2sigma/dxdQ2 and d2sigma/dQ2dy in pb (CC, neutrino: +xF3)."""
    X_LO = sf_grid_nodes()[0][0]        # the interpolation grid's lower edge

    def sf(name, x, q2):
        x, q2 = np.broadcast_arrays(np.asarray(x, float), np.asarray(q2, float))
        # CLIP TO THE GRID'S OWN LOWER EDGE, not a literal -- this was
        # np.clip(x, 0.002, XMAX), the 1 TeV grid's lower node.  See the same
        # fix in yadism_calc.make_sigma_red: at 4 TeV it silently evaluated
        # every x below 0.002 AT 0.002 and cost 2.1% on sigma.
        v = splines[name].ev(np.log(np.clip(x, X_LO, XMAX)), np.log(q2))
        return np.where(x > XMAX, 0.0, v)

    def dsig_dxdq2(x, q2, y):
        """d2sigma/dxdQ2 in pb/GeV2 (x, q2, y broadcastable arrays/scalars)."""
        yp = 1.0 + (1.0 - y) ** 2
        ym = 1.0 - (1.0 - y) ** 2
        comb = (yp * sf("F2", x, q2) - y * y * sf("FL", x, q2)
                + ym * sf("F3", x, q2))
        prop = (MW2 / (MW2 + q2)) ** 2
        return GF * GF / (4.0 * math.pi * x) * prop * comb * GEV2_PB

    def dsig_dq2dy(q2, y):
        """d2sigma/dQ2dy in pb/GeV2; q2, y broadcastable. Jacobian x/y."""
        x = q2 / (y * TWO_KP)
        return dsig_dxdq2(x, q2, y) * x / y
    return dsig_dxdq2, dsig_dq2dy


# ------------------------------------------------------------ quadrature grids
def gl_panels(edges, nsub, order, log=False):
    """Gauss-Legendre nodes/weights on each [edges[i], edges[i+1]] split into
    nsub sub-panels (in log space if log=True). Returns (nodes, weights, ibin)."""
    t, w = np.polynomial.legendre.leggauss(order)
    nodes, wgts, ibin = [], [], []
    for i in range(len(edges) - 1):
        a, b = edges[i], edges[i + 1]
        subs = np.geomspace(a, b, nsub + 1) if log else np.linspace(a, b, nsub + 1)
        for j in range(nsub):
            lo, hi = subs[j], subs[j + 1]
            if log:
                la, lb = math.log(lo), math.log(hi)
                xx = np.exp(0.5 * (la + lb) + 0.5 * (lb - la) * t)
                ww = 0.5 * (lb - la) * w * xx
            else:
                xx = 0.5 * (lo + hi) + 0.5 * (hi - lo) * t
                ww = 0.5 * (hi - lo) * w
            nodes.append(xx); wgts.append(ww); ibin.append(np.full(order, i))
    return (np.concatenate(nodes), np.concatenate(wgts),
            np.concatenate(ibin, dtype=int))


def gl_interval(a, b, order, log=False):
    t, w = np.polynomial.legendre.leggauss(order)
    if log:
        la, lb = math.log(a), math.log(b)
        xx = np.exp(0.5 * (la + lb) + 0.5 * (lb - la) * t)
        return xx, 0.5 * (lb - la) * w * xx
    return 0.5 * (a + b) + 0.5 * (b - a) * t, 0.5 * (b - a) * w


def theta_of(q2, y):
    """Lab muon scattering angle; massless incoming nu, massive outgoing mu."""
    ep = EBEAM * (1.0 - y)
    k = EBEAM                                    # massless neutrino
    kp = np.sqrt(np.maximum(ep * ep - M_MU**2, 0.0))
    c = (2.0 * EBEAM * ep - M_MU**2 - q2) / (2.0 * k * kp)
    return np.arccos(np.clip(c, -1.0, 1.0))


# ----------------------------------------------------------------- integrators
def integrate_q2y(dsig_dq2dy, hist_accum):
    """Master integration on (Q2, y): sigma_fid, means, Q2/y/Emu/nu hists."""
    q2n, q2w, q2i = gl_panels(BINS["Q2"], 3, 8, log=True)
    q2o, q2ow = gl_interval(BINS["Q2"][-1], Q2_KIN_MAX, 8, log=True)
    q2_all = np.concatenate([q2n, q2o])
    q2w_all = np.concatenate([q2w, q2ow])
    q2i_all = np.concatenate([q2i, np.full(len(q2o), -1)])  # -1: Q2 overflow

    # ends at Y_HI -- see yadism_calc.y_lattice for the strip it used to drop
    y_lat = y_lattice(Y_MIN, Y_HI)
    tg, wg = np.polynomial.legendre.leggauss(4)

    sigma = 0.0
    sums = {"Q2": 0.0, "xbj": 0.0, "y": 0.0}
    for key in ("Q2", "y", "Emu", "nu"):
        hist_accum.setdefault(key, np.zeros(len(BINS[key]) - 1))

    emu_edges, nu_edges, y_edges = BINS["Emu"], BINS["nu"], BINS["y"]
    for q2, wq, iq in zip(q2_all, q2w_all, q2i_all):
        ymin_eff = max(Y_MIN, q2 / (XMAX * TWO_KP), y_min_from_w2(q2, TWO_KP))
        if ymin_eff >= Y_HI:
            continue
        lat = y_lat[y_lat > ymin_eff]
        edges = np.concatenate([[ymin_eff], lat]) if lat.size else \
            np.array([ymin_eff, Y_HI])
        yy = (0.5 * (edges[:-1] + edges[1:])[:, None]
              + 0.5 * np.diff(edges)[:, None] * tg).ravel()
        wy = (0.5 * np.diff(edges)[:, None] * wg).ravel()

        f = dsig_dq2dy(q2, yy)
        w = wq * wy * f
        sigma += w.sum()
        sums["Q2"] += q2 * w.sum()
        sums["y"] += (yy * w).sum()
        sums["xbj"] += (q2 / (yy * TWO_KP) * w).sum()

        if iq >= 0:
            hist_accum["Q2"][iq] += w.sum()
        ib = np.searchsorted(y_edges, yy) - 1
        np.add.at(hist_accum["y"], np.clip(ib, 0, len(y_edges) - 2), w)
        emu = EBEAM * (1.0 - yy)
        ib = np.searchsorted(emu_edges, emu) - 1
        np.add.at(hist_accum["Emu"], np.clip(ib, 0, len(emu_edges) - 2), w)
        nu = EBEAM * yy
        ib = np.searchsorted(nu_edges, nu) - 1
        np.add.at(hist_accum["nu"], np.clip(ib, 0, len(nu_edges) - 2), w)
    means = {k: s / sigma for k, s in sums.items()}
    return sigma, means


def integrate_xq2(dsig_dxdq2, hist_accum):
    """xbj histogram on (x, Q2) with panels aligned to the xbj bin edges."""
    xn, xw, xi = gl_panels(BINS["xbj"], 3, 8, log=True)
    hist_accum.setdefault("xbj", np.zeros(len(BINS["xbj"]) - 1))
    sigma = 0.0
    for x, wx, ix in zip(xn, xw, xi):
        q2lo = max(Q2_MIN, Y_MIN * x * TWO_KP, q2_min_from_w2(x))
        q2hi = Y_HI * x * TWO_KP
        if q2hi <= q2lo or x > XMAX:
            continue
        qq, wq = gl_panels(np.geomspace(q2lo, q2hi, 21), 1, 4, log=True)[:2]
        yy = qq / (x * TWO_KP)
        f = dsig_dxdq2(np.full(len(qq), x), qq, yy)
        s = (wq * f).sum() * wx
        hist_accum["xbj"][ix] += s
        sigma += s
    return sigma


def integrate_ytheta(dsig_dq2dy, hist_accum):
    """theta histogram on (y, theta) with panels aligned to the theta bins."""
    th_edges = BINS["theta"]
    hist_accum.setdefault("theta", np.zeros(len(th_edges) - 1))
    yn, yw = gl_panels(y_lattice(Y_MIN, Y_HI), 1, 4, log=False)[:2]
    k = EBEAM
    sigma = 0.0
    for y, wy in zip(yn, yw):
        ep = EBEAM * (1.0 - y)
        kp = math.sqrt(ep * ep - M_MU**2)
        q2max = min(y * TWO_KP * XMAX, q2_max_from_w2(y, TWO_KP),
                    2.0 * EBEAM * ep - M_MU**2 + 2.0 * k * kp)
        thmin = float(theta_of(Q2_MIN, y))
        thmax = float(theta_of(q2max, y))
        for i in range(len(th_edges) - 1):
            a, b = max(th_edges[i], thmin), min(th_edges[i + 1], thmax)
            if b <= a:
                continue
            tt, wt = gl_interval(a, b, 6, log=False)
            q2 = 2.0 * (EBEAM * ep - k * kp * np.cos(tt)) - M_MU**2
            q2 = np.clip(q2, Q2_MIN, q2max)
            jac = 2.0 * k * kp * np.sin(tt)         # dQ2/dtheta
            f = dsig_dq2dy(q2, y)
            s = wy * (wt * f * jac).sum()
            hist_accum["theta"][i] += s
            sigma += s
    return sigma


# -------------------------------------------------------- LO closure check
def ckm2_masked():
    """CKM^2 sums as used by yadism's nf = 5 ZM CC treatment (top row removed).

    Returns (w_down[d,s,b], w_up[u,c]) : the coefficient multiplying each
    down-type quark and each up-type ANTIquark in F2^{nu p}/2x."""
    m = np.array([float(v) for v in CKM_STR.split()]).reshape(3, 3) ** 2
    w_down = m[:2, :].sum(axis=0)     # d, s, b : sum over produced u, c
    w_up = m[:2, :].sum(axis=1)       # u, c: sum over produced d, s, b
    return w_down, w_up


def manual_lo_sigma(pdf, per_channel=False):
    """Direct (spline-free) LO parton-model CC sigma_fid from LHAPDF, in pb.

    If per_channel, returns a dict {pid: sigma} of the incoming-parton
    breakdown (pid > 0: down-type quark, pid < 0: up-type antiquark)."""
    w_down, w_up = ckm2_masked()
    chan = {1: w_down[0], 3: w_down[1], 5: w_down[2],
            -2: w_up[0], -4: w_up[1]}

    q2n, q2w = gl_panels(np.geomspace(Q2_MIN, Q2_KIN_MAX, 40), 1, 6, log=True)[:2]
    tg, wg = np.polynomial.legendre.leggauss(6)
    tot = {p: 0.0 for p in chan}
    for q2, wq in zip(q2n, q2w):
        ymin_eff = max(Y_MIN, q2 / (XMAX * TWO_KP), y_min_from_w2(q2, TWO_KP))
        if ymin_eff >= Y_HI:
            continue
        sub = np.linspace(ymin_eff, Y_HI, 31)
        yy = (0.5 * (sub[:-1] + sub[1:])[:, None]
              + 0.5 * np.diff(sub)[:, None] * tg).ravel()
        wy = (0.5 * np.diff(sub)[:, None] * wg).ravel()
        xx = q2 / (yy * TWO_KP)
        yp = 1.0 + (1.0 - yy) ** 2
        ym = 1.0 - (1.0 - yy) ** 2
        prop = (MW2 / (MW2 + q2)) ** 2
        pref = GF * GF / (4.0 * math.pi * xx) * prop * GEV2_PB * xx / yy
        for pid, cw in chan.items():
            f = np.array([pdf.xfxQ2(pid, x, q2) if x < 1 else 0.0 for x in xx])
            # F2 = 2x f, xF3 = +2x f (quark) or -2x f (antiquark)
            s3 = 1.0 if pid > 0 else -1.0
            comb = cw * 2.0 * f * (yp + s3 * ym)
            tot[pid] += wq * (wy * pref * comb).sum()
    if per_channel:
        return tot
    return sum(tot.values())


# ------------------------------------------------------------------------ main
def main():
    order = (sys.argv[1] if len(sys.argv) > 1 else "lo").lower()
    if order not in ORDERS:
        sys.exit(f"unknown order '{order}' (use lo|nlo|nnlo)")
    pto = ORDERS[order]

    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDFSET, 0)

    print(f"[yadism_cc] running yadism CC (nu_mu) at {order.upper()} "
          f"(PTO={pto}) ...")
    splines = structure_function_splines(pto, pdf)
    dsig_dxdq2, dsig_dq2dy = make_sigma_red(splines)

    hists = {}
    sigma1, means = integrate_q2y(dsig_dq2dy, hists)
    sigma2 = integrate_xq2(dsig_dxdq2, hists)
    sigma3 = integrate_ytheta(dsig_dq2dy, hists)
    print(f"[yadism_cc] sigma_fid = {sigma1:.4f} pb "
          f"(cross-checks: x-Q2 {sigma2:.4f}, y-theta {sigma3:.4f})")

    if order == "lo":
        ch = manual_lo_sigma(pdf, per_channel=True)
        direct = sum(ch.values())
        print(f"[yadism_cc] LO closure: direct LHAPDF CC = {direct:.4f} pb "
              f"(yadism/direct = {sigma1/direct:.5f})")
        names = {1: "d", 3: "s", 5: "b", -2: "ubar", -4: "cbar"}
        for pid in (1, 3, 5, -2, -4):
            print(f"              channel {names[pid]:>4}: {ch[pid]:8.4f} pb "
                  f"({100*ch[pid]/direct:5.2f}%)")

    err = max(abs(sigma1 - sigma2), abs(sigma1 - sigma3), 2e-4 * sigma1)

    gen = "yadism" if order == "lo" else f"yadism_{order}"
    out = {
        "generator": gen,
        "label": f"YADISM CC {order.upper()}",
        "n_parsed": None,
        "n_fiducial": None,
        "sigma_fid_pb": sigma1,
        "sigma_fid_err_pb": err,
        "means": {k: float(v) for k, v in means.items()},
        "hists": {},
    }
    for key in ("Q2", "xbj", "y", "Emu", "theta", "nu"):
        edges = np.asarray(BINS[key], dtype=float)
        widths = np.diff(edges)
        out["hists"][key] = {
            "edges": edges.tolist(),
            "dsig": (hists[key] / widths).tolist(),
            "err": np.zeros(len(widths)).tolist(),
        }

    os.makedirs(OUTDIR, exist_ok=True)
    # per-energy result key; untagged at the 1 TeV anchor
    # The selection belongs in the name for the same reason it does on the
    # muon side (see analysis/yadism_calc.py): without it, a run under any
    # non-inclusive selection overwrites the inclusive reference in silence.
    # sel_suffix() is "" for the inclusive selection, so nothing existing moves.
    # THE TARGET IS IN THE NAME for the same reason the selection is -- a
    # neutron cross-section is 1.8x the proton one on this current and would
    # otherwise overwrite it in silence -- and in the same PLACE, before the
    # energy tag, so that one lookup rule finds every file.
    ofn = (f"{OUTDIR}/"
           f"{beams_mod.at_energy(f'histos_{gen}' + sel_suffix() + target.suffix() + target.tmc_suffix(), EBEAM)}"
           f"{selection.q2_suffix()}.json")
    with open(ofn, "w") as f:
        stamp(out, current="nu")
        prov = target.provenance(PDFSET)
        if prov:              # absent at proton / TMC off: earlier stays byte-identical
            out.update(prov)
        json.dump(out, f)
    print(f"[yadism_cc] wrote {ofn}")


if __name__ == "__main__":
    main()
