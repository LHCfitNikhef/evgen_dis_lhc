#!/usr/bin/env python3
"""YADISM inclusive (analytic) NC muon-DIS reference for the generator benchmark.

Computes d^2sigma/dQ2dy for NC (gamma + gamma-Z + Z) DIS of a 1 TeV mu- on a
proton at rest from yadism structure functions F2/FL/xF3 convolved with
NNPDF40_nnlo_as_01180, integrates it over the fiducial region
Q2 > 4 GeV2, 0.2 < y < 0.9, and writes results/histos_yadism[_<order>].json
in the same schema as analysis/analyze.py (identical binning, imported from
analyze.BINS).

Conventions (matched to the Sherpa/Pythia benchmark setup):
  * fixed alpha_em = 1/132.119 (enters only the overall 2 pi alpha^2 prefactor;
    yadism's NC coupling weights carry no alpha at all: the Z admixture enters
    via kappa = Q2/(Q2+MZ^2) / (4 sin2tw (1-sin2tw)) exactly as in the
    generators' Gmu-scheme couplings)
  * sin2thetaW = 0.22291, MZ = 91.1876 GeV, no Z width (t-channel, spacelike)
  * muF = muR = Q (XIF = XIR = 1)
  * ZM-VFNS: massless charm and bottom coefficient functions, NNPDF4.0 NNLO PDF
  * alpha_s(muR) taken from the LHAPDF set itself at every order
  * TMC from $BENCH_TMC (target.tmc(); default 0 = none), massless
    cross-section formula

Cross-section formula (massless leptons/hadron in the flux and ME, exactly the
parton-model formula the LO generators implement):

  d2sigma/dx dQ2 = 2 pi alpha^2 Y+ / (x Q^4) * sigma_red
  sigma_red      = F2 - y^2/Y+ FL + (Y-/Y+) xF3        (+ sign: mu^- beam)
  Y+- = 1 +- (1-y)^2

with x = Q2 / (y * 2 k.P), 2 k.P = 2 E_mu M_P = 1876.544 GeV^2 (exact for a
target at rest, all masses kept in the invariants).

Lab-frame observables (target at rest, exact to all orders in masses):
  E_mu' = (1 - y) * E_beam            (from P.k' = (1-y) P.k)
  nu    = y * E_beam
  theta: Q2 = 2 (E E' - |k||k'| cos th) - 2 m_mu^2,  |k| = sqrt(E^2 - m_mu^2)

Usage:  yadism_calc.py [lo|nlo|nnlo]     (default lo)
"""
import json
import math
import os
import sys

import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import (BINS, M_MU, M_P, Q2_MIN, W2_MIN, Y_MAX,  # noqa: E402
                     Y_MIN, sel_suffix)
import target  # noqa: E402  -- the one target nucleon
import beams as beams_mod  # noqa: E402
import selection  # noqa: E402  -- the one Q2 floor
from beams import Beams  # noqa: E402
from runmeta import stamp, write_result  # noqa: E402  -- the one (selection, energy) stamp

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ------------------------------------------------------------------ constants
# GeV, lab.  From the environment so the analytic reference follows the scan:
#   BENCH_ENERGY=400 analysis/yadism_calc.py lo
# Everything below that depends on the beam (TWO_KP, Q2_KIN_MAX, the BINS
# imported from analyze) follows from this one value.
EBEAM = float(os.environ.get("BENCH_ENERGY") or beams_mod.ANCHOR_ENERGY)
ALPHA = 1.0 / 132.119                # fixed (Gmu-scheme value, as in Sherpa/Pythia)
S2TW = 0.22291
MZ = 91.1876
# 2 k.P comes from analysis/beams.py, the one place that turns
# "E GeV lepton on a proton at rest" into each code's beam inputs.
# Numerically identical to the literal it replaces at 1 TeV; the point
# is that an energy scan changes ONE line rather than fifteen cards.
TWO_KP = Beams("mu", EBEAM).two_kP   # 2 k.P = 1876.544 GeV^2 (exact, target at rest)

# THE KINEMATIC CEILING ON y.  The scattered muon carries E' = E (1 - y), and
# it cannot carry less than its own rest mass, so y <= 1 - m_mu/E.  Every
# selection the benchmark publishes stops at y = 0.9 and never comes near it,
# which is why this was not needed until a region without a y window appeared:
# dis2506 runs to y = 1, and integrate_ytheta then evaluates
# sqrt(E'^2 - m_mu^2) on a negative argument and dies with "math domain
# error".  At 1 TeV the ceiling is 0.99989 and at 10 GeV it is 0.98943, so
# this is a no-op wherever Y_MAX is 0.9 and only bites where the region
# genuinely reaches the endpoint.
Y_HI = min(Y_MAX, 1.0 - M_MU / EBEAM)
GEV2_PB = 3.893793e8                 # (1 GeV^-2) in pb
XMAX = 0.999                         # SF grid upper edge; F2(x>XMAX) ~ 0
Q2_KIN_MAX = Y_HI * TWO_KP * XMAX    # fiducial kinematic limit (x < XMAX)
# THE TARGET NUCLEON, and with it the PDF set (CONVENTIONS.md rule 2b: the muon
# side is owed whatever the neutrino side gets).  BENCH_TARGET=n swaps to the
# isospin-mirrored NNPDF4.0 and stamps "_n" on the filename; the default is the
# proton, so every existing result is bit-for-bit what it was.
PDFSET = target.pdfset()

ORDERS = {"lo": 0, "nlo": 1, "nnlo": 2}

# Beyond this the LO closure is treated as a failure, not a wobble.
# 400 GeV and 1 TeV both close to 1.00000; 4 TeV does not.
CLOSURE_TOL = 5e-3


# ------------------------------------------------------------------- yadism
def install_nan_audit():
    """Make yadism SAY which structure-function points it zeroes as NaN.

    >>> yadism's "CRITICAL Some NaNs are encountered and set to zero!" IS NOT
    A NaN REPORT (2026-09-13). <<<  Runner.replace_nans_with_0 (yadism
    0.13.11, runner.py:218) logs that line once for every key of the output
    that is NOT an observable -- xgrid, polynomial_degree, is_log, pids and
    projectilePID, hence always exactly five times -- and then zeroes any non-finite
    coefficient in the real observables WITHOUT logging anything.  So the
    warning appears on every run whatever happens, and a genuine NaN inside
    the fiducial region would be set to zero in silence, biasing the reference
    low.  This wrapper counts the non-finite entries per observable before the
    replacement and prints them with their (x, Q2), so every log says
    "[nan-audit] none" or names the points.  Installed once, idempotently,
    by make_cards -- every yadism run here (direct, charm, PineAPPL grid)
    builds its card through this function or yadism_cc_calc's, which calls
    it too.
    """
    import yadism.runner as yr
    if getattr(yr.Runner.replace_nans_with_0, "_bench_audit", False):
        return
    from yadism import observable_name
    orig = yr.Runner.replace_nans_with_0

    def audited(self, out):
        bad = []
        for name, points in out.items():
            if not observable_name.ObservableName.is_valid(name):
                continue
            for p in points:
                if p is None:
                    continue
                n = sum(int((~np.isfinite(v[i])).sum())
                        for v in p.orders.values() for i in range(2))
                if n:
                    bad.append((name, p.x, p.Q2, n))
        th = getattr(out, "theory", None) or {}
        ob = getattr(out, "observables", None) or {}
        tag = (f"FNS={th.get('FNS')} PTO={th.get('PTO')} TMC={th.get('TMC')} "
               f"prDIS={ob.get('prDIS')}")
        if not bad:
            print(f"[nan-audit] none: no non-finite coefficient ({tag})",
                  flush=True)
        else:
            print(f"[nan-audit] {len(bad)} point(s) with non-finite "
                  f"coefficients, set to zero by yadism ({tag}):", flush=True)
            for name, x, q2, n in bad:
                print(f"[nan-audit]   {name} x={x:.6g} Q2={q2:.6g} "
                      f"entries={n}", flush=True)
        return orig(self, out)
    audited._bench_audit = True
    yr.Runner.replace_nans_with_0 = audited


def make_cards(pto, x_nodes, q2_nodes):
    install_nan_audit()
    theory = {
        "ID": 0, "PTO": pto, "PTODIS": None, "FNS": "ZM-VFNS", "NfFF": 5,
        "DAMP": 0, "IC": 0, "IB": 0,
        "ModEv": "EXA", "ModSV": None, "XIR": 1.0, "XIF": 1.0,
        "fact_to_ren_scale_ratio": 1.0, "RenScaleVar": True, "FactScaleVar": True,
        "MaxNfAs": 5, "MaxNfPdf": 5,
        "Q0": 1.65, "nf0": 4, "alphas": 0.1180, "Qref": MZ, "nfref": 5,
        "QED": 0, "alphaqed": ALPHA, "Qedref": 1.777,
        "SIN2TW": S2TW, "MZ": MZ, "MW": 80.377, "GF": 1.1663787e-05,
        "CKM": ("0.97428 0.22530 0.003470 0.22520 0.97345 0.041000 "
                "0.00862 0.04030 0.999152"),
        "mc": 1.51, "Qmc": 1.51, "kcThr": 1.0,
        "mb": 4.92, "Qmb": 4.92, "kbThr": 1.0,
        "mt": 172.5, "Qmt": 172.5, "ktThr": 1.0,
        # TARGET-MASS CORRECTIONS from the one knob, target.tmc() ($BENCH_TMC,
        # default 0 = off, which is what every pre-result used).  Every
        # other card builder imports this one or yadism_cc_calc's, and every
        # cache key and filename downstream reads the same function.
        "MP": M_P, "TMC": target.tmc(), "n3lo_cf_variation": 0,
        "EScaleVar": 1,
    }
    kins = [{"x": float(x), "Q2": float(q2), "y": 0.5}
            for q2 in q2_nodes for x in x_nodes]
    xgrid = np.concatenate([np.geomspace(1e-4, 0.1, 40, endpoint=False),
                            np.linspace(0.1, 1.0, 41)]).tolist()
    observables = {
        "observables": {"F2_total": kins, "FL_total": kins, "F3_total": kins},
        "prDIS": "NC",
        "ProjectileDIS": "electron",   # mu-: identical NC couplings, mass not used
        "TargetDIS": "proton",
        "PolarizationDIS": 0.0,
        "PropagatorCorrection": 0.0,
        "NCPositivityCharge": None,
        "interpolation_xgrid": xgrid,
        "interpolation_polynomial_degree": 4,
        "interpolation_is_log": True,
    }
    return theory, observables


def sf_grid_nodes():
    """(x_nodes, q2_nodes) for the structure-function interpolation grid.

    THE GRID MUST FOLLOW THE BEAM.  These nodes were literals tuned to 1 TeV,
    where the fiducial region needs x >= 2.4e-3 and Q2 <= 1687 and so sits
    just inside them.  At 4 TeV it needs x >= 5.9e-4 and Q2 <= 6748, both
    OUTSIDE -- and RectBivariateSpline EXTRAPOLATES silently rather than
    raising, which cost 2.1% on sigma with no warning at all (the LO closure
    read 0.97876 instead of 1.00000; that closure is the only reason it was
    caught).

    x_min scales as 1/E and Q2_max as E, since both follow 2k.P.  The margins
    over the kinematic limits are therefore preserved exactly, and r = 1 at the
    anchor reproduces the original literals bit for bit.
    """
    r = EBEAM / 1000.0
    x_nodes = np.unique(np.concatenate([
        np.geomspace(0.002 / r, 0.5, 48, endpoint=False),
        1.0 - np.geomspace(0.5, 1.0 - XMAX, 18),
    ]))
    q2_nodes = np.geomspace(3.9, 1750.0 * r, 42)
    return x_nodes, q2_nodes


def structure_function_splines(pto, pdf):
    """Run yadism on an (x, Q2) grid, return cubic splines in (ln x, ln Q2)."""
    import yadism

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
    """sigma_red(Q2, y) and the differential cross sections (pb)."""
    X_LO = sf_grid_nodes()[0][0]        # the interpolation grid's lower edge

    def sf(name, x, q2):
        x, q2 = np.broadcast_arrays(np.asarray(x, float), np.asarray(q2, float))
        # CLIP TO THE GRID'S OWN LOWER EDGE, not a literal.  This was
        # np.clip(x, 0.002, XMAX) -- the 1 TeV grid's lower node.  At 1 TeV the
        # fiducial region only reaches x = 2.37e-3, just above it, so the clip
        # never fired; at 4 TeV x goes down to 5.9e-4 and everything below
        # 0.002 was evaluated AT 0.002 instead.  F2 rises as x falls, so the
        # clipped value is too small and sigma came out 2.1% low -- the whole
        # 4 TeV closure failure, and invisible at the anchor.
        v = splines[name].ev(np.log(np.clip(x, X_LO, XMAX)), np.log(q2))
        return np.where(x > XMAX, 0.0, v)

    def dsig_dxdq2(x, q2, y):
        """d2sigma/dxdQ2 in pb/GeV2 (x, q2, y broadcastable arrays/scalars)."""
        yp = 1.0 + (1.0 - y) ** 2
        ym = 1.0 - (1.0 - y) ** 2
        sred = (sf("F2", x, q2) - y * y / yp * sf("FL", x, q2)
                + ym / yp * sf("F3", x, q2))
        return 2.0 * math.pi * ALPHA**2 * yp / (x * q2 * q2) * sred * GEV2_PB

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
    """Lab muon scattering angle, exact with muon mass (analyze.py convention)."""
    ep = EBEAM * (1.0 - y)
    k = math.sqrt(EBEAM**2 - M_MU**2)
    kp = np.sqrt(ep * ep - M_MU**2)
    c = (2.0 * EBEAM * ep - 2.0 * M_MU**2 - q2) / (2.0 * k * kp)
    return np.arccos(np.clip(c, -1.0, 1.0))


# ----------------------------------------------------------------- integrators

# ------------------------------------------------------------ the W^2 bound
#
# W^2 = m_p^2 + 2 P.q - Q^2 and 2 P.q = y * (2 k.P), so the hadronic-mass cut
# is EXACTLY a bound on the (Q2, y) rectangle the integrators already work in:
#
#     W^2 > W2_MIN   <=>   y > (W2_MIN - m_p^2 + Q2) / (2 k.P)
#                    <=>   Q2 < m_p^2 + y (2 k.P) - W2_MIN
#
# -- the same shape as the x < XMAX kinematic limit already imposed here, so
# it costs a term in a max() rather than a change of variables.  W2_MIN is
# None for every selection the benchmark publishes; it is set only by the
# foreign arXiv:2506.13889 region (see analysis/selection.py).
def y_min_from_w2(q2, two_kp):
    """Lower y limit imposed by W2 > W2_MIN at this Q2 (-inf if no cut).

    two_kp IS PASSED IN, not read from this module.  The neutrino integrators
    import these helpers, and a helper that closed over the muon module's
    TWO_KP would silently apply the wrong beam invariant on the other current
    -- a cross-module coupling of exactly the kind that has gone unnoticed
    here before.
    """
    if W2_MIN is None:
        return -np.inf
    return (W2_MIN - M_P * M_P + q2) / two_kp


def q2_max_from_w2(y, two_kp):
    """Upper Q2 limit imposed by W2 > W2_MIN at this y (+inf if no cut)."""
    if W2_MIN is None:
        return np.inf
    return M_P * M_P + y * two_kp - W2_MIN


def q2_min_from_w2(x):
    """Lower Q2 limit imposed by W2 > W2_MIN at this x (-inf if no cut).

    At fixed x, W^2 = m_p^2 + Q^2 (1-x)/x, so the cut becomes a FLOOR on Q2
    that rises with x -- it is the large-x corner it removes.
    """
    if W2_MIN is None:
        return -np.inf
    if x >= 1.0:
        return np.inf
    return (W2_MIN - M_P * M_P) * x / (1.0 - x)


def y_lattice(y_lo, y_hi):
    """The 0.005 y lattice of the integrators, ENDING AT y_hi.

    >>> IT DID NOT END THERE WHEN y_hi IS OFF THE LATTICE (found 2026-09-13).
    np.arange(Y_MIN, Y_HI + 1e-9, 0.005) stops at the last multiple of 0.005
    below Y_HI, and the panels were built from that array alone -- so with
    no y window, where Y_HI is the kinematic ceiling 1 - m_mu/E (0.99989 at
    1 TeV), the strip 0.995 < y < Y_HI was never integrated.  Nothing failed:
    the q4w3 LO closure read 0.99549 (CC, 400 GeV), 0.99610 (CC, 4 TeV) and
    0.99765 (NC, 4 TeV), and the x-Q2 cross-check, which has no lattice,
    agreed with the direct LHAPDF integral to 1e-4.  Every no-y-cut region
    (ally, dis2506, sidis_*, q4w3) was short by that strip.  With a y window
    ending on the lattice (0.9) the appended edge is not added, so every
    published inclusive-region number is bit-for-bit unchanged.
    """
    lat = np.round(np.arange(y_lo, y_hi + 1e-9, 0.005), 10)
    if y_hi - lat[-1] > 1e-9:
        lat = np.append(lat, y_hi)
    return lat


def integrate_q2y(dsig_dq2dy, hist_accum):
    """Master integration on (Q2, y): sigma_fid, means, Q2/y/Emu/nu hists.

    Q2 panels aligned to the Q2 bin edges (+ overflow up to the kinematic
    limit); y panels on a 0.005 lattice, which contains the bin edges of
    y (0.025), Emu and nu (both 0.02 in y) exactly.
    """
    q2n, q2w, q2i = gl_panels(BINS["Q2"], 3, 8, log=True)
    q2o, q2ow = gl_interval(BINS["Q2"][-1], Q2_KIN_MAX, 8, log=True)
    q2_all = np.concatenate([q2n, q2o])
    q2w_all = np.concatenate([q2w, q2ow])
    q2i_all = np.concatenate([q2i, np.full(len(q2o), -1)])  # -1: Q2 overflow

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
        # clip the y lattice
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
    k = math.sqrt(EBEAM**2 - M_MU**2)
    sigma = 0.0
    for y, wy in zip(yn, yw):
        ep = EBEAM * (1.0 - y)
        kp = math.sqrt(ep * ep - M_MU**2)
        q2max = min(y * TWO_KP * XMAX, q2_max_from_w2(y, TWO_KP),
                    2.0 * EBEAM * ep - 2.0 * M_MU**2
                    - 2.0 * k * kp * (-1.0))
        thmin = float(theta_of(Q2_MIN, y))
        thmax = float(theta_of(q2max, y))
        for i in range(len(th_edges) - 1):
            a, b = max(th_edges[i], thmin), min(th_edges[i + 1], thmax)
            if b <= a:
                continue
            tt, wt = gl_interval(a, b, 6, log=False)
            q2 = 2.0 * (EBEAM * ep - k * kp * np.cos(tt)) - 2.0 * M_MU**2
            q2 = np.clip(q2, Q2_MIN, q2max)
            jac = 2.0 * k * kp * np.sin(tt)         # dQ2/dtheta
            f = dsig_dq2dy(q2, y)
            s = wy * (wt * f * jac).sum()
            hist_accum["theta"][i] += s
            sigma += s
    return sigma


# -------------------------------------------------------- LO closure check
def manual_lo_sigma(pdf, em_only=False):
    """Direct (spline-free) LO parton-model sigma_fid from LHAPDF, in pb."""
    def kappa(q2):
        return q2 / (q2 + MZ * MZ) / (4.0 * S2TW * (1.0 - S2TW))

    el, vl, al = -1.0, -0.5 + 2.0 * S2TW, -0.5
    w2c, w3c = {}, {}
    for pid in (1, 2, 3, 4, 5):
        eq = 2.0 / 3.0 if pid % 2 == 0 else -1.0 / 3.0
        t3 = 0.5 if pid % 2 == 0 else -0.5
        vq, aq = t3 - 2.0 * eq * S2TW, t3
        w2c[pid] = (eq, 2.0 * el * eq * vl * vq, (vl**2 + al**2) * (vq**2 + aq**2))
        w3c[pid] = (2.0 * el * al * eq * aq, 2.0 * vl * al * 2.0 * vq * aq)

    q2n, q2w = gl_panels(np.geomspace(Q2_MIN, Q2_KIN_MAX, 40), 1, 6, log=True)[:2]
    tg, wg = np.polynomial.legendre.leggauss(6)
    sigma = 0.0
    for q2, wq in zip(q2n, q2w):
        kap = kappa(q2)
        ymin_eff = max(Y_MIN, q2 / (XMAX * TWO_KP), y_min_from_w2(q2, TWO_KP))
        if ymin_eff >= Y_HI:
            continue
        sub = np.linspace(ymin_eff, Y_HI, 31)
        yy = (0.5 * (sub[:-1] + sub[1:])[:, None]
              + 0.5 * np.diff(sub)[:, None] * tg).ravel()
        wy = (0.5 * np.diff(sub)[:, None] * wg).ravel()
        xx = q2 / (yy * TWO_KP)
        f2 = np.zeros(len(yy))
        xf3 = np.zeros(len(yy))
        for pid, (eq, c1, c2) in w2c.items():
            a2 = eq * eq + c1 * kap + c2 * kap * kap
            if em_only:
                a2 = eq * eq
            a3 = 0.0 if em_only else w3c[pid][0] * kap + w3c[pid][1] * kap * kap
            qv = np.array([pdf.xfxQ2(pid, x, q2) if x < 1 else 0.0 for x in xx])
            qbv = np.array([pdf.xfxQ2(-pid, x, q2) if x < 1 else 0.0 for x in xx])
            f2 += a2 * (qv + qbv)
            xf3 += a3 * (qv - qbv)
        yp = 1.0 + (1.0 - yy) ** 2
        ym = 1.0 - (1.0 - yy) ** 2
        dxdy = (2.0 * math.pi * ALPHA**2 * yp / (xx * q2 * q2)
                * (f2 + ym / yp * xf3) * GEV2_PB)
        sigma += wq * (wy * dxdy * xx / yy).sum()
    return sigma


# ------------------------------------------------------------------------ main
def main():
    order = (sys.argv[1] if len(sys.argv) > 1 else "lo").lower()
    if order not in ORDERS:
        sys.exit(f"unknown order '{order}' (use lo|nlo|nnlo)")
    pto = ORDERS[order]

    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDFSET, 0)

    print(f"[yadism_calc] running yadism at {order.upper()} (PTO={pto}) ...")
    splines = structure_function_splines(pto, pdf)
    dsig_dxdq2, dsig_dq2dy = make_sigma_red(splines)

    hists = {}
    sigma1, means = integrate_q2y(dsig_dq2dy, hists)
    sigma2 = integrate_xq2(dsig_dxdq2, hists)
    sigma3 = integrate_ytheta(dsig_dq2dy, hists)
    print(f"[yadism_calc] sigma_fid = {sigma1/1e3:.4f} nb "
          f"(cross-checks: x-Q2 {sigma2/1e3:.4f}, y-theta {sigma3/1e3:.4f})")

    closure = None
    # The direct LHAPDF integral is the bare parton model: with target-mass
    # corrections on, a "failed" closure would be the TMC, not a bug.
    if order == "lo" and target.tmc() != 0:
        print("[yadism_calc] LO closure SKIPPED: TMC is on, and the direct "
              "parton-model integral carries none")
    if order == "lo" and target.tmc() == 0:
        direct = manual_lo_sigma(pdf)
        gamma = manual_lo_sigma(pdf, em_only=True)
        closure = sigma1 / direct
        print(f"[yadism_calc] LO closure: direct LHAPDF NC = {direct/1e3:.4f} nb "
              f"(yadism/direct = {closure:.5f}), gamma-only = "
              f"{gamma/1e3:.4f} nb")
        # The closure is the ONLY thing that catches a wrong reference: yadism
        # returns a perfectly plausible number either way.  Say so loudly and
        # record it, rather than leaving one line in a log to be scrolled past.
        if abs(closure - 1.0) > CLOSURE_TOL:
            print("*" * 72)
            print(f"[yadism_calc] LO CLOSURE FAILS at E = {EBEAM:g} GeV: "
                  f"yadism/direct = {closure:.5f}")
            print("  The spline-free parton-model integral and yadism disagree by "
                  f"{100*abs(closure-1):.2f}%.")
            print("  DO NOT USE THIS AS A VALIDATED REFERENCE; the result is")
            print("  written with \"lo_closure_ok\": false so downstream code can see it.")
            print("*" * 72)

    # numerical error: spread of the three independent integrations + 2e-4 floor
    err = max(abs(sigma1 - sigma2), abs(sigma1 - sigma3), 2e-4 * sigma1)

    gen = "yadism" if order == "lo" else f"yadism_{order}"
    # SELECTION BEFORE ENERGY, which is the order analyze.result_path() uses
    # ("histos_<gen>[_<selection>][_<energy>]") and therefore the order
    # faser_rates._result() looks in.  Building it the other way round gave
    # histos_yadism_nlo_400GeV_dis2506.json against the generators'
    # histos_powheg_dis2506_400GeV.json, and one lookup rule cannot find both.
    # sel_suffix() is "" for the inclusive selection, so every existing file
    # keeps its name.
    # _tmc after the target (target.tmc_suffix(), "" when TMC is off), so the
    # lookup rule is histos_<gen>[_<sel>][_<target>][_tmc][_<energy>].
    gen = beams_mod.at_energy(gen + sel_suffix() + target.suffix()
                              + target.tmc_suffix(), EBEAM)
    out = {
        "generator": gen,
        "label": f"YADISM {order.upper()}",
        "n_parsed": None,
        "n_fiducial": None,
        "sigma_fid_pb": sigma1,
        "sigma_fid_err_pb": err,
        "means": {k: float(v) for k, v in means.items()},
        "lo_closure": None if closure is None else float(closure),
        "lo_closure_ok": (None if closure is None
                          else bool(abs(closure - 1.0) <= CLOSURE_TOL)),
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

    # THE SELECTION MUST BE IN THE NAME.  It was not, and the consequence is
    # the anchor-overwrite trap this repo has sprung before: running this
    # module under any non-inclusive selection wrote its answer straight over
    # results/histos_yadism_nlo.json -- the inclusive 1 TeV reference that the
    # whole benchmark is normalised against -- with no warning, because the
    # energy tag and the Q2 floor were in the filename and the selection was
    # not.  Caught on 2026-09-03 by a dis2506 test run that replaced 35.1 nb
    # with 173.2 nb; git had the file, so nothing was lost.
    #
    # sel_suffix() returns "" for the inclusive selection, so every existing
    # file keeps exactly the name it has.
    # The target joins the name BEFORE the energy tag, for the same reason
    # the selection does and in the same place: the lookup rule is
    # "histos_<gen>[_<selection>][_<target>][_<energy>]" and one rule has to
    # find every file.  target.suffix() is "" for the proton, so nothing that
    # exists moves.  See analysis/target.py.
    ofn = f"{BASE}/results/histos_{gen}{selection.q2_suffix()}.json"
    stamp(out, current="mu")
    prov = target.provenance(PDFSET)
    if prov:                  # absent at proton / TMC off: earlier stays byte-identical
        out.update(prov)
    write_result(ofn, out)
    print(f"[yadism_calc] wrote {ofn}")


if __name__ == "__main__":
    main()
