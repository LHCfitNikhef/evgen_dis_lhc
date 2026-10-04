#!/usr/bin/env python3
"""The impact of NUCLEAR PDFs on the inclusive DIS distributions, region.

The redo of the earlier study (analysis/npdf_impact.py, deleted with every the earlier production
result in ea61871) for the final benchmark region, and the input to appendix
paper plot A2 (analysis/paper_plots/ppA2_npdf_impact.py).

WHAT CHANGED FROM the earlier production, AND WHY.  earlier reweighted POWHEG-V2 events member by
member.  Those samples are gone, and the muon NLO generator (POWHEG-RES)
has no reweighting, so a POWHEG route would cover one current only -- the
asymmetry CONVENTIONS.md rule 2b forbids.  This version is ANALYTIC: YADISM at NLO
in ZM-VFNS with target-mass corrections (mode 3), exactly the calculation of
the references histos_yadism_nlo_q4w3_W_tmc.json, evaluated from the stored
PineAPPL grids (analysis/pineappl_grids.py) with each PDF member in turn.  A
grid is PDF- and region-independent, so a member costs one convolution plus
the fiducial integration, about half a second.  Every observable here is
leptonic, so nothing is lost by having no hadrons.

THE CALCULATION IS THE REFERENCE'S OWN, NOT A COPY OF IT.  The splines, the
dsigma formula and the fiducial integration are the calculator modules'
functions (make_sigma_red, integrate_q2y, integrate_xq2), the same path
mhou_diff.py takes.  CLOSURE: member 0 of the baseline must reproduce the
tracked reference -- sigma_fid and every x, Q2 and y bin -- to 1e-6, or the
run dies.  Measured: 4e-16.

CONVERGENCE, so that no bin-to-bin wiggle can be numerical (user,
2026-09-21, against the noisy reweighted earlier figure).  Two tests:
  * the fiducial QUADRATURE -- every run re-integrates member 0 of each band
    set with four more Gauss-Legendre points per rule, every panel halved and
    the y lattice halved (quadrature_check, stored in the JSON): 3e-5 per bin,
    1e-6 on the ratios;
  * the STRUCTURE-FUNCTION NODES -- done once, off-line, on 2026-09-21: grids
    rebuilt on a node grid with a midpoint (in log) inserted between every
    pair of x and Q2 nodes moved sigma_fid by < 1e-4, any bin holding 0.3%
    of the rate by < 2e-3 and any RATIO to the baseline by < 1.3e-3 (the
    worst being EPPS21 at the largest x), against nuclear effects of 3-40%.
    Not repeated per run: the dense grids take ~20 min and 330 MB.

THE REGION: Q2 > 4 GeV2 and W > 3 GeV, no y cut, 1 TeV, per nucleon of
tungsten (74 p + 110 n).  Fixed here rather than read from the environment,
because this study exists for one region; the BENCH_* knobs are set
unconditionally before any benchmark module is imported.

THE SETS.
  * baseline: NNPDF40_nnlo_as_01180_W184free, the free-nucleon average of
    tungsten (tools/make_isoscalar_pdf.py) -- nuclear effects assumed to
    vanish, which is what every reference uses; 100 replicas.
  * nNNPDF30_nlo_as_0118_A184_Z74 (200 MC replicas), EPPS21nlo_CT18Anlo_W184
    (Hessian, 106 eigenvector members, 90% CL) and nCTEQ15HQ_FullNuc_184_74
    (symmetric Hessian, 38 members, 90% CL; arXiv:2204.09982).  nCTEQ15HQ, not
    nCTEQ15 (user, 2026-09-30): the version with the heavy-quark data, which
    is also the one arXiv:2509.00144 uses.  These are ALREADY the average bound
    nucleon of tungsten, not the bound proton: no isospin is applied on top.
    That is CHECKED on every run with the valence sum rules (below).
  * each nuclear set's own free-nucleon baseline, central member only
    (nNNPDF30_nlo_as_0118_p_W184free, CT18ANLO_W184free,
    nCTEQ15HQ_1_1_W184free), so the nuclear modification proper can be
    separated from the difference between proton fits.
  * the free PROTON (NNPDF40_nnlo_as_01180), central only, recorded as the
    isospin-only line, as earlier did.  Not drawn.

UNCERTAINTIES through LHAPDF's own prescription, at 68% CL, exactly as the earlier production:
set.uncertainty(values, 68.27) knows each set's error type and confidence
level, so EPPS21's and nCTEQ15's 90% CL Hessian eigenvectors are rescaled to
68% (asymmetric, err_plus/err_minus) and the replica sets give their standard
deviation.  Writing the 1.645 by hand is how a 90% band ends up drawn as 68%.

alpha_s is each set's own (alphasQ2 of the member), as in every reference.

Usage:  analysis/npdf_impact.py [--current mu|nu] [--jobs N]
Writes: results/npdf_impact_q4w3_W_tmc.json (muon NC)
        results_nu/npdf_impact_q4w3_W_tmc.json (neutrino CC)
"""
import math
import os
import sys
import time

# THE REGION, fixed.  Unconditional assignment: a default on these would let a
# stray $BENCH_SELECTION from the calling shell compute another region under
# this file's name (CONVENTIONS.md rule 1).
os.environ["BENCH_SELECTION"] = "q4w3"
os.environ["BENCH_TARGET"] = "W"
os.environ["BENCH_TMC"] = "3"
os.environ["BENCH_ENERGY"] = "1000"

import numpy as np                                            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

PTO = 1
ENERGY = 1000.0
OBS = ("xbj", "Q2", "y")
NAMES = ["F2_total", "FL_total", "F3_total"]
CL68 = 100.0 * math.erf(1.0 / math.sqrt(2.0))
CLOSURE_TOL = 1e-6

BASELINE = "NNPDF40_nnlo_as_01180_W184free"
PROTON = "NNPDF40_nnlo_as_01180"
NUCLEAR = ["nNNPDF30_nlo_as_0118_A184_Z74", "EPPS21nlo_CT18Anlo_W184",
           "nCTEQ15HQ_FullNuc_184_74"]
OWN_BASELINE = {"nNNPDF30_nlo_as_0118_A184_Z74":
                "nNNPDF30_nlo_as_0118_p_W184free",
                "EPPS21nlo_CT18Anlo_W184": "CT18ANLO_W184free",
                "nCTEQ15HQ_FullNuc_184_74": "nCTEQ15HQ_1_1_W184free"}
# the sets whose full member list is evaluated; the rest are central only
WITH_BAND = [BASELINE] + NUCLEAR
CENTRAL_ONLY = [PROTON] + list(OWN_BASELINE.values())

REF = {"mu": f"{BASE}/results/histos_yadism_nlo_q4w3_W_tmc.json",
       "nu": f"{BASE}/results_nu/histos_yadism_nlo_q4w3_W_tmc.json"}
OUT = {"mu": f"{BASE}/results/npdf_impact_q4w3_W_tmc.json",
       "nu": f"{BASE}/results_nu/npdf_impact_q4w3_W_tmc.json"}


def die(msg):
    sys.exit(f"npdf_impact: {msg}")


# --------------------------------------------------------------- per member
_CACHE = {}


def _setup(current):
    """(calc module, grids, x nodes, Q2 nodes) for one current, cached."""
    if current in _CACHE:
        return _CACHE[current]
    import importlib
    import lhapdf
    lhapdf.setVerbosity(0)
    import target                                             # noqa: F401
    import pineappl_grids as pg
    from pineappl.grid import Grid
    calc = importlib.import_module("yadism_cc_calc" if current == "nu"
                                   else "yadism_calc")
    from yadism_calc import sf_grid_nodes
    xn, qn = sf_grid_nodes()
    key = pg.grid_key(current, ENERGY, PTO, xn, qn, "ZM-VFNS", 5)
    if not key.endswith("_tmc3"):
        die(f"grid key {key} is not the TMC-3 grid")
    grids = {}
    for n in NAMES:
        f = pg.grid_file(key, n)
        if not os.path.exists(f):
            die(f"no grid {os.path.relpath(f, BASE)}: build it with\n  "
                f"BENCH_TMC=3 analysis/pineappl_grids.py build {current} "
                f"1000 1")
        grids[n] = Grid.read(f)
    _CACHE[current] = (calc, grids, xn, qn, key)
    return _CACHE[current]


def evaluate(args):
    """sigma_fid, mean x, and the x, Q2, y histograms for one PDF member.

    The structure-function values come from the grid (pineappl_grids.convolve
    with the Grid read once per process instead of once per call); from there
    on it is the calculator module's own path.
    """
    current, setname, member = args
    import lhapdf
    from scipy.interpolate import RectBivariateSpline
    lhapdf.setVerbosity(0)
    calc, grids, xn, qn, _key = _setup(current)
    pdf = lhapdf.mkPDF(setname, member)
    sp = {}
    for n in NAMES:
        g = grids[n]
        v = np.asarray(g.convolve(
            pdg_convs=g.convolutions,
            xfxs=[lambda pid, x, q2: pdf.xfxQ2(pid, x, q2)],
            alphas=lambda q2: pdf.alphasQ2(q2), xi=None))
        sp[n[:2]] = RectBivariateSpline(
            np.log(xn), np.log(qn), v.reshape(len(qn), len(xn)).T,
            kx=3, ky=3)
    dsig_dxdq2, dsig_dq2dy = calc.make_sigma_red(sp)
    h = {}
    sig, means = calc.integrate_q2y(dsig_dq2dy, h)
    calc.integrate_xq2(dsig_dxdq2, h)
    return (setname, member, float(sig), float(means["xbj"]),
            float(means["Q2"]), {k: np.asarray(h[k], float) for k in OBS})


# ------------------------------------------------- quadrature convergence
QUAD_SETS = [BASELINE] + NUCLEAR
QUAD_MIN_SHARE = 1e-3


class _FineQuadrature:
    """Temporarily refine every quadrature the fiducial integration uses:
    each Gauss-Legendre rule gets four more points, every panel of
    gl_panels is split in two, and the 0.005 y lattice is halved.  Applied
    in THIS process only (the pool workers never see it), and undone on exit.
    """

    def __init__(self, calc):
        self.calc = calc

    def __enter__(self):
        import numpy.polynomial.legendre as leg
        c = self.calc
        self.saved = (leg.leggauss, c.gl_panels, c.y_lattice)
        orig_lg, orig_yl = leg.leggauss, c.y_lattice
        leg.leggauss = lambda n: orig_lg(n + 4)

        def gl_panels(edges, nsub, order, log=False):
            t, w = leg.leggauss(order)
            nodes, wgts, ibin = [], [], []
            for i in range(len(edges) - 1):
                a, b = edges[i], edges[i + 1]
                subs = (np.geomspace(a, b, 2 * nsub + 1) if log
                        else np.linspace(a, b, 2 * nsub + 1))
                for lo, hi in zip(subs[:-1], subs[1:]):
                    if log:
                        la, lb = math.log(lo), math.log(hi)
                        xx = np.exp(0.5 * (la + lb) + 0.5 * (lb - la) * t)
                        ww = 0.5 * (lb - la) * w * xx
                    else:
                        xx = 0.5 * (lo + hi) + 0.5 * (hi - lo) * t
                        ww = 0.5 * (hi - lo) * w
                    nodes.append(xx)
                    wgts.append(ww)
                    ibin.append(np.full(len(t), i))
            return (np.concatenate(nodes), np.concatenate(wgts),
                    np.concatenate(ibin).astype(int))

        def y_lattice(a, b):
            lat = orig_yl(a, b)
            return np.sort(np.concatenate([lat, 0.5 * (lat[:-1] + lat[1:])]))

        c.gl_panels, c.y_lattice = gl_panels, y_lattice
        return self

    def __exit__(self, *exc):
        import numpy.polynomial.legendre as leg
        leg.leggauss, self.calc.gl_panels, self.calc.y_lattice = self.saved
        return False


def quadrature_check(current, by):
    """Is the fiducial integration converged in every bin that matters?

    Member 0 of the baseline and of each nuclear set is re-integrated with the
    refined quadrature and compared bin by bin with the production numbers,
    over bins holding at least 0.1% of the baseline rate: both the absolute
    bin content and the RATIO to the baseline, which is what the figure
    draws.  Measured 2026-09-21: bins to 3e-5, ratios to 1e-6 -- far below
    the smallest nuclear effect shown.
    """
    calc = _setup(current)[0]
    fine = {}
    with _FineQuadrature(calc):
        for name in QUAD_SETS:
            fine[name] = evaluate((current, name, 0))
    base_c = by[BASELINE][0][3]
    share = {k: base_c[k] / base_c[k].sum() for k in OBS}
    out = {"refinement": ("Gauss-Legendre +4 points, panels split in two, "
                          "y lattice halved"),
           "min_share": QUAD_MIN_SHARE, "bin": {}, "ratio": {}, "sigma": {}}
    for name in QUAD_SETS:
        cs, ch = by[name][0][0], by[name][0][3]
        fs, fh = fine[name][2], fine[name][5]
        out["sigma"][name] = abs(fs / cs - 1.0)
        out["bin"][name] = max(float(np.max(np.abs(
            fh[k][share[k] >= QUAD_MIN_SHARE]
            / ch[k][share[k] >= QUAD_MIN_SHARE] - 1.0))) for k in OBS)
        out["ratio"][name] = max(float(np.max(np.abs(
            (fh[k] / fine[BASELINE][5][k])[share[k] >= QUAD_MIN_SHARE]
            / (ch[k] / base_c[k])[share[k] >= QUAD_MIN_SHARE] - 1.0)))
            for k in OBS)
    print(f"  quadrature: worst bin {max(out['bin'].values()):.1e}, worst "
          f"ratio {max(out['ratio'].values()):.1e} (refined vs production)")
    if max(out["bin"].values()) > 1e-3:
        die("the fiducial integration is not converged to 1e-3 in some bin")
    return out


# ------------------------------------------------------- target composition
def valence_integrals(name):
    """(int u_v dx, int d_v dx) of member 0 at Q = 10 GeV."""
    import lhapdf
    lhapdf.setVerbosity(0)
    p = lhapdf.mkPDF(name, 0)
    xs = np.logspace(-6, 0, 3000)
    return tuple(float(np.trapezoid(np.array(
        [(p.xfxQ(pid, x, 10.0) - p.xfxQ(-pid, x, 10.0)) / x for x in xs]),
        xs)) for pid in (2, 1))


def check_target_composition():
    """Every set must describe the average nucleon of tungsten -- except the
    free proton, which is recorded as the isospin-only line.

    THE TEST IS THE VALENCE SUM RULE, which pins the composition exactly:
    (2, 1) for u_v, d_v in a proton grid, ((2Z+N)/A, (Z+2N)/A) = (1.402,
    1.598) per nucleon of tungsten.  If a future release shipped bound-PROTON
    grids under the same name, the isospin mix would be missing and the
    neutrino rate wrong by ten per cent with every curve looking fine.  The
    tolerance is 3% (the replica average of an MC set does not satisfy the
    sum rule exactly, and the quadrature stops at x = 1e-6); the two
    hypotheses are 0.6 apart.
    """
    import target
    want_w = ((2 * target.Z_W + target.N_W) / target.A_W,
              (target.Z_W + 2 * target.N_W) / target.A_W)
    want_p = (2.0, 1.0)
    out = {}
    for name in WITH_BAND + CENTRAL_ONLY:
        got = valence_integrals(name)
        dw = max(abs(g - w) for g, w in zip(got, want_w))
        dp = max(abs(g - w) for g, w in zip(got, want_p))
        kind = "W" if dw < dp else "p"
        expect = "p" if name == PROTON else "W"
        if kind != expect or min(dw, dp) > 0.03 * 2.0:
            die(f"{name}: valence integrals u_v = {got[0]:.3f}, d_v = "
                f"{got[1]:.3f} -- expected the "
                f"{'free proton' if expect == 'p' else 'average nucleon of W'}"
                f" ({want_p if expect == 'p' else want_w}).")
        out[name] = {"u_v": got[0], "d_v": got[1], "is": kind,
                     "deviation": min(dw, dp)}
        print(f"  {name:36s} u_v {got[0]:.4f}  d_v {got[1]:.4f}  -> "
              f"{'free proton' if kind == 'p' else 'average nucleon of W'}")
    return {"expected_W": list(want_w), "expected_p": list(want_p),
            "sets": out}


# ------------------------------------------------------------------ bands
def band(pset, values):
    """(central, err_plus, err_minus) at 68% CL through LHAPDF."""
    if len(values) < 2:
        return values[0], 0.0, 0.0
    u = pset.uncertainty(list(values), CL68)
    return u.central, u.errplus, u.errminus


def run(current, jobs, composition):
    import lhapdf
    import json
    lhapdf.setVerbosity(0)
    import beams
    import selection
    import target
    from analyze import BINS
    from runmeta import write_result

    print(f"\n=== {current}: nuclear PDFs, YADISM NLO ZM-VFNS, TMC 3, q4w3, "
          f"W per nucleon, 1 TeV ===")
    sizes = {n: lhapdf.getPDFSet(n).size for n in WITH_BAND}
    tasks = [(current, n, m) for n in WITH_BAND for m in range(sizes[n])]
    tasks += [(current, n, 0) for n in CENTRAL_ONLY]
    t0 = time.time()
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(jobs) as pool:
            rows = pool.map(evaluate, tasks, chunksize=4)
    else:
        rows = [evaluate(t) for t in tasks]
    print(f"  {len(tasks)} members in {time.time()-t0:.0f} s "
          f"({jobs} process{'es' if jobs > 1 else ''})")
    by = {}
    for setname, m, sig, mx, mq2, h in rows:
        by.setdefault(setname, {})[m] = (sig, mx, mq2, h)

    edges = {k: np.asarray(BINS[k], float) for k in OBS}
    widths = {k: np.diff(edges[k]) for k in OBS}

    # ---- closure: member 0 of the baseline IS the tracked reference
    with open(REF[current]) as f:
        ref = json.load(f)
    if (ref.get("pdf_set") != BASELINE or ref.get("tmc") != 3
            or ref.get("target") != "W"
            or ref["selection"]["name"] != "q4w3"):
        die(f"{REF[current]} is not the q4w3 / W / TMC-3 reference")
    s0, _mx, _mq, h0 = by[BASELINE][0]
    dev = {"sigma_fid": abs(s0 / ref["sigma_fid_pb"] - 1.0)}
    for k in OBS:
        r = np.asarray(ref["hists"][k]["dsig"], float)
        if not np.allclose(ref["hists"][k]["edges"], edges[k]):
            die(f"the {k} binning differs from the reference's")
        m = r != 0
        dev[k] = float(np.max(np.abs(h0[k][m] / widths[k][m] / r[m] - 1.0)))
    worst = max(dev.values())
    print(f"  closure vs {os.path.relpath(REF[current], BASE)}: worst "
          f"{worst:.2e}  ({', '.join(f'{k} {v:.1e}' for k, v in dev.items())})")
    if worst > CLOSURE_TOL:
        die(f"member 0 of {BASELINE} does not reproduce the tracked "
            f"reference ({worst:.2e} > {CLOSURE_TOL:g}): this is not the "
            f"reference's calculation, and a few-per-cent nuclear effect on "
            f"top of it would look fine.")

    quad = quadrature_check(current, by)

    # ---- bands
    sets_out = {}
    for name in WITH_BAND + CENTRAL_ONLY:
        pset = lhapdf.getPDFSet(name)
        mem = by[name]
        idx = sorted(mem)
        c, ep, em = band(pset, [mem[i][0] for i in idx])
        rec = {"members": len(idx), "error_type": pset.errorType,
               "conf_level_native": pset.errorConfLevel,
               "integrated": {"sigma_fid_pb": c, "err_plus_pb": ep,
                              "err_minus_pb": em,
                              "member0_pb": mem[0][0],
                              "mean_x_member0": mem[0][1],
                              "mean_q2_member0": mem[0][2]},
               "observables": {}}
        for k in OBS:
            cen, up, dn = [], [], []
            for j in range(len(widths[k])):
                a, b_, cc = band(pset, [mem[i][3][k][j] for i in idx])
                cen.append(a / widths[k][j])
                up.append(b_ / widths[k][j])
                dn.append(cc / widths[k][j])
            rec["observables"][k] = {"dsig": cen, "err_plus": up,
                                     "err_minus": dn,
                                     "dsig_member0":
                                     (mem[0][3][k] / widths[k]).tolist()}
        sets_out[name] = rec

    base = sets_out[BASELINE]["integrated"]["sigma_fid_pb"]
    ratios = {}
    print(f"  sigma_fid, 68% CL, and the ratio to {BASELINE}:")
    for name in WITH_BAND + CENTRAL_ONLY:
        v = sets_out[name]["integrated"]
        ratios[name] = {"to_baseline": v["sigma_fid_pb"] / base,
                        "rel_plus": v["err_plus_pb"] / v["sigma_fid_pb"],
                        "rel_minus": v["err_minus_pb"] / v["sigma_fid_pb"]}
        print(f"    {name:36s} {v['sigma_fid_pb']:12.6g} pb  "
              f"+{100*ratios[name]['rel_plus']:.2f}% "
              f"-{100*ratios[name]['rel_minus']:.2f}%   ratio "
              f"{ratios[name]['to_baseline']:.4f}   <x> "
              f"{v['mean_x_member0']:.4f}")
    for name in NUCLEAR:
        ob = sets_out[OWN_BASELINE[name]]["integrated"]["sigma_fid_pb"]
        rr = sets_out[name]["integrated"]["sigma_fid_pb"] / ob
        ratios[name]["own_baseline"] = OWN_BASELINE[name]
        ratios[name]["to_own_baseline"] = rr
        print(f"    {name:36s} nuclear modification alone: "
              f"{100*(rr-1):+.2f}%")

    sel = selection.get()
    out = {
        "what": ("impact of nuclear PDFs on the inclusive DIS distributions: "
                 "YADISM NLO in ZM-VFNS with target-mass corrections, "
                 "evaluated per PDF member from the PineAPPL grids"),
        "current": current,
        "process": "nu_mu CC" if current == "nu" else "mu- NC",
        "calculation": {"code": "YADISM", "order": "NLO", "pto": PTO,
                        "scheme": "ZM-VFNS", "tmc": target.tmc(),
                        "alphas": "each member's own",
                        "grid_key": _setup(current)[4]},
        "energy_gev": ENERGY, "beam_tag": beams.Beams(current, ENERGY).tag,
        "selection": sel.as_dict(),
        "target": {"name": "W", "Z": target.Z_W, "N": target.N_W,
                   "A": target.A_W, "per": "nucleon"},
        "baseline_set": BASELINE, "proton_set": PROTON,
        "nuclear_sets": NUCLEAR, "own_baseline": OWN_BASELINE,
        "confidence_level_percent": CL68,
        "band_prescription": ("LHAPDF PDFSet.uncertainty at 68.27% CL: "
                              "Hessian 90% CL sets rescaled, asymmetric; "
                              "replica sets their standard deviation"),
        "closure": {"reference": os.path.relpath(REF[current], BASE),
                    "max_rel_dev": worst, "per_observable": dev,
                    "tolerance": CLOSURE_TOL},
        "target_composition": composition,
        "quadrature_check": quad,
        "binning": {k: {"edges": edges[k].tolist(),
                        "log": k in ("xbj", "Q2")} for k in OBS},
        "sets": sets_out, "ratios": ratios,
    }
    write_result(OUT[current], out)
    print(f"  wrote {os.path.relpath(OUT[current], BASE)}")


def main():
    argv = sys.argv[1:]
    currents = ["mu", "nu"]
    if "--current" in argv:
        currents = [argv[argv.index("--current") + 1]]
        if currents[0] not in ("mu", "nu"):
            die("--current takes mu or nu")
    jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 6
    print("valence sum rules (the target each set describes):")
    composition = check_target_composition()
    for c in currents:
        run(c, jobs, composition)


if __name__ == "__main__":
    main()
