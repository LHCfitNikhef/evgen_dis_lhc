#!/usr/bin/env python3
"""YADISM analytic LO CHARM-PRODUCTION reference for the CC NEUTRINO benchmark.

CC counterpart of yadism_charm_calc.py.  In the ZM-VFNS at LO, charm
production in CC neutrino DIS is CKM-suppressed quark flavour conversion,

    W+ + d/s/b -> c        with weights |Vcd|^2, |Vcs|^2, |Vcb|^2,

i.e. the produced-charm piece of the CC structure functions.  Dominated by
s -> c off the strange sea (|Vcs|^2 = 0.9476).  As on the muon side this is
NOT W-gluon fusion (O(alpha_s), charm-pair production).

Everything -- cards, integrators, fiducial region, EW inputs, binning -- is
imported unchanged from yadism_cc_calc.py; only the structure functions are
restricted to the charm-production channel.

Two INDEPENDENT calculations, compared bin by bin:
  (a) yadism F2_charm / FL_charm / F3_charm at PTO = 0, ZM-VFNS, prDIS: CC;
  (b) a direct LHAPDF parton-model integral over CKM-row-c weighted d, s, b.
NOTE their definitions can differ by the incoming-cbar channel
(W+ cbar -> dbar/sbar/bbar), which yadism's "charm" heavyness may include;
the ratio (a)/(b) printed below quantifies it.  The MC tag (outgoing hard
charm) corresponds to definition (b), which is what is written out.

Usage:  yadism_cc_charm_calc.py   -> results_nu/histos_yadism_charm.json
"""
import json
import hashlib
import os
import sys

import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import BINS, sel_suffix  # noqa: E402
import beams as beams_mod  # noqa: E402
import selection  # noqa: E402  -- the one Q2 floor
from yadism_calc import sf_grid_nodes  # noqa: E402
from yadism_cc_calc import (ALPHA, CKM_STR, EBEAM, GEV2_PB, GF, MW2, OUTDIR,  # noqa: E402
                            PDFSET, TWO_KP, XMAX, integrate_q2y,
                            integrate_xq2, integrate_ytheta, make_cards,
                            make_sigma_red)
import math  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp


def ckm2_charm_row():
    """|V_cd|^2, |V_cs|^2, |V_cb|^2 : the weights of d, s, b -> c."""
    m = np.array([float(v) for v in CKM_STR.split()]).reshape(3, 3) ** 2
    return {1: m[1, 0], 3: m[1, 1], 5: m[1, 2]}


# --------------------------------------------------------------- (a) yadism
CC_CACHE = f"{OUTDIR}/fonll_cache"


def _cc_grid(pto, pdf, x_nodes, q2_nodes, fns, nfff):
    """One yadism CC run for the charm flavour, as raw [x, Q2] grids."""
    import yadism

    theory, observables = make_cards(pto, x_nodes, q2_nodes)
    theory["FNS"] = fns
    theory["NfFF"] = nfff
    kins = observables["observables"]["F2_total"]
    names = ["F2_charm", "FL_charm", "F3_charm"]
    observables["observables"] = {n: kins for n in names}
    out = yadism.run_yadism(theory, observables)
    res = out.apply_pdf_alphas_alphaqed_xir_xif(
        pdf, lambda muR: pdf.alphasQ(muR), lambda _muR: ALPHA, 1.0, 1.0)
    return {n[:2]: np.array([p["result"] for p in res[n]]).reshape(
        len(q2_nodes), len(x_nodes)).T for n in names}


def cc_fonll_splines(pto, pdf, damp=False, ffns_only=False):
    """FONLL general-mass CC charm, same construction as the NC side.

        F_FONLL = F^(nl+1) + f_thr * F^(d),   F^(d) = F^(nl) - F^(nl,0)
                =    ZM    + f_thr * (FFNS  -  FFN0)

    Damping multiplies the MASS CORRECTION only (see yadism_charm_calc.py for
    why that grouping and not the one the yadism docs imply, and for the
    docs' (1 - Q2/m2) typo).  The three component runs are cached, so the
    damped variant is free once the undamped one has run.
    """
    from yadism_charm_calc import cache_tag, damping_factor, grid_components

    # one definition, in yadism_calc.sf_grid_nodes()
    x_nodes, q2_nodes = sf_grid_nodes()
    if ffns_only:
        # THE FFNS n_f = 3 PIECE, WITH THE CHARM PDF SET TO ZERO (2026-09-19):
        # POWHEG-V2mc's process.  qmass > 0 hard-sets st_nlight = 3 and keeps
        # only d/s -> c and g -> c dbar/sbar (nu-DIS-master/init_processes.f).
        # yadism's FONLL-FFNS "charm", given a PDF WITH charm (NNPDF4.0 fits
        # it), ALSO convolves the charm PDF -- measured: its LO F2 is 22-83%
        # above 2 xi s(xi), and exactly equal to it once charm is zeroed.
        # So the grid is convolved with a charm-free copy of the same set,
        # never read from the FONLL cache.  alpha_s is still the set's own
        # (V2mc runs n_f = 3; worth -0.3% on V2mc, measured by reweighting).
        class _NoCharm:
            def __init__(self, base):
                self.base = base
            def xfxQ2(self, pid, x, q2):
                return 0.0 if abs(pid) == 4 else self.base.xfxQ2(pid, x, q2)
            def alphasQ2(self, q2):
                return self.base.alphasQ2(q2)
            def alphasQ(self, q):
                return self.base.alphasQ(q)
            def set(self):
                return self.base.set()
        _zm, ffns, _f0 = grid_components("nu", pto, _NoCharm(pdf), "charm",
                                         x_nodes, q2_nodes)
        lx, lq = np.log(x_nodes), np.log(q2_nodes)
        return {k: RectBivariateSpline(lx, lq, ffns[k], kx=3, ky=3)
                for k in ("F2", "FL", "F3")}
    os.makedirs(CC_CACHE, exist_ok=True)
    # >>> THE GRID IS PART OF THE KEY -- see yadism_charm_calc.py. <<<
    # Keyed on pto alone, the 400 GeV and 4 TeV runs reused the 1 TeV cache
    # and combined 1 TeV structure functions with their own kinematics,
    # giving FONLL/ZM = 1.208 at 400 GeV and 0.721 at 4 TeV when mass effects
    # must VANISH at high Q2.
    grid_id = hashlib.sha1(
        np.concatenate([np.asarray(x_nodes, dtype=float),
                        np.asarray(q2_nodes, dtype=float)]).tobytes()
    ).hexdigest()[:12]
    # THE PDF AND THE TMC MODE ARE IN THE KEY (2026-09-13) -- these arrays
    # are PDF-convolved; see yadism_charm_calc._fonll_components.  "" at the
    # proton with TMC off, so the pre-cache files keep their names.
    fn = f"{CC_CACHE}/cc_charm_pto{pto}_{grid_id}{cache_tag(pdf)}.npz"
    if os.path.exists(fn):
        z = np.load(fn)
        print(f"      [cache] reusing {os.path.basename(fn)}", flush=True)
        zm, ffns, ffn0 = ({k: z[f"{t}_{k}"] for k in ("F2", "FL", "F3")}
                          for t in ("zm", "ffns", "ffn0"))
    else:
        # the massive runs through the shared PineAPPL grids, as on the NC
        # side (yadism_charm_calc.grid_components)
        zm, ffns, ffn0 = grid_components("nu", pto, pdf, "charm", x_nodes,
                                         q2_nodes)
        np.savez_compressed(
            fn, **{f"zm_{k}": v for k, v in zm.items()},
            **{f"ffns_{k}": v for k, v in ffns.items()},
            **{f"ffn0_{k}": v for k, v in ffn0.items()})

    lx, lq = np.log(x_nodes), np.log(q2_nodes)
    f_thr = damping_factor(q2_nodes) if damp else np.ones(len(q2_nodes))
    if damp:
        print(f"      damping f_thr: {f_thr[0]:.4f} at Q2={q2_nodes[0]:.2f}"
              f" ... {f_thr[-1]:.4f} at Q2={q2_nodes[-1]:.0f}", flush=True)
    return {k: RectBivariateSpline(
        lx, lq, zm[k] + f_thr[None, :] * (ffns[k] - ffn0[k]), kx=3, ky=3)
        for k in ("F2", "FL", "F3")}


def charm_structure_function_splines(pdf, pto=0):
    """yadism CC single-flavour (charm) F2/FL/F3 at PTO on the shared grid.

    CAREFUL -- this is NOT the benchmark's MC-tag definition.  yadism's
    charm-flavour observable also carries the INCOMING-charm channel
    (c -> s, c -> d), which the tag (outgoing charm, s/d -> c) excludes; at
    LO the two differ by 7.24% (0.69092 vs 0.64429 pb).  From O(alpha_s) it
    additionally carries the charm share of the gluon terms.
    """
    import yadism

    # one definition, in yadism_calc.sf_grid_nodes()
    x_nodes, q2_nodes = sf_grid_nodes()
    theory, observables = make_cards(pto, x_nodes, q2_nodes)
    kins = observables["observables"]["F2_total"]
    names = ["F2_charm", "FL_charm", "F3_charm"]
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
def direct_charm_dsig(pdf):
    """LO parton-model d2sigma for W+ + d/s/b -> c, straight out of LHAPDF.

    Quark (not antiquark) channels only -> +xF3, so the y dependence is flat:
    Y+ F2 + Y- xF3 = 2 * F2.  Same differential interfaces as
    yadism_cc_calc.make_sigma_red so the shared integrators apply unchanged.
    """
    chan = ckm2_charm_row()

    def qsum(x, q2):
        x = np.asarray(x, float)
        q2 = np.asarray(q2, float)
        flat = np.zeros(x.size)
        xr, qr = x.ravel(), q2.ravel()
        for pid, cw in chan.items():
            flat += cw * np.array(
                [pdf.xfxQ2(pid, float(xx), float(qq)) if xx < 1.0 else 0.0
                 for xx, qq in zip(xr, qr)])
        return flat.reshape(x.shape)

    def dsig_dxdq2(x, q2, y):
        x, q2, y = np.broadcast_arrays(np.asarray(x, float),
                                       np.asarray(q2, float),
                                       np.asarray(y, float))
        f2 = 2.0 * qsum(x, q2)   # F2 = 2x sum w q  (the x is inside xfxQ2)
        # comb = Y+ F2 + Y- xF3 with FL = 0, xF3 = +F2 (pure quark channels)
        #      = (Y+ + Y-) F2 = 2 F2  -> flat in y
        prop = (MW2 / (MW2 + q2)) ** 2
        out = (GF * GF / (4.0 * math.pi * x) * prop
               * 2.0 * f2 * GEV2_PB)
        return np.where(x > XMAX, 0.0, out)

    def dsig_dq2dy(q2, y):
        q2, y = np.broadcast_arrays(np.asarray(q2, float), np.asarray(y, float))
        x = q2 / (y * TWO_KP)
        return dsig_dxdq2(x, q2, y) * x / y

    return dsig_dxdq2, dsig_dq2dy


def run(dsig_dxdq2, dsig_dq2dy):
    hists = {}
    s1, means = integrate_q2y(dsig_dq2dy, hists)
    s2 = integrate_xq2(dsig_dxdq2, hists)
    s3 = integrate_ytheta(dsig_dq2dy, hists)
    return s1, s2, s3, means, hists


ORDERS = {"lo": 0, "nlo": 1, "nnlo": 2}


def main():
    order = (sys.argv[1] if len(sys.argv) > 1 else "lo").lower()
    if order not in ORDERS:
        sys.exit(f"unknown order '{order}' (use {'|'.join(ORDERS)})")
    pto = ORDERS[order]
    scheme = (sys.argv[2] if len(sys.argv) > 2 else "zm").lower()
    if scheme not in ("zm", "fonll", "fonll_damp", "ffns3"):
        sys.exit("scheme must be zm, fonll, fonll_damp or ffns3")
    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(PDFSET, 0)

    print("[cc_charm] (b) direct LHAPDF d/s/b -> c integral ...", flush=True)
    dxq_b, dqy_b = direct_charm_dsig(pdf)
    sb1, sb2, sb3, means_b, hists_b = run(dxq_b, dqy_b)
    print(f"[cc_charm] (b) direct   sigma_fid = {sb1:.5f} pb "
          f"(cross-checks: x-Q2 {sb2:.5f}, y-theta {sb3:.5f})")

    print(f"[cc_charm] (a) yadism {order.upper()} (PTO={pto}), "
          f"F2/FL/F3_charm, prDIS CC ...", flush=True)
    try:
        splines = (cc_fonll_splines(pto, pdf,
                                    damp=(scheme == 'fonll_damp'),
                                    ffns_only=(scheme == 'ffns3'))
                   if scheme.startswith('fonll') or scheme == 'ffns3' else
                   charm_structure_function_splines(pdf, pto))
        dxq_a, dqy_a = make_sigma_red(splines)
        sa1, sa2, sa3, means_a, hists_a = run(dxq_a, dqy_a)
        print(f"[cc_charm] (a) yadism   sigma_fid = {sa1:.5f} pb "
              f"(cross-checks: x-Q2 {sa2:.5f}, y-theta {sa3:.5f})")
        print(f"[cc_charm] ratio (a)/(b) = {sa1/sb1:.5f}  (any excess over 1 "
              "is the incoming-cbar channel in yadism's charm definition)")
    except Exception as e:                                    # noqa: BLE001
        print(f"[cc_charm] yadism CC charm failed ({e}); "
              "using the direct integral alone")
        sa1 = None

    # definition (b) -- produced charm from d/s/b, matching the MC tag -- is
    # what the benchmark compares against
    # DEFINITION (user, 2026-08-20): on the YADISM side charm production is
    # ALL CONTRIBUTIONS PROPORTIONAL TO V_cd OR V_cs -- which is exactly
    # yadism's charm-flavour observable -- while on the generator side it is
    # all events with at least one charm quark in the FINAL state.
    #
    # The two are not the same object.  The channel nu cbar -> mu- sbar
    # carries |V_cs|^2, so it belongs to the YADISM side, but its charm is in
    # the INITIAL state so the generator tag drops it.  Measured directly:
    # that channel alone is 0.04679 pb against the yadism-minus-direct gap of
    # 0.69092 - 0.64429 = 0.04663 pb, agreeing to 0.3%.  So the definitions
    # differ by +7.3% AT LO, not only at NNLO.  It is also an ANTIQUARK
    # channel (xF3 = -F2, a (1-y)^2 shape against the flat quark channels), so
    # it tilts the y distribution as well as the rate.
    #
    # The user's definition is adopted because it is the one that EXISTS above
    # LO: the direct d/s/b -> c integral is a LO parton-model construction
    # with no NLO counterpart, whereas the observable runs at any PTO.  The
    # direct integral is still computed at LO and reported as a cross-check,
    # and the offset is stated on the report page.
    if sa1 is None:
        sys.exit("[cc_charm] yadism CC charm failed -- it is the published "
                 "definition, so there is no fallback")
    s1, means, hists = sa1, means_a, hists_a
    err = max(abs(sa1 - sa2), abs(sa1 - sa3), 2e-4 * sa1)
    if pto == 0:
        print(f"[cc_charm] LO cross-check: direct d/s/b -> c integral "
              f"{sb1:.5f} pb = final-state-charm definition; published "
              f"(V_cd/V_cs terms) {sa1:.5f} pb; offset "
              f"{100*(sa1/sb1 - 1):+.2f}% = the incoming-cbar channel")
    print(f"[cc_charm] means: <Q2> = {means['Q2']:.3f} GeV2, "
          f"<x> = {means['xbj']:.5f}, <y> = {means['y']:.4f}")

    suffix = "" if order == "lo" else f"_{order}"
    if scheme == "fonll":
        suffix += "_fonll"
    elif scheme == "fonll_damp":
        suffix += "_fonll_damp"
    elif scheme == "ffns3":
        suffix += "_ffns3"
    defn = ("all contributions proportional to V_cd or V_cs (yadism "
            "charm-flavour observable); the generator tag is final-state "
            "charm, which excludes the incoming-cbar channel, +7.3% at LO")
    out = {
        "generator": f"yadism_charm{suffix}",
        "label": (f"YADISM CC {order.upper()} (charm"
                  + {"fonll": ", FONLL", "fonll_damp": ", FONLL",
                     "ffns3": ", FFNS n_f = 3"}
                  .get(scheme, "") + ")"),
        "scheme": {"fonll": "FONLL general-mass (undamped)",
                   "fonll_damp": "FONLL general-mass, damped",
                   "ffns3": "FFNS n_f = 3: massive charm produced from "
                            "s/d/g, charm PDF set to zero (POWHEG-V2mc's "
                            "process)"
                   }.get(scheme, "ZM-VFNS"),
        "pto": pto,
        "definition": defn,
        "n_parsed": None,
        "n_fiducial": None,
        "sigma_fid_pb": s1,
        "sigma_fid_err_pb": err,
        "sigma_fid_pb_yadism_charm_obs": sa1,
        "sigma_fid_pb_direct_finalstate_charm": sb1,
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
    # SELECTION, TARGET AND TMC IN THE NAME, as on the NC side -- see
    # yadism_charm_calc.main.  All "" at the defaults.
    import target
    ofn = (f"{OUTDIR}/"
           f"{beams_mod.at_energy(f'histos_yadism_charm{suffix}' + sel_suffix() + target.suffix() + target.tmc_suffix(), EBEAM)}"
           f"{selection.q2_suffix()}.json")
    with open(ofn, "w") as f:
        stamp(out, current="nu")
        prov = target.provenance(PDFSET)
        if prov:              # absent at proton / TMC off: earlier stays byte-identical
            out.update(prov)
        json.dump(out, f)
    print(f"[cc_charm] wrote {ofn}")


if __name__ == "__main__":
    main()
