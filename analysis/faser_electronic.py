#!/usr/bin/env python3
"""Our predictions against FASER's ELECTRONIC-detector muon-neutrino results.

The sub-tab the user asked for on 2026-09-08.  Two measurements, both made by
reconstructing the muon in the FASER spectrometer while the FASERnu tungsten
is only the target: CERN-FASER-CONF-2026-005 (186 fb^-1) and arXiv:2412.03186
(65.6 fb^-1).  See data/faser_electronic/README.md for what is transcribed and
what is digitised.

>>> THE OBSERVABLE IS -L/E_nu, AND THAT IS WHY NO RAPIDITY INFORMATION IS
NEEDED. <<<  L is the LEPTON NUMBER, so the SIGN of the axis is the neutrino
charge and its magnitude is 1/E: one variable carrying energy and charge
together, which our one-dimensional per-charge flux files predict directly.
The note's other figures ARE differential in rapidity and are not reproduced;
data/faser_flux has no rapidity information and there is no way to fake it.

THREE COMPARISONS, IN INCREASING ORDER OF HOW MUCH DETECTOR IS IN THE WAY.

  1. THE FLUX-WEIGHTED CROSS-SECTION PER NUCLEON, bin by bin, against
     Table III of the PRL.  No detector, no geometry, no luminosity -- the
     cleanest comparison anywhere in this benchmark, and it is possible only
     because FASER published that column as numbers.
  2. THE ACCEPTANCE of the muon cuts, against Table II of the PRL.  FASER's
     spectrometer selection is a muon above 100 GeV within 25 mrad of the
     axis, which is EXACTLY this benchmark's Tier S -- so the fraction of our
     events passing it is directly comparable with their alpha.  This is the
     one place the benchmark's own selection can be checked against an
     experiment's.
  3. THE EVENT YIELDS, against Fig. 5 of the note and Table I of the PRL.
     This one needs the geometry, and carries the approximation named below.

>>> THE ONE APPROXIMATION, AND IT IS IN THE YIELDS ALONE. <<<  The measured
fiducial volume is a cylinder of 100 mm radius; the flux files are a COUNT
through their authors' 25 x 25 cm aperture.  Converting between them assumes
the flux is UNIFORM across the aperture, which it is not -- it peaks on axis,
so this UNDERESTIMATES the count through a central cylinder.  The factor is
geometric, A_fid / A_flux, and is reported so a reader can undo it.  The size
of the miss is measured rather than guessed: the PRL publishes its own
simulated flux per bin (Table III), so our folded flux is compared with it
directly and the ratio IS the non-uniformity.

Usage: analysis/faser_electronic.py -> results_nu/faser_electronic.json

It reads results_nu/faser_powheg_rates.json, the ubexcess-corrected ladder
with genuine neutrons (analysis/faser_powheg_rates.py); it is what paper plot
12b draws.
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import faser_rates as fr                                     # noqa: E402
import faser_dimuon_cutflow as fdc                           # noqa: E402
from faser_genie_rates import NUCLEI_PER_CM2, Z_W, N_W       # noqa: E402

DATA = f"{BASE}/data/faser_electronic/faser_electronic.json"
LADDER = f"{BASE}/results_nu/faser_powheg_rates.json"
OUT = f"{BASE}/results_nu/faser_electronic.json"

# The measured fiducial volume: a cylinder of 100 mm radius about the
# spectrometer axis (both papers).  >>> THE FLUX IS COUNTED THROUGH THAT
# CYLINDER ITSELF (2026-10-07), not taken as the face count times an area
# ratio: the flux peaks on the line of sight, and inside r < 100 mm it is 2%
# (below 100 GeV) to 25% (above 1.5 TeV) above the face average -- about 6%
# of the rate, invisibly. <<<  data/faser_flux_2025/*_fid_r100.txt, written by
# tools/make_flux_fid_spectra.py.  The column density stays that of the face
# (NUCLEI_PER_CM2): the tungsten is equally thick everywhere.
R_FID_CM = 10.0
A_FID_CM2 = math.pi * R_FID_CM ** 2
FID_APERTURE = "fid_r100"
GEOM = 1.0                     # the count already is through the cylinder

# The seven scale points POWHEG-V2's reweighting carries, in the ladder's own
# order; index 0 is the nominal.
N_SCALE = 7

# >>> SHERPA NLO, ADDED FOR FIG. 6.2 (user, 2026-09-30). <<<  Its FASER ladder
# (tools/sherpa_faser_ladder.sh, Q2 > 4 and no y cut -- the region POWHEG-V2's
# "ally" fraction selects) was built for the emulsion, whose 200 GeV lepton
# needs nothing below 200 GeV and nothing much above 4 TeV.  FASER's
# electronic bins start at 100 GeV and the > 1 TeV bin runs with the flux to
# 6.8 TeV, so the 100 and 6800 GeV points were added for this figure: without
# them the fold ZEROES the cross-section below the first point and freezes it
# above the last, which would gut the two 100-300 GeV bins silently.
# Protons are 2212 beams; neutrons the genuine-2112 re-runs.
SHERPA_ES = [100, 200, 400, 700, 1000, 2000, 4000, 6800]


def sherpa_ladder(cur, tgt):
    """(E, sigma_pb) for one Sherpa arm, the delivered cross-section with the
    neutrino spin factor (faser_emulsion_shapes.sherpa_total_pb, the one
    source for reading it).  Every point must exist: a missing end point is
    not an interpolation, it is a zero."""
    import paths
    from faser_emulsion_shapes import sherpa_total_pb
    pre = "V2_" if tgt == "n" else ""
    sig, miss = [], []
    for e in SHERPA_ES:
        d = f"{paths.SHERPA_RUNS}/{pre}NuDIS_NLO_FASER_{cur}_{tgt}_E{e}"
        v = sherpa_total_pb(d)
        if not v:
            miss.append(os.path.basename(d))
        sig.append(v)
    if miss:
        raise SystemExit("Sherpa ladder incomplete, refusing to fold: "
                         + " ".join(miss))
    return np.array(SHERPA_ES, float), np.array(sig, float)


HERWIG_ES = [30, 60, 100, 200, 300, 400, 700, 1000, 2000, 4000, 6800]


def herwig_ladder(cur, tgt):
    """(E, sigma_pb) for one Herwig FASER-ladder arm (tools/herwig_faser_ladder.sh,
    2026-10-04): positive minus negative half, every point required."""
    from faser_emulsion_shapes import herwig_ladder_sigma
    sig, miss = [], []
    for e in HERWIG_ES:
        v = herwig_ladder_sigma(cur, tgt, e)
        if v is None:
            miss.append(f"{cur}_{tgt}_E{e}")
        sig.append(v[0] if v else None)
    if miss:
        raise SystemExit("Herwig ladder incomplete, refusing to fold: " + " ".join(miss))
    return np.array(HERWIG_ES, float), np.array(sig, float)


def ladder(cur, tgt, key="ally"):
    """(E, sigma_pb per nucleon) and the seven scale variants, for one arm."""
    with open(LADDER) as f:
        d = json.load(f)
    pts = d["ladder"][f"{cur}_{tgt}"]
    e = np.array([p["E_GeV"] for p in pts], float)
    tot = np.array([p["sigma_pb"] for p in pts], float)
    fid = np.array([p["fid_fraction"][key] for p in pts], float)
    sc = None
    if all(p.get("sigma_fid_by_scale_pb") for p in pts):
        sc = np.array([p["sigma_fid_by_scale_pb"][key] for p in pts], float)
    return e, tot, fid, sc


# >>> THE Q2 < 4 GeV2 PIECE, FROM GENIE (user, 2026-10-04). <<<  POWHEG-V2 and
# Sherpa are generated above Q2 = 4 GeV2, while FASER's cross-section and
# yields are the total charged-current ones.  To compare like with like the
# region below the floor is added to both from GENIE: at each ladder energy
# the weighted fraction of GENIE events (all CC channels, default tune) with
# Q2 < 4 GeV2, from the per-event tables of the emulsion ladder
# (genie/fdladder_job_<cur>_<tgt>_E<E>_1/faserdata_events.npz), times
# GENIE's total from its splines.  The q2 > 4 rows are kept beside the
# extended ones in the output (keys without "_lowq2").
Q2_FLOOR = 4.0


def genie_lowq2_sigma(cur, tgt, es):
    """GENIE's sigma(Q2 < 4 GeV2) per nucleon [pb] at energies `es`."""
    from faser_emulsion_shapes import LADDER_ES
    pid = 14 if cur == "nu" else -14
    ee, fr_low = [], []
    for e in LADDER_ES:
        f = f"{BASE}/genie/fdladder_job_{cur}_{tgt}_E{e}_1/faserdata_events.npz"
        if not os.path.exists(f):
            raise SystemExit(f"missing {os.path.relpath(f, BASE)}: cannot "
                             "estimate the Q2 < 4 GeV2 piece")
        z = np.load(f, allow_pickle=False)
        w, q2 = z["w"].astype(float), z["q2"].astype(float)
        ee.append(float(e))
        fr_low.append(float(w[q2 < Q2_FLOOR].sum() / w.sum()))
    ee, fr_low = np.array(ee), np.array(fr_low)
    frac = np.interp(np.log(es), np.log(ee), fr_low)
    tot = fdc.genie_sigma_tot(pid, tgt, np.asarray(es, float))
    return frac * tot, dict(zip([float(x) for x in ee], fr_low.tolist()))


def _interp_sigma(e_grid, es, sig):
    """log-log in energy, ZERO below the ladder and frozen above it -- the
    same convention faser_powheg_rates.rate and faser_dimuon_cutflow.convolve
    use, so the three cannot disagree about the tails."""
    with np.errstate(divide="ignore"):
        v = np.exp(np.interp(np.log(e_grid), np.log(es), np.log(sig),
                             left=-np.inf, right=np.log(sig[-1])))
    return np.where(e_grid < es[0], 0.0, v)


def fold(pid, sig_p, sig_n, es, e_lo, e_hi, lumi_fb):
    """(flux, <sigma> per nucleon, events) restricted to one energy window.

    The flux files are a COUNT at 150 fb^-1 through the authors' aperture; the
    events are that count times sigma times the column density, and then times
    the geometric factor for the fiducial cylinder.
    """
    e, phi = fr.flux(str(pid), aperture=FID_APERTURE)
    phi = phi * (lumi_fb / fr.LUMI_FLUX_NU)
    m = (e >= e_lo) & ((e < e_hi) if e_hi else np.ones_like(e, bool))
    if not m.any():
        return 0.0, float("nan"), 0.0
    sp = _interp_sigma(e, es, sig_p)
    sn = _interp_sigma(e, es, sig_n)
    per_nucleon = (Z_W * sp + N_W * sn) / (Z_W + N_W)
    phi_b = float(phi[m].sum())
    sig_mean = (float((phi[m] * per_nucleon[m]).sum()) / phi_b
                if phi_b > 0 else float("nan"))
    n = float(((Z_W * sp[m] + N_W * sn[m]) * phi[m]).sum()
              * 1e-36 * NUCLEI_PER_CM2) * GEOM
    return phi_b, sig_mean, n


def _scale_arms(charge, j, sc_p, sc_n, sc_pb, sc_nb):
    """(pid, sigma_p, sigma_n) at the j-th scale point, per charge."""
    out = []
    if charge in ("nu", "both"):
        out.append((14, sc_p[:, j], sc_n[:, j]))
    if charge in ("nubar", "both"):
        out.append((-14, sc_pb[:, j], sc_nb[:, j]))
    return out


def main():
    with open(DATA) as f:
        dat = json.load(f)
    bins = dat["bins"]
    prl, conf = dat["prl"], dat["conf"]

    es, tot_p, fid_p, sc_p = ladder("nu", "p", "ally")   # sc_* : 7 scales
    _e, tot_n, fid_n, sc_n = ladder("nu", "n", "ally")
    _e, _t, fs_p, _s = ladder("nu", "p", "faser_s")
    _e, _t, fs_n, _s = ladder("nu", "n", "faser_s")
    _e, _t, ka_p, _s = ladder("nu", "p", "faser_s_ally")
    _e, _t, ka_n, _s = ladder("nu", "n", "faser_s_ally")
    es_b, tot_pb_, fid_pb_, sc_pb_ = ladder("nubar", "p", "ally")
    _e, tot_nb, fid_nb, sc_nb = ladder("nubar", "n", "ally")
    _e, _t, fs_pb, _s = ladder("nubar", "p", "faser_s")
    _e, _t, fs_nb, _s = ladder("nubar", "n", "faser_s")
    _e, _t, ka_pb, _s = ladder("nubar", "p", "faser_s_ally")
    _e, _t, ka_nb, _s = ladder("nubar", "n", "faser_s_ally")

    def arms(charge):
        if charge == "nu":
            return [(14, tot_p * fid_p, tot_n * fid_n,
                     tot_p * fs_p, tot_n * fs_n,
                     tot_p * ka_p, tot_n * ka_n)]
        if charge == "nubar":
            return [(-14, tot_pb_ * fid_pb_, tot_nb * fid_nb,
                     tot_pb_ * fs_pb, tot_nb * fs_nb,
                     tot_pb_ * ka_pb, tot_nb * ka_nb)]
        return arms("nu") + arms("nubar")

    # >>> GENIE'S TOTAL, WHICH IS WHAT FASER'S OWN SIMULATION USES. <<<  Its
    # splines carry the whole charged-current cross-section with no Q2 floor,
    # so it is the one row here that is like for like with FASER's published
    # number -- and how closely it lands on it is a check of our machinery
    # rather than of the physics.  POWHEG-V2's ladder is generated above
    # Q2 > 4 GeV2 and is low by that part, by construction.
    def genie_arms(charge):
        out = []
        if charge in ("nu", "both"):
            out.append(14)
        if charge in ("nubar", "both"):
            out.append(-14)
        return out

    sh = {(c, t): sherpa_ladder(c, t) for c in ("nu", "nubar")
          for t in ("p", "n")}
    low = {(c, t): genie_lowq2_sigma(c, t, es) for c in ("nu", "nubar")
           for t in ("p", "n")}
    low_sh = {(c, t): genie_lowq2_sigma(c, t, sh[(c, t)][0])[0]
              for c in ("nu", "nubar") for t in ("p", "n")}
    hw = {(c, t): herwig_ladder(c, t) for c in ("nu", "nubar") for t in ("p", "n")}
    low_hw = {(c, t): genie_lowq2_sigma(c, t, hw[(c, t)][0])[0]
              for c in ("nu", "nubar") for t in ("p", "n")}

    def ext_arms(charge):
        """POWHEG-V2 Q2 > 4 plus GENIE's Q2 < 4 piece, per charge."""
        out = []
        if charge in ("nu", "both"):
            out.append((14, tot_p * fid_p + low[("nu", "p")][0],
                        tot_n * fid_n + low[("nu", "n")][0]))
        if charge in ("nubar", "both"):
            out.append((-14, tot_pb_ * fid_pb_ + low[("nubar", "p")][0],
                        tot_nb * fid_nb + low[("nubar", "n")][0]))
        return out

    rows = []
    for i, b in enumerate(bins):
        r = {"bin": i + 1, "charge": b["charge"],
             "e_lo": b["e_lo"], "e_hi": b["e_hi"]}
        for lumi, tag in ((conf["lumi_fb"], "conf"), (prl["lumi_fb"], "prl")):
            phi = sig_w = 0.0
            n_all = n_s = n_k = 0.0
            for pid, sp, sn, sps, sns, spk, snk in arms(b["charge"]):
                p, sm, n = fold(pid, sp, sn, es, b["e_lo"], b["e_hi"], lumi)
                _p2, _sm2, ns = fold(pid, sps, sns, es,
                                     b["e_lo"], b["e_hi"], lumi)
                _p3, _sm3, nk = fold(pid, spk, snk, es,
                                     b["e_lo"], b["e_hi"], lumi)
                phi += p
                sig_w += p * sm if p > 0 and np.isfinite(sm) else 0.0
                n_all += n
                n_s += ns
                n_k += nk
            r[f"n_fid_{tag}"] = n_all
            r[f"n_tiers_{tag}"] = n_s
            r[f"n_muoncuts_{tag}"] = n_k
            if tag == "prl":
                r["flux_1e6_fb_cm2"] = (phi / (A_FID_CM2
                                               * prl["lumi_fb"]) / 1e6)
                r["sigma_1e38_cm2_per_nucleon"] = (
                    (sig_w / phi) * 1e2 if phi > 0 else float("nan"))
        # >>> TWO ACCEPTANCES, AND THE DIFFERENCE BETWEEN THEM IS THE POINT.
        # <<<  `acceptance_pct` is the KINEMATIC one -- FASER's muon cuts and
        # nothing else -- which is what their Table II alpha is a fraction of.
        # `tier_s_pct` is this benchmark's own Tier S, which carries the
        # 0.2 < y < 0.9 window too; at 100-300 GeV that window is the dominant
        # loss (E' > 100 GeV at E = 200 already means y < 0.5, so y > 0.2
        # removes the stiffest muons) and comparing it with their alpha would
        # be comparing two different things.
        # GENIE, from its own splines: total CC per nucleon, no Q2 floor
        gphi = gsig = gn = 0.0
        for pid in genie_arms(b["charge"]):
            sp = fdc.genie_sigma_tot(pid, "p", es)
            sn = fdc.genie_sigma_tot(pid, "n", es)
            p, sm, n = fold(pid, sp, sn, es, b["e_lo"], b["e_hi"],
                            conf["lumi_fb"])
            gphi += p
            gsig += p * sm if p > 0 and np.isfinite(sm) else 0.0
            gn += n
        r["genie_sigma_1e38_cm2_per_nucleon"] = (
            (gsig / gphi) * 1e2 if gphi > 0 else float("nan"))
        r["genie_n_fid_conf"] = gn
        # Sherpa NLO, folded exactly as GENIE is
        sphi = ssig = sn_ = 0.0
        for pid in genie_arms(b["charge"]):
            c = "nu" if pid > 0 else "nubar"
            e_sh, sp = sh[(c, "p")]
            _e, snn = sh[(c, "n")]
            p, sm, n = fold(pid, sp, snn, e_sh, b["e_lo"], b["e_hi"],
                            conf["lumi_fb"])
            sphi += p
            ssig += p * sm if p > 0 and np.isfinite(sm) else 0.0
            sn_ += n
        r["sherpa_sigma_1e38_cm2_per_nucleon"] = (
            (ssig / sphi) * 1e2 if sphi > 0 else float("nan"))
        r["sherpa_n_fid_conf"] = sn_
        # ... and both extended below Q2 = 4 GeV2 with GENIE's piece
        xphi = xsig = xn = xs_sig = xs_n = 0.0
        for pid, sp, sn in ext_arms(b["charge"]):
            p, sm, n = fold(pid, sp, sn, es, b["e_lo"], b["e_hi"],
                            conf["lumi_fb"])
            xphi += p
            xsig += p * sm if p > 0 and np.isfinite(sm) else 0.0
            xn += n
        for pid in genie_arms(b["charge"]):
            c = "nu" if pid > 0 else "nubar"
            e_sh, sp = sh[(c, "p")]
            _e, snn = sh[(c, "n")]
            p, sm, n = fold(pid, sp + low_sh[(c, "p")], snn + low_sh[(c, "n")],
                            e_sh, b["e_lo"], b["e_hi"], conf["lumi_fb"])
            xs_sig += p * sm if p > 0 and np.isfinite(sm) else 0.0
            xs_n += n
        # Herwig NLO, folded as Sherpa, raw and extended below Q2 = 4
        hphi = hsig = hn_ = hxsig = hxn = 0.0
        for pid in genie_arms(b["charge"]):
            c = "nu" if pid > 0 else "nubar"
            e_hw, sp = hw[(c, "p")]
            _e, snn = hw[(c, "n")]
            p, sm, n = fold(pid, sp, snn, e_hw, b["e_lo"], b["e_hi"], conf["lumi_fb"])
            hphi += p
            hsig += p * sm if p > 0 and np.isfinite(sm) else 0.0
            hn_ += n
            p2, sm2, n2 = fold(pid, sp + low_hw[(c, "p")], snn + low_hw[(c, "n")],
                               e_hw, b["e_lo"], b["e_hi"], conf["lumi_fb"])
            hxsig += p2 * sm2 if p2 > 0 and np.isfinite(sm2) else 0.0
            hxn += n2
        r["herwig_sigma_1e38_cm2_per_nucleon"] = (hsig / hphi) * 1e2 if hphi > 0 else float("nan")
        r["herwig_n_fid_conf"] = hn_
        r["herwig_sigma_1e38_cm2_per_nucleon_lowq2"] = (hxsig / hphi) * 1e2 if hphi > 0 else float("nan")
        r["herwig_n_fid_conf_lowq2"] = hxn
        r["sigma_1e38_cm2_per_nucleon_lowq2"] = (
            (xsig / xphi) * 1e2 if xphi > 0 else float("nan"))
        r["n_fid_conf_lowq2"] = xn
        r["sherpa_sigma_1e38_cm2_per_nucleon_lowq2"] = (
            (xs_sig / sphi) * 1e2 if sphi > 0 else float("nan"))
        r["sherpa_n_fid_conf_lowq2"] = xs_n
        r["acceptance_pct"] = (100.0 * r["n_muoncuts_conf"] / r["n_fid_conf"]
                               if r["n_fid_conf"] > 0 else float("nan"))
        r["tier_s_pct"] = (100.0 * r["n_tiers_conf"] / r["n_fid_conf"]
                           if r["n_fid_conf"] > 0 else float("nan"))
        # >>> THE 7-POINT SCALE BAND ON POWHEG-V2, which the user asked for
        # explicitly. <<<  POWHEG's own reweighting of every ladder point
        # gives the fiducial cross-section at each of the seven scales
        # (tools/powheg_v2_faser_ladder_rwgt.sh), so the whole fold is redone
        # at each and the envelope taken -- rather than a band on the
        # cross-section being propagated by hand, which would lose the
        # energy dependence the fold is sensitive to.
        if sc_p is not None:
            band_n, band_s, band_x = [], [], []
            for j in range(N_SCALE):
                nn = ss = 0.0
                nx = 0.0
                for pid, sp, sn in _scale_arms(b["charge"], j, sc_p, sc_n,
                                               sc_pb_, sc_nb):
                    p2, sm2, n2 = fold(pid, sp, sn, es,
                                       b["e_lo"], b["e_hi"], conf["lumi_fb"])
                    nn += n2
                    ss += p2 * sm2 if p2 > 0 and np.isfinite(sm2) else 0.0
                    c = "nu" if pid > 0 else "nubar"
                    _p3, _s3, n3 = fold(pid, sp + low[(c, "p")][0],
                                        sn + low[(c, "n")][0], es,
                                        b["e_lo"], b["e_hi"], conf["lumi_fb"])
                    nx += n3
                band_n.append(nn)
                band_s.append(ss)
                band_x.append(nx)
            r["n_fid_conf_by_scale"] = band_n
            r["n_fid_conf_scale_hi"] = max(band_n)
            r["n_fid_conf_scale_lo"] = min(band_n)
            r["n_fid_conf_lowq2_scale_hi"] = max(band_x)
            r["n_fid_conf_lowq2_scale_lo"] = min(band_x)
        rows.append(r)

    # ------------------------------------------------------------- report
    out = {"what": ("this benchmark against FASER's electronic-detector "
                    "muon-neutrino measurements, in bins of -L/E_nu"),
           "geometry": {
               "fiducial_radius_cm": R_FID_CM,
               "fiducial_area_cm2": A_FID_CM2,
               "flux_aperture_cm2": A_FID_CM2,
               "factor": GEOM,
               "note": ("the flux files are a COUNT through their authors' "
                        "aperture, so the fiducial cylinder takes this "
                        "geometric share of it.  That assumes the flux is "
                        "UNIFORM across the aperture and it is not -- it "
                        "peaks on axis -- so the yields are an "
                        "UNDERESTIMATE.  The size of the miss is measured "
                        "against the PRL's own simulated flux, below")},
           "cross_section_note": ("our POWHEG-V2 ladder is generated above "
                                  "Q2 > 4 GeV2, where FASER's quoted "
                                  "cross-section is the total; 4-8% of the "
                                  "neutrino CC rate lies below that floor "
                                  "over this energy range (see the non-DIS "
                                  "appendix), so our sigma is low by about "
                                  "that much and by construction"),
           "selection_note": ("Tier S is E' > 100 GeV and theta < 25 mrad, "
                              "which is FASER's spectrometer selection; the "
                              "acceptance column compares our fraction with "
                              "their Table II"),
           "lumi_conf_fb": conf["lumi_fb"], "lumi_prl_fb": prl["lumi_fb"],
           "ladder": os.path.basename(LADDER),
           "lowq2_note": ("keys ending _lowq2 add GENIE's Q2 < 4 GeV2 piece "
                          "(weighted event fraction x spline total) to "
                          "POWHEG-V2 and Sherpa, so they are total CC like "
                          "FASER's numbers"),
           "genie_lowq2_fraction": {f"{c}_{t}": low[(c, t)][1]
                                    for c in ("nu", "nubar")
                                    for t in ("p", "n")},
           "rows": rows}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)

    print("bin   E range      POWHEG <s>   GENIE <s>  FASER sim   P/F   G/F"
          "    our flux/FASER   accept  FASER")
    for i, (r, b) in enumerate(zip(rows, bins)):
        s_f = prl["sigma_1e38_cm2_per_nucleon"][i]
        f_f = prl["flux_1e6_fb_cm2"][i]
        rng = (f"{b['e_lo']:.0f}-{b['e_hi']:.0f}" if b["e_hi"]
               else f">{b['e_lo']:.0f}")
        print(f"{i+1} {b['charge']:5s} {rng:>10s} "
              f"{r['sigma_1e38_cm2_per_nucleon']:9.1f} "
              f"{r['genie_sigma_1e38_cm2_per_nucleon']:9.1f} {s_f:9.1f} "
              f"{r['sigma_1e38_cm2_per_nucleon']/s_f:6.3f} "
              f"{r['genie_sigma_1e38_cm2_per_nucleon']/s_f:6.3f} "
              f"      {r['flux_1e6_fb_cm2']/f_f:6.3f}   "
              f"{r['acceptance_pct']:6.1f} {prl['acceptance'][i]:6.1f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
