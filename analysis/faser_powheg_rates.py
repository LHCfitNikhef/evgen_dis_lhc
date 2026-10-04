#!/usr/bin/env python3
"""Step B of the neutrino-rate reproduction: POWHEG-V2 fiducial rates at
FASERnu, both neutrino flavours and both charges, on tungsten.

WHAT IT IS.  sigma_fid(E) -- the benchmark's fiducial region, Q2 > 4 GeV2
and 0.2 < y < 0.9 (user, 2026-09-03: "we still need Q2 > 4 GeV2 else the MC
output is not trustable") -- for nu_mu and nubar_mu on a proton and on a
neutron, from the POWHEG-V2 ladder of tools/powheg_v2_faser_ladder.sh:
at each energy the run's integrated cross-section (pwg-stat.dat) times the
fiducial fraction of its Les Houches events.  The fraction is exact for
these cuts: Q2 and y are built from the lepton momenta, which the shower
does not touch with QED radiation off, so no shower is needed.  A closure
against the showered productions at 400, 1000 and 4000 GeV is printed.

Tungsten is 74 p + 110 n per nucleus; nu_e uses the nu_mu cross-sections
(flavour universality; the lepton mass is irrelevant here).  The flux is
the vendored per-flavour, per-charge files scaled to 250 fb^-1 and the
target is the 1.1 t of tungsten behind the 25 x 30 cm face, as in
analysis/faser_genie_rates.py, so the two steps differ ONLY in the
cross-section: GENIE's total CC against POWHEG-V2's NLO fiducial region.
Their ratio is the fiducial acceptance under the FASER flux.

Usage: analysis/faser_powheg_rates.py [--rwgt-ladder]
Writes results_nu/faser_powheg_rates.json and prints the tables.

The fold is on the ladder of tools/faser_emulsion_ladder_ubexcess.sh -- ten
energies, ubexcess_correct 1, genuine 2112 neutron beams and 100000 events a
point.  Those Les Houches files carry no <rwgt> block, so the seven-point
band is IMPORTED as a RELATIVE shift per (arm, energy, region) from the
earlier, REWEIGHTED ladder ($POWHEG_V2/ladder-faser, 2026-09-08):
sigma_fid x sigma_fid_rwgt(scale j) / sigma_fid_rwgt -- the band is a
property of the scale choice, not of the unweighting bound.

--rwgt-ladder folds that reweighted ladder itself and writes
results_nu/faser_powheg_rates_rwgt.json, the file the band is read from.
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import beams  # noqa: E402
import faser_rates as fr  # noqa: E402
import paths  # noqa: E402
import selection  # noqa: E402
from faser_genie_rates import NUCLEI_PER_CM2, Z_W, N_W, REF  # noqa: E402

LADDER = f"{paths.POWHEG_V2}/ladder-faser-v2"          # the production ladder
LADDER_RWGT = f"{paths.POWHEG_V2}/ladder-faser"        # the reweighted one
OUT = f"{BASE}/results_nu/faser_powheg_rates.json"
OUT_RWGT = f"{BASE}/results_nu/faser_powheg_rates_rwgt.json"
# The benchmark's own region and the all-y one, evaluated together on every
# ladder point: the Les Houches pass is the cost and the cuts are free.
# >>> TIER S IS IN THE LADDER (2026-09-08), for the electronic-detector
# comparison. <<<  FASER's spectrometer analyses select a muon above 100 GeV
# within 25 mrad of the axis, and that is exactly what this benchmark calls
# Tier S -- so the ladder has to carry it if our prediction is to be compared
# with theirs bin by bin.  It costs nothing: the fraction is computed from the
# SAME Les Houches pass as the other two.
SELS = {name: selection.SELECTIONS[name]
        for name in ("inclusive", "ally", "faser_s")}

# >>> AND THE MUON CUTS WITH NO y WINDOW, WHICH IS THE ONE THAT MATCHES
# FASER'S ACCEPTANCE. <<<  Their Table II alpha is the fraction of ALL
# charged-current interactions in the fiducial volume whose muon passes
# p > 100 GeV and theta < 25 mrad -- there is no y cut anywhere in it.  Tier S
# carries this benchmark's 0.2 < y < 0.9 as well, and at 100-300 GeV that
# window is the DOMINANT loss: E' > 100 GeV at E = 200 GeV already means
# y < 0.5, so y > 0.2 removes the stiffest muons in the bin.  Comparing Tier S
# against their alpha therefore compares two different things and comes out a
# factor of two low.  This region is Tier S minus the y window, and it is
# defined here rather than in selection.py because nothing else needs it: it
# is not a fiducial region of this benchmark, it is FASER's acceptance.
SELS["faser_s_ally"] = selection.Selection(
    "faser_s_ally", selection.Q2_FLOOR, 0.0, 1.0,
    label="Tier S muon cuts, no y window (FASER's acceptance)",
    e_lep_min=100.0, theta_max=0.025,
    note="E' > 100 GeV and theta < 25 mrad and nothing else, so that the "
         "fraction passing is comparable with Table II of arXiv:2412.03186.")

# The proton mass, from the one place the benchmark keeps it.
M_P = beams.M_P


def dot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


def total_pb(rundir):
    for line in open(os.path.join(rundir, "pwg-stat.dat"), errors="replace"):
        if "total (btilde+remnants) cross section in pb" in line:
            v = line.split("pb")[1].split()[0]
            return float(v.replace("D", "E").replace("d", "e"))
    raise SystemExit(f"no total in {rundir}/pwg-stat.dat")


SCALE_IDS = ["1001", "1002", "1003", "1004", "1005", "1006", "1007"]
WGT_RE = re.compile(r"<wgt\s+id=['\"]([^'\"]+)['\"]\s*>\s*([0-9eEdD.+-]+)\s*<")


def fiducial_fraction(lhe, e_beam=None):
    """sum w(fid) / sum w(all) over the Les Houches events, PER REGION.

    Both regions come out of ONE pass, because they are two predicates on the
    same two invariants and the file is the expensive part: the benchmark's
    own (Q2 > 4, 0.2 < y < 0.9) and the all-y region the arXiv:2402.13318
    comparison needs (user, 2026-09-04).  Returns {name: fraction}, n,
    {name: [fraction at each of the 7 scale points]} -- the last from the
    <rwgt> block that tools/powheg_v2_faser_ladder_rwgt.sh adds, and None
    when the file carries none.  Each scale point is normalised to the
    NOMINAL total, so tot * fraction_i is that scale's fiducial
    cross-section and the envelope over i is the MHOU (user, 2026-09-04).

    >>> A SELECTION WITH LEPTON CUTS NEEDS THE LAB FRAME, AND THE LES HOUCHES
    FILE IS NOT IN IT. <<<  Tier S asks for E' > 100 GeV and theta < 25 mrad
    -- FASER's own spectrometer cuts -- and POWHEG-V2 writes the neutrino side
    in the lepton-proton centre of mass.  Both are rebuilt from invariants,
    exactly as analyze.py does: E_lab = p.P_proton / M_p, and the lab angle
    from p_T (unchanged by a boost along z) together with that energy.  The
    incoming record in the file is the PARTON, whose normalisation is
    arbitrary, so the proton is reconstructed from the beam energy instead --
    which is why `e_beam` has to be passed in.  Without it a lepton cut would
    be applied to centre-of-mass numbers and would look perfectly reasonable.
    """
    P_prot = None
    if e_beam is not None:
        ss = M_P * M_P + 2.0 * M_P * float(e_beam)
        rs = math.sqrt(ss)
        P_prot = (0.5 * (ss + M_P * M_P) / rs, 0.0, 0.0,
                  -0.5 * (ss - M_P * M_P) / rs)
    sw = 0.0
    sfid = {name: 0.0 for name in SELS}
    sfid_sc = {name: [0.0] * len(SCALE_IDS) for name in SELS}
    have_sc, n_missing = False, 0
    n = 0
    with open(lhe) as f:
        inev = False
        for line in f:
            if line.startswith("<event"):
                inev, head, parts, wg = True, None, [], {}
                continue
            if not inev:
                continue
            if line.startswith("<wgt"):
                m = WGT_RE.search(line)
                if m:
                    wg[m.group(1)] = float(m.group(2).replace("D", "E"))
                continue
            if line.startswith("</event"):
                inev = False
                w = float(head[2])
                k = P = lep = None
                for p in parts:
                    pid, st = int(p[0]), int(p[1])
                    p4 = (float(p[9]), float(p[6]), float(p[7]), float(p[8]))
                    if st == -1 and abs(pid) in (12, 14, 16):
                        k = p4
                    elif st == -1:
                        # the incoming parton: the proton is along it, and
                        # y = P.q / P.k does not depend on its normalisation
                        P = p4
                    elif st == 1 and abs(pid) in (11, 13, 15):
                        lep = p4
                sw += w
                n += 1
                if k is None or P is None or lep is None:
                    continue
                q = (k[0]-lep[0], k[1]-lep[1], k[2]-lep[2], k[3]-lep[3])
                Q2 = -dot(q, q)
                y = dot(P, q) / dot(P, k)
                e_lab = theta = None
                if P_prot is not None:
                    e_lab = dot(lep, P_prot) / M_P
                    pt = math.hypot(lep[1], lep[2])
                    m2 = max(dot(lep, lep), 0.0)
                    pz2 = max(max(e_lab * e_lab - m2, 0.0) - pt * pt, 0.0)
                    pz = math.sqrt(pz2)
                    theta = math.atan2(pt, pz) if pz > 0 else math.pi / 2
                for name, sel in SELS.items():
                    if not (Q2 > sel.q2_min and sel.y_min < y < sel.y_max):
                        continue
                    if sel.e_lep_min is not None or sel.theta_max is not None:
                        if e_lab is None:
                            raise SystemExit(
                                f"{name} has lepton cuts but no beam energy "
                                f"was given for {lhe} -- refusing to apply "
                                f"them in the wrong frame")
                        if not sel.passes_lepton(e_lab, theta):
                            continue
                    if True:
                        sfid[name] += w
                        if wg:
                            have_sc = True
                            if all(i in wg for i in SCALE_IDS):
                                for j, i in enumerate(SCALE_IDS):
                                    sfid_sc[name][j] += wg[i]
                            else:
                                n_missing += 1
                continue
            if line.startswith("<") or line.startswith("#"):
                continue
            t = line.split()
            if head is None:
                head = t
            elif len(t) >= 10:
                parts.append(t)
    if n_missing:
        # a truncated <rwgt> block is an error, not a smaller band
        raise SystemExit(f"{lhe}: {n_missing} events carry an incomplete "
                         f"<rwgt> block")
    by_scale = ({name: [v / sw for v in vals] for name, vals in sfid_sc.items()}
                if have_sc else None)
    return {name: v / sw for name, v in sfid.items()}, n, by_scale


def _rwgt_scale_ratios():
    """{(cur, tgt, E): {region: [sigma_fid(scale j) / sigma_fid]}} from the
    reweighted ladder's own weights -- the band the production fold imports."""
    with open(OUT_RWGT) as f:
        d = json.load(f)
    out = {}
    for arm, pts in d["ladder"].items():
        cur, tgt = arm.split("_")
        for p in pts:
            sc = p.get("sigma_fid_by_scale_pb")
            if not sc:
                raise SystemExit(f"{OUT}: {arm} E={p['E_GeV']} has no scale "
                                 f"band to import")
            out[(cur, tgt, p["E_GeV"])] = {
                # a region empty at this energy (Tier S at 30 GeV) has no
                # band to import; the point contributes nothing either way
                name: ([v / p["sigma_fid_pb"][name] for v in vals]
                       if p["sigma_fid_pb"][name] > 0 else [1.0] * len(vals))
                for name, vals in sc.items()}
    return out


def ladder(cur, tgt, base=LADDER, import_band=None):
    pts = []
    for d in sorted(glob.glob(f"{base}/{cur}_{tgt}_E*")):
        lhe = os.path.join(d, "pwgevents.lhe")
        if not os.path.exists(lhe):
            continue
        # the reweighted twin carries the scale weights; same events
        rw = os.path.join(d, "rwgt", "pwgevents-rwgt.lhe")
        if os.path.exists(rw):
            lhe = rw
        e = float(os.path.basename(d).split("_E")[1])
        tot = total_pb(d)
        frac, n, by_scale = fiducial_fraction(lhe, e_beam=e)
        sfs = ({name: [tot * f for f in fs] for name, fs in by_scale.items()}
               if by_scale else None)
        if import_band is not None:
            if by_scale:
                raise SystemExit(f"{lhe} carries its own scale weights -- "
                                 f"refusing to overwrite them with the reweighted ladder's")
            rel = import_band.get((cur, tgt, e))
            if rel is None:
                raise SystemExit(f"no reweighted-ladder scale band at {cur}_{tgt}_E{e:g}")
            sfs = {name: [tot * frac[name] * r for r in rel[name]]
                   for name in frac}
        pts.append((e, tot, frac,
                    {name: tot * f for name, f in frac.items()}, n, sfs))
    pts.sort()
    return pts


def rate(pid, e_knots, sig_pb, per_cm2):
    e, phi = fr.flux(str(pid))
    phi = phi * fr.LUMI_FB / 150.0
    s = np.exp(np.interp(np.log(e), np.log(e_knots), np.log(sig_pb),
                         left=-np.inf, right=np.log(sig_pb[-1])))
    return float(np.sum(phi * s * 1e-36 * per_cm2))


def main():
    rwgt = "--rwgt-ladder" in sys.argv[1:]
    base, out_path = (LADDER_RWGT, OUT_RWGT) if rwgt else (LADDER, OUT)
    band = None if rwgt else _rwgt_scale_ratios()
    out = {"reference": REF,
           "selections": {name: sel.as_dict() for name, sel in SELS.items()},
           "selection": "Q2 > 4 GeV2, 0.2 < y < 0.9",
           "ladder": {}, "flavours": {}}
    lad = {}
    for cur in ("nu", "nubar"):
        for tgt in ("p", "n"):
            pts = ladder(cur, tgt, base, band)
            if len(pts) < 3:
                raise SystemExit(f"ladder {cur}_{tgt} has {len(pts)} points -- "
                                 f"run tools/powheg_v2_faser_ladder.sh")
            lad[(cur, tgt)] = pts
            out["ladder"][f"{cur}_{tgt}"] = [
                {"E_GeV": e, "sigma_pb": t, "fid_fraction": f,
                 "sigma_fid_pb": sf, "n_lhe": n,
                 "sigma_fid_by_scale_pb": sfs}
                for e, t, f, sf, n, sfs in pts]
            print(f"{cur} on {tgt}:")
            for e, t, f, sf, n, sfs in pts:
                print(f"  E = {e:7.0f} GeV  sigma = {t:9.4f} pb  "
                      + "  ".join(f"{name} {f[name]:.4f} -> {sf[name]:9.4f} pb"
                                  for name in SELS))
    # closure against the showered nu_mu-on-proton productions
    print("closure of the LHE-fraction ladder against the showered productions "
          "(nu_mu on p):")
    closure = {}
    for e in beams.ENERGIES:
        tag = beams.Beams("nu", e).tag
        suf = "" if e == beams.ANCHOR_ENERGY else f"_{tag}"
        p = f"{BASE}/results_nu/histos_powheg_nu{suf}.json"
        if not os.path.exists(p):
            continue
        prod = json.load(open(p))["sigma_fid_pb"]
        pts = lad[("nu", "p")]
        es = np.array([x[0] for x in pts])
        sf = np.array([x[3]["inclusive"] for x in pts])
        v = float(np.exp(np.interp(np.log(e), np.log(es), np.log(sf))))
        closure[str(e)] = {"ladder_pb": v, "production_pb": prod, "ratio": v / prod}
        print(f"  {e:6.0f} GeV: ladder {v:.4f} pb, production {prod:.4f} pb, "
              f"ratio {v/prod:.4f}")
    out["closure"] = closure
    # rates, in BOTH regions: the benchmark's own and the all-y one that the
    # arXiv:2402.13318 comparison is made on
    tot = {name: {"nue": 0.0, "numu": 0.0} for name in SELS}
    # the same, at each of the 7 scale points -- one index across every
    # flavour, charge and nucleon, so the envelope of the SUM is taken, not
    # the sum of envelopes
    have_scale = all(x[5] is not None for pts in lad.values() for x in pts)
    tot_sc = {name: {"nue": [0.0] * len(SCALE_IDS), "numu": [0.0] * len(SCALE_IDS)}
              for name in SELS}
    print(f"{'flavour':8s} " + " ".join(
        f"{'sig_W(1TeV) ' + name:>26s} {'N ' + name:>12s}" for name in SELS))
    for pid in (12, -12, 14, -14):
        cur = "nu" if pid > 0 else "nubar"
        pp, pn = lad[(cur, "p")], lad[(cur, "n")]
        es = np.array([x[0] for x in pp])
        key = "nue" if abs(pid) == 12 else "numu"
        rec, line = {}, f"{pid:8d} "
        for name in SELS:
            sw = (Z_W * np.array([x[3][name] for x in pp])
                  + N_W * np.array([x[3][name] for x in pn]))
            n_ev = rate(pid, es, sw, NUCLEI_PER_CM2)
            tot[name][key] += n_ev
            s1 = float(np.exp(np.interp(np.log(1000.0), np.log(es), np.log(sw))))
            rec[f"events_fid_{name}"] = n_ev
            rec[f"sigma_W_1TeV_pb_{name}"] = s1
            line += f" {s1:26.2f} {n_ev:12.1f}"
            if have_scale:
                by_sc = []
                for j in range(len(SCALE_IDS)):
                    sw_j = (Z_W * np.array([x[5][name][j] for x in pp])
                            + N_W * np.array([x[5][name][j] for x in pn]))
                    v = rate(pid, es, sw_j, NUCLEI_PER_CM2)
                    by_sc.append(v)
                    tot_sc[name][key][j] += v
                rec[f"events_fid_{name}_by_scale"] = by_sc
        # the benchmark region keeps the historical key names
        rec["events_fid"] = rec["events_fid_inclusive"]
        rec["sigma_W_1TeV_pb"] = rec["sigma_W_1TeV_pb_inclusive"]
        out["flavours"][str(pid)] = rec
        print(line)
    # >>> THE NEUTRON/PROTON RATIO, CHECKED RATHER THAN ASSUMED. <<<
    # faser_rates.tungsten_factor() builds tungsten for EVERY generator out of
    # one sigma_n/sigma_p ratio, taken from the analytic reference, on the
    # assumption that the ratio is generator-independent.  POWHEG-V2 is the
    # one generator here that computes both nucleons itself, so it can test
    # that assumption instead of inheriting it -- and it is the same fiducial
    # region, so the two ratios are directly comparable.
    print()
    print("sigma_n/sigma_p, POWHEG-V2's own against the analytic reference "
          "the other generators inherit:")
    try:
        e_p, s_p = fr.shape("nu", "inclusive", "p")
        e_n, s_n = fr.shape("nu", "inclusive", "n")
        have_ref = True
    except SystemExit as exc:
        print(f"  (no neutron ladder yet: {exc})")
        have_ref = False
    npr = []
    if have_ref:
        pp, pn = lad[("nu", "p")], lad[("nu", "n")]
        for (e, _t, _f, sfp, _n, _s), (_e2, _t2, _f2, sfn, _n2, _s2) in zip(pp, pn):
            r_pw = sfn["inclusive"] / sfp["inclusive"]
            r_ya = float(fr._loglog(np.array([e]), e_n, s_n)[0]
                         / fr._loglog(np.array([e]), e_p, s_p)[0])
            npr.append({"E_GeV": e, "powheg": r_pw, "yadism": r_ya,
                        "ratio": r_pw / r_ya})
            print(f"  E = {e:7.0f} GeV  POWHEG-V2 {r_pw:.4f}  "
                  f"YADISM {r_ya:.4f}  POWHEG/YADISM {r_pw/r_ya:.4f}")
    out["neutron_over_proton"] = npr

    gen = {}
    pg = f"{BASE}/results_nu/faser_genie_rates.json"
    if os.path.exists(pg):
        gen = json.load(open(pg))
    print()
    for key, lab in (("nue", "nu_e+bar"), ("numu", "nu_mu+bar")):
        c = REF[key][0]
        ti, ta = tot["inclusive"][key], tot["ally"][key]
        line = (f"{lab:10s} POWHEG-V2 fiducial {ti:8.0f}  all-y {ta:8.0f}  "
                f"paper (total CC) {c:6.0f}   fid/paper {ti/c:.3f}")
        d = {"events_fid": ti, "events_fid_ally": ta, "paper": c,
             "fid_over_paper": ti / c, "ally_over_paper": ta / c}
        if have_scale:
            # 7-point mu_R, mu_F envelope from POWHEG-V2's own reweighting of
            # every ladder point (tools/powheg_v2_faser_ladder_rwgt.sh)
            for name, suf in (("inclusive", ""), ("ally", "_ally")):
                v = tot_sc[name][key]
                d[f"events_fid{suf}_by_scale"] = v
                d[f"events_fid{suf}_scale_lo"] = min(v)
                d[f"events_fid{suf}_scale_hi"] = max(v)
            d["scale_band"] = ("7-point mu_R, mu_F envelope, POWHEG-V2 "
                               "reweighting of the ladder")
            line += (f"   scale band all-y "
                     f"[{min(tot_sc['ally'][key]):.0f}, "
                     f"{max(tot_sc['ally'][key]):.0f}]")
        if gen:
            g = gen[key]["events"]
            line += (f"   GENIE total {g:6.0f}   fid/GENIE {ti/g:.3f}  "
                     f"all-y/GENIE {ta/g:.3f}")
            d["genie_total"] = g
            d["fid_over_genie"] = ti / g
            d["ally_over_genie"] = ta / g
        out[key] = d
        print(line)
    out["ladder_dir"] = os.path.basename(base)
    if not rwgt:
        out["scale_band_note"] = ("imported from the reweighted ladder's "
                                  "weights as a relative shift per arm, "
                                  "energy and region; the production LHE "
                                  "carry no <rwgt> block")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {os.path.relpath(out_path, BASE)}")


if __name__ == "__main__":
    main()
