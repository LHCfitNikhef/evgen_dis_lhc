#!/usr/bin/env python3
"""Where FASERnu's charged-current events sit in the (Q2, W) plane, from GENIE.

WHAT THIS ANSWERS.  The NLO generators of the benchmark are run above
Q2 = 4 GeV2 and W = 3 GeV, while a FASER measurement that sees only the muon
(the electronic detector) collects every CC event.  The question, raised by
F. Kling in his reading of the paper (App. B, 2026-10-07), is what fraction
of the events lies in each kinematic region:

  DIS            Q2 > 4 GeV2 and W > 2 GeV, split into the benchmark region
                 (W > 3 GeV) and the strip 2 < W < 3 GeV
  soft DIS       Q2 < 4 GeV2 and W > 2 GeV (the shallow-inelastic region)
  non-DIS        W < 2 GeV: resonances, quasi-elastic, and GENIE's
                 non-resonant inelastic background there

The same partition, made with Pythia alone, is Fig. 7.11 of the FPF white
paper, arXiv:2203.05090; Q > 2 GeV there is Q2 > 4 GeV2 here.

INPUTS.  The per-event tables of GENIE's emulsion ladder (default tune, all
CC channels, no cut at generation), genie/fdladder_job_<cur>_<p|n>_E<E>_1/
faserdata_events.npz, written by analysis/faser_emulsion_shapes.py: beam
energy, Q2 and y of every event, from the outgoing lepton in the target rest
frame.  W is rebuilt from the same lepton, W2 = M2 + 2 M y E - Q2, which is
how GENIE's gst ntuple defines it and how an experiment would.  Each nucleon
is weighted by its own spline total (faser_dimuon_cutflow.genie_sigma_tot),
and tungsten per nucleon is (74 p + 110 n)/184.

THE FLUX-WEIGHTED COLUMN is the share of the FASERnu interaction rate,
sum_E phi(E) sigma(E) f(E) / sum_E phi(E) sigma(E), with the vendored Run 3
flux (EPOS-LHC light + POWHEG charm, analysis/faser_rates.flux, E > 10 GeV), the fractions interpolated in
log E between ladder points and frozen below the lowest (30 GeV).

CROSS-CHECK.  At 100 GeV, 400 GeV and 1 TeV the nu_mu fractions are compared
with the dedicated non-DIS samples behind Fig. B.1
(results_nu/genie_nondis_diff_W*.json), an independent GENIE production.

Output: results_nu/genie_kinematic_regions.json
Usage:  analysis/genie_kinematic_regions.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import beams  # noqa: E402
import faser_rates as fr  # noqa: E402
from faser_dimuon_cutflow import genie_sigma_tot  # noqa: E402
from faser_emulsion_shapes import LADDER_ES  # noqa: E402

M_P = beams.M_P
Q2_CUT, W_DIS, W_BENCH = 4.0, 2.0, 3.0
Z_W, N_W = 74, 110
E_LO = 10.0
REGIONS = ["dis_benchmark", "dis_w2to3", "soft_dis", "nondis"]
PIDS = {"nu": 14, "nubar": -14}


def fractions(cur, tgt, e):
    f = f"{BASE}/genie/fdladder_job_{cur}_{tgt}_E{e}_1/faserdata_events.npz"
    z = np.load(f, allow_pickle=False)
    w, q2, y, en = (z[k].astype(float) for k in ("w", "q2", "y", "enu"))
    # the neutron tables quote E as (k.P)/m_p, i.e. 1.0014 x the nominal
    # energy; W below uses only the invariants y (k.P) = P.q and Q2, so the
    # convention changes W2 by m_n^2 - m_p^2 = 0.002 GeV2, which is nothing
    if not np.allclose(en, float(e), rtol=2e-3):
        raise SystemExit(f"{f}: beam energy is not {e} GeV")
    W = np.sqrt(np.maximum(M_P ** 2 + 2.0 * M_P * y * en - q2, 0.0))
    hi_q, hi_w = q2 > Q2_CUT, W > W_DIS
    sel = {"dis_benchmark": hi_q & (W > W_BENCH),
           "dis_w2to3": hi_q & hi_w & (W <= W_BENCH),
           "soft_dis": ~hi_q & hi_w,
           "nondis": ~hi_w}
    tot = w.sum()
    return {r: float(w[s].sum() / tot) for r, s in sel.items()}, len(w)


def tungsten(cur):
    """{E: {region: fraction}} per nucleon on tungsten, and sigma_W(E)."""
    es = np.array([float(e) for e in LADDER_ES])
    out, sig_w, nev = {}, [], {}
    for e in LADDER_ES:
        fp, n_p = fractions(cur, "p", e)
        fn, n_n = fractions(cur, "n", e)
        sp = genie_sigma_tot(PIDS[cur], "p", np.array([float(e)]))[0]
        sn = genie_sigma_tot(PIDS[cur], "n", np.array([float(e)]))[0]
        s = Z_W * sp + N_W * sn
        out[float(e)] = {r: (Z_W * sp * fp[r] + N_W * sn * fn[r]) / s
                         for r in REGIONS}
        sig_w.append(s / (Z_W + N_W))
        nev[float(e)] = n_p + n_n
    return es, out, np.array(sig_w), nev


def flux_weighted(cur, es, fr_e, sig_w):
    e, phi = fr.flux(str(PIDS[cur]))
    keep = e >= E_LO
    e, phi = e[keep], phi[keep]
    sig = np.exp(np.interp(np.log(e), np.log(es), np.log(sig_w)))
    den = np.sum(phi * sig)
    res = {}
    for r in REGIONS:
        f = np.interp(np.log(e), np.log(es), [fr_e[x][r] for x in es])
        res[r] = float(np.sum(phi * sig * f) / den)
    res["mean_energy_gev"] = float(np.sum(phi * sig * e) / den)
    return res


def crosscheck():
    """nu_mu fractions from the Fig. B.1 samples (gst W, all channels)."""
    out = {}
    for e, tag in ((100.0, "_100GeV"), (400.0, "_400GeV"), (1000.0, "")):
        fn = f"{BASE}/results_nu/genie_nondis_diff_W{tag}.json"
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        sig, reg = d["sigma_pb"], d["regions"]
        groups = [g for g in ("DIS", "DIS charm", "RES", "QEL", "DFR")]
        tot = sum(sig[g] for g in groups)
        def share(key):
            return sum(sig[g] * reg[g][key] for g in groups) / tot
        fid, q2w2, w2 = share("fiducial"), share("q2_above_floor_and_w_above_2"), share("w_above_2")
        out[e] = {"dis_benchmark": fid, "dis_w2to3": q2w2 - fid,
                  "soft_dis": w2 - q2w2, "nondis": 1.0 - w2}
    return out


def main():
    res = {"regions": {"dis_benchmark": "Q2 > 4 GeV2, W > 3 GeV",
                       "dis_w2to3": "Q2 > 4 GeV2, 2 < W < 3 GeV",
                       "soft_dis": "Q2 < 4 GeV2, W > 2 GeV",
                       "nondis": "W < 2 GeV"},
           "target": "tungsten per nucleon, (74 p + 110 n)/184",
           "generator": "GENIE default tune, all CC channels (emulsion ladder)",
           "flux": "data/faser_flux_2025 EPOS-LHC light + POWHEG charm (Run 3), E > 10 GeV, "
                   "fractions frozen below the lowest ladder energy"}
    xc = crosscheck()
    for cur in ("nu", "nubar"):
        es, fr_e, sig_w, nev = tungsten(cur)
        res[cur] = {"by_energy": {f"{e:g}": fr_e[e] for e in es},
                    "events_per_energy": {f"{e:g}": nev[e] for e in es},
                    "flux_weighted": flux_weighted(cur, es, fr_e, sig_w)}
        print(f"\n{cur}: tungsten per nucleon, % of the CC cross-section")
        print(f"{'E [GeV]':>10} " + " ".join(f"{r:>14}" for r in REGIONS))
        for e in es:
            print(f"{e:10g} " + " ".join(f"{100*fr_e[e][r]:14.2f}" for r in REGIONS))
        fw = res[cur]["flux_weighted"]
        print(f"{'FASERnu':>10} " + " ".join(f"{100*fw[r]:14.2f}" for r in REGIONS)
              + f"   <E> = {fw['mean_energy_gev']:.0f} GeV")
        if cur == "nu":
            res["crosscheck_nondis_samples"] = {f"{e:g}": v for e, v in xc.items()}
            for e, v in xc.items():
                print(f"{'B.1 ' + format(e, 'g'):>10} "
                      + " ".join(f"{100*v[r]:14.2f}" for r in REGIONS))
    fn = f"{BASE}/results_nu/genie_kinematic_regions.json"
    with open(fn, "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nwrote {os.path.relpath(fn, BASE)}")


if __name__ == "__main__":
    main()
