#!/usr/bin/env python3
"""POWHEG-V2 SIDIS in the benchmark region alone: the numbers an analytic
NLO calculation is benchmarked against.

User, 2026-09-21: "they will benchmark the NLO calculation with ours. Can you
add a subtab to this SIDIS study tab the expected event yields, obtained with
POWHEG-V2 using the benchmark settings, where only the Q^2 > 4 GeV2 and
W > 3 GeV cuts are applied, but not the FASER Tier E ones?  Then I can use the
efficiencies that you already computed to transform their NLO analytic
calculation based on fragmentation functions to the FASER event yields."

So the region here is EXACTLY the denominator of faser_sidis_efficiency.py:
events with Q2 > 4 GeV2 and W > 3 GeV and nothing else, counting ALL charged
hadrons of the final state with no cut on the hadron, z = E_h/nu in the
target rest frame, tungsten per nucleon (74 p + 110 n, no nuclear effects).
No hadron selection efficiency either: that 0.6 lives in eps.  Two views:

  fixed_energy   dsigma/dz per nucleon at every energy of the SIDIS ladder,
                 in pb per z bin -- sigma_region x the per-event count in the
                 bin.  This is the like-for-like comparison with a
                 fragmentation-function calculation, before any flux.
  flux_folded    dN/dz at the luminosity of the pion study, over the FASERnu
                 nu_mu and nu_e fluxes, through the SAME faser_pions
                 machinery as the efficiency factors' denominator -- so
                 eps(z) x these yields reproduces the FASER yields exactly.

Species: pi, K (charge summed) and pip, pim, Kp, Km.

The fixed-energy MC error is approximate: sigma x sqrt((<n^2> - <n>^2)/N)
from the stored second moment, ignoring the (near-uniform) POWHEG weights.

Output: results_nu/faser_sidis_region_yields.json
Usage:  analysis/faser_sidis_region_yields.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import faser_pions as fp                                     # noqa: E402

CURRENT = "nu"
KEY = "powheg_nu"
REGION = "sidis_w"          # Q2 > 4, W > 3, nothing else
SPECIES = ("pi", "K", "pip", "pim", "Kp", "Km")
FLAVOURS = (("nu_mu", "14"), ("nu_e", "12"))
OUT = f"{BASE}/results_nu/faser_sidis_region_yields.json"


def fixed_energy(sp):
    out = {}
    for e in sorted(sp):
        sel = sp[e]["selections"].get(REGION)
        if sel is None:
            continue
        sig, n = sel["sigma_pb"], sel["n_selected"]
        rec = {"sigma_region_pb": sig, "sigma_region_err_pb": sel["sigma_err_pb"],
               "n_events": n}
        for h in SPECIES:
            m1 = np.array(sel[f"{h}_all_z"])
            m2 = np.array(sel[f"{h}_all_z2"])
            err = sig * np.sqrt(np.clip(m2 - m1 ** 2, 0.0, None) / max(n, 1))
            rec[h] = {"dsigma_pb": (sig * m1).tolist(),
                      "dsigma_err_pb": err.tolist(),
                      "multiplicity": float(sel[f"{h}_all_sum"])}
        out[f"{e:g}"] = rec
    return out


def flux_folded(sp, pid):
    e_pts, n_ev, how = fp.flux_weights(CURRENT, REGION, KEY, sp, pid=pid)
    rec = {"events_region": float(n_ev.sum()), "normalisation": how,
           "flux_energies_gev": np.asarray(e_pts).tolist(),
           "events_per_flux_point": np.asarray(n_ev).tolist()}
    measured = sorted(sp)
    rec["events_below_lowest_sample"] = float(n_ev[e_pts < measured[0]].sum())
    rec["events_above_highest_sample"] = float(n_ev[e_pts > measured[-1]].sum())
    for h in SPECIES:
        shape = fp.interp_spectrum(sp, REGION, f"{h}_all_z", e_pts)
        rec[h] = (n_ev[:, None] * shape).sum(axis=0).tolist()
    return rec


def main():
    sp = fp.load_spectra(CURRENT, KEY)
    if not sp:
        raise SystemExit(f"no {KEY} spectra")
    out = {
        "what": ("POWHEG-V2 + Pythia 8 SIDIS in the benchmark region alone "
                 "(Q2 > 4 GeV2, W > 3 GeV, no other cut): all charged "
                 "hadrons, no cut on the hadron, no selection efficiency"),
        "generator": KEY, "region": fp.REGION_ALIAS.get(REGION),
        "z": "E_h / nu, both in the target rest frame",
        "target": "tungsten per nucleon, (74 sigma_p + 110 sigma_n) / 184, no nuclear effects",
        "z_edges": fp.Z_EDGES.tolist(),
        "lumi_fb": fp.LUMI_FB,
        "sample_energies_gev": [float(e) for e in sorted(sp)],
        "fixed_energy": fixed_energy(sp),
        "flux_folded": {fl: flux_folded(sp, pid) for fl, pid in FLAVOURS},
    }
    with open(OUT, "w") as f:
        json.dump(out, f)
    print(f"written {OUT}")
    z = np.array(out["z_edges"])
    for fl, r in out["flux_folded"].items():
        cut = z[:-1] >= 0.1 - 1e-12
        print(f"{fl}: {r['events_region']:.0f} region events at "
              f"{out['lumi_fb']:g} fb^-1 ({r['events_below_lowest_sample'] / r['events_region']:.1%} "
              f"below the lowest sample energy); pi {sum(r['pi']):.0f} "
              f"(z>0.1: {np.array(r['pi'])[cut].sum():.0f}), "
              f"K {sum(r['K']):.0f} (z>0.1: {np.array(r['K'])[cut].sum():.0f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
