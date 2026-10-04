#!/usr/bin/env python3
"""All charged pions in z, or only the hardest one per event?

User, 2026-09-21: "I don't understand how the yields in 'Charged pions above
z = 0.1' can be larger than the number of events?  For a SIDIS measurement,
one takes the hardest pion in the event right?"

The SIDIS yields of this study count EVERY charged pion in the z bin, which
is the observable collinear factorisation describes:
    dsigma^h/dz = sum_q sigma_q (x) D_q^h(z),
with D_q^h a number density (its integral is the multiplicity, and
sum_h int z D_q^h dz = 1).  It is also what HERMES and COMPASS publish
(multiplicities M^h = dN^h/dz / N_DIS count every hadron).  The HARDEST
pion per event is a different observable, not given by a fragmentation
function.  This script measures how different, on the 1 TeV POWHEG-V2
anchor samples (the only event files kept): per z bin, all charged pions,
the hardest charged pion of each event, and the number of pions above
z = 0.1 per event, in the benchmark region (Q2 > 4, W > 3, no other cut).
Proton and neutron are mixed with the tungsten event weight of the pion
pass (74 p + 110 n).

Output: results_nu/sidis_leading_vs_all.json
Usage:  analysis/sidis_leading_vs_all.py [--max-events N]   (per nucleon)
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import faser_pions as fp                                     # noqa: E402
import analyze_nu                                            # noqa: E402
from analyze import M_P, dot, dis_invariants, lab_energy     # noqa: E402

OUT = f"{BASE}/results_nu/sidis_leading_vs_all.json"
ENERGY = 1000.0
REGION = "sidis_w"
Z = fp.Z_EDGES


def one_nucleon(t, max_events):
    analyze_nu.NU_BEAM_PID = 14
    sel = fp.SELS[REGION]
    nz = len(Z) - 1
    all_z, lead_z = np.zeros(nz), np.zeros(nz)
    mult = np.zeros(6)                    # events with 0..5+ pions above z = 0.1
    sumw = 0.0
    zsum_max = 0.0
    n = 0
    for fn in fp.sample_files("nu", "powheg_nu", ENERGY, t):
        for w, k, P, parts, _d, _h in analyze_nu.parse_hepmc3(fn):
            if n >= max_events:
                break
            n += 1
            w = 1.0 if w is None else w
            best, best_e = None, -1.0
            for pid, p in parts:
                if pid == 13:
                    e = lab_energy(p, P)
                    if e > best_e:
                        best, best_e = p, e
            if best is None:
                continue
            q2, y, xbj, kP = dis_invariants(k, P, best)
            nu_had = kP / M_P - best_e
            w2 = M_P * M_P + 2.0 * M_P * nu_had - q2
            if not sel.passes(q2, y, w2=w2, x=xbj) or nu_had <= 0:
                continue
            zs = [dot(p, P) / M_P / nu_had for pid, p in parts if abs(pid) == 211]
            sumw += w
            zsum_max = max(zsum_max, sum(zs))
            for z in zs:
                i = np.searchsorted(Z, z, side="right") - 1
                if 0 <= i < nz:
                    all_z[i] += w
            if zs:
                i = np.searchsorted(Z, max(zs), side="right") - 1
                if 0 <= i < nz:
                    lead_z[i] += w
            mult[min(sum(1 for z in zs if z >= fp.Z_TAG), 5)] += w
        if n >= max_events:
            break
    print(f"  {t}: {n} events parsed, {sumw:.4g} weight in the region")
    return {"per_event_all": (all_z / sumw).tolist(),
            "per_event_leading": (lead_z / sumw).tolist(),
            "n_above_ztag": (mult / sumw).tolist(), "n_parsed": n,
            "max_sum_z_pions": zsum_max}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-events", type=int, default=200000)
    a = ap.parse_args()
    with open(fp.spectra_path("nu", "powheg_nu", ENERGY, "W")) as f:
        wp = json.load(f)["selections"][REGION]["weight_p"]
    rec = {t: one_nucleon(t, a.max_events) for t in ("p", "n")}
    mix = {k: (wp * np.array(rec["p"][k]) + (1 - wp) * np.array(rec["n"][k])).tolist()
           for k in ("per_event_all", "per_event_leading", "n_above_ztag")}
    out = {"what": "all charged pions vs the hardest charged pion per event, per region event",
           "energy_gev": ENERGY, "region": REGION, "z_edges": Z.tolist(),
           "z_tag": fp.Z_TAG, "weight_p": wp, "tungsten": mix, "per_nucleon": rec}
    with open(OUT, "w") as f:
        json.dump(out, f)
    al, le = np.array(mix["per_event_all"]), np.array(mix["per_event_leading"])
    cut = Z[:-1] >= fp.Z_TAG - 1e-9
    print(f"written {OUT}\nper region event at 1 TeV, tungsten:")
    print(f"  all pions above z=0.1:      {al[cut].sum():.3f}")
    print(f"  hardest pion above z=0.1:   {le[cut].sum():.3f}  (= fraction of events with one)")
    print(f"  events with 0/1/2/3/4/5+ pions above z=0.1: "
          + " ".join(f"{x:.3f}" for x in mix["n_above_ztag"]))
    print(f"  {'z bin':>11} {'all':>8} {'hardest':>8} {'ratio':>6}")
    for i in range(len(Z) - 1):
        r = le[i] / al[i] if al[i] > 0 else float("nan")
        print(f"  {Z[i]:4.2f}-{Z[i+1]:4.2f} {al[i]:8.4f} {le[i]:8.4f} {r:6.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
