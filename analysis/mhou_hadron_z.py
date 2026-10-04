#!/usr/bin/env python3
"""The 7-point scale band on the charged-hadron z spectrum of the SIDIS yields.

WHAT THIS IS FOR (user, 2026-09-30): the inclusive-hadron yields figure
(paper plot 15) is normalised to POWHEG-V2 and needs "the corresponding MHOU
band, which seems to be missing from the plot".  analysis/mhou_hadron.py
bands PER-EVENT observables (one value per event: a multiplicity, an angle);
the z spectrum is PER HADRON, several entries per event, so it needs the
per-event hadron counts the extraction pass records and nothing in
mhou_hadron.py carries.

THE INPUT is the per-event sidecar analysis/faser_pions.py writes for
POWHEG-V2 at 1 TeV (pions_events_<key>_q4w3_<t>.npz): for every event passing
the yields region (sidis_e == q4w3_faser_e), its file, lhe_index, weight and
the number of charged hadrons inside the emulsion acceptance in each z bin.
Those events are joined to the seven scale weights of their own job by
lhe_index, with mhou_hadron's closure test and its shifted-join control, and
the band is the envelope of the seven z histograms.

WHAT THE BAND CONTAINS.  Each nucleon's seven curves are scaled by ONE
common factor so the nominal integral equals that nucleon's tracked Tier E
cross-section; the variations therefore carry the scale dependence of the
rate and of the hadron spectrum together, as the yields do.  Tungsten is
(74 p + 110 n)/184.  Only the base job of each nucleon is reweighted
(mhou_hadron.py says why that is enough: the band is a ratio of the same
events), so the band is measured on fewer events than the central value.

AT 1 TeV, APPLIED AS A RELATIVE BAND TO THE FLUX-FOLDED YIELDS, the same
approximation as the emulsion figure's fixed-energy piece: the scale weights
exist at the 1 TeV anchor only.

Usage: analysis/mhou_hadron_z.py -> results_nu/mhou_hadron_z_q4w3_faser_e_W.json
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import faser_pions as fp                                     # noqa: E402
from mhou_hadron import NOMINAL_ID, join_scale_weights, sample_layout_key  # noqa: E402

KEY = "powheg_nu"
ENERGY = 1000.0
REGION = "q4w3_faser_e"          # sidis_e, the yields region
SPECIES = "h_emul"
Z, A = 74, 184
OUT = f"{BASE}/results_nu/mhou_hadron_z_{REGION}_W.json"


def nucleon(t):
    side = fp.events_sidecar_path("nu", KEY, ENERGY, t)
    if not os.path.exists(side):
        sys.exit(f"no per-event sidecar {side} -- run faser_pions.py extract "
                 f"--current nu --key {KEY} --energy 1000 --target {t}")
    s = np.load(side, allow_pickle=False)
    files = [str(f) for f in s["files"]]
    job = f"powheg/v2nu_{t}_job_1/events.hepmc"      # relative to BASE
    if job not in files:
        sys.exit(f"{job} is not among the sidecar's files")
    m = s["file_idx"] == files.index(job)
    w = s["w"][m].astype(float)
    idx = s["lhe_index"][m].astype(int)
    cnt = s[SPECIES][m].astype(float)
    wz = np.load(os.path.join(BASE, "results_nu/mhou_weights",
                              f"scale_nu_{t}_job1.npz"))
    ids = [str(x) for x in wz["ids"]]
    F, worst = join_scale_weights(w, idx, wz["weights"], ids, job)
    curves = np.stack([(w * F[k])[:, None] * cnt for k in range(len(ids))]
                      ).sum(axis=1)
    tot = np.array([(w * F[k]).sum() for k in range(len(ids))])
    with open(f"{BASE}/results_nu/histos_{sample_layout_key(KEY, t, REGION)}.json") as f:
        sig = json.load(f)["sigma_fid_pb"]
    scale = sig / tot[ids.index(NOMINAL_ID)]
    print(f"  {t}: {len(w)} selected events joined, closure {worst:.2g}, "
          f"sigma(Tier E) {sig:.6g} pb")
    return ids, curves * scale, tot * scale, len(w), worst


def main():
    ids_p, cp, tp, np_, wp = nucleon("p")
    ids_n, cn, tn, nn, wn = nucleon("n")
    if ids_p != ids_n:
        sys.exit("p and n carry different weight ids")
    ids = ids_p
    inom = ids.index(NOMINAL_ID)
    curves = (Z * cp + (A - Z) * cn) / A
    tot = (Z * tp + (A - Z) * tn) / A
    base = curves[inom]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(base != 0, curves / np.where(base != 0, base, 1.0), 1.0) - 1.0
    others = [k for k in range(len(ids)) if k != inom]
    zero = np.zeros(len(base))
    edges = fp.Z_EDGES
    lo_z = int(np.argmax(edges[:-1] >= fp.Z_TAG - 1e-12))
    zsum = curves[:, lo_z:].sum(axis=1)
    out = {"what": "7-point scale band on the charged-hadron z spectrum "
                   "(all stable charged hadrons inside tan(theta) < 0.5) in "
                   "the Tier E yields region, POWHEG-V2 at 1 TeV, tungsten "
                   "per nucleon (74 p + 110 n)/184",
           "species": SPECIES, "region": REGION, "energy_gev": ENERGY,
           "ids": ids, "nominal_id": NOMINAL_ID,
           "n_events_band": {"p": np_, "n": nn},
           "join_closure": max(wp, wn), "z_edges": edges.tolist(),
           "curves": curves.tolist(),
           "rel_hi": np.max(np.vstack([r[others], zero]), axis=0).tolist(),
           "rel_lo": np.min(np.vstack([r[others], zero]), axis=0).tolist(),
           "zcut_rel_hi": float(max(0.0, max(zsum[k] / zsum[inom] - 1 for k in others))),
           "zcut_rel_lo": float(min(0.0, min(zsum[k] / zsum[inom] - 1 for k in others))),
           "sigma_rel_hi": float(max(tot[k] / tot[inom] - 1 for k in others)),
           "sigma_rel_lo": float(min(tot[k] / tot[inom] - 1 for k in others))}
    with open(OUT, "w") as f:
        json.dump(out, f)
    print(f"hadrons above z = {fp.Z_TAG}: +{100*out['zcut_rel_hi']:.2f} / "
          f"{100*out['zcut_rel_lo']:.2f}%;  Tier E rate: "
          f"+{100*out['sigma_rel_hi']:.2f} / {100*out['sigma_rel_lo']:.2f}%")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
