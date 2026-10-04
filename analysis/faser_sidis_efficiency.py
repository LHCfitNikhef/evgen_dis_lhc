#!/usr/bin/env python3
"""The FASER Tier E efficiency factors, bin by bin in z.

WHY THIS EXISTS (user, 2026-09-21).  "They can compute the SIDIS cross
section with the same Q2 > 4 GeV and W > 3 GeV cuts as we do, but no other
cuts can be applied.  So I would like to evaluate the efficiency factors that
they need to apply to relate their calculation to the actual event yields at
FASER."  So this is the bridge between a theory calculation in the benchmark
region and an event yield in FASERnu:

    dN_FASER/dz  =  eps(z)  x  dN_region/dz

with dN_region/dz the SAME quantity a calculation delivers -- ALL charged
pions (or kaons) of the final state, no angular or energy cut on the hadron,
in events with Q2 > 4 GeV2 and W > 3 GeV and NOTHING ELSE, binned in
z = E_h/nu with both energies in the TARGET REST FRAME.

>>> eps FACTORISES INTO THREE PIECES, AND THEY ARE REPORTED SEPARATELY <<<
because a collaborator may want to replace one of them:

  eps_tier(z)  the FASER Tier E event selection, nested in the region: a
               scattered lepton above 200 GeV at more than 5 mrad, at least
               5 charged tracks above 1 GeV inside tan(theta) < 0.5 of which
               4 inside 0.1, and Delta phi > pi/2 between the lepton and the
               hadronic system.  It is z DEPENDENT even though the cut is not
               made on the hadron: Tier E prefers high-nu events, whose pion
               spectrum is softer in z.
  eps_acc(z)   the emulsion track acceptance ON THE HADRON, tan(theta) < 0.5
               and E_h > 1 GeV in the lab, measured inside Tier E events.
  eps_sel      a flat hadron selection efficiency, faser_pions.
               HADRON_SELECTION_EFF = 0.6 (user, 2026-09-21).

    eps(z) = eps_tier(z) * eps_acc(z) * eps_sel

THE ENERGY DEPENDENCE IS THE INTERESTING PART, and it is why both a
per-energy and a flux-averaged table are written.  Tier E asks for E' > 200
GeV, so its efficiency runs from 0.03 at 300 GeV to 0.48 at 1 TeV; a
calculation done at one beam energy needs eps at THAT energy, while a
calculation already folded over the FASER flux needs the flux-averaged one.
The flux average is over the same flux and the same cross-section ladder the
yields use, separately for nu_mu and for nu_e -- the cross-section is
flavour blind, so the two differ only through the flux.

WHAT IS NOT IN eps.  Nuclear corrections (the target is tungsten by isospin,
74 p + 110 n, with nuclear effects assumed to vanish, on both sides of the
ratio); detector smearing in z; and any efficiency of the vertex or of the
lepton identification beyond the Tier E kinematics above.

Inputs: results_nu/pions_<gen>_q4w3_W[_<E>].json (tools/faser_pions_passes.sh)
Output: results_nu/faser_sidis_efficiency.json
Usage:  analysis/faser_sidis_efficiency.py
"""
import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import faser_pions as fp                                     # noqa: E402

CURRENT = "nu"
REGION_NUM = "sidis_e"      # Tier E nested in the region
REGION_DEN = "sidis_w"      # the region alone: Q2 > 4, W > 3, nothing else
SPECIES = ("pi", "K")
FLAVOURS = (("nu_mu", "14"), ("nu_e", "12"))


def _per_energy(sp, key):
    """eps(z, E) for one generator, at every energy of the SIDIS ladder."""
    out = {}
    for e in sorted(sp):
        sels = sp[e]["selections"]
        if REGION_NUM not in sels or REGION_DEN not in sels:
            continue
        num_sig = sels[REGION_NUM]["sigma_pb"]
        den_sig = sels[REGION_DEN]["sigma_pb"]
        rec = {"sigma_region_pb": den_sig, "sigma_tier_e_pb": num_sig,
               "tier_e_event_efficiency": num_sig / den_sig if den_sig else None}
        for h in SPECIES:
            all_e = np.array(sels[REGION_NUM][f"{h}_all_z"])
            emu_e = np.array(sels[REGION_NUM][f"{h}_emul_z"])
            all_w = np.array(sels[REGION_DEN][f"{h}_all_z"])
            with np.errstate(divide="ignore", invalid="ignore"):
                # the rate ratio, all hadrons: the event selection alone
                tier = np.where(all_w > 0, num_sig * all_e / (den_sig * all_w),
                                np.nan)
                acc = np.where(all_e > 0, emu_e / all_e, np.nan)
            rec[h] = {"eps_tier": tier.tolist(), "eps_acc": acc.tolist(),
                      "eps_total": (tier * acc
                                    * fp.HADRON_SELECTION_EFF).tolist()}
        out[f"{e:g}"] = rec
    return out


def _flux_averaged(sp, key, pid):
    """eps(z) with the ladder folded over the FASER flux, one flavour.

    Both the numerator and the denominator ride the SAME rate machinery --
    faser_pions.flux_weights, which normalises a subset through its enclosing
    region -- so the ratio is the efficiency of one region inside the other
    and not two normalisations divided.  That is the trap recorded in
    faser_pions: two regions normalised by different routes are not
    comparable even when each is right on its own.
    """
    e_num, n_num, how_num = fp.flux_weights(CURRENT, REGION_NUM, key, sp, pid=pid)
    e_den, n_den, how_den = fp.flux_weights(CURRENT, REGION_DEN, key, sp, pid=pid)
    if n_num is None or n_den is None:
        return None
    rec = {"events_tier_e": float(n_num.sum()),
           "events_region": float(n_den.sum()),
           "tier_e_event_efficiency": float(n_num.sum() / n_den.sum()),
           "normalisation": {"tier_e": how_num, "region": how_den}}
    for h in SPECIES:
        emu_e = fp.interp_spectrum(sp, REGION_NUM, f"{h}_emul_z", e_num)
        all_e = fp.interp_spectrum(sp, REGION_NUM, f"{h}_all_z", e_num)
        all_w = fp.interp_spectrum(sp, REGION_DEN, f"{h}_all_z", e_den)
        y_emul_e = (n_num[:, None] * emu_e).sum(axis=0)
        y_all_e = (n_num[:, None] * all_e).sum(axis=0)
        y_all_w = (n_den[:, None] * all_w).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            tier = np.where(y_all_w > 0, y_all_e / y_all_w, np.nan)
            acc = np.where(y_all_e > 0, y_emul_e / y_all_e, np.nan)
        rec[h] = {"eps_tier": tier.tolist(), "eps_acc": acc.tolist(),
                  "eps_total": (tier * acc * fp.HADRON_SELECTION_EFF).tolist(),
                  "yield_region": y_all_w.tolist(),
                  "yield_faser": (y_emul_e * fp.HADRON_SELECTION_EFF).tolist()}
    return rec


def _generator_mean(fa):
    """The flux-averaged factors averaged over the generators (user,
    2026-09-21: "the efficiency [is] computed from the average between the
    five generators considered (the band is kept as it is, namely the
    envelope)").

    eps_total is averaged AS ITSELF, not rebuilt from averaged pieces: the
    mean of a product is not the product of the means, and eps_total is the
    number a collaborator multiplies by.  eps_tier and eps_acc are averaged
    too, for the table, so their product with eps_sel differs from eps_total
    at the per-mille level.  A generator whose bin is empty (NaN) drops out
    of that bin; n_generators records how many entered.
    """
    keys = sorted(fa)
    rec = {"generators": keys,
           "tier_e_event_efficiency": float(np.mean(
               [fa[k]["tier_e_event_efficiency"] for k in keys]))}
    for h in SPECIES:
        r = {}
        for field in ("eps_tier", "eps_acc", "eps_total"):
            a = np.array([fa[k][h][field] for k in keys], dtype=float)
            with warnings.catch_warnings():
                # an all-NaN bin (beyond every generator's reach) stays NaN
                warnings.simplefilter("ignore", RuntimeWarning)
                r[field] = np.nanmean(a, axis=0).tolist()
                if field == "eps_total":
                    r["eps_total_min"] = np.nanmin(a, axis=0).tolist()
                    r["eps_total_max"] = np.nanmax(a, axis=0).tolist()
                    r["n_generators"] = np.isfinite(a).sum(axis=0).tolist()
        rec[h] = r
    return rec


def main():
    gens = [(lab, key) for lab, key in fp.GENERATORS[CURRENT]]
    out = {
        "what": ("FASER Tier E efficiency factors in z: what multiplies a "
                 "calculation done in the benchmark region alone to give the "
                 "FASERnu event yield"),
        "definition": ("dN_FASER/dz = eps(z) x dN_region/dz, with "
                       "dN_region/dz counting ALL charged hadrons of the "
                       "final state (no angular or energy cut) in events "
                       "with Q2 > 4 GeV2 and W > 3 GeV and no other cut, "
                       "z = E_h/nu in the target rest frame"),
        "factorisation": "eps = eps_tier(z) * eps_acc(z) * eps_sel",
        "eps_sel": fp.HADRON_SELECTION_EFF,
        "region": fp.REGION_ALIAS.get(REGION_DEN),
        "region_cuts": "Q2 > 4 GeV2, W > 3 GeV, no cut on y or x",
        "tier_e": ("E'_lepton > 200 GeV, theta_lepton > 5 mrad, at least 5 "
                   "charged tracks with E > 1 GeV inside tan(theta) < 0.5 of "
                   "which at least 4 inside tan(theta) < 0.1, and "
                   "Delta phi > pi/2 between the lepton and the hadronic "
                   "system"),
        "emulsion_acceptance": (f"tan(theta) < {fp.TAN_EMULSION:g} and "
                                f"E_h > {fp.TRACK_E_MIN:g} GeV in the lab"),
        "target": "tungsten per nucleon, 74 p + 110 n (nuclear effects assumed to vanish)",
        "lumi_fb": fp.LUMI_FB,
        "z_edges": fp.Z_EDGES.tolist(),
        "energies": list(fp.ENERGIES),
        "species": list(SPECIES),
        "generators": {key: lab for lab, key in gens},
        "per_energy": {}, "flux_averaged": {},
    }
    for lab, key in gens:
        sp = fp.load_spectra(CURRENT, key)
        if not sp:
            print(f"{key}: no spectra yet")
            continue
        out["per_energy"][key] = _per_energy(sp, key)
        for fl, pid in FLAVOURS:
            rec = _flux_averaged(sp, key, pid)
            if rec:
                out["flux_averaged"].setdefault(fl, {})[key] = rec
    out["flux_averaged_mean"] = {fl: _generator_mean(fa)
                                 for fl, fa in out["flux_averaged"].items()}
    p = f"{BASE}/results_nu/faser_sidis_efficiency.json"
    with open(p, "w") as f:
        json.dump(out, f)
    print(f"written {p}\n")
    report(out)
    return 0


def report(d):
    edges = np.array(d["z_edges"])
    for fl in ("nu_mu", "nu_e"):
        fa = d["flux_averaged"].get(fl, {})
        mean = d.get("flux_averaged_mean", {}).get(fl)
        if not mean:
            continue
        print(f"=== {fl}: Tier E efficiency in z, mean of {len(fa)} "
              f"generators (event efficiency "
              f"{mean['tier_e_event_efficiency']:.4f}) ===")
        print(f"{'z bin':>13} {'eps_tier':>9} {'eps_acc':>8} {'eps':>8} "
              f"{'spread':>16}")
        for h in d["species"]:
            print(f"  -- {h}")
            for i in range(len(edges) - 1):
                vals = [fa[k][h]["eps_total"][i] for k in fa
                        if fa[k][h]["eps_total"][i] == fa[k][h]["eps_total"][i]]
                if not vals:
                    continue
                r = mean[h]
                print(f"  [{edges[i]:4.2f},{edges[i+1]:4.2f}) "
                      f"{r['eps_tier'][i]:9.4f} {r['eps_acc'][i]:8.4f} "
                      f"{r['eps_total'][i]:8.4f} "
                      f"{min(vals):8.4f}-{max(vals):.4f}")
        print()


if __name__ == "__main__":
    sys.exit(main())
