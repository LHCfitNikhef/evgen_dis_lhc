#!/usr/bin/env python3
"""Charged-pion yields in z at FASERnu, for comparison with arXiv:2504.05376.

WHAT THIS IS FOR (user, 2026-09-07): "compare the predictions with our
pipeline for some of the distributions shown in arXiv:2504.05376 ... produce
event yields predictions for the z distribution.  Impose W > 3 GeV and
z > 0.1, in addition to the FASER Tier-E cuts ... produce event yields both
for electron and muon DIS ... predictions for these event rates should be
produced with POWHEG-V2 which is our reference MC program."

Revised the same day, and the revisions ARE the region definitions:

  "we should not use the 0.2 < y < 0.9 cut here, we already have a W cut and
   this is sufficient"                       -> NO y WINDOW anywhere below.
  "we don't cut in x, it is fine to account for this when comparing with
   their numbers but for the final event rate predictions we should not cut
   in x"                                     -> NO x CUT on any yield.
  "the comparison with their numbers should be made at the same energy"
                                             -> a 300 GeV sample, generated.
  "compute which is the percentage of muon-neutrino inclusive DIS events
   which have a pion in the final state"     -> pion_presence() below.

THE PAPER.  Bonino, Gehrmann, Loechner, Schoenwald and Stagnitto, "Identified
Hadron Production in Deeply Inelastic Neutrino-Nucleon Scattering"
(arXiv:2504.05376): the first NNLO QCD calculation of charged pion production
in (anti-)neutrino-induced SIDIS.  Their observable is the MULTIPLICITY

    dM^h/dz = [ d^3 sigma^h / dx dy dz ] / [ d^2 sigma / dx dy ]

-- pions per DIS event per unit z -- in the region x > 0.1 and W > 3 GeV, on
a PROTON, and their Fig. 4 is the FPF-type prediction at E_nu ~ 300 GeV.

>>> THE COMPARISON IS MADE AT THEIR ENERGY, NOT AT OURS. <<<  300 GeV is not
one of the benchmark's energies, so samples were generated for it -- POWHEG-V2
(powheg/cards/POWHEG-V2/nu300GeV-prod, written by
tools/make_energy_cards.py --extra 300) and both GENIE tunes -- and this
module reads them directly.  300 GeV is a COMPARISON POINT and deliberately
not a fourth benchmark energy: it is absent from beams.ENERGIES, so no
cross-section ladder and no FASER rate sees it.  Where a generator has no
300 GeV sample the record says which energy was used instead, rather than
quietly comparing 400 GeV with 300.

>>> THE V2 REGION (user, 2026-09-18). <<<  Every spectrum comes from the
paper-plots samples -- tungsten by isospin, proton and neutron generated
separately, Q2 > 4 GeV2 and W > 3 GeV at generation, NO y window anywhere --
so all five neutrino generators can serve these regions, Herwig and Sherpa
included (their earlier samples carried the window as a generation cut and were
refused here; the guard that refused them, on the measured fraction of a
sample below y = 0.2, stays).  The comparison with the calculation is made
on the PROTON sample alone, as the calculation is on a proton, at its own
300 GeV, a point of the SIDIS ladder (beams.SIDIS_ENERGIES) generated for
every generator.

THE TARGET IS TUNGSTEN, 74 protons and 110 neutrons per nucleus, assembled
per region from the two nucleon passes (faser_pions.combine).

>>> THE STATISTICAL ERROR IS NOT sqrt(N_pions). <<<  Pions in a z bin are not
independent counts: several come from the same event.  With the number of DIS
events Poisson about N_ev and n(b) pions per event in bin b, the bin total is
compound Poisson with mean N_ev <n> and variance N_ev <n^2>, so the relative
error is sqrt(<n^2>/N_ev)/<n> and exceeds the naive 1/sqrt(N_pions) by
sqrt(<n^2>/<n>).  faser_pions.py accumulates that second moment per bin and
this module uses it; the naive number is kept beside it so the size of the
clustering correction stays visible.

Usage:
  analysis/faser_sidis.py            -> results_nu/faser_sidis.json + tables
Requires the tungsten per-event spectra of analysis/faser_pions.py
(tools/faser_pions_passes.sh: extract p and n, then combine).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import beams                                                 # noqa: E402
import faser_pions as fp                                     # noqa: E402
import selection as selmod                                   # noqa: E402

LUMI_FB = fp.LUMI_FB                # 300 fb^-1, the pion tab's luminosity
Z_MIN = fp.Z_TAG                    # 0.1, the hadron-level cut asked for
PAPER = "arXiv:2504.05376"
PAPER_ENERGY = 300.0                # their Fig. 4 beam energy

# The regions, in the order the report shows them.  The chain is nested, so
# each row is the one above plus one cut and every efficiency is meaningful.
REGIONS = [
    ("sidis_w", "Q2 > 4 GeV2, W > 3 GeV (the benchmark region)"),
    ("sidis_e", "+ FASER Tier E"),
    ("sidis_ex", "+ x > 0.1 (diagnostic only)"),
    ("sidis_paper", "Q2 > 4 + W > 3 + x > 0.1, no tier"),
]
YIELD_REGION = "sidis_e"            # what the yields are quoted for
PAPER_REGION = "sidis_paper"        # what the comparison is made in
OUR_REGION = "sidis_w"              # the same without their x floor

FLAVOURS = {
    "nu": [("nu_mu", "14", r"\nu_\mu CC"), ("nu_e", "12", r"\nu_e CC")],
    "mu": [("mu", None, r"\mu^- NC")],
}

# "h" is every stable charged hadron (faser_pions.SPECIES), the observable
# of the inclusive-hadron yields figure (user, 2026-09-30).
SPECIES = ["pi", "pip", "pim", "K", "Kp", "Km", "h"]
REF = {"nu": "powheg_nu", "mu": "powheg"}

# A sample with fewer than this fraction of its events below y = 0.2 was
# generated inside the y window.  0.1% is far below the 19-26% POWHEG and
# GENIE actually have and far above the exactly zero Herwig and Sherpa have,
# so the gate never has to be tuned.
Y_COVERAGE_MIN = 1e-3

REF_DIR = f"{BASE}/data/sidis_bonino2504"
REF_CURVES = [("fig4_LO", "LO"), ("fig4_NLO", "NLO"), ("fig4_NNLO", "NNLO")]


# --------------------------------------------------------- published curves
def load_reference(name):
    """One published curve: z bin edges, central value and scale band."""
    p = f"{REF_DIR}/theory_{name}.dat"
    if not os.path.exists(p):
        return None
    rows = np.loadtxt(p)
    d = {"file": os.path.relpath(p, BASE),
         "z_lo": rows[:, 0].tolist(), "z_hi": rows[:, 2].tolist(),
         "dmdz": rows[:, 3].tolist()}
    if rows.shape[1] >= 6:
        d["dmdz_min"] = rows[:, 4].tolist()
        d["dmdz_max"] = rows[:, 5].tolist()
    return d


def _align(ref):
    """Index into OUR z bins for each of their bins.  Raises if a bin edge
    does not land on ours: a silently rebinned comparison is exactly the
    failure mode CONVENTIONS.md rule 2 is about."""
    edges = np.array(fp.Z_EDGES)
    idx = []
    for lo, hi in zip(ref["z_lo"], ref["z_hi"]):
        j = int(np.argmin(np.abs(edges[:-1] - lo)))
        if not (np.isclose(edges[j], lo) and np.isclose(edges[j + 1], hi)):
            raise SystemExit(
                f"the published bin [{lo:g}, {hi:g}] is not one of ours "
                f"({edges[j]:g}, {edges[j + 1]:g}) -- rebin explicitly")
        idx.append(j)
    return np.array(idx)


# -------------------------------------------------------------- the spectra
def z_bin_mask(edges, z_min):
    """Bins entirely at or above z_min.  Z_EDGES steps by 0.05, so z = 0.1 is
    an edge and the cut is exact; the assert keeps it that way."""
    e = np.asarray(edges)
    assert np.any(np.isclose(e, z_min)), f"z = {z_min} is not a bin edge of {e}"
    return e[:-1] >= z_min - 1e-12


def load_one(current, key, energy, t="p"):
    """The spectra of one nucleon at one energy, or None: the comparison
    with the calculation is made on the proton, as the calculation is."""
    p = fp.spectra_path(current, key, energy, t)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def y_coverage(spectra):
    """The smallest below-y=0.2 fraction over a generator's energies and
    nucleons (the tungsten record carries one entry per nucleon)."""
    fr = []
    for s in spectra.values():
        yc = s.get("y_coverage", {})
        if "frac_below_y02" in yc:
            fr.append(yc["frac_below_y02"])
        else:
            fr.extend(v.get("frac_below_y02", 0.0) for v in yc.values())
    return min(fr) if fr else 0.0


def multiplicity(spectra, region):
    """dM/dz per DIS event at each computed energy: flux independent, so NOT
    convolved.  Every pion counted, with no angular acceptance, because the
    calculation this is compared with has no detector in it."""
    width = np.diff(np.array(fp.Z_EDGES))
    mask = z_bin_mask(np.array(fp.Z_EDGES), Z_MIN)
    out = {}
    for e in sorted(spectra):
        s = spectra[e]["selections"][region]
        m = {"n_selected": s["n_selected"], "sum_w": s["sum_w"]}
        for sp in fp.species_in(spectra):
            per_ev = np.array(s[f"{sp}_all_z"])
            m[f"{sp}_dmdz"] = (per_ev / width).tolist()
            m[f"{sp}_total"] = float(per_ev.sum())
            m[f"{sp}_total_zcut"] = float(per_ev[mask].sum())
        out[f"{e:g}"] = m
    return out


def pion_presence(spectra, region, e_pts, n):
    """The fraction of DIS events that HAVE a charged pion (user, 2026-09-07).

    Three numbers, because "has a pion" has three useful meanings: any
    charged pion at all; one above z = 0.1; and one above z = 0.1 inside the
    emulsion track acceptance, which is the only one a measurement sees.
    Each is 1 - P(0 pions), read off the per-event multiplicity histogram,
    and each is given per energy AND flux averaged with the same weights the
    yields use -- the flux-averaged one is the number for FASER.
    """
    out = {"per_energy": {}}
    es = sorted(spectra)
    kinds = {"any_z": "pi_all_n", "zcut": "pi_all_nz",
             "zcut_emulsion": "pi_emul_nz"}
    curves = {}
    for kind, hist in kinds.items():
        vals = []
        for e in es:
            f0 = spectra[e]["selections"][region][hist][0]
            vals.append(1.0 - f0)
            out["per_energy"].setdefault(f"{e:g}", {})[kind] = 1.0 - f0
        curves[kind] = np.array(vals)
    le = np.log(np.array(es))
    x = np.clip(np.log(e_pts), le[0], le[-1])
    tot = n.sum()
    for kind, v in curves.items():
        f = np.interp(x, le, v)
        out[kind] = float((n * f).sum() / tot) if tot else 0.0
    out["note"] = ("1 - P(no charged pion). any_z: any charged pion in the "
                   f"event; zcut: one above z = {Z_MIN:g}; zcut_emulsion: one "
                   f"above z = {Z_MIN:g} inside tan(theta) < "
                   f"{fp.TAN_EMULSION:g}. The top-level values are flux "
                   "averaged with the same weights as the yields.")
    return out


def one_generator(current, key, label, pid, spectra, nubar=None):
    """Every region's yields, errors and multiplicities for one generator.

    nubar = (key, spectra) of the generator's ANTINEUTRINO row (user,
    2026-10-04: the yields are nu + nubar, as FASER's selection takes both).
    Its yields, with the antineutrino flux of the same flavour, are ADDED to
    the neutrino ones and kept separately as *_nubar; the multiplicities and
    pion-presence fractions stay the neutrino ones."""
    edges = np.array(fp.Z_EDGES)
    mask = z_bin_mask(edges, Z_MIN)
    rec = {"label": label, "energies": sorted(spectra), "regions": {}}
    if nubar is not None:
        rec["nubar_key"], rec["nubar_energies"] = nubar[0], sorted(nubar[1])
    for region, rlabel in REGIONS:
        e_pts, n, how = fp.flux_weights(current, region, key, spectra, pid=pid)
        if n is None:
            continue
        nb = None
        if nubar is not None:
            eb, nbn, _how = fp.flux_weights("nubar", region, nubar[0], nubar[1],
                                            pid="-" + pid)
            if nbn is not None:
                nb = (eb, nbn)
        r = {"label": rlabel,
             "events": float(n.sum()) + (float(nb[1].sum()) if nb else 0.0),
             "events_nu": float(n.sum()),
             "events_nubar": float(nb[1].sum()) if nb else 0.0,
             "normalisation": how,
             "selection": selmod.SELECTIONS[region].as_dict(),
             "v2_selection": fp.REGION_ALIAS.get(region),
             "sigma_pb": {f"{e:g}": spectra[e]["selections"][region]["sigma_pb"]
                          for e in sorted(spectra)}}
        for sp in fp.species_in(spectra):
            # yields: pions inside the emulsion track acceptance, which is
            # where FASER would count them
            m = fp.interp_spectrum(spectra, region, f"{sp}_emul_z", e_pts)
            s2 = fp.interp_spectrum(spectra, region, f"{sp}_emul_z2", e_pts)
            y = (n[:, None] * m).sum(axis=0)
            var = (n[:, None] * s2).sum(axis=0)
            if nb:
                mb = fp.interp_spectrum(nubar[1], region, f"{sp}_emul_z", nb[0])
                s2b = fp.interp_spectrum(nubar[1], region, f"{sp}_emul_z2", nb[0])
                yb = (nb[1][:, None] * mb).sum(axis=0)
                y = y + yb
                var = var + (nb[1][:, None] * s2b).sum(axis=0)
                r[f"{sp}_total_zcut_nubar"] = float(fp.hadron_yield(
                    np.maximum(yb, 0.0))[mask].sum())
            # MC@NLO negative weights (Sherpa) can turn a nearly empty bin at
            # the top of z negative -- a tail fluctuation, not a yield.  Such
            # bins are clipped at zero and counted, so the figure shows an
            # empty bin rather than a NaN error.
            neg = int((y < 0).sum())
            if neg:
                r.setdefault("negative_bins_clipped", {})[sp] = neg
            y = np.maximum(y, 0.0)
            var = np.maximum(var, 0.0)
            # THE HADRON SELECTION EFFICIENCY (user, 2026-09-21), applied
            # once, in faser_pions: it scales the yield and DECORRELATES the
            # hadrons within an event, so the variance is not simply scaled.
            y, var = fp.hadron_yield(y, var)
            r[f"{sp}_z"] = y.tolist()
            r[f"{sp}_z_err"] = np.sqrt(var).tolist()
            r[f"{sp}_z_err_poisson"] = np.sqrt(np.maximum(y, 0.0)).tolist()
            r[f"{sp}_total"] = float(y.sum())
            r[f"{sp}_total_zcut"] = float(y[mask].sum())
            r[f"{sp}_total_zcut_err"] = float(np.sqrt(var[mask].sum()))
            r[f"{sp}_per_event_zcut"] = (float(y[mask].sum() / r["events"])
                                         if r["events"] else 0.0)
        r["pion_presence"] = pion_presence(spectra, region, e_pts, n)
        r["multiplicity"] = multiplicity(spectra, region)
        rec["regions"][region] = r
    return rec


# ------------------------------------------------------------ the comparison
def comparison(out, spectra_at_paper_energy):
    """Our pion multiplicity in the paper's region against their Fig. 4.

    Made at THEIR beam energy: the samples handed in here are the 300 GeV
    ones where they exist, and each generator's record says which energy it
    actually used.  The remaining differences from their setup are our
    Q2 > 4 GeV2 floor, and that ours is a showered generator rather than a
    fixed-order calculation -- which is the point of the comparison.
    """
    nnlo = out["reference_curves"].get("NNLO")
    if not nnlo:
        return None
    idx = _align(nnlo)
    ref = np.array(nnlo["dmdz"])
    zlo = np.array(nnlo["z_lo"])
    zw = np.array(nnlo["z_hi"]) - zlo
    keep = zlo >= Z_MIN - 1e-12
    hi = zlo >= 0.2
    ref_int = float((ref * zw)[keep].sum())
    res = {"note": "our dM(pi+)/dz over the published NNLO curve of Fig. 4, "
                   f"both at E_nu = {PAPER_ENERGY:g} GeV; their region "
                   "(x > 0.1, W > 3 GeV) with our Q2 > 4 GeV2 floor on top; "
                   "both on a proton (the proton sample at 300 GeV)",
           "energy_gev": PAPER_ENERGY,
           "z_lo": nnlo["z_lo"], "z_hi": nnlo["z_hi"],
           "reference_dmdz": nnlo["dmdz"],
           "reference_total_zcut": ref_int,
           "reference_total_zcut_note": f"pi+ per DIS event with z > {Z_MIN:g}, "
                                        "integrating their published NNLO curve",
           "generators": {}}
    for key, (d, e_used) in spectra_at_paper_energy.items():
        mult = multiplicity({e_used: d}, PAPER_REGION)[f"{e_used:g}"]
        nox = multiplicity({e_used: d}, OUR_REGION)[f"{e_used:g}"]
        ours = np.array(mult["pip_dmdz"])[idx]
        ours_nox = np.array(nox["pip_dmdz"])[idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref > 0, ours / ref, np.nan)
        our_int = float((ours * zw)[keep].sum())
        res["generators"][key] = {
            "label": d.get("label") or key,
            "energy_used_gev": e_used,
            "matched_energy": abs(e_used - PAPER_ENERGY) < 1e-6,
            "n_selected": mult["n_selected"],
            "dmdz": ours.tolist(), "ratio": ratio.tolist(),
            "mean_ratio_z_gt_0p2": float(np.nanmean(ratio[hi])),
            "total_zcut": our_int,
            "ratio_total_zcut": our_int / ref_int if ref_int else float("nan"),
            # the same multiplicity on OUR region, i.e. without their x floor,
            # so the cost of that cut to a multiplicity is a measured number
            "total_zcut_no_x": float((ours_nox * zw)[keep].sum()),
        }
    return res


# ------------------------------------------------------------------- build
def build():
    out = {"lumi_fb": LUMI_FB, "hadron_selection_eff": fp.HADRON_SELECTION_EFF,
           "z_edges": fp.Z_EDGES.tolist(), "z_min": Z_MIN,
           "paper": PAPER, "paper_energy_gev": PAPER_ENERGY,
           "paper_title": "Identified Hadron Production in Deeply Inelastic "
                          "Neutrino-Nucleon Scattering",
           "paper_region": "x > 0.1, W > 3 GeV, on a proton; their Fig. 4 is "
                           f"at E_nu ~ {PAPER_ENERGY:g} GeV",
           "regions": [{"key": k, "label": v} for k, v in REGIONS],
           "yield_region": YIELD_REGION, "paper_region_key": PAPER_REGION,
           "reference": REF, "region": fp.REGION, "energies": list(fp.ENERGIES),
           "target": "tungsten, 74 p + 110 n per nucleus, from the proton "
                     "and neutron samples (faser_pions.combine)",
           "y_coverage_min": Y_COVERAGE_MIN,
           "currents": {}}
    paper_spectra = {}
    for current in ("nu", "mu"):
        cur = {"flavours": {}}
        for fname, pid, flabel in FLAVOURS[current]:
            gens, excluded = {}, {}
            for label, key in fp.GENERATORS[current]:
                sp = fp.load_spectra(current, key)
                if not sp:
                    excluded[key] = f"{label}: no spectra on disk"
                    continue
                missing = [e for e in sp if sp[e].get("target") != "W"]
                if missing:
                    raise SystemExit(
                        f"{current} {key}: the spectra at {missing} are not the "
                        "tungsten combination -- re-run tools/faser_pions_passes.sh")
                cov = y_coverage(sp)
                if cov < Y_COVERAGE_MIN:
                    excluded[key] = (
                        f"{label}: {100 * cov:.3f}% of the sample lies below "
                        "y = 0.2, so it was generated inside the 0.2 < y < 0.9 "
                        "window and cannot serve a region without one")
                    continue
                nubar = None
                if current == "nu" and key in fp.NUBAR_OF:
                    sp_nb = fp.load_spectra("nubar", fp.NUBAR_OF[key])
                    if sp_nb:
                        nubar = (fp.NUBAR_OF[key], sp_nb)
                    else:
                        excluded[key] = (f"{label}: no antineutrino spectra yet "
                                         "-- shown as neutrino only")
                rec = one_generator(current, key, label, pid, sp, nubar=nubar)
                if not rec["regions"]:
                    excluded[key] = (f"{label}: no cross-section ladder for the "
                                     "benchmark region, so no rate can be formed")
                    continue
                gens[key] = rec
                if current == "nu" and fname == "nu_mu":
                    at = load_one(current, key, PAPER_ENERGY, "p")
                    if at is not None:
                        at["label"] = label
                        paper_spectra[key] = (at, PAPER_ENERGY)
                    else:
                        d = load_one(current, key, min(sp), "p") or dict(sp[min(sp)])
                        d["label"] = label
                        paper_spectra[key] = (d, min(sp))
            cur["flavours"][fname] = {"pid": pid, "label": flabel,
                                      "generators": gens, "excluded": excluded}
        out["currents"][current] = cur
    out["reference_curves"] = {}
    for name, lab in REF_CURVES:
        c = load_reference(name)
        if c:
            c["order"] = lab
            out["reference_curves"][lab] = c
    out["comparison"] = comparison(out, paper_spectra)
    p = f"{BASE}/results_nu/faser_sidis.json"
    with open(p, "w") as f:
        json.dump(out, f)
    print(f"written {p}")
    return out


# ------------------------------------------------------------------ tables
def report(d):
    zc = d["z_min"]
    for current in ("nu", "mu"):
        for fname, fl in d["currents"][current]["flavours"].items():
            print(f"\n=== {current} / {fname} ({fl['label']}), "
                  f"{d['lumi_fb']:g} fb^-1, z > {zc:g} ===")
            for key, why in fl["excluded"].items():
                print(f"  EXCLUDED {key}: {why}")
            print(f"{'generator':24s} {'region':34s} {'events':>9s} "
                  f"{'pi+':>9s} {'pi-':>9s} {'+-':>7s} {'pi+/pi-':>8s} "
                  f"{'has pi':>7s}")
            for key, g in fl["generators"].items():
                for region, _lab in REGIONS:
                    r = g["regions"].get(region)
                    if not r:
                        continue
                    pp, pm = r["pip_total_zcut"], r["pim_total_zcut"]
                    print(f"{g['label'][:24]:24s} {r['label'][:34]:34s} "
                          f"{r['events']:9.0f} {pp:9.0f} {pm:9.0f} "
                          f"{r['pi_total_zcut_err']:7.0f} "
                          f"{(pp / pm if pm else float('nan')):8.2f} "
                          f"{100 * r['pion_presence']['zcut_emulsion']:6.1f}%")
    c = d.get("comparison")
    if c:
        print(f"\n=== dM(pi+)/dz over the published NNLO Fig. 4 curve "
              f"({PAPER}), E_nu = {c['energy_gev']:g} GeV ===")
        print(f"published integral above z = {zc:g}: "
              f"{c['reference_total_zcut']:.4f} pi+ per DIS event")
        for key, g in c["generators"].items():
            tag = "" if g["matched_energy"] else \
                  f"  [!! at {g['energy_used_gev']:g} GeV, NOT matched]"
            print(f"{g['label'][:26]:26s} {g['total_zcut']:.4f} "
                  f"({g['ratio_total_zcut']:.3f})   mean ratio z > 0.2: "
                  f"{g['mean_ratio_z_gt_0p2']:.2f}   without their x cut "
                  f"{g['total_zcut_no_x']:.4f}{tag}")
    print("\n=== fraction of DIS events with a charged pion ===")
    for current in ("nu", "mu"):
        fl = next(iter(d["currents"][current]["flavours"].values()))
        for key, g in fl["generators"].items():
            for region in ("sidis_w", "sidis_e"):
                r = g["regions"].get(region)
                if not r:
                    continue
                pp = r["pion_presence"]
                print(f"{current:3s} {g['label'][:24]:24s} {region:12s} "
                      f"any z {100 * pp['any_z']:6.2f}%   "
                      f"z > {zc:g} {100 * pp['zcut']:6.2f}%   "
                      f"in the emulsion {100 * pp['zcut_emulsion']:6.2f}%")


def main():
    report(build())


if __name__ == "__main__":
    main()
