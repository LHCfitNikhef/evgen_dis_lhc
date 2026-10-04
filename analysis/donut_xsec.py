#!/usr/bin/env python3
"""DONUT's tau-neutrino cross-section against our GENIE, on the same target.

WHAT DONUT MEASURED.  Writing the charged-current cross-section per nucleon as

    sigma_nu_l(E) = sigma_const_l * E * K(E),      l = e, mu, tau

with K the kinematic suppression from the tau-lepton mass (their Fig. 15, and
safely 1 for e and mu over their energy range), DONUT quotes

    sigma_const_tau = 0.72 +- 0.24 (stat) +- 0.36 (syst)  x 1e-38 cm2/GeV

(their Eq. 16), to be compared with 0.51 x 1e-38 for the average of nu_mu and
nu_bar_mu, which is the number they compare it against and which assumes the
equal neutrino and antineutrino fluxes their beam has to 1.05 +- 0.13.  The tau
mass is DIVIDED OUT of their number, so it is directly comparable to ours.

>>> AND OURS IS A nu_mu CALCULATION, WHICH IS THE POINT (user, 2026-09-10): <<<
"neutrino cross-sections are flavour independent, so if we have interaction
cross-sections for muon neutrinos, we can also use them for tau neutrinos".
Exactly so, once the tau mass is out of the way -- and DONUT has already taken
it out.  What is compared here is therefore our sigma_CC per nucleon on IRON,
averaged over neutrino and antineutrino, against both of their numbers.

THE SUM CONVENTION IS CHECKED, NOT ASSUMED.  A gmkspl file for a NUCLEUS holds
one spline per (bound nucleon, channel), and whether those values are per
nucleus or per nucleon is exactly the kind of factor-of-56 that would leave a
plausible-looking curve.  So the summed result is required to land within a
factor of two of the world's isoscalar value, 0.51 x 1e-38 cm2/GeV, which no
wrong convention can do.

Usage: analysis/donut_xsec.py
Writes results_nu/donut_xsec.json
"""
import json
import math
import os
import re
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SPLDIR = f"{BASE}/genie/splines/donut"
RESULTS = f"{BASE}/results_nu"
A_FE = 56
EMAX = 400

# DONUT's own numbers, arXiv:0711.0728 Eq. (16) and Section IX B, in
# 1e-38 cm2/GeV.  Transcribed once, here.
DONUT_TAU = 0.72
DONUT_TAU_STAT = 0.24
DONUT_TAU_SYST = 0.36
DONUT_REF_NUMU = 0.51           # the nu_mu + nubar_mu average they compare to

GEV2_TO_CM2 = 0.3893793656e-27  # 1 GeV^-2 in cm^2


def die(msg):
    sys.exit(f"donut_xsec: {msg}")


def cc_total(pid):
    """(E [GeV], sigma_CC [1e-38 cm2] per NUCLEON) for one beam on Fe-56.

    Only the charged-current splines are summed: the file is built with the
    Default generator list, which also holds every neutral-current channel,
    and DONUT's sigma_const is a charged-current quantity.
    """
    fn = f"{SPLDIR}/{pid}_Fe56_all_e{EMAX}.xml"
    if not os.path.exists(fn):
        die(f"missing {fn} -- run tools/donut_genie.sh splines")
    es = np.logspace(0, math.log10(EMAX), 400)
    tot = np.zeros_like(es)
    n_cc = n_all = 0
    for m in re.finditer(r'<spline name="([^"]+)" nknots="(\d+)">(.*?)</spline>',
                         open(fn).read(), flags=re.S):
        name, body = m.group(1), m.group(3)
        n_all += 1
        if "Weak[CC]" not in name:
            continue
        n_cc += 1
        e = np.array([float(v) for v in re.findall(r"<E>\s*([0-9.eE+-]+)", body)])
        x = np.array([float(v) for v in re.findall(r"<xsec>\s*([0-9.eE+-]+)", body)])
        tot = tot + np.interp(es, e, x, left=0.0, right=x[-1])
    if not n_cc:
        die(f"{fn} holds no charged-current spline (of {n_all})")
    # GENIE tabulates in GeV^-2; 1e-38 cm2 is DONUT's unit
    return es, tot * GEV2_TO_CM2 / 1e-38 / A_FE, n_cc, n_all


def main():
    out = {"what": ("GENIE's charged-current cross-section per nucleon on "
                    "iron against DONUT's tau-neutrino measurement"),
           "target": "Fe-56", "unit": "1e-38 cm2",
           "donut": {"sigma_const_tau": DONUT_TAU,
                     "stat": DONUT_TAU_STAT, "syst": DONUT_TAU_SYST,
                     "reference_numu_average": DONUT_REF_NUMU,
                     "source": "arXiv:0711.0728 Eq. (16) and Section IX B"}}
    es = None
    sig = {}
    for pid, name in ((14, "nu"), (-14, "nubar")):
        e, s, n_cc, n_all = cc_total(pid)
        es = e
        sig[name] = s
        print(f"  {name:6s}: {n_cc} charged-current splines of {n_all}")
    # THE BEAM IS HALF ANTINEUTRINO, measured at 1.05 +- 0.13 and taken as
    # equal by DONUT; the same assumption is made here so the two numbers are
    # the same object.
    avg = 0.5 * (sig["nu"] + sig["nubar"])
    const = avg / es

    # the flux average over DONUT's own interacting spectrum.  FIG. 2 is
    # proportional to flux x sigma, so the FLUX is FIG. 2 divided by sigma and
    # the flux-weighted mean of sigma_const is sum(FIG2) / sum(FIG2 / const).
    fl = json.load(open(f"{BASE}/data/donut/flux.json"))
    edges = np.array(fl["e_edges_gev"])
    ctr = 0.5 * (edges[:-1] + edges[1:])
    w = sum(np.array(v) for v in fl["spectra"].values())
    keep = (ctr >= es[0]) & (ctr <= es[-1]) & (w > 0)
    c_at = np.interp(ctr[keep], es, const)
    flux_avg = float(w[keep].sum() / (w[keep] / c_at).sum())

    # the check that no wrong sum convention can pass
    at100 = float(np.interp(100.0, es, const))
    if not 0.25 < at100 < 1.02:
        die(f"sigma_const comes out {at100:.3g} x 1e-38 cm2/GeV at 100 GeV, "
            f"against the world's isoscalar {DONUT_REF_NUMU}.  The splines "
            f"have been summed per nucleon where they are per nucleus, or "
            f"the other way round.")

    out["energies_gev"] = es.tolist()
    out["sigma_const_nu"] = (sig["nu"] / es).tolist()
    out["sigma_const_nubar"] = (sig["nubar"] / es).tolist()
    out["sigma_const_avg"] = const.tolist()
    out["sigma_const_at_100gev"] = at100
    out["sigma_const_flux_averaged"] = flux_avg
    out["ratio_to_donut_tau"] = flux_avg / DONUT_TAU
    out["ratio_to_reference"] = flux_avg / DONUT_REF_NUMU
    err = math.hypot(DONUT_TAU_STAT, DONUT_TAU_SYST)
    out["pull_vs_donut_tau"] = (flux_avg - DONUT_TAU) / err

    print(f"donut_xsec: GENIE on Fe, nu/nubar averaged")
    print(f"  sigma_const at 100 GeV      = {at100:.3f} x 1e-38 cm2/GeV")
    print(f"  flux-averaged over FIG. 2   = {flux_avg:.3f}")
    print(f"  DONUT nu_tau                = {DONUT_TAU:.2f} "
          f"+- {DONUT_TAU_STAT:.2f} +- {DONUT_TAU_SYST:.2f}  "
          f"-> ours/theirs = {flux_avg/DONUT_TAU:.2f}, "
          f"{out['pull_vs_donut_tau']:+.2f} sigma")
    print(f"  their nu_mu reference       = {DONUT_REF_NUMU:.2f}"
          f"  -> ours/theirs = {flux_avg/DONUT_REF_NUMU:.3f}")
    p = f"{RESULTS}/donut_xsec.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"  wrote {os.path.relpath(p, BASE)}")


if __name__ == "__main__":
    main()
