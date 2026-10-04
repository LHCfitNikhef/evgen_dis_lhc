#!/usr/bin/env python3
"""Pion and kaon fragmentation functions for the FASER SIDIS study.

User, 2026-09-21: "do a deep dive in the literature and collect
determinations of hadron fragmentation functions (pions and kaons) that we
could use for this study?  I think FASERnu cannot tell apart pi^+ from pi^-
(no magnet) so no need to have FFs which separate hadron charges. ... add to
the report a tab where you compare all these hadron fragmentation functions
(LHAPDF format, likely) for pion and kaon production in the z region
relevant for FASER."

So everything here is CHARGE SUMMED, D^{pi+ + pi-}_i and D^{K+ + K-}_i.
Three storage conventions meet, and each is turned into that sum here:
  "sum"   the set IS the charge sum (NNFF1.0 *sum, MAPFF1.0 *sum, NPC23 *sum)
  "plus"  the set is h+ only (JAM): D^{h+ + h-}_q = D^{h+}_q + D^{h+}_qbar,
          by charge conjugation
The sum is formed MEMBER BY MEMBER, and the uncertainty is then taken on the
sum with the set's own error prescription (LHAPDF PDFSet.uncertainty: the
standard deviation for replicas, the Hessian formula for eigenvectors), at
68% CL.

Flavours shown: u, d, s, c, g (for the charge sum q and qbar are equal), and
a LO illustration of the combination nu_mu CC on tungsten actually probes --
the struck-quark weights after W+ absorption at the region's typical
kinematics (x = 0.23, y = 0.5, POWHEG-V2's means at 300 GeV), with
NNPDF4.0 NNLO:
    u    <- d (|Vud|^2), s (|Vus|^2)
    c    <- d (|Vcd|^2), s (|Vcs|^2)
    dbar <- ubar (|Vud|^2) x (1-y)^2,  sbar <- ubar (|Vus|^2) x (1-y)^2
It is an illustration of which FF the measurement is sensitive to, not a
cross-section.

Output: results_nu/ff_comparison.json
Usage:  analysis/ff_comparison.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

OUT = f"{BASE}/results_nu/ff_comparison.json"
Q2S = (10.0, 100.0)
Z = np.round(np.linspace(0.05, 0.90, 35), 4)
FLAVOURS = ("nucc", "u", "d", "s", "c", "g")
PID = {"u": 2, "d": 1, "s": 3, "c": 4, "g": 21}

# The families.  meta: collaboration label, arXiv, data, order, error type.
# Picked for the comparison: every LHAPDF family with identified charged
# pions and kaons, at every order it is released at.  NPC23's NLO entry is
# its lowQ set (arXiv:2502.17837, QMin 1.3 GeV): the main NLO set of
# arXiv:2407.04422 starts at Q = 4 GeV, above the region's Q = 2 GeV.
SETS = {
    "pi": [
        ("NNFF10_PIsum_nlo", "NNFF1.0", "nlo", "sum"),
        ("MAPFF10NLOPIsum", "MAPFF1.0", "nlo", "sum"),
        ("NPC23_lowQ_PIsum_nlo", "NPC23", "nlo", "sum"),
        ("JAM19FF_pion_nlo", "JAM19", "nlo", "plus"),
        ("JAM20-SIDIS_FF_pion_nlo", "JAM20-SIDIS", "nlo", "plus"),
        ("JAM24_FF_pion_nlo", "JAM24", "nlo", "plus"),
        ("NNFF10_PIsum_nnlo", "NNFF1.0", "nnlo", "sum"),
        ("MAPFF10NNLOPIsum", "MAPFF1.0", "nnlo", "sum"),
        ("NPC23_PIsum_nnlo", "NPC23", "nnlo", "sum"),
    ],
    "K": [
        ("NNFF10_KAsum_nlo", "NNFF1.0", "nlo", "sum"),
        ("MAPFF10NLOKAsum", "MAPFF1.0", "nlo", "sum"),
        ("NPC23_lowQ_KAsum_nlo", "NPC23", "nlo", "sum"),
        ("JAM19FF_kaon_nlo", "JAM19", "nlo", "plus"),
        ("JAM20-SIDIS_FF_kaon_nlo", "JAM20-SIDIS", "nlo", "plus"),
        ("JAM24_FF_kaon_nlo", "JAM24", "nlo", "plus"),
        ("NNFF10_KAsum_nnlo", "NNFF1.0", "nnlo", "sum"),
        ("MAPFF10NNLOKAsum", "MAPFF1.0", "nnlo", "sum"),
        ("NPC23_KAsum_nnlo", "NPC23", "nnlo", "sum"),
    ],
}
REFERENCE = {"nlo": "MAPFF1.0", "nnlo": "MAPFF1.0"}

# the families, for the report's table (data as the papers describe them)
FAMILIES = {
    "NNFF1.0": {"arxiv": "1706.07049", "who": "Bertone, Carrazza, Hartland, Nocera, Rojo (NNPDF)",
                "data": "SIA", "orders": "LO, NLO, NNLO", "errors": "100 MC replicas",
                "method": "neural network"},
    "MAPFF1.0": {"arxiv": "2105.08725, 2204.10331", "who": "Abdul Khalek, Bertone, Khoudli, Nocera (MAP)",
                 "data": "SIA + SIDIS (HERMES, COMPASS)", "orders": "NLO, NNLO",
                 "errors": "200 MC replicas", "method": "neural network"},
    "NPC23": {"arxiv": "2401.02781, 2407.04422, 2502.17837",
              "who": "Gao, Liu, Shen, Xing, Zhao, Zhou",
              "data": "SIA + SIDIS + pp (incl. LHC)", "orders": "NLO, NNLO",
              "errors": "Hessian, 68% CL", "method": "parametric"},
    "JAM19": {"arxiv": "1905.03788", "who": "Sato, Andres, Ethier, Melnitchouk (JAM)",
              "data": "SIA + SIDIS, simultaneous with PDFs", "orders": "NLO",
              "errors": "MC replicas", "method": "parametric, Bayesian MC"},
    "JAM20-SIDIS": {"arxiv": "2101.04664", "who": "Moffat, Melnitchouk, Rogers, Sato (JAM)",
                    "data": "SIA + SIDIS, simultaneous with PDFs", "orders": "NLO",
                    "errors": "MC replicas", "method": "parametric, Bayesian MC"},
    "JAM24": {"arxiv": "2501.00665", "who": "Anderson, Melnitchouk, Sato (JAM)",
              "data": "SIA + SIDIS + LHC W+c, simultaneous with PDFs", "orders": "NLO",
              "errors": "MC replicas", "method": "parametric, Bayesian MC"},
}

# the struck-quark weights of the LO illustration
X_NUCC, Y_NUCC, Q2_NUCC_NOTE = 0.23, 0.5, "x = 0.23, y = 0.5"
VUD2, VUS2, VCD2, VCS2 = 0.9484, 0.0504, 0.0484, 0.9474


def nucc_weights(q2):
    """LO weights of the fragmenting (charge-summed) flavour for nu_mu CC on
    tungsten per nucleon, NNPDF4.0 NNLO at x = X_NUCC."""
    import lhapdf
    p = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
    f = {i: p.xfxQ2(i, X_NUCC, q2) for i in (1, 2, 3, -1, -2, -3)}
    zp, zn = 74.0 / 184.0, 110.0 / 184.0
    d = zp * f[1] + zn * f[2]          # isospin: d_n = u_p
    ub = zp * f[-2] + zn * f[-1]
    s = f[3]
    a = (1.0 - Y_NUCC) ** 2
    w = {"u": d * VUD2 + s * VUS2, "c": d * VCD2 + s * VCS2,
         "d": ub * VUD2 * a, "s": ub * VUS2 * a}
    t = sum(w.values())
    return {k: v / t for k, v in w.items()}


def evaluate(name, conv, q2, weights):
    import lhapdf
    pset = lhapdf.getPDFSet(name)
    members = pset.mkPDFs()
    out = {}
    for fl in FLAVOURS:
        vals = np.zeros((len(members), len(Z)))
        for m, p in enumerate(members):
            for j, z in enumerate(Z):
                def dsum(i):
                    if i == 21:
                        return p.xfxQ2(21, z, q2) * (2.0 if conv == "plus" else 1.0)
                    if conv == "plus":
                        return p.xfxQ2(i, z, q2) + p.xfxQ2(-i, z, q2)
                    return p.xfxQ2(i, z, q2)
                if fl == "nucc":
                    vals[m, j] = sum(w * dsum(PID[k]) for k, w in weights.items())
                else:
                    vals[m, j] = dsum(PID[fl])
        cen, lo, hi = [], [], []
        for j in range(len(Z)):
            u = pset.uncertainty(list(vals[:, j]), 68.27)
            cen.append(u.central)
            lo.append(u.central - u.errminus)
            hi.append(u.central + u.errplus)
        out[fl] = {"central": cen, "lo": lo, "hi": hi}
    return out


def main():
    import lhapdf
    lhapdf.setVerbosity(0)
    res = {"what": ("charge-summed pion and kaon fragmentation functions "
                    "z D_i^{h+ + h-}(z, Q2), central and 68% CL"),
           "z": Z.tolist(), "q2": list(Q2S), "flavours": list(FLAVOURS),
           "z_window_faser": [0.1, 0.8], "reference": REFERENCE,
           "families": FAMILIES,
           "nucc_note": ("LO struck-quark weights, nu_mu CC on tungsten per nucleon, "
                         f"{Q2_NUCC_NOTE}, NNPDF4.0 NNLO"),
           "nucc_weights": {f"{q2:g}": nucc_weights(q2) for q2 in Q2S},
           "sets": {}}
    for h, sets in SETS.items():
        for name, fam, order, conv in sets:
            info = lhapdf.getPDFSet(name)
            rec = {"hadron": h, "family": fam, "order": order, "storage": conv,
                   "members": info.size, "error_type": info.errorType,
                   "qmin": float(info.get_entry("QMin")),
                   "xmin": float(info.get_entry("XMin")),
                   "values": {}}
            for q2 in Q2S:
                rec["values"][f"{q2:g}"] = evaluate(name, conv, q2,
                                                    res["nucc_weights"][f"{q2:g}"])
            res["sets"][name] = rec
            v = rec["values"]["10"]["u"]["central"]
            print(f"{name:26s} {info.size:4d} members; zD_u(z=0.5, Q2=10) = "
                  f"{np.interp(0.5, Z, v):.4f}")
    with open(OUT, "w") as f:
        json.dump(res, f)
    print(f"written {OUT}")


if __name__ == "__main__":
    main()
