#!/usr/bin/env python3
"""The FONLL INCLUSIVE reference: light in ZM-VFNS, heavy in FONLL.

WHY THIS EXISTS (user, 2026-08-27).  The charm FRACTION is
sigma_charm / sigma_inclusive, and once the charm numerator moved to FONLL the
denominator had to move with it.  It did not, so the figures were plotting

    sigma_charm(FONLL, damped) / sigma_inclusive(ZM-VFNS)

which mixes two schemes in one ratio.  The fraction has to be formed inside a
single scheme or it means nothing.

WHAT "FONLL INCLUSIVE" IS.  In FONLL only the HEAVY-quark part of the
structure function changes scheme; the light-quark contributions are the
ZM-VFNS ones and are untouched.  So

    F_incl^FONLL  =  F_light^ZM  +  F_heavy^FONLL
                  =  (F_incl^ZM - F_heavy^ZM)  +  F_heavy^FONLL

and the same identity holds after integration over the fiducial region and
bin by bin in any distribution, because every term is a cross-section built
from the same grid with the same cuts.  That is what this script forms, so it
needs no yadism run of its own -- the three inputs already exist.

WHAT IT DOES NOT DO.  Only CHARM is moved to FONLL: no bottom FONLL reference
is computed in this benchmark, so the bottom contribution stays in ZM inside
`F_light` above.  At FASER kinematics bottom is a per-mille effect on the
inclusive rate, but it IS an approximation and is recorded in the output as
`heavy_in_fonll: ["charm"]` rather than left implicit.

INPUTS MUST BE FROM THE SAME GRID.  The identity is exact only if all three
were computed on one x/Q2 grid at one beam energy.  That is checked, not
assumed: the FONLL cache silently mixed a 1 TeV grid into the 400 GeV and
4 TeV results until 2026-08-27, and the resulting numbers looked plausible
enough to reach the report.

Usage:
  BENCH_ENERGY=1000 analysis/yadism_fonll_inclusive.py [nlo|nnlo]
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")

import beams  # noqa: E402

import selection  # noqa: E402  -- the one Q2 floor
import target  # noqa: E402  -- the target and the TMC mode

ENERGY = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
TAG = beams.Beams("mu", ENERGY).tag
# The energy tag AND the Q2 floor: this script reads ZM results and writes
# FONLL ones, so it has to name the same region on both sides or it would
# silently combine a Q2 > 4 input into a Q2 > 11 output.
#
# AND THE SELECTION, THE TARGET AND THE TMC MODE (2026-09-13), for the same
# reason and in the order the calculators write them,
# histos_<gen>[_<sel>][_<target>][_tmc][_<energy>][_q2min<N>].  Without them a
# q4w3 tungsten TMC-on run would have read the inclusive proton inputs and
# overwritten the inclusive proton FONLL reference.  All "" at the defaults.
SEL = selection.get()
PRE = (("" if SEL.name == selection.DEFAULT else f"_{SEL.name}")
       + target.suffix() + target.tmc_suffix())
SUFFIX = ("" if ENERGY == beams.ANCHOR_ENERGY else f"_{TAG}") \
    + selection.q2_suffix()
# the fields that must agree between the three inputs, when present
STAMPS = ("selection", "energy_gev", "pdf_set", "target", "tmc")


def load(resdir, name):
    p = f"{resdir}/histos_{name}{PRE}{SUFFIX}.json"
    if not os.path.exists(p):
        return None, p
    with open(p) as f:
        return json.load(f), p


def combine(resdir, order):
    """F_incl^ZM - F_charm^ZM + F_charm^FONLL, as a full result record."""
    o = "" if order == "lo" else f"_{order}"
    incl, p1 = load(resdir, f"yadism{o}")
    czm, p2 = load(resdir, f"yadism_charm{o}")
    cfo, p3 = load(resdir, f"yadism_charm{o}_fonll_damp")
    for d, p in ((incl, p1), (czm, p2), (cfo, p3)):
        if d is None:
            return None, f"missing input {os.path.relpath(p, BASE)}"

    # The three inputs must be the SAME calculation apart from the flavour:
    # one region, one beam, one PDF, one TMC mode.  Their names already say
    # so; the stamps are checked as well, since a name is only as right as
    # the run that wrote it.
    # An ABSENT pdf_set/target/tmc means the default (target.provenance
    # writes none at proton / TMC off), so it is compared as that value; an
    # absent selection or energy is an older writer and is not compared.
    dflt = {"pdf_set": target.PDFSET_P, "target": target.DEFAULT, "tmc": 0}
    for k in STAMPS:
        vals = [d.get(k, dflt.get(k)) for d in (incl, czm, cfo)]
        seen = {json.dumps(v, sort_keys=True) for v in vals if v is not None}
        if len(seen) > 1:
            return None, f"inputs disagree on '{k}': {vals}"

    # The identity is exact only on a common grid; the bin edges are the
    # cheapest witness to that, and a mismatch means the inputs are not
    # comparable however plausible their cross-sections look.
    for key in incl.get("hists", {}):
        for other, lab in ((czm, "ZM charm"), (cfo, "FONLL charm")):
            if key not in other.get("hists", {}):
                return None, f"{lab} has no '{key}' histogram"
            if other["hists"][key]["edges"] != incl["hists"][key]["edges"]:
                return None, (f"bin edges differ between the inclusive and "
                              f"{lab} for '{key}' -- these were not computed "
                              f"on the same grid")

    sig = incl["sigma_fid_pb"] - czm["sigma_fid_pb"] + cfo["sigma_fid_pb"]
    out = {
        "generator": f"yadism{o}_fonll_damp",
        "label": "YADISM, FONLL",
        "n_parsed": incl.get("n_parsed"),
        "n_fiducial": incl.get("n_fiducial"),
        "sigma_fid_pb": sig,
        # the analytic inputs carry an integration error, not a statistical
        # one; they add in quadrature since the three runs are independent
        "sigma_fid_err_pb": (incl.get("sigma_fid_err_pb", 0.0) ** 2
                             + czm.get("sigma_fid_err_pb", 0.0) ** 2
                             + cfo.get("sigma_fid_err_pb", 0.0) ** 2) ** 0.5,
        "scheme": "FONLL inclusive = light ZM-VFNS + charm FONLL",
        "heavy_in_fonll": ["charm"],
        "built_from": {"inclusive_zm": os.path.basename(p1),
                       "charm_zm": os.path.basename(p2),
                       "charm_fonll_damp": os.path.basename(p3)},
        "energy_gev": ENERGY,
        "beam_tag": TAG,
        "hists": {},
    }
    # stamped only away from the defaults, so a earlier re-run is byte-identical
    prov = target.provenance(incl.get("pdf_set"))
    if prov or SEL.name != selection.DEFAULT:
        out["selection"] = incl.get("selection", SEL.as_dict())
        out.update(prov or {})
    for key, h in incl.get("hists", {}).items():
        z, f = czm["hists"][key], cfo["hists"][key]
        out["hists"][key] = {
            "edges": h["edges"],
            "dsig": [a - b + c for a, b, c in
                     zip(h["dsig"], z["dsig"], f["dsig"])],
            "err": [(a*a + b*b + c*c) ** 0.5 for a, b, c in
                    zip(h.get("err", [0]*len(h["dsig"])),
                        z.get("err", [0]*len(z["dsig"])),
                        f.get("err", [0]*len(f["dsig"])))],
        }
    # means are re-formed from the combined distribution rather than copied:
    # a mean is not additive, and carrying the ZM one over would quietly
    # describe a different spectrum from the one in `hists`.
    means = {}
    for key in ("Q2", "xbj", "y"):
        h = out["hists"].get(key)
        if not h:
            continue
        edges, dsig = h["edges"], h["dsig"]
        centres = [(edges[i] + edges[i+1]) / 2 for i in range(len(dsig))]
        widths = [edges[i+1] - edges[i] for i in range(len(dsig))]
        tot = sum(d*w for d, w in zip(dsig, widths))
        if tot:
            means[key] = sum(c*d*w for c, d, w
                             in zip(centres, dsig, widths)) / tot
    out["means"] = means
    return out, None


def main():
    order = (sys.argv[1] if len(sys.argv) > 1 else "nlo").lower()
    if order not in ("nlo", "nnlo"):
        sys.exit("FONLL exists at NLO and NNLO only -- the massive "
                 "contribution gamma* g -> c cbar first enters at "
                 "O(alpha_s), so there is no LO counterpart.")
    wrote = 0
    for resdir, cur in ((f"{BASE}/results", "mu"),
                        (f"{BASE}/results_nu", "nu")):
        rec, why = combine(resdir, order)
        if rec is None:
            print(f"  [{cur} {order}] skipped: {why}")
            continue
        path = f"{resdir}/histos_yadism_{order}_fonll_damp{PRE}{SUFFIX}.json"
        with open(path, "w") as f:
            json.dump(rec, f)
        zm = json.load(open(f"{resdir}/histos_yadism_{order}{PRE}{SUFFIX}.json"))
        print(f"  [{cur} {order} {TAG}] inclusive FONLL "
              f"{rec['sigma_fid_pb']:.4f} pb vs ZM "
              f"{zm['sigma_fid_pb']:.4f} pb "
              f"({100*(rec['sigma_fid_pb']/zm['sigma_fid_pb']-1):+.2f}%)")
        print(f"      wrote {os.path.relpath(path, BASE)}")
        wrote += 1
    return 0 if wrote else 1


if __name__ == "__main__":
    sys.exit(main())
