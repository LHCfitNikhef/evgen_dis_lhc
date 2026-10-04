#!/usr/bin/env python3
"""PDF and scale uncertainties on the charm fraction, from POWHEG-V2 weights.

WHAT THIS IS FOR (user, 2026-08-28).  The PDF-dependence and MHOU studies on
this page are ANALYTIC -- YADISM structure functions convoluted with each PDF
member and each scale point.  POWHEG-V2 answers the same two questions on the
same process by a completely different route: reweighting the Les Houches
events it already generated.  The goal is not the central values, which are
different objects, but the RELATIVE uncertainties: if the two agree, then
"the PDF choice beats the missing higher orders" rests on two independent
calculations rather than on one.

THE FRACTION'S BAND IS NOT THE TWO BANDS COMBINED.  Numerator and denominator
move together under a PDF variation, and strongly, since they share the same
quark densities.  The ratio is therefore formed PER WEIGHT ID and the
uncertainty taken of the ratio -- doing it the other way round roughly doubles
the band.  The same argument applies to the scale envelope.

THE JOIN IS POSITIONAL, and that is safe for a stated reason: every
reweighting pass reads the same pristine pwgevents.lhe, and POWHEG-V2's
reweighting is a per-event operation that neither drops nor reorders events.
It is CHECKED, not assumed -- a positional join that has slipped produces
perfectly plausible numbers attached to the wrong events, which is how the
first LHE-to-HepMC join here failed.

THE CHARM TAG IS PARTON LEVEL -- an outgoing charm quark in the LHE -- because
what it is compared against, YADISM's F2_charm, is itself parton level.

ON THE SCHEME.  POWHEG-V2 is the massive-charm code (`powheg-cmass`), but
every production card in this benchmark sets `qmass 0d0` to meet the
massless-charm convention, so these samples are ZM-VFNS.  That is a run
setting, not a limitation of the code.

Usage:
  BENCH_ENERGY=1000 analysis/powheg_pdf_members.py
Writes results_nu/powheg_pdf_members_<tag>.json
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import beams                                               # noqa: E402
import paths                                               # noqa: E402
import selection as selection_mod                          # noqa: E402
from analyze import M_P, dis_invariants, lab_energy, lepton_theta  # noqa: E402
from powheg_nlo_uncertainty import read_init, parse_events, lhe_path  # noqa: E402

CL68 = 100.0 * math.erf(1.0 / math.sqrt(2.0))
SCALE_IDS = ["1001", "1002", "1003", "1004", "1005", "1006", "1007"]
NOMINAL_ID = "1001"


def rwgt_dir(energy):
    tag = beams.Beams("nu", energy).tag
    return os.path.join(paths.POWHEG_V2, f"rwgt-nu{tag}")


def main():
    energy = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
    sel = selection_mod.get()
    d = rwgt_dir(energy)
    npz = os.path.join(d, "weights.npz")
    lhe = os.path.join(d, "pwgevents.lhe")
    for p in (npz, lhe):
        if not os.path.exists(p):
            sys.exit(f"missing {p} -- run the reweighting for this energy first")

    z = np.load(npz, allow_pickle=False)
    ids = [str(x) for x in z["ids"]]
    W = z["weights"]                     # (n_events, n_weights)
    col = {i: k for k, i in enumerate(ids)}
    with open(f"{BASE}/powheg/reweight/rwl_members.json") as f:
        rwl = json.load(f)

    print(f"powheg_pdf_members: nu {energy:g} GeV, selection {sel.name}")
    print(f"  {W.shape[0]} events x {W.shape[1]} weights")

    # ---- stream the LHE and accumulate, positionally joined to W ----
    n_ev = n_fid = n_charm = 0
    nw = W.shape[1]
    s_fid = np.zeros(nw)
    s_charm = np.zeros(nw)
    with open(lhe, errors="replace") as fh:
        beam_ids, k, P = read_init(fh, "nu")
        e_lab_beam = (k[0]*P[0] - k[3]*P[3]) / M_P
        if abs(e_lab_beam/energy - 1.0) > 1e-3:
            sys.exit(f"rebuilt beam energy {e_lab_beam:.4f} GeV disagrees "
                     f"with the requested {energy:g} GeV")
        for lep, _w, _nom, has_charm in parse_events(fh, 13):
            if n_ev >= W.shape[0]:
                # the weights cover a CAPPED prefix of a longer sample; stop
                # at the cap rather than sliding the join off its end
                break
            row = W[n_ev]
            n_ev += 1
            Q2, y, xbj, kP = dis_invariants(k, P, lep)
            if Q2 < sel.q2_min or y < sel.y_min or y > sel.y_max:
                continue
            e_lab = lab_energy(lep, P)
            if not sel.passes_lepton(e_lab, lepton_theta(e_lab, lep)):
                continue
            n_fid += 1
            s_fid += row
            if has_charm:
                n_charm += 1
                s_charm += row
    if n_ev != W.shape[0]:
        sys.exit(f"the LHE yielded {n_ev} events but the weight array has "
                 f"{W.shape[0]} rows -- the positional join would be "
                 f"misaligned; refusing")
    print(f"  {n_fid} fiducial, {n_charm} charm-tagged "
          f"({100.0*n_charm/max(n_fid,1):.2f}%)")

    def frac(i):
        c = col.get(i)
        return None if c is None or s_fid[c] == 0 else s_charm[c] / s_fid[c]

    out = {"current": "nu", "energy_gev": energy, "selection": sel.name,
           "beam_tag": beams.Beams("nu", energy).tag,
           "n_events": int(n_ev), "n_fiducial": int(n_fid),
           "n_charm": int(n_charm),
           "confidence_level_percent": CL68,
           "source": os.path.relpath(npz, BASE), "sets": {},
           # THE PER-MEMBER INCLUSIVE WEIGHT SUM, normalised to the nominal
           # member.  It is the DENOMINATOR any selection-level fraction
           # needs, and it can only be formed here: this is the pass that
           # sees every fiducial event with every member weight.  A tier's
           # own module sees a SUBSET, so its "sum over all events" is a sum
           # over the subset and cannot normalise anything.
           #
           # Stored as a ratio to the nominal so it is a pure number, free of
           # the sample size and of the integrator normalisation -- which is
           # exactly what cancels when a fraction is formed member by member,
           # keeping the correlation between numerator and denominator (the
           # same rule pdf_dependence.py follows for the charm fraction).
           "fid_sums_by_id": {}}

    _nomcol = col.get(NOMINAL_ID)
    if _nomcol is not None and s_fid[_nomcol]:
        out["fid_sums_by_id"] = {
            str(i): float(s_fid[c] / s_fid[_nomcol]) for i, c in col.items()}

    # ---- scale envelope, formed on the fraction ----
    fn = frac(NOMINAL_ID)
    fs = [frac(i) for i in SCALE_IDS]
    fs = [v for v in fs if v]
    if fn and fs:
        out["scale"] = {"central": fn, "min": min(fs), "max": max(fs),
                        "rel_up": max(fs)/fn - 1.0,
                        "rel_down": 1.0 - min(fs)/fn,
                        "rel_half_width": 0.5*(max(fs)-min(fs))/fn,
                        "points": len(fs)}
        print(f"  scale envelope on the fraction: {100*fn:.3f}% "
              f"+{100*(max(fs)/fn-1):.2f}% -{100*(1-min(fs)/fn):.2f}% (rel)")

    # ---- PDF bands, per set, through LHAPDF's own prescription ----
    import lhapdf
    lhapdf.setVerbosity(0)
    for name, rec in rwl["sets"].items():
        wids = [i for i in rec["weight_ids"] if i in col]
        vals = [frac(i) for i in wids]
        vals = [v for v in vals if v]
        if len(vals) < 2:
            print(f"  {name:24s} -- no usable members in the weight file")
            continue
        pset = lhapdf.getPDFSet(name)
        try:
            u = pset.uncertainty(list(vals), CL68)
            c, ep, em = u.central, u.errplus, u.errminus
        except Exception as exc:                              # noqa: BLE001
            print(f"    !! LHAPDF uncertainty() failed ({exc}); "
                  f"symmetric Hessian instead")
            c = vals[0]
            ep = em = math.sqrt(sum((v-c)**2 for v in vals[1:]))
        out["sets"][name] = {"fraction": c, "err_plus": ep, "err_minus": em,
                             "members_used": len(vals),
                             "rel_plus": ep/c if c else None,
                             "rel_minus": em/c if c else None,
                             "error_type": pset.errorType}
        print(f"  {name:24s} f = {100*c:.3f} +{100*ep:.3f} -{100*em:.3f} %"
              f"   (rel +{100*ep/c:.2f}% -{100*em/c:.2f}%, "
              f"{len(vals)} members)")

    tag = out["beam_tag"]
    suffix = "" if energy == beams.ANCHOR_ENERGY else f"_{tag}"
    fn_out = f"{BASE}/results_nu/powheg_pdf_members{suffix}.json"
    with open(fn_out, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {fn_out}")


if __name__ == "__main__":
    main()
