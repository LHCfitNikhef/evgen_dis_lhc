#!/usr/bin/env python3
"""Assemble Herwig's full NLO result from its positive and negative halves.

WHY THIS EXISTS.  Herwig's DIS matrix elements do not deliver an NLO cross
section in one run.  DISBase computes an NLO weight per phase-space point and
then splits it (MatrixElement/DIS/DISBase.cc::NLOWeight):

    Contribution 1  "PositiveNLO"   ->  max(0,  w)
    Contribution 2  "NegativeNLO"   ->  max(0, -w)

so each run is a positive-weight sample of one half, and the physical answer
is the DIFFERENCE.  Herwig's own defaults set Contribution 1 on
PowhegMEDISNC / PowhegMEDISCC, which is why a run "just works" and quietly
reports only the positive half.  Measured here on muon NC:

    positive 36.221 nb, negative 1.228 nb, difference 34.994 nb
    = 0.9969 x YADISM NLO, where the positive half alone reads 1.0318.

The charm channel is where it really bites, because that is where the NLO
correction is large and negative.

WHAT SUBTRACTS AND HOW.  Cross sections and differential distributions
subtract bin by bin, with errors added in quadrature (the two samples are
independent runs).  Means and fractions do NOT: they are ratios, so each is
first multiplied back up by its own sigma, subtracted, and divided by the
combined sigma.

Usage: combine_nlo.py <pos.json> <neg.json> <out.json> ["label"]
"""
import json
import math
import os
import sys
import sys as _sys, os as _os  # noqa: E402
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "analysis"))
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp

# The fields that say WHICH beam and WHICH fiducial region a result describes.
# They are properties of the two input halves, not of the environment the
# combination happens to run in -- see inherit_stamp().
STAMP_KEYS = ("energy_gev", "beam_tag", "two_kP", "selection")
# Paper plots results also carry the NUCLEON they were generated on.  the earlier production
# results have no such key, so both halves read None and the check is a no-op
# there; on it refuses to subtract a neutron negative half from a proton
# positive one.
STAMP_KEYS = STAMP_KEYS + ("target",)


def inherit_stamp(out, pos, neg, out_path):
    """Carry the (selection, energy) stamp from the HALVES, and check them.

    >>> THIS USED TO BE `stamp(out, current="mu")` AND IT WAS WRONG TWICE. <<<

    1. `stamp()` derives the energy from $BENCH_ENERGY, but this script takes
       its data from ARGV.  Run as `combine_nlo.py <pos_400GeV> <neg_400GeV>
       <out_400GeV>` from a shell with no BENCH_ENERGY set, it wrote the
       correct 400 GeV numbers into a file stamped `energy_gev: 1000.0`.  All
       six `*_full_{400GeV,4TeV}` results carried the anchor's stamp.
       Nothing errored, and the numbers were right -- but the stamp is the
       benchmark's tripwire for an energy mix-up, so a wrong stamp does not
       merely mislabel a result, it DISARMS THE ALARM that guards every other
       result.  A file's metadata must come from the same place as its data.
    2. `current="mu"` was hardcoded, so the NEUTRINO combinations were stamped
       as muon beams.  Harmless today only because `beam_tag` and `two_kP`
       happen not to depend on the lepton; it would not stay harmless.

    Both go away by INHERITING rather than re-deriving: the halves were each
    stamped by the analysis run that produced them, from the environment that
    actually generated them.

    And because the two halves are separate runs read from separate argv
    paths, this is also the only place that can catch a positive half being
    subtracted by a negative half from a DIFFERENT energy or selection -- the
    pooled-two-energies failure that has already happened once in this
    benchmark.  It is a hard error.
    """
    for k in STAMP_KEYS:
        pv, nv = pos.get(k), neg.get(k)
        if pv != nv:
            sys.exit(f"HALVES DISAGREE on {k!r}: positive has {pv!r}, "
                     f"negative has {nv!r}. Refusing to subtract a "
                     f"{sys.argv[2]} from a {sys.argv[1]} -- one of them is "
                     f"the wrong sample.")
    if pos.get("energy_gev") is None:
        # a pre-scan half with no stamp: fall back, but say so, because the
        # fallback reads the environment and that is what caused the bug.
        print("WARNING: the input halves carry no energy_gev stamp; falling "
              "back to $BENCH_ENERGY / the anchor. Check the output.")
        return stamp(out)
    for k in STAMP_KEYS:
        out[k] = pos[k]
    return out

# keys that are sigma-weighted ratios: rebuild as (f_pos*s_pos - f_neg*s_neg)/s
RATIO_KEYS = ("mean_n_d_ch", "frac_ge1_d_ch", "frac_ge2_d_ch",
              "mean_n_d_0", "frac_ge1_d_0", "frac_ge2_d_0",
              "mean_n_d_s", "frac_ge1_d_s", "frac_ge2_d_s",
              "mean_dphi", "frac_no_hadron")


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    pos = json.load(open(sys.argv[1]))
    neg = json.load(open(sys.argv[2]))
    out_path = sys.argv[3]
    label = sys.argv[4] if len(sys.argv) > 4 else \
        pos["label"].replace(" (POWHEG, negative)", " (POWHEG)")

    s_pos, s_neg = pos["sigma_fid_pb"], neg["sigma_fid_pb"]
    sigma = s_pos - s_neg
    err = math.hypot(pos.get("sigma_fid_err_pb", 0.0),
                     neg.get("sigma_fid_err_pb", 0.0))

    out = {k: v for k, v in pos.items()
           if k not in ("hists", "means", "label", "inputs")}
    out.update({
        "label": label,
        "sigma_fid_pb": sigma,
        "sigma_fid_err_pb": err,
        "nlo_positive_pb": s_pos,
        "nlo_negative_pb": s_neg,
        "n_parsed": pos.get("n_parsed", 0) + neg.get("n_parsed", 0),
        "n_fiducial": pos.get("n_fiducial", 0) + neg.get("n_fiducial", 0),
        "combination": ("full NLO = Contribution 1 (positive) minus "
                        "Contribution 2 (negative); see combine_nlo.py"),
        "inputs": {"positive": pos.get("inputs"), "negative": neg.get("inputs")},
    })

    # A QUANTITY PRESENT IN ONE HALF ONLY IS AN ERROR, NOT A POSITIVE-HALF
    # RESULT (pipeline review 2026-09-19, finding C6): falling back to the
    # positive half silently published an UNSUBTRACTED number -- for muon
    # charm the negative half is 25-56% of the positive one.
    missing = []
    for k in RATIO_KEYS:
        if k in pos and k in neg and pos[k] is not None and neg[k] is not None:
            out[k] = (pos[k] * s_pos - neg[k] * s_neg) / sigma if sigma else 0.0
        elif pos.get(k) is not None:
            missing.append(f"ratio key {k}")

    means = {}
    for k, v in (pos.get("means") or {}).items():
        w = (neg.get("means") or {}).get(k)
        if v is None:
            continue
        if w is None:
            missing.append(f"mean {k}")
            continue
        means[k] = (v * s_pos - w * s_neg) / sigma if sigma else v
    out["means"] = means

    hists = {}
    for key, hp in (pos.get("hists") or {}).items():
        hn = (neg.get("hists") or {}).get(key)
        if hn is None:
            missing.append(f"histogram {key}")
            continue
        hists[key] = {
            "edges": hp["edges"],
            "dsig": [a - b for a, b in zip(hp["dsig"], hn["dsig"])],
            "err": [math.hypot(a, b) for a, b in zip(hp["err"], hn["err"])],
        }
    out["hists"] = hists
    if missing:
        sys.exit(f"{label}: the negative half lacks {len(missing)} quantity(ies) "
                 f"the positive half has ({', '.join(missing[:6])}"
                 f"{' ...' if len(missing) > 6 else ''}); a positive-only value "
                 f"would be published unsubtracted.  Re-analyse the negative "
                 f"half with the same selection and tag.  Nothing written.")

    inherit_stamp(out, pos, neg, out_path)
    with open(out_path, "w") as f:
        json.dump(out, f)
    print(f"{label}: {sigma:.5g} +- {err:.3g} pb "
          f"(positive {s_pos:.5g} - negative {s_neg:.5g}, "
          f"negative is {100*s_neg/s_pos:.2f}% of it)")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
