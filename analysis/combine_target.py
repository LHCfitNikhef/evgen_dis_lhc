#!/usr/bin/env python3
"""Tungsten per nucleon from a proton and a neutron result (paper plots).

The generators run each nucleon SEPARATELY (user, 2026-09-13), so the
per-nucleon tungsten result is assembled here:

    sigma_W = (Z sigma_p + N sigma_n) / A,        Z = 74, N = 110, A = 184

Cross-sections and differential distributions add bin by bin with those
weights, errors in quadrature (independent samples).  Means and fractions are
RATIOS, so each is multiplied back up by its own sigma first -- the same rule
herwig7/combine_nlo.py applies to the positive and negative NLO halves.

THE HALVES MUST DESCRIBE THE SAME THING.  Energy, selection, current and
generator are checked to agree, and each input must carry the target it
claims to be (stamped by the analysis run, never inferred from the name).
A proton result paired with a neutron result from another energy is exactly
the pooled-two-samples failure this benchmark has already met once.

Usage: combine_target.py <gen> <mu|nu> [energy ...] [--region R]
       gen is the earlier result key (powheg, sherpa, herwig_nlo_powheg_full,
       powheg_nu, sherpa_nlo, herwig_nlo_full); energies default to
       beams.BENCH_ENERGIES.  Reads histos_<gen>_q4w3_{p,n}[_TAG].json, writes
       histos_<gen>_q4w3_W[_TAG].json beside them.  --region names another
       selection nested in q4w3 (q4w5, q4w3_faser_e, ...; sample_layout).
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
BASE = os.path.dirname(HERE)

import beams                                                # noqa: E402
import target                                               # noqa: E402
import sample_layout                                             # noqa: E402

STAMP_KEYS = ("energy_gev", "beam_tag", "two_kP", "selection")
RATIO_KEYS = ("mean_n_d_ch", "frac_ge1_d_ch", "frac_ge2_d_ch",
              "mean_n_d_0", "frac_ge1_d_0", "frac_ge2_d_0",
              "mean_n_d_s", "frac_ge1_d_s", "frac_ge2_d_s",
              "mean_dphi", "frac_no_hadron")


def path(gen, current, t, energy, region=sample_layout.REGION):
    resdir = "results" if current == "mu" else "results_nu"
    key = beams.at_energy(sample_layout.result_key(gen, t, region), energy)
    return f"{BASE}/{resdir}/histos_{key}.json"


def combine(p, n, gen, region=sample_layout.REGION):
    for k in STAMP_KEYS:
        if p.get(k) != n.get(k):
            sys.exit(f"{gen}: p and n disagree on {k!r} "
                     f"({p.get(k)!r} vs {n.get(k)!r}) -- refusing to combine")
    for d, t in ((p, "p"), (n, "n")):
        if d.get("target") != t:
            sys.exit(f"{gen}: the {t} input is stamped target="
                     f"{d.get('target')!r} -- refusing to combine")
        if (d.get("selection") or {}).get("name") != region:
            sys.exit(f"{gen}: the {t} input is not a {region} result")
    fz, fn = target.Z_W / target.A_W, target.N_W / target.A_W
    sp, sn = p["sigma_fid_pb"], n["sigma_fid_pb"]
    sigma = fz * sp + fn * sn
    err = math.hypot(fz * p.get("sigma_fid_err_pb", 0.0),
                     fn * n.get("sigma_fid_err_pb", 0.0))
    out = {k: v for k, v in p.items()
           if k not in ("hists", "means", "inputs", "tag_stats")}
    out.update({
        "target": "W",
        "per_nucleon": True,
        "Z": target.Z_W, "N": target.N_W,
        "pdfset_nucleons": {"p": p.get("pdfset"), "n": n.get("pdfset")},
        "sigma_fid_pb": sigma,
        "sigma_fid_err_pb": err,
        "sigma_fid_p_pb": sp,
        "sigma_fid_n_pb": sn,
        "n_over_p": sn / sp if sp else None,
        "n_parsed": p.get("n_parsed", 0) + n.get("n_parsed", 0),
        "n_fiducial": p.get("n_fiducial", 0) + n.get("n_fiducial", 0),
        "integrator_closure": {"p": p.get("integrator_closure"),
                               "n": n.get("integrator_closure")},
        "combination": "(74 p + 110 n)/184 per nucleon; "
                       "analysis/combine_target.py",
        "inputs": {"p": p.get("inputs"), "n": n.get("inputs")},
    })
    if p.get("stat_insufficient") or n.get("stat_insufficient"):
        out["stat_insufficient"] = True
        out["stat_insufficient_why"] = "; ".join(
            f"{t}: {d.get('stat_insufficient_why')}" for d, t in
            ((p, "p"), (n, "n")) if d.get("stat_insufficient"))
    for k in RATIO_KEYS:
        if p.get(k) is not None and n.get(k) is not None and sigma:
            out[k] = (fz * sp * p[k] + fn * sn * n[k]) / sigma
    means = {}
    for k, v in (p.get("means") or {}).items():
        w = (n.get("means") or {}).get(k)
        if v is not None and w is not None and sigma:
            means[k] = (fz * sp * v + fn * sn * w) / sigma
    out["means"] = means
    hists = {}
    for key, hp in (p.get("hists") or {}).items():
        hn = (n.get("hists") or {}).get(key)
        if hn is None or hn["edges"] != hp["edges"]:
            continue                 # a histogram on one nucleon only is not W
        hists[key] = {
            "edges": hp["edges"],
            "dsig": [fz * a + fn * b for a, b in zip(hp["dsig"], hn["dsig"])],
            "err": [math.hypot(fz * a, fn * b)
                    for a, b in zip(hp["err"], hn["err"])],
        }
    out["hists"] = hists
    return out


def main():
    if len(sys.argv) < 3 or sys.argv[2] not in ("mu", "nu"):
        sys.exit(__doc__)
    gen, current = sys.argv[1], sys.argv[2]
    args = sys.argv[3:]
    region = sample_layout.REGION
    if "--region" in args:
        i = args.index("--region")
        region = args[i + 1]
        del args[i:i + 2]
    es = [float(a) for a in args] or list(beams.BENCH_ENERGIES)
    rc = 0
    for e in es:
        fp, fn = (path(gen, current, t, e, region) for t in ("p", "n"))
        if not (os.path.exists(fp) and os.path.exists(fn)):
            print(f"{gen} {e:g} GeV: missing "
                  + " and ".join(os.path.relpath(f, BASE)
                                 for f in (fp, fn) if not os.path.exists(f)))
            rc = 1
            continue
        with open(fp) as f:
            p = json.load(f)
        with open(fn) as f:
            n = json.load(f)
        out = combine(p, n, gen, region)
        ofn = path(gen, current, "W", e, region)
        with open(ofn, "w") as f:
            json.dump(out, f)
        print(f"{gen} {e:g} GeV: W {out['sigma_fid_pb']:.6g} "
              f"+- {out['sigma_fid_err_pb']:.3g} pb/nucleon "
              f"(p {p['sigma_fid_pb']:.6g}, n {n['sigma_fid_pb']:.6g}, "
              f"n/p {out['n_over_p'] if out['n_over_p'] is None else round(out['n_over_p'], 4)})"
              f" -> {os.path.relpath(ofn, BASE)}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
