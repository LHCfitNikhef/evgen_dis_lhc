#!/usr/bin/env python3
"""Prove the selected-subset reweighting is the SAME calculation, not a cheaper one.

THE CLAIM.  A PDF band on a selection is the spread, across members, of the
sum of member weights over the SELECTED events, normalised by a
member-independent constant.  Events failing the selection therefore cannot
influence it, so reweighting only the selected events must return identical
numbers.

WHY IT IS WORTH PROVING RATHER THAN ARGUING.  The saving is enormous -- the
400 GeV dimuon band goes from a measured 32 hours to seconds -- and an
argument that is subtly wrong would show up as a plausible band, which is the
failure mode this benchmark keeps producing.  4 TeV is where BOTH routes are
affordable, so it is the case that can be checked, and it is checked here.

Measured on first run: worst difference 3.1e-13 percentage points, i.e. the
two routes agree to double precision.

Usage: check_subset_reweight.py [energy] [selection]; non-zero on failure.
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))
import beams                                                # noqa: E402
import paths                                                # noqa: E402

CL68 = 100.0 * math.erf(1.0 / math.sqrt(2.0))
TOL = 1e-9          # far above the 3e-13 measured, far below anything physical


def main():
    energy = float(sys.argv[1]) if len(sys.argv) > 1 else 4000.0
    sel = sys.argv[2] if len(sys.argv) > 2 else "faser_dimuon"
    tag = beams.Beams("nu", energy).tag
    full_p = f"{BASE}/results_nu/powheg_selmembers_{sel}_{tag}.json"
    sub_d = f"{paths.POWHEG_V2}/rwgtsel-nu{tag}-{sel}"
    for p in (full_p, f"{sub_d}/weights.npz",
              f"{sub_d}/pwgevents-subset.index.json"):
        if not os.path.exists(p):
            print(f"SKIP -- cannot run: missing {p}")
            return 77           # a skip, not a pass (tools/run_checks.sh)
    with open(full_p) as f:
        full = json.load(f)
    z = np.load(f"{sub_d}/weights.npz", allow_pickle=False)
    ids = [str(x) for x in z["ids"]]
    W = z["weights"]
    col = {i: k for k, i in enumerate(ids)}
    with open(f"{sub_d}/pwgevents-subset.index.json") as f:
        man = json.load(f)
    with open(f"{BASE}/powheg/reweight/rwl_members.json") as f:
        rwl = json.load(f)

    if man["rows"] != full["n_selected"]:
        print(f"FAIL: subset has {man['rows']} rows but the selection kept "
              f"{full['n_selected']} events")
        return 1
    import lhapdf
    lhapdf.setVerbosity(0)
    S = W.sum(axis=0)
    worst, n = 0.0, 0
    for name, rec in rwl["sets"].items():
        w = [i for i in rec["weight_ids"] if i in col]
        if len(w) < 2 or name not in full["sets"]:
            continue
        u = lhapdf.getPDFSet(name).uncertainty([S[col[i]] for i in w], CL68)
        sub = 100.0 * u.errplus / u.central
        ful = 100.0 * full["sets"][name]["rel_plus"]
        worst = max(worst, abs(sub - ful))
        n += 1
        print(f"  {name:24s} full {ful:8.4f}%   subset {sub:8.4f}%   "
              f"{sub-ful:+.2e}")
    if not n:
        print("SKIP -- cannot run: no sets in common")
        return 77
    print(f"\n  worst difference over {n} set(s): {worst:.2e} pp")
    if worst > TOL:
        print(f"FAILED: the subset route is NOT reproducing the full-sample "
              f"band (tolerance {TOL:.0e}).")
        return 1
    print("  the selected-subset reweighting is the same calculation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
