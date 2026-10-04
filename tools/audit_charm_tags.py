#!/usr/bin/env python3
"""Audit every charm definition in the benchmark, and what each figure uses.

TODO N3d, at the user's instruction (2026-08-29): "revise globally charm tags
in all generators, we need to take stock we are doing things kosher."

THERE ARE FOUR DEFINITIONS OF "CHARM PRODUCTION" IN THIS PROJECT and they are
not interchangeable:

  out    charm OUT of the hard process        `<gen>_charm`
  any    charm ANYWHERE in the hard process   `<gen>_charmany`
  in     charm-INITIATED (LO, ZM, NC only)    `<gen>_charmin`
  final  a charm HADRON in the final state    `<gen>_charmfinal`
         -- THE BENCHMARK DEFINITION (user, 2026-08-21)

and a fifth that belongs to the reference rather than to a generator:

  YADISM F2_charm, a PARTON-LEVEL structure function.  On the charged current
  at LO in ZM-VFNS this is the CKM-weighted V_cd/V_cs contribution, which
  COUNTS nu cbar -> mu- sbar (charm-initiated, no outgoing charm) and MISSES
  the beam-remnant charm the final-state tag sees.

WHAT THIS TOOL CHECKS
  1. every published charm result records WHICH definition produced it;
  2. everything the report and the paper plots consume uses ONE definition,
     the benchmark's `final`;
  3. the size of the gap between definitions, per generator and per current,
     is measured and printed -- because it is not small and not uniform.

WHAT IT FOUND ON ITS FIRST RUN, and why the audit was worth doing:

  * the two currents stamped the SAME provenance under DIFFERENT KEYS --
    `hard_tag_mode` on the muon side, `charm_tag_mode` on the neutrino side --
    so no cross-current audit could see both without knowing both names;
  * on the NEUTRINO current the benchmark tag exceeds the hard-process tag by
    6.5% (Herwig), 7.9% (Pythia) and 11.9% (Sherpa) -- a spread BETWEEN
    generators of 5.4 points, comparable to the physics differences the charm
    figures quote;
  * on the MUON current the same gap is 0.2-3.3%, an order of magnitude
    smaller.  The currents are NOT equally sensitive to the definition.

Usage: audit_charm_tags.py [--verbose]; non-zero if a check fails.
"""
import glob
import json
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# `tag_mode` is the name now written by both currents; the other
# two are the legacy keys, read so results stored before
# 2026-08-29 still audit.
MODE_KEYS = ("tag_mode", "charm_tag_mode", "hard_tag_mode")


def mode_of(j):
    for k in MODE_KEYS:
        if j.get(k):
            return j[k]
    return None


def sigma(d, key):
    p = f"{BASE}/{d}/histos_{key}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f).get("sigma_fid_pb")


def consumed_keys():
    """Charm result keys the report and the paper plots actually read."""
    src = ""
    for p in [f"{BASE}/analysis/make_report.py"] + sorted(
            glob.glob(f"{BASE}/analysis/paper_plots/*.py")):
        with open(p) as f:
            src += f.read()
    return {k for k in re.findall(r'"([a-z0-9_]*charm[a-z0-9_]*)"', src)
            if os.path.exists(f"{BASE}/results/histos_{k}.json")
            or os.path.exists(f"{BASE}/results_nu/histos_{k}.json")}


def main():
    verbose = "--verbose" in sys.argv
    bad = []

    # ---- 1. every result records its definition ------------------------
    unstamped = []
    for d in ("results", "results_nu"):
        for f in sorted(glob.glob(f"{BASE}/{d}/histos_*charm*.json")):
            b = os.path.basename(f)[7:-5]
            if b.startswith("yadism"):
                continue            # analytic: no tag, by nature
            with open(f) as fh:
                j = json.load(fh)
            if mode_of(j) is None:
                unstamped.append(f"{d}/{b}")
    print(f"1. provenance: {len(unstamped)} charm result(s) do not record "
          f"which definition produced them")
    if unstamped:
        for u in unstamped[:8]:
            print(f"     {u}")
        bad.append("unstamped results")

    # ---- 2. one definition reaches the page ----------------------------
    keys = consumed_keys()
    modes = {}
    for k in sorted(keys):
        if k.startswith("yadism"):
            continue
        for d in ("results", "results_nu"):
            p = f"{BASE}/{d}/histos_{k}.json"
            if os.path.exists(p):
                with open(p) as f:
                    modes.setdefault(mode_of(json.load(f)), []).append(f"{d}/{k}")
    print(f"\n2. what the report and paper plots consume, by definition:")
    for m, ks in sorted(modes.items(), key=lambda kv: str(kv[0])):
        print(f"     {str(m):7s} : {len(ks):3d} key(s)"
              + ("   " + ", ".join(ks[:4]) if verbose else ""))
    non_final = {m: ks for m, ks in modes.items() if m != "final"}
    if non_final:
        print("     !! a published figure or table uses a definition other "
              "than the benchmark's `final`:")
        for m, ks in non_final.items():
            for k in ks:
                print(f"        {k}  [{m}]")
        bad.append("mixed definitions on the page")

    # ---- 3. how far apart the definitions are --------------------------
    print("\n3. the gap between definitions, benchmark `final` over hard-process"
          " `out`:")
    print(f"     {'generator':12s} {'muon NC':>12s} {'neutrino CC':>14s}")
    gaps = {"results": [], "results_nu": []}
    for gen in ("pythia", "herwig", "sherpa_lo", "powheg"):
        cells = []
        for d in ("results", "results_nu"):
            k = "powheg_nu" if (gen == "powheg" and d == "results_nu") else gen
            o, fi = sigma(d, k + "_charm"), sigma(d, k + "_charmfinal")
            if o and fi:
                gaps[d].append(fi / o - 1.0)
                cells.append(f"{100*(fi/o-1):+7.1f}%")
            else:
                cells.append("      --")
        print(f"     {gen:12s} {cells[0]:>12s} {cells[1]:>14s}")
    for d, lab in (("results", "muon NC"), ("results_nu", "neutrino CC")):
        if gaps[d]:
            lo, hi = 100 * min(gaps[d]), 100 * max(gaps[d])
            print(f"     {lab}: {lo:+.1f}% to {hi:+.1f}%, spread "
                  f"{hi-lo:.1f} points between generators")

    print("\n4. the reference is a FIFTH definition and is NOT like-for-like:")
    print("     YADISM F2_charm is parton level; on the CC at LO in ZM-VFNS it")
    print("     is the CKM-weighted V_cd/V_cs contribution, which COUNTS")
    print("     charm-initiated scattering the final-state tag cannot see and")
    print("     MISSES beam-remnant charm the tag does see.")
    for d, lab in (("results", "muon NC"), ("results_nu", "neutrino CC")):
        y, p = sigma(d, "yadism_charm"), sigma(d, "pythia_charmfinal")
        if y and p:
            print(f"     {lab:12s} YADISM LO {y:10.5g}  vs Pythia `final` "
                  f"{p:10.5g}   ({100*(p/y-1):+.1f}%)")

    print()
    if bad:
        print("AUDIT FAILED: " + "; ".join(bad))
        return 1
    print("audit passed: every published charm number records its definition, "
          "and they all use the benchmark's `final`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
