#!/usr/bin/env python3
"""Every cross-energy figure that exists must appear in the report.

WHY (user, 2026-08-29): "double check that no results / plots / comparisons
are missing there".  Twice in this project a figure was built, committed and
never rendered -- cmp_unc_methods.png and cmp_pdf_dependence_mc_E.png, both of
them things the user had explicitly asked for days earlier.  Neither failed;
they simply were not referenced, so nothing complained.

The check is on the RENDERED PAGE, not on the source: a figure can be
referenced by code that never runs, and asking the HTML is the only way to
know it arrived.  Per-energy and per-selection figures are excluded -- there
are hundreds and they are generated in bulk by one driver -- so this covers
the cross-energy and study figures, which are the ones added by hand and
therefore the ones that get forgotten.

Usage: check_report_complete.py; exits non-zero if a figure is missing.
"""
import glob
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = f"{BASE}/results/report.html"

# figures generated in bulk per (energy, selection); not this tool's business
BULK = re.compile(r"_(400GeV|4TeV|faser_s|faser_e|faser_dimuon)\.png$")
# the per-view comparison figures, likewise
PERVIEW = re.compile(r"^cmp_(shw|qed|lo|nlo|orders|shower)[a-z0-9_]*_"
                     r"(Q2|y|x|nu|W|Emu|Elead|nch|nch1|dphi|Ed_|nd_|theta)")

# DELIBERATELY NOT SHOWN, each with the decision behind it.  A figure absent
# on purpose and a figure absent by oversight look identical on disk, so the
# difference is written down rather than left to memory.
EXCLUDED = {
    "cmp_me_": "matrix-element-level views, removed from the report by user "
               "decision on 2026-08-27: with QED radiation off the shower "
               "cannot touch the scattered lepton, so the six DIS-kinematics "
               "observables are identical before and after showering "
               "(measured: Pythia <Q2> 15.258 vs 15.263, Sherpa 15.123 vs "
               "15.095) and the sections measured the same thing twice. "
               "The figures are still generated; only the display was dropped.",
    "cmp_sigma_vs_E.png": "superseded 2026-08-27: the figure was split into "
                          "_lo and _nlo on user request, and both are shown.",
    "cmp_faser_electronic.png": "superseded in the report 2026-09-21 by paper "
                                "plot 12b (pp_faser_electronic.png, current "
                                "POWHEG-V2 ladder); kept because the user's "
                                "FASER meeting deck (talks/) still reads it.",
}


def excluded(name):
    for k, why in EXCLUDED.items():
        if name.startswith(k) or name == k:
            return why
    return None


def main():
    man = f"{BASE}/results/report_figures.json"
    if not os.path.exists(man):
        print("no figure manifest -- rebuild the report first "
              "(analysis/make_report.py writes it)")
        return 0
    import json
    with open(man) as f:
        embedded = {os.path.basename(p) for p in json.load(f)}

    missing, skipped = [], []
    for d in ("results", "results_nu"):
        for p in sorted(glob.glob(f"{BASE}/{d}/*.png")):
            b = os.path.basename(p)
            if BULK.search(b) or PERVIEW.match(b):
                continue
            if b.startswith("paper_"):
                continue            # paper plots have their own tab and check
            if excluded(b):
                skipped.append(b)
                continue
            if b not in embedded:
                missing.append(f"{d}/{b}")
    if skipped:
        seen = sorted({excluded(b)[:60] for b in skipped})
        print(f"{len(skipped)} figure(s) deliberately not shown:")
        for w in seen:
            print(f"   - {w}...")
        print()
    print(f"cross-energy and study figures on disk that are NOT in the "
          f"rendered report: {len(missing)}")
    for m in missing:
        print(f"   {m}")
    if missing:
        print("\nEither reference them from make_report.py, or delete them if "
              "they are superseded -- a figure that exists and is never shown "
              "is how work goes missing.")
        return 1
    print("every cross-energy figure reaches the page")
    return 0


if __name__ == "__main__":
    sys.exit(main())
