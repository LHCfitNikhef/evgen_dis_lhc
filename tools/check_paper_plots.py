#!/usr/bin/env python3
"""Verify every paper plot's physics message against its own numbers.

USER INSTRUCTION, 2026-08-28: each paper plot carries the physics message it
is meant to convey, and the message must be CHECKED rather than believed --
"you should automatically check that my interpretation is correct, in the
light of what is shown in the plot and the project memory".

So each script in analysis/paper_plots/ declares CLAIMS: small, executable
statements that read the SAME result files the figure plots and return True
only if the message holds.  This runs them all.

WHY THIS MATTERS HERE SPECIFICALLY.  Numbers in this project have moved under
their own prose repeatedly -- a Sherpa caveat outlived its fix by two days, a
paper table kept pre-CKM numbers whose conclusion had inverted, and a report
blurb quoted a spread that had changed.  Every one of those was written by
someone reading a correct figure at the time.  A claim that re-derives itself
from the results on every build is the only version that cannot go stale.

A CLAIM THAT RESTATES A HAND-TYPED NUMBER CHECKS NOTHING.  If the message says
"CT18 lies 13% below", the claim must recompute that 13% from the result JSON,
not compare 13 to a literal 13 written beside it.

Usage: check_paper_plots.py [--dir paper_plots] [slug ...]
       exits non-zero if a claim fails.  --dir selects the figure directory
       under analysis/ (default paper_plots, the final region; the the earlier production
       directory analysis/paper_plots/ was deleted on 2026-09-19).
"""
import importlib.util
import os
import sys
import traceback

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PPDIR = os.path.join(BASE, "analysis", "paper_plots")
sys.path.insert(0, os.path.join(BASE, "analysis"))


def load(path):
    name = os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(f"paperplot_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    global PPDIR
    args = sys.argv[1:]
    if args[:1] == ["--dir"]:
        PPDIR = os.path.join(BASE, "analysis", args[1])
        args = args[2:]
    only = set(args)
    if not os.path.isdir(PPDIR):
        print(f"no {os.path.relpath(PPDIR, BASE)}/ directory")
        return 0
    scripts = sorted(f for f in os.listdir(PPDIR)
                     if f.endswith(".py") and not f.startswith("_"))
    if not scripts:
        print("no paper plots defined yet -- nothing to check.\n"
              f"  Add one as {os.path.relpath(PPDIR, BASE)}/<slug>.py.")
        return 0
    nfail = ntot = nskip = 0
    for fn in scripts:
        slug = os.path.splitext(fn)[0]
        if only and slug not in only:
            continue
        path = os.path.join(PPDIR, fn)
        try:
            mod = load(path)
        except Exception:                                     # noqa: BLE001
            print(f"\n{fn}: FAILED TO IMPORT")
            traceback.print_exc()
            nfail += 1
            continue
        claims = getattr(mod, "CLAIMS", None)
        msg = getattr(mod, "MESSAGE", "")
        print(f"\n{fn}")
        if not msg.strip():
            print("   !! no MESSAGE -- every paper plot must say what it "
                  "teaches")
            nfail += 1
        if not claims:
            print("   !! no CLAIMS -- a physics message that nothing checks "
                  "is exactly what this tool exists to prevent")
            nfail += 1
            continue
        # A claim marked "private" needs an input that is not part of the
        # public release (CONVENTIONS.md); without it the claim is SKIPPED,
        # never silently passed, and says so.
        avail = getattr(mod, "private_available", lambda: True)
        have_private = bool(avail())
        for c in claims:
            ntot += 1
            what = c.get("what", "(unnamed claim)")
            if c.get("private") and not have_private:
                print(f"   SKIP   {what}   [needs data not in the public release]")
                nskip += 1
                continue
            try:
                ok = bool(c["check"]())
                detail = c.get("detail")
                extra = f"   [{detail()}]" if callable(detail) else ""
            except Exception as exc:                          # noqa: BLE001
                print(f"   ERROR  {what}\n          {type(exc).__name__}: "
                      f"{exc}")
                nfail += 1
                continue
            print(f"   {'ok  ' if ok else 'FAIL'}   {what}{extra}")
            if not ok:
                nfail += 1
    print()
    if nfail:
        print(f"{nfail} problem(s) across {ntot} claim(s). The physics "
              f"message and the figure disagree -- fix one of them.")
        return 1
    skipped = f" ({nskip} skipped: their inputs are not public)" if nskip else ""
    print(f"all {ntot - nskip} claim(s) hold{skipped}: every physics message "
          f"matches its figure's own numbers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
