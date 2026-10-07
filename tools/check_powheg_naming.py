#!/usr/bin/env python3
"""NEVER A BARE "POWHEG" IN READER-FACING TEXT (user rule, restated 2026-08-28).

Two related but DIFFERENT codes carry the name in this benchmark:

    POWHEG-RES   powheg-dis-main/DIS_v        the muon NC NLO entry
    POWHEG-V2    powheg-cmass/nu-DIS-master   the neutrino CC entry

They differ in ways that change results, not just in provenance: POWHEG-RES
hardcodes a diagonal CKM matrix and so has no charm decomposition on the
charged current, while POWHEG-V2 is the massive-charm code -- that is what
`cmass` means and why it is the neutrino entry -- though every production card
here runs it at qmass = 0 to meet the benchmark's massless-charm convention.
A sentence that says only "POWHEG" therefore does not identify a calculation,
in exactly the way "YADISM" without a scheme does not (CLAUDE.md rule 3).

WHAT IS NOT A VIOLATION.  "POWHEG" is also the name of a MATCHING METHOD, and
in that sense it is correct unqualified -- Herwig's POWHEG mode, "the POWHEG
hardest emission", "POWHEG-BOX".  Those senses are allowed below.  The rule is
about naming the GENERATOR ENTRY.

Checks the rendered report and the paper sources, which is what a reader sees;
the code's own comments and identifiers are a different audience, like the
capitals rule.

Usage: check_powheg_naming.py [page.html] ; reads paper/*.tex (and the page, if
given); exits non-zero on a violation.
"""
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# senses in which a bare POWHEG is correct: the matching method, the framework
# name, and the LaTeX macros that expand to a qualified name
ALLOWED = re.compile(
    # the matching METHOD, in every phrasing it appears in here.  "(POWHEG)"
    # and "(POWHEG, charm)" are the annotation on a Herwig row -- Herwig's
    # POWHEG mode -- and are the single most common correct bare use.
    r"\(POWHEG[,)]|"
    r"POWHEG(?:-BOX|\s+(?:method|matching|hardest|emission|mode|Hooks|hooks|"
    r"weights|rwl|LHE|box|implementation|family))|"
    r"(?:its|the|LO,|NLO,|through|independent)\s+POWHEG|"
    # a CITATION of the POWHEG method papers names the method, not a code
    r"POWHEG~?\\cite|"
    # a FILE NAME is a name, not a claim about which code: faser_powheg_rates.py
    r"faser_powheg_rates|powheg_v2_faser|powheg_v2/|"
    # Herwig's row of the generator table: "NLO ({\sc POWHEG})" is its
    # matching method, the LaTeX spelling of the "(POWHEG)" annotation above
    r"\(\{\\sc POWHEG\}\)|"
    # \powhegres and \powhegtwo are fine.  THE BARE \powheg MACRO IS ALSO
    # FINE IN THE PAPER (user ruling, 2026-09-21): "when I refer to \powheg
    # without specification, I mean RES for mu and V2 for nu" -- the paper's
    # generator section assigns each code to one current, so the current
    # identifies the code.  This is the paper's convention only: a bare
    # "POWHEG" typed as text, and anything on the report page, still fails.
    r"\\powheg(?:res|two(?:mc)?)?(?![a-zA-Z])|"
    r"(?:Herwig[^.]{0,40}|angular-ordered[^.]{0,20})POWHEG|"
    r"POWHEG[^ ]{0,3}(?:MEDIS|DIS_v)|"
    # THE FORWARD-CHARM POWHEG IS A THIRD CODE, and naming it is what this
    # rule wants rather than what it forbids: arXiv:2309.12793 generates
    # charm-hadron production with POWHEG's heavy-quark process showered by
    # Pythia 8.3, which is neither of the two DIS codes.  "POWHEG + Pythia 8.3"
    # is the flux authors' own name for it and identifies it unambiguously.
    r"POWHEG \+ Pythia|POWHEG charm|charm flux[^.]{0,30}POWHEG",
    re.I)
BARE = re.compile(r"POWHEG(?!-(?:V2|RES|BOX))", re.I)


def strip_html(t):
    t = re.sub(r"<script.*?</script>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<style.*?</style>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    return t


def scan(text, label):
    bad = []
    for n, line in enumerate(text.splitlines(), 1):
        # macro DEFINITIONS, and LaTeX COMMENTS, are not what a reader sees --
        # this file's own rule-explaining comments would otherwise trip it
        if "newcommand" in line or line.lstrip().startswith("%"):
            continue
        for m in BARE.finditer(line):
            seg = line[max(0, m.start() - 45):m.start() + 45]
            if ALLOWED.search(seg):
                continue
            bad.append((n, seg.strip()))
    return bad


def main():
    targets = []
    rep = sys.argv[1] if len(sys.argv) > 1 else None
    if rep and os.path.exists(rep):
        targets.append((rep, strip_html(open(rep, errors="replace").read())))
    pdir = f"{BASE}/paper"
    for fn in sorted(os.listdir(pdir)) if os.path.isdir(pdir) else []:
        if fn.endswith(".tex"):
            p = os.path.join(pdir, fn)
            targets.append((p, open(p, errors="replace").read()))
    if not targets:
        print("nothing to check (no report and no paper sources)")
        return 0
    total = 0
    for path, text in targets:
        bad = scan(text, path)
        if bad:
            print(f"\n{os.path.relpath(path, BASE)}: "
                  f"{len(bad)} unqualified mention(s)")
            for n, seg in bad[:12]:
                print(f"   line {n}: ...{seg}...")
            total += len(bad)
    if total:
        print(f"\n{total} bare \"POWHEG\" mention(s). Say POWHEG-V2 or "
              f"POWHEG-RES -- they are different codes and the difference "
              f"changes results.")
        return 1
    print("POWHEG naming: every mention says -V2 or -RES, or is the "
          "matching method")
    return 0


if __name__ == "__main__":
    sys.exit(main())
