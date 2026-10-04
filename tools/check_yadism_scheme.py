#!/usr/bin/env python3
"""Refuse a bare "YADISM" in the report: the scheme must always be named.

WHY (user rule, 2026-08-27).  YADISM computes structure functions in whichever
flavour scheme it is asked for, and this benchmark uses TWO of them:

    ZM-VFNS   massless charm and bottom -- the benchmark's own convention,
              and what every inclusive YADISM curve on the page is
    FONLL     charm mass effects restored in the heavy-quark part -- what
              every CHARM reference on the page is, since 2026-08-27
    FFNS      n_f = 3 with the charm PDF set to zero -- POWHEG-V2mc's own
              scheme, drawn beside it on paper plots 5 and 6 (2026-09-19)

They differ by 7-13% on charm and by up to 0.9% on the inclusive rate, and the
difference has a sign that changes with energy and with current.  So "YADISM
NLO" on its own does not identify a calculation: a reader cannot tell which of
two curves differing by 13% is meant, and neither can an author six months
later.  The first survey of this found 294 bare mentions against 1 qualified.

WHAT COUNTS AS QUALIFIED.  A scheme word -- ZM-VFNS, ZM, FONLL -- within a
short window either side of the mention, so that
"YADISM CC NLO, FONLL" and "the FONLL YADISM reference" both pass while
"YADISM CC NLO" alone does not.  The window is deliberately generous: this is
a prompt to say which scheme, not a grammar checker.

WHERE IT LOOKS.  The RENDERED page, like tools/check_report_caps.py, because
one source string can render on twelve tabs and the count that matters is what
a reader meets.

Usage: tools/check_yadism_scheme.py [path/to/report.html]
Exit status 1 if any mention is unqualified.
"""
import collections
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = f"{BASE}/results/report.html"
sys.path.insert(0, os.path.join(BASE, "tools"))

# "either scheme" / "both schemes" / "same scheme" are explicit statements
# about the scheme rather than silence about it, and are as unambiguous as
# naming one.  "same scheme" earns its place in a RATIO context: a table whose
# rows each divide by their own reference cannot name a single scheme without
# being wrong for some row, and saying the two sides match is the true
# statement.
# FFNS joined on 2026-09-19: the fixed-flavour n_f = 3 curve (charm PDF set to
# zero) drawn beside POWHEG-V2mc on paper plots 5 and 6.
SCHEME = (r"(?:ZM-VFNS|ZM-VFN|ZM|FONLL|FFNS|fixed-flavour|either scheme|both schemes|"
          r"same(?: order and)? scheme|same-scheme)")
# how far either side of "YADISM" a scheme word may sit and still count
WINDOW = 60


def prose(html):
    """The text a reader actually sees -- same reduction as check_report_caps."""
    h = re.sub(r"data:image/[^\"']+", "", html)
    h = re.sub(r"<script.*?</script>|<style.*?</style>", " ", h, flags=re.S)
    h = re.sub(r"<[^>]+>", " ", h)
    for ent, ch in (("&mdash;", "-"), ("&ndash;", "-"), ("&nbsp;", " "),
                    ("&amp;", "&"), ("&plusmn;", "+-"), ("&times;", "x"),
                    ("&rsquo;", "'"), ("&ldquo;", '"'), ("&rdquo;", '"')):
        h = h.replace(ent, ch)
    return re.sub(r"\s+", " ", h)


def unqualified(text):
    """Every YADISM mention with no scheme word nearby."""
    out = []
    for m in re.finditer(r"YADISM", text):
        lo = max(0, m.start() - WINDOW)
        ctx = text[lo:m.end() + WINDOW]
        if re.search(SCHEME, ctx):
            continue
        out.append(text[m.start():m.end() + 55].strip())
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(path):
        sys.exit(f"no report at {path} -- build it with "
                 f"analysis/make_report.py first")
    text = prose(open(path, encoding="utf-8", errors="replace").read())
    total = len(re.findall(r"YADISM", text))
    bad = unqualified(text)
    if not bad:
        print(f"YADISM scheme: all {total} mention(s) name ZM-VFNS or FONLL")
        return 0
    print(f"*** {len(bad)} of {total} YADISM mention(s) do not name a scheme "
          f"***", file=sys.stderr)
    for t, n in collections.Counter(bad).most_common(25):
        print(f"  x{n:<4} {t[:74]}", file=sys.stderr)
    print("\nYADISM computes in whichever scheme it is asked for, and this "
          "benchmark uses two:\nZM-VFNS for the inclusive curves, FONLL for "
          "the charm references. They differ by\n7-13% on charm, so a bare "
          "'YADISM' does not identify a calculation.\nName the scheme "
          "wherever YADISM is mentioned (user rule, 2026-08-27).",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
