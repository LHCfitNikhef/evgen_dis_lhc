#!/usr/bin/env python3
"""Refuse a report whose JavaScript does not parse, or whose tabs and figures dangle.

WHY.  The page is one self-contained file whose figures and tabs are wired up
by a single inline script: the figures are filled in from an IMG map, and the
tab buttons are bound to their panels.  If that script has a syntax error,
NOTHING runs -- every figure is blank and no tab can be opened -- yet the HTML
is still written, its size looks normal and every prose checker passes, since
they read the text.  On 2026-09-30 a bulk text clean-up stripped every "()"
from the code, including make_report.py's script
("addEventListener('click',  => select(t))"), and the page shipped dead.

WHAT IT CHECKS, on the rendered page:
  * the inline script parses (`node --check`; SKIPPED with a note if node is
    not installed -- the structural checks below still run);
  * every <img data-img="iN"> has an entry in the IMG map;
  * every tab's aria-controls names an element that exists.

Usage: tools/check_report_js.py [path/to/report.html]
Exit status 1 on any failure.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = f"{BASE}/results/report.html"


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.isfile(path):
        print(f"FAIL: {path} does not exist")
        return 1
    page = open(path, encoding="utf-8").read()
    bad = []

    scripts = re.findall(r"<script>(.*?)</script>", page, re.S)
    if not scripts:
        bad.append("no inline <script> -- figures and tabs cannot work")
    node = shutil.which("node")
    if node and scripts:
        with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as fh:
            fh.write("\n".join(scripts))
            js = fh.name
        try:
            r = subprocess.run([node, "--check", js], capture_output=True, text=True)
            if r.returncode:
                bad.append("inline script does not parse:\n" + r.stderr.strip())
        finally:
            os.unlink(js)
    elif not node:
        print("note: node not installed, script syntax not checked")

    keys = set(re.findall(r'data-img="([^"]+)"', page))
    have = set(re.findall(r'^"(i\d+)":', "\n".join(scripts), re.M))
    for k in sorted(keys - have):
        bad.append(f"figure {k} has no entry in the IMG map")

    ids = set(re.findall(r'\bid="([^"]+)"', page))
    for a in re.findall(r'aria-controls="([^"]+)"', page):
        if a not in ids:
            bad.append(f"tab controls missing panel '{a}'")

    for b in bad:
        print("FAIL:", b)
    if not bad:
        print(f"ok: script parses, {len(keys)} figures and "
              f"{page.count('aria-controls=')} tabs resolve")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
