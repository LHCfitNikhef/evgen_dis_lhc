#!/usr/bin/env python3
"""Catch figure labels that LaTeX will silently mangle.

WHY THIS EXISTS.  With `text.usetex` on (analysis/plotstyle.py), a bare `%` in
a label starts a LaTeX COMMENT and everything after it on that line is
DISCARDED.  The figure still renders, nothing errors, and the label is simply
wrong:

    r"$\\sigma_{\\rm charm}/\\sigma_{\\rm inclusive}$  [%]"   ->  "... ["
    "bands are PDF uncertainties at 68% CL"                 ->  "... at 68"

Both were live in results_nu/cmp_pdf_dependence_charmfrac_E.png on 2026-08-27
and were spotted only by looking at the picture.  That is the CONVENTIONS.md rule 2
failure mode applied to figures: the run completed, the output looked
plausible, and nothing anywhere returned non-zero.

`_` and `&` and `#` are the same hazard with a louder failure -- they usually
raise instead of truncating -- but they are checked here too so the diagnosis
arrives before the traceback does.

THE CONTRACT: every user-visible string reaching matplotlib goes through
plotstyle.tex(), which escapes the specials outside math mode and passes
$...$ segments through untouched.  This checks that contract statically, so a
new label cannot quietly skip it.

Usage: tools/check_figure_labels.py
Exit status 1 if any label bypasses tex() while containing a special.
"""
import ast
import glob
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# the calls that put text on a figure
TEXT_FUNCS = {"set_xlabel", "set_ylabel", "set_title", "suptitle", "text",
              "annotate"}
# specials that are fatal or silently destructive in LaTeX text mode
SPECIALS = "%_&#"
MATH = re.compile(r"\$[^$]*\$")


def literal_parts(node):
    """Every string literal reachable in this expression, flattened.

    Handles the two shapes these scripts use: adjacent-string concatenation
    (which the parser has already folded into one Constant) and explicit
    `a + b` joins across lines.
    """
    out = []
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        out.append(node.value)
    elif isinstance(node, ast.JoinedStr):
        for v in node.values:
            out.extend(literal_parts(v))
    elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        out.extend(literal_parts(node.left))
        out.extend(literal_parts(node.right))
    elif isinstance(node, ast.FormattedValue):
        pass                       # a runtime value; nothing static to check
    return out


def guarded(node):
    """True if this expression is wrapped in tex(), wholly or piecewise."""
    if isinstance(node, ast.Call):
        f = node.func
        return (isinstance(f, ast.Name) and f.id == "tex") or \
               (isinstance(f, ast.Attribute) and f.attr == "tex")
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        # "tex(a) + '\n' + tex(b)": every literal-bearing side must be guarded
        return all(guarded(s) or not literal_parts(s)
                   for s in (node.left, node.right))
    return False


def offending(text):
    """The specials in `text` that sit OUTSIDE math mode."""
    prose = MATH.sub("", text)
    return sorted({c for c in SPECIALS if c in prose})


def check(path):
    tree = ast.parse(open(path).read(), filename=path)
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else \
            (fn.id if isinstance(fn, ast.Name) else None)
        if name not in TEXT_FUNCS or not node.args:
            continue
        # ax.text(x, y, s) puts the string third; the rest put it first
        arg = node.args[2] if (name in ("text", "annotate")
                               and len(node.args) >= 3) else node.args[0]
        if guarded(arg):
            continue
        for lit in literal_parts(arg):
            chars = offending(lit)
            if chars:
                bad.append((node.lineno, name, "".join(chars), lit[:70]))
    return bad


def main():
    total = 0
    for path in sorted(glob.glob(f"{BASE}/analysis/*.py")):
        if os.path.basename(path) == "plotstyle.py":
            continue
        bad = check(path)
        if not bad:
            continue
        total += len(bad)
        print(f"\n{os.path.relpath(path, BASE)}:", file=sys.stderr)
        for lineno, fn, chars, lit in bad:
            print(f"  line {lineno}: {fn}() has bare {chars!r} outside math "
                  f"mode and does not go through tex()\n"
                  f"    {lit!r}", file=sys.stderr)
    if total:
        print(f"\n{total} label(s) will be mangled by LaTeX. A bare '%' "
              f"COMMENTS OUT the rest of\nthe string and the figure still "
              f"renders, so this does not show up as an error.\n"
              f"Wrap the string in plotstyle.tex().", file=sys.stderr)
        return 1
    print("figure labels: every text call with a LaTeX special goes through "
          "tex()")
    return 0


if __name__ == "__main__":
    sys.exit(main())
