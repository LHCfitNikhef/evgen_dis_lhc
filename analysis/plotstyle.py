#!/usr/bin/env python3
"""ONE plotting style for every figure in the benchmark.

Adopted 2026-08-27 at the user's request, from the style they use for
publication figures (faser/Sec6-Fig1-NSBIg-bin-vs-unbin-N7-ENWM2.py).  It is
the standard PDF-fitting / HEP-paper look, and the report and the paper should
not be typeset differently from each other.

WHAT IT SETS, and why each one matters:

  * LaTeX text rendering with Helvetica, so figure text is typeset by the same
    engine as paper/, rather than matplotlib's own approximation of it.
  * TICKS POINTING INWARD, ON ALL FOUR SIDES, major length 7 / minor 4.  This
    is the single change that most separates a publication figure from a
    default matplotlib one.
  * A very light grid (lw=0.1) that guides the eye without competing with the
    data.
  * Explicit font sizes: axis labels 16, ylabel 15, ticks 12, legend 12-14.
  * subplots_adjust rather than tight_layout -- see fused_gridspec() below.

WHY THERE IS NO SILENT FALLBACK.  usetex needs a working LaTeX with helvet and
type1cm, plus dvipng.  If it is missing, the honest options are to fail or to
render differently -- and rendering differently is the worse one: the figures
would silently change appearance between this machine and the Nikhef cluster,
and nobody would notice until two versions of the same plot appeared side by
side.  So a missing LaTeX is a hard error naming what to install, and
BENCH_NO_USETEX=1 is the deliberate, explicit opt-out for a machine that
cannot have it.

Usage, at the top of a plotting script and BEFORE importing pyplot:

    import plotstyle
    plotstyle.apply()
    import matplotlib.pyplot as plt
"""
import os
import re
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib import rc                                  # noqa: E402

USETEX = os.environ.get("BENCH_NO_USETEX") != "1"

# Font sizes, named once so every script agrees and a change is one edit.
FS_XLABEL = 16
FS_YLABEL = 15
FS_TICKS = 12
FS_LEGEND = 12
FS_TITLE = 16
FS_SUPTITLE = 18
FS_PANEL_TITLE = 14

_applied = False


def _usetex_available():
    """(ok, reason).  Probed once, cheaply, before matplotlib tries it."""
    for exe in ("latex", "dvipng"):
        if shutil.which(exe) is None:
            return False, f"{exe} is not on PATH"
    try:
        out = subprocess.run(["kpsewhich", "helvet.sty", "type1cm.sty"],
                             capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError) as e:
        return False, f"kpsewhich failed ({e.__class__.__name__})"
    found = [ln for ln in out.stdout.split("\n") if ln.strip()]
    if len(found) < 2:
        return False, "helvet.sty or type1cm.sty is missing from the TeX tree"
    return True, ""


def apply():
    """Install the style.  Call once, before pyplot is imported."""
    global _applied
    if _applied:
        return
    if USETEX:
        ok, why = _usetex_available()
        if not ok:
            sys.exit(
                f"plotstyle: LaTeX rendering is not available -- {why}.\n"
                f"  Every figure in this benchmark is typeset with LaTeX so "
                f"the report and paper/ match.\n"
                f"  Install a TeX distribution with helvet and type1cm (and "
                f"dvipng), or set\n"
                f"  BENCH_NO_USETEX=1 to render with matplotlib's own mathtext "
                f"instead -- which\n"
                f"  produces VISIBLY DIFFERENT figures, so do not mix the two "
                f"in one report.")
        rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
        rc("text", usetex=True)
    else:
        rc("font", **{"family": "sans-serif"})
        rc("text", usetex=False)

    matplotlib.rcParams.update({
        "axes.labelsize": FS_YLABEL,
        "axes.titlesize": FS_TITLE,
        "legend.fontsize": FS_LEGEND,
        "xtick.labelsize": FS_TICKS,
        "ytick.labelsize": FS_TICKS,
        # inward ticks on all four sides -- the defining feature of the style
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 7,
        "xtick.minor.size": 4,
        "ytick.major.size": 7,
        "ytick.minor.size": 4,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.linewidth": 1.0,
        "grid.linewidth": 0.1,
        "grid.color": "black",
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "figure.dpi": 150,
        "savefig.dpi": 150,
    })
    _applied = True


def ticks(ax, labelbottom=True):
    """The tick treatment, applied per axis.

    rcParams covers most of it, but direction/length have to be re-applied to
    an axis whose ticks were configured after creation, and `labelbottom` is
    how the top panel of a fused pair hands its labels to the one below.
    """
    ax.tick_params(which="both", direction="in", labelsize=FS_TICKS,
                   right=True, top=True, labelbottom=labelbottom)
    ax.tick_params(which="major", length=7)
    ax.tick_params(which="minor", length=4)
    ax.grid(True, lw=0.1)


# --------------------------------------------------------------- label TeX
# Existing labels are a mix: some are plain prose ("Sherpa 3.0.5 MC@NLO"),
# some already carry math ("charm seen as $D^{\\pm}$").  Blanket-wrapping
# either one breaks the other, so the math segments are preserved verbatim and
# only the prose between them is wrapped in \rm.
_MATH = re.compile(r"(\$[^$]*\$)")
# LaTeX specials that must be escaped in TEXT mode
_ESCAPE = {"_": r"\_", "%": r"\%", "&": r"\&", "#": r"\#",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
           "^": r"\textasciicircum{}"}


def tex(s):
    """Plain-or-mixed label -> a usetex-safe string in the house style.

    PROSE IS LEFT IN TEXT MODE, not wrapped in ${\rm ...}$.  Both render as
    upright roman and look identical, but inside math mode LaTeX retypesets
    punctuation: "POWHEG-V2" comes out as "POWHEG - V2" with a math MINUS and
    its spacing, and "MC@NLO*" grows a centred asterisk operator.  Text mode
    gives the same face with the punctuation intact.

    Math segments written by the caller ($D^{\pm}$, $E_{\nu}$) are passed
    through untouched, so a label may mix the two freely.

    Returns the input unchanged when usetex is off, so a BENCH_NO_USETEX run
    still produces readable (if differently typeset) figures.
    """
    if not USETEX or s is None or s == "":
        return s
    out = []
    for part in _MATH.split(s):
        if part.startswith("$") and part.endswith("$"):
            out.append(part)                     # already math: leave alone
        elif part:
            out.append("".join(_ESCAPE.get(c, c) for c in part))
    return "".join(out)


def fused_gridspec(fig, nrows_top=2, height_ratios=(2.0, 1.0), **kw):
    """A main panel with its ratio panel FUSED beneath it: no vertical gap.

    hspace=0.0 and a shared x-axis, so the ratio panel reads as part of the
    same plot rather than a second figure stacked under it.

    NB the caller must NOT then call tight_layout(): it overrides hspace and
    reopens the gap.  Use fig.subplots_adjust(...) with explicit margins --
    the reference script carries the same warning for the same reason.
    """
    from matplotlib import gridspec
    kw.setdefault("hspace", 0.0)
    return gridspec.GridSpec(nrows_top, 1, figure=fig,
                             height_ratios=list(height_ratios), **kw)


# THE SCALE BAND SHARES ITS CURVE'S LEGEND ENTRY (user, 2026-10-04: "whenever
# possible, integrate the 'NLO MHOU' legend into the corresponding legend for
# the curve of the tool used to evaluate the associated central value").
MHOU_SUFFIX = " w. MHOU"   # "w. MHOU", not "± MHOU" (user, 2026-10-04)


def merge_bands(handles, labels, merges):
    """Fold each band's legend entry into the curve(s) it belongs to.

    `merges` maps a band's label to the labels of the curves whose central
    value it surrounds, e.g. {tex("NLO MHOU"): [tex("YADISM (ZM)")]}.  Each
    such curve's key becomes (band, curve) -- the line drawn over the band,
    matplotlib's default for a tuple handle -- and its label gains
    MHOU_SUFFIX; the band's own entry is dropped.  A band whose curves are
    all absent keeps its entry, so nothing silently disappears.
    """
    h, lab = list(handles), list(labels)
    for band, curves in merges.items():
        if band not in lab:
            continue
        i = lab.index(band)
        bh, hit = h[i], False
        for c in curves:
            if c in lab:
                j = lab.index(c)
                h[j] = (bh, h[j])
                # tex() works segment by segment, so tex(a) + tex(b) is
                # tex(a + b): the suffix can be appended to a tex'd label
                lab[j] = c + tex(MHOU_SUFFIX)
                hit = True
        if hit:
            k = lab.index(band)
            del h[k], lab[k]
    return h, lab

