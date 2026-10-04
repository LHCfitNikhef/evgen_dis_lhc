#!/usr/bin/env python3
"""Splice the benchmark's DIS cut and scale into a generated MG5 process dir.

    install_hooks.py <PROC dir>

MG5 writes SubProcesses/dummy_fct.f with stub versions of `dummy_cuts` and
`user_dynamical_scale`; this replaces exactly those two routines with the ones
in dis_hooks.f and leaves the rest of the file as MG5 wrote it, so a change of
MG5 version keeps its own stubs for everything else instead of inheriting a
frozen copy from this repository.

THE CUTS COME FROM analysis/selection.py, not from a number typed here: the
fiducial region is defined once, and a sample generated to a different Q2 or y
than the analysis applies would be short by an amount nothing reports.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "analysis"))
import selection  # noqa: E402

SEL = selection.SELECTIONS["inclusive"]   # the fiducial region: Q2 and y only

# THE GENERATION REGION IS SLIGHTLY WIDER THAN THE FIDUCIAL ONE, ON PURPOSE.
#
# A sample generated EXACTLY on the fiducial boundary can only lose events
# across it -- there is nothing outside to come back the other way -- so any
# reshuffling of the kinematics between generation and analysis turns into a
# one-sided cross-section deficit.  The first MG5 run was made that way and
# lost 1.9% of its events, which reads as a 2% deficit and is an artefact of
# the setup rather than a statement about the generator.
#
# WHAT THE MARGIN DOES AND DOES NOT BUY, measured on this sample.  At PARTON
# level it makes the normalisation exact: the fiducial fraction is a measured
# number (91.3% muon, 97.1% neutrino) multiplying the generated cross section,
# and the result reproduces the analytic leading-order reference to 1.0000 and
# 0.9999.  It does NOT recover what SHOWERING moves, because that migration is
# one-way: 52 events of 3000 crossed out of the region and NONE crossed in.
#
# THE CAUSE IS THE HEAVY-QUARK MASSES, and it is worth stating because it is
# invisible in any cross-section: the benchmark's convention is massless charm
# and bottom, MG5 writes them massless in the Les Houches file, and Pythia
# puts every parton back on ITS OWN mass shell as it reads one -- 1.5 GeV for
# charm, 4.8 for bottom.  Forcing that mass on a massless-generated event
# reshuffles the whole hard process, so a charm- or bottom-initiated event
# comes back with a DIFFERENT reconstructed Q2: 9.6% of muon events move by
# more than 1%, by -14% on average in that tail.  At the Q2 = 4 edge, where
# dsigma/dQ2 is largest, that costs about 1.5% of the fiducial rate.
# Re-showering the same file with `4:m0 = 0` and `5:m0 = 0` removes the shift
# COMPLETELY -- zero crossings, mean shift -5e-5 GeV2 -- which is what
# identifies the cause; ISR, the beam remnants, the primordial kT and
# hadronisation were each switched off in turn and changed nothing.
#
# THE SAMPLE IS STILL SHOWERED WITH PYTHIA'S OWN MASSES, deliberately: every
# other showered sample in the benchmark is, so this one has to be for the
# comparison to be like-for-like.  It shows up as the showered MG5 row sitting
# 1.5% below the parton-level one, next to Pythia's own leading-order sample
# at 0.984 -- the same offset, and Pythia's massless-quark diagnostic
# (pythia8/dis_mu1TeV_me_massless.cmnd) recovers the same 1.3%.
Q2_MARGIN = 0.95        # generate from 0.95 x the fiducial Q2 floor
Y_MARGIN = 0.01         # and 0.01 outside each y edge

# A fixed-form Fortran routine: header line, then everything up to the "end"
# that closes it (at the start of a line, no continuation).
ROUTINES = ("dummy_cuts", "user_dynamical_scale")


def split_routines(text):
    """{name: source} for every routine in a fixed-form Fortran file."""
    head = re.compile(
        r"^ {6}\s*(?:(?:logical|double precision|real\*8|integer)\s+)?"
        r"(?:function|subroutine)\s+([A-Za-z_]\w*)", re.I | re.M)
    out, marks = {}, [(m.start(), m.group(1).lower()) for m in head.finditer(text)]
    for i, (pos, name) in enumerate(marks):
        stop = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        out[name] = (pos, stop)
    return out


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    proc = sys.argv[1]
    target = os.path.join(proc, "SubProcesses", "dummy_fct.f")
    if not os.path.exists(target):
        sys.exit(f"no {target} -- is {proc} an MG5 process directory?")

    hooks = open(os.path.join(HERE, "dis_hooks.f")).read()
    q2_gen = Q2_MARGIN * float(SEL.q2_min)
    y_lo = max(0.0, float(SEL.y_min) - Y_MARGIN)
    y_hi = min(1.0, float(SEL.y_max) + Y_MARGIN)
    hooks = (hooks.replace("@Q2MIN@", repr(q2_gen))
                  .replace("@YMIN@", repr(y_lo))
                  .replace("@YMAX@", repr(y_hi)))
    ours = split_routines(hooks)
    missing = [r for r in ROUTINES if r not in ours]
    if missing:
        sys.exit(f"dis_hooks.f does not define {missing}")

    text = open(target).read()
    if "evgen-benchmark" in text:
        print(f"  hooks already installed in {target}")
        return
    theirs = split_routines(text)
    missing = [r for r in ROUTINES if r not in theirs]
    if missing:
        sys.exit(f"{target} has no stub for {missing} -- MG5 layout changed")

    # replace back to front so the earlier offsets stay valid
    for name in sorted(ROUTINES, key=lambda n: theirs[n][0], reverse=True):
        a, b = theirs[name]
        src = hooks[ours[name][0]:ours[name][1]].rstrip() + "\n\n"
        text = text[:a] + src + text[b:]

    banner = ("C     evgen-benchmark: dummy_cuts and user_dynamical_scale below are\n"
              "C     INSTALLED by mg5/install_hooks.py from mg5/dis_hooks.f.  Edit that\n"
              "C     file, not this one -- this copy is regenerated and gitignored.\n")
    open(target, "w").write(banner + text)
    print(f"  installed DIS hooks in {target} "
          f"(generation: Q2 > {q2_gen:g}, {y_lo:g} < y < {y_hi:g}; "
          f"fiducial: Q2 > {SEL.q2_min:g}, "
          f"{SEL.y_min:g} < y < {SEL.y_max:g})")


if __name__ == "__main__":
    main()
