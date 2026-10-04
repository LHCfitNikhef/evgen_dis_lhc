#!/usr/bin/env python3
"""The TARGET nucleon, and how tungsten is built out of two of them.

ONE PLACE, for the same reason selection.py is one place.  The benchmark's
process is a lepton on a PROTON at rest, and that is unchanged: every
published cross-section here is a proton cross-section.  But the FASER rate
predictions are for a TUNGSTEN target, and tungsten is 74 protons and 110
neutrons.

>>> WHY THIS MODULE EXISTS (user, 2026-09-04). <<<  "For the sake of the
'Predictions for FASER' studies we have to generate cross-sections also for
neutron target, so that the W nucleus is correctly reproduced (for neutrino
DIS xsections on proton and neutron are quite different)."  They are: at
1 TeV, GENIE's total CC on a neutron is 1.80x its proton value for nu_mu and
0.64x for nubar_mu, because nu + d -> mu- + u runs on the valence d quark and
a neutron has two of them.  Counting all 184 nucleons as protons, which is
what faser_rates.py did between 2026-09-02 and today, therefore UNDERSTATES
the neutrino rate and OVERSTATES the antineutrino one, and the two do not
cancel because the flux is not charge-symmetric.

>>> THE NEUTRON IS THE ISOSPIN MIRROR OF THE PROTON, NOT A NEW FIT. <<<
u <-> d and ubar <-> dbar in NNPDF4.0 member 0, built once by
tools/make_neutron_pdf.py as the LHAPDF set NNPDF40_nnlo_as_01180_n (id
339900) under data/pdfs/lhapdf/.  So a neutron cross-section is the same
calculation with a different PDF set, generator-independent, and the isospin
relation is exact up to QED, which this benchmark switches off.

>>> A NEUTRON RESULT CARRIES "_n" IN ITS FILENAME. <<<  Same rule as the
energy tag, the selection suffix and the Q2 floor, and for the same reason:
a neutron cross-section is a perfectly plausible number that would silently
replace the proton one otherwise.  suffix() is "" for the proton, so nothing
that already exists moves.

Pick the target from the environment, like the selection:

    BENCH_TARGET=n analysis/yadism_cc_calc.py nlo
"""
import os

# THE NEUTRON SET LIVES IN THIS REPOSITORY, so the search path has to include
# it before anything calls lhapdf.  Done here rather than in each driver
# because every consumer imports this module and none of them should have to
# remember: an LHAPDF set that is not found does not fall back, it aborts, but
# a caller that forgot would abort in the middle of a ladder rather than at
# its start.  Prepended, never replacing, so the system sets stay reachable.
_PDFDIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "data", "pdfs", "lhapdf")
if _PDFDIR not in os.environ.get("LHAPDF_DATA_PATH", "").split(":"):
    os.environ["LHAPDF_DATA_PATH"] = (
        _PDFDIR + ":" + os.environ.get("LHAPDF_DATA_PATH", ""))

# Tungsten, the FASERnu target.  A = Z + N is the mass number the rate
# convolution uses; the elemental mixture is 74/110 to better than the
# accuracy of anything downstream (natural W is Z = 74, <A> = 183.84).
Z_W, N_W = 74, 110
A_W = Z_W + N_W
A_W_GMOL = 183.84

PDFSET_P = "NNPDF40_nnlo_as_01180"
PDFSET_N = "NNPDF40_nnlo_as_01180_n"      # tools/make_neutron_pdf.py, id 339900
# The free-nucleon average (74p + 110n)/184, tools/make_isoscalar_pdf.py, id
# 341000.  Structure functions are LINEAR in the PDF, so a calculation with
# this set on a proton target IS the tungsten cross-section per nucleon -- the
# analytic references use it directly.  THE GENERATORS DO NOT: they run p and
# n separately (user, 2026-09-13), because on a 2212 beam this set would leave
# a proton remnant under a PDF that is not a proton's.
PDFSET_W = "NNPDF40_nnlo_as_01180_W184free"

# "W" is a pseudo-target: per-nucleon tungsten.  For the analytic side it
# selects PDFSET_W; for a generator it can only be the COMBINATION of a p and
# an n result (per_nucleon below), never a run of its own.
TARGETS = ("p", "n", "W")
DEFAULT = "p"


def current(name=None):
    """The target: the argument, then $BENCH_TARGET, then the proton."""
    t = name or os.environ.get("BENCH_TARGET") or DEFAULT
    if t not in TARGETS:
        raise SystemExit(f"BENCH_TARGET must be one of {TARGETS}, got {t!r}")
    return t


def pdfset(name=None):
    """The LHAPDF set for one target."""
    return {"p": PDFSET_P, "n": PDFSET_N, "W": PDFSET_W}[current(name)]


def suffix(name=None):
    """Legacy filename suffix: "" for the proton, "_n" / "_W" otherwise.

    The proton stays untagged so nothing published before moves.
    """
    t = current(name)
    return "" if t == "p" else f"_{t}"


def tag(name=None):
    """Filename suffix, ALWAYS explicit: "_p", "_n" or "_W".

    The final-region results (selection q4w3) exist on three targets side by
    side, so an untagged proton file would be the one name a reader could
    mistake for the tungsten result the figure wants.
    """
    return f"_{current(name)}"


# ------------------------------------------------ target-mass corrections --
# >>> BENCH_TMC IS THE ONE TMC KNOB (user, 2026-09-13). <<<  0 = off, which
# is what every result before was computed with and stays the default, so
# nothing published moves; 1 / 2 / 3 are yadism's own modes (esf/tmc.py):
#   1  "APFEL"        F2/FL with the h2 integral, the g2 term dropped
#   2  "approximate"  the integrals replaced by one evaluation at xi
#                     (Schienbein et al. eq. 4)
#   3  "exact"        Georgi-Politzer with both h2 and g2 (F3 with h3)
# 3 IS THE ONE THE BENCHMARK USES ("_tmc"), measured on the 1 TeV NLO grid
# before choosing: 2 is off by a factor 2.4 on FL at Q2 = 4, x = 0.5 (TMC
# shift +679% against +284% exact), 1 is within 1% of 3 on FL and exact on
# F3, and 3 costs 2.7x a TMC-off run on the massless structure functions and
# 1.5x on the massive ones -- affordable, so there is no reason to take less.
#
# It lives here because a target-mass correction is a property of the target
# (yadism's MP card entry is the nucleon mass), and every card builder,
# filename and cache key that must know it already imports this module.  It is
# READ AT CALL TIME, like current(), so a driver that loops over settings in
# one process cannot be served a stale value.
TMC_MODES = (0, 1, 2, 3)
TMC_CHOSEN = 3


def tmc():
    """The yadism TMC mode: $BENCH_TMC, else 0 (off)."""
    v = os.environ.get("BENCH_TMC") or "0"
    try:
        m = int(v)
    except ValueError:
        m = -1
    if m not in TMC_MODES:
        raise SystemExit(f"BENCH_TMC must be one of {TMC_MODES}, got {v!r}")
    return m


def tmc_suffix():
    """Filename suffix: "" off, "_tmc" for the chosen mode, "_tmc<N>" else.

    A TMC-on reference is a plausible number a few per cent from the TMC-off
    one, so it can never share its name -- the same rule as the target, the
    selection and the energy tag.
    """
    m = tmc()
    return "" if m == 0 else ("_tmc" if m == TMC_CHOSEN else f"_tmc{m}")


def provenance(pdfset_name=None):
    """{pdf_set, target, tmc} to stamp into a result, or None at the defaults.

    NONE AT THE DEFAULTS (proton, TMC off) so that every pre-result, when
    re-derived, stays byte-identical; anything else carries the stamp, so a
    tungsten or TMC-on JSON says what it is once loaded, not only in its name.
    """
    t, m = current(), tmc()
    if t == DEFAULT and m == 0 and pdfset_name in (None, PDFSET_P):
        return None
    return {"pdf_set": pdfset_name or pdfset(), "target": t, "tmc": m}


def combine(sigma_p, sigma_n):
    """The tungsten cross-section per NUCLEUS from the two nucleon ones."""
    return Z_W * sigma_p + N_W * sigma_n


def per_nucleon(sigma_p, sigma_n):
    """The tungsten cross-section per NUCLEON -- what a column density in
    nucleons/cm^2 multiplies."""
    return (Z_W * sigma_p + N_W * sigma_n) / A_W
