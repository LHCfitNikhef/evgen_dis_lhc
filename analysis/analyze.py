#!/usr/bin/env python3
"""Common analysis for the muon-DIS generator benchmark.

Parses HepMC3 Ascii files from any generator, computes frame-invariant DIS
observables, and writes binned histograms (dsigma/dO in pb per unit observable)
plus the fiducial cross-section to a JSON file in results/.

Fiducial region: Q2 > 4 GeV2, 0.2 < y < 0.9 (all invariants built from the
beam records of each event, so the generation frame is irrelevant).

Usage:
  analyze.py pythia   # reads pythia8/job_*/events.hepmc + *_xsec.json
  analyze.py sherpa   # reads the Sherpa Runs/MuonDIS_NLO job_*/evtfull
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- config.sh is the one source of every path
import beams  # noqa: E402
import selection  # noqa: E402
import sample_layout  # noqa: E402  -- where the paper-plots samples live
from runmeta import stamp  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHERPA_RUNS = paths.SHERPA_RUNS
POWHEG_RES = paths.POWHEG_RES
POWHEG_V2 = paths.POWHEG_V2
# from the integration logs of the validated runs (pb)
SHERPA_SIGMA_PB = 35019.0        # MC@NLO: BVI 34782.3 + RS 236.7
# massless-charm convention (2026-08-13): c/b in the 93 container
SHERPA_LO_SIGMA_PB = 38097.5     # LO Comix, +- 34.6 (integ_massless.log)
SHERPA_LO_ME_SIGMA_PB = 38032.8  # same, ME-level run dir, +- 34.5
# POWHEG-DIS NLO, generation region Q2 > 2.25 (Qmin 1.5), full y: combined
# stage-2 integrator value from DIS_v/parallel-mu1TeV-wide (LHE <init>),
# +- 284. The generation cut sits below the analysis Q2 > 4 so radiation
# migrates both ways across the fiducial boundary (a Q2 > 4 generation cut
# biased sigma_fid low by ~1%).
POWHEG_SIGMA_PB = 216553.0
# POWHEG-V2 run on the MUON NC (cards/POWHEG-V2/mu1TeV-NC).  A third
# independent NLO matching on the muon side.  From prod-mu1TeV's pwg-stat.dat,
# "total (btilde+remnants)".
#
# RE-INTEGRATED 2026-08-27 at 8x (ncall2 4e5 x itmx2 4, seed 137, grids
# rebuilt) as the decisive test of TODO item 5b:
#
#     2026-08-26, seed 12    215685.94 +- 3047.47   (1.41%)
#     2026-08-27, seed 137   214675.62 +-  939.63   (0.44%)   <- this one
#     POWHEG-RES, 10 seeds   216553    +-  284      (0.13%)
#
# The two V2 integrations are INDEPENDENT (different seed, rebuilt grids) and
# agree at 0.33 sigma, so the first was imprecise rather than wrong.  Against
# POWHEG-RES on the SAME wide generation region the gap is now -0.867%, which
# at the sharpened error is 1.9 sigma -- it was "0.40%, inside one sigma"
# only because the old error was 3x larger.
#
# NOTE WHAT THIS DOES NOT FIX.  The number here is a NORMALISATION, and
# sigma_fid scales with it through the ratio estimator, so the fiducial row
# moves by exactly this -0.468% -- from 1.0370 to about 1.0321 against YADISM.
# The 3.7% excess is therefore NOT an integrator error.  It survives into the
# fiducial region, where V2 and POWHEG-RES sit 4.8% apart while their
# inclusive integrals agree to 0.9%, i.e. it lives in the ACCEPTANCE into
# Q2 > 4, 0.2 < y < 0.9 rather than in the total rate.
POWHEG_V2_MU_SIGMA_PB = 214675.62
POWHEG_V2_MU_SIGMA_ERR_PB = 939.63

M_P = 0.938272
M_MU = 0.10566

# --------------------------------------------------------- selection + energy
# The fiducial region and the beam energy are the two things this analysis is
# run once per.  Both come from a module that owns them (selection.py,
# beams.py) rather than being constants here, because the energy scan and the
# FASER selection are the same parameterisation problem.
#
#   BENCH_SELECTION=inclusive BENCH_ENERGY=400 analysis/analyze.py pythia
#
# Q2_MIN/Y_MIN/Y_MAX are kept as names because analyze_nu.py and
# me_nlo_samples.py import them; they are now just the active selection's.
SELECTION = selection.get()
TRACK_E_MIN = selection.TRACK_E_MIN
Q2_MIN, Y_MIN, Y_MAX = SELECTION.q2_min, SELECTION.y_min, SELECTION.y_max
# None for every benchmark selection; a number only for the foreign
# arXiv:2506.13889 region used by the FASER-tab muon reproduction.
W2_MIN = SELECTION.w2_min

ENERGY = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
BEAMS = beams.Beams("mu", ENERGY)
ETAG = BEAMS.tag


def at_energy(name):
    """Per-energy directory/result name -- see beams.at_energy()."""
    return beams.at_energy(name, ENERGY)


def job_files(parent, base, fname):
    """Files at <parent>/<base>_<N>/<fname>, N an INTEGER job number.

    NOT a plain glob of "<base>_*".  At the anchor `base` is "job", and the
    energy scan's own directories are job_400GeV_N and job_4TeV_N -- which
    "job_*" matches perfectly happily, silently folding three beam energies
    into one result.  It did exactly that: the 1 TeV Pythia result went from
    2,000,000 parsed events to 6,000,000 and sigma_fid rose 10.6%, with no
    error anywhere.  Only a before/after diff of the anchor results
    caught it.

    So the job number is matched exactly, and a directory whose suffix is not
    a bare integer is not a job of this sample.
    """
    pat = re.compile(rf"^{re.escape(base)}_\d+$")
    out = []
    for d in sorted(glob.glob(f"{parent}/{base}_*")):
        if os.path.isdir(d) and pat.match(os.path.basename(d)):
            f = os.path.join(d, fname)
            if os.path.exists(f):
                out.append(f)
    return sorted(out)


# ----------------------------------------------- integrator-file parsers
# Each of the four normalisation readers below needs its file parsed, and so
# does check_anchor_literals().  Split out so there is one parser per format
# rather than one per (format x caller) -- the same reason the DIS kinematics
# were pulled out of the two analysers.
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def _sherpa_last_xs(log):
    """Last '<process> : X pb +- ( Y pb' line -- a single-process integration.

    Sherpa colours its output, so ANSI has to come off before matching.  A raw
    grep for this pattern finds nothing at all, which is misleading enough to
    be worth stating.
    """
    pat = re.compile(r":\s*([0-9.eE+-]+)\s*pb\s*\+-\s*\(\s*([0-9.eE+-]+)\s*pb")
    hit = None
    with open(log, errors="replace") as f:
        for line in f:
            m = pat.search(_ANSI.sub("", line))
            if m:
                hit = m
    return None if hit is None else (float(hit.group(1)), float(hit.group(2)))


def _sherpa_bvi_rs(log):
    """BVI + RS, which is what an MC@NLO total is.  Both or nothing."""
    pat = re.compile(r"\((BVI|RS)\)\s*:\s*([0-9.eE+-]+)\s*pb\s*\+-\s*"
                     r"\(\s*([0-9.eE+-]+)\s*pb")
    parts = {}
    with open(log, errors="replace") as f:
        for line in f:
            m = pat.search(_ANSI.sub("", line))
            if m:
                parts[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    if {"BVI", "RS"} - set(parts):
        return None
    return parts


def _powheg_res_grand_total(f):
    """pos - |neg| from a stage-3 statistics file."""
    pos = neg = None
    for line in open(f, errors="replace"):
        m = re.search(r"grand total pos\.\s+weights:\s+([0-9.eE+-]+)", line)
        if m:
            pos = float(m.group(1))
        m = re.search(r"grand total \|neg\.\|\s+weights:\s+([0-9.eE+-]+)", line)
        if m:
            neg = float(m.group(1))
    return None if pos is None or neg is None else (pos, neg)


def _powheg_v2_total(f):
    """'total (btilde+remnants) cross section in pb' from pwg-stat.dat."""
    for line in open(f, errors="replace"):
        m = re.search(r"total \(btilde\+remnants\) cross section in pb\s+"
                      r"([0-9.eEdD+-]+)", line)
        if m:
            return float(m.group(1).replace("D", "E").replace("d", "e"))
    return None


# ------------------------------------------- the anchor literals, and a check
# EVERY reader below returns a hardcoded literal at the anchor energy and only
# parses a file at the others.  That guard is deliberate and right: the anchor
# is what every published number was produced from, and a re-integration must
# not be able to move it by itself.
#
# WHAT IT IS NOT ALLOWED TO BE IS SILENT.  Until 2026-08-27 nothing ever
# compared a literal against the file it came from, so re-integrating the
# anchor left the literal pointing at the superseded number and the new sample
# was normalised against it without a word.  That happened on the POWHEG-V2
# muon row the same day (TODO item 5b): the literal had to be edited by hand,
# and had it not been, the analysis would have printed a perfectly plausible
# cross-section built on the old integrator.
#
# So the literal still WINS -- the file never silently overrides it -- but the
# two are compared and a disagreement is fatal, with both numbers named.
#
# THE VERIFICATION FILE IS NOT ALWAYS THE FILE THE NON-ANCHOR PATH READS, and
# Sherpa LO is why.  MuonDIS_LO/ holds four integration logs; integ.log there
# is the OLD MASSIVE run and reads 34228.1 pb against the massless
# convention's 38097.5, an 11% difference.  Other energies have only one log,
# so their integ.log is the right file -- but the anchor's is
# integ_massless.log, and a check that read integ.log would "discover" an 11%
# discrepancy that is really a wrong-file error.
ANCHOR_LITERAL_TOL = 2e-3


def _anchor_sources():
    """(label, literal, path, parse) for every anchor literal, or None."""
    return [
        ("Sherpa LO (muon NC)", SHERPA_LO_SIGMA_PB,
         f"{SHERPA_RUNS}/MuonDIS_LO/integ_massless.log",
         lambda f: (_sherpa_last_xs(f) or (None,))[0]),
        ("Sherpa MC@NLO (muon NC)", SHERPA_SIGMA_PB,
         f"{SHERPA_RUNS}/MuonDIS_NLO/integ.log",
         lambda f: (lambda p: None if p is None
                    else sum(v for v, _ in p.values()))(_sherpa_bvi_rs(f))),
        ("POWHEG-RES (muon NC)", POWHEG_SIGMA_PB,
         f"{POWHEG_RES}/parallel-mu1TeV-wide/pwg-0001-st3-stat.dat",
         lambda f: (lambda g: None if g is None else g[0] - g[1])(
             _powheg_res_grand_total(f))),
        ("POWHEG-V2 (muon NC cross-variant)", POWHEG_V2_MU_SIGMA_PB,
         f"{POWHEG_V2}/prod-mu1TeV/pwg-stat.dat",
         _powheg_v2_total),
    ]


def check_anchor_literals(verbose=False):
    """Verify every anchor literal still matches the file it was taken from.

    A missing file is NOT an error -- samples are pruned and install trees get
    rebuilt, and the literal is precisely the record that survives that.  A
    file that is present and DISAGREES is fatal.
    """
    bad, checked = [], 0
    for label, lit, path, parse in _anchor_sources():
        if not os.path.exists(path):
            if verbose:
                print(f"  [anchor] {label}: {lit:.6g} pb "
                      f"(file gone, literal is the record)")
            continue
        try:
            got = parse(path)
        except (OSError, ValueError):
            got = None
        if got is None:
            if verbose:
                print(f"  [anchor] {label}: {lit:.6g} pb (unparseable file)")
            continue
        checked += 1
        rel = got/lit - 1.0
        if abs(rel) > ANCHOR_LITERAL_TOL:
            bad.append(f"{label}\n     literal {lit:.6g} pb (published)\n"
                       f"     file    {got:.6g} pb  ({100*rel:+.3f}%)\n"
                       f"     {path}")
        elif verbose:
            print(f"  [anchor] {label}: {lit:.6g} pb, file agrees "
                  f"({100*rel:+.4f}%)")
    if bad:
        raise SystemExit(
            "ANCHOR LITERAL DISAGREES WITH ITS INTEGRATION FILE\n\n  "
            + "\n\n  ".join(bad)
            + "\n\n  The literal is what every published 1 TeV number was "
              "normalised with, so it wins\n  until someone changes it "
              "deliberately.  If the re-integration is intended, edit the\n"
              "  literal in analysis/analyze.py and re-run every result that "
              "depends on it.\n  If it is not, the run directory has been "
              "overwritten and the file is the wrong one.")
    return checked


def sherpa_lo_sigma_pb():
    """Sherpa LO integrator cross-section for the configured energy, in pb.

    THE ANCHOR KEEPS ITS LITERAL.  SHERPA_LO_SIGMA_PB above comes from
    integ_massless.log -- the massless-convention integration -- and NOT from
    the MuonDIS_LO/integ.log sitting beside it, which is the older massive run
    and reads 34228.1 pb.  Parsing the wrong one of those two would quietly
    renormalise every published 1 TeV Sherpa number by 11%.

    Other energies have only one integration, done with the current card, so
    their integ.log is the right file and is read here.  A hardcoded 1 TeV
    literal was what the delivery-closure gate caught: at 400 GeV it compared
    32.31 nb of delivered events against the 1 TeV 38.0975 nb and reported a
    15% "failure" that was purely the wrong reference.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return SHERPA_LO_SIGMA_PB
    log = f"{SHERPA_RUNS}/{at_energy('MuonDIS_LO')}/integ.log"
    if not os.path.exists(log):
        raise SystemExit(f"no Sherpa integration log at {log} -- integrate first")
    hit = _sherpa_last_xs(log)
    if hit is None:
        raise SystemExit(f"could not find the integrated cross-section in {log}")
    print(f"  [sherpa] integrator {hit[0]/1e3:.4f} nb "
          f"+- {hit[1]/1e3:.4f} (from {os.path.basename(log)})")
    return hit[0]


def sherpa_nlo_sigma_pb():
    """Sherpa MC@NLO integrator cross-section for the configured energy, in pb.

    NOT one number but the SUM OF TWO PROCESSES.  Sherpa integrates MC@NLO as
    BVI (Born + virtual + integrated subtraction) and RS (real minus
    subtraction) separately, and the physical total is their sum.  Taking the
    last matching line -- which is what the LO reader does, correctly, for a
    single-process integration -- would return the RS piece alone.  RS is
    -1043.6 pb at 400 GeV, i.e. NEGATIVE, so that mistake would not merely be
    wrong, it would report a negative cross-section.

    The anchor keeps its literal SHERPA_SIGMA_PB (= BVI 34782.3 + RS 236.7),
    so re-running 1 TeV cannot move a published number.

    This branch had the 1 TeV LITERAL that the LO branch was already fixed for:
    sigma_sample_pb was SHERPA_SIGMA_PB at every energy, so a 400 GeV sample
    would have been normalised with the 1 TeV cross-section.  The
    delivery-closure gate would have caught it -- it caught exactly this on the
    LO side -- but it should never have been reachable.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return SHERPA_SIGMA_PB
    log = f"{SHERPA_RUNS}/{at_energy('MuonDIS_NLO')}/integ.log"
    if not os.path.exists(log):
        raise SystemExit(f"no Sherpa NLO integration log at {log} -- "
                         f"integrate first")
    parts = _sherpa_bvi_rs(log)
    missing = {"BVI", "RS"} if parts is None else set()
    if missing:
        raise SystemExit(f"{log}: no {'/'.join(sorted(missing))} cross-section "
                         f"found -- MC@NLO needs BOTH pieces, and their sum is "
                         f"the answer")
    tot = sum(v for v, _ in parts.values())
    print(f"  [sherpa NLO] integrator {tot/1e3:.4f} nb "
          f"= BVI {parts['BVI'][0]/1e3:.4f} + RS {parts['RS'][0]/1e3:.4f} "
          f"(from {os.path.basename(log)})")
    return tot


# --- Sherpa SHOWER-MODEL arms ----------------------------------------------
# The NLO half of the shower study (sherpa/run_sherpa_shower.sh).  Each arm is
# the SAME card and the SAME matrix element with SHOWER_GENERATOR overridden on
# the command line, run in its own directory <base>_<shower>.
#
# EACH ARM READS ITS OWN integ.log, never the baseline's.  MC@NLO's subtraction
# is built from the shower's splitting kernels, so a Dire sample normalised to
# a CSS integration would be internally inconsistent while looking entirely
# reasonable -- the failure mode CONVENTIONS.md rule 2 is about.  There is no anchor
# literal here for the same reason: the arms are new, so nothing is published
# that a literal would have to protect.
def sherpa_arm_sigma_pb(rundir, nlo, with_err=False):
    """Integrated cross-section in pb for one shower arm, from its own log.

    `with_err` returns (sigma, error).  The neutrino side needs it: its Sherpa
    generation cuts EQUAL the fiducial cuts, so an inclusive result is pinned
    straight to the integrator and carries the integrator's own error rather
    than a weight-sum error.
    """
    log = f"{SHERPA_RUNS}/{rundir}/integ.log"
    if not os.path.exists(log):
        raise SystemExit(f"no integration log at {log} -- run "
                         f"sherpa/run_sherpa_shower.sh for this arm first")
    if nlo:
        parts = _sherpa_bvi_rs(log)
        if parts is None or {"BVI", "RS"} - set(parts):
            raise SystemExit(f"{log}: MC@NLO needs BOTH the BVI and RS pieces "
                             f"and their sum is the answer")
        tot = sum(v for v, _ in parts.values())
        err = sum(e ** 2 for _, e in parts.values()) ** 0.5
        print(f"  [{rundir}] integrator {tot/1e3:.4f} nb = "
              f"BVI {parts['BVI'][0]/1e3:.4f} + RS {parts['RS'][0]/1e3:.4f}")
        return (tot, err) if with_err else tot
    hit = _sherpa_last_xs(log)
    if hit is None:
        raise SystemExit(f"could not find the integrated cross-section in {log}")
    print(f"  [{rundir}] integrator {hit[0]/1e3:.4f} nb +- {hit[1]/1e3:.4f}")
    return (hit[0], hit[1]) if with_err else hit[0]


# key -> (run directory, is NLO, label)
SHERPA_SHOWER_ARMS = {
    "sherpa_dire":    ("MuonDIS_NLO_Dire", True,  "Sherpa 3.0.5 MC@NLO, Dire"),
    "sherpa_lo_dire": ("MuonDIS_LO_Dire",  False, "Sherpa 3.0.5 LO, Dire"),
}


def powheg_res_sigma_pb():
    """POWHEG-RES integrated cross-section for the configured energy, in pb.

    Read from the run's own stage-3 statistics file as

        grand total pos. weights  -  grand total |neg.| weights

    which at 1 TeV gives 229448.73 - 12895.33 = 216553.4 pb, reproducing the
    POWHEG_SIGMA_PB literal above to the digit.  The anchor keeps that literal
    so nothing published can move; other energies are read.

    Seed 0001's file is used, matching how the anchor number was obtained --
    each seed integrates the full phase space, so it is a complete estimate
    rather than a partial sum.

    Without this the delivery-closure gate compares a 4 TeV sample against the
    1 TeV integrator and reports +51%, which is exactly what it did.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return POWHEG_SIGMA_PB
    tag = BEAMS.tag
    f = f"{POWHEG_RES}/parallel-mu{tag}-wide/pwg-0001-st3-stat.dat"
    if not os.path.exists(f):
        raise SystemExit(f"no POWHEG statistics at {f} -- run "
                         f"powheg/run_powheg_res.sh {ENERGY:g} first")
    got = _powheg_res_grand_total(f)
    if got is None:
        raise SystemExit(f"could not read the grand totals from {f}")
    pos, neg = got
    sig = pos - neg
    print(f"  [powheg] integrator {sig/1e3:.4f} nb "
          f"(pos {pos/1e3:.4f} - |neg| {neg/1e3:.4f}, from {os.path.basename(f)})")
    return sig


def powheg_v2_mu_sigma_pb():
    """POWHEG-V2 muon NC integrated cross-section for this energy, in pb.

    Read from that run's own pwg-stat.dat, "total (btilde+remnants)" -- the
    single-run POWHEG-V2 shape, not POWHEG-RES's per-seed stage-3 files.

    The anchor keeps its literal, as every other entry here does, so a
    re-integration cannot move a published number without someone editing this
    file.  Unlike the older literals it currently equals what the stat file
    says, because this sample is new and has no history the file has lost.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return POWHEG_V2_MU_SIGMA_PB
    tag = BEAMS.tag
    f = f"{POWHEG_V2}/prod-mu{tag}/pwg-stat.dat"
    if not os.path.exists(f):
        raise SystemExit(f"no POWHEG-V2 statistics at {f} -- run "
                         f"powheg/run_powheg_v2.sh --current mu {ENERGY:g} "
                         f"first")
    v = _powheg_v2_total(f)
    if v is None:
        raise SystemExit(f"could not read the total cross section from {f}")
    print(f"  [powheg-v2-mu] integrator {v/1e3:.4f} nb (from pwg-stat.dat)")
    return v


# Where each POWHEG entry's LHE files live, for lhe_events_offered().  Stated
# once, next to the cross-section readers that look in the same directories --
# a wrong directory here would silently give the closure gate the wrong
# denominator, which is the very failure this sidecar exists to expose.
POWHEG_LHE_DIR = {
    "powheg":    lambda: f"{POWHEG_RES}/parallel-mu{BEAMS.tag}-wide",
    "powheg_v2": lambda: f"{POWHEG_V2}/prod-mu{BEAMS.tag}",
    # paper plots: one integration directory per nucleon, from sample_layout
    "v2_powheg": lambda: sample_layout.powheg_rundir(
        "mu", sample_layout.require()[1], ENERGY),
}

# The QED study's shower arms (powheg/qed/, run_powheg_shower.sh --qed ARM).
# They re-shower the SAME LHE files as `powheg`, so they share its integrator
# cross-section, its LHE directory and therefore its closure gate -- which is
# the point: matrix element, matching, PDF and scale are held fixed and only
# the radiation changes.
#
# THE `off` ARM IS NOT LISTED because it is not a separate sample: re-showering
# with QED off reproduces powheg/job_N BYTE-IDENTICALLY (logs/qed_closure_off
# .log), so the baseline arm IS the published `powheg` result and generating a
# copy of it would be 3.5 GB of duplicate.
QED_ARMS = {
    "fsr":    "QED FSR (lepton line)",
    "fsrisr": "QED FSR+ISR (lepton line)",
    "full":   "QED FSR+ISR + quark line",
}


def qed_arm_of(pbase):
    """The arm name for a `powheg_qed_<arm>` key, or None.

    An UNKNOWN arm raises rather than returning None: `powheg_qed_typo` would
    otherwise fall through the dispatch to the generic weighted branch and be
    normalised by a rule that does not apply to it.
    """
    if not pbase.startswith("powheg_qed_"):
        return None
    arm = pbase[len("powheg_qed_"):]
    if arm not in QED_ARMS:
        sys.exit(f"no such QED arm {arm!r} -- powheg/qed/ defines "
                 f"{', '.join(sorted(QED_ARMS))} (and 'off', which is the "
                 f"published `powheg` sample itself)")
    return arm


def fmt_sigma(pb):
    """A fiducial cross-section in whichever unit keeps it readable.

    The benchmark spans nb (muon NC inclusive, ~38 nb) to well below a pb
    (the dimuon tier), so no single fixed unit and precision serves them all.
    """
    if pb >= 1e3:
        return f"{pb/1e3:.4g} nb"
    if pb >= 1.0:
        return f"{pb:.4g} pb"
    return f"{pb*1e3:.4g} fb"


# A selection can be so tight that the sample cannot measure it.  That is a
# legitimate outcome, but it must not be presented as a cross-section.
STAT_REL_ERR_MAX = 0.5
# A cross-section built from a handful of events is not a measurement whatever
# its formal error says.  GENIE's muon dimuon points (3 and 4 selected events)
# slipped through on the relative-error test alone, because the binomial error
# on a tiny selected fraction lands right at the 50% boundary.
MIN_SELECTED_EVENTS = 10


def flag_low_stats(out, gen):
    """Mark a result the sample is too small to support, and say why.

    Two ways this shows up, and the first is not obviously a failure at all:

    * A NEGATIVE sigma_fid.  POWHEG carries negative-weight events, as any NLO
      subtraction does, and they cancel against positive ones only in the
      limit of decent statistics.  POWHEG's 400 GeV dimuon point selects FOUR
      events and sums to -1568 fb.  Nothing errored, and "-1568 fb" would have
      gone into a table looking like a measurement.
    * A relative statistical error above STAT_REL_ERR_MAX, i.e. a number whose
      own error bar makes it consistent with anything.

    Either way the result is still written -- "this sample cannot measure this
    selection" is information, and deleting it would leave a hole that looks
    like an unrun job -- but it is stamped so the report can render it as such
    rather than plotting it.
    """
    sig, err = out.get("sigma_fid_pb"), out.get("sigma_fid_err_pb")
    n_sel = out.get("n_fiducial") or 0
    why = None
    if n_sel == 0:
        # NOT a statistics failure and NOT negative-weight cancellation: no
        # event passed at all.  That is its own outcome, already announced by
        # the EMPTY SELECTION path, and it has a real physical reading (GENIE's
        # EM DIS has no charm channel, so the dimuon tier is empty for it by
        # construction).  Attributing it to weight cancellation, as the first
        # version of this guard did, invents a mechanism that is not there.
        pass
    elif sig is not None and sig < 0:
        why = (f"sigma_fid = {sig:.4g} pb is NEGATIVE from {n_sel} selected "
               f"events -- negative-weight cancellation with too few events")
    elif sig is not None and sig == 0:
        why = (f"sigma_fid is exactly zero despite {n_sel} selected events -- "
               f"weights cancelled exactly, which needs looking at")
    elif 0 < n_sel < MIN_SELECTED_EVENTS:
        why = (f"only {n_sel} selected events (< {MIN_SELECTED_EVENTS}) -- too "
               f"few for a cross-section whatever the formal error says")
    elif sig and err and abs(err / sig) > STAT_REL_ERR_MAX:
        why = (f"relative statistical error {abs(err/sig):.0%} from {n_sel} "
               f"selected events exceeds {STAT_REL_ERR_MAX:.0%}")
    if why:
        out["stat_insufficient"] = True
        out["stat_insufficient_why"] = why
        print(f"  !! {gen}: SAMPLE CANNOT MEASURE THIS SELECTION -- {why}.")
        print(f"     Result written and stamped stat_insufficient; the report "
              f"shows it as such rather than quoting it.")
    return out


def sel_suffix():
    """Result-name suffix for the active selection.

    EMPTY for the inclusive selection, so every published result keeps its
    original name; a `_<name>` suffix for the FASER tiers.  Without this a
    tier run would overwrite the inclusive result under the same filename --
    the identical trap that energy tagging had to fix twice, once when
    job_400GeV_* was matched by job_*, and once when analyze_nu wrote every
    energy to the anchor name.
    """
    return "" if SELECTION.name == selection.DEFAULT else f"_{SELECTION.name}"


def result_path(gen):
    """results/histos_<gen>[_<selection>][_<energy tag>][_q2min<N>].json.

    A key ("v2_<gen>") is written as histos_<gen>_q4w3_<p|n>[_<tag>].json
    instead (analysis/sample_layout.py): the region and the NUCLEON are part of
    the name, and sample_layout.require() refuses a run that did not set both.
    """
    if gen.startswith("v2_"):
        _sel, t = sample_layout.require()
        return (f"{BASE}/results/histos_"
                f"{at_energy(sample_layout.result_key(gen[3:], t, _sel.name))}.json")
    return (f"{BASE}/results/histos_{at_energy(gen + sel_suffix())}"
            f"{selection.q2_suffix()}.json")

# charged hadrons that are stable under the ctau > 10 mm convention
CHARGED_HADRONS = {211, 321, 2212, 3112, 3222, 3312, 3334}

# The neutrinos, which carry energy out of the event invisibly and are
# therefore excluded from the hadronic energy sum.  Written as a set rather
# than a chain of comparisons because it is read once per final-state
# particle of every event in the benchmark.
NEUTRINO_PIDS = {12, 14, 16}

# DeltaPhi sentinel: events with no charged hadron passing the N_ch selection
# have no minimum |dphi| at all. They are stored as this value, which sits
# outside every bin edge, so np.histogram drops them (see dphi_summary()).
NO_DPHI = -1.0

# ------------------------------------------------------------------ D mesons
# ctau(D+-) = 312 um and ctau(D0) = 123 um are both far below the 10 mm
# stability convention of this benchmark, so a D meson is NEVER a status-1
# final-state particle in these samples: it always decays. They are therefore
# counted at the point of PRODUCTION, from the decayed (status-2) records --
# every generator here writes them with status 2 -- and de-duplicated by
# dedup_dmesons() below. A filter on final-state particles returns zero.
D_CH_PID = 411        # D+ / D-
D_0_PID = 421         # D0 / D0bar
D_S_PID = 431         # D_s+ / D_s-: NOT folded into the charged-D observable
                      # (tracked separately, reported in "means")
D_PIDS_SIGNED = {D_CH_PID, -D_CH_PID, D_0_PID, -D_0_PID, D_S_PID, -D_S_PID}


# ------------------------------------------------------- final-state charm
# The user's definition of charm production (2026-08-21): an event counts if,
# AFTER shower and hadronisation, it contains at least one charm quark in the
# final state.  Confined charm is observable only as charm HADRONS, so that is
# what is tested for -- every charm quark leaving the shower ends up inside
# one.  D+-, D0, Ds, all their excited states, charmed baryons (Lambda_c,
# Xi_c, Omega_c, ...) and charmonium all count; so does a charm diquark, which
# only ever appears on the way to a charmed baryon.
#
# WHY THIS IS NOT THE HARD-PROCESS TAG.  Three things move between the two:
#   * shower g -> c cbar, absent at hard-process level;
#   * the beam remnant, which carries the partner c of an initial-state cbar
#     (the proton has zero net charm), so the cbar-initiated CC channel
#     nu cbar -> mu- sbar DOES have final-state charm even though its hard
#     process has none outgoing;
#   * charm that the hard process makes but the hadronisation loses -- nothing
#     does that, so this direction is empty.
# The tag is therefore a superset of the hard-process "charm out" tag and is
# closer to, but still not identical with, YADISM's V_cd/V_cs definition.
#
# The test is on the PDG quark-content digits rather than a hand-written list:
# generators write charm states this benchmark never anticipated (Herwig in
# particular emits high excitations), and a missing code would silently
# UNDERCOUNT.  Nuclear codes (10LZZZAAAI, e.g. the 1000010010 proton target)
# have to be excluded first -- their digits are not quark content.
def is_heavy_hadron(pdg, nq):
    """True if this PDG code is a hadron (or diquark) containing flavour nq.

    nq = 4 for charm, 5 for bottom.  Only the quark-content digits are read;
    the units digit is 2J+1, not a flavour (Delta++ = 2224 is not charm).
    """
    a = abs(int(pdg))
    if a < 100 or a >= 1000000000:
        return False        # quarks/leptons/bosons/strings/clusters; nuclei
    return nq in ((a // 10) % 10, (a // 100) % 10, (a // 1000) % 10)


def is_charm_hadron(pdg):
    """True if this PDG code is a hadron (or diquark) containing charm."""
    return is_heavy_hadron(pdg, 4)


# ------------------------------------------------------- delivery closure gate
# A weighted sample's sigma_fid is pinned to the INTEGRATOR's cross-section
# while its shapes come from the delivered events.  That combination is
# self-consistent only if the events actually reproduce the integrator.  When
# they do not -- because the generator silently discarded a kinematically
# biased subset -- sigma_fid still looks perfect by construction while every
# distribution is wrong.  That is exactly how Sherpa's RESPECT_MASSIVE_FLAG
# default hid for months: 8.9% of muon-NC events (23% of neutrino-CC ones)
# were thrown away, the charm channel was gutted 8x, and the only visible
# symptom was an event-level cross-section 25 sigma below the integrator.
#
# So: never trust the pinned number silently.  Always compute the delivered
# cross-section, compare, print, and persist the comparison.
CLOSURE_WARN = 0.005    # 0.5%: investigate
CLOSURE_FAIL = 0.02     # 2%: refuse to write a result


# ------------------------------------------------- declared closure deficits
# A sample may deliver measurably less than its integrator for a reason that is
# UNDERSTOOD and VERIFIED not to bias sigma_fid.  The gate must still refuse
# everything else, so such a case is DECLARED HERE, per sample, with its size,
# its cause and the evidence -- never waved through.
#
# WHY NOT ALLOW_BAD_CLOSURE=1.  That env var suppresses the gate for EVERY
# sample in the run, which is the wrong granularity twice over: it hides
# unrelated problems in the same pass, and it leaves no record of which sample
# was excused or why.  A declaration is per sample, is version controlled, and
# is printed on every run.
#
# WHY EACH ENTRY CARRIES A CEILING.  `max_dev` is the largest deviation the
# declaration covers, set just above what was measured.  A deficit that GROWS
# past it fails the gate again -- which is the point: the known cause is
# excused, a worsening of it is not.
DECLARED_DEFICITS = {
    ("powheg_nu", "400GeV"): {
        "max_dev": 0.025,          # measured 0.0216 on 2026-08-27
        "cause": "Pythia cannot put ~2.2% of POWHEG-V2's LHE final states on "
                 "shell ('ProcessContainer::constructProcess: setting mass "
                 "failed', 21527 of 1000000), the mirror of POWHEG's own "
                 "'lhefinitemasses: Cannot reshuffle'. Both codes struggle "
                 "with the same near-threshold configurations, and the rate "
                 "grows as the beam energy falls: 1.60% at 4 TeV, 1.88% at "
                 "1 TeV, 2.22% at 400 GeV.",
        "evidence": "sigma_fid is unaffected: on 20x the statistics "
                    "(48891 -> 977809 events) it moves -0.09%, i.e. "
                    "-0.16 sigma (1.2393 -> 1.2382 pb). The ratio estimator "
                    "cancels a proportional loss.",
    },
    # THE POWHEG-LHE SAMPLE SHOWERED BY HERWIG (2026-09-08).  Herwig discards
    # events whose shower or remnant it cannot reconstruct -- 3.6% of the
    # muon sample by count -- and each discard consumes one more Les Houches
    # event, so the delivered cross-section falls short of POWHEG's
    # integrator.  This is NOT a defect of this arm: it is the same discard
    # herwig7/herwig_xsec.py is written around, and the benchmark's standing
    # answer to it is the one used here, normalising to the integrated
    # cross-section so that the total is restored.
    #
    # >>> WHAT THAT REPAIRS AND WHAT IT DOES NOT. <<<  The fiducial rate is
    # fine: 35.136 +- 0.185 nb against the Pythia arm's 34.785 +- 0.040 on
    # THE SAME Les Houches events, +1.0% and 1.9 sigma, so the ratio
    # estimator absorbs the loss to within the statistics.  The SHAPES are
    # not: the discards are flavour dependent (d 0.1%, ubar 3.5%, s 7.4%,
    # cbar 17.4% -- measured in NEUTRINO_NOTES.md and quoted in
    # herwig_xsec.py) and they concentrate at high x, so Herwig's charm and
    # its highest-x bins are depleted by more than the average.  Any figure
    # using this arm's charm or high-x shape says so.
    # >>> PAPER PLOTS V2, POWHEG-V2 + Pythia 8: NO DECLARATION, BECAUSE THE
    # LOSS IS CURED (2026-09-14). <<<  Entries for ("v2_powheg_nu", "400GeV")
    # and "700GeV" stood here from 2026-09-13: Pythia dropped 2.2% / 2.1% of
    # the LHE events ("setting mass failed"), 99.8% of them charm, so the
    # inclusive q4w3 bias was +0.04..+0.07% and charm-tagged sigma -9..-10%,
    # declared inclusive_only (powheg/production/lost_events.py).  The samples were
    # re-showered with LesHouches:matchInOut = off (powheg/powheg_nu_v2.cmnd,
    # powheg/production/reshower_nu_matchinout.sh): on the 1 TeV proton batch 1 the
    # loss fell from 2.02% to 0.04%, the charm-tagged bias from -8.21% to
    # -0.00% and the inclusive one from +0.050% to +0.024%, and what is still
    # lost sits mostly outside the region (median W2 = 7).  analyze_nu.py
    # refuses a  POWHEG-V2 job whose shower did not run with the switch off.
    ("powheg_hw", "1TeV"): {
        "max_dev": 0.03,           # measured 0.0225 on 2026-09-08
        "cause": "Herwig discards 3.6% of these events when its shower or "
                 "remnant reconstruction fails, and each discard consumes "
                 "another Les Houches event; the delivered cross-section is "
                 "correspondingly short of POWHEG's integrator.",
        "evidence": "sigma_fid on the SAME LHE events is 35.136 +- 0.185 nb "
                    "against the Pythia arm's 34.785 +- 0.040, +1.0% and "
                    "1.9 sigma, so the ratio estimator absorbs the loss in "
                    "the rate. Shapes are NOT repaired: the discards are "
                    "flavour dependent and concentrate at high x (see "
                    "herwig7/herwig_xsec.py).",
    },
    # The charged-current half of the same arm.  Herwig discards MORE here
    # than on the muon side -- 5.2% against 3.6% -- and the closure is
    # correspondingly further out, but the fiducial rate closes BETTER:
    # 3.0744 +- 0.0102 pb against the Pythia arm's 3.0647 +- 0.0039 on the
    # same Les Houches events, +0.32% and 0.9 sigma.
    ("powheg_nu_hw", "1TeV"): {
        "max_dev": 0.06,           # measured 0.0502 on 2026-09-08
        "cause": "Herwig discards 5.2% of these events when its shower or "
                 "remnant reconstruction fails, and each discard consumes "
                 "another Les Houches event.  Same mechanism as the muon "
                 "arm above, larger here.",
        "evidence": "sigma_fid on the SAME LHE events is 3.0744 +- 0.0102 pb "
                    "against the Pythia arm's 3.0647 +- 0.0039, +0.32% and "
                    "0.9 sigma. Shapes are NOT repaired: mean N(D+-) is "
                    "0.0464 against 0.0566, an 18% charm depletion that is "
                    "the flavour-dependent discard of "
                    "herwig7/herwig_xsec.py, not hadronisation.",
    },
    # THE SIDIS LADDER BELOW 300 GeV, POWHEG-V2 ONLY (2026-09-21,
    # beams.SIDIS_ENERGIES_LOW).  Two understood deviations, each measured
    # and each negligible in the one number these points feed -- the
    # flux-folded region-only SIDIS yield, where E < 35 GeV is 0.8% and
    # 35-75 GeV is 3.5% of the nu_mu rate.
    ("v2_powheg_nu", "20GeV"): {
        "max_dev": 0.035,          # measured -0.0079 (p), -0.0297 (n)
        "cause": "Pythia cannot shower 2.1% (p) / 3.4% (n) of the LHE events "
                 "near threshold (BeamRemnants 'no momentum left', "
                 "ministring fragmentation failures): at 20 GeV sqrt(s) is "
                 "6.2 GeV and a quarter of the offered events have W < 3 GeV. "
                 "matchInOut = off is on; this is what is left.",
        "evidence": "powheg/production/lost_events.py: 80-82% of the lost events have "
                    "W2 < 9, i.e. lie outside the region, so the "
                    "delivered-event fiducial fraction is biased HIGH, "
                    "+1.40% +- 0.02% (p) and +2.42% +- 0.03% (n); lost "
                    "heavy flavour 78% (p) / 28% (n). In the flux-folded "
                    "yield that is +0.02%.",
    },
    ("v2_powheg_nu", "50GeV"): {
        "max_dev": 0.05,           # measured +0.0416 (p, 2 batches), -0.0017 (n)
        "cause": "POWHEG-V2's INTEGRATOR is low on the proton at 50 GeV: "
                 "0.18496 +- 0.00067 pb, against event-weight means of "
                 "0.1890, 0.1898, 0.1892, 0.1901 pb (+- 0.0004-0.0005) in "
                 "four spike-free batches, i.e. 2.5% and ~6 of its own "
                 "sigma -- the integrator-error question (decision C2) at "
                 "one more point. Batch 1 also holds one event of weight "
                 "988.8 against a mean of 0.19 (the unweighting-excess "
                 "correction where the bound is poor), which is why that "
                 "point runs to 7 batches: the spike is 0.7% of the weight.",
        "evidence": "The neutron at 50 GeV closes (-0.17%), as do 100 and "
                    "200 GeV (|dev| <= 0.6%). The pinned normalisation "
                    "carries the 2.5% into sigma_fid at 50 GeV proton; in "
                    "the flux-folded yield that is ~0.04%.",
    },
    # >>> PAPER PLOT 9, THE HERWIG ARM: NO DECLARATION, THE LOSS IS CURED
    # (2026-09-19). <<<  Entries ("v2_powheg_hw", "1TeV") at 3% and
    # ("v2_powheg_nu_hw", "1TeV") at 6% stood here from 2026-09-15: Herwig
    # dropped 3.7% / 3.2% (mu p / n) and 5.2% / 2.6% (nu) of the LHE events
    # when it could not put the remnant on shell, for closures of -2.30% /
    # -1.73% and -4.94% / -2.29%, and the drops were charm/sea/high x.
    # LesHouchesHandler:MaxEventErrorRetries 100 (patches/thepeg-2.3.0-event-
    # error-retry.diff, herwig7/LHE-{mu,nu}.in) now re-showers the SAME LHE
    # event instead; regenerated arm, q4w3: closures +0.44% / +0.04% (mu p/n)
    # and -0.11% / -0.39% (nu p/n), survival 99.4%.  herwig_xsec.read_lhe_job
    # refuses a arm job run without the switch.
}


def declared_deficit(gen):
    """The declaration covering this sample, or None.

    Keyed on the BASE generator and the beam tag, so a charm variant of a
    declared sample inherits the declaration -- it is the same events --
    UNLESS the declaration says `inclusive_only`: then the loss is known to
    be concentrated in the charm channel, and excusing a charm result with
    evidence gathered on the inclusive rate would be exactly the hidden bias
    this gate exists to stop.
    """
    base = gen.split("_charm")[0]
    d = DECLARED_DEFICITS.get((base, BEAMS.tag))
    if d and d.get("inclusive_only") and "_charm" in gen:
        return None
    return d


def input_manifest(files, sig_pb=None, counts=None, tag=""):
    """Record exactly which files produced a result, and flag mixed vintages.

    The production scripts mkdir -p predictable job_* directories and the
    analyses glob every match, so re-running a subset silently merges new
    output with old.  Nothing in the result recorded which files it came from.
    This writes that provenance into the JSON, and warns when the inputs look
    like they came from more than one production: a spread of modification
    times wider than a day is the cheap, reliable tripwire.
    """
    entries, mtimes = [], []
    for i, fn in enumerate(files):
        try:
            st = os.stat(fn)
        except OSError:
            continue
        mtimes.append(st.st_mtime)
        e = {"path": os.path.relpath(fn, BASE), "bytes": st.st_size,
             "mtime": st.st_mtime}
        if sig_pb is not None and i < len(sig_pb):
            e["sigma_gen_pb"] = sig_pb[i]
        if counts is not None and i < len(counts):
            e["n_accepted"] = counts[i]
        entries.append(e)

    man = {"n_files": len(entries), "files": entries}
    if mtimes:
        span_h = (max(mtimes) - min(mtimes)) / 3600.0
        man["mtime_span_hours"] = span_h
        if span_h > 24.0:
            print(f"WARNING{tag}: input files span {span_h/24:.1f} DAYS of "
                  f"modification time -- this result may mix jobs from "
                  f"different productions.  Check the job directories.")
    return man


def pooled_sigma_pb(sig_pb, counts, tag=""):
    """Combine per-job generated cross-sections into one number.

    Each job contributes n_j events that each carry weight sigma_j / n_j, so
    the pooled cross-section is the ACCEPTED-COUNT-WEIGHTED mean, not the plain
    mean: sum_j(n_j sigma_j) / sum_j(n_j).  The two coincide only for jobs of
    equal size, and the production drivers permit aborts and partial reruns
    (Pythia in particular has an abort budget), so equal sizes cannot be
    assumed.  Also warns when the jobs disagree more than their spread should
    allow, which is the signature of mixed-vintage job directories.
    """
    n_tot = sum(counts)
    if not n_tot:
        return 0.0
    sigma = sum(n * s for n, s in zip(counts, sig_pb)) / n_tot
    if len(sig_pb) > 1:
        spread = (max(sig_pb) - min(sig_pb)) / sigma
        if spread > 0.05:
            print(f"WARNING{tag}: per-job sigma_gen spans {100*spread:.1f}% "
                  f"({min(sig_pb)/1e3:.4f} to {max(sig_pb)/1e3:.4f} nb) -- "
                  f"are these jobs from the same production?")
        if len(set(counts)) > 1:
            print(f"NOTE{tag}: unequal job sizes {sorted(set(counts))}; "
                  f"using the count-weighted mean")
    return sigma


def lhe_events_offered_for(rundir, files):
    """As lhe_events_offered, but only the seeds these job directories used.

    An arm that showers a SUBSET of the seeds must be closed against that
    subset.  The seed is the trailing number of the job directory --
    qedfsr_job_7 showered pwgevents-0007.lhe -- so it is recovered from the
    inputs rather than passed in, and a directory whose name does not end in a
    number is skipped loudly by returning None, which makes the closure
    "NOT CHECKED" instead of wrong.
    """
    seeds = []
    for f in files:
        base = os.path.basename(os.path.dirname(f))
        n = base.rsplit("_", 1)[-1]
        if not n.isdigit():
            return None
        seeds.append(int(n))
    total, found = 0, False
    for i in sorted(set(seeds)):
        side = os.path.join(rundir, "pwgevents-%04d.lhe.nevents" % i)
        if not os.path.exists(side):
            return None
        with open(side) as fh:
            total += int(fh.read().split()[0])
        found = True
    return total if found else None


def lhe_events_offered(*rundirs):
    """How many LHE events the shower was OFFERED, from the .nevents sidecars.

    powheg/strip_lhe_nan.py writes one beside every LHE it passes over, and
    both shower drivers run it, so the number is recorded at the moment it is
    still knowable.  Returns None if no sidecar is present -- an older sample,
    which then falls back to the delivered count with a printed note rather
    than silently using the wrong denominator.

    A run directory's ADDED BATCHES count too.  powheg/add_v2_batches.sh
    raises POWHEG-V2 statistics by generating into "<rundir>-b<N>" siblings
    with the same cached grids and a different iseed, because raising numevts
    in place would regenerate the anchor LHE and break the lhe_index join.
    Those batches shower into the same pooled sample, so their events are
    OFFERED to the shower just as the parent's are, and leaving them out makes
    the denominator too small by exactly the factor the statistics grew.
    That is not hypothetical: it is what this gate caught on 2026-09-03,
    reporting +554% when 980156 delivered events were divided by the parent
    directory's 150000 -- the same shape as the +194% the note on
    POWHEG_NU_LHE_DIR records from the previous time the sample grew.
    """
    total, found = 0, False
    for d in rundirs:
        for dd in [d] + sorted(glob.glob(f"{d}-b*")):
            for f in sorted(glob.glob(f"{dd}/*.nevents")):
                with open(f) as fh:
                    total += int(fh.read().strip())
                found = True
    return total if found else None


def delivered_sigma_pb(run_stats, n_events, n_read=None):
    """The cross-section the DELIVERED events actually represent, in pb.

    Sherpa: sigma = sum(Weight) / sum(NTrials).  An accepted event may have
    taken several trials (mean ~4 here), so the mean weight overestimates the
    cross-section by exactly that factor -- it is not a bias, it is the wrong
    normalisation.  Generators writing a bare weight and no trial count
    (POWHEG here) fall back to the mean weight.

    !!! FOR AN LHE-DRIVEN SAMPLE THE DENOMINATOR IS THE NUMBER OF EVENTS THE
    !!! SHOWER WAS OFFERED, NOT THE NUMBER THAT CAME OUT.  (2026-08-26.)
    Pythia answers a `partonLevel failed; try again` on LHE input by
    DISCARDING that event and reading the next one.  Nothing counts those:
    they are not in the HepMC and they are not in main_powheg's abort tally.
    Measured on a 5000-event probe: POWHEG-RES offered 5000 and delivered
    5000, POWHEG-V2 muon NC offered 5000 and delivered 4861 + 11 aborts, with
    105 `partonLevel ... try again`.

    Dividing by what came OUT makes the loss cancel out of its own check --
    the numerator and the denominator both shrink -- which is why this went
    unnoticed.  Measured at 1 TeV:

                              sum(w)/N_out   sum(w)/N_offered   dropped
        POWHEG-RES mu (pub)      1.01082          1.00702           0
        POWHEG-RES nu            1.00166          1.00130           0
        POWHEG-V2  mu            1.02732          1.00070        1187
        POWHEG-V2  nu (pub)      1.00080          0.98199         929

    The V2 muon sample looked 2.7% biased and is not: its dropped events come
    in near-cancelling +/- pairs, so sum(w) survives while N falls 2.6%.  The
    published V2 neutrino sample looked fine and is genuinely 1.8% short of
    the integrator, because its dropped events took their weight with them.
    Both readings were wrong, in opposite directions, for the same reason.

    No sigma_fid moves: those use the ratio estimator pinned to the
    integrator, which is insensitive to a proportional loss.  What changes is
    only what this gate is able to see.
    """
    sum_w = run_stats.get("sum_w", 0.0)
    sum_nt = run_stats.get("sum_ntrials")
    n_wt = run_stats.get("n_with_trials", 0)
    if sum_nt and n_wt == n_events:
        return sum_w / sum_nt, "sum(w)/sum(ntrials)"
    if sum_nt:                      # partial: some events carried no trials
        return None, (f"trial counts on only {n_wt}/{n_events} events")
    if n_read:
        lost = n_read - n_events
        how = f"sum(w)/{n_read} offered"
        if lost > 0:
            how += f", {lost} lost in the shower ({100*lost/n_read:.2f}%)"
        return sum_w / n_read, how
    return (sum_w / n_events if n_events else 0.0),\
        "mean weight over DELIVERED events -- no .nevents sidecar"


def herwig_lhe_delivered_pb(files):
    """(sigma_delivered_pb, n_written, n_attempted) for a POWHEG-LHE sample
    showered by Herwig, read from each job's own .out file.

    ONE SOURCE FOR A NORMALISATION RULE (CONVENTIONS.md rule 2): the parsing lives
    in herwig7/herwig_xsec.py with the rest of Herwig's cross-section
    bookkeeping, and this is only the glob from events.hepmc paths to .out
    paths.
    """
    sys.path.insert(0, os.path.join(BASE, "herwig7"))
    from herwig_xsec import combine_lhe_jobs      # noqa: E402
    outs = []
    for fn in files:
        outs.extend(glob.glob(os.path.join(os.path.dirname(fn), "*-S*.out")))
    if not outs:
        return None, 0, 0
    return combine_lhe_jobs(outs)


def closure_check(gen, delivered_pb, integrator_pb, w, out, how="", note=""):
    """Compare the event-level cross-section against the pinned integrator one.

    Records the result in `out` and returns the ratio.  Raises SystemExit past
    CLOSURE_FAIL unless ALLOW_BAD_CLOSURE=1 is set in the environment (for
    deliberately analysing a known-broken sample).

    note: a known, physical reason for a deviation (e.g. MC@NLO's S/H-event
    split), printed alongside a warning so that an expected offset does not
    read as a defect -- and so nobody learns to ignore this warning.
    """
    if delivered_pb is None:
        print(f"{gen}: delivery closure NOT CHECKED -- {how}")
        out["integrator_closure"] = None
        return None
    ratio = delivered_pb / integrator_pb if integrator_pb else float("nan")
    # A NON-FINITE CLOSURE IS A FAILURE, NOT A PASS (pipeline review
    # 2026-09-19, finding C3).  A zero/None integrator or a NaN weight gives
    # ratio = nan, every comparison below is then False, and the result used
    # to be written as clean.
    if not (math.isfinite(ratio) and math.isfinite(delivered_pb)):
        if os.environ.get("ALLOW_BAD_CLOSURE") != "1":
            sys.exit(f"{gen}: FAILED -- non-finite delivery closure "
                     f"(delivered {delivered_pb!r}, integrator "
                     f"{integrator_pb!r}, {how}).  A NaN weight or a missing "
                     f"integrator; nothing is written.")
    # Deliberately no "n sigma" here.  For the sum(w)/sum(ntrials) estimator
    # the dominant fluctuation is in the TRIAL COUNTS, not the weights, so the
    # naive std(w)/sqrt(N) understates the error by orders of magnitude and
    # would flag a harmless 0.05% deviation as "21 sigma".  The percentage
    # deviation against the bands below is the honest gate; the integrator's
    # own quoted precision is ~0.1% for these runs.
    out["sigma_delivered_pb"] = delivered_pb
    out["sigma_delivered_method"] = how
    out["integrator_closure"] = ratio

    # nb for the ~38 nb muon runs, pb for the ~3 pb neutrino ones
    scale, unit = ((1e3, "nb") if min(delivered_pb, integrator_pb) >= 1e3
                   else (1.0, "pb"))
    msg = (f"{gen}: delivery closure = {ratio:.4f} "
           f"({100*(ratio-1):+.2f}%; events {delivered_pb/scale:.4f} {unit} "
           f"[{how}] vs integrator {integrator_pb/scale:.4f} {unit})")
    dev = abs(ratio - 1.0)
    decl = declared_deficit(gen)
    covered = decl is not None and dev <= decl["max_dev"]
    if covered:
        out["closure_declared_deficit"] = {
            "max_dev": decl["max_dev"],
            "cause": decl["cause"],
            "evidence": decl["evidence"],
        }
    if dev > CLOSURE_FAIL and not covered \
            and os.environ.get("ALLOW_BAD_CLOSURE") != "1":
        extra = ""
        if decl is not None:
            extra = (f"\n{gen}: a deficit IS declared for this sample, but "
                     f"only up to {100*decl['max_dev']:.1f}% and this is "
                     f"{100*dev:.2f}%.\n  A declared cause that is GETTING "
                     f"WORSE is not covered -- re-check it before widening "
                     f"the declaration.")
        sys.exit(msg + f"\n{gen}: FAILED -- the delivered events do not "
                 f"reproduce the integrator to {100*CLOSURE_FAIL:.0f}%, so the "
                 f"pinned normalisation is hiding a biased sample.  Check the "
                 f"job log for 'Massive PS flavours' and for a nonzero "
                 f"'From \"Jet_Evolution:CSS\"' discard rate (see the Sherpa "
                 f"trap notes)." + extra +
                 f"\n  If this deficit is understood AND verified not to bias "
                 f"sigma_fid, declare it in\n  analyze.DECLARED_DEFICITS with "
                 f"its size, cause and evidence.  ALLOW_BAD_CLOSURE=1 exists "
                 f"for a\n  one-off diagnostic and suppresses the gate for "
                 f"EVERY sample in the run.")
    if covered and dev > CLOSURE_FAIL:
        # printed in full every time: a declaration that scrolls past unread
        # is the same as no gate at all
        print(msg + f"  <-- DECLARED DEFICIT (up to "
                    f"{100*decl['max_dev']:.1f}%)")
        print(f"     cause:    {decl['cause']}")
        print(f"     evidence: {decl['evidence']}")
    elif dev > CLOSURE_WARN:
        tail = f"  <-- outside the 0.5% closure band"
        tail += f": {note}" if note else " -- WARNING"
        print(msg + tail)
    else:
        print(msg)
    return ratio


def dedup_dmesons(dcand, vin):
    """Collapse the record copies of each D meson to a single entry.

    Generators write a particle several times along its own history (a
    "self-chain": D+ -> D+ -> D+), and we must count it once. But a D from a
    D* -> D pi decay is a genuinely new meson that must count, and a D from a
    B decay likewise. The rule that separates the two: keep the LAST copy of
    each self-chain, i.e. drop any D that has a same-PDG daughter, and keep
    the one that actually decays to non-D children.

    dcand: list of (pdg, p4, particle_id, parent_field) for every D-flavoured
    record in the event, any status. `parent_field` is the HepMC3 ASCII P-line
    field that is either a production-vertex id (negative; mothers are then
    the vertex's incoming particles, from `vin`) or a mother particle id
    (positive). Returns [(pdg, p4)] for the surviving mesons.
    """
    if not dcand:
        return []
    byid = {rec[2]: i for i, rec in enumerate(dcand)}
    superseded = set()
    for pdg, _, _, parent in dcand:
        if parent > 0:
            mothers = (parent,)
        else:
            s = vin.get(parent)
            if not s:
                continue
            mothers = s.split(",")
        for m in mothers:
            j = byid.get(int(m))
            # a same-PDG mother is this particle's earlier copy, not its parent
            if j is not None and dcand[j][0] == pdg:
                superseded.add(j)
    return [(rec[0], rec[1]) for i, rec in enumerate(dcand)
            if i not in superseded]


def weakly_decaying_charm(chad, vin):
    """The PROMPT, WEAKLY-DECAYING charm hadrons of an event.

    WHY THIS IS NOT `dedup_dmesons` (2026-09-01).  That function answers "how
    many D mesons did this event make", and it keeps a D that came from a D*
    because that is a second, genuinely different meson.  This one answers a
    different question -- "what charm object would a detector reconstruct" --
    and there the D* and the D it decays to are ONE object seen twice.  So the
    rule here is stricter: drop any charm hadron that has a charm DAUGHTER of
    any species, not merely a same-PDG copy of itself.  What survives is the
    charm hadron at the end of every charm chain, which is the one that
    travels far enough to leave a displaced vertex.

    chad: [(pdg, p4, id, parent)] for every prompt charm-hadron record, in
    file order; the b -> c cascade is already excluded by the caller.
    vin: vertex id -> raw incoming-id string, as parse_hepmc3 collects it.
    Returns [(pdg, p4)].
    """
    if not chad:
        return []
    byid = {rec[2]: i for i, rec in enumerate(chad)}
    superseded = set()
    for _pdg, _p4, _pid, parent in chad:
        mothers = ((parent,) if parent > 0
                   else tuple(vin.get(parent, "").split(",")))
        for m in mothers:
            if m == "":
                continue
            j = byid.get(int(m))
            if j is not None:
                superseded.add(j)
    return [(rec[0], rec[1]) for i, rec in enumerate(chad)
            if i not in superseded]


D_HIST_KEYS = ("nd_ch", "nd_0", "Ed_ch", "Ed_0")


def dmeson_summary(obs, w):
    """Mean D multiplicities and the fraction of sigma with at least one D.

    Also sets "d_supported": a sample in which not a single D meson appears
    cannot express these observables at all -- e.g. the GENIE HepMC3 files are
    written by genie/gtohepmc3.cc with only kIStStableFinalState (status 1)
    and beam records, so every decayed particle, and hence every D, is absent
    from the file. Writing all-zero histograms for such a sample would draw a
    "no charm" curve that is really "no information", so main() drops the
    D histograms instead and the plotting code skips the sample.
    """
    w = np.asarray(w, dtype=float)
    sw = float(w.sum())
    out, n_tot = {}, 0
    for key, tag in (("nd_ch", "d_ch"), ("nd_0", "d_0"), ("nd_s", "d_s")):
        if key not in obs or len(obs[key]) != len(w):
            continue
        v = np.asarray(obs[key], dtype=float)
        n_tot += int(v.sum())
        out[f"mean_n_{tag}"] = float(np.average(v, weights=w)) if sw else 0.0
        out[f"frac_ge1_{tag}"] = float(w[v > 0].sum() / sw) if sw else 0.0
        out[f"frac_ge2_{tag}"] = float(w[v > 1].sum() / sw) if sw else 0.0
    out["n_d_records"] = n_tot
    out["d_supported"] = n_tot > 0
    return out

# ---------------------------------------------------------------- binnings
def bins_for(e_lab):
    """Histogram edges at one beam energy.

    THE ENERGY-LIKE AXES MUST SCALE WITH THE BEAM.  They were literals tuned to
    1 TeV, and np.histogram DISCARDS out-of-range entries without a word: at
    4 TeV, nu runs 800-3600 GeV against edges of 200-900, so nearly every event
    would vanish from the spectrum -- while sigma_fid, which is summed from the
    weights and not from the histogram, stayed perfectly correct.  Exactly the
    silent-failure shape CONVENTIONS.md rule 2 is about.

    Scalings, all exact at the 1 TeV anchor so no published plot moves:
      nu    = y * E        -> [0.2 E, 0.9 E], i.e. the y range itself
      Emu   = (1 - y) E    -> [0.1 E, 0.8 E]
      Elead, Ed_*          -> upper edge 0.9 E (the nu endpoint)
      Q2max = y * 2k.P     -> proportional to E
      xmin  = Q2min/(y 2k.P) -> inversely proportional to E
      theta ~ sqrt(Q2/(E E')) -> inversely proportional to E
    The shape-only axes (y, nch, dphi, N(D)) are energy-independent.

    >>> AND THE y-LIKE AXES MUST FOLLOW THE SELECTION'S y WINDOW (2026-09-14).
    <<<  They were literals for 0.2 < y < 0.9, which is right for every
    selection that carries that window and wrong for every one that does not.
    In the no-y-cut final region (q4w3) the generator E_lep' histogram held
    35% of the muon sigma_fid and 70% of the neutrino one -- np.histogram drops
    the rest -- while the YADISM calculators CLIP out-of-range entries into
    the edge bins, so the reference's first y bin came out 24x too high.
    Neither errors.  The ranges are now [y_min, y_max] of the active
    selection at the SAME bin widths (dy = 0.025, 20 GeV at 1 TeV), which
    reproduces the old edges exactly for a 0.2-0.9 window:
      y     = [y_min, y_max]
      nu    = [y_min E, y_max E]     (Ehad likewise)
      Emu   = [(1 - y_max) E, (1 - y_min) E]
    The hadron-energy axes (Elead, Ed_*) keep their 0.9 E endpoint: they feed
    published earlier no-y-cut figures, and nothing needs them moved.
    """
    r = e_lab / 1000.0                                  # 1 at the anchor
    y0, y1 = SELECTION.y_min, SELECTION.y_max
    ny = int(round((y1 - y0) / 0.025)) + 1              # 29 for 0.2-0.9
    ne = int(round((y1 - y0) / 0.020)) + 1              # 36 for 0.2-0.9
    return {
    # THE FIRST Q2 EDGE IS THE FIDUCIAL FLOOR ITSELF, read from the active
    # selection rather than from the environment.  It used to read
    # BENCH_Q2MIN directly, which moved the BINNING without moving the CUT --
    # so a Q2 > 11 run histogrammed a Q2 > 4 sample from 11 upwards and looked
    # entirely reasonable.  selection.Q2_FLOOR is now the one place the floor
    # is set, and every result computed above the default carries the value in
    # its filename (selection.q2_suffix) so the two can never be confused.
    "Q2":      np.geomspace(SELECTION.q2_min, 1600.0 * r, 27),   # GeV2
    "xbj":     np.geomspace(2e-3 / r, 1.0, 27),
    "y":       np.linspace(y0, y1, ny),
    # Total hadronic energy in the lab, on the same grid as `nu` because on a
    # FREE PROTON the two are the same number up to the proton mass -- which
    # makes their ratio a closure test here and a measure of what the nucleus
    # absorbs once these samples are run on one.
    "Ehad":    np.linspace(y0 * 1000.0 * r, y1 * 1000.0 * r, ne),
    # (1000 - 1000 y), not (1 - y) 1000: the latter is 99.99999999999997
    "Emu":     np.linspace((1000.0 - y1 * 1000.0) * r,
                           (1000.0 - y0 * 1000.0) * r, ne),  # GeV, lab
    "theta":   np.geomspace(1e-3 / r, 4e-1 / r, 27),   # rad, lab
    # width-2 bins: each contains one odd N (even N is ~1% by charge conserv.)
    "nch":     np.arange(-0.5, 31.5, 2.0),
    "nch1":    np.arange(-0.5, 31.5, 2.0),             # E_lab > 1 GeV
    # >>> THE NUCLEAR OBSERVABLES (user, 2026-09-08). <<<  Added for the
    # "Comparison with Neutrino Generators" study, where GENIE is set beside
    # the dedicated neutrino generators.  Fig. 7.20 of the FPF whitepaper
    # (arXiv:2203.05090) compares GENIE, NEUT and NuWro at 1 TeV on tungsten
    # and its finding is about NUCLEON multiplicity -- GENIE's HEDIS/CSMS mode
    # sits low because it has no final-state interactions -- and nothing in
    # this benchmark counted nucleons at all: `nch` counts CHARGED hadrons
    # (protons among them, indistinguishably) and `nd_*` counts D mesons.
    #
    # ANTIPARTICLES ARE COUNTED WITH PARTICLES.  A knocked-out nucleon and an
    # antinucleon from the string are different physics, but the quantity the
    # nuclear generators differ on is the number of baryons leaving the
    # vertex, and splitting the count would halve the statistics of the
    # channel that matters.  Said here because it is a choice.
    "nprot":   np.arange(-0.5, 24.5, 1.0),             # |pdg| = 2212
    "nneut":   np.arange(-0.5, 24.5, 1.0),             # |pdg| = 2112
    "npi":     np.arange(-0.5, 30.5, 1.0),             # |pdg| = 211
    # THE MULTIPLICITY THE EMULSION ACTUALLY COUNTS (2026-09-08).  FASER's
    # vertex selection keeps tracks with tan(theta) <= 0.5 (arXiv:2403.12520),
    # so a multiplicity over ALL charged hadrons is not the measured one --
    # it counts tracks that leave the emulsion's angular acceptance.  This is
    # the same E_lab > 1 GeV threshold as nch1 with that acceptance applied,
    # and it is a NEW key rather than a change to nch1, which several
    # published figures use.
    # >>> INTEGER BINS, unlike nch and nch1 above (user, 2026-09-08: "for the
    # N_ch, maybe use integer bins since this is a discrete quantity"). <<<
    # The width-2 grid exists because on the FULL charged-hadron multiplicity
    # the target's charge makes N odd -- every even bin would be ~1% -- so
    # pairing them hides a comb that is bookkeeping rather than physics.  That
    # argument does NOT carry over here: nch05 counts only hadrons above 1 GeV
    # inside tan(theta) < 0.5, and the two cuts break the parity, so the even
    # bins are populated and the distribution is smooth on an integer grid.
    "nch05":   np.arange(-0.5, 31.5, 1.0),             # E_lab > 1, tan < 0.5
    "nu":      np.linspace(y0 * 1000.0 * r, y1 * 1000.0 * r, ne),  # GeV
    "Elead":   np.geomspace(1.0, 900.0 * r, 27),       # GeV, leading ch. hadron
    # min |phi_h - phi_mu| over the charged hadrons of the N_ch selection;
    # 30 uniform bins of width pi/30 = 0.105 rad
    "dphi":    np.linspace(0.0, math.pi, 31),          # rad, lab
    # >>> THE BACK-TO-BACK VARIABLE FASER CUTS ON, which is NOT `dphi`. <<<
    # The emulsion analysis requires Delta phi > pi/2 between the lepton and
    # the SUM of the other tracks in the vertex; `dphi` above is the minimum
    # over INDIVIDUAL hadrons, a different quantity that happens to have the
    # same name in ordinary speech.  Tier E has cut on the summed one since it
    # was written (selection.passes_hadrons), but it was never histogrammed,
    # so the observable the measurement is defined by could not be drawn.
    # Finer than dphi because Tier E keeps only the upper half of the range.
    "dphix":   np.linspace(0.0, math.pi, 49),          # rad, lab
    # D mesons at production (see dedup_dmesons). Width-1 integer bins so the
    # N = 0 bin -- which holds ~95% of sigma -- is visible on its own.
    "nd_ch":   np.arange(-0.5, 6.5, 1.0),              # N(D+-) per event, 0..5
    "nd_0":    np.arange(-0.5, 6.5, 1.0),              # N(D0/D0bar), 0..5
    # lab energy of the LEADING D of each species, one entry per event (the
    # per-event form keeps every consumer of BINS length-safe; see NOTE below).
    # Range from just under threshold (m_D = 1.87) to the Elead endpoint.
    "Ed_ch":   np.geomspace(2.0, 900.0 * r, 17),       # GeV, lab
    "Ed_0":    np.geomspace(2.0, 900.0 * r, 17),       # GeV, lab
    }


BINS = bins_for(ENERGY)
# NOTE on the D energy spectra: every script that consumes BINS histograms
# obs[key] against the per-EVENT weight array, so a key holding one entry per
# D meson (rather than per event) would raise a length mismatch in scripts
# owned elsewhere (herwig7/make_histos_{nlo,me}.py). Ed_ch/Ed_0 are therefore
# the leading-D energy per event, exactly parallel to the existing Elead. As
# with Elead the "no such meson" case needs no sentinel: the value stays 0.0
# and falls below the first bin edge, so np.histogram drops it.
LOGBINS = {"Q2", "xbj", "theta", "Elead"}


def dphi_summary(obs, w):
    """Book-keeping for the DeltaPhi observable: how many fiducial events had
    no charged hadron at all (so DeltaPhi is undefined and not histogrammed),
    and the mean of the ones that were filled."""
    v = np.asarray(obs["dphi"], dtype=float)
    w = np.asarray(w, dtype=float)
    miss = v < 0.0
    sw = float(w.sum())
    return {"n_no_hadron": int(miss.sum()),
            "frac_no_hadron": float(w[miss].sum() / sw) if sw else 0.0,
            "mean_dphi": (float(np.average(v[~miss], weights=w[~miss]))
                          if (~miss).any() else float("nan"))}


def dot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


# --------------------------------------------------- DIS kinematics, ONE copy
# These three were written out inline, IDENTICALLY, in analyze.py and
# analyze_nu.py, and powheg_nlo_uncertainty.py would have been the third copy.
# CONVENTIONS.md rule 2 names exactly this: Herwig's cross-section bookkeeping was
# duplicated and the copies diverged.  A fiducial region that drifts between
# two analysers is invisible in any single number, because each one is
# self-consistent.
#
# They are split at the point the Q2/y cut is applied rather than returning
# everything at once, so the ORDER of operations is unchanged: theta costs a
# hypot, a sqrt and an asin, and the invariant cut rejects most events before
# it is ever reached.
def lab_energy(p, P):
    """Scattered-lepton energy in the PROTON REST FRAME, computed covariantly.

    The event record is not in that frame -- POWHEG generates a symmetric
    massless pair at sqrt(s)/2 and Pythia puts the beams back on shell -- so
    this cannot read an energy component.  dot(p, P)/M_P is frame-independent
    and gives the lab energy the benchmark's cuts are specified in.
    """
    return dot(p, P) / M_P


def dis_invariants(k, P, lep):
    """(Q2, y, x_Bjorken, k.P) from the two beams and the scattered lepton.

    y is built from the PROTON momentum, never from the struck parton's.
    Herwig's SimpleDISCut builds it from x_parton instead, which is why a
    generation-level y cut there does not cut the y an analysis reconstructs
    (see CONVENTIONS.md rule 2).
    """
    kP = dot(k, P)
    q = (k[0]-lep[0], k[1]-lep[1], k[2]-lep[2], k[3]-lep[3])
    Q2 = -dot(q, q)
    Pq = dot(P, q)
    return Q2, Pq / kP, Q2 / (2.0*Pq), kP


def lepton_theta(e_lab, lep):
    """Scattered-lepton polar angle about the beam axis.

    p_T is invariant under the boost between the generation frame and the lab
    because that boost is along z, and e_lab is covariant, so this is correct
    in either frame without an explicit boost.  It cannot distinguish forward
    from backward, which is immaterial for the > 100 GeV muons it is asked
    about.
    """
    pt = math.hypot(lep[1], lep[2])
    pmag = math.sqrt(max(e_lab*e_lab - M_MU*M_MU, 1e-12))
    return math.asin(min(pt/pmag, 1.0))


def track_tan(p, e_lab, P, k):
    """(p_T, tan theta_lab) of a particle about the beam, in the target rest frame.

    p_T is invariant under the boost along the beam between the generation
    frame and the lab, and |p_z,lab| = sqrt(E_lab^2 - m^2 - p_T^2) follows from
    the covariant lab energy -- the construction every Tier E count has always
    used, kept so that no forward track moves.  WHAT IT CANNOT GIVE IS THE SIGN
    of p_z,lab, and until 2026-09-19 it was taken positive: a soft hadron going
    BACKWARD in the lab, target fragmentation, counted as a forward track
    (review S1; <= 0.11% of the Tier E rate).  The sign now comes from the
    covariant
        p_z,lab = (E_lab E_k - p.k) / |k_lab|,   E_k = k.P / M_P,
    and a particle that is not forward gets tan = +inf, i.e. no track.
    """
    pt = math.hypot(p[1], p[2])
    e_k = dot(k, P) / M_P
    if e_lab * e_k - dot(p, k) <= 0.0:
        return pt, math.inf
    m2 = max(dot(p, p), 0.0)
    pz = math.sqrt(max(max(e_lab * e_lab - m2, 0.0) - pt * pt, 0.0))
    return pt, (pt / pz if pz > 0 else math.inf)


def reconstruct_beams(mu3):
    """Beams for events without status-4 records (Sherpa BEAM_REMNANTS: false):
    the beam muon is the hardest status-3 muon (no QED radiation, no lepton
    PDF), and the proton beam is back-to-back with it in the generation
    (c.m.) frame with E fixed by the proton mass."""
    e_p = math.sqrt(mu3[1]**2 + mu3[2]**2 + mu3[3]**2 + M_P*M_P)
    return mu3, (e_p, -mu3[1], -mu3[2], -mu3[3])


# ------------------------------------------------------- hard-process flavour
# HepMC3 keeps Pythia's own status codes for non-final particles, so the hard
# 2 -> 2 scattering is identifiable in the record: status 21 = incoming to the
# hard process (after backward evolution of the shower, i.e. the flavour that
# actually enters the matrix element), status 23 = outgoing from it. This is
# "particle level before showering and hadronisation" in the sense required for
# a channel tag: it is the ME flavour, unaffected by what the shower does
# afterwards. The two incoming records are the lepton and the struck parton;
# hard_partons() drops the leptons and returns the QCD flavours.
HARD_IN_STATUS = 21
HARD_OUT_STATUS = 23
from hardproc import (CC_BOSONS, NC_BOSONS, LEPTON_PIDS,  # noqa: E402
                      find_hard_vertex, hard_partons)


def parse_hepmc3(fname, beam_pid=13, run_stats=None, hard_graph=False):
    """Yield (weight, mu_beam, p_beam, particles, dmesons, hard) per event;
    momenta are (E,px,py,pz). particles = final-state (pid, p4) list; dmesons =
    de-duplicated (pid, p4) list of D mesons counted at production; hard =
    {"in": [pdg...], "out": [pdg...]} of the hard-process records (statuses 21
    and 23; empty for generators that do not write them).
    beam_pid selects the status-4 lepton beam record (13 = mu-, 14 = nu_mu).

    run_stats: optional dict, accumulated across files.  Sherpa's weight block
    is "Weight | MEWeight | WeightNormalisation | NTrials | ...", and its
    cross-section is sum(Weight) / sum(NTrials) -- NOT the mean weight, since
    an event can take several trials to be accepted (mean ~4 here).  The
    trial counts are collected so that the delivery-closure gate can compare
    the delivered cross-section against the integrator's.  Generators that
    write a bare single weight (POWHEG here) contribute no trials, and the
    gate falls back to the mean weight for them."""
    w = None
    ntrials = None
    mu_beam = p_beam = mu3 = None
    pall = {}           # particle id -> (pdg, status), only if hard_graph
    vout = {}           # vertex id -> [outgoing particle ids]
    vgraph = {}         # vertex id -> [incoming particle ids]
    parts = []
    dcand = []          # every D-flavoured record, any status
    chad = []           # every prompt charm-HADRON record (hard_graph only)
    vin = {}            # vertex id -> raw "id,id,..." of its incoming particles
    hard_in, hard_out = [], []
    # Heavy flavour in the FINAL STATE of the event.  fs_heavy holds 4 when
    # the event has PROMPT charm and 5 when it has bottom at all.
    #
    # A hadron-level sample contributes charm/bottom HADRONS (any status --
    # they always decay under the ctau > 10 mm convention); an ME-level sample
    # has no hadronisation, so there a status-1 heavy QUARK is what "heavy
    # quark in the final state" means.  Both are collected here.
    #
    # THE b -> c CASCADE MUST BE EXCLUDED.  A bottom hadron always decays, and
    # it decays to charm, so a naive "is there a charm hadron" test makes
    # every b event a charm-production event.  Ancestry is therefore tracked:
    # b_tainted holds the ids of bottom hadrons and everything descended from
    # them, and a charm hadron counts as prompt only if it is not in that set.
    # A crude "any b vetoes the event" rule is NOT good enough -- in the
    # POWHEG-V2 neutrino sample 17% of the genuine charm events also contain a
    # spectator B, and vetoing them threw away real signal.
    #
    # The forward pass is exact because HepMC3 Ascii writes a vertex before
    # the particles it produces (verified on all five generators here), so a
    # particle's mothers are always already classified.  The bookkeeping costs
    # nothing until the first bottom hadron appears, and most events have none.
    fs_heavy = set()
    b_tainted = set()       # particle ids: a bottom hadron or its descendant
    v_tainted = set()       # vertex ids with a tainted incoming particle
    have_event = False

    def finish():
        hard = {"in": hard_in, "out": hard_out, "fs_heavy": frozenset(fs_heavy),
                "lhe_index": lhe_index[0]}
        if hard_graph:
            hard["charm_final"] = weakly_decaying_charm(chad, vin)
        if hard_graph and not hard_in and not hard_out:
            # Sherpa / Herwig: no Pythia status codes, walk the graph instead
            bosons = CC_BOSONS if abs(beam_pid) in (12, 14, 16) else NC_BOSONS
            hin, hout = find_hard_vertex(pall, vgraph, vout, beam_pid, bosons)
            hard = {"in": hin, "out": hout,
                    "fs_heavy": frozenset(fs_heavy),
                    "charm_final": weakly_decaying_charm(chad, vin),
                    "lhe_index": lhe_index[0]}
        if mu_beam is not None and p_beam is not None:
            return (w, mu_beam, p_beam, parts, dedup_dmesons(dcand, vin), hard)
        if mu_beam is None and p_beam is None and mu3 is not None:
            return ((w,) + reconstruct_beams(mu3)
                    + (parts, dedup_dmesons(dcand, vin), hard))
        return None

    lhe_index = [None]
    with open(fname) as f:
        for line in f:
            c0 = line[0:2]
            if c0 == "E ":
                if have_event:
                    ev = finish()
                    if ev is not None:
                        yield ev
                w, mu_beam, p_beam, mu3 = None, None, None, None
                ntrials = None
                pall, vout, vgraph = {}, {}, {}
                parts, dcand, vin = [], [], {}
                chad = []
                hard_in, hard_out = [], []
                fs_heavy, b_tainted, v_tainted = set(), set(), set()
                lhe_index = [None]
                have_event = True
            elif c0 == "A " and line.startswith("A 0 lhe_index "):
                # main_powheg's Les Houches event number: the join key for
                # every reweighting (memory: lhe-index-and-wrong-lhe), so it
                # travels with the event rather than being re-derived
                lhe_index[0] = int(line.split()[3])
            elif c0 == "W " and w is None:
                f_ = line.split()
                try:
                    w = float(f_[1])
                except ValueError:
                    pass          # the "W Weight|EXTRA__..." header line
                else:
                    if run_stats is not None:
                        run_stats["sum_w"] = run_stats.get("sum_w", 0.0) + w
                        # sum_w2 was accumulated here and read NOWHERE -- the
                        # only difference tools/check_parser_drift.py found
                        # between the two copies of this function. Removed
                        # rather than mirrored into analyze_nu.py: the closure
                        # gate deliberately quotes no "n sigma" (see
                        # closure_check), so there is nothing for it to feed.
                        if len(f_) > 4:       # Sherpa: NTrials is column 4
                            try:
                                ntrials = float(f_[4])
                            except ValueError:
                                ntrials = None
                            else:
                                run_stats["sum_ntrials"] = (
                                    run_stats.get("sum_ntrials", 0.0) + ntrials)
                                run_stats["n_with_trials"] = (
                                    run_stats.get("n_with_trials", 0) + 1)
            elif c0 == "P ":
                c = line.split()
                pid, status = int(c[3]), int(c[9])
                if hard_graph:
                    pall[int(c[1])] = (pid, status)
                    par = int(c[2])
                    if par < 0:            # produced at a vertex
                        vout.setdefault(par, []).append(int(c[1]))
                if status == 1 or status == 4:
                    p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                    if status == 4:
                        if pid == beam_pid:
                            mu_beam = p
                        elif pid == 2212 or pid == 2112:
                            # 2112: the FASER-rate ladders scatter off a free
                            # NEUTRON too (tools/genie_dimuon_ladder.sh); the
                            # kinematics use M_P, 0.14% off, immaterial there
                            p_beam = p
                    else:
                        parts.append((pid, p))
                elif status == 3 and pid == 13:
                    p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                    if mu3 is None or p[0] > mu3[0]:
                        mu3 = p
                elif status == HARD_IN_STATUS:
                    hard_in.append(pid)
                elif status == HARD_OUT_STATUS:
                    hard_out.append(pid)
                # deliberately outside the status branches: charm/bottom
                # hadrons are always decayed, never status 1
                if is_heavy_hadron(pid, 5) or (status == 1 and abs(pid) == 5):
                    fs_heavy.add(5)
                    b_tainted.add(int(c[1]))
                elif b_tainted:
                    par_ = int(c[2])
                    if (par_ in v_tainted) if par_ < 0 else (par_ in b_tainted):
                        b_tainted.add(int(c[1]))
                    elif 4 not in fs_heavy and (
                            is_heavy_hadron(pid, 4) or
                            (status == 1 and abs(pid) == 4)):
                        fs_heavy.add(4)
                elif 4 not in fs_heavy and (is_heavy_hadron(pid, 4) or
                                            (status == 1 and abs(pid) == 4)):
                    fs_heavy.add(4)
                # D mesons are always decayed (status 2), never final state,
                # so this test is deliberately outside the status branches
                if pid in D_PIDS_SIGNED:
                    dcand.append(
                        (pid, (float(c[7]), float(c[4]), float(c[5]),
                               float(c[6])), int(c[1]), int(c[2])))
                # EVERY prompt charm hadron, with its momentum, for the charm
                # ANATOMY (2026-09-01): which charm the hadron-level tag adds
                # over the hard-process one, and whether that charm is soft.
                # Opt-in via hard_graph so the production path is untouched.
                if hard_graph and is_heavy_hadron(pid, 4) \
                        and int(c[1]) not in b_tainted:
                    chad.append(
                        (pid, (float(c[7]), float(c[4]), float(c[5]),
                               float(c[6])), int(c[1]), int(c[2])))
            elif c0 == "V ":
                # "V <id> <status> [<in>,<in>] @ x y z t"; the incoming list is
                # kept as a raw string and only split when a D needs its
                # mothers, which is a few per mille of vertices
                c = line.split()
                raw = c[3][1:-1]
                vin[int(c[1])] = raw
                if b_tainted and any(int(t) in b_tainted
                                     for t in raw.split(",") if t):
                    v_tainted.add(int(c[1]))
                if hard_graph:
                    vgraph[int(c[1])] = [int(t) for t in raw.split(",") if t]
        if have_event:
            ev = finish()
            if ev is not None:
                yield ev


def n_opposite_sign_muons(parts, P, lep_out, sel):
    """Muons of charge OPPOSITE to the primary that pass the SAME cuts as it.

    The one implementation for both currents -- analyze_nu.py imports this
    rather than keeping a copy, because the two event loops have already
    drifted apart twice on exactly this kind of addition (the DeltaPhi
    observable, which the neutrino side simply lacked, and the sigma_fid
    pinning fix, which had to be made twice).  See CONVENTIONS.md rule 2b.

    The kinematics are built exactly as the primary lepton's are, a few lines
    below: p_T from the generation frame, which a boost along z leaves alone,
    and the lab energy from the invariant p.P/M_p.  So theta is the LAB angle,
    which is the one the cut means.  It inherits the primary's small-angle
    convention too -- asin(p_T/|p|) cannot tell a forward muon from a backward
    one -- which is immaterial for the > 100 GeV muons this is ever asked about.
    """
    n = 0
    for pid, p in parts:
        if pid == -lep_out:
            e2 = dot(p, P) / M_P
            pt2 = math.hypot(p[1], p[2])
            pm2 = math.sqrt(max(e2 * e2 - M_MU * M_MU, 1e-12))
            if sel.passes_lepton(e2, math.asin(min(pt2 / pm2, 1.0))):
                n += 1
    return n


# SHERPA'S SHOWER CHARM IS SUBTRACTED (user decision, 2026-09-03).  The
# benchmark's RESPECT_MASSIVE_FLAG: true, which Sherpa needs so its ME -> shower
# interface does not discard 9-23% of the events, makes charm MASSLESS in its
# parton shower too ("Massive PS flavours: (none)"), so g -> c cbar splits with
# no charm-mass threshold: 12-16x the rate of the massive showers Pythia and
# Herwig run (analysis/charm_shower_mass_check.py).  Those events inflate
# Sherpa's charm cross-section by 3-5% at 1 TeV and 9% at 4 TeV, and its
# D-meson yield by twice that.  So for Sherpa, and only Sherpa, the published
# charm tag is "a prompt charm hadron in the final state AND charm in the hard
# record" (hard process or initial state, i.e. the anatomy's hard + remnant
# pieces), the D-meson observables drop the D mesons of shower-only charm
# events, and the dimuon tier vetoes them.  The raw final-state number is kept
# beside the subtracted one in every result (sigma_fid_pb_raw_final).
SHOWER_CHARM_SUBTRACTED_PREFIXES = ("sherpa", "v2_sherpa")


def subtracts_shower_charm(gen):
    return gen.startswith(SHOWER_CHARM_SUBTRACTED_PREFIXES)


def stamp_shower_charm(out, tag_stats, w, subtract):
    """Record whether the shower charm was subtracted, and the raw number.

    Every Sherpa result carries `shower_charm_subtracted`; a charm-tagged one
    also carries `sigma_fid_pb_raw_final`, the cross-section the plain
    final-state tag would have given, so the size of the subtraction is
    visible on every result rather than only in the anatomy."""
    out["shower_charm_subtracted"] = bool(subtract)
    if (subtract and tag_stats and "sum_w_final_raw" in tag_stats
            and out.get("sigma_fid_pb") and len(w) and float(w.sum())):
        out["sigma_fid_pb_raw_final"] = (out["sigma_fid_pb"]
                                         * tag_stats["sum_w_final_raw"]
                                         / float(w.sum()))
    return out


def analyze(files, beam_pid=13, lep_out=13, hard_flavour=None, tag_stats=None,
            run_stats=None, hard_tag_mode="out", subtract_shower_charm=False):
    """beam_pid/lep_out: incoming and scattered lepton PDG codes
    (13/13 for muon NC DIS; 14/13 for nu_mu CC; 14/14 for nu_mu NC).

    hard_flavour: if set (e.g. 4), keep only events whose HARD PROCESS has that
    |PDG| among its incoming partons -- i.e. tag the channel by the flavour of
    the matrix element, before any showering or hadronisation. n_events and
    sum_w_all still count ALL parsed events, so the caller's normalisation
    (sigma_generated / N_generated) is unchanged and the returned histograms are
    the channel's absolute contribution to sigma_fid.
    tag_stats: optional dict, filled with hard-process book-keeping.
    run_stats: optional dict, filled with sum_w / sum_ntrials over
    ALL parsed events, for the delivery-closure gate (see closure_check)."""
    obs = {k: [] for k in BINS}
    obs["nd_s"] = []      # D_s multiplicity: reported as a mean, not binned
    # THE LES HOUCHES EVENT NUMBER OF EVERY SURVIVING EVENT.  Not a histogram
    # axis -- it is the join key that lets a scale- or PDF-reweighted LHE be
    # married to the SHOWERED events, which is the only way to put a band on
    # an observable the LHE cannot compute (a hadron multiplicity, a Delta phi
    # against the hadron system).  -1 where the sample carries none, as GENIE
    # and Sherpa do.  Same precedent as nd_s above: an obs key that no BINS
    # entry consumes.
    obs["lhe_index"] = []
    weights = []
    n_events = 0          # all parsed events (before fiducial cuts)
    sum_w_all = 0.0
    if tag_stats is not None:
        tag_stats.setdefault("n_no_hard_record", 0)
        tag_stats.setdefault("n_tagged", 0)
        tag_stats.setdefault("n_in_out_mismatch", 0)
        tag_stats.setdefault("mismatch_examples", [])
        tag_stats.setdefault("in_flavours", {})
    for fn in files:
        print(f"  parsing {fn}", flush=True)
        want_hard = (hard_flavour is not None or tag_stats is not None
                     or subtract_shower_charm)
        for w, k, P, parts, dmes, hard in parse_hepmc3(
                fn, beam_pid=beam_pid, run_stats=run_stats,
                hard_graph=want_hard):
            if w is None:
                w = 1.0
            n_events += 1
            sum_w_all += w
            # shower-only charm: a prompt charm hadron in the event, none in
            # the hard record (see SHOWER_CHARM_SUBTRACTED_PREFIXES)
            shower_only = False
            raw_final_event = False
            if subtract_shower_charm:
                _hin = hard_partons(hard["in"])
                _hout = hard_partons(hard["out"])
                _fs = hard.get("fs_heavy", frozenset())
                shower_only = (4 in _fs and not any(abs(q) == 4 for q in _hin)
                               and not any(abs(q) == 4 for q in _hout))
            if hard_flavour is not None or tag_stats is not None:
                hin = hard_partons(hard["in"])
                hout = hard_partons(hard["out"])
                if tag_stats is not None:
                    if not hin:
                        tag_stats["n_no_hard_record"] += 1
                    for p in hin:
                        key = str(abs(p))
                        tag_stats["in_flavours"][key] = \
                            tag_stats["in_flavours"].get(key, 0) + 1
                if hard_flavour is not None:
                    # TAG MODE.  The benchmark definition of heavy-quark
                    # production is "at least one such quark in the FINAL
                    # state".  At LO the t-channel exchange conserves flavour,
                    # so incoming and outgoing tags coincide exactly -- but at
                    # NLO they do NOT: gamma*g -> c cbar has a GLUON incoming
                    # and charm only outgoing, and an incoming-flavour tag
                    # silently drops it.  Measured on POWHEG-RES: charm is
                    # 1.33% of incoming partons against 13.67% gluon, so the
                    # incoming tag captured only the charm-initiated piece.
                    # FINAL-STATE HEAVY FLAVOUR, with the b -> c cascade
                    # already removed by the ancestry pass in parse_hepmc3():
                    # "4 in fs_heavy" means PROMPT charm.  How that was
                    # diagnosed: on POWHEG-RES the naive test added 1338
                    # events beyond the hard-process tag yet pushed
                    # sigma_charm DOWN, because 1219 of them were b-INITIATED
                    # and carried net NEGATIVE weight (-0.166 nb).
                    fs_heavy = (hard.get("fs_heavy", frozenset())
                                if isinstance(hard, dict) else frozenset())
                    tagged_final = hard_flavour in fs_heavy
                    if hard_tag_mode == "in":
                        keep = any(abs(p) == hard_flavour for p in hin)
                    elif hard_tag_mode == "any":
                        keep = (any(abs(p) == hard_flavour for p in hin) or
                                any(abs(p) == hard_flavour for p in hout))
                    elif hard_tag_mode == "final":
                        # THE BENCHMARK DEFINITION (user, 2026-08-21): the
                        # heavy quark is in the final state AFTER shower and
                        # hadronisation.  Unlike every mode above, this needs
                        # nothing from the hard-process record at all --
                        # except for Sherpa, whose shower-only charm is
                        # subtracted (SHOWER_CHARM_SUBTRACTED_PREFIXES).
                        keep = tagged_final and not (
                            subtract_shower_charm and hard_flavour == 4
                            and shower_only)
                        # the raw (unsubtracted) tag keeps the shower-only
                        # events too; counted after the selection cuts below
                        raw_final_event = tagged_final and not keep
                        if raw_final_event:
                            keep = True         # dropped after the cuts, below
                    else:                       # "out"
                        keep = any(abs(p) == hard_flavour for p in hout) or \
                            any(abs(p) == hard_flavour for p, _ in parts)
                    if tag_stats is not None:
                        tag_stats["n_tagged_final"] = (
                            tag_stats.get("n_tagged_final", 0)
                            + int(tagged_final))
                        tag_stats["n_charm_with_bottom"] = (
                            tag_stats.get("n_charm_with_bottom", 0)
                            + int(4 in fs_heavy and 5 in fs_heavy))
                        tag_stats["n_tagged_out"] = (
                            tag_stats.get("n_tagged_out", 0)
                            + int(any(abs(p) == hard_flavour for p in hout) or
                                  any(abs(p) == hard_flavour
                                      for p, _ in parts)))
                    if not keep:
                        continue
                    if tag_stats is not None:
                        tag_stats["n_tagged"] += 1
                        # LO t-channel gamma*/Z exchange conserves flavour, so
                        # the outgoing hard parton must be the same charm quark.
                        # A hard-process parton that neither showers nor
                        # hadronises (ME-level sample) can be written as a
                        # status-1 record instead of status 23, so final-state
                        # partons count as a match too.
                        ok = any(abs(p) == hard_flavour for p in hout) or \
                            any(abs(p) == hard_flavour for p, _ in parts)
                        if not ok and hard_tag_mode != "final":
                            # not a mismatch under the hadron-level tag: an
                            # event tagged on a shower g -> c cbar has no
                            # hard-process charm by construction
                            tag_stats["n_in_out_mismatch"] += 1
                            if len(tag_stats["mismatch_examples"]) < 20:
                                tag_stats["mismatch_examples"].append(
                                    {"file": os.path.basename(
                                        os.path.dirname(fn)),
                                     "in": hard["in"], "out": hard["out"]})
            # hardest scattered lepton by invariant lab energy
            best, best_elab = None, -1.0
            for pid, p in parts:
                if pid == lep_out:
                    elab = lab_energy(p, P)
                    if elab > best_elab:
                        best, best_elab = p, elab
            if best is None:
                continue
            Q2, y, xbj, kP = dis_invariants(k, P, best)
            if Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX:
                continue
            # W^2 = m_p^2 + 2 P.q - Q^2, and P.q = y * (k.P).  Only a selection
            # that sets w2_min pays for this; for every benchmark selection
            # W2_MIN is None and the branch is one comparison against None.
            # It is NOT folded into the line above because that line is the
            # published region and is deliberately left in the form it has
            # always had.
            if W2_MIN is not None and M_P*M_P + 2.0*y*kP - Q2 < W2_MIN:
                continue
            e_lab = best_elab
            theta = lepton_theta(e_lab, best)
            # FASER tiers add scattered-lepton energy and angle requirements
            # on top of the leptonic invariants.  No-ops for the inclusive
            # selection, whose tier fields are all None.
            if not SELECTION.passes_lepton(e_lab, theta):
                continue
            # The DIMUON tier: a second, opposite-sign muon passing those same
            # cuts.  Guarded on SELECTION.dimuon so that no other selection
            # pays for the extra scan or can change behaviour because of it.
            if SELECTION.dimuon and not SELECTION.passes_dimuon(
                    n_opposite_sign_muons(parts, P, lep_out, SELECTION)):
                continue
            if shower_only:
                # Sherpa: the charm hadrons of this event exist only because
                # its shower splits g -> c cbar massless; they are not
                # counted, and neither is a dimuon they decay into
                if SELECTION.dimuon:
                    continue
                dmes = []
            ebeam_lab = kP / M_P

            # lab-frame azimuth of the scattered muon; the lab is reached from
            # the generation frame by a boost along the beam (z) axis, which
            # leaves phi -- and hence every dphi below -- unchanged
            phi_mu = math.atan2(best[2], best[1])

            nch = nch1 = 0
            elead = 0.0
            dphi_min = math.inf
            # FASER Tier E: charged-hadron angular counts, and the azimuth of
            # the summed hadron system.
            #
            # THE ANGLE MUST BE THE LAB ANGLE, and theta is NOT invariant under
            # the boost from the generation frame (Sherpa and Herwig generate
            # in the c.m.).  It is rebuilt from invariants only: p_T is
            # unchanged by a boost along z, the lab energy is p.P/M_p as
            # everywhere else here, and the mass is p.p, so
            #     |p|_lab = sqrt(E_lab^2 - m^2),  |p_z,lab| = sqrt(|p|^2 - p_T^2),
            # with the SIGN of p_z,lab from the covariant form: track_tan().
            # Azimuth needs no such care: a z-boost leaves phi alone, which is
            # why the existing dphi uses the generation-frame components.
            n05 = n01 = 0
            nch05 = 0
            sum_px = sum_py = 0.0
            # THE NUCLEAR OBSERVABLES.  Counted over EVERY final-state
            # particle and not only the charged hadrons, since a neutron is
            # neither.  Ehad sums the lab energy of everything that is not a
            # neutrino and then removes the primary lepton, so it is the
            # energy a calorimeter would call hadronic: photons from pi0
            # included (pi0 decays under the ctau > 10 mm convention), the
            # scattered lepton excluded, secondary leptons from heavy-flavour
            # decay included -- which is what a detector sees.
            nprot = nneut = npi = 0
            ehad = 0.0
            for pid, p in parts:
                a = abs(pid)
                if a == 2212:
                    nprot += 1
                elif a == 2112:
                    nneut += 1
                elif a == 211:
                    npi += 1
                if a not in NEUTRINO_PIDS:
                    ehad += dot(p, P) / M_P
                if abs(pid) in CHARGED_HADRONS:
                    elab_h = dot(p, P) / M_P
                    nch += 1
                    if elab_h > 1.0:
                        nch1 += 1
                    if elab_h > elead:
                        elead = elab_h
                    _pt, tan_h = track_tan(p, elab_h, P, k)
                    # a TRACK needs E > TRACK_E_MIN (1 GeV), as in FASER's
                    # emulsion reconstruction; Tier E had no threshold until
                    # 2026-09-19 (user: "sign fix and the threshold applied
                    # together").  nch05 always had it.
                    if tan_h < 0.5 and elab_h > TRACK_E_MIN:
                        n05 += 1
                        nch05 += 1
                        if tan_h < 0.1:
                            n01 += 1
                    sum_px += p[1]
                    sum_py += p[2]
                    # phi_h - phi_mu wrapped into (-pi, pi], then |.|
                    d = abs(math.remainder(math.atan2(p[2], p[1]) - phi_mu,
                                           2.0*math.pi))
                    if d < dphi_min:
                        dphi_min = d

            # Delta phi between the lepton and the SUMMED charged hadron
            # system: the emulsion analysis' back-to-back topology cut.  Not
            # the dphi_min above, which is the closest single hadron.
            if sum_px or sum_py:
                dphi_x = abs(math.remainder(
                    math.atan2(sum_py, sum_px) - phi_mu, 2.0 * math.pi))
            else:
                dphi_x = None
            if not SELECTION.passes_hadrons(n05, n01, dphi_x):
                continue

            # D mesons at production; energies in the lab via the invariant
            # p.P/M_p, as for every other hadron here
            nd_ch = nd_0 = nd_s = 0
            ed_ch = ed_0 = 0.0
            for pdg, p in dmes:
                a = abs(pdg)
                if a == D_CH_PID:
                    nd_ch += 1
                    e_d = dot(p, P) / M_P
                    if e_d > ed_ch:
                        ed_ch = e_d
                elif a == D_0_PID:
                    nd_0 += 1
                    e_d = dot(p, P) / M_P
                    if e_d > ed_0:
                        ed_0 = e_d
                else:
                    nd_s += 1

            vals = {"Q2": Q2, "xbj": xbj, "y": y, "Emu": e_lab,
                    "theta": theta, "nch": nch, "nch1": nch1,
                    "nch05": nch05,
                    "nu": ebeam_lab - e_lab, "Elead": elead,
                    "dphi": dphi_min if nch else NO_DPHI,
                    # NO_DPHI when there is no charged hadron at all: the
                    # sentinel sits below every bin edge, so np.histogram
                    # drops it rather than piling those events into bin 0.
                    "dphix": dphi_x if dphi_x is not None else NO_DPHI,
                    "nd_ch": nd_ch, "nd_0": nd_0, "nd_s": nd_s,
                    "Ed_ch": ed_ch, "Ed_0": ed_0,
                    "nprot": nprot, "nneut": nneut, "npi": npi,
                    # the primary lepton is removed here rather than skipped
                    # in the loop: it is identified by being the HIGHEST-energy
                    # lep_out of the event, which is not known until the loop
                    # that finds it has finished
                    "Ehad": ehad - e_lab}
            # >>> THE SUBTRACTED EVENT MUST NOT REACH `obs` EITHER. <<<
            # This block used to sit AFTER the loop that appends to obs, so a
            # Sherpa shower-only charm event was recorded in every observable
            # list and then skipped for `weights`.  The arrays then had
            # different lengths AND, worse, were MISALIGNED from the first
            # subtracted event on: obs[i] no longer belonged to weights[i].
            # np.average refuses shapes that differ, which is the only reason
            # this was ever visible -- it made the charm-tagged Sherpa results
            # fail to be written rather than be written wrong, and the failure
            # was swallowed by a driver that greps its output.  Found on
            # 2026-09-08 while asking why sherpa_dire_charmfinal had no Tier S
            # or Tier E result; introduced with the subtraction itself
            # (24407c1, 2026-09-03).  CONVENTIONS.md rule 2, once more.
            if raw_final_event:
                # a Sherpa shower-only charm event under the final tag: counted
                # in the raw sum, not in the published result
                if tag_stats is not None:
                    tag_stats["sum_w_final_raw"] = tag_stats.get("sum_w_final_raw", 0.0) + w
                continue
            for kk, v in vals.items():
                obs[kk].append(v)
            li = hard.get("lhe_index") if isinstance(hard, dict) else None
            obs["lhe_index"].append(-1 if li is None else int(li))
            weights.append(w)
            if tag_stats is not None and (raw_final_event or
                                          (hard_flavour == 4 and hard_tag_mode == "final")):
                # every kept charm event plus the subtracted shower-only ones
                tag_stats["sum_w_final_raw"] = tag_stats.get("sum_w_final_raw", 0.0) + w
    return obs, np.asarray(weights), n_events, sum_w_all


def main():
    # Before anything is normalised: every anchor literal must still
    # agree with the integration file it was taken from.  Cheap (four
    # small files), and it is the check whose absence let a
    # re-integrated anchor be normalised against the superseded
    # number on 2026-08-27.
    check_anchor_literals()
    gen = sys.argv[1] if len(sys.argv) > 1 else "pythia"
    # not every branch reads per-job sidecars; the manifest tolerates empty
    sig_pb, job_counts = [], []
    # Heavy-flavour tag modes, selected by the suffix on the sample name:
    #   <gen>_charmfinal  -- the BENCHMARK DEFINITION (user, 2026-08-21): the
    #                        quark is in the final state after shower and
    #                        hadronisation, i.e. a charm HADRON is present
    #                        (or, in an ME-level sample with hadronisation
    #                        off, a status-1 charm quark)
    #   <gen>_charmin     -- the LO charm-INITIATED definition
    #   <gen>_charmany    -- either side of the hard process
    #   <gen>_charm       -- out of the hard process (the original tag)
    # Same suffixes for "bottom".
    hard_tag_mode = "out"
    for _suf, _mode in (("final", "final"), ("in", "in"), ("any", "any")):
        for _fl in ("charm", "bottom"):
            if gen.endswith(_fl + _suf):
                gen = gen[:-len(_suf)]
                hard_tag_mode = _mode
                break
        if hard_tag_mode != "out":
            break
    beam_pid, lep_out = 13, 13
    hard_flavour, tag_stats = None, None
    # hard-process flavour tags: <sample>_charm / <sample>_bottom select the
    # charm- / bottom-INITIATED channel (gamma*/Z + q -> q off the heavy-quark
    # PDF), tagged on the matrix element before shower and hadronisation
    HARD_TAGS = {"charm": 4, "bottom": 5}
    # Pythia samples: <sample>[_<flavour tag>]. "pythia_me_massless" is the
    # MASSLESS-QUARK DIAGNOSTIC (pythia8/dis_mu1TeV_me_massless.cmnd, all five
    # quark masses set to 0) that puts Pythia's phase space in the same
    # zero-mass scheme as the YADISM reference; see that card's header.
    # The shower-model arms (pythia8/shower/, run_pythia.sh --shower ARM) slot
    # in here because they ARE Pythia LO samples: same card, same seed, same
    # unweighted normalisation, only PartonShowers:model differs.  Done at LO
    # deliberately -- on the POWHEG samples the veto hooks are a SimpleShower
    # feature and Vincia and Dire deliver essentially nothing
    # (logs/shower_probe_powheg.log).
    PYTHIA_JOBDIRS = {"pythia": "job", "pythia_me": "me_job",
                      "pythia_me_massless": "me_ml_job",
                      "pythia_shw_simple": "shwsimple_job",
                      "pythia_shw_vincia": "shwvincia_job",
                      "pythia_shw_dire": "shwdire_job",
                      # DIAGNOSTIC: simple shower with dipoleRecoil OFF, to
                      # separate the shower algorithm from the DIS recoil
                      # treatment.  See pythia8/shower/norecoil.cmnd.
                      "pythia_shw_norecoil": "shwnorecoil_job"}
    SHOWER_LABEL = {"pythia_shw_simple": "Pythia 8.311 LO, simple shower",
                    "pythia_shw_vincia": "Pythia 8.311 LO, Vincia",
                    "pythia_shw_dire": "Pythia 8.311 LO, Dire",
                    "pythia_shw_norecoil":
                        "Pythia 8.311 LO, simple shower, no dipole recoil"}
    pbase, ptag = gen, None
    if "_" in gen and gen.rsplit("_", 1)[-1] in HARD_TAGS:
        pbase, ptag = gen.rsplit("_", 1)
    if pbase in PYTHIA_JOBDIRS:
        jdir = PYTHIA_JOBDIRS[pbase]
        files = job_files(f"{BASE}/pythia8", at_energy(jdir), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/pythia8", at_energy(jdir), "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        quals = ["ME"] if "_me" in pbase else []
        if pbase.endswith("_massless"):
            quals.append("massless quarks")
        tag = ptag
        if pbase in SHOWER_LABEL:
            label = SHOWER_LABEL[pbase]
        if tag in HARD_TAGS:
            # heavy-quark-INITIATED scattering gamma*/Z + Q -> Q off the
            # heavy-quark PDF, tagged on the hard process (status-21/23
            # records) before shower and hadronisation. NOT photon-gluon
            # fusion, which is O(alpha_s) and absent at LO.
            hard_flavour, tag_stats = HARD_TAGS[tag], {}
            quals.append(tag)
        label = "Pythia 8.311 LO" + (f" ({', '.join(quals)})" if quals else "")
    elif pbase == "genie":
        # unweighted like pythia: per-job events.hepmc + events_xsec.json,
        # sigma_gen_mb evaluated from the gmkspl spline at 1 TeV
        files = job_files(f"{BASE}/genie", at_energy('job'), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/genie", at_energy('job'), "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        label = "GENIE G18_02a (EMDIS)"
        # GENIE has NO electromagnetic/neutral-current DIS charm generator:
        # EventGeneratorListAssembler offers only DIS-CC-CHARM (Aivazis with
        # slow rescaling), so the EMDIS list has no charm matrix element, and
        # AGKY's charm model produces none in fragmentation either.  The tag
        # is accepted anyway so the muon charm table can carry a MEASURED zero
        # rather than an assumed one.
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase in ("genie_nnpdf", "genie_nucc", "genie_nunc"):
        # further unweighted GENIE samples, same per-job layout as "genie"
        #
        # KEYED ON pbase, NOT gen: a charm-tagged request arrives here as
        # "genie_nnpdf_charm" and the old `gen in (...)` test rejected it with
        # "unknown generator", so the muon NNPDF row had no charm entry at any
        # energy while its GRV98 twin did.  The tag is accepted for the same
        # reason it is on the plain genie branch -- GENIE's EMDIS list has no
        # charm matrix element and AGKY makes none in fragmentation, so the
        # muon charm table carries a measured zero rather than a gap.
        jdir, label, (beam_pid, lep_out) = {
            "genie_nnpdf": ("nnpdf_job", "GENIE G18_02a (EMDIS, NNPDF4.0)",
                            (13, 13)),
            "genie_nucc":  ("nucc_job", "GENIE HEDIS BGR18 NLO (nu CC)",
                            (14, 13)),
            "genie_nunc":  ("nunc_job", "GENIE HEDIS BGR18 NLO (nu NC)",
                            (14, 14)),
        }[pbase]
        files = job_files(f"{BASE}/genie", at_energy(jdir), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/genie", at_energy(jdir), "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase in ("v2_genie", "v2_genie_nnpdf"):
        # >>> PAPER PLOTS V2: GENIE G18_02a EMDIS on ONE nucleon (user,
        # 2026-09-14, for pp03/pp04). <<<  the earlier production's genie / genie_nnpdf rows
        # with the target as the only change (genie/production/genie_job.sh mu_grv /
        # mu_nnpdf).  Unweighted, sigma_gen from this nucleon's own spline;
        # the build's EM Q2 > 4 floor is the only generation cut, so the q4w3
        # W cut is applied here.  The unweighted normaliser below is keyed on
        # the name -- "v2_genie" is listed there explicitly.
        _sel, _t = sample_layout.require()
        _k = pbase[3:]
        files = sample_layout.genie_files("mu", _k, _t, ENERGY, "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in sample_layout.genie_files("mu", _k, _t, ENERGY,
                                       "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"] * 1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        if not files:
            sys.exit(f"no complete GENIE {_k} jobs for {_t} at {ENERGY:g} GeV")
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        label = ("GENIE G18_02a (EMDIS)" if _k == "genie"
                 else "GENIE G18_02a (EMDIS, NNPDF4.0)")
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "v2_powheg":
        # >>> PAPER PLOTS V2 (user, 2026-09-13): POWHEG-RES + Pythia 8 in the
        # final region, ONE NUCLEON per run (PRODUCTION.md). <<<  Generated
        # at Born Q2 > 2.25 with no y cut -- WIDER than q4w3, deliberately:
        # the LHE writer's on-shell reshuffle moves the lepton, and a
        # generation cut at 4 measured 0.4% low (powheg/production/).  So this is the
        # wider-generation recipe, integrator x fiducial weight fraction, and
        # the integrator is read from THIS nucleon's own stage-3 file -- no
        # literal, since nothing is published that one would protect.
        # A neutron run has a genuine 2112 beam: POWHEG and Pythia each apply
        # the isospin swap to the proton set themselves.
        _sel, _t = sample_layout.require()
        files = job_files(f"{BASE}/powheg",
                          at_energy(sample_layout.powheg_job_base("mu", _t)),
                          "events.hepmc")
        _f = f"{sample_layout.powheg_rundir('mu', _t, ENERGY)}/pwg-0001-st3-stat.dat"
        _got = _powheg_res_grand_total(_f) if os.path.exists(_f) else None
        if _got is None:
            sys.exit(f"no POWHEG-RES grand totals in {_f}")
        sigma_sample_pb = _got[0] - _got[1]
        print(f"  [powheg {_t}] integrator {sigma_sample_pb/1e3:.4f} nb "
              f"(from {_f})")
        n_gen_tot = None
        label = "POWHEG-RES NLO + Pythia8"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif (pbase.startswith("v2_powheg_")
          and pbase[len("v2_powheg_"):] in sample_layout.ARMS):
        # >>> PAPER PLOTS 9-11: THE SHOWER ARMS (user, 2026-09-14). <<<
        # The SAME 1 TeV POWHEG-RES Les Houches files as `v2_powheg`, on the
        # arm SUBSET of its jobs (sample_layout.arm_jobnums): `sub` is the Pythia
        # baseline on those jobs, `qedfsr/qedfsrisr/qedfull` the QED re-showers
        # (powheg/production/qed_arms.sh), `hw` Herwig 7 (powheg/production/herwig_arm.sh).
        # Same integrator as v2_powheg -- each seed is an independent estimate
        # of the whole integral -- so the recipe is integrator x fiducial
        # weight fraction over the subset, closed against the subset's own
        # offered events (Pythia) or Herwig's own table (hw).
        _sel, _t = sample_layout.require()
        _arm = pbase[len("v2_powheg_"):]
        files = sample_layout.arm_files("mu", _arm, _t)
        _f = f"{sample_layout.powheg_rundir('mu', _t, ENERGY)}/pwg-0001-st3-stat.dat"
        _got = _powheg_res_grand_total(_f) if os.path.exists(_f) else None
        if _got is None:
            sys.exit(f"no POWHEG-RES grand totals in {_f}")
        sigma_sample_pb = _got[0] - _got[1]
        n_gen_tot = None
        label = {"sub": "POWHEG-RES NLO + Pythia8",
                 "qedfsr": f"POWHEG-RES NLO + Pythia8, {QED_ARMS['fsr']}",
                 "qedfsrisr": f"POWHEG-RES NLO + Pythia8, {QED_ARMS['fsrisr']}",
                 "qedfull": f"POWHEG-RES NLO + Pythia8, {QED_ARMS['full']}",
                 "hw": "POWHEG-RES NLO + Herwig7"}[_arm]
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "v2_sherpa_dire":
        # >>> PAPER PLOT 9: Sherpa MC@NLO with the DIRE shower, one
        # nucleon, 1 TeV (tools/sherpa_production.sh with SHOWER=Dire). <<<  The same
        # card as v2_sherpa run with SHOWER_GENERATOR Dire, integrated on its
        # own (the MC@NLO subtraction follows the shower), normalised to its
        # own BVI + RS exactly as v2_sherpa is.
        _sel, _t = sample_layout.require()
        _rd = sample_layout.sherpa_rundir("mu", _t, ENERGY) + "_Dire"
        files = job_files(_rd, "job", "evtfull")
        sigma_sample_pb = sherpa_arm_sigma_pb(os.path.relpath(_rd, SHERPA_RUNS),
                                              True)
        n_gen_tot = None
        label = "Sherpa 3.0.5 MC@NLO, Dire"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "v2_sherpa":
        # >>> PAPER PLOTS V2: Sherpa MC@NLO in the final region, one nucleon.
        # Generated at Q2 > 4 and INEL y > the W = 3 edge at the Q2 floor,
        # which is looser than W > 3 above it, so q4w3 is a SUBSET and the
        # generic weighted branch below takes the weight fraction (its subset
        # test fires on the selection name).  Normalised to BVI + RS from the
        # run's own integ.log, like every other Sherpa arm.
        _sel, _t = sample_layout.require()
        _rd = sample_layout.sherpa_rundir("mu", _t, ENERGY)
        files = job_files(_rd, "job", "evtfull")
        sigma_sample_pb = sherpa_arm_sigma_pb(os.path.relpath(_rd, SHERPA_RUNS),
                                              True)
        n_gen_tot = None
        label = "Sherpa 3.0.5 MC@NLO"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "sherpa":
        files = job_files(f"{SHERPA_RUNS}/{at_energy('MuonDIS_NLO')}", "job", "evtfull")
        sigma_sample_pb = sherpa_nlo_sigma_pb()
        n_gen_tot = None                            # weighted: use sum of weights
        label = "Sherpa 3.0.5 MC@NLO"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase in SHERPA_SHOWER_ARMS:
        rundir, is_nlo, lab = SHERPA_SHOWER_ARMS[pbase]
        files = job_files(f"{SHERPA_RUNS}/{rundir}", "job", "evtfull")
        sigma_sample_pb = sherpa_arm_sigma_pb(rundir, is_nlo)
        n_gen_tot = None                            # weighted: sum of weights
        label = lab
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "sherpa_lo":
        files = job_files(f"{SHERPA_RUNS}/{at_energy('MuonDIS_LO')}", "job", "evtfull")
        sigma_sample_pb = sherpa_lo_sigma_pb()
        n_gen_tot = None                            # weighted: use sum of weights
        label = "Sherpa 3.0.5 LO"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "mg5_ps":
        # MG5_aMC LO matrix element, showered and hadronised by Pythia 8
        # (mg5/run_mg5.sh).  Unweighted like the Pythia and GENIE samples:
        # per-job events.hepmc + events_xsec.json, the cross-section being
        # MG5's own integrated weight, which Pythia reports back unchanged
        # through the Les Houches file.
        #
        # THE SAMPLE IS GENERATED INSIDE THE FIDUCIAL REGION: Q2 and y are cut
        # at generation by mg5/dis_hooks.f, since MG5's run card has no DIS
        # cuts at all.  So n_fiducial/n_parsed is 1 by construction, and a
        # number materially below it would mean the shower is moving the
        # leptonic invariants -- which, with QED radiation off and the dipole
        # recoil on, it must not.
        #
        # ITS PARTNER IS THE PYTHIA LO ROW, not the parton-level `mg5_me` one:
        # the two carry the same shower, the same hadronisation, the same
        # parton distribution and the same scale, and differ in the matrix
        # element and the phase-space generation alone.
        files = job_files(f"{BASE}/mg5", at_energy('ps_job'), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/mg5", at_energy('ps_job'),
                            "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            # THE MEAN DELIVERED WEIGHT, NOT Pythia's sigma_gen.  The events
            # are unweighted and every Les Houches weight IS MG5's integrated
            # cross-section, so their mean is that number exactly.  Pythia's
            # own sigmaGen divides the same sum by the events it TRIED, which
            # includes the handful whose hard process it could not construct
            # (30 in 200000 here, "setting mass failed"), while the fiducial
            # fraction below is taken over the events DELIVERED -- mixing the
            # two denominators would lose that fraction twice.
            sig_pb.append(d["sum_weights_pb"] / d["n_accepted"])
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        label = "MG5_aMC 3.7.2 LO + Pythia8"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "sherpa_lo_me":
        # same process/cuts as MuonDIS_LO, shower+hadronisation off
        files = job_files(f"{SHERPA_RUNS}/{at_energy('MuonDIS_LO_ME')}", "job", "evtfull")
        sigma_sample_pb = SHERPA_LO_ME_SIGMA_PB
        n_gen_tot = None                            # weighted: use sum of weights
        label = "Sherpa 3.0.5 LO (ME)"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "powheg":
        # POWHEG-DIS NLO + Pythia 8.311 shower (powheg/run_powheg_shower.sh);
        # weighted +-sigma_abs events; generation region (Q2 > 2.25, full y)
        # is wider than the fiducial cuts
        files = job_files(f"{BASE}/powheg", at_energy('job'), "events.hepmc")
        sigma_sample_pb = powheg_res_sigma_pb()
        n_gen_tot = None                            # weighted: use sum of weights
        label = "POWHEG-RES NLO + Pythia8"
        if ptag in HARD_TAGS:
            # NLO, so the tag catches BOTH charm-initiated scattering
            # and gamma*g -> c cbar; that is the "at least one charm in
            # the final state" definition the benchmark uses.
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "powheg_hw":
        # >>> THE SAME POWHEG-RES LHE FILES, SHOWERED BY HERWIG 7. <<<
        # powheg/run_powheg_herwig.sh, card herwig7/LHE-mu.in.  Everything
        # about the normalisation is identical to `powheg` above -- same
        # integrator cross-section, same generation region -- because it is
        # the same file: the matrix element, the matching, the parton
        # distribution and the scale are held fixed by construction and the
        # only difference is the shower and the hadronisation model.  That is
        # what makes this the one comparison in the benchmark where the
        # hadronisation model is the thing that changes, rather than being
        # confounded with the generator (user, 2026-09-08).
        #
        # ONE CAVEAT TRAVELS WITH IT, and it is the same one Herwig's own NLO
        # sample carries: Herwig DISCARDS about 3.6% of these events when its
        # shower reconstruction fails, and the discards sit at the bottom of
        # the Q2 range rather than being a random subset (the measurement is
        # in herwig7/DIS-mu-POWHEG.in).  Normalising to the POWHEG integrator
        # restores the total rate and spreads the missing events over every
        # bin.  Each job records what it was offered and what it delivered in
        # delivery.json, so the fraction is on the record rather than inferred.
        files = job_files(f"{BASE}/powheg", at_energy('hw_job'),
                          "events.hepmc")
        sigma_sample_pb = powheg_res_sigma_pb()
        n_gen_tot = None                            # weighted: sum of weights
        label = "POWHEG-RES NLO + Herwig7"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif qed_arm_of(pbase):
        # A QED arm: the SAME LHE files as `powheg` above, re-showered with
        # photon radiation switched on in stages.  Everything about the
        # normalisation is therefore identical to `powheg` -- same integrator
        # cross-section, same generation region, same closure -- and the only
        # difference between this sample and the published one is the shower.
        arm = qed_arm_of(pbase)
        files = job_files(f"{BASE}/powheg", at_energy(f"qed{arm}_job"),
                          "events.hepmc")
        sigma_sample_pb = powheg_res_sigma_pb()
        n_gen_tot = None                            # weighted: sum of weights
        label = f"POWHEG-RES NLO + Pythia8, {QED_ARMS[arm]}"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    elif pbase == "powheg_v2":
        # POWHEG-V2 (powheg-cmass/nu-DIS-master) run on the MUON NEUTRAL
        # CURRENT -- a THIRD independent NLO matching on the muon side, beside
        # POWHEG-RES and Herwig.  Card powheg/cards/POWHEG-V2/mu1TeV-NC,
        # showered through powheg/main_powheg with powheg_mu1TeV.cmnd.
        #
        # Generation region (Born Q2 > 2.25 via q2cut, full y) is WIDER than
        # the fiducial cuts, exactly as for POWHEG-RES, so the same recipe
        # applies.  q2cut is MANDATORY here and absent from the CC cards: NC
        # has a photon pole and POWHEG's default 5d-2 silently returns 15
        # microbarn.
        #
        # This sample loses ~2.6% of its events inside Pythia (see
        # delivered_sigma_pb); the loss is +/- balanced, so the ratio
        # estimator below is unaffected and only the closure denominator
        # cares.
        files = job_files(f"{BASE}/powheg", at_energy('v2mu_job'),
                          "events.hepmc")
        sigma_sample_pb = powheg_v2_mu_sigma_pb()
        n_gen_tot = None                            # weighted: use sum of weights
        label = "POWHEG-V2 NLO + Pythia8"
        if ptag in HARD_TAGS:
            hard_flavour, tag_stats = HARD_TAGS[ptag], {}
            label += f" ({ptag})"
    else:
        sys.exit(f"unknown generator {gen}")
    if not files:
        sys.exit(f"no input files found for {gen}")

    run_stats = {}
    subtract = subtracts_shower_charm(gen)
    if subtract and tag_stats is None and hard_flavour is None:
        tag_stats = {}          # the raw final-state count travels with it
    obs, w, n_events, sum_w_all = analyze(files, beam_pid=beam_pid,
                                          run_stats=run_stats,
                                          lep_out=lep_out,
                                          hard_flavour=hard_flavour,
                                          tag_stats=tag_stats,
                                          hard_tag_mode=hard_tag_mode,
                                          subtract_shower_charm=subtract)
    print(f"{gen}: parsed {n_events} events, {len(w)} in fiducial region")

    # The tag suffix was stripped so the sample branches could match; put it
    # back before anything is written, or a "_charmfinal" run would overwrite
    # the "_charm" result file.
    if hard_tag_mode != "out":
        gen = gen + hard_tag_mode

    out = {"generator": gen, "label": label, "n_parsed": n_events,
           # ONE NAME FOR ONE FIELD (2026-08-29).  This was `hard_tag_mode`
           # here and `charm_tag_mode` in analyze_nu -- the same quantity under
           # two names, so a cross-current audit read one and concluded the
           # other side recorded nothing.  Both are written now: `tag_mode` is
           # the name to use, the other two are kept so already-stored results
           # and anything reading them keep working.
           "tag_mode": hard_tag_mode,
           "hard_tag_mode": hard_tag_mode,
           "charm_tag_mode": hard_tag_mode,
           "n_fiducial": len(w)}
    if gen.startswith("v2_"):
        # the nucleon is STAMPED, not only named: combine_target.py refuses
        # a p/n pair whose stamps do not say p and n
        out["target"] = sample_layout.require()[1]
    # Which fiducial region and which beam produced this.  Stamped so a plot
    # cannot silently combine two selections or two energies -- the same reason
    # input_manifest() records the files.  runmeta.stamp() is the one
    # definition, shared with analyze_nu, me_nlo_samples, the yadism_*
    # calculators and the Herwig histogrammers.
    stamp(out)
    out["inputs"] = input_manifest(files, sig_pb if sig_pb else None,
                                   job_counts if job_counts else None,
                                   f" [{gen}]")
    if tag_stats is not None:
        tot = sum(tag_stats["in_flavours"].values())
        print(f"{gen}: hard-process incoming flavours (all parsed events): "
              + ", ".join(f"|{k}|: {v} ({100*v/tot:.2f}%)"
                          for k, v in sorted(tag_stats["in_flavours"].items(),
                                             key=lambda kv: int(kv[0]))))
        print(f"{gen}: {tag_stats['n_no_hard_record']} events with no hard-"
              f"process parton record; {tag_stats['n_tagged']} tagged "
              f"|PDG| = {hard_flavour}, of which "
              f"{tag_stats['n_in_out_mismatch']} WITHOUT the same flavour "
              f"outgoing")
        print(f"{gen}: tag comparison on the same sample -- hard-out "
              f"{tag_stats.get('n_tagged_out', 0)}, final-state "
              f"{tag_stats.get('n_tagged_final', 0)} of {n_events} parsed; "
              f"{tag_stats.get('n_charm_with_bottom', 0)} events vetoed as "
              f"charm-from-bottom")
        for ex in tag_stats["mismatch_examples"]:
            print(f"    mismatch example: {ex}")
        out["tag_stats"] = tag_stats
        out["hard_flavour"] = hard_flavour
    if (pbase.startswith("pythia") or pbase.startswith("genie")
            or pbase.startswith("v2_genie") or pbase == "mg5_ps"):
        # unweighted: dsigma = (N_bin / N_generated) * sigma_generated / width
        #
        # mg5_ps IS NAMED HERE, and it has to be: the test below is on the
        # sample name, so a new unweighted sample that matches neither prefix
        # falls through to the WEIGHTED branch at the bottom, which pins
        # sigma_fid to the integrator and therefore reports the generated
        # cross-section however many events the selection actually keeps.
        # That is what it did on the first MG5 + Pythia run -- 38.076 nb with
        # 2% of the events outside the region, quoted to the digit.
        norm = sigma_sample_pb / n_gen_tot
        out["sigma_fid_pb"] = norm * len(w)
        out["sigma_fid_err_pb"] = norm * math.sqrt(len(w))
    elif (pbase in ("powheg", "powheg_v2", "powheg_hw", "v2_powheg")
          or qed_arm_of(pbase)
          or (pbase.startswith("v2_powheg_")
              and pbase[len("v2_powheg_"):] in sample_layout.ARMS)):
        # weighted, generation region WIDER than fiducial: sigma_fid is the
        # integrator sigma times the fiducial weight fraction
        norm = sigma_sample_pb / sum_w_all
        out["sigma_fid_pb"] = norm * float(w.sum())
        out["sigma_fid_err_pb"] = norm * float(np.sqrt((w * w).sum()))
        out["mean_weight_pb"] = sum_w_all / n_events
        if pbase in ("powheg_hw", "v2_powheg_hw"):
            # THE HERWIG ARM CANNOT USE THE MEAN-WEIGHT CLOSURE.  ThePEG
            # writes HepMC weights normalised to the largest weight in the
            # file, so they carry no absolute scale: the mean of them is a
            # pure number and comparing it to a cross-section in picobarns
            # reports "-100%", which says nothing about the sample.  Herwig's
            # own end-of-run table quotes the cross-section it delivered, and
            # that IS an independent measurement of the pinned normalisation
            # -- see herwig7/herwig_xsec.py:read_lhe_job.
            dlv, n_w, n_a = herwig_lhe_delivered_pb(files)
            out["n_lhe_offered"] = n_a
            out["herwig_shower_survival"] = (n_w / n_a) if n_a else None
            closure_check(gen, dlv, sigma_sample_pb, w, out,
                          "Herwig's own delivered sigma over the Les Houches "
                          "events consumed")
        else:
            # a QED arm reads the SAME LHE directory as `powheg`, by
            # construction
            lhe_key = "powheg" if qed_arm_of(pbase) else pbase
            if pbase.startswith("v2_powheg_"):
                # a arm or its baseline: the subset's own offered events
                n_read = sample_layout.powheg_offered("mu", sample_layout.require()[1],
                                                 ENERGY, files)
            elif qed_arm_of(pbase):
                # >>> BUT ONLY THE SEEDS IT ACTUALLY SHOWERED. <<<  The QED
                # arms were run on twelve seeds; the muon POWHEG-RES sample
                # has since grown to a hundred and four, and counting every
                # LHE in the directory made the denominator eight times too
                # big -- a closure of -87.6% on a sample that had not moved.
                # Found on 2026-09-08 when the arms were first analysed under
                # a selection since the sample grew.  The job directory names
                # carry the seed (qedfsr_job_7 <- pwgevents-0007.lhe), so the
                # seeds are read off the inputs rather than assumed.
                n_read = lhe_events_offered_for(
                    POWHEG_LHE_DIR[lhe_key](), files)
            else:
                n_read = lhe_events_offered(POWHEG_LHE_DIR[lhe_key]())
            out["n_lhe_offered"] = n_read
            dlv, how = delivered_sigma_pb(run_stats, n_events, n_read)
            closure_check(gen, dlv, sigma_sample_pb, w, out, how)
    else:
        # weighted, generation cuts == fiducial cuts: shape from weights,
        # normalisation pinned to the integrator's sigma.  MANDATORY: check
        # that the delivered events actually reproduce that sigma before
        # trusting it (see closure_check).
        norm = sigma_sample_pb / sum_w_all
        sw, sw2 = float(w.sum()), float((w * w).sum())
        # A SELECTION NARROWER THAN THE GENERATION CUTS breaks the pinning in
        # exactly the same way a channel tag does.  "Generation cuts ==
        # fiducial cuts" is true only for the INCLUSIVE selection, which is
        # what these samples were generated for; a FASER tier keeps a subset,
        # so the integrator total is no longer the answer for it.
        #
        # Pinning it anyway is not a small error, it is a silent one: Sherpa's
        # Tier S came out at 38.097 nb, its inclusive value to the digit, while
        # Pythia moved from 37.5 to 111.2.  That is CONVENTIONS.md rule 2's original
        # failure mode returning in a new place.
        # A RAISED Q2 FLOOR IS A SUBSET TOO, and it does not change the
        # selection's NAME -- so testing the name alone let a Q2 > 11 run
        # report the Q2 > 4 integrator cross-section to the digit, with the
        # event count correctly down to 39%.  Found on the first such run
        # (2026-08-31); it is the Tier S pinning bug in a new place, which is
        # the third time this exact shape has appeared.
        subset = (hard_flavour is not None
                  or SELECTION.name != selection.DEFAULT
                  or SELECTION.q2_min != selection.Q2_FLOOR_DEFAULT)
        if subset:
            # sum_w_all still counts every parsed event (see analyze()'s
            # docstring), so the absolute contribution of the surviving subset
            # is the integrator sigma times the surviving weight FRACTION --
            # the same recipe as the wider-generation POWHEG branch.
            out["sigma_fid_pb"] = norm * sw
            out["sigma_fid_err_pb"] = norm * math.sqrt(sw2)
            out["selected_weight_fraction"] = sw / sum_w_all
            if hard_flavour is not None:
                out["tagged_weight_fraction"] = sw / sum_w_all
        else:
            out["sigma_fid_pb"] = sigma_sample_pb
            # the quoted error is the MC error on the weighted sum, not
            # sigma/sqrt(N): the events carry unequal weights, so N
            # unit-weight events is the wrong statistics.
            out["sigma_fid_err_pb"] = sigma_sample_pb * math.sqrt(sw2) / sw
        out["n_eff"] = sw * sw / sw2 if sw2 else 0.0
        out["mean_weight_pb"] = sum_w_all / n_events
        dlv, how = delivered_sigma_pb(run_stats, n_events)
        closure_check(gen, dlv, sigma_sample_pb, w, out, how)
    # AUTOSCALE THE UNIT.  The dimuon tier's cross-sections are ~1e-4 nb, which
    # a fixed ".3f nb" renders as "0.000 nb" -- indistinguishable on sight from
    # the genuinely empty selection printed a few lines below, and Pythia's
    # 400 GeV dimuon result (10 selected events from 2,000,000) did exactly
    # that.  A result that small is a physics statement about the selection,
    # not a null, so it has to be legible as a number.
    print(f"{gen}: sigma_fid = {fmt_sigma(out['sigma_fid_pb'])}"
          f"  ({out.get('n_fiducial', 0)} events selected)")
    flag_low_stats(out, gen)

    # An EMPTY selection is a legitimate result, not a crash.  GENIE's EMDIS
    # list has no charm matrix element at all (there is only DIS-CC-CHARM), so
    # "genie_charmfinal" selects zero events out of 200k -- a measured zero the
    # charm table should be able to show.  Every weighted mean below divides by
    # sum(w), so write the result and stop rather than raising ZeroDivisionError.
    if len(w) == 0 or not np.any(w):
        out["empty_selection"] = True
        out["means"] = {}
        out["hists"] = {}
        os.makedirs(f"{BASE}/results", exist_ok=True)
        ofn = result_path(gen)
        with open(ofn, "w") as f:
            json.dump(stamp_shower_charm(out, tag_stats, w, subtract), f)
        print(f"{gen}: EMPTY SELECTION -- 0 events pass; wrote {ofn} "
              f"with sigma_fid = 0 and no histograms")
        return

    out["means"] = {kk: float(np.average(np.asarray(obs[kk], dtype=float),
                                         weights=w))
                    for kk in ("nch", "nch1", "nch05", "Q2", "xbj", "y")}
    dsum = dphi_summary(obs, w)
    out.update(dsum)
    out["means"]["dphi"] = dsum["mean_dphi"]
    print(f"{gen}: <dphi_min> = {dsum['mean_dphi']:.4f} rad, "
          f"{dsum['n_no_hadron']} events ({100*dsum['frac_no_hadron']:.3f}% "
          f"of sigma) with no charged hadron -> DeltaPhi unfilled")
    msum = dmeson_summary(obs, w)
    out.update(msum)
    out["means"].update({"nd_ch": msum.get("mean_n_d_ch"),
                         "nd_0": msum.get("mean_n_d_0"),
                         "nd_s": msum.get("mean_n_d_s")})
    print(f"{gen}: <N(D+-)> = {msum.get('mean_n_d_ch', 0):.4f}, "
          f"<N(D0)> = {msum.get('mean_n_d_0', 0):.4f}, "
          f"<N(Ds)> = {msum.get('mean_n_d_s', 0):.4f}; "
          f"sigma fraction with >=1: D+- {100*msum.get('frac_ge1_d_ch', 0):.2f}%"
          f", D0 {100*msum.get('frac_ge1_d_0', 0):.2f}%")

    hists = {}
    for kk, edges in BINS.items():
        v = np.asarray(obs[kk], dtype=float)
        cnt, _ = np.histogram(v, bins=edges, weights=w)
        err2, _ = np.histogram(v, bins=edges, weights=w*w)
        widths = np.diff(edges)
        hists[kk] = {"edges": edges.tolist(),
                     "dsig": (cnt*norm/widths).tolist(),
                     "err": (np.sqrt(err2)*norm/widths).tolist()}
    if not msum.get("d_supported"):
        for kk in D_HIST_KEYS:
            hists.pop(kk, None)
        print(f"{gen}: NO D mesons in the sample -- D histograms omitted "
              f"(the plots skip this generator rather than showing a zero)")
    out["hists"] = hists

    os.makedirs(f"{BASE}/results", exist_ok=True)
    ofn = result_path(gen)
    with open(ofn, "w") as f:
        json.dump(stamp_shower_charm(out, tag_stats, w, subtract), f)
    print(f"wrote {ofn}")


if __name__ == "__main__":
    main()
