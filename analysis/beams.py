#!/usr/bin/env python3
"""Beam kinematics for the benchmark, at any lepton energy.

ONE PLACE where "1 TeV mu- on a proton at rest" turns into the numbers each
generator wants.  Every generator states the beams differently, and the
conversions are easy to get subtly wrong -- which is exactly what happened
(see CHECK below).  Nothing here changes an existing card; `beams.py check`
compares the committed cards against the formula and reports.

WHAT EACH GENERATOR WANTS
  Herwig / Sherpa   the two beam energies in the CENTRE-OF-MASS frame,
                    (E_A, E_B), which are unequal because both beams are
                    massive.
  POWHEG-RES        the LAB energies with `fixed_target 1`: ebeam1 = E,
                    ebeam2 = m_p.
  POWHEG-V2         a symmetric pair, ebeam1 = ebeam2 = sqrt(s)/2.
  MG5_aMC           the c.m. pair, like Herwig/Sherpa.
  GENIE             the LAB energy, `gevgen -e E`.
  YADISM            only ever needs 2 k.P = 2 E m_p.

Usage:
  beams.py                 show every configured energy
  beams.py check           compare the committed cards against the formula
  beams.py energies        the configured energies, space separated
  beams.py tag mu 400      the short tag for one beam        -> 400GeV
  beams.py name job 400    the per-energy directory name     -> job_400GeV
"""
import math
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Masses: the same values analyze.py uses, so the analysis and the generation
# cannot drift apart.
M_P = 0.938272
M_MU = 0.10566

# The lepton of each current, and its mass.  A neutrino is exactly massless
# here; that is not an approximation but the benchmark's convention.
LEPTONS = {"mu": ("mu-", M_MU), "nu": ("nu_mu", 0.0)}


class Beams:
    """Beam kinematics for one lepton of one energy on a proton at rest."""

    def __init__(self, current, energy_gev):
        if current not in LEPTONS:
            raise ValueError(f"current must be one of {list(LEPTONS)}")
        self.current = current
        self.name, self.m_lep = LEPTONS[current]
        self.e_lab = float(energy_gev)

    # --- invariants ---------------------------------------------------------
    @property
    def two_kP(self):
        """2 k.P, the only combination the analytic reference needs."""
        return 2.0 * self.e_lab * M_P

    @property
    def s(self):
        return self.m_lep**2 + M_P**2 + self.two_kP

    @property
    def sqrt_s(self):
        return math.sqrt(self.s)

    # --- per-generator forms -----------------------------------------------
    @property
    def cm_energies(self):
        """(E_lepton, E_proton) in the c.m. frame -- Herwig, Sherpa, MG5."""
        rs = self.sqrt_s
        return ((self.s + self.m_lep**2 - M_P**2) / (2.0 * rs),
                (self.s + M_P**2 - self.m_lep**2) / (2.0 * rs))

    @property
    def powheg_v2_ebeam(self):
        """The symmetric pair POWHEG-V2 wants: sqrt(s)/2."""
        return self.sqrt_s / 2.0

    @property
    def tag(self):
        """Short label for filenames and result keys: 500GeV, 1TeV, 2TeV."""
        e = self.e_lab
        if e >= 1000.0 and abs(e / 1000.0 - round(e / 1000.0)) < 1e-9:
            return f"{round(e / 1000.0):g}TeV"
        return f"{e:g}GeV"

    def summary(self):
        a, b = self.cm_energies
        return (f"{self.name} {self.e_lab:g} GeV on p at rest  [{self.tag}]\n"
                f"    s          = {self.s:.5f} GeV^2\n"
                f"    sqrt(s)    = {self.sqrt_s:.6f} GeV\n"
                f"    2 k.P      = {self.two_kP:.4f} GeV^2\n"
                f"    c.m. beams = ({a:.6f}, {b:.6f}) GeV\n"
                f"    POWHEG-V2  = {self.powheg_v2_ebeam:.11f} GeV each")


# The energies the benchmark is configured for (user-set 2026-08-24).  1 TeV is
# the published anchor that every existing sample sits at; 400 GeV and 4 TeV
# bracket it by a factor 2.5 either way, spanning the FASER flux.
ENERGIES = (400.0, 1000.0, 4000.0)

# THE PAPER-PLOTS-V2 LADDER (user, 2026-09-13): the final region on tungsten
# adds 700 GeV and 2 TeV.  A SEPARATE tuple, not a wider ENERGIES: 22 scripts
# loop over ENERGIES and would start asking the earlier samples for energies they
# were never generated at.  Differential distributions are shown at 1 TeV
# only, so the other four points carry the integrated rate and need less
# statistics (PRODUCTION.md).
BENCH_ENERGIES = (400.0, 700.0, 1000.0, 2000.0, 4000.0)

# THE SIDIS LADDER (user, 2026-09-18: the pion study moved to the region):
# the energies plus 300 GeV.  The pion yields at FASER are the Tier E rate
# folded over the flux, and Tier E asks for a scattered lepton above 200 GeV,
# so its efficiency falls from 0.15 at 400 GeV to 0.04 at 300 and to EXACTLY
# zero at and below TIER_E_ZERO_ENERGY -- while a third of the flux-weighted
# rate of the region sits below 400 GeV.  Holding the 400 GeV efficiency
# constant below it, as the earlier study did, overestimates the yield by several
# per cent; a measured 300 GeV point and the exact zero fix that.  300 GeV is
# a SIDIS-ladder point only: not in BENCH_ENERGIES, so no cross-section figure
# and no FASER rate outside the pion study reads it.
SIDIS_ENERGIES = (300.0,) + BENCH_ENERGIES
TIER_E_ZERO_ENERGY = 200.0      # E' > 200 GeV cannot pass at E_beam <= 200 GeV

# THE SIDIS LADDER BELOW 300 GeV, POWHEG-V2 ONLY (2026-09-21).  The Tier E
# yield does not need it (Tier E is exactly zero at <= 200 GeV), but the
# REGION-ONLY yield does -- Q2 > 4 and W > 3 alone, the rate a collaborator's
# analytic NLO calculation is benchmarked against, and the denominator of the
# Tier E efficiency factors -- and a quarter of the nu_mu region rate sits
# below 300 GeV, where it used to be a log-log extrapolation of the YADISM
# shape with the z spectra frozen at 300 GeV.  POWHEG-V2 ONLY: it is the
# benchmark's reference for the region-only yield, and the other generators
# keep the extrapolation.  Muon side NOT extended: the SIDIS study draws the
# neutrino current alone (CONVENTIONS.md 2b carve-out).  A SEPARATE tuple so the
# card makers of Herwig, Sherpa and GENIE, which loop over
# SIDIS_ENERGIES, do not start producing cards for points never run.
SIDIS_ENERGIES_LOW = (20.0, 50.0, 100.0, 200.0)

# The anchor: samples at this energy predate the scan and keep their original
# untagged paths and result keys, so nothing published moves.
ANCHOR_ENERGY = 1000.0


def configured():
    return [Beams(c, e) for e in ENERGIES for c in ("mu", "nu")]


def at_energy(name, energy):
    """Job-dir / run-dir / result name for one energy.

    THE ANCHOR KEEPS ITS ORIGINAL UNTAGGED NAME, so every sample, card and
    result that existed before the scan stays exactly where it is and nothing
    published moves; other energies get a `_<tag>` suffix.

    This lives here, not in analyze.py or in the run scripts, because all
    three must agree: a job written to a tagged directory and read back from
    an untagged glob would silently analyse an empty set, or worse, mix two
    energies into one histogram.
    """
    return (name if float(energy) == ANCHOR_ENERGY
            else f"{name}_{Beams('mu', energy).tag}")


# --------------------------------------------------------------- card check
# (file, regex capturing the numbers, what they mean)
CARD_CHECKS = [
    ("herwig7/DIS-mu-POWHEG.in", "mu", 1000.0, "cm",
     r"BeamEMaxA\s+([0-9.]+)\*GeV.*?BeamEMaxB\s+([0-9.]+)\*GeV"),
    ("herwig7/DIS-nu-POWHEG.in", "nu", 1000.0, "cm",
     r"BeamEMaxA\s+([0-9.]+)\*GeV.*?BeamEMaxB\s+([0-9.]+)\*GeV"),
    ("sherpa/Runs/MuonDIS_NLO/Sherpa.yaml", "mu", 1000.0, "cm",
     r"BEAM_ENERGIES:\s*\[\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\]"),
    ("sherpa/Runs/NuDIS_NLO/Sherpa.yaml", "nu", 1000.0, "cm",
     r"BEAM_ENERGIES:\s*\[\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\]"),
    ("powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input", "nu", 1000.0, "v2",
     r"ebeam1\s+([0-9.]+)d0.*?ebeam2\s+([0-9.]+)d0"),
    ("powheg/cards/POWHEG-RES/mu1TeV-wide/powheg.input", "mu", 1000.0, "lab",
     r"ebeam1\s+([0-9.]+)d0.*?ebeam2\s+([0-9.]+)d0"),
]

TOL = 1e-4          # relative; the cards quote 6 significant figures


def check():
    bad = 0
    for rel, current, energy, kind, pat in CARD_CHECKS:
        path = os.path.join(BASE, rel)
        if not os.path.exists(path):
            print(f"  ?  {rel}: not found")
            continue
        m = re.search(pat, open(path).read(), re.S)
        if not m:
            print(f"  ?  {rel}: no beam-energy line matched")
            continue
        got = (float(m.group(1)), float(m.group(2)))
        b = Beams(current, energy)
        want = {"cm": b.cm_energies,
                "v2": (b.powheg_v2_ebeam, b.powheg_v2_ebeam),
                "lab": (b.e_lab, M_P)}[kind]
        dev = max(abs(g - w) / w for g, w in zip(got, want) if w)
        ok = dev < TOL
        bad += not ok
        print(f"  {'OK ' if ok else 'BAD'} {rel}\n"
              f"        card    ({got[0]:.6f}, {got[1]:.6f})\n"
              f"        formula ({want[0]:.6f}, {want[1]:.6f})   "
              f"max deviation {dev:.3%}")
    print(f"\n{len(CARD_CHECKS) - bad}/{len(CARD_CHECKS)} cards match the "
          f"formula at {TOL:.0e} relative.")
    if bad:
        print("A mismatch is NOT auto-corrected: the committed cards are what\n"
              "produced the published samples, so changing them silently would\n"
              "invalidate every number on the page. Decide, then regenerate.")
        print(KNOWN_MISMATCH)
    return bad


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        sys.exit(1 if check() else 0)
    # Shell entry points, so run scripts share this module's naming rules
    # instead of reimplementing them:
    #   beams.py tag  mu 400    -> 400GeV
    #   beams.py name job 400   -> job_400GeV     (job at the anchor)
    #   beams.py energies       -> 400 1000 4000
    if len(sys.argv) == 4 and sys.argv[1] == "tag":
        print(Beams(sys.argv[2], float(sys.argv[3])).tag)
        return
    if len(sys.argv) == 4 and sys.argv[1] == "name":
        print(at_energy(sys.argv[2], float(sys.argv[3])))
        return
    if len(sys.argv) == 2 and sys.argv[1] == "energies":
        print(" ".join(f"{e:g}" for e in ENERGIES))
        return
    for b in configured():
        print(b.summary())
        print()


if __name__ == "__main__":
    main()


KNOWN_MISMATCH = """
KNOWN, as of 2026-08-21, and MEASURED rather than estimated:

  * The MUON cards (Herwig DIS-mu-*.in, Sherpa MuonDIS_*/Sherpa.yaml) use
    c.m. beam energies 21.6438 / 21.6641, i.e. s = 1875.574 and
    2 k.P = 1874.68 -- as if the muon were at 998.94 GeV rather than 1000.
    The muon YADISM reference uses 2 k.P = 1876.544 exactly, so the muon
    generator samples and the muon analytic reference sit at slightly
    different energies.  Worth -0.0187% on sigma_fid, measured by running
    yadism_calc.manual_lo_sigma at both values (38.11128 -> 38.10415 nb).
    That is an order of magnitude below the tightest closure quoted on the
    page (0.9996), so it changes no conclusion -- but it is a real
    inconsistency and every muon number carries it.

  * POWHEG-RES uses m_p = 0.938 rather than 0.938272 (0.029% in s).
    NOT AN ERROR, and do not "fix" it (established 2026-08-25).  POWHEG
    writes ebeam2 into the LHE <init> line as beam 2's ENERGY, and Pythia's
    proton mass is 0.93827.  Below that, Pythia reads the target as AT REST --
    the same trick as `Beams:eB = 0.` in the Pythia card.  Raised to 0.938272
    it sits just ABOVE m_p, Pythia tries to build a moving proton, and every
    event dies with "ProcessContainer::constructProcess: setting mass failed"
    (12000 aborts, 0 events per job).  The 0.029% is the price of a
    fixed target that Pythia will actually shower.

  * The NEUTRINO cards are exact on all three generators.

Nothing has been regenerated on account of this: the 1 TeV samples cost a
great deal of effort and the bias is negligible.  New energies take correct
values from this module, which does leave 1 TeV very slightly inconsistent
with the rest of any future scan -- that is the trade, recorded here so the
choice is visible rather than inherited.
"""
