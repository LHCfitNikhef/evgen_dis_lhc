#!/usr/bin/env python3
"""Boost a POWHEG DIS Les Houches file from the fixed-target lab frame into the
lepton-proton centre-of-mass frame, so that Herwig 7 can read it.

WHY THIS EXISTS, and it is not a matter of taste.  POWHEG runs these samples
with `fixed_target 1`: the Les Houches init line declares beam B as a proton of
energy 0.938 GeV, which is its MASS, i.e. a proton AT REST WITH ZERO MOMENTUM.
Pythia does not care -- it reads the two incoming legs out of each event.
ThePEG builds parton bins from momentum fractions along the beam direction, and
a beam with no momentum makes that a division by zero: the run dies at
`saverun` with

    While writing object ... of class ThePEG::Cuts:
    Tried to write a NaN or Inf double to a persistent stream.

AND IT WRITES A .run FILE ANYWAY, which then fails to read back with a message
about a corrupted file that names neither the beam nor the cause.  Two hours of
this is why the paragraph is long.

WHAT THE BOOST IS.  A single boost along z into the frame where the incoming
lepton and the proton have equal and opposite momenta.  It is exact and it
changes no physics: every invariant -- Q2, y, x, the event weight, SCALUP --
is untouched, and `analysis/analyze.py` reconstructs the lab frame covariantly
from the beams in the event record, which is how the Herwig and Sherpa samples
(generated in this very frame) have always been analysed.  The boosted file is
therefore the SAME EVENTS, which is the whole point of the comparison it
serves: one matrix element, two showers and two hadronisation models.

THE NEW BEAM ENERGIES ARE NOT sqrt(s)/2.  They are written so that each beam
carries exactly the common momentum p*, with its own mass: E_A = sqrt(p*^2 +
m_lepton^2) and E_B = sqrt(p*^2 + m_proton^2).  Declaring sqrt(s)/2 for both
instead would put the light-cone fraction of the incoming lepton a part in 10^5
ABOVE ONE, and x > 1 is rejected downstream rather than warned about.

Usage:  tools/lhe_boost_cm.py <in.lhe> <out.lhe>
"""
import math
import os
import sys

M_P = 0.9382720813          # GeV, proton
LEPTONS = {11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16}


def beta_of(e_l, pz_l, m_target):
    """Boost velocity taking (lepton, target at rest) to their c.m. frame."""
    return pz_l / (e_l + m_target)


def boost_z(p, beta):
    """(E, px, py, pz) boosted along z by beta."""
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    e, px, py, pz = p
    return (g * (e - beta * pz), px, py, g * (pz - beta * e))


def read_init(line):
    f = line.split()
    return [int(f[0]), int(f[1]), float(f[2]), float(f[3])] + f[4:]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    lines = open(src).readlines()

    # ---------------------------------------------------------- the header
    try:
        i_init = next(k for k, l in enumerate(lines)
                      if l.strip().startswith("<init>"))
    except StopIteration:
        sys.exit(f"{src}: no <init> block -- not a Les Houches file")
    init_fields = lines[i_init + 1].split()
    if len(init_fields) < 10:
        sys.exit(f"{src}: malformed init line")
    # NOT EVERY POWHEG DIS FILE NEEDS THIS.  POWHEG-RES runs the muon side
    # with `fixed_target 1` and writes the proton at rest; POWHEG-V2 writes
    # the neutrino side already in the lepton-proton c.m. frame, with both
    # beams at 21.66 GeV.  A file of the second kind is passed through as a
    # SYMLINK, so one call site serves both currents and the log says which
    # happened.  The test is on beam B's energy against the proton mass, and
    # it is exact enough: the two cases differ by a factor of twenty.
    if abs(float(init_fields[3]) - M_P) > 1e-3:
        if float(init_fields[3]) < 2.0 * M_P:
            sys.exit(f"{src}: beam B energy {init_fields[3]} GeV is neither "
                     f"the proton mass nor a collinear beam.  Refusing rather "
                     f"than guessing which frame this is.")
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(src), dst)
        print(f"{os.path.basename(src)}: beams already collinear at "
              f"{float(init_fields[2]):.4f} / {float(init_fields[3]):.4f} GeV "
              f"-- no boost needed, symlinked")
        return

    # ------------------------------------- the boost, from the first event
    beta = lep_ref = ebeam_a = ebeam_b = None
    for k, l in enumerate(lines):
        if not l.strip().startswith("<event>"):
            continue
        npart = int(lines[k + 1].split()[0])
        for j in range(npart):
            f = lines[k + 2 + j].split()
            if int(f[1]) == -1 and int(f[0]) in LEPTONS:
                e, pz, m_lep = float(f[9]), float(f[8]), float(f[10])
                beta = beta_of(e, pz, M_P)
                lep_ref = e
                pstar = abs(boost_z((e, 0.0, 0.0, pz), beta)[3])
                ebeam_a = math.hypot(pstar, m_lep)
                ebeam_b = math.hypot(pstar, M_P)
                break
        break
    if beta is None:
        sys.exit(f"{src}: no incoming lepton found in the first event")

    # ------------------------------------------------------------- rewrite
    n_ev = 0
    with open(dst, "w") as out:
        k, n = 0, len(lines)
        while k < n:
            l = lines[k]
            if k == i_init + 1:
                init_fields[2] = "%.10E" % ebeam_a
                init_fields[3] = "%.10E" % ebeam_b
                out.write(" " + " ".join(init_fields) + "\n")
                k += 1
                continue
            if not l.strip().startswith("<event>"):
                out.write(l)
                k += 1
                continue
            # an event: the tag, its header line, then exactly NUP particles.
            # A file cut short mid-event is a TRUNCATED SAMPLE, which is the
            # failure mode this whole repository is organised around -- it
            # parses, so it must be refused here by name.
            if k + 1 >= n:
                sys.exit(f"{src}: file ends inside an <event> tag -- the "
                         f"sample is truncated")
            out.write(l)
            out.write(lines[k + 1])
            npart = int(lines[k + 1].split()[0])
            if k + 2 + npart > n:
                sys.exit(f"{src}: event {n_ev + 1} declares {npart} particles "
                         f"but the file ends first -- the sample is truncated")
            for j in range(npart):
                f = lines[k + 2 + j].split()
                if len(f) < 13:
                    sys.exit(f"{src}: event {n_ev + 1}, particle {j + 1} has "
                             f"{len(f)} fields, not 13 -- truncated or not a "
                             f"Les Houches particle line")
                if int(f[1]) == -1 and int(f[0]) in LEPTONS:
                    if abs(float(f[9]) - lep_ref) > 1e-6 * lep_ref:
                        sys.exit(f"{src}: incoming lepton energy {f[9]} "
                                 f"differs from the first event's {lep_ref} "
                                 f"-- one sample must have one beam")
                b = boost_z((float(f[9]), float(f[6]), float(f[7]),
                             float(f[8])), beta)
                f[6] = "%.9E" % b[1]
                f[7] = "%.9E" % b[2]
                f[8] = "%.9E" % b[3]
                f[9] = "%.9E" % b[0]
                out.write(" " + " ".join(f) + "\n")
            k += 2 + npart
            n_ev += 1
    print(f"{os.path.basename(src)}: {n_ev} events boosted, "
          f"beta = {beta:.12f}, E_A = {ebeam_a:.6f}, E_B = {ebeam_b:.6f} GeV "
          f"(sqrt(s) = {ebeam_a + ebeam_b:.4f} GeV)")


if __name__ == "__main__":
    main()
