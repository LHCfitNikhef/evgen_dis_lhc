#!/usr/bin/env python3
"""Emit the Herwig 7 POWHEG-NLO cards of the production (PRODUCTION.md).

    herwig7/production/make_cards.py            # write herwig7/production/cards/*.in
    herwig7/production/make_cards.py --check    # report what would change

The earlier production cards DIS-{mu,nu}-POWHEG.in / -PWGNEG.in are the templates
and are NEVER touched (earlier must stay reproducible).  Every card is one of
them with exactly six things changed, each by an ASSERTED substitution (a
pattern that is not there aborts the script -- Herwig trap 8, a card that
silently stayed LO):

  1. BEAM ENERGIES, from analysis/beams.py.  The earlier muon card sits at
     E = 998.94 GeV (beams.KNOWN_MISMATCH); is exact at every energy.
  2. TARGET.  p: unchanged.  n: a GENUINE NEUTRON BEAM, /Herwig/Particles/n0,
     with the PROTON set NNPDF40_nnlo_as_01180 -- see NEUTRON below.
  3. CUTS, the final region: Q2 > 4, W > 3, NO y cut -- see CUTS below.
  4. saverun name, so no the earlier production .run file can be clobbered.
  5. EventHandler:MaxEventErrorRetries -- see RETRIES below.
  6. RemnantDecayer:DISRemnantOption NoLepton -- see RETRIES below.

Everything else -- POWHEG ME, Contribution 1/2, the polarised neutrino beam,
EW scheme, QED splittings removed, ctau > 10 mm, IntrinsicPtGaussian 0.5 GeV,
HepMC3 output, MaxErrors -- is inherited byte for byte.  `diff` any card
against its template to see it.

NEUTRON (settled 2026-09-13 from the source and then by a run)
==============================================================
ThePEG::LHAPDF (ThePEG-2.3.0/PDF/LHAPDF6.cc, xfx/xfvx/xfsx) ISOSPIN-CONJUGATES
THE SET ITSELF when the beam is a neutron: for particle id 2112 a request for
u returns the set's d, ubar returns dbar, and so on; canHandleParticle accepts
|id| = 2212 or 2112.  So a neutron beam must be given the PROTON set.  Giving
it NNPDF40_nnlo_as_01180_n swaps twice and silently generates a PROTON.  The
remaining pieces accept a neutron as they are: n0 is a ThePEG::BeamParticleData
in Herwig's defaults (baryons.in), HwRemDecayer builds the valence content from
the PDG code (2112 -> d d u, so the remnant is a neutron's), and DISBase and
the shower's Sudakov both evaluate the PDF through the beam particle, so the
same isospin swap reaches the NLO weight and the backward evolution.
The validation arms below measure it (herwig7/production/validate_neutron.sh).

The neutron beam is set up with the SAME invariant k.P = E_lab * M_P as the
proton one, M_P = beams.M_P.  That is the benchmark's per-nucleon convention:
the analytic reference (YADISM on the isoscalar set) and target.per_nucleon
both take 2k.P = 2 E M_P for either nucleon, and an analysis that recovers
E_lab = p.P / M_P from the event's target record gets exactly E back.  The
c.m. energies are then solved with the neutron's own mass (Herwig's 0.9395654
GeV), so ThePEG does not have to reshuffle the beams on shell.  Physically it
is a neutron at rest hit by E * M_P/m_n, 0.14% below E -- stated, not hidden.

CUTS
====
The fiducial region is applied EXACTLY in the analysis (selection q4w3,
W2 = M_P^2 + y 2kP - Q2 > 9, from the FINAL-STATE lepton).  At generation
Herwig is cut as tightly as it allows WITHOUT EVER BEING TIGHTER:

  LeptonicCut                  MinQ2 3.5, MaxQ2 1e6, Miny 0, Maxy 1, MinW2 5
  {Neutral,Charged}CurrentCut  MinQ2 3.5, MaxQ2 1e6, Miny 0, Maxy 1

so sigma_fid = sigma_gen x (fiducial fraction) -- the cuts do NOT coincide
and an analysis must not assume they do.

Maxy 1 on ThePEG's SimpleDISCut is REQUIRED, not cosmetic: left at the
template's 0.9 it silently re-imposes the y window under the new name.

WHY LOOSER THAN THE REGION, ON BOTH VARIABLES (measured 2026-09-13).

1. The cut's own W2.  LeptonicDISCut builds W2 = y_lep * SMax - Q2 with
   y_lep = (q0 + qz)/(2 E_lepton) in the partonic frame, dropping the hadron
   mass.  A 4M-point scan of the exact ThePEG kinematics per energy, current
   and nucleon (massive muon and nucleon, parton = x times the hadron's
   light-cone momentum, SMax = 4 E_A E_B) gives, for Q2 > 4 and analysis
   W2 > 9: min of the cut's W2 = 8.126-8.159 GeV2 in all 20 cases (the
   hadron mass, 0.88 GeV2); y_lep in (0, 1) and SimpleDISCut's y < 1
   everywhere.  So at MATRIX-ELEMENT level 8 would do, by 0.13 GeV2.

2. Herwig MOVES THE SCATTERED LEPTON.  "Warning had to adjust the momentum of
   the non-colour connected final-state, e.g. the scattered lepton in DIS":
   6% of muon and 2.5% of neutrino events leave the shower with a lepton
   different from the hard process's, Q2 shifted by up to +-20% and W2 raised
   by up to 30 GeV2.  Measured on 50k-event proton samples generated LOOSE
   (Q2 > 2.8, MinW2 5) at 1 TeV and 400 GeV, comparing the hard-process lepton
   (status 11, parent of the final one) with the final lepton:

     of events in the 4 < Q2 < 5 bin    generated at Q2 <= 4   <= 3.5   <= 3
         mu 1 TeV / 400 GeV                  0.97 / 1.06%     0.15 / 0.28%
         nu 1 TeV / 400 GeV                  1.11 / 1.13%     0.32 / 0.38%
     of events in the 9 < W2 < 12 bin   generated at W2_hard <= 8.38 (the most
         a MinW2 7.5 can remove): 2.8 / 3.6% (mu), 3.7 / 1.2% (nu);
         generated at W2_hard <= 5.88 (the most MinW2 5 can remove): 0 / 0.

   A Q2 cut coincident with the region (as earlier had) therefore leaves a ~1%
   hole in the first Q2 bin, and MinW2 7.5 a ~3% hole in the first W bin --
   silently, since the in-flow is simply never generated.  MinQ2 3.5 cuts the
   Q2 hole to <= 0.4% of the edge bin (0.04-0.08% of sigma_fid, below the edge
   bin's own statistical error) at the price of 11-12% extra muon events
   (1-2% neutrino); 3.0 would cost 27-29%.  MinW2 5 closes the W hole and
   costs 0.5% of events.  NNPDF4.0's Q0^2 = 2.72 GeV2 is below 3.5, so muF = Q
   stays on the grid.

RETRIES (2026-09-19)
====================
Herwig cannot put the beam remnant on shell after the forced g -> q qbar
splitting of a SEA quark at high x and low W (HwRemDecayer::setRemMasses).
Unpatched ThePEG then throws the whole event away and samples a NEW point,
while the analysis normalises to the attempted cross-section -- so the failed
points' rate is handed to the ones that survive.  The failures are charm, s
and sea at high x: at 400 GeV on a proton the delivered hard-charm fraction in
Q2 > 4, W > 5 was 14.45% (nu) / 4.49% (mu) against 18.49% / 5.29% at matrix-
element level, and the charm figure read 0.74-0.98 of the reference.
patches/thepeg-2.3.0-event-error-retry.diff adds
EventHandler:MaxEventErrorRetries: the event is rebuilt on the SAME point up
to that many times before it is dropped.  With 100: 18.63% / 5.60%, 96-99% of
failures recovered, survival 99.4-99.9%.  A card carrying the switch cannot be
read by an unpatched ThePEG (no such interface), which is the point: the
fingerprint changes and no pre-fix sample can be mistaken for a fixed one.

NoLepton (same day, user-approved).  Where the remnant cannot take the
recoil, Herwig's default (DISRemnantOption Default) MOVES THE SCATTERED
LEPTON instead (the "had to adjust the momentum of the non-colour connected
final-state" warning; 6% of muon events, see CUTS).  Every observable here
is built from that lepton, and for high-x charm the move is large: with
retries alone the muon charm x spectrum (1 TeV p, pos - neg) was
1.37/1.58/1.25 x the matrix-element level at x = 0.09-0.19 and 0.5-0.8
above 0.24 -- a migration, not a loss.  NoLepton makes such an event fail
instead, and the retry rebuilds it; the bump goes (0.97/1.01/1.11).  What
never succeeds without moving the lepton is dropped after the retries (0.3%
of positive- and 1.5% of negative-half muon events, all at high x, ~2% of
the charm in the region), the smaller evil for a lepton-defined benchmark.

Card inventory (herwig7/production/cards/):
  V2-<cur>-PWG[NEG]-<t>[-<TAG>].in     production, t = p | n, TAG absent at 1 TeV
  V2VAL-<cur>-PWG-<arm>.in            1 TeV neutron validation arms:
      pn   proton beam  + NNPDF40_nnlo_as_01180_n  (the fallback route)
      nn   neutron beam + NNPDF40_nnlo_as_01180_n  (the double swap: must
           reproduce the PROTON, which is the proof ThePEG swaps)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "analysis"))

import beams   # noqa: E402  -- the one source of beam kinematics
import target  # noqa: E402  -- the one source of the PDF set names

CHECK = "--check" in sys.argv
CARDDIR = os.path.join(HERE, "cards")

# Herwig 7.3.0's own masses (src/defaults/baryons.in, leptons.in).  Used only
# to put each beam ON SHELL in ThePEG's eyes; the invariant k.P is beams.M_P.
HW_MASS = {"p+": 0.93827208816, "n0": 0.9395654205}
HW_MMU = 0.1056583755

MIN_Q2 = 3.5          # GeV2 -- see CUTS in the docstring
MIN_W2 = 5.0          # GeV2 -- see CUTS in the docstring
RETRIES = 100         # see RETRIES in the docstring
DIS_REMNANT = "NoLepton"   # never move the scattered lepton; see RETRIES
ENERGIES = beams.SIDIS_ENERGIES     #  + the 300 GeV SIDIS-ladder point
HALVES = {"pos": ("DIS-{cur}-POWHEG.in", "DIS-{cur}-PWG", "PWG"),
          "neg": ("DIS-{cur}-PWGNEG.in", "DIS-{cur}-PWGNEG", "PWGNEG")}


def cm_energies(cur, energy, beam_b):
    """(E_lepton, E_nucleon) in the c.m. frame at k.P = E * beams.M_P."""
    if beam_b == "p+":
        # the proton: exactly beams.py's numbers, so `beams.py check` style
        # comparisons hold for as they do for the earlier production
        return beams.Beams(cur, energy).cm_energies
    m_a = HW_MMU if cur == "mu" else 0.0
    m_b = HW_MASS[beam_b]
    s = m_a**2 + m_b**2 + 2.0 * energy * beams.M_P
    rs = math.sqrt(s)
    return ((s + m_a**2 - m_b**2) / (2.0 * rs),
            (s + m_b**2 - m_a**2) / (2.0 * rs))


def sub(text, old, new, what):
    """Replace exactly one occurrence; abort if the template has moved."""
    n = text.count(old)
    assert n == 1, f"{what}: expected 1 occurrence of {old!r}, found {n}"
    return text.replace(old, new)


FASER_ENERGIES = (30.0, 60.0, 100.0, 200.0, 300.0, 400.0, 700.0, 1000.0,
                  2000.0, 4000.0, 6800.0)


def faser_stem(cur, half, tgt, energy):
    """FL-<cur>-PWG[NEG]-<t>-<tag>: the FASER-ladder cards, always tagged."""
    return f"FL-{cur}-{HALVES[half][2]}-{tgt}-{beams.Beams('nu', energy).tag}"


def stem(cur, half, tgt, energy):
    s = f"V2-{cur}-{HALVES[half][2]}-{tgt}"
    return s if energy == beams.ANCHOR_ENERGY else \
        f"{s}-{beams.Beams(cur, energy).tag}"


def make_card(cur, half, energy, beam_b, pdfset, runname, label,
              region="q4w3"):
    """region "q4w3": the production (Q2 > 3.5, MinW2 5).  region "faser":
    the FASER ladder of paper plots 12/12b and the antineutrino SIDIS
    (2026-10-04): Q2 > 4 AT GENERATION, as the POWHEG-V2 and Sherpa FASER
    ladders have it, and no W cut beyond MinW2 0.3 (the cut's partonic W2).
    cur "nubar" is the neutrino template with the beam turned round."""
    beam_cur = "nu" if cur == "nubar" else cur
    min_q2, min_w2 = (MIN_Q2, MIN_W2) if region == "q4w3" else (4.0, 0.3)
    src_name, save_stem, _ = (x.format(cur=beam_cur) for x in HALVES[half])
    text = open(os.path.join(REPO, "herwig7", src_name)).read()
    simple = "NeutralCurrentCut" if cur == "mu" else "ChargedCurrentCut"
    if cur == "nubar":
        # the template defines nu_mubar_L (polarisation +1) already
        text = sub(text, "set EventHandler:BeamA /Herwig/Particles/nu_mu\n",
                   "set EventHandler:BeamA /Herwig/Particles/nu_mubar\n", "beam A")
        text = sub(text, "set EventHandler:BeamA /Herwig/Particles/nu_mu_L\n",
                   "set EventHandler:BeamA /Herwig/Particles/nu_mubar_L\n",
                   "polarised beam A")

    # 4. saverun
    text = sub(text, f"saverun {save_stem} EventGenerator",
               f"saverun {runname} EventGenerator", "saverun")

    # 1. beam energies (the only two BeamEMax lines in the card)
    ea, eb = cm_energies(beam_cur, energy, beam_b)
    n = 0
    out = []
    for line in text.splitlines(keepends=True):
        ls = line.strip()
        if ls.startswith("set Luminosity:BeamEMaxA"):
            line = f"set Luminosity:BeamEMaxA {ea:.6f}*GeV\n"
            n += 1
        elif ls.startswith("set Luminosity:BeamEMaxB"):
            line = f"set Luminosity:BeamEMaxB {eb:.6f}*GeV\n"
            n += 1
        out.append(line)
    assert n == 2, f"{src_name}: expected 2 BeamEMax lines, replaced {n}"
    text = "".join(out)

    # 2. target and PDF
    if beam_b == "n0":
        text = sub(text, "set EventHandler:BeamB /Herwig/Particles/p+\n",
                   "set EventHandler:BeamB /Herwig/Particles/n0\n", "beam B")
        note = ("# v2: the NEUTRON beam carries the PROTON set; ThePEG::LHAPDF\n"
                "# swaps u<->d for id 2112 itself (LHAPDF6.cc::xfx)\n"
                if pdfset == target.PDFSET_P else
                "# VALIDATION: neutron beam AND neutron set = a double isospin\n"
                "# swap, expected to reproduce the PROTON\n")
        text = sub(text, "set /Herwig/Particles/p+:PDF NNPDF40\n",
                   "set /Herwig/Particles/p+:PDF NNPDF40\n" + note +
                   "set /Herwig/Particles/n0:PDF NNPDF40\n", "n0 PDF")
    if pdfset != target.PDFSET_P:
        text = sub(text, f"set NNPDF40:PDFName {target.PDFSET_P}\n",
                   f"set NNPDF40:PDFName {pdfset}\n", "PDFName")

    # 3. cuts: y window out, W2 in, Q2 floor lowered, on BOTH cut objects
    text = sub(text, "set /Herwig/Cuts/LeptonicCut:MinQ2 4.0*GeV2\n",
               f"set /Herwig/Cuts/LeptonicCut:MinQ2 {min_q2}*GeV2\n",
               "LeptonicCut MinQ2")
    text = sub(text, f"set /Herwig/Cuts/{simple}:MinQ2 4.0*GeV2\n",
               f"set /Herwig/Cuts/{simple}:MinQ2 {min_q2}*GeV2\n",
               f"{simple} MinQ2")
    text = sub(text, "set /Herwig/Cuts/LeptonicCut:Miny 0.2\n",
               "set /Herwig/Cuts/LeptonicCut:Miny 0.0\n", "LeptonicCut Miny")
    text = sub(text, "set /Herwig/Cuts/LeptonicCut:Maxy 0.9\n",
               "set /Herwig/Cuts/LeptonicCut:Maxy 1.0\n"
               f"set /Herwig/Cuts/LeptonicCut:MinW2 {min_w2}*GeV2\n",
               "LeptonicCut Maxy")
    text = sub(text, f"set /Herwig/Cuts/{simple}:Miny 0.2\n",
               f"set /Herwig/Cuts/{simple}:Miny 0.0\n", f"{simple} Miny")
    text = sub(text, f"set /Herwig/Cuts/{simple}:Maxy 0.9\n",
               f"set /Herwig/Cuts/{simple}:Maxy 1.0\n", f"{simple} Maxy")

    # 5. rebuild an event whose remnant fails on the SAME point (RETRIES)
    text = sub(text, f"saverun {runname} EventGenerator",
               "# v2: a failed remnant/shower is rebuilt on the SAME phase-space\n"
               "# point, never replaced by a new one (make_cards.py RETRIES)\n"
               "set /Herwig/EventHandlers/EventHandler:MaxEventErrorRetries "
               f"{RETRIES}\n"
               "# v2: never move the scattered lepton to rescue the remnant; the\n"
               "# event fails and is retried instead (make_cards.py RETRIES)\n"
               "set /Herwig/Partons/RemnantDecayer:DISRemnantOption "
               f"{DIS_REMNANT}\n"
               f"saverun {runname} EventGenerator", "retries")

    banner = (
        f"# GENERATED by herwig7/production/make_cards.py from herwig7/{src_name}"
        f" -- edit the generator, not this file.\n"
        f"# v2 production (PRODUCTION.md): {label}\n"
        f"#   E_lab = {energy:g} GeV, k.P = E_lab * {beams.M_P} GeV2, "
        f"beams (c.m.) = [{ea:.6f}, {eb:.6f}] GeV\n"
        f"#   target beam {beam_b}, PDF {pdfset}\n"
        f"#   {'region faser; ' if region == 'faser' else ''}generation cuts: Q2 > {min_q2} GeV2, NO y cut (Miny 0 / Maxy 1 "
        f"on BOTH cut objects),\n"
        f"#   LeptonicCut MinW2 {min_w2} GeV2 -- deliberately LOOSER than the "
        f"region, see\n"
        f"#   make_cards.py.  The exact Q2 > 4, W > 3 GeV is applied in the "
        f"analysis.\n"
        f"#   Comments below that\n"
        f"#   mention 0.2 < y < 0.9, p+ or 1 TeV describe the template.\n")
    return banner + text


def write(path, text):
    rel = os.path.relpath(path, REPO)
    old = open(path).read() if os.path.exists(path) else None
    if old == text:
        print(f"  ok        {rel}")
        return 0
    print(f"  {'would write' if CHECK else ('rewrite  ' if old else 'create   ')} {rel}")
    if not CHECK:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(text)
    return 1


def main():
    changed = 0
    for energy in ENERGIES:
        for cur in ("mu", "nu"):
            for tgt, beam_b in (("p", "p+"), ("n", "n0")):
                for half in ("pos", "neg"):
                    rn = stem(cur, half, tgt, energy)
                    changed += write(
                        os.path.join(CARDDIR, rn + ".in"),
                        make_card(cur, half, energy, beam_b, target.PDFSET_P,
                                  rn, f"{cur} {half} target {tgt}"))
    # THE FASER LADDER (2026-10-04): nu and nubar on p and n at the energies
    # of the emulsion/electronic folds and of the SIDIS ladder
    for energy in FASER_ENERGIES:
        for cur in ("nu", "nubar"):
            for tgt, beam_b in (("p", "p+"), ("n", "n0")):
                for half in ("pos", "neg"):
                    rn = faser_stem(cur, half, tgt, energy)
                    changed += write(
                        os.path.join(CARDDIR, rn + ".in"),
                        make_card(cur, half, energy, beam_b, target.PDFSET_P,
                                  rn, f"FASER ladder {cur} {half} target {tgt}",
                                  region="faser"))
    for cur in ("mu", "nu"):
        for arm, beam_b in (("pn", "p+"), ("nn", "n0")):
            rn = f"V2VAL-{cur}-PWG-{arm}"
            changed += write(
                os.path.join(CARDDIR, rn + ".in"),
                make_card(cur, "pos", beams.ANCHOR_ENERGY, beam_b,
                          target.PDFSET_N, rn,
                          f"VALIDATION ONLY, {cur} pos, beam {beam_b} with "
                          f"the neutron set"))
    print(f"\n{changed} card(s) {'would change' if CHECK else 'written'}.")


if __name__ == "__main__":
    main()
