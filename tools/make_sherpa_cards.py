#!/usr/bin/env python3
"""Emit the Sherpa MC@NLO cards of the paper-plots production.

    tools/make_sherpa_cards.py            # write the cards
    tools/make_sherpa_cards.py --check    # report what would change

THE SPEC IS PRODUCTION.md: Q2 > 4 GeV2 and W > 3 GeV with NO y cut, on a
proton AND on a neutron, at beams.BENCH_ENERGIES, both currents (rule 2b).  One
card per (current, target, energy), 20 in all, written to
sherpa/Runs/V2_{MuonDIS,NuDIS}_NLO_<t>[_TAG]/Sherpa.yaml and linked into
$SHERPA_RUNS by tools/link_cards.sh / tools/sherpa_production.sh.  <t> is always
explicit; the 1 TeV name is untagged (beams.at_energy).

The templates are the earlier production cards, which are READ ONLY here:
    muon NC     sherpa/Runs/MuonDIS_NLO/Sherpa.yaml
    neutrino CC sherpa/Runs/NuDIS_NLO_ckm3/Sherpa.yaml   (full CKM, patched V)
Everything that is not a beam, a PDF or a selector is kept character for
character: RESPECT_MASSIVE_FLAG, KIN_SCHEME, the integration targets, the
massless-charm convention, CKM Order 3 on the neutrino side.

>>> THE NEUTRON IS A GENUINE NEUTRON BEAM (2112), RUN FROM A PATCHED COPY. <<<
(user decision 2026-09-13: "patched copy for neutron runs")
  * PDF/LHAPDF/LHAPDF_CPP_Interface.C applies NO isospin transformation (it
    only flips an antiparticle bunch), so a 2112 beam is given the neutron
    set NNPDF40_nnlo_as_01180_n.  The proton set on a 2112 beam would be a
    proton PDF with a neutron remnant.
  * SHERPA/Initialization/Initialization_Handler.C defaults PDF_LIBRARY only
    for a 2212 bunch; without it a 2112 beam stops LOUDLY ("PDF ... does not
    exist in any of the loaded libraries for n bunch").  So the neutron cards
    name it.
  * STOCK 3.0.5 SEGFAULTS on the first event of a 2112 beam:
    REMNANTS/Main/Hadron_Remnant.C, RemnantFlavour(), sets its `taken` flag
    after the first constituent whatever it is, which for u d d builds kf
    21101 (u struck) or 201 (d struck).  patches/sherpa-hadron-remnant-
    neutron.diff is the one-line fix; tools/build_sherpa_npatch.sh builds it
    into a separate install copy ($SHERPA_NPATCH, config.sh), and
    tools/sherpa_production.sh runs every neutron dir from that copy and checks, in
    each job, the library actually mapped and the remnant in every event.
  * Evidence (PRODUCTION.md, Status, Sherpa): at 1 TeV the 2112-beam and
    the 2212-beam + _n integrations agree to 3e-6; with the patch a 2112 beam
    gives valence (u,d) = (1,2) and diquarks 2101/1103 in every event, and
    proton events at a fixed seed are byte-identical with and without it.
  * The neutron keeps Sherpa's own mass (0.939566).  Its c.m. beam energies
    are chosen so that k.P = E M_P EXACTLY, the invariant the analysis
    (analyze.py: E_lab = p.P / M_P, W2 = M_P^2 + 2P.q - Q2) and the YADISM
    reference (2kP = 2 E M_P) are built on.  So a neutron event and a proton
    event at the same nominal E share every DIS invariant; only s differs, by
    m_n^2 - m_p^2 = 0.0024 GeV2.
NEUTRON_BEAM = False reverts to the tools/sherpa_faser_ladder.sh method (a
2212 beam carrying the neutron PDF: right cross-section, proton remnant).

>>> THE W CUT AT GENERATION IS AN INEL y FLOOR. <<<  Sherpa has no W
selector (PHASIC++/Selectors: DIS_Selector.C has INEL and Q2 only), so the
card carries y > (9 - M_P^2 + 4) / 2kP, the W = 3 edge at the Q2 floor, and
no upper bound (y <= 1).  Above the floor W > 3 needs a larger y, so the
floor is looser than the region everywhere else; the analysis applies the
exact W cut.
  NOTE Sherpa's INEL y is not the invariant P.q/P.k: IINEL_Selector::Trigger
  builds it as 1 - (E'/E)(1 + cos theta)/2 from the two lepton momenta, which
  is exact only for a massless target along -z.  With the massive target of
  the c.m. frame it comes out LARGER than the invariant y, by 4e-5 (4 TeV) to
  4e-4 (400 GeV) relative at the Q2 = 4 edge and by up to 1.2e-3 at
  Q2 = 500 (solved exactly for the edge kinematics, PRODUCTION.md).  Larger
  means the selector passes every event whose invariant y is above the floor,
  so the formula value is never tighter than the fiducial region.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "analysis"))

import beams   # noqa: E402

CHECK = "--check" in sys.argv

# Sherpa's own beam masses (BEAM/Main/Beam_Parameters.C, SHERPA/SoftPhysics/
# Hadron_Init.C).  The PROTON value is beams.M_P; the neutron's is Sherpa's.
M_N_SHERPA = 0.939566
M_MU_SHERPA = 0.105658      # Sherpa's default muon mass

TEMPLATES = {"mu": "sherpa/Runs/MuonDIS_NLO/Sherpa.yaml",
             "nu": "sherpa/Runs/NuDIS_NLO_ckm3/Sherpa.yaml"}
BEAM_PID = {"mu": 13, "nu": 14}
BASE = {"mu": "V2_MuonDIS_NLO", "nu": "V2_NuDIS_NLO"}
# See the docstring.  True since the patched copy exists (2026-09-13); the
# neutron dirs MUST then be run from $SHERPA_NPATCH (tools/sherpa_production.sh does).
NEUTRON_BEAM = True
PDF = {"p": "NNPDF40_nnlo_as_01180", "n": "NNPDF40_nnlo_as_01180_n"}

W2_MIN = 9.0
Q2_MIN = 4.0


def run_name(cur, tgt, energy):
    return beams.at_energy(f"{BASE[cur]}_{tgt}", energy)


def target_beam(tgt, neutron_beam):
    """(PDG code, Sherpa mass) of the target beam."""
    if tgt == "n" and neutron_beam:
        return 2112, M_N_SHERPA
    return 2212, beams.M_P


def cm_energies(cur, tgt, energy, neutron_beam=NEUTRON_BEAM):
    """c.m. beam energies with k.P = E * M_P exactly, whatever the target mass.

    For the proton this is beams.Beams.cm_energies.  (Sherpa treats the muon
    beam as MASSLESS -- its particle table has Massive 0 -- so the generated
    k.P sits 1.5e-5 (400 GeV) to 1.5e-6 (4 TeV) above E M_P; negligible.)
    """
    ml = M_MU_SHERPA if cur == "mu" else 0.0
    mt = target_beam(tgt, neutron_beam)[1]
    kP = energy * beams.M_P
    s = ml * ml + mt * mt + 2.0 * kP
    rs = math.sqrt(s)
    return (s + ml * ml - mt * mt) / (2.0 * rs), (s + mt * mt - ml * ml) / (2.0 * rs)


def y_floor(energy):
    """The W = 3 GeV edge at Q2 = 4, as an invariant y: W2 = M_P^2 + y 2kP - Q2.

    Rounded DOWN at 1e-9, so the value written to the card is never above
    the formula (a rounded-up floor would be tighter than the region).
    """
    y = (W2_MIN - beams.M_P ** 2 + Q2_MIN) / (2.0 * energy * beams.M_P)
    return math.floor(y * 1e9) / 1e9


def card(cur, tgt, energy, neutron_beam=NEUTRON_BEAM):
    t = open(os.path.join(REPO, TEMPLATES[cur])).read()
    pid, lep = BEAM_PID[cur], 13
    tpid, tmass = target_beam(tgt, neutron_beam)
    nbeam = tpid == 2112
    ea, eb = cm_energies(cur, tgt, energy, neutron_beam)
    ymin = y_floor(energy)
    pdf = PDF[tgt]

    # the template's own header describes the earlier region and a 1 TeV proton;
    # it is replaced, not kept, so no comment survives to describe a cut or a
    # target that is not in this card
    body = t[t.index("\nBEAMS:") + 1:]
    head = (
        f"# PAPER PLOTS v2 -- GENERATED by tools/make_sherpa_cards.py from\n"
        f"# {TEMPLATES[cur]}; edit the generator, not this file.\n"
        f"#\n"
        f"# {'Neutral-current muon' if cur == 'mu' else 'Charged-current neutrino'}"
        f" DIS at NLO QCD (MC@NLO), {energy:g} GeV "
        f"{'mu-' if cur == 'mu' else 'nu_mu'} on a\n"
        f"# {'PROTON' if tgt == 'p' else 'NEUTRON'} at rest"
        f" ({tpid} beam, {pdf}).\n"
        f"# Region (PRODUCTION.md): Q2 > 4 GeV2 and W > 3 GeV, NO y cut.\n"
        f"# Collider mode in the c.m. frame, beam energies chosen so that\n"
        f"# 2 k.P = 2 E M_P = {2 * energy * beams.M_P:.6f} GeV2 exactly"
        f" (target mass {tmass}).\n"
        f"# Generation cuts: Q2 > 4 and Sherpa INEL y > {ymin:.9f}, the W = 3\n"
        f"# edge at the Q2 floor; looser than W > 3 everywhere above it.  The\n"
        f"# analysis applies the exact W cut (selection q4w3).\n")
    if tgt == "n" and not nbeam:
        head += (
            "#\n"
            "# THE NEUTRON AS A PROTON BEAM CARRYING THE ISOSPIN-MIRROR PDF: the\n"
            "# cross-section is the neutron's, the beam REMNANT is a proton's.\n"
            "# (NEUTRON_BEAM = False in tools/make_sherpa_cards.py.)\n")
    if nbeam:
        head += (
            "#\n"
            "# A GENUINE NEUTRON BEAM.  Sherpa's LHAPDF interface does no isospin\n"
            "# swap, so the neutron PDF set is given explicitly; PDF_LIBRARY is\n"
            "# named because Sherpa only defaults it for a 2212 bunch, and the\n"
            "# remnant is built from the 2112 code (u d d).  STOCK SHERPA 3.0.5\n"
            "# SEGFAULTS ON THIS BEAM: run it only from the patched copy\n"
            "# $SHERPA_NPATCH (tools/build_sherpa_npatch.sh, tools/sherpa_production.sh).\n")
    head += "#\n# ---- below: the v1 template, beams/PDF/selectors replaced ----\n"

    body = re.sub(r"^BEAMS: .*$", f"BEAMS: [{pid}, {tpid}]", body,
                  count=1, flags=re.M)
    body = re.sub(r"^BEAM_ENERGIES: .*$",
                  f"BEAM_ENERGIES: [{ea:.7f}, {eb:.7f}]", body, count=1, flags=re.M)
    pdfline = f"PDF_SET: [None, {pdf}]"
    if nbeam:
        pdfline = "PDF_LIBRARY: [None, LHAPDFSherpa]\n" + pdfline
    body = re.sub(r"^PDF_SET: .*$", pdfline, body, count=1, flags=re.M)
    body = re.sub(r"^MPI_PDF_SET: .*$", f"MPI_PDF_SET: {pdf}", body,
                  count=1, flags=re.M)

    # the selector block is rewritten whole: the template's comments there
    # describe the y window this production removes
    i = body.index("\nSELECTORS:")
    body = body[:i + 1] + (
        "SELECTORS:\n"
        "# Q2 > 4 GeV2, and the W > 3 GeV edge at Q2 = 4 as a y floor (Sherpa\n"
        "# has no W selector).  No upper y bound.\n"
        f"- [Q2, {pid}, {lep}, 4, 1e12]\n"
        f"- [INEL, {pid}, {lep}, {ymin:.9f}, 1.0]\n")

    for must in (f"BEAMS: [{pid}, {tpid}]", f"PDF_SET: [None, {pdf}]",
                 f"MPI_PDF_SET: {pdf}", "RESPECT_MASSIVE_FLAG: true",
                 f"- {pid} 93 -> {lep} 93:", "NLO_Mode: MC@NLO"):
        assert must in body, f"{cur}/{tgt}/{energy}: card lost '{must}'"
    assert "0.2, 0.9" not in body
    if cur == "nu":
        assert "Order: 3" in body, "the neutrino card lost CKM Order 3"
    return head + body


# ------------------------------------------------------------ validation
# THE NEUTRON-BEAM CHECK (PRODUCTION.md, Status, Sherpa).  Two arms at
# 1 TeV that differ ONLY in the target beam, both with the neutron PDF:
#   V2check_NuDIS_NLO_nbeam       2112 beam  (what the production uses)
#   V2check_NuDIS_NLO_pbeam_npdf  2212 beam  (tools/sherpa_faser_ladder.sh)
# Q2 > 4 only, no INEL at all, i.e. exactly the ladder's region, so the
# ladder's own NuDIS_NLO_FASER_nu_n_E1000 integration is a third number.
# The two cross-sections must agree within integration errors; the showered
# events of the first must carry a neutron remnant (no uu diquark, dd present).
def check_card(beam):
    text = card("nu", "n", 1000.0, neutron_beam=(beam == "nbeam"))
    text = re.sub(r"^- \[INEL, .*\]\n", "", text, flags=re.M)
    text = re.sub(r"^# Region \(PRODUCTION.md\).*\n", "", text, flags=re.M)
    text = re.sub(r"^# Generation cuts: .*\n.*\n.*\n",
                  "# Generation cut: Q2 > 4 ONLY (validation arm, no INEL).\n",
                  text, flags=re.M)
    text = re.sub(r"^# Q2 > 4 GeV2, and the W > 3 GeV edge.*\n.*\n",
                  "# Q2 > 4 GeV2 only.\n", text, flags=re.M)
    text = ("# VALIDATION ARM, not production: the neutron-beam check, Q2 > 4 only\n"
            "# (no INEL), the region of tools/sherpa_faser_ladder.sh.\n" + text)
    assert "- [INEL" not in text and "NNPDF40_nnlo_as_01180_n" in text
    return text


def write(path, text):
    rel = os.path.relpath(path, REPO)
    old = open(path).read() if os.path.exists(path) else None
    if old == text:
        print(f"  ok       {rel}")
        return
    print(f"  {'would write' if CHECK else ('rewrite  ' if old else 'create   ')} {rel}")
    if not CHECK:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(text)


def main():
    for cur in ("mu", "nu"):
        for tgt in ("p", "n"):
            for e in beams.SIDIS_ENERGIES:   #  + the 300 GeV SIDIS point
                d = os.path.join(REPO, "sherpa", "Runs", run_name(cur, tgt, e))
                write(os.path.join(d, "Sherpa.yaml"), card(cur, tgt, e))
    for beam in ("nbeam", "pbeam_npdf"):
        d = os.path.join(REPO, "sherpa", "Runs", f"V2check_NuDIS_NLO_{beam}")
        write(os.path.join(d, "Sherpa.yaml"), check_card(beam))
    if "--table" in sys.argv:
        for e in beams.BENCH_ENERGIES:
            print(f"  E = {e:6g} GeV   2kP = {2*e*beams.M_P:10.4f}   "
                  f"INEL y_min = {y_floor(e):.9f}")


if __name__ == "__main__":
    main()
