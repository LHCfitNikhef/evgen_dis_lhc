#!/usr/bin/env python3
"""Refuse all-caps emphasis in the report's prose.

WHY.  CONVENTIONS.md rule 5: report prose is written in sentence case, and emphasis
is <b>, not capitals.  "THE TWO CURRENTS ARE LIMITED BY DIFFERENT THINGS" and
"the SAME-ORDER YADISM reference" both shipped before anyone noticed, because
capitals read as deliberate rather than as a mistake -- nothing errors, the
sentence is still true, and it is only jarring once someone reads the page as
a reader rather than as its author.

WHAT IT DOES NOT FLAG.  Acronyms and names are not shouting: YADISM, FONLL,
NLO, CKM, GENIE, ZM-VFNS, HEDIS keep their capitals, and so does anything in
ALLOWED below.

SINGLE WORDS ARE FLAGGED TOO, since 2026-08-27.  They did not used to be, on
the theory that too many single capitalised words are legitimate -- and the
report duly shipped "the two studies above are ANALYTIC", "requires BOTH muons
above 100 GeV", "the error on the SHIFT is what matters" and six more, none of
which the two-word rule could see.  The fix is not a cleverer rule but a real
word list: a single all-caps word is shouting unless it is a known acronym or
name.  That means ALLOWED has to be kept up to date, and a NEW acronym will
raise a false positive once, which is the right way round -- a false positive
costs one line in this file, a false negative ships.

WHERE IT LOOKS.  The RENDERED page, not the source.  The emphatic capitals in
CONVENTIONS.md, in code comments and in commit messages are for a different
audience and are deliberate; only what a reader of the report sees is in
scope, and rendering is the only way to tell those apart reliably.

Usage: tools/check_report_caps.py [path/to/report.html]
Exit status 1 if the report shouts.
"""
import collections
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT = f"{BASE}/results/report.html"

# capitals that are names, units or acronyms rather than emphasis
ALLOWED = {
    # the neutrino-telescope report (results_telescopes/, 2026-10-01): an
    # IceCube event selection, a KM3NeT detector, the evolution equations
    "ESTES", "ORCA", "ARCA", "DGLAP", "HKKM2014",
    # a file format, on the SIDIS tab's pointer to share/faser-sidis/
    "CSV",
    # fragmentation-function families, experiments and the process names of
    # the "Fragmentation functions" tab (2026-09-21)
    "SIA", "JAM", "JAM20", "HERMES", "COMPASS", "HAPS", "MAP", "KKP",
    "NPC23",
    # FASER's own report numbers, cited on the "Comparison with FASER"
    # sub-tab: CERN-EP-2024-079 and CERN-FASER-CONF-2026-002.
    "CERN", "CONF", "CERN-FASER-CONF", "CERN-EP", "FASER-CONF",
    # HL-LHC, in the Run 3 + Run 4 projection note; and README, which names
    # data/faser_emulsion/README.md, the file that records the digitisation.
    "HL", "LHC", "HL-LHC", "README",
    # NEUT is a neutrino generator (T2K/Super-K); "III" numbers a table of
    # arXiv:2412.03186 that the electronic-detector sub-tab cites.
    "NEUT", "III",
    # ZM-VFN is the spelling the paper plots use (user, 2026-08-30); ZM-VFNS
    # is the spelling the rest of the report uses.  Both are the same scheme
    # and both are acronyms, so both are allowed -- and "VFN" appears on its
    # own because the checker splits on the hyphen.
    "ZM-VFN", "VFN", "VFNS",
    # experiments whose data the "Comparison with data" tab uses
    "NOMAD", "HERA", "NMC", "EMC", "CHORUS", "NuTeV", "CCFR",
    # the simulation the FASER muon flux is made with (2026-09-02)
    "FLUKA",
    # the Forward Physics Facility, named on the SIDIS comparison (2026-09-07)
    "FPF",
    # the forward hadronic-interaction models, named where the neutrino FLUX
    # is discussed and nowhere else (2026-09-04)
    "SIBYLL", "DPMJET", "EPOS", "QGSJET",
    # GENIE's process names, on the "Impact of non-DIS processes" tab
    # (2026-09-04): quasi-elastic, diffractive, meson-exchange, coherent,
    # and the electromagnetic current
    "QEL", "DFR", "MEC", "COH", "EM",
    "YADISM", "GENIE", "POWHEG", "FONLL", "HEDIS", "EMDIS", "AGKY", "AHADIC",
    "NLO", "LO", "NNLO", "NNPDF", "GRV98", "CT18", "MSHT20", "ATLAS", "BGR18",
    "DIS", "CC", "NC", "MC", "PDF", "CKM", "QED", "QCD", "EW", "IR", "UV",
    "ZM", "VFNS", "FFNS", "FFN0", "RES", "V2", "LHE", "HTML", "JSON", "CL",
    "TeV", "GeV", "MeV", "FASER", "SM", "ME", "PS", "II", "I", "A", "B",
    "S", "E", "D", "K", "Q", "W", "Z", "X", "Y", "N",
    # names and acronyms that appear as LONE all-caps words in the prose
    "BYPDF", "APFEL", "PDG", "LHAPDF", "MHOU", "FSR", "ISR", "AGKY",
    "AHADIC", "EMDIS", "HEDIS", "SCALUP", "MPI", "QCD", "QED", "CTEQ",
    "ABMP", "HERA", "SIDIS", "VFNS", "FFNS", "POWHEG", "YADISM", "GENIE",
    "FONLL", "NNPDF", "MSHT", "GRV", "BGR", "FASER", "NLO", "NNLO", "CKM",
    "HTML", "JSON", "LHE", "TODO",
    # Sherpa's Catani-Seymour shower module, and the accelerator
    "CSS", "LHC", "DIM", "SHERPA", "ABMP16", "ABMP",
    # DONUT is the Fermilab tau-neutrino experiment (E872) and LEPTO its
    # own event generator, both names (2026-09-10)
    "DONUT", "LEPTO",
    # final-state interactions, an acronym the DONUT figure's legend uses and
    # the caption beside it defines (2026-09-10)
    "FSI",
}
RUN = re.compile(r"\b([A-Z][A-Z0-9]{1,}(?:[ -][A-Z][A-Z0-9]{1,})+)\b")
# a lone all-caps word of three or more letters.  Two-letter words are left
# alone: too many are units, symbols or initials, and the run rule already
# covers them when they appear together.
LONE = re.compile(r"\b([A-Z]{3,})\b")


def prose(html):
    """The text a reader actually sees."""
    h = re.sub(r"data:image/[^\"']+", "", html)
    h = re.sub(r"<script.*?</script>|<style.*?</style>", " ", h, flags=re.S)
    h = re.sub(r"<[^>]+>", " ", h)
    for ent, ch in (("&mdash;", "-"), ("&ndash;", "-"), ("&nbsp;", " "),
                    ("&amp;", "&"), ("&plusmn;", "+-"), ("&times;", "x")):
        h = h.replace(ent, ch)
    return h


def shouted(text):
    out = []
    for m in RUN.finditer(text):
        words = re.split(r"[ -]", m.group(1))
        if all(w in ALLOWED for w in words):
            continue
        out.append(m.group(1))
    # single words, checked against the word list rather than against a shape
    for m in LONE.finditer(text):
        if m.group(1) not in ALLOWED:
            out.append(m.group(1))
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(path):
        sys.exit(f"no report at {path} -- build it with "
                 f"analysis/make_report.py first")
    hits = shouted(prose(open(path, encoding="utf-8", errors="replace").read()))
    if not hits:
        print("report prose: no all-caps emphasis")
        return 0
    print(f"*** {len(hits)} all-caps phrase(s) in the report prose ***",
          file=sys.stderr)
    for w, n in collections.Counter(hits).most_common():
        print(f"  x{n:<4} {w}", file=sys.stderr)
    print("\nReport prose is sentence case; emphasis is <b>, not capitals "
          "(CONVENTIONS.md rule 5).\nIf one of these is a name or an acronym, add "
          "it to ALLOWED in this file.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
