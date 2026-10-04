#!/usr/bin/env python3
"""Write the GENIE config overlay for the LO neutrino CC DIS run.

Two files, both REGENERATED from the installed GENIE config on every run so
they cannot drift from it: an EventGeneratorListAssembler.xml carrying one
added generator list, and a CommonParam.xml with one raised parameter.

1. THE GENERATOR LIST: CCDISCHARM = DIS-CC + DIS-CC-CHARM

WHY THIS EXISTS.  GENIE's stock EventGeneratorListAssembler.xml has no
param_set that runs the classic Bodek-Yang DIS-CC generator TOGETHER with
DIS-CC-CHARM, and the two are not optional halves of each other:
QPMDISPXSec::XSec ends with

    // Subtract the inclusive charm production cross section
    xsec = TMath::Max(0., xsec - xsec_charm);

so `--event-generator-list CCDIS` alone delivers CC DIS with the charm piece
REMOVED -- a silently charm-free sample, and an inclusive cross-section short
by that amount.  The stock `CC` list has both, but drags in QEL/RES/COH/MEC/DFR
as well, which is not the DIS benchmark process.

2. THE VALIDITY CEILING: GVLD-Emax

GENIE's stock CommonParam.xml sets GVLD-Emax = 1000 GeV, the declared upper
end of the validity range of the classic (non-HEDIS) generators.  The
benchmark beam is EXACTLY 1 TeV, so gmkspl truncates every spline at the beam
energy -- and a GENIE spline evaluated AT its last knot returns zero, because
Spline::Evaluate asks for the knot above and FindClosestKnot reads past the
array (genie/README.md records the same trap on the muon side).  The result is
not an error: gevgen finds a total cross-section of zero, cannot select an
interaction, and spins forever printing "Could not select interaction".  So the
ceiling is raised to 5 TeV here, exactly as genie/config/nnpdf/CommonParam.xml
already does for the NNPDF muon run.

That ceiling is GENIE's own statement that this model is not validated at
1 TeV.  Raising it does not make it validated -- GENIE enters this benchmark
as-is (genie/README.md) and the LO CC row is a stress test of the Bodek-Yang
model far outside its tuning region, not a closure test.

Usage: make_nucc_lo_overlay.py [outdir]  (default: <repo>/genie/config/nu_cc_lo)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis"))
import paths  # noqa: E402  -- config.sh is the one source of every path

# $GENIE still wins when genie_setup.sh has been sourced; otherwise config.sh
# decides, rather than a second copy of the default living here.
GENIE = os.environ.get("GENIE") or paths.GENIE_DIR
STOCK = os.path.join(GENIE, "config", "EventGeneratorListAssembler.xml")
COMMON = os.path.join(GENIE, "config", "CommonParam.xml")

NEW_SET = """
  <!-- ADDED by evgen-benchmark/genie/make_nucc_lo_overlay.py.  Do not edit here:
       this file is regenerated from the installed GENIE config every run.
       Classic (Bodek-Yang, GRV98LO) charged-current DIS, complete.  The
       inclusive generator has its charm piece subtracted internally, so the
       charm generator must run alongside it or the sample has no charm and
       the rate is short.
       NB: a double hyphen is illegal inside an XML comment and libxml2
       rejects the WHOLE file, which GENIE reports only as "No Configuration
       available for .../CCDISCHARM". -->
  <param_set name="CCDISCHARM">
     <param type="int" name="NGenerators">    2                                 </param>
     <param type="alg" name="Generator-0 ">   genie::EventGenerator/DIS-CC       </param>
     <param type="alg" name="Generator-1 ">   genie::EventGenerator/DIS-CC-CHARM </param>
  </param_set>

</alg_conf>"""


GVLD_EMAX_GEV = 5000.0


def write_generator_list(outdir):
    with open(STOCK) as f:
        text = f.read()
    if "CCDISCHARM" in text:
        sys.exit(f"{STOCK} already defines CCDISCHARM -- upstream changed, "
                 f"drop this overlay instead of shadowing it.")
    tail = "</alg_conf>"
    if text.count(tail) != 1:
        sys.exit(f"unexpected structure in {STOCK}: {text.count(tail)} "
                 f"closing tags")
    out = os.path.join(outdir, "EventGeneratorListAssembler.xml")
    with open(out, "w") as f:
        f.write(text.replace(tail, NEW_SET.lstrip("\n")))
    print(f"wrote {out} (from {STOCK})")


def write_common_param(outdir):
    with open(COMMON) as f:
        text = f.read()
    pat = re.compile(r'(<param\s+type="double"\s+name="GVLD-Emax">\s*)'
                     r'([0-9.eE+-]+)(\s*</param>)')
    m = pat.search(text)
    if not m:
        sys.exit(f"no GVLD-Emax parameter in {COMMON} -- upstream changed, "
                 f"re-check the validity-ceiling trap before running.")
    if float(m.group(2)) > GVLD_EMAX_GEV:
        sys.exit(f"{COMMON} already allows {m.group(2)} GeV -- drop this "
                 f"overlay instead of lowering the ceiling.")
    out = os.path.join(outdir, "CommonParam.xml")
    with open(out, "w") as f:
        f.write(pat.sub(rf"\g<1>{GVLD_EMAX_GEV:.3f}\g<3>", text, count=1))
    print(f"wrote {out} (from {COMMON}, "
          f"GVLD-Emax {m.group(2)} -> {GVLD_EMAX_GEV:.3f} GeV)")


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "config", "nu_cc_lo")
    os.makedirs(outdir, exist_ok=True)
    write_generator_list(outdir)
    write_common_param(outdir)


if __name__ == "__main__":
    main()
