#!/usr/bin/env python3
"""Write the GENIE config overlay for the non-DIS study's dedicated runs.

Two files, both REGENERATED from the installed GENIE config on every run so
they cannot drift from it (the same pattern, and the same CommonParam
function, as make_nucc_lo_overlay.py):

1. EventGeneratorListAssembler.xml with two added lists, each the stock
   default list MINUS its DIS generators:

     CCNONDIS = QEL-CC + RES-CC + DFR-CC + QEL-CC-CHARM + QEL-CC-LAMBDA
     EMNONDIS = QEL-EM + RES-EM

   WHY.  The non-DIS channels are 0.1-1% of the rate at FASER energies
   (analysis/genie_nondis.py), so a fully inclusive run of a million events
   holds a few thousand of them.  A run with ONLY those channels gives their
   kinematic shapes with full statistics; each shape is then normalised to
   the channel's own spline cross-section, which is what gevgen would have
   done anyway.  The inclusive run (stock lists CC and EM) is still made and
   is the closure: its per-channel event fractions must reproduce the
   spline fractions.

   The nuclear-only channels (COH, MEC) are left out: on a free nucleon they
   produce no spline (genie/splines/faser/*.xml has none) and a list that
   names them would make gevgen look for one.

2. CommonParam.xml with GVLD-Emax raised to 5 TeV, so that gevgen at the
   1 TeV benchmark point is not AT the declared ceiling (the spline trap of
   make_nucc_lo_overlay.py).

Usage: make_nondis_overlay.py [outdir]  (default: <repo>/genie/config/nondis)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_nucc_lo_overlay as base  # noqa: E402

NEW_SETS = """
  <!-- ADDED by evgen-benchmark/genie/make_nondis_overlay.py.  Do not edit
       here: this file is regenerated from the installed GENIE config every
       run.  The default lists WITHOUT their DIS generators, for the non-DIS
       study's high-statistics shape runs.  (No double hyphen in this
       comment: libxml2 would reject the whole file.) -->
  <param_set name="CCNONDIS">
     <param type="int" name="NGenerators">    5                                  </param>
     <param type="alg" name="Generator-0 ">   genie::EventGenerator/QEL-CC        </param>
     <param type="alg" name="Generator-1 ">   genie::EventGenerator/RES-CC        </param>
     <param type="alg" name="Generator-2 ">   genie::EventGenerator/DFR-CC        </param>
     <param type="alg" name="Generator-3 ">   genie::EventGenerator/QEL-CC-CHARM  </param>
     <param type="alg" name="Generator-4 ">   genie::EventGenerator/QEL-CC-LAMBDA </param>
  </param_set>
  <param_set name="EMNONDIS">
     <param type="int" name="NGenerators">    2                                  </param>
     <param type="alg" name="Generator-0 ">   genie::EventGenerator/QEL-EM        </param>
     <param type="alg" name="Generator-1 ">   genie::EventGenerator/RES-EM        </param>
  </param_set>

</alg_conf>"""


def write_generator_list(outdir):
    with open(base.STOCK) as f:
        text = f.read()
    for name in ("CCNONDIS", "EMNONDIS"):
        if name in text:
            sys.exit(f"{base.STOCK} already defines {name} -- upstream "
                     f"changed, drop this overlay instead of shadowing it.")
    tail = "</alg_conf>"
    if text.count(tail) != 1:
        sys.exit(f"unexpected structure in {base.STOCK}")
    out = os.path.join(outdir, "EventGeneratorListAssembler.xml")
    with open(out, "w") as f:
        f.write(text.replace(tail, NEW_SETS.lstrip("\n")))
    print(f"wrote {out} (from {base.STOCK})")


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "config", "nondis")
    os.makedirs(outdir, exist_ok=True)
    write_generator_list(outdir)
    base.write_common_param(outdir)


if __name__ == "__main__":
    main()
