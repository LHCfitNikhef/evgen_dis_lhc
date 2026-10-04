#!/usr/bin/env python3
"""Paper plot A1b: the non-DIS channels of GENIE's default tune at
E = 100 GeV, the companion of ppA1_genie_nondis.py (1 TeV).

WHY (user, 2026-09-30): "Can we maybe add the same plots for E_l = 100 GeV?
These non-DIS effects can be larger in this case?"  They are, in the total
cross-section: the explicit non-DIS channels carry a few per cent of the
neutrino rate at 100 GeV against 0.3-0.7% at the benchmark energies.  The
figure shows how much of that survives the benchmark region, Q2 > 4 GeV2 and
W > 3 GeV, with the same panels, the same drawing code and the same
conventions as ppA1 (DIS and DIS charm drawn as one, a measured zero drawn
as one); only the beam energy differs.

Inputs: results{,_nu}/genie_nondis_diff_W_100GeV.json
        (genie/run_genie_nondis.sh <cur> 100 1000000 200000 8 <p|n>, then
         analysis/genie_nondis_diff.py <cur> 100).
Usage: analysis/paper_plots/ppA1b_genie_nondis_100GeV.py
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ppA1_genie_nondis as a1                              # noqa: E402

SLUG = "ppA1b_genie_nondis_100GeV"
V1_NUMBER = None
IN_PAPER = True
TITLE = "Impact of non-DIS processes in GENIE at 100 GeV"
OUTPUT = "pp_genie_nondis_100GeV.png"
RESULTS = "results_nu"
ENERGY = 100.0

CAPTION = ("As the 1 TeV figure, for E<sub>&mu;</sub> = E<sub>&nu;</sub> = "
           "100 GeV: GENIE G18_02a on tungsten per nucleon, inside Q<sup>2</sup> "
           "&gt; 4 GeV<sup>2</sup> and W &gt; 3 GeV, with the non-DIS share "
           "per bin beneath each distribution.")

MESSAGE = """At 100 GeV the explicit non-DIS channels of GENIE's default tune carry 2.6% of the neutrino charged-current cross-section and 1.2% of the muon one, against 0.3-0.7% at 1 TeV. The benchmark region removes them all the same: resonances and quasi-elastic scattering stay below W = 3 GeV, so inside Q2 > 4 GeV2 and W > 3 GeV only neutrino diffractive scattering survives, at 0.008% of the fiducial rate (0.004% at 1 TeV), and the muon share is exactly zero. What grows at low energy is GENIE's DIS channel below the Q2 floor: 17% of it at 100 GeV, against 3-6% at the benchmark energies."""


def _load(cur):
    with open(a1._file(cur, ENERGY)) as f:
        return json.load(f)


def _share_fid(cur):
    return _load(cur)["regions"]["share_of_fiducial"]


def _share_tot(cur):
    s = _load(cur)["sigma_pb"]
    return 1.0 - (s["DIS"] + s["DIS charm"]) / s["total"]


CLAIMS = [
    {"what": "the non-DIS channels carry 2.6% (neutrino) and 1.2% (muon) of "
             "the TOTAL cross-section at 100 GeV, against 0.3-0.7% at 1 TeV",
     "check": lambda: (abs(100 * _share_tot("nu") - 2.59) < 0.01
                       and abs(100 * _share_tot("mu") - 1.22) < 0.01),
     "detail": lambda: f"nu {100*_share_tot('nu'):.3f}%, mu {100*_share_tot('mu'):.3f}%"},
    {"what": "inside the region only neutrino DFR survives, 0.008% of the "
             "fiducial rate (0.004% at 1 TeV); RES and QEL are zero on both "
             "currents and the muon non-DIS share is exactly zero",
     "check": lambda: (abs(100 * _share_fid("nu")["non-DIS"] - 0.0076) < 0.0005
                       and _share_fid("nu")["RES"] == 0.0 and _share_fid("nu")["QEL"] == 0.0
                       and _share_fid("mu")["non-DIS"] == 0.0),
     "detail": lambda: (f"nu DFR {100*_share_fid('nu')['DFR']:.4f}%, "
                        f"mu {100*_share_fid('mu')['non-DIS']:.4f}%")},
    {"what": "17% of GENIE's neutrino DIS channel lies below Q2 = 4 GeV2 at "
             "100 GeV (3-6% at 400 GeV and 1 TeV)",
     "check": lambda: abs(100 * _load("nu")["dis_channel"]["fraction_below_q2_floor"] - 16.9) < 0.1,
     "detail": lambda: f"{100*_load('nu')['dis_channel']['fraction_below_q2_floor']:.2f}%"},
    {"what": "the inputs are the tungsten combination at 100 GeV, in the "
             "benchmark region",
     "check": lambda: all(_load(c)["target"] == "W" and _load(c)["w_min"] == 3.0
                          and abs(_load(c)["energy_gev"] - ENERGY) < 1e-6
                          for c in ("nu", "mu")),
     "detail": lambda: ", ".join(f"{c} {_load(c)['region']}" for c in ("nu", "mu"))},
]


def main():
    a1.main(energy=ENERGY, output=OUTPUT)


if __name__ == "__main__":
    main()
