#!/usr/bin/env python3
"""Neutrino spectra through the electronic-detector fiducial cylinder.

FASER's electronic-detector measurements (arXiv:2412.03186, 2026 update) select
interactions inside a cylinder of 100 mm radius about the spectrometer axis.
The vendored spectra of data/faser_flux_2025/ are counts through the whole
25 x 30 cm tungsten face, and the flux is NOT flat over it: the high-energy
neutrinos are collimated about the line of sight, so the count inside the
cylinder is 2% (10-100 GeV) to 25% (> 1.5 TeV) above the face count scaled by
area.  Scaling by area under-predicts the electronic-detector rate by about 6%
overall, rising with energy -- silently, since the cross-section is unaffected.

This script counts the neutrinos themselves, from the per-neutrino text files
of $FASER_DATA/fluxes/2025/text/events/minus/ (written by
faser_format/flux_to_text.py; one row per neutrino, columns w_per_fb x_mm y_mm
px py pz E pdg parent, FASER global frame, the crossing angle already in the
positions), inside r < 100 mm of x = y = 0 (the spectrometer axis), and writes
spectra in exactly the format of the face files:

    data/faser_flux_2025/<name>_fid_r100.txt

Run once, after sourcing config.sh (it needs only numpy and pandas):

    tools/make_flux_fid_spectra.py
"""
import datetime
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
FASER_DATA = os.environ.get("FASER_DATA")
if not FASER_DATA:
    sys.exit("FASER_DATA not set: source config.sh first")
SRC = os.path.join(FASER_DATA, "fluxes", "2025", "text", "events", "minus")
OUT = os.path.join(REPO, "data", "faser_flux_2025")
NAMES = ("EPOSLHC_light", "Powheg_pythia_charm_central",
         "Powheg_pythia_charm_min", "Powheg_pythia_charm_max")
R_MM, X0, Y0 = 100.0, 0.0, 0.0
EDGES = np.logspace(1, 4, 31)          # the face files' bins: 10 per decade
PDGS = (12, -12, 14, -14, 16, -16)


def main():
    for name in NAMES:
        fn = os.path.join(SRC, f"{name}.txt.gz")
        t = pd.read_csv(fn, sep=r"\s+", comment="#", header=None,
                        names="w x y px py pz E pdg parent".split())
        cyl = (t.x - X0) ** 2 + (t.y - Y0) ** 2 < R_MM ** 2
        cols = [EDGES[:-1], EDGES[1:]]
        for pdg in PDGS:
            s = t[cyl & (t.pdg == pdg)]
            n = np.histogram(s.E, EDGES, weights=s.w)[0]
            e2 = np.histogram(s.E, EDGES, weights=s.w ** 2)[0]
            cols += [n, np.sqrt(e2)]
        out = os.path.join(OUT, f"{name}_fid_r100.txt")
        hdr = (f"Kling 2025 FASER neutrino flux, {name}.root, minus 160 urad "
               f"crossing angle, through r < {R_MM:g} mm of (x, y) = ({X0:g}, "
               f"{Y0:g}) mm at z = 0 (the electronic-detector fiducial "
               f"cylinder).\nNeutrinos per fb^-1 in each energy bin (NOT per "
               f"GeV); err = MC statistical error, sqrt(sum w^2).  Written "
               f"{datetime.date.today()} by tools/make_flux_fid_spectra.py.\n"
               "E_lo_GeV E_hi_GeV N_nue err_nue N_nuebar err_nuebar N_numu "
               "err_numu N_numubar err_numubar N_nutau err_nutau N_nutaubar "
               "err_nutaubar")
        np.savetxt(out, np.column_stack(cols), fmt="%.6e", header=hdr)
        print(f"wrote {os.path.relpath(out, REPO)}")


if __name__ == "__main__":
    main()
