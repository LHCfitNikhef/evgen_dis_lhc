#!/usr/bin/env python3
"""Comparison plots for the muon-DIS generator benchmark.

Reads results/histos_<gen>.json (from analyze.py) and produces, for each
view and each observable, an overlay of dsigma/dO with a ratio panel
(each generator over the view's reference).

Views (see VIEWS below):
  lo         LO generators + shower/hadronisation, vs analytic YADISM LO (ZM-VFNS)
  nlo        the two NLO-matched samples, vs analytic YADISM NLO (ZM-VFNS)
  orders     one LO MC (Pythia) vs one NLO MC (POWHEG) vs the ZM-VFNS YADISM LO/NLO/NNLO ladder
  me_lo      matrix-element level, LO samples vs analytic YADISM LO (ZM-VFNS)
  me_orders  matrix-element level, Pythia ME vs the ZM-VFNS YADISM LO/NLO/NNLO ladder

The two me_* views mirror the hadron-level lo / orders views on the ME-level
samples (muon-side observables only). There is no ME-level counterpart of the
hadron-level "nlo" view: no NLO-matched ME sample exists (see ME_NLO_MISSING).

Usage: plot_compare.py [view ...]   (no argument = all views)
       plot_compare.py me           (alias for both me_* views)
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")
# the house style, applied BEFORE pyplot is imported (analysis/plotstyle.py)
import plotstyle  # noqa: E402
plotstyle.apply()
from plotstyle import tex  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
RESULTS = f"{BASE}/results"

# --- energy scan -----------------------------------------------------------
# Results and FIGURES are both keyed by beam energy.  The anchor keeps its
# original untagged names so nothing published moves; other energies get a
# suffix, from the one helper that owns the rule.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import beams  # noqa: E402
import selection  # noqa: E402
import labels  # noqa: E402  -- the one display-name table

ENERGY = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)


def at_e(name):
    """Per-energy result/figure name -- see beams.at_energy()."""
    return beams.at_energy(name, ENERGY)


# --- FASER selection tiers -------------------------------------------------
# Results and FIGURES are keyed by (selection, energy), in that order, exactly
# as analyze.result_path() writes them: histos_<gen>[_<sel>][_<energy>].json.
# The inclusive selection adds no suffix, so every published figure keeps its
# name.  sel_suffix() is imported rather than re-derived: duplicating a naming
# rule is how job_400GeV_* came to be matched by job_*.
from analyze import sel_suffix, SELECTION  # noqa: E402


def at_sel_e(name):
    """Per-(selection, energy) result/figure name."""
    return at_e(name + sel_suffix())

COLORS = {"sherpa": "#2a78d6", "pythia": "#eb6834", "sherpa_lo": "#3fa66a",
          "herwig": "#8e5bd0", "genie": "#c2317b", "powheg": "#0f9b9b",
          "genie_nnpdf": "#e88bb8",
          # MG5_aMC, at parton level and showered by Pythia: one colour, since
          # no view carries both
          "mg5_me": "#b8860b", "mg5_ps": "#b8860b",
          "pythia_me": "#eb6834", "sherpa_lo_me": "#3fa66a",
          "herwig_me": "#8e5bd0", "sherpa_nlo_me": "#2a78d6",
          "powheg_lhe": "#0f9b9b",
          "yadism": "#33383f", "yadism_nlo": "#2a78d6",
          "yadism_nnlo": "#7fb3e8",
          # charm view: the two Pythia charm samples are two shades of the
          # Pythia orange, the analytic charm reference keeps the YADISM ink
          "pythia_me_charm": "#eb6834", "pythia_charm": "#a33d12",
          "sherpa_lo_charm": "#2a78d6", "herwig_charm": "#7a6ba8",
          "yadism_charm": "#33383f",
          "yadism_charm_nlo": "#6b7280", "yadism_charm_nnlo": "#9aa1ab",
          # neutrino-side NLO samples: keep Sherpa blue / POWHEG teal
          "sherpa_nlo": "#2a78d6", "powheg_nu": "#0f9b9b",
          # Herwig NLO keeps the Herwig purple in both currents
          "herwig_nlo_powheg": "#8e5bd0",
          "herwig_nlo_powheg_full": "#8e5bd0", "herwig_nlo_full": "#8e5bd0",
          # CC charm, final-state-charm tag (2026-08-21): each generator keeps
          # its inclusive colour so the charm view reads against the LO/NLO
          # views without a second colour key
          "sherpa_lo_charmfinal": "#3fa66a", "pythia_charmfinal": "#eb6834",
          "herwig_charmfinal": "#8e5bd0", "genie_lo_charmfinal": "#c2317b",
          "sherpa_nlo_charmfinal": "#2a78d6",
          "powheg_nu_charmfinal": "#0f9b9b",
          "herwig_nlo_charmfinal": "#6a3fb0",
          "genie_charmfinal": "#e88bb8",
          "genie_lo": "#c2317b",
          "powheg_charm": "#0f9b9b",
          # muon NC charm, final-state tag: each keeps its inclusive colour
          "sherpa_lo_charmfinal": "#3fa66a", "herwig_charmfinal": "#8e5bd0",
          "pythia_charmfinal": "#a33d12", "pythia_me_charmfinal": "#eb6834",
          "sherpa_charmfinal": "#2a78d6", "powheg_charmfinal": "#0f9b9b",
          "herwig_nlo_powheg_charmfinal": "#6a3fb0",
          "herwig_nlo_powheg_full_charmfinal": "#6a3fb0",
          "herwig_nlo_full_charmfinal": "#6a3fb0",
          "yadism_charm_nlo_fonll_damp": "#c2317b",
          "yadism_charm_nnlo_fonll_damp": "#e88bb8"}

# short names for the ratio-panel axis label (the json labels carry versions)
SHORT = {"sherpa_lo": "Sherpa LO", "sherpa": "Sherpa MC@NLO",
         "pythia": "Pythia LO", "herwig": "Herwig LO", "genie": "GENIE",
         "genie_nnpdf": "GENIE NNPDF4.0",
         "powheg": "POWHEG-RES NLO", "powheg_lhe": "POWHEG-RES NLO (LHE)",
         "yadism": "YADISM LO (ZM-VFNS)",
         "sherpa_lo_me": "Sherpa LO (ME)", "pythia_me": "Pythia LO (ME)",
         "herwig_me": "Herwig LO (ME)", "sherpa_nlo_me": "Sherpa NLO (ME)",
         "yadism_nlo": "YADISM NLO (ZM-VFNS)",
         "yadism_nnlo": "YADISM NNLO (ZM-VFNS)",
         "sherpa_lo_charm": "Sherpa LO (charm)",
         "herwig_charm": "Herwig LO (charm)",
         "pythia_me_charm": "Pythia LO (ME, charm)",
         "pythia_charm": "Pythia LO (charm)",
         "yadism_charm": "YADISM LO (ZM-VFNS, charm)",
         "yadism_charm_nlo": "YADISM NLO (ZM-VFNS, charm)",
         "yadism_charm_nnlo": "YADISM NNLO (ZM-VFNS, charm)",
         # neutrino-side NLO samples (results_nu/)
         "sherpa_nlo": "Sherpa MC@NLO", "powheg_nu": "POWHEG-V2 NLO",
         # the QED study's reference arms -- the ratio ylabel is set from
         # this table, and the CC baseline's stored label is long enough to
         # run the whole height of the panel
         "powheg_res_nu": "POWHEG-RES, QED off",
         "pythia_shw_simple": "Pythia, simple shower",
         "sherpa_lo_dire": "Sherpa LO, Dire",
         "sherpa_dire": "Sherpa MC@NLO, Dire",
         "sherpa_nlo_dire": "Sherpa MC@NLO, Dire",
         "herwig_nlo_powheg": "Herwig NLO (POWHEG)",
         "herwig_nlo_powheg_full": "Herwig NLO (POWHEG)",
         "herwig_nlo_full": "Herwig NLO (POWHEG)",
         "herwig_nlo_powheg_full_charmfinal": "Herwig NLO (charm)",
         "herwig_nlo_full_charmfinal": "Herwig NLO (charm)",
         "genie_lo": "GENIE (GRV98LO)",
         "sherpa_lo_charmfinal": "Sherpa LO (charm)",
         "pythia_charmfinal": "Pythia LO (charm)",
         "herwig_charmfinal": "Herwig LO (charm)",
         "genie_lo_charmfinal": "GENIE (GRV98LO, charm)",
         "sherpa_nlo_charmfinal": "Sherpa MC@NLO (charm)",
         "powheg_nu_charmfinal": "POWHEG-V2 NLO (charm)",
         "herwig_nlo_charmfinal": "Herwig NLO (charm)",
         "genie_charmfinal": "GENIE (HEDIS, charm)",
         "powheg_charm": "POWHEG-RES NLO (charm)",
         "sherpa_lo_charmfinal": "Sherpa LO (charm)",
         "herwig_charmfinal": "Herwig LO (charm)",
         "pythia_charmfinal": "Pythia LO (charm)",
         "pythia_me_charmfinal": "Pythia LO (ME, charm)",
         "sherpa_charmfinal": "Sherpa MC@NLO (charm)",
         "powheg_charmfinal": "POWHEG-RES NLO (charm)",
         "herwig_nlo_powheg_charmfinal": "Herwig NLO (charm)",
         "yadism_charm_nlo_fonll_damp": "YADISM NLO (charm, FONLL)",
         "yadism_charm_nnlo_fonll_damp": "YADISM NNLO (charm, FONLL)"}

# muon-side observables: the only ones the analytic YADISM curves provide
DIS_OBS = ("Q2", "xbj", "y", "Emu", "theta", "nu")

# samples that have NOT been validated against the analytic reference and must
# not be read as part of the closure. Their legend entry is marked, and
# make_report.py flags them in the table and blurbs. GENIE sits 17% below
# YADISM LO. The genie_nnpdf variant (identical config, GRV98 LO ->
# NNPDF4.0 NNLO inside BYPDF) recovers +7% of rate (31.67 -> 34.03 nb), so
# ~60% of the deficit is the PDF; the remaining ~11% is Bodek-Yang model
# content (x-shift + higher-twist corrections tuned at few GeV, photon-only
# exchange), not a configuration problem.
# GENIE is NOT flagged "preliminary".  It uses a different physical model for
# the hard scattering (Bodek-Yang with GRV98 LO PDFs, photon exchange only,
# its own EW parameters), so it is not expected to reproduce the analytic
# reference and its offset is a RESULT, not a caveat about data quality.
# The model differences are stated in the report blurbs instead.
PRELIMINARY = set()

# Each view: which generators to overlay, which one the ratio panel divides
# by (must be an MC sample -- YADISM has no N_ch / E_lead histograms), the
# output filename prefix, and any per-view colour overrides.
VIEWS = {
    # -------------------------------------------------- SHOWER, Sherpa
    # THE NLO HALF of the shower study, and its LO control.  Same generator,
    # same matrix element, same matching, same PDF, same scales AND THE SAME
    # CLUSTER HADRONISATION -- only SHOWER_GENERATOR differs.  That last point
    # is what the Pythia LO study cannot claim: Vincia and Dire ship their own
    # tunes, so switching model there switches tune too.
    #
    # THE PAIR IS THE POINT.  The question is whether an NLO-matched sample is
    # stiffer against the choice of shower than a LO one, so the LO arms are
    # run with the SAME two showers in the SAME generator.  Comparing the
    # Pythia LO spread against a Sherpa NLO spread would answer nothing.
    "shower_sherpa_lo": dict(
        gens=["sherpa_lo", "sherpa_lo_dire"],
        ref="sherpa_lo", prefix="cmp_shwsherpa_lo_",
        subtitle="shower model at LO: Sherpa CSS (published) vs Dire",
        colors={"sherpa_lo": "#33383f", "sherpa_lo_dire": "#eb6834"},
    ),
    "shower_sherpa_nlo": dict(
        gens=["sherpa", "sherpa_dire"],
        ref="sherpa", prefix="cmp_shwsherpa_nlo_",
        subtitle="shower model at NLO: Sherpa MC@NLO with CSS (published) vs Dire",
        colors={"sherpa": "#33383f", "sherpa_dire": "#2b6cb0"},
    ),
    # ------------------------------------------------------------- SHOWER
    # THE SAME 250000 LO EVENTS IN EVERY ARM, same card and same seed, with
    # only PartonShowers:model changed.  The reference is the simple shower,
    # which is BYTE-IDENTICAL to pythia8/job_1 (logs/shower_closure_simple
    # .log), so the ratio panels are maximally correlated.
    #
    # DONE AT LO, NOT ON THE POWHEG SAMPLES.  The POWHEG veto is a
    # SimpleShower feature that Vincia and Dire do not expose, and both arms
    # deliver essentially nothing there (logs/shower_probe_powheg.log).  At LO
    # there is no matching veto, so the shower is the only thing that differs.
    #
    # WHAT TO EXPECT, and what it means if the expectation fails: with QED off
    # and dipoleRecoil on, the shower does not touch the SCATTERED LEPTON, so
    # Q2, x and y should be essentially identical across the arms.  The
    # differences belong in the HADRONIC observables -- multiplicity, leading
    # hadron, the D-meson panels.  A shift in Q2 or y would mean the recoil
    # scheme is moving the lepton, which would be a finding in itself.
    "shower": dict(
        gens=["pythia_shw_simple", "pythia_shw_vincia", "pythia_shw_dire"],
        ref="pythia_shw_simple", prefix="cmp_shower_",
        subtitle=("parton-shower model on identical LO events: simple "
                  "(published), Vincia antenna, Dire dipole"),
        colors={"pythia_shw_simple": "#33383f",
                "pythia_shw_vincia": "#eb6834",
                "pythia_shw_dire": "#3fa66a"},
    ),
    # ---------------------------------------------------------------- QED
    # THE SAME LHE EVENTS IN EVERY ARM, re-showered with photon radiation
    # switched on in stages, so the matrix element, the matching, the PDF and
    # the scale are held fixed and the ratio panels show the radiation alone.
    # The reference is the QED-OFF arm, which is the published `powheg` sample
    # itself: re-showering with QED off reproduces it byte-identically
    # (logs/qed_closure_off.log), which is what makes this study
    # self-validating.
    #
    # Q2, x and y are reconstructed from the SCATTERED LEPTON, so a photon
    # radiated before reconstruction moves the event in the fiducial variables
    # themselves.  The quark-line arm does not touch those and shows up in the
    # hadronic panels instead.
    "qed": dict(
        gens=["powheg", "powheg_qed_fsr", "powheg_qed_fsrisr",
              "powheg_qed_full"],
        ref="powheg", prefix="cmp_qed_",
        subtitle=("QED radiation on the same POWHEG-RES NLO events: off "
                  "(published), lepton FSR, lepton FSR+ISR, and quark line"),
        colors={"powheg": "#33383f", "powheg_qed_fsr": "#eb6834",
                "powheg_qed_fsrisr": "#2b6cb0", "powheg_qed_full": "#3fa66a"},
    ),
    "lo": dict(
        # mg5_ps is the MG5_aMC matrix element showered and hadronised by
        # Pythia 8.  Its partner in this view is the Pythia row: same shower,
        # same hadronisation, same parton distribution and same scale, so the
        # pair isolates the matrix element and the phase-space generation.
        gens=["sherpa_lo", "pythia", "herwig", "genie", "genie_nnpdf",
              "mg5_ps", "yadism"],
        ref="sherpa_lo", prefix="cmp_lo_",
        subtitle="LO matrix elements + parton shower + hadronisation",
    ),
    "nlo": dict(
        # herwig_nlo_powheg, NOT herwig_nlo: the Matchbox MC@NLO sample has
        # N_eff = 6-10% and max|w| ~ 1131, so its rate is quotable but its
        # shapes are not.  The PowhegMEDISNC sample carries UNIT WEIGHTS
        # (N_eff = 100%, zero negative), so it belongs in the plots.
        gens=["sherpa", "powheg", "herwig_nlo_powheg_full", "yadism_nlo"],
        ref="sherpa", prefix="cmp_nlo_",
        subtitle="NLO-matched samples",
        # yadism_nlo shares Sherpa's blue in the combined palette
        colors={"yadism_nlo": "#33383f"},
    ),
    "orders": dict(
        gens=["pythia", "powheg", "yadism", "yadism_nlo", "yadism_nnlo"],
        ref="powheg", prefix="cmp_orders_",
        subtitle="perturbative order",
        # theory ladder as a grey ramp darkening with order; MCs keep their colours
        colors={"yadism": "#9aa0a6", "yadism_nlo": "#5f5e5a",
                "yadism_nnlo": "#1f1f1e"},
    ),
    # ME-level views mirror the hadron-level ones as far as the samples allow;
    # the analytic curves are the reference here, so they stay solid
    "me_lo": dict(
        gens=["sherpa_lo_me", "pythia_me", "herwig_me", "yadism"],
        ref="sherpa_lo_me", prefix="cmp_me_lo_",
        subtitle="matrix-element level, LO samples (no shower / hadronisation)",
        obs=DIS_OBS, yadism_dashed=False,
    ),
    "me_nlo": dict(
        # the fixed-order sample is the ratio reference: it is the validated,
        # shower-independent NLO calculation, so POWHEG LHE (NLO + first
        # emission, Sudakov-resummed) is read against it and not vice versa
        gens=["sherpa_nlo_me", "powheg_lhe", "yadism_nlo"],
        ref="sherpa_nlo_me", prefix="cmp_me_nlo_",
        subtitle=("matrix-element level, NLO: fixed-order and POWHEG-RES LHE "
                  "(NLO + first emission) vs analytic YADISM NLO in ZM-VFNS"),
        obs=DIS_OBS, yadism_dashed=False,
        colors={"yadism_nlo": "#33383f"},
    ),
    # CHARM PRODUCTION: charm-INITIATED scattering gamma*/Z + c -> c off the
    # charm PDF, i.e. a charm quark in BOTH the initial and the final state of
    # the hard process. NOT photon-gluon fusion gamma* g -> c cbar, which is
    # O(alpha_s) and belongs to charm-PAIR production.
    # The ratio reference is the ME-level Pythia sample: it is the
    # apples-to-apples match for the analytic calculation (particle level
    # before shower and hadronisation) and, being an MC sample, it carries the
    # per-event statistics the ratio band needs. The hadron-level charm-tagged
    # sample rides along to show what the shower does to the tag.
    # Split into three sub-views (2026-08-21) so the charm tab reads like the
    # inclusive ones: LO generators / NLO generators / perturbative order.
    # It used to be a single plot carrying the LO generators AND the whole
    # LO-NLO-NNLO analytic ladder, which made the generator comparison and
    # the order comparison fight for the same axes.
    "charm": dict(
        gens=["sherpa_lo_charmfinal", "herwig_charmfinal",
              "pythia_me_charmfinal", "pythia_charmfinal", "yadism_charm"],
        ref="pythia_me_charmfinal", prefix="cmp_charm_",
        subtitle=("charm production at LO: $\\geq 1$ prompt charm quark in the "
                  "final state"),
        obs=DIS_OBS, yadism_dashed=False,
    ),
    "charm_nlo": dict(
        gens=["sherpa_charmfinal", "powheg_charmfinal",
              "herwig_nlo_powheg_full_charmfinal", "yadism_charm_nlo",
              "yadism_charm_nlo_fonll_damp"],
        ref="yadism_charm_nlo_fonll_damp", prefix="cmp_charm_nlo_",
        subtitle=("charm production at NLO: $\\geq 1$ prompt charm quark in the "
                  "final state"),
        obs=DIS_OBS, yadism_dashed=False,
        colors={"yadism_charm_nlo": "#33383f",
                "yadism_charm_nlo_fonll_damp": "#c2317b"},
    ),
    "charm_orders": dict(
        gens=["pythia_charmfinal", "powheg_charmfinal", "yadism_charm",
              "yadism_charm_nlo", "yadism_charm_nnlo",
              "yadism_charm_nlo_fonll_damp"],
        # THE LADDER IS SHOWN AGAINST THE MASS-CORRECT SCHEME.  The point of
        # this view is the ZM non-convergence, and a ratio to one of the ZM
        # rungs shows it only relative to itself; against FONLL each
        # rung's distance from the physical answer is what is plotted.
        ref="yadism_charm_nlo_fonll_damp", prefix="cmp_charm_orders_",
        subtitle=("perturbative order, charm \u2014 the ZM-VFNS charm series "
                  "does not converge at these kinematics"),
        obs=DIS_OBS, yadism_dashed=False,
        colors={"yadism_charm": "#9aa0a6", "yadism_charm_nlo": "#5f5e5a",
                "yadism_charm_nnlo": "#1f1f1e"},
    ),
    "me_orders": dict(
        gens=["pythia_me", "yadism", "yadism_nlo", "yadism_nnlo"],
        ref="pythia_me", prefix="cmp_me_orders_",
        subtitle="matrix-element level, perturbative order",
        obs=DIS_OBS, yadism_dashed=False,
        colors={"yadism": "#9aa0a6", "yadism_nlo": "#5f5e5a",
                "yadism_nnlo": "#1f1f1e"},
    ),
}

# one name per hadron-level view group, for convenience on the command line
ALIASES = {"me": ["me_lo", "me_orders"],
           "hadron": ["lo", "nlo", "orders"]}

# why there is no ME counterpart of the hadron-level "nlo" view -- quoted by
# make_report.py so the page says it rather than silently dropping the tab
ME_NLO_MISSING = (
    "No NLO-matched matrix-element sample exists: neither Sherpa MC@NLO nor "
    "POWHEG-RES was run with shower and hadronisation switched off, so there "
    "is nothing to compare at this order.")

INK, MUTED, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e9e8e4", "#fcfcfb"

AXES = {
    "Q2":    (r"$Q^2$  [GeV$^2$]", r"$d\sigma/dQ^2$  [pb/GeV$^2$]", "log", "log"),
    "xbj":   (r"$x_{\rm Bj}$", r"$d\sigma/dx$  [pb]", "log", "log"),
    "y":     (r"$y$", r"$d\sigma/dy$  [pb]", "linear", "linear"),
    "Emu":   (r"$E_{\mu}'$ (lab)  [GeV]", r"$d\sigma/dE_{\mu}'$  [pb/GeV]",
              "linear", "linear"),
    "theta": (r"$\theta_{\mu}$ (lab)  [rad]", r"$d\sigma/d\theta_{\mu}$  [pb/rad]",
              "log", "log"),
    "nch":   (r"$N_{\rm ch}$ (charged hadrons)", r"$d\sigma/dN_{\rm ch}$  [pb]",
              "linear", "log"),
    "nch1":  (r"$N_{\rm ch}$ ($E_{\rm lab}>1$ GeV)", r"$d\sigma/dN_{\rm ch}$  [pb]",
              "linear", "log"),
    "nu":    (r"$\nu = E_h$  [GeV]", r"$d\sigma/d\nu$  [pb/GeV]",
              "linear", "linear"),
    "Elead": (r"$E$ leading charged hadron (lab)  [GeV]", r"$d\sigma/dE$  [pb/GeV]",
              "log", "log"),
    "dphi":  (r"$\Delta\phi_{\min}(\mu, h^{\pm})$  [rad]",
              r"$d\sigma/d\Delta\phi$  [pb/rad]", "linear", "log"),
    "nd_ch": (r"$N(D^{\pm})$ per event", r"$d\sigma/dN(D^{\pm})$  [pb]",
              "linear", "log"),
    "nd_0":  (r"$N(D^0/\bar{D}^0)$ per event",
              r"$d\sigma/dN(D^0)$  [pb]", "linear", "log"),
    "Ed_ch": (r"$E$ leading $D^{\pm}$ (lab)  [GeV]",
              r"$d\sigma/dE_{D^{\pm}}$  [pb/GeV]", "log", "log"),
    "Ed_0":  (r"$E$ leading $D^0/\bar{D}^0$ (lab)  [GeV]",
              r"$d\sigma/dE_{D^0}$  [pb/GeV]", "log", "log"),
}

# observables that only exist where there is hadronisation: the ME-level and
# analytic entries have no such histogram. plot_obs() already skips a key that
# is missing from a sample's "hists", and the me_* views restrict themselves to
# DIS_OBS, so these drop out of those views without further plumbing.
HADRONISATION_OBS = ("nch", "nch1", "Elead", "dphi",
                     "nd_ch", "nd_0", "Ed_ch", "Ed_0")

# ratio-panel SCALE per observable. The limits are no longer fixed: they are
# computed per figure from the actual spread (see ratio_limits), because the
# DIS-kinematics observables agree at the few-% level while the hadronisation
# observables differ by orders of magnitude -- one hard-coded range for either
# leaves most of the panel empty.
# The D-meson observables join the hadronisation set: the multiplicities span
# the 0 bin (~95% of sigma) against tails that differ by large factors between
# generators, and the charm channel itself differs by ~10% in rate, so a linear
# ratio axis cannot hold both ends.
RATIO_LOG = HADRONISATION_OBS

# a ratio bin only constrains the axis if it is statistically meaningful;
# beyond this relative error it is noise and would blow the range up
RATIO_MAX_RELERR = 0.5
# never shrink below this much total span, so a flat ratio keeps a sane axis
RATIO_MIN_SPAN_LIN = 0.10      # +-5% around 1
RATIO_MIN_SPAN_LOG = 1.15      # multiplicative


def ratio_limits(series, scale, band=None):
    """Auto-range a ratio panel to the natural variability of the ratios.

    `series` is a list of (ratio, ratio_err) arrays, one per non-reference
    generator; `band` is the reference's own relative error per bin. A bin
    only constrains the axis if BOTH sides of the ratio are statistically
    meaningful -- in the sparse tails a single-event reference bin produces a
    huge but meaningless ratio that would flatten everything else. (The large
    ratios that survive this cut are real: Herwig's N_ch tail genuinely runs
    an order of magnitude above Sherpa's.) 1.0 is always inside the window.
    """
    lo, hi = [], []
    ref_ok = None
    if band is not None:
        band = np.asarray(band, dtype=float)
        ref_ok = np.isfinite(band) & (band <= RATIO_MAX_RELERR) & (band > 0)
    for r, re in series:
        r = np.asarray(r, dtype=float)
        re = np.asarray(re, dtype=float)
        ok = np.isfinite(r) & (r > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            good = ok & np.isfinite(re) & (re <= RATIO_MAX_RELERR*np.abs(r))
        if ref_ok is not None and ref_ok.shape == good.shape:
            good &= ref_ok
        use = r[good] if good.any() else r[ok]
        if use.size:
            lo.append(use.min())
            hi.append(use.max())
    if not lo:
        return (0.5, 1.5) if scale == "linear" else (0.1, 10.0)
    lo, hi = min(lo + [1.0]), max(hi + [1.0])

    if scale == "log":
        span = max(hi/lo, RATIO_MIN_SPAN_LOG)
        pad = span**0.06
        mid = (lo*hi)**0.5
        lo, hi = min(lo, mid/span**0.5)/pad, max(hi, mid*span**0.5)*pad
    else:
        span = max(hi - lo, RATIO_MIN_SPAN_LIN)
        pad = 0.06*span
        mid = 0.5*(lo + hi)
        lo, hi = min(lo, mid - span/2) - pad, max(hi, mid + span/2) + pad
    return lo, hi

def beam_label(e_lab=None):
    """"1 TeV" / "400 GeV", from the energy actually being plotted."""
    e = ENERGY if e_lab is None else e_lab
    return f"{e/1000:g} TeV" if e >= 1000 else f"{e:g} GeV"


def region_label(sel=None):
    """The fiducial region ACTUALLY in force, as a title fragment.

    Built from selection.py rather than written out, for the same reason the
    result names are: it is the definition that cut the events.
    """
    sel = sel or SELECTION
    parts = [f"$Q^2>{sel.q2_min:g}$ GeV$^2$",
             f"${sel.y_min:g}<y<{sel.y_max:g}$"]
    if sel.e_lep_min:
        parts.append(f"$E'>{sel.e_lep_min:g}$ GeV")
    if sel.theta_max:
        parts.append(f"$\\theta'<{1e3 * sel.theta_max:g}$ mrad")
    if sel.theta_min:
        parts.append(f"$\\theta'>{1e3 * sel.theta_min:g}$ mrad")
    if sel.has_hadronic_cuts:
        parts.append("hadronic vertex cuts")
    if sel.dimuon:
        parts.append("opposite-sign 2nd $\\mu$")
    return ", ".join(parts)


def make_title(beam="NC muon DIS", lepton=r"$\mu^-$"):
    """>>> THE FIGURE TITLE MUST TRACK THE FIGURE. <<<

    This was a hardcoded string reading "1 TeV ... $Q^2>4$ GeV$^2$,
    $0.2<y<0.9$" on EVERY figure the module produced -- so a 400 GeV Tier E
    plot was titled 1 TeV and announced the inclusive region, while the data
    in it was neither.  Over a thousand figures carried it.  The axes were
    right and only the title lied, which is why nothing looked wrong.

    Same failure as the report blurbs quoting the anchor on all twelve tabs,
    and the fix is the same: derive it from the energy and selection actually
    in force instead of writing it down once.
    """
    head = f"{beam}, {beam_label()} {lepton} on p at rest"
    if SELECTION.name == "inclusive":
        # one line, and byte-identical to what the inclusive figures have
        # always carried, so nothing already published moves
        return f"{head},  {region_label()},  NNPDF4.0 NNLO"
    # A tier region does not fit on one line, and letting it run wide does not
    # merely look bad: the figures are saved with bbox_inches="tight", so an
    # over-long title STRETCHES THE WHOLE CANVAS and the plot comes out squat
    # and twice the width of its inclusive counterpart. Break it deliberately
    # rather than at a character count -- the region string is full of $...$
    # and a naive wrap would split a LaTeX group.
    return (f"{head},  NNPDF4.0 NNLO\n"
            f"{SELECTION.label}:  {region_label()}")


TITLE = make_title()


# ---------------------------------------------------------------- rebinning
# ADAPTIVE TAIL MERGING (user, 2026-08-27): every bin drawn should carry
# roughly the same Monte Carlo error, so bins are merged until they do and the
# drawn bins end up unequal in width.
#
# WHY AT PLOT TIME AND NOT IN analyze.py's BINS.  Merging here is arithmetic on
# the stored histogram -- no event file is touched, so it costs nothing and is
# reversible.  Changing BINS instead would mean re-parsing 120 GB, and every
# sample would have to be re-parsed TOGETHER: a generator left on the old edges
# would silently stop being comparable.
#
# WHY THE BULK IS LEFT ALONE.  The muon NLO sample is already at 2-3% relative
# error in the median bin and only its tails are bad (4 of 24 Q2 bins above
# 20%, worst 58%).  A uniform rebin would throw away resolution across the
# well-measured bulk to fix four bins.
#
# THE GROUPING IS DECIDED ONCE, FROM THE RATIO REFERENCE, and applied to every
# curve in the view.  It has to be: the ratio panel divides curves bin by bin,
# so two curves on different edges cannot be divided at all.  The reference is
# the right one to drive it because it is the denominator of every ratio, so
# its noise propagates into all of them.
# 2% RATHER THAN 10% (user, 2026-08-28).  10% fixed the 100%-error tails but
# left the bulk untouched, and the bulk is where the ratio panels look ragged:
# the visible bin-to-bin scatter in a ratio is the REFERENCE's statistics,
# shared by every curve, which is why the analytic YADISM line wiggles in step
# with the Monte Carlo ones.  Measured on the muon NLO E'mu panel, where the
# reference carries 2.98% per bin:
#
#     target   bins   ref err/bin   ratio scatter
#      10%      35       2.98%          4.14%
#       3%      26       2.42%          3.64%
#       2%      13       1.75%          2.72%
#     1.5%       8       1.43%          1.33%
#
# 2% keeps thirteen bins across 100-800 GeV -- enough to read the shape -- and
# roughly halves the scatter.  Tighter than that buys smoothness by throwing
# away resolution.
MERGE_TARGET = float(os.environ.get("BENCH_REBIN_TARGET", "0.02"))


def merge_groups(edges, y, e, target):
    """Index groups whose merged relative error is <= target, greedily.

    Walks left to right accumulating bins until the running group is precise
    enough, then starts a new one.  A trailing group that never gets there is
    folded into its predecessor rather than left under-determined -- otherwise
    the very last bin, which is exactly the one this exists for, would survive
    unmerged.
    """
    w = np.diff(edges)
    groups, start = [], 0
    acc_i, acc_e2 = 0.0, 0.0
    for k in range(len(y)):
        acc_i += y[k] * w[k]
        acc_e2 += (e[k] * w[k]) ** 2
        good = acc_i > 0 and math.sqrt(acc_e2) <= target * acc_i
        if good:
            groups.append((start, k + 1))
            start, acc_i, acc_e2 = k + 1, 0.0, 0.0
    if start < len(y):                      # trailing remainder
        if groups:
            groups[-1] = (groups[-1][0], len(y))
        else:
            groups = [(0, len(y))]
    return groups


def apply_merge(edges, y, e, groups):
    """Collapse a density histogram onto `groups`, conserving the integral."""
    w = np.diff(edges)
    ne = [edges[groups[0][0]]]
    ny, nev = [], []
    for a, b in groups:
        width = w[a:b].sum()
        integral = float((y[a:b] * w[a:b]).sum())
        abserr = float(np.sqrt(((e[a:b] * w[a:b]) ** 2).sum()))
        ne.append(edges[b])
        ny.append(integral / width if width else 0.0)
        nev.append(abserr / width if width else 0.0)
    return np.array(ne), np.array(ny), np.array(nev)


# MULTIPLICITY HISTOGRAMS GET THREE BINS: 0, 1, 2+ (user, 2026-08-28).
# The analysis books these 0..5, but bins 2 and above hold a per-mille tail
# that is mostly empty and, drawn on a log axis, spends most of the width on
# nothing.  Collapsing is done HERE rather than in analyze.py because the
# binning lives in the stored histograms: rebooking it there would mean
# re-parsing every HepMC sample to change a picture.
COUNT_2PLUS = {"nd_ch", "nd_0", "nd_s"}


def collapse_2plus(edges, y, e):
    """Collapse a unit-bin multiplicity onto 0, 1, 2+ on UNIT-WIDTH bins.

    The merged bin carries the INTEGRAL -- the cross-section for two or more
    -- on a bin of width one, not a density spread over the four unit bins it
    came from.  Going through the generic density merge instead would divide
    that integral by four and draw the 2+ bar a factor of four low while every
    number behind it stayed right, which on a log axis reads as a real
    suppression rather than as a plotting convention.
    """
    y = np.asarray(y, dtype=float)
    e = np.asarray(e, dtype=float)
    w = np.diff(np.asarray(edges, dtype=float))
    iy, ie = y * w, e * w                       # per-bin integrals
    ny = np.array([iy[0], iy[1], float(iy[2:].sum())])
    ne = np.array([ie[0], ie[1],
                   float(np.sqrt((ie[2:] ** 2).sum()))])
    return np.array([-0.5, 0.5, 1.5, 2.5]), ny, ne


def _count_ticks(active, *axs):
    """Label the three multiplicity bins 0, 1, 2+ and nothing between them.

    The last bin is drawn at x = 2 but MEANS two-or-more, so it has to say so:
    left to matplotlib's own ticker it reads "2", and a reader would take the
    distribution to stop there rather than to have been summed into it.
    """
    if not active:
        return
    for a in axs:
        if a is None:
            continue
        a.set_xticks([0, 1, 2])
        a.set_xticklabels(["0", "1", "2+"])
        a.set_xticks([], minor=True)


def style_axes(*axs):
    for a in axs:
        a.set_facecolor(SURFACE)
        a.grid(axis="y", color=GRID, lw=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_color(MUTED)
        a.tick_params(colors=MUTED, labelcolor=INK, which="both")


def plot_obs(key, data, view, fname):
    xlab, ylab, xs, ys = AXES[key]
    colors = {**COLORS, **view.get("colors", {})}
    dashed = view.get("yadism_dashed", True)
    # >>> A RATIO PANEL NEEDS SOMETHING TO DIVIDE. <<<
    # Under the FASER tiers a view can come down to a single surviving sample
    # -- the analytic references have no tier counterpart, and not every
    # generator was run at every tier and energy.  The panel then plotted that
    # sample against ITSELF: a flat line at one inside its own error band,
    # which reads as a comparison and is not one.  Count first, and drop the
    # panel when there is nothing to compare.
    ref_g = view["ref"]
    others = [g for g in view["gens"]
              if g != ref_g and g in data and key in data[g].get("hists", {})]
    with_ratio = bool(others)
    if with_ratio:
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(7.2, 6.2), facecolor=SURFACE, sharex=True,
            gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.0})
        style_axes(ax, axr)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(7.2, 5.0), facecolor=SURFACE)
        axr = None
        style_axes(ax)
    ref = data[ref_g]
    if key not in ref["hists"]:
        # the ratio reference has no such histogram (e.g. a sample whose
        # HepMC3 records cannot express the D mesons) -- nothing to divide by
        plt.close(fig)
        print(f"  (skipping {key}: reference {ref_g} has no such histogram)")
        return
    edges0 = np.asarray(ref["hists"][key]["edges"])
    ref_y = np.asarray(ref["hists"][key]["dsig"])
    ref_e = np.asarray(ref["hists"][key]["err"])
    # decide the grouping ONCE, from the reference, and keep the original
    # edges so every other curve can be collapsed onto the same one
    # a fixed 0/1/2+ collapse takes precedence over the adaptive merge: the
    # grouping is prescribed, so there is nothing for the error target to
    # decide, and letting both run would re-merge the three bins into fewer
    count_mode = key in COUNT_2PLUS and len(ref_y) > 3
    if count_mode:
        groups = None
        edges, ref_y, ref_e = collapse_2plus(edges0, ref_y, ref_e)
    elif MERGE_TARGET > 0 and len(ref_y) > 2:
        groups = merge_groups(edges0, ref_y, ref_e, MERGE_TARGET)
        edges, ref_y, ref_e = apply_merge(edges0, ref_y, ref_e, groups)
    else:
        groups, edges = None, edges0
    mids = 0.5*(edges[:-1] + edges[1:])

    ratios = []
    for g in view["gens"]:
        if g not in data:
            continue
        d = data[g]
        if key not in d["hists"]:
            continue
        yv = np.asarray(d["hists"][key]["dsig"])
        ev = np.asarray(d["hists"][key]["err"])
        if count_mode:
            _, yv, ev = collapse_2plus(edges0, yv, ev)
        elif groups is not None:
            _, yv, ev = apply_merge(edges0, yv, ev, groups)
        col = colors[g]
        ls = "--" if g.startswith("yadism") and dashed else "-"
        # labels.display(), not the raw JSON label: the legend and the
        # table row beside it must name the sample identically, and GENIE
        # appears four times under two modules and two tunes.
        lab = labels.display(d["label"]) + \
            ("  [PRELIMINARY]" if g in PRELIMINARY else "")
        ax.stairs(yv, edges, color=col, lw=2, ls=ls, label=lab)
        ax.errorbar(mids, yv, yerr=ev, fmt="none", ecolor=col,
                    elinewidth=1, capsize=0, alpha=0.7)
        if g != ref_g:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(ref_y > 0, yv/ref_y, np.nan)
                re = np.where(ref_y > 0, ev/ref_y, np.nan)
            if axr is not None:
                axr.stairs(r, edges, color=col, lw=2, ls=ls)
                axr.errorbar(mids, r, yerr=re, fmt="none", ecolor=col,
                             elinewidth=1, alpha=0.7)
                ratios.append((r, re))
    # reference stat band around 1
    with np.errstate(divide="ignore", invalid="ignore"):
        band = np.where(ref_y > 0, ref_e/ref_y, 0.0)
    if axr is not None:
        axr.fill_between(mids, 1.0-band, 1.0+band, step="mid",
                         color=colors[ref_g], alpha=0.25, lw=0)
        axr.axhline(1.0, color=MUTED, lw=0.8, ls="--")

    ax.set_xscale(xs)
    ax.set_yscale(ys)
    ax.set_ylabel(tex(ylab), color=INK)
    _h, _l = ax.get_legend_handles_labels()
    ax.legend(_h, [tex(x) for x in _l], frameon=True,
              labelcolor=INK, fontsize=plotstyle.FS_LEGEND)
    ax.set_title(tex(TITLE) + "\n" + tex(view['subtitle']),
                 color=INK, fontsize=plotstyle.FS_PANEL_TITLE)
    if axr is not None:
        rscale = "log" if key in RATIO_LOG else "linear"
        axr.set_yscale(rscale)
        axr.set_ylim(*ratio_limits(ratios, rscale, band))
        axr.set_ylabel(tex("ratio to") + "\n"
                       + tex(SHORT.get(ref_g, ref['label'])),
                       color=INK, fontsize=9)
        axr.set_xlabel(tex(xlab), color=INK)
        axr.set_xlim(edges[0], edges[-1])
        _count_ticks(count_mode, ax, axr)
    else:
        ax.set_xlabel(tex(xlab), color=INK)
        ax.set_xlim(edges[0], edges[-1])
        _count_ticks(count_mode, ax)
        ax.text(0.5, -0.16,
                tex("only one sample survives this selection "
                    "-- no ratio panel"),
                transform=ax.transAxes, ha="center",
                va="top", fontsize=9, style="italic", color=MUTED)
    fig.savefig(fname, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {os.path.basename(fname)}")


# how each tier cut reads in a summary header, in the order selection.py
# applies them.  Built from as_dict() so a cut added there cannot go unprinted.
CUT_LABELS = [("q2_min", "Q2 > {:g} GeV2"), ("y_min", "y > {:g}"),
              ("y_max", "y < {:g}"), ("e_lep_min", "E' > {:g} GeV"),
              ("theta_max", "theta' < {:g} rad"),
              ("theta_min", "theta' > {:g} rad"),
              ("n05_min", "N_ch(tan theta < 0.5) >= {:g}"),
              ("n01_min", "N_ch(tan theta < 0.1) >= {:g}"),
              ("dphi_x_min", "Dphi(l', X) > {:g} rad")]


def absent_under_selection(view, data):
    """Generators in this view with no result under the ACTIVE selection.

    Split into the two reasons, because they mean opposite things to a reader.
    An ANALYTIC row (YADISM) cannot have a FASER tier applied to it at all: the
    tiers cut on hadron-level quantities -- charged multiplicities, final-state
    angles, a second muon -- and a structure-function calculation has no events
    to cut.  Its absence from a tier table is by construction.  Anything else
    missing is simply not run yet, which is a to-do.

    Without this the two look identical: a reader moving from the inclusive tab
    to a tier tab just sees the YADISM row vanish.
    """
    analytic, unrun = [], []
    for g in view["gens"]:
        if g in data:
            continue
        try:                       # the INCLUSIVE result says which kind it is
            with open(f"{RESULTS}/histos_{at_e(g)}.json") as f:
                inc = json.load(f)
        except FileNotFoundError:
            unrun.append(g)
            continue
        (analytic if inc.get("n_fiducial") is None else unrun).append(
            inc.get("label", g))
    return analytic, unrun


def cut_line(sel):
    d = sel.as_dict()
    parts = [fmt.format(d[k]) for k, fmt in CUT_LABELS if k in d]
    unknown = set(d) - {k for k, _ in CUT_LABELS} - {"name", "label", "note"}
    if unknown:
        # a new cut in selection.py with no label here would otherwise be
        # applied to the numbers but absent from the header describing them
        parts.append("UNLABELLED CUTS: " + ", ".join(sorted(unknown)))
    return ", ".join(parts)


def load(gens):
    data = {}
    for g in gens:
        if g in data:
            continue
        try:
            # read the per-energy file, but key the dict by the PLAIN
            # generator name so the VIEWS table needs no changes
            with open(f"{RESULTS}/histos_{at_sel_e(g)}.json") as f:
                data[g] = json.load(f)
        except FileNotFoundError:
            print(f"(skipping {g}: no histos_{at_sel_e(g)}.json yet)")
    return data



def pick_ref(view, data):
    """The ratio denominator, falling back when the designated one is absent.

    A view used to be dropped WHOLE when its reference was missing, even
    though every other sample in it was present.  That cost real figures: the
    neutrino NLO view disappeared at 400 GeV and 4 TeV because Sherpa MC@NLO
    exists only at 1 TeV, taking POWHEG-V2 and Herwig NLO down with it, and
    the muon NLO view disappeared under the FASER tiers for the same reason.
    Which sample the ratio panel divides by is a PRESENTATION choice, not a
    precondition for drawing the comparison, and the panel labels its own
    denominator, so a substitution is self-documenting on the figure.

    A USABLE reference has to satisfy two things, and PRESENCE ON DISK IS
    NEITHER OF THEM:
      * it must be an MC sample -- the analytic YADISM entries carry no N_ch
        or E_lead histograms, so normalising to one would empty those panels;
      * it must have selected events and carry histograms.  A result file can
        exist and be EMPTY, and that is not hypothetical: Tier E cuts on the
        charged hadronic final state, which a matrix-element-level sample does
        not have at all, so every `*_me_*_faser_e` result is a legitimate zero.
        Testing only `in data` accepted one of those as a reference and then
        dropped all six observables one by one for want of histograms, which
        looks like six unrelated failures rather than one selection that
        cannot be applied to that sample.

    Returns (reference_key, substituted?) or (None, False) if the view has no
    usable reference -- which IS a reason to skip, and is reported as one.
    """
    def usable(g):
        d = data.get(g)
        if not d or not d.get("hists"):
            return False
        # TEST WHAT THIS VIEW ACTUALLY PLOTS, not whether the candidate is an
        # MC sample.  The old rule was `bool(d["n_fiducial"])`, and an ANALYTIC
        # entry has n_fiducial = None by construction -- so designating YADISM
        # as a ratio denominator was silently overridden, and the charm views
        # went on normalising to Sherpa while the code said FONLL.  It printed
        # a substitution line that scrolled past in a 1667-figure log; only
        # looking at the figure caught it.
        #
        # What the MC-only rule was really protecting is real: an analytic
        # entry carries no N_ch or E_lead, so it would empty those panels.
        # Requiring a histogram for every observable THIS view plots covers
        # that exactly, and lets an analytic reference serve a view that plots
        # only the six DIS observables -- which is what the charm views do.
        if not all(o in d["hists"] for o in view.get("obs", ())):
            return False
        # An MC sample with zero selected events is an empty result and must
        # not become a denominator; an analytic one legitimately has no events.
        analytic = d.get("n_parsed") is None and d.get("n_fiducial") is None
        return analytic or bool(d.get("n_fiducial"))

    if usable(view["ref"]):
        return view["ref"], False
    for g in view["gens"]:
        if usable(g):
            return g, True
    return None, False


def main():
    names = [x for a in (sys.argv[1:] or list(VIEWS))
             for x in ALIASES.get(a, [a])]
    for n in names:
        if n not in VIEWS:
            sys.exit(f"unknown view {n} (have: {', '.join(VIEWS)}"
                     f"; aliases: {', '.join(ALIASES)})")

    data = load([g for n in names for g in VIEWS[n]["gens"]])
    for g, d in data.items():
        s, e = d["sigma_fid_pb"], d["sigma_fid_err_pb"]
        nev = d.get("n_fiducial")
        tag = f"{nev} events" if nev else "analytic"
        print(f"{d['label']:24s} sigma_fid = {s/1e3:8.3f} "
              f"+- {e/1e3:.3f} nb   ({tag})")

    for n in names:
        view = VIEWS[n]
        ref_g, substituted = pick_ref(view, data)
        if ref_g is None:
            print(f"(skipping view {n}: no MC sample present to normalise to)")
            continue
        if substituted:
            print(f"(view {n}: reference {view['ref']} absent here -- ratio "
                  f"panels normalised to {ref_g} instead)")
            view = dict(view, ref=ref_g)
        print(f"\n[{n}] {view['subtitle']}")
        for key in view.get("obs") or AXES:
            plot_obs(key, data, view,
                     f"{RESULTS}/{at_sel_e(view['prefix'] + key)}.png")


if __name__ == "__main__":
    main()
