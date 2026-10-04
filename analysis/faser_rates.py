#!/usr/bin/env python3
"""Expected FASER event rates: this benchmark's cross-sections x the LHC flux.

WHAT THIS IS FOR (user, 2026-09-01): "convolute the DIS cross-sections (for
the FASER selection) with the incoming FASER muon and neutrino fluxes to make
an estimate for the expected event rates."  The neutrino side landed first;
the muon side (user, 2026-09-02) is here now, together with a comparison of
the neutrino rates against the published FASER predictions.

    N = integral dE  Phi(E) x sigma_fid(E) x T

  Phi     the flux through the detector in LHC Run 3, from the flux authors'
          own files -- NOT ours.  `data/faser_flux/` for neutrinos (150 fb^-1)
          and `data/faser_muon_flux/` for muons (250 fb^-1); each directory's
          README carries the provenance and the authors' caveats.
  sigma   THE BENCHMARK'S OWN fiducial cross-section, per generator, per
          selection tier.  This is the factor this project contributes.
  T       target nucleons per cm^2 -- measured for neutrinos, and already
          inside the muon flux normalisation (see below).

>>> THE FLUX IS A SPECTRUM, NOT THREE ENERGIES, AND THAT COST A FACTOR 1.7.
<<<  The first version of this module gave each of three energy bins the
cross-section of the benchmark energy representing it.  But the flux-averaged
neutrino energy in the lowest bin is 171 GeV, not the 400 GeV of its
representative, and sigma_fid is very nearly linear in E on the charged
current -- so that prescription overestimated the rate by 68%, invisibly,
because every generator was overestimated by the same factor and the spread
the page is about was untouched.  The muon side would have been worse:
sigma_fid(10 GeV) is a TENTH of sigma_fid(400 GeV), where a power law through
the benchmark energies says a half, because the Q2 > 4 GeV2 and 0.2 < y < 0.9
fiducial region closes as the energy falls.

So the convolution now runs over the flux's own energy points, against a
CONTINUOUS sigma_fid(E):

    sigma_g(E) = sigma_yadism_nlo(E) x R_g(E)

with the shape from the analytic reference, measured on the benchmark's own
fiducial region at 16 energies from 10 GeV to 6 TeV
(`tools/yadism_energy_ladder.sh`), and R_g -- an O(1) ratio, 0.98 to 1.09 --
interpolated in log E between each generator's three computed points and held
constant outside them.  Every generator therefore still enters through the
cross-sections the benchmark actually computed at 400, 1000 and 4000 GeV;
only the interpolation between and below them is new.  The three bins survive
as a display of WHERE the rate comes from, and the old bin-representative
total is kept in the JSON as `total_binrep` so the size of the bias stays on
the record rather than in a commit message.

>>> THE NEUTRINO TARGET IS CALIBRATED AGAINST THE FLUX AUTHORS' OWN
    INTERACTION COUNTS, NOT GUESSED FROM A DETECTOR DRAWING. <<<  The flux
    repository ships both the neutrinos passing through the detector and the
    charged-current interactions in it, and their ratio is sigma_tot(E) x T.
    Dividing by a standard neutrino-nucleon cross-section therefore MEASURES
    the column density the published FASER numbers were made with,
    1.23e27 nucleons/cm^2, which agrees with 1.1 tonnes of tungsten over a
    25 x 25 cm face to about 20%.

>>> THE MUON SIDE NEEDS NO SUCH CALIBRATION, WHICH IS WHAT HELD IT UP. <<<
    Eq. (2.1) of arXiv:2506.13889 defines the muon flux as
    f(x) = n_T L_T dN/dx in pb^-1, so the tungsten column density -- nucleon
    density times the 50 cm of FASERnu that paper uses -- is already in the
    number, and the rate is a plain integral of f against sigma in pb.

>>> THE TARGET IS TUNGSTEN, AND TUNGSTEN IS 74 PROTONS AND 110 NEUTRONS
    (user, 2026-09-04). <<<  "For the sake of the 'Predictions for FASER'
    studies we have to generate cross-sections also for neutron target, so
    that the W nucleus is correctly reproduced (for neutrino DIS xsections on
    proton and neutron are quite different)."

    They are: on the charged current the neutron is the LARGER by 79% at
    1 TeV, because nu + d -> mu- + u runs on the valence d quark and a neutron
    has two of them, so sigma_W per nucleon is 1.47 sigma_p.  The prescription
    this module carried from 2026-09-02 until then -- all 184 nucleons counted
    as protons, "neglecting any differences between protons and neutrons and
    nuclear effects, for the sake of the comparison" -- therefore understated
    every neutrino rate on the page by a third, in a way no single number
    showed, since every generator was understated by the same factor.

    The neutron cross-section is the same calculation with the isospin-
    mirrored PDF set (analysis/target.py, tools/make_neutron_pdf.py), and the
    ratio it defines is applied to every generator: see tungsten_factor() for
    what that assumes and how the assumption is checked.  The superseded
    all-proton number is kept in the JSON as `total_all_protons`, and the
    protons-alone one as `total_protons_only`.

    WHAT IS STILL NEGLECTED, and deliberately: nuclear effects.  Shadowing,
    anti-shadowing and the EMC effect are a few per cent at these x and Q2 and
    would need a nuclear PDF set, which this benchmark does not carry.

Usage:
  analysis/faser_rates.py [--current nu|mu] [--selection inclusive|faser_s|...]
  analysis/faser_rates.py --compare        # against the published predictions
Writes results_nu/faser_rates_<current>_<selection>.json and prints the table.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                          # noqa: E402

import beams                                                # noqa: E402
import labels                                               # noqa: E402
import paths                                                # noqa: E402
import selection as _sel                                    # noqa: E402
import target as _tgt                                       # noqa: E402

# Every non-default selection's filename suffix, taken FROM the selection
# module rather than listed here, so a new selection cannot be added without
# this knowing about it.  Used to stop the inclusive ladder glob from picking
# up another region's files.
SELECTION_SUFFIXES = tuple(f"_{n}" for n in _sel.SELECTIONS
                           if n != _sel.DEFAULT)

BASE = paths.REPO
FLUXDIR = f"{BASE}/data/faser_flux"
MUFLUXDIR = f"{BASE}/data/faser_muon_flux"

# THE ONE LUMINOSITY EVERY NUMBER ON THE PAGE IS QUOTED AT: 250 fb^-1, the
# total FASERnu expects by the end of Run 3, and the normalisation of both
# published predictions this module compares against.  The neutrino flux
# files are 150 fb^-1 and are scaled; the muon flux files are already 250.
LUMI_FB = 250.0
LUMI_FLUX_NU = 150.0
LUMI_FLUX_MU = 250.0

# Tungsten.  A is the MASS number: the user's prescription counts all A
# nucleons of the nucleus as protons.  Z is kept only for the proton-only
# number retained in the JSON.
A_W = 184.0
Z_W = 74.0

# The three beam energies the benchmark has, and the bin each represents.
# LOGARITHMIC boundaries, so a bin edge sits at the geometric mean of the two
# energies it separates: 632 GeV between 400 and 1000, 2000 GeV between 1000
# and 4000.  These now only GROUP the rate for display -- the cross-section
# is evaluated at the flux's own energies, not at E_REP.
E_REP = [400.0, 1000.0, 4000.0]
E_EDGES = [0.0, float(np.sqrt(400.0 * 1000.0)),
           float(np.sqrt(1000.0 * 4000.0)), np.inf]

# >>> THE COLUMN DENSITY BELONGS TO THE FLUX FILE, NOT TO THE DETECTOR WE
#     WANT TO PREDICT (2026-09-04).  GETTING THAT WRONG IS A PURE, SILENT
#     NORMALISATION ERROR, AND IT WAS HERE FOR THREE DAYS. <<<
#
# The vendored flux is a COUNT of neutrinos through the aperture its authors
# defined, not a flux density, so the rate is
#
#     N = Phi_file(E) x sigma(E) x M / A_file
#
# where A_file is the aperture the COUNT was made through and M the target
# mass we want.  A_file cancels against the detector face only if the two are
# the same, and they are not: the files are those of arXiv:2105.08270, whose
# FASERnu is |x|, |y| < 12.5 cm (25 x 25 cm) holding 1.2 tonnes, while
# arXiv:2402.13318 -- the prediction this module is compared against -- has
# the Run 3 detector, 25 x 30 cm holding 1.1 tonnes.  Dividing the Run 3 mass
# by the Run 3 FACE, as this module did, describes a detector of
# 1.1 t x 625/750 = 0.92 t behind the file's aperture: 17% of the rate, gone,
# with every cross-section on the page still perfectly correct.
#
# THE IDENTIFICATION IS NOT A GUESS.  Summed over energy, the vendored
# FASER_CCint_* files reproduce Table I of arXiv:2105.08270 to four figures --
# 1710.6 against 1710 for nu_e + nubar_e, 5782.3 against 5782 for the muon
# flavour, 40.5 against 40.5 for the tau -- so they ARE that paper's
# combination row, made with its 1.2 t over its 25 x 25 cm aperture.  The
# arithmetic is checked every run by check_flux_provenance() below.
#
# The flux is nearly flat over the aperture ("the neutrino flux is almost
# constant throughout the detector's cross sectional area", 2105.08270 Sec.
# IV), which is what lets a 25 x 30 cm detector be predicted from a 25 x 25 cm
# count at all; the residual is an acceptance shape, not a factor.
N_A = 6.02214076e23
A_W_GMOL = 183.84                     # tungsten, molar mass
FLUX_APERTURE_CM2 = 25.0 * 25.0       # 2105.08270: |x|, |y| < 12.5 cm
FLUX_TARGET_MASS_G = 1.2e6            # 2105.08270: 1.2 tonnes of tungsten
FASERNU_MASS_G = 1.1e6                # 2402.13318: 1.1 t, 25 x 30 x 80 cm

# A standard isoscalar neutrino-nucleon charged-current cross-section.  It is
# NO LONGER used to set the column density -- sigma/E is not constant over the
# flux, and the fit returned 1.229e27 against the truth of 1.157e27, 6% high
# and in the wrong direction to be noticed.  Kept for the diagnostic print.
SIGMA_CC_PER_GEV = 0.677e-38          # cm^2 per GeV per nucleon, nu
SIGMA_CC_PER_GEV_BAR = 0.334e-38      # cm^2 per GeV per nucleon, nubar

PB_TO_CM2 = 1.0e-36

# sqrt(s) of the pp collision the muon flux grid is expressed in: its variable
# is x_mu = 2 E_mu / sqrt(s_pp), so x_mu = 1 is 6.8 TeV.
SQRT_S_PP = 13600.0

FLAVOURS = [("14", r"$\nu_\mu$", SIGMA_CC_PER_GEV),
            ("-14", r"$\bar\nu_\mu$", SIGMA_CC_PER_GEV_BAR),
            ("12", r"$\nu_e$", SIGMA_CC_PER_GEV),
            ("-12", r"$\bar\nu_e$", SIGMA_CC_PER_GEV_BAR),
            ("16", r"$\nu_\tau$", SIGMA_CC_PER_GEV),
            ("-16", r"$\bar\nu_\tau$", SIGMA_CC_PER_GEV_BAR)]

# The generators to report, in the order the report shows them.  Both GENIE
# rows on each current, per CONVENTIONS.md rule 1b.
GENERATORS_NU = [
    ("POWHEG-V2", "powheg_nu"),
    ("Sherpa MC@NLO", "sherpa_nlo"),
    ("Herwig 7 (POWHEG)", "herwig_nlo_full"),
    (labels.GENIE_NU, "genie_lo"),
    (labels.GENIE_NU_NNPDF, "genie_nnpdf"),
    (labels.GENIE_NU_HEDIS, "genie"),
    ("YADISM NLO (FONLL)", "yadism_nlo_fonll_damp"),
]

# The muon counterpart, same shape and the same order of matchings
# (CONVENTIONS.md rule 2b).  There is no HEDIS row: HEDIS is a neutrino module.
GENERATORS_MU = [
    ("POWHEG-RES", "powheg"),
    ("Sherpa MC@NLO", "sherpa"),
    ("Herwig 7 (POWHEG)", "herwig_nlo_powheg_full"),
    (labels.GENIE_MU, "genie"),
    (labels.GENIE_MU_NNPDF, "genie_nnpdf"),
    ("YADISM NLO (FONLL)", "yadism_nlo_fonll_damp"),
]

# The analytic reference row, and the fallback used on a region where the
# FONLL assembly does not exist.  Named once, read by rates() and the report.
FONLL_KEY = "yadism_nlo_fonll_damp"

CURRENTS = {
    "nu": {"dir": "results_nu", "generators": GENERATORS_NU,
           "beam": r"$\nu_\mu$", "process": "nu_mu + p -> mu- + X (CC)"},
    "mu": {"dir": "results", "generators": GENERATORS_MU,
           "beam": r"$\mu^-$", "process": "mu- + p -> mu- + X (NC)"},
}

# ---------------------------------------------------------------- references
# The published predictions this benchmark's rates are set beside.  Numbers
# transcribed from the papers, with the table they come from, so a reader can
# check them and a checker can test the prose against them.
REF_NU = {
    "arxiv": "2402.13318",
    "title": "Neutrino Rate Predictions for FASER",
    "table": "Table I, combination row",
    "lumi_fb": 250.0,
    "detector": "FASERnu, 1.1 t of tungsten",
    "quantity": "CC interactions, ALL Q2 and y",
    "numu_plus_numubar": 8507.0,
    "numu_plus_numubar_err": (992.0, 962.0),      # (+, -)
    "nue_plus_nuebar": 1675.0,
    "nutau_plus_nutaubar": 28.0,
}
REF_MU = {
    "arxiv": "2506.13889",
    "title": "Deep-Inelastic Scattering at TeV Energies with LHC Muons",
    "table": "Table 2.1, NNPDF4.0 row",
    "lumi_fb": 250.0,
    "detector": "FASERnu, 25 x 30 cm, 50 cm of tungsten",
    "dis_cuts": 2.7e5,          # Q > 1.65 GeV, W > 2 GeV -- and NOTHING else
    "dis_plus_fiducial": 1.7e5,  # + E'_mu > 100 GeV, n_tracks >= 3
    # >>> ONE FACTOR SEPARATES BOTH COLUMNS FROM OURS, AND IT IS A TARGET
    #     COLUMN DENSITY (2026-09-04). <<<  The muon flux is normalised as
    #     n_T L_T dN/dx, so the target is inside it, and Sec. 2 of that paper
    #     says L_T = 50 cm -- 965 g/cm^2.  The full FASERnu is 1.1 t behind
    #     25 x 30 cm, 1467 g/cm^2, a factor 1.520 more.  Our rates times that
    #     factor reproduce Table 2.1 to 5% in BOTH columns; without it they
    #     are 0.63 and 0.68 of it.  See analysis/faser_muflux_check.py for the
    #     full argument, including why the DIS column cannot carry the
    #     E' > 100 GeV cut.  THE FACTOR IS NOT APPLIED TO ANY RATE ON THIS
    #     PAGE: our muon numbers use the flux exactly as released, which is
    #     the 50 cm the paper describes.  It is reported beside the ratio so
    #     that the comparison says what it measures.
    "full_column_factor": (1.1e6 / (25.0 * 30.0)) / (19.3 * 50.0),
}


# ------------------------------------------------------------------ the flux
def flux(pid, cc=False):
    """(E [GeV], N) from the vendored neutrino files; N is for 150 fb^-1."""
    name = f"FASER_CCint_{pid}.txt" if cc else f"FASER_{pid}.txt"
    d = np.loadtxt(f"{FLUXDIR}/{name}")
    return d[:, 0], d[:, 1]


def column_density(mass_g=FASERNU_MASS_G):
    """Target nucleons per cm^2 to pair with the vendored neutrino flux.

    The mass we want, divided by the aperture the FLUX was counted through --
    see the block on FLUX_APERTURE_CM2 above for why the detector face is the
    wrong denominator.  Default: FASERnu as arXiv:2402.13318 describes it,
    1.1 tonnes of tungsten, giving 1.061e27 nucleons/cm^2.
    """
    return mass_g / A_W_GMOL * N_A * A_W / FLUX_APERTURE_CM2


def column_density_sigma_fit(pid="14"):
    """The SUPERSEDED estimate, kept so its size stays on the record.

    The flux authors' CC count over their flux is sigma_tot(E) x T, so
    dividing by a cross-section returns T -- but only if the cross-section is
    right, and sigma = 0.677e-38 E cm^2 is the low-energy linear form, which
    the W propagator bends over well inside the flux.  Fitted over 10 GeV to
    2 TeV it returns 1.229e27 against the file's true 1.157e27.
    """
    e, phi = flux(pid)
    _e2, cc = flux(pid, cc=True)
    s = SIGMA_CC_PER_GEV if not pid.startswith("-") else SIGMA_CC_PER_GEV_BAR
    m = (e > 10.0) & (e < 2000.0) & (phi > 0)
    t = (cc[m] / phi[m]) / (s * e[m])
    return float(np.mean(t)), float(np.std(t))


# Table I of arXiv:2105.08270, "Combination (all)" and "Combination (w/o
# DPMJET)", at 150 fb^-1 over 1.2 t: the row the vendored CCint files ARE, and
# the row that shows what the DPMJET average costs.  Transcribed, then checked
# against the files themselves by check_flux_provenance().
REF_FLUX_2105 = {
    "arxiv": "2105.08270", "table": "Table I, FASERnu", "lumi_fb": 150.0,
    "target_mass_g": FLUX_TARGET_MASS_G, "aperture_cm2": FLUX_APERTURE_CM2,
    "combination_all": {"nue": 1710.0, "numu": 5782.0, "nutau": 40.5},
    "combination_no_dpmjet": {"nue": 1128.0, "numu": 5346.0, "nutau": 21.6},
}


def check_flux_provenance(tol=1e-3):
    """The vendored CCint files ARE Table I of arXiv:2105.08270.

    Four figures on three flavours is not a coincidence, and it is what fixes
    the aperture and the target mass the files were made with -- the two
    numbers the whole neutrino normalisation hangs on.  Checked rather than
    asserted in a comment, because if a future vendoring replaces the files
    with a different configuration every rate on the page moves silently.
    """
    got, want, bad = {}, REF_FLUX_2105["combination_all"], []
    for key, pids in (("nue", (12, -12)), ("numu", (14, -14)),
                      ("nutau", (16, -16))):
        got[key] = float(sum(flux(str(p), cc=True)[1].sum() for p in pids))
        if abs(got[key] / want[key] - 1.0) > tol:
            bad.append(f"{key}: files {got[key]:.1f} vs Table I {want[key]:g}")
    return got, bad


def muon_flux(geometry="25x30"):
    """(E [GeV], x f(x) summed over mu- and mu+) from the LHAPDF grid.

    The grid is `x f(x)` on 22 nodes in x_mu = 2 E_mu / sqrt(s_pp), as every
    LHAPDF grid is, and f itself is n_T L_T dN/dx in pb^-1 -- the target is
    already in it (Eq. 2.1 of 2506.13889).  Muon and antimuon are stored
    separately and summed: FASERnu muon DIS is photon exchange at
    Q2 << m_Z^2, which does not know the lepton charge.
    """
    var = {"25x30": "var2", "r9": "var1"}[geometry]
    p = f"{MUFLUXDIR}/muon_flux_FASERv_Run3_{var}_0000.dat"
    with open(p) as f:
        block = f.read().split("---")[1].strip().split("\n")
    xs = np.array([float(v) for v in block[0].split()])
    qs = np.array([float(v) for v in block[1].split()])
    fl = [int(v) for v in block[2].split()]
    vals = np.array([[float(v) for v in ln.split()]
                     for ln in block[3:] if ln.strip()])
    vals = vals.reshape(len(xs), len(qs), len(fl))
    xf = vals[:, 0, fl.index(13)] + vals[:, 0, fl.index(-13)]
    return xs * SQRT_S_PP / 2.0, xf


# ------------------------------------------------------- sigma_fid(E), fine
def _result(direc, key, selection, energy):
    """One result JSON, or None if it does not exist / is under-populated."""
    tag = beams.Beams("mu", energy).tag
    esuf = "" if energy == beams.ANCHOR_ENERGY else f"_{tag}"
    ssuf = "" if selection == "inclusive" else f"_{selection}"
    p = f"{BASE}/{direc}/histos_{key}{ssuf}{esuf}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    if d.get("stat_insufficient"):
        return None
    return d


def sigma_fid(key, energy, selection, current="nu"):
    """Fiducial cross-section in pb at one beam energy, or None."""
    d = _result(CURRENTS[current]["dir"], key, selection, energy)
    return None if d is None else d.get("sigma_fid_pb")


# CHARM PRODUCTION AS AN APPLICATION OF ITS OWN (user, 2026-09-07: the
# "Charm production in NC DIS" tab).  The same convolution with two
# substitutions: every generator key carries the benchmark's final-state
# charm tag (`<key>_charmfinal`), and the analytic shape is the FONLL charm
# ladder (`histos_yadism_charm_nlo_fonll_damp_<E>`, tools/yadism_charm_ladder.sh)
# rather than the inclusive one -- the charm fraction rises from a few per
# cent at 100 GeV to 14% at 4 TeV, so a charm rate cannot ride the
# inclusive shape.  The reference row is FONLL charm at NLO.  The tungsten
# factor is ONE: charm-initiated and gluon-initiated production see isoscalar
# partons, so a neutron and a proton give the same charm cross-section up to
# the light-quark share of the NLO real emission, and the neutron charm
# ladder has not been computed (stated in the JSON, not assumed silently).
CHARM_STEM = "histos_yadism_charm_nlo_fonll_damp"
CHARM_REF_KEY = "yadism_charm_nlo_fonll_damp"


def charm_key(key):
    """The charm-tagged result key of a generator key (the reference is
    already the charm calculation)."""
    return key if key == CHARM_REF_KEY else key + "_charmfinal"


def shape(current, sel="inclusive", tgt="p", charm=False):
    """The analytic sigma_fid(E) ladder: (E, sigma_pb), sorted.

    YADISM NLO in ZM-VFNS -- the benchmark's own inclusive convention -- on
    the SAME fiducial region as everything else, at the ladder energies plus
    the three benchmark energies.  ZM rather than FONLL because what is
    wanted here is the SHAPE in E; the two schemes differ by under a per cent
    on the inclusive rate and that difference is carried by R_g below, which
    is formed against this same ZM curve.
    """
    direc = CURRENTS[current]["dir"]
    es, ss = [], []
    import glob
    import re

    # ONLY a bare energy tag: `_fonll_damp_400GeV` and `_q2min11` are other
    # calculations and other fiducial regions, and float() on the first of
    # them would abort the whole convolution.
    #
    # THE LADDER IS PER REGION.  A non-inclusive selection has its own ladder
    # under its own suffix, and must NOT fall back to the inclusive one: the
    # 0.2 < y < 0.9 window closes as the beam energy falls, so the shape in E
    # is a property of the region and borrowing it would hide that difference
    # inside R_g, where it is then extrapolated below the lowest computed
    # point.  Selection comes BEFORE the energy tag, matching
    # analyze.result_path() and _result() above.
    #
    # A NESTED TIER RIDES THE INCLUSIVE SHAPE ON PURPOSE.  Tiers S, E and the
    # dimuon signal all keep Q2 > 4 and 0.2 < y < 0.9 and only add cuts, so
    # R_g against the inclusive curve is a genuine efficiency and no separate
    # ladder is wanted or computed.  A region that is NOT nested -- dis2506 is
    # LARGER than the inclusive region, not smaller -- must have its own, and
    # gets it by having its own ladder files on disk.  So: use the region's
    # own ladder when it exists, and fall back to the inclusive one when it
    # does not.  The fallback is what every tier has always used.
    ssuf = "" if sel == "inclusive" else f"_{sel}"
    # the fallback test is asked of THIS target's files, so a region that has
    # a proton ladder and no neutron one falls back on both nucleons together
    # rather than pairing a region's protons with the inclusive neutrons
    tsuf = _tgt.suffix(tgt)
    base_stem = CHARM_STEM if charm else "histos_yadism_nlo"
    if ssuf and not glob.glob(f"{BASE}/{direc}/{base_stem}{ssuf}{tsuf}*.json"):
        ssuf = ""
    # THE TARGET NEVER FALLS BACK.  A missing neutron ladder must abort, not
    # quietly return the proton curve: the two differ by 79% on the charged
    # current, and a tungsten rate silently built out of protons is the exact
    # shape of bug CONVENTIONS.md rule 2 is about.  target.suffix() is "" for the
    # proton, so the existing lookup is unchanged.
    stem = f"{base_stem}{ssuf}{tsuf}"
    pat = re.compile(r"^_(\d+(?:\.\d+)?)(GeV|TeV)$")
    for p in glob.glob(f"{BASE}/{direc}/{stem}*.json"):
        name = os.path.basename(p)[len(stem):-len(".json")]
        m = pat.match(name)
        if name == "":
            e = beams.ANCHOR_ENERGY
        elif m:
            e = float(m.group(1)) * (1000.0 if m.group(2) == "TeV" else 1.0)
        else:
            continue                       # _fonll_damp, _q2min11, ...
        # An INCLUSIVE ladder must not swallow the other regions' files:
        # "histos_yadism_nlo_dis2506.json" strips to "_dis2506", which the
        # energy pattern rejects, but "histos_yadism_nlo_dis2506_400GeV.json"
        # would strip to "_dis2506_400GeV" -- also rejected.  Belt and braces,
        # since a future suffix might not be so lucky.
        if not ssuf and any(k in os.path.basename(p)
                            for k in SELECTION_SUFFIXES):
            continue
        if not tsuf and "_n_" in os.path.basename(p):
            continue                       # a neutron file under the proton stem
        with open(p) as f:
            ss.append(json.load(f)["sigma_fid_pb"])
        es.append(e)
    if len(es) < 4:
        raise SystemExit(
            f"[faser_rates] only {len(es)} points on the {current} sigma(E) "
            f"ladder for target {tgt!r}, selection {sel!r}"
            f"{' (charm)' if charm else ''} -- run "
            + (f"tools/yadism_charm_ladder.sh {current}" if charm else
               f"BENCH_TARGET={tgt} tools/yadism_ladder_region.sh {sel} {current}")
            + " first.")
    o = np.argsort(es)
    return np.asarray(es)[o], np.asarray(ss)[o]


def _loglog(x, xp, fp):
    """Interpolate in log-log, extrapolating with the END SLOPE.

    Used for the analytic shape only, where the end slope is measured from
    two nearby ladder points.  The generator ratios use a CLAMPED version:
    see gen_sigma.
    """
    lx, lxp, lfp = np.log(x), np.log(xp), np.log(fp)
    out = np.interp(lx, lxp, lfp)
    lo = lx < lxp[0]
    if lo.any():
        p = (lfp[1] - lfp[0]) / (lxp[1] - lxp[0])
        out[lo] = lfp[0] + p * (lx[lo] - lxp[0])
    hi = lx > lxp[-1]
    if hi.any():
        p = (lfp[-1] - lfp[-2]) / (lxp[-1] - lxp[-2])
        out[hi] = lfp[-1] + p * (lx[hi] - lxp[-1])
    return np.exp(out)


def tungsten_factor(current, selection):
    """sigma_W per NUCLEON over sigma per PROTON, as a callable of E.

    >>> TUNGSTEN IS 74 PROTONS AND 110 NEUTRONS, AND ON THE CHARGED CURRENT
        THAT IS NOT A DETAIL (user, 2026-09-04). <<<  "For the sake of the
    'Predictions for FASER' studies we have to generate cross-sections also
    for neutron target, so that the W nucleus is correctly reproduced (for
    neutrino DIS xsections on proton and neutron are quite different)."  At
    1 TeV, YADISM NLO on the benchmark's region gives 5.4611 pb on a neutron
    against 3.0464 pb on a proton -- a factor 1.79, because nu + d -> mu- + u
    runs on the valence d quark and a neutron has two of them.  So

        (74 sigma_p + 110 sigma_n) / 184 = 1.47 sigma_p

    and counting all 184 nucleons as protons, which is what this module did
    from 2026-09-02 until today, understated the neutrino rate by a THIRD.
    The muon neutral current is milder -- the photon couples to charge, so
    the neutron is the smaller by the u/d charge ratio rather than the larger
    -- but it is the same size of correction in the other direction, and
    rule 2b means both are done.

    THE FACTOR IS TAKEN FROM THE ANALYTIC REFERENCE and applied to every
    generator, which IS an approximation and is the one to argue with: it
    assumes sigma_n/sigma_p is generator-independent.  That is a statement
    about isospin and about the PDF, not about showering or hadronisation, so
    it should hold far better than the generator spread itself -- and on the
    neutrino side it is CHECKED rather than assumed, because POWHEG-V2 has its
    own proton and neutron ladders (tools/powheg_v2_faser_ladder.sh) and
    analysis/faser_powheg_rates.py compares the two ratios.

    The neutron ladder is per REGION, like the proton one, for the same
    reason: the y window closes with falling beam energy and it does not close
    identically on the two nucleons.
    """
    # >>> BOTH NUCLEONS MUST BE READ ON THE SAME REGION, AND shape()'S OWN
    #     PER-CALL FALLBACK IS NOT ENOUGH TO GUARANTEE IT. <<<  A region with
    #     a proton ladder and no neutron one made shape() keep the region for
    #     the proton and fall back to `inclusive` for the neutron, and the
    #     "ratio" was then (inclusive neutron)/(this region's proton) -- 1.9
    #     on the muon dis2506 region, which took that comparison from 0.642 to
    #     0.340 with nothing erroring.  So the region is resolved ONCE, here,
    #     from the NEUTRON's files, and both calls are given the same one.
    import glob
    direc = CURRENTS[current]["dir"]
    sel = selection
    if sel != "inclusive" and not glob.glob(
            f"{BASE}/{direc}/histos_yadism_nlo_{sel}_n*.json"):
        sel = "inclusive"
    e_p, s_p = shape(current, sel, "p")
    e_n, s_n = shape(current, sel, "n")

    def factor(e):
        e = np.atleast_1d(np.asarray(e, dtype=float))
        r = _loglog(e, e_n, s_n) / _loglog(e, e_p, s_p)
        return (_tgt.Z_W + _tgt.N_W * r) / _tgt.A_W
    return factor


def gen_sigma(current, key, selection, charm=False):
    """sigma_g(E) for one generator, as a callable, or None if it has none.

    sigma_g(E) = sigma_yadism_nlo(E) x R_g(E), with R_g interpolated in log E
    between the benchmark's computed energies and HELD CONSTANT outside them
    -- the shape is measured, the ratio is not extrapolated.  A generator
    with a single computed energy (Sherpa under the FASER tiers) gets a
    constant R_g, which is exactly the assumption its one point supports.
    """
    e_lad, s_lad = shape(current, selection, charm=charm)
    rkey = charm_key(key) if charm else key
    pts = [(e, sigma_fid(rkey, e, selection, current)) for e in E_REP]
    pts = [(e, s) for e, s in pts if s]
    if not pts:
        return None, []
    e_g = np.array([e for e, _ in pts])
    r_g = np.array([s / float(_loglog(np.array([e]), e_lad, s_lad)[0])
                    for e, s in pts])

    def sig(e):
        e = np.atleast_1d(np.asarray(e, dtype=float))
        base = _loglog(e, e_lad, s_lad)
        if len(e_g) == 1:
            return base * r_g[0]
        return base * np.exp(np.interp(np.log(e), np.log(e_g), np.log(r_g)))
    return sig, [e for e, _ in pts]


# --------------------------------------------------------------- the rates
def flux_weights(current, pid="14", geometry="25x30"):
    """(E points, flux per point at LUMI_FB, pb -> events factor).

    The one place the flux enters a rate: rates() below and the hadron
    spectra of faser_pions.py both multiply these by sigma_fid(E) x the
    tungsten factor.  For neutrinos the points are the flux file's own and
    the factor is the column density; for muons the LHAPDF grid is sampled
    on 2000 log-spaced points and the target is already inside the flux.
    """
    if current == "nu":
        e_pts, phi = flux(pid)
        phi = phi * (LUMI_FB / LUMI_FLUX_NU)
        _got, bad = check_flux_provenance()
        if bad:
            raise SystemExit("[faser_rates] the vendored flux files no longer "
                             "reproduce Table I of arXiv:2105.08270, so the "
                             "aperture and target mass the column density "
                             "assumes are not theirs: " + "; ".join(bad))
        return e_pts, phi, PB_TO_CM2 * column_density()   # pb -> events per nucleon
    e_nodes, xf = muon_flux(geometry)
    e_pts = np.geomspace(e_nodes[0], e_nodes[-1], 2000)
    dlnx = np.gradient(np.log(e_pts))
    phi = np.exp(np.interp(np.log(e_pts), np.log(e_nodes),
                           np.log(np.maximum(xf, 1e-30)))) * dlnx
    phi[np.interp(np.log(e_pts), np.log(e_nodes), xf) <= 0] = 0.0
    return e_pts, phi, 1.0                              # the target is in the flux


def _bins_of(e):
    """Index of the display bin each energy falls in."""
    return np.searchsorted(np.array(E_EDGES[1:-1]), e, side="right")


def rates(current="nu", selection="inclusive", pid="14", geometry="25x30",
          charm=False):
    """Expected number of fiducial events per generator, LHC Run 3.
    charm=True: charm-production events (see CHARM_STEM above)."""
    cfg = CURRENTS[current]
    # TUNGSTEN IS 74 p + 110 n (user, 2026-09-04).  The factor is the tungsten
    # cross-section per nucleon over the proton one, measured on this region
    # by the analytic reference -- see tungsten_factor().  It replaces the
    # 2026-09-02 prescription of counting all 184 nucleons as protons, whose
    # answer is kept in the JSON as `total_all_protons` so the size of the
    # change stays on the record rather than in a commit message.
    if charm:
        def wfac(e):                       # isoscalar: see CHARM_STEM above
            return np.ones_like(np.atleast_1d(np.asarray(e, dtype=float)))
    else:
        wfac = tungsten_factor(current, selection)
    # THE ANALYTIC ROW IS FONLL WHERE FONLL EXISTS, AND ZM WHERE IT DOES NOT.
    # The benchmark's convention is FONLL wherever charm appears, and the
    # inclusive FONLL reference is ASSEMBLED from three runs (ZM inclusive, ZM
    # charm, FONLL charm) rather than computed directly -- so on a foreign
    # region that has only the ZM ladder there is no FONLL row to show.
    # Falling back to the ZM curve, which is the one the whole convolution
    # rides anyway, is honest; falling back SILENTLY is not, so the key and
    # the label are both written into the JSON and the report prints them.
    # The fallback only happens when the ZM row can actually be filled: under
    # a FASER tier NEITHER analytic row exists, because the tier cuts live in
    # the event analysis, and swapping one empty row for another empty row
    # under a different name would be a change of label with no change of
    # content -- which is worse than leaving it alone.
    gens = list(cfg["generators"])
    if charm:
        gens = [(lab, key) for lab, key in gens if key != FONLL_KEY]
        gens.append(("YADISM NLO (FONLL) charm", CHARM_REF_KEY))
    if (not charm and not any(sigma_fid(FONLL_KEY, e, selection, current) for e in E_REP)
            and all(sigma_fid("yadism_nlo", e, selection, current)
                    for e in E_REP)):
        gens = [g for g in gens if g[1] != FONLL_KEY]
        gens.append(("YADISM NLO (ZM-VFNS)", "yadism_nlo"))
    out = {"current": current, "selection": selection, "lumi_fb": LUMI_FB,
           "charm": charm,
           "charm_note": ("charm-tagged rates: every generator key carries the "
                          "final-state charm tag, the shape is the FONLL charm "
                          "ladder, and the tungsten factor is ONE (isoscalar "
                          "charm and gluon; no neutron charm ladder computed)"
                          if charm else None),
           "e_edges": [None if not np.isfinite(x) else x for x in E_EDGES],
           "e_rep": E_REP,
           "target": f"tungsten, {_tgt.Z_W} p + {_tgt.N_W} n per nucleus; "
                     "the neutron cross-section is the isospin-mirrored PDF "
                     "(analysis/target.py)",
           "tungsten_factor_at_e_rep": [float(wfac(e)[0]) for e in E_REP],
           "reference_key": gens[-1][1], "reference_label": gens[-1][0],
           "generators": {}}

    e_pts, phi, norm = flux_weights(current, pid, geometry)
    if current == "nu":
        t = column_density()
        t_fit, t_fit_rms = column_density_sigma_fit(pid)
        out.update({"pid": pid, "column_density_per_cm2": t,
                    "column_density_target_mass_g": FASERNU_MASS_G,
                    "column_density_aperture_cm2": FLUX_APERTURE_CM2,
                    "column_density_sigma_fit": t_fit,
                    "column_density_sigma_fit_rms": t_fit_rms,
                    "flux_note": f"vendored at {LUMI_FLUX_NU:g} fb^-1, "
                                 f"scaled to {LUMI_FB:g}"})
    else:
        # integrate f(x) sigma dx = (x f) sigma dln x on a fine grid
        out.update({"geometry": geometry,
                    "flux_note": f"vendored at {LUMI_FLUX_MU:g} fb^-1, "
                                 "target column density included"})

    ib = _bins_of(e_pts)
    out["flux_per_bin"] = [float(phi[ib == i].sum()) for i in range(3)]
    out["flux_total"] = float(phi.sum())
    out["flux_mean_energy_per_bin"] = [
        float((phi[ib == i] * e_pts[ib == i]).sum() / phi[ib == i].sum())
        if phi[ib == i].sum() else None for i in range(3)]

    for label, key in gens:
        sig, e_have = gen_sigma(current, key, selection, charm=charm)
        rkey = charm_key(key) if charm else key
        if sig is None:
            out["generators"][key] = {"label": label, "per_bin": [None] * 3,
                                      "total": None, "energies": []}
            continue
        sig_p = sig(e_pts)
        n = phi * sig_p * wfac(e_pts) * norm
        # A BIN IS ONLY FILLED IF THE GENERATOR COMPUTED ITS OWN ENERGY.
        # Sherpa's FASER-tier results exist at 1 TeV alone, and a tier
        # acceptance is strongly energy dependent -- 99% at 1 TeV against 45%
        # at 400 GeV on the neutrino side -- so carrying its single ratio
        # across the whole spectrum does not produce a conservative estimate,
        # it produces a wrong one, larger than every generator that HAS the
        # samples.  A missing bin stays missing, and without all three there
        # is no total.
        have = [sigma_fid(rkey, e, selection, current) is not None
                for e in E_REP]
        per_bin = [float(n[ib == i].sum()) if have[i] else None
                   for i in range(3)]
        # what the superseded bin-representative prescription would have said
        binrep = 0.0
        for i, e in enumerate(E_REP):
            s = sigma_fid(rkey, e, selection, current)
            if s:
                binrep += (float(phi[ib == i].sum()) * s
                           * float(wfac(e)[0]) * norm)
        total = float(n.sum()) if all(have) else None
        out["generators"][key] = {
            "label": label, "per_bin": per_bin, "total": total,
            "total_binrep": binrep if all(have) else None,
            "energies": e_have,
            # the SUPERSEDED prescription, all 184 nucleons counted as
            # protons, kept so the size of the tungsten correction is visible
            "total_all_protons": (None if not all(have)
                                  else float((phi * sig_p * norm).sum())),
            # and the protons alone, which is what the benchmark computes
            "total_protons_only": (None if not all(have)
                                   else float((phi * sig_p * norm).sum())
                                   * Z_W / A_W)}
    return out


# ----------------------------------------------------- published comparison
def compare_nu(selection="inclusive", pid="14"):
    """Our neutrino rates against arXiv:2402.13318, Table I.

    THREE THINGS HAVE TO BE MADE COMMENSURATE, and each is a factor:

      luminosity  their 250 fb^-1 against the flux files' 150 -- done in
                  `rates`, so every number on the page is at 250;
      target      their tungsten against our nucleon -- 74 p + 110 n, also
                  already in `rates` (see tungsten_factor);
      flavour     their nu_mu + nubar_mu against our nu_mu.  The split is
                  taken from the flux authors' OWN interaction counts, the
                  same files the column density is calibrated on, so it is
                  not a number we invented either.

    What is left over after all three is NOT a discrepancy: their number is
    every CC interaction and ours is the fiducial region Q2 > 4 GeV2,
    0.2 < y < 0.9.  The ratio is therefore an ACCEPTANCE, and the honest way
    to read this comparison is as a measurement of it.

    AND IT IS AN ACCEPTANCE DOMINATED BY THE y WINDOW, WHICH IS WHY THERE IS
    NOW A SECOND FUNCTION (user, 2026-09-04): compare_nu_like_for_like() runs
    the same convolution on the `ally` region -- the benchmark's Q2 floor with
    NO y window -- against a reference of our own that is measured the same
    way.  Read that one for the physics; read this one for the size of the
    benchmark's fiducial region under the FASER flux.
    """
    r = rates("nu", selection, pid)
    _e, cc14 = flux("14", cc=True)
    _e, cc14b = flux("-14", cc=True)
    f_numu = float(cc14.sum() / (cc14.sum() + cc14b.sum()))
    scale = LUMI_FB / LUMI_FLUX_NU
    files_numu = float(cc14.sum() * scale)
    files_both = float((cc14.sum() + cc14b.sum()) * scale)
    ref_numu = REF_NU["numu_plus_numubar"] * f_numu
    out = {"selection": selection, "lumi_fb": LUMI_FB, "reference": REF_NU,
           "numu_fraction_of_cc": f_numu,
           "reference_numu_only": ref_numu,
           "flux_files_cc_numu": files_numu,
           "flux_files_cc_numu_plus_bar": files_both,
           "flux_files_over_paper": files_both / REF_NU["numu_plus_numubar"],
           "generators": {}}
    for key, g in r["generators"].items():
        if g["total"] is None:
            continue
        out["generators"][key] = {
            "label": g["label"], "total": g["total"],
            "over_reference_numu": g["total"] / ref_numu,
            "over_flux_files_cc": g["total"] / files_numu,
            "total_binrep": g.get("total_binrep")}
    return out


def compare_nu_like_for_like():
    """Our neutrino rate against arXiv:2402.13318 WITHOUT the y window.

    WHY THIS EXISTS (user, 2026-09-04).  "In comparing with the numbers of
    2402.13318, we should not be imposing the 0.2 < y < 0.9 cut at all, but
    rather going fully inclusive.  We can still keep the Q2 > 4 GeV2 cut to
    make sure MCs are well behaved."  compare_nu() divides our fiducial rate
    by their number and the two are measured on different regions: theirs is
    every CC interaction at every Q2 and y, ours had a y window that discards
    30% of the rate above Q2 > 4.  That ratio measures the window.  Here the
    same convolution runs on the `ally` region, so what is left is the Q2 > 4
    floor alone -- and that floor is MEASURED rather than assumed, because
    GENIE supplies both sides of it.

    THE CHAIN, EACH LINK CHECKABLE SEPARATELY:

      1. GENIE total CC x this flux x this target reproduces Table I to 1.03
         on nu_mu + nubar_mu (analysis/faser_genie_rates.py, which also shows
         the cross-section itself closing on the flux authors' own to ~1%).
      2. The reference here is that same GENIE total-CC number for nu_mu
         ALONE -- our own, from the splines -- so no flavour split has to be
         invented and the flux, the luminosity and the tungsten are literally
         the same arrays on both sides.
      3. Each generator's `ally` rate over it is then the Q2 > 4 ACCEPTANCE of
         the total charged-current rate under the FASER flux, and nothing
         else.  GENIE appears on both sides, so its own ratio is that
         acceptance measured within one generator -- the number the others
         are read against.

    WHAT CANNOT APPEAR HERE, and it is a property of the samples rather than
    an oversight: Sherpa and Herwig apply 0.2 < y < 0.9 at GENERATION level on
    the neutrino side as they do on the muon side (fiducial/parsed 0.99998 and
    0.99974), so they hold no events outside the window at all and cannot be
    asked for this region.  POWHEG-V2, Pythia and all three GENIE rows can:
    their fiducial fractions are 0.665, 0.707 and 0.67-0.69.
    """
    r = rates("nu", "ally")
    ref = f"{BASE}/results_nu/faser_genie_rates.json"
    if not os.path.exists(ref):
        raise SystemExit("[faser_rates] run analysis/faser_genie_rates.py "
                         "first -- it supplies the total-CC reference")
    with open(ref) as f:
        g = json.load(f)
    numu_total = g["flavours"]["14"]["events"]
    out = {"selection": "ally", "lumi_fb": LUMI_FB, "reference": REF_NU,
           "region": "Q2 > 4 GeV2, no y window; the reference is TOTAL CC",
           "reference_numu_total_cc": numu_total,
           "reference_source": "analysis/faser_genie_rates.py, GENIE "
                               "G18_02a total CC, nu_mu only, same flux and "
                               "same tungsten",
           "genie_total_cc_over_paper": g["numu"]["ours_over_paper"],
           "benchmark_region_for_contrast": compare_nu("inclusive"),
           "generators": {}, "absent": {
               "sherpa_nlo": "0.2 < y < 0.9 applied at generation level",
               "herwig_nlo_full": "0.2 < y < 0.9 applied at generation level"}}
    for key, gg in r["generators"].items():
        if gg["total"] is None:
            continue
        out["generators"][key] = {
            "label": gg["label"], "total": gg["total"],
            "q2_acceptance": gg["total"] / numu_total,
            "total_all_protons": gg.get("total_all_protons")}
    return out


def compare_mu(selection="inclusive"):
    """Our muon rates against arXiv:2506.13889, Table 2.1.

    Commensurate already in luminosity (both 250 fb^-1), in target (the flux
    normalisation carries n_T L_T over tungsten NUCLEONS) and in beam (both
    sum mu- and mu+, which photon exchange cannot tell apart).  What differs
    is the SELECTION, and it differs a lot: that paper's DIS cuts are
    Q > 1.65 GeV and W > 2 GeV, and its fiducial cuts add E'_mu > 100 GeV and
    at least three charged tracks, while this benchmark's fiducial region is
    Q2 > 4 GeV2 and 0.2 < y < 0.9.  The y window is the severe one -- muon
    DIS at FASER is dominated by small y -- so our rate is expected to be a
    fraction of theirs, and the ratio measures that fraction.
    """
    r = rates("mu", selection)
    out = {"selection": selection, "lumi_fb": LUMI_FB, "reference": REF_MU,
           "generators": {}}
    for key, g in r["generators"].items():
        if g["total"] is None:
            continue
        out["generators"][key] = {
            "label": g["label"], "total": g["total"],
            "over_reference_dis": g["total"] / REF_MU["dis_cuts"],
            "over_reference_dis_full_column": (g["total"]
                                              * REF_MU["full_column_factor"]
                                              / REF_MU["dis_cuts"]),
            "over_reference_fiducial": g["total"] / REF_MU["dis_plus_fiducial"],
            "total_binrep": g.get("total_binrep")}
    return out


def compare_mu_like_for_like():
    """Our muon rate against arXiv:2506.13889 ON THEIR OWN FIDUCIAL REGION.

    WHY THIS EXISTS (user, 2026-09-03).  compare_mu() divides OUR fiducial
    rate by THEIR number, and the two are measured on different regions: ours
    is Q2 > 4 GeV2 with 0.2 < y < 0.9, theirs is Q > 1.65 GeV with W > 2 GeV.
    The y window alone discards 66% of the muon rate above Q2 > 4, so that
    ratio has been reading 0.146 and measuring mostly the region mismatch.
    Here the same convolution is run on the dis2506 selection -- their region,
    no y window -- so the ratio measures the physics instead.  The benchmark's
    own fiducial region is untouched; this is the one place dis2506 is read.

    ONLY POWHEG-RES CAN ANSWER IT, and that is a property of the samples, not
    an oversight:

      Sherpa and Herwig apply 0.2 < y < 0.9 at GENERATION level (their
      fiducial/parsed ratios are 1.000 and 0.998), so they hold no events
      outside it at all.

      Pythia's muon cards set PhaseSpace:Q2Min = 4, and GENIE's sample starts
      at Q2 = 3.5; this region begins at Q2 = 2.72.  Asked for it anyway they
      return 114.0 nb and 105.3 nb against POWHEG's 168.0 -- entirely
      plausible numbers that are silently measuring a smaller region, with
      0.25% and 0.00% of their rate in the first Q2 bin against POWHEG's
      25.6%.  Their results were computed, inspected and DELETED rather than
      shown, and this is the "say so in place" that CONVENTIONS.md rule 1b asks
      for when GENIE genuinely cannot appear.

    THE NUMBER CARRIES TWO KNOWN BIASES, both toward under-counting:

      POWHEG-RES generates at Q2 > 2.25 and this region starts at 2.72, with
      25.6% of its rate in the first bin -- so events that should migrate up
      from below the generation cut are missing.  Measured size: POWHEG sits
      at 0.970 +- 0.001 of the analytic reference here against 0.985-0.998 on
      the benchmark region, i.e. about two points low, consistently at all
      three energies.

      The analytic ladder this rides is less reliable in this region than in
      the benchmark's: yadism's three independent integrators agree to six
      figures on the benchmark region and spread by 0.72% at 1 TeV rising to
      8.1% at 10 GeV here, because the region reaches x -> 0.99 where the
      calculation carries no target-mass corrections.  R_g divides out the
      absolute normalisation at the three computed energies, so what survives
      is the shape error between them, ~0.3% over 100-400 GeV.

    The independent check is the flux-mode reproduction
    (analysis/faser_muflux_check.py), which generates POWHEG-RES with the
    FASERnu muon flux as the beam and needs no analytic ladder at all: it
    gives 168800 events on the same selection against this route's 164200,
    agreeing to 2.8%.

    AND THE RATIO IS NOT THE WHOLE STORY EITHER (2026-09-04): both this route
    and the flux-mode one land at 0.92-0.95 of Table 2.1 once the target
    column is taken from the 50 cm the muon flux normalisation is stated to
    carry to the full FASERnu, a factor REF_MU["full_column_factor"] = 1.520.
    That factor is NOT applied to any rate here -- see the note on REF_MU and
    analysis/faser_muflux_check.py for the argument.
    """
    r = rates("mu", "dis2506")
    out = {"selection": "dis2506", "lumi_fb": LUMI_FB, "reference": REF_MU,
           "region": "Q > 1.65 GeV, W > 2 GeV (arXiv:2506.13889), no y window",
           "benchmark_region_for_contrast": compare_mu("inclusive"),
           "generators": {}, "absent": {
               "sherpa": "0.2 < y < 0.9 applied at generation level",
               "sherpa_lo": "0.2 < y < 0.9 applied at generation level",
               "herwig": "0.2 < y < 0.9 applied at generation level",
               "herwig_nlo_powheg": "0.2 < y < 0.9 applied at generation level",
               "pythia": "PhaseSpace:Q2Min = 4 at generation level; region "
                         "starts at Q2 = 2.72",
               "genie": "sample starts at Q2 = 3.5; region starts at 2.72"}}
    for key, g in r["generators"].items():
        if g["total"] is None:
            continue
        out["generators"][key] = {
            "label": g["label"], "total": g["total"],
            "over_reference_dis": g["total"] / REF_MU["dis_cuts"],
            "over_reference_dis_full_column": (g["total"]
                                              * REF_MU["full_column_factor"]
                                              / REF_MU["dis_cuts"]),
            "total_binrep": g.get("total_binrep")}
    return out


# ------------------------------------------------------------------- output
def _print(r):
    cur = r["current"]
    print(f"current {cur}, selection {r['selection']}, {LUMI_FB:g} fb^-1, "
          f"tungsten ({_tgt.Z_W} p + {_tgt.N_W} n; sigma_W/nucleon over "
          "sigma_p = "
          + ", ".join(f"{w:.3f}" for w in r["tungsten_factor_at_e_rep"])
          + " at 400/1000/4000 GeV)")
    if cur == "nu":
        print(f"  target column density "
              f"{r['column_density_per_cm2']:.3e} nucleons/cm2 = "
              f"{r['column_density_target_mass_g']/1e6:.2g} t of tungsten "
              f"behind the flux files' own "
              f"{r['column_density_aperture_cm2']:.0f} cm2 aperture "
              f"(the superseded sigma fit said "
              f"{r['column_density_sigma_fit']:.3e})")
    print("  flux per bin: "
          + ", ".join(f"{x:.4g}" for x in r["flux_per_bin"])
          + f"   (total {r['flux_total']:.4g})")
    print("  flux-averaged energy per bin: "
          + ", ".join("n/a" if x is None else f"{x:.0f} GeV"
                      for x in r["flux_mean_energy_per_bin"]))
    print(f"{'generator':28s} {'<632':>10s} {'632-2000':>10s} {'>2000':>10s}"
          f" {'total':>10s} {'binrep':>10s}")
    for key, g in r["generators"].items():
        cells = " ".join("       n/a" if v is None else f"{v:10.1f}"
                         for v in g["per_bin"])
        tot = "n/a" if g["total"] is None else f"{g['total']:10.1f}"
        br = ("n/a" if not g.get("total_binrep")
              else f"{g['total_binrep']:10.1f}")
        print(f"{g['label']:28s} {cells} {tot} {br}")


def main():
    argv = sys.argv[1:]
    if "--compare" in argv:
        cn, cm = compare_nu(), compare_mu()
        print("=== neutrinos vs arXiv:2402.13318 Table I "
              f"({REF_NU['numu_plus_numubar']:.0f} nu_mu + nubar_mu CC "
              f"at {REF_NU['lumi_fb']:g} fb^-1)")
        print(f"  nu_mu is {cn['numu_fraction_of_cc']:.3f} of the CC rate "
              f"-> reference for nu_mu alone {cn['reference_numu_only']:.0f}")
        print(f"  the vendored flux files' own CC count is "
              f"{cn['flux_files_cc_numu_plus_bar']:.0f}, "
              f"{cn['flux_files_over_paper']:.2f}x the paper's")
        for key, g in cn["generators"].items():
            print(f"    {g['label']:28s} {g['total']:8.0f}  "
                  f"{g['over_reference_numu']:.3f} of the reference, "
                  f"{g['over_flux_files_cc']:.3f} of the flux files' CC")
        print("=== muons vs arXiv:2506.13889 Table 2.1 "
              f"({REF_MU['dis_cuts']:.2g} with DIS cuts, "
              f"{REF_MU['dis_plus_fiducial']:.2g} with fiducial cuts)")
        for key, g in cm["generators"].items():
            print(f"    {g['label']:28s} {g['total']:8.0f}  "
                  f"{g['over_reference_dis']:.3f} of the DIS-cut rate")
        # THE SAME MUON COMPARISON ON THEIR REGION, which is the one that
        # measures physics rather than the difference between two fiducial
        # definitions.  See compare_mu_like_for_like() for why only
        # POWHEG-RES can answer it and what biases the answer carries.
        cl = compare_mu_like_for_like()
        print("=== muons vs the SAME reference on ITS OWN region "
              "(Q > 1.65 GeV, W > 2 GeV, no y window)")
        print(f"    (x{REF_MU['full_column_factor']:.3f} takes the flux's "
              f"50 cm column to the full FASERnu one -- see "
              f"analysis/faser_muflux_check.py)")
        for key, g in cl["generators"].items():
            print(f"    {g['label']:28s} {g['total']:8.0f}  "
                  f"{g['over_reference_dis']:.3f} of the DIS-cut rate, "
                  f"{g['over_reference_dis_full_column']:.3f} on the full "
                  f"column")
        for key, why in sorted(cl["absent"].items()):
            print(f"    {key:28s}    absent: {why}")
        # AND THE NEUTRINO COUNTERPART OF THE SAME MOVE (rule 2b): the
        # comparison on a region with no y window, against a reference
        # measured the same way.  See compare_nu_like_for_like().
        nl = compare_nu_like_for_like()
        print("=== neutrinos vs the TOTAL CC rate on a region with no y "
              "window (Q2 > 4 GeV2, all y), nu_mu alone")
        print(f"  reference: our own GENIE total CC, "
              f"{nl['reference_numu_total_cc']:.0f} events, which is "
              f"{nl['genie_total_cc_over_paper']:.3f} of "
              f"{REF_NU['arxiv']} Table I for nu_mu + nubar_mu")
        for key, g in nl["generators"].items():
            print(f"    {g['label']:28s} {g['total']:8.0f}  "
                  f"Q2 > 4 acceptance {g['q2_acceptance']:.3f}")
        for key, why in sorted(nl["absent"].items()):
            print(f"    {key:28s}    absent: {why}")
        p = f"{BASE}/results_nu/faser_rates_comparison.json"
        with open(p, "w") as f:
            json.dump({"nu": cn, "mu": cm, "mu_like_for_like": cl,
                       "nu_like_for_like": nl}, f, indent=1)
        print("wrote", p)
        return

    cur = "nu"
    sel = "inclusive"
    if "--current" in argv:
        cur = argv[argv.index("--current") + 1]
    if "--selection" in argv:
        sel = argv[argv.index("--selection") + 1]
    charm = "--charm" in argv
    r = rates(cur, sel, charm=charm)
    _print(r)
    p = f"{BASE}/results_nu/faser_rates_{cur}{'_charm' if charm else ''}_{sel}.json"
    with open(p, "w") as f:
        json.dump(r, f, indent=1)
    print("wrote", p)


if __name__ == "__main__":
    main()
