#!/usr/bin/env python3
"""Single-inclusive charged pion (and kaon) production at FASERnu.

WHAT THIS IS FOR (user, 2026-09-07): "FASER event rate predictions, using
the emulsion detector, for single pion production" -- the TODO entry that
asked for it says why: "inclusive single pion and kaon production ... to test
QCD models of inclusive hadron production, fragmentation functions".  So
this is SINGLE-INCLUSIVE hadron production in DIS, l N -> l' h X, the
observable fragmentation functions are fitted to, and not resonant
single-pion production (which the non-DIS study measured to be per mille of
the rate and below 0.02% of the fiducial one).

>>> THE V2 REGION, AND NOTHING ELSE (user, 2026-09-18). <<<  "Move the pion
generator slide to the no-y-cut region, redo simulation and plots for SIDIS
without any y cut, only the fiducial cuts."  So every spectrum here is
read from the paper-plots samples -- tungsten by isospin, each nucleon
generated separately (sample_layout), Q2 > 4 GeV2 and W > 3 GeV at generation, no
y window -- and the four regions below are the region and its subsets:

    sidis_w      Q2 > 4, W > 3                     == q4w3          (the region)
    sidis_e      + FASER Tier E                    == q4w3_faser_e  (the YIELDS)
    sidis_ex     + x > 0.1                          diagnostic only
    sidis_paper  Q2 > 4, W > 3, x > 0.1, no tier    the arXiv:2504.05376 region

The first two are the selections under the names the SIDIS study has
used since 2026-09-07 (selection.py states the equivalence); the JSON records
both names.  The earlier regions with a y window (`inclusive`, `faser_e`, `ally`)
are gone from here: a sample is generated with W > 3 GeV, so an all-y
region read from it would be silently truncated, which is exactly the trap
this benchmark keeps meeting.

WHAT IS COMPUTED.  For every DIS event of a sample, every charged pion
(|PDG| = 211) and charged kaon (321) in the final state -- both are STABLE
under the benchmark's ctau > 10 mm convention, so these are the primary
hadrons the emulsion sees at the vertex -- and for each its lab energy
E_h = p.P/M_p, its energy fraction z = E_h / nu with nu = E_beam - E'_lepton
the energy of the hadronic system, its transverse momentum relative to the
beam axis, and whether it lies inside tan(theta) < 0.5, the angular
acceptance of the emulsion tracks used by the Tier E selection.  Spectra
are filled PER SELECTED EVENT, so a sample's absolute normalisation never
enters the extraction: it comes from the cross-section ladders below.

THE SPECIES.  "pi" and "K" are the charge sums; the four signed ones carry
the flavour message of neutrino SIDIS: nu p -> mu- pi+ X runs on the valence
d -> u transition and is "favoured", nu p -> mu- pi- X needs a sea flavour
combination and is "unfavoured", and a charge-summed spectrum averages that
away.

THE TARGET is assembled the way every result is (analysis/
combine_target.py): a per-event spectrum on tungsten is the proton and
neutron spectra weighted by 74 sigma_p and 110 sigma_n OF THE REGION, so a
region where the two nucleons differ (Tier E asks for a hadronic system) is
weighted by what each nucleon contributes to it, and the region's tungsten
cross-section per nucleon, (74 sigma_p + 110 sigma_n)/184, comes out of the
same combination.  sigma_p and sigma_n of the region are the tracked
histos_<gen>_q4w3_{p,n}[_TAG].json; a subset's is that times the fraction
this pass measures, so the regions nest by construction.

THE RATE.  Pions per bin at FASER in Run 3:

    N_h(bin) = L x integral dE  Phi(E) sigma^W_region(E) T  x  n_h(bin | E)

with sigma^W_region(E) the region's tungsten cross-section per nucleon at
the energies of the SIDIS ladder (beams.SIDIS_ENERGIES: 300, 400, 700,
1000, 2000, 4000 GeV; POWHEG-V2 also at 20, 50, 100 and 200 GeV,
beams.SIDIS_ENERGIES_LOW, for the region-only yield) carried across the flux on the YADISM (ZM-VFNS, TMC)
shape of the region -- the generator's ratio to it interpolated in log E
between the ladder points and held constant outside -- and n_h(bin | E) the
per-event spectrum, interpolated in log E and held constant outside.  The
Tier E regions are the exception below the ladder: their efficiency is
EXACTLY zero at and below 200 GeV (the tier asks for a scattered lepton
above 200 GeV), and between 200 and 300 GeV it is taken linear in log E
from that zero to the measured 300 GeV point.  The 300 GeV point exists for
this reason: a third of the flux-weighted rate of the region sits below
400 GeV, where the tier's efficiency falls from 0.15 to 0.04.  Quoted at
300 fb^-1, the full Run 3, as the dimuon cut flow is.

BOTH CURRENTS (CONVENTIONS.md rule 2b), GENIE ALWAYS (rule 1b): the neutrino CC
samples of POWHEG-V2, Herwig 7 (POWHEG matching, positive half), Sherpa
MC@NLO and the two GENIE tunes, and the muon NC ones of POWHEG-RES, Herwig,
Sherpa and GENIE.  What is DRAWN is the neutrino current alone (user,
2026-09-07, the one carve-out from rule 2b); the muon numbers are computed
and written beside them.

Usage:
  faser_pions.py extract --current nu|mu --key <key> --target p|n --energy E
      -> results{,_nu}/pions_<key>_q4w3_<t>[_TAG].json   (one nucleon)
  faser_pions.py combine --current nu|mu --key <key> --energy E
      -> results{,_nu}/pions_<key>_q4w3_W[_TAG].json     (tungsten per nucleon)
  faser_pions.py rates
      -> results_nu/faser_pions.json  (pions and kaons at FASER, both currents)
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import analyze                                               # noqa: E402
import analyze_nu                                            # noqa: E402
import beams                                                 # noqa: E402
import faser_rates as fr                                     # noqa: E402
import paths                                                 # noqa: E402
import selection as selmod                                   # noqa: E402
import target                                                # noqa: E402
import sample_layout                                              # noqa: E402
from analyze import (M_P, CHARGED_HADRONS, dot, dis_invariants,  # noqa: E402
                     lab_energy, lepton_theta, job_files)

TRACK_E_MIN = selmod.TRACK_E_MIN   # GeV: a Tier E track (2026-09-19)

# >>> THE HADRON SELECTION EFFICIENCY (user, 2026-09-21). <<<  "I think a
# selection efficiency of 60% would be realistic in both cases" -- pions AND
# kaons.  Until then every yield here was a GENERATED rate inside the
# acceptance: the only cuts in the chain are geometric and kinematic (Tier E
# on the event; tan theta < 0.5 and E_h > TRACK_E_MIN on the hadron), with no
# reconstruction or identification factor anywhere.
#
# IT MULTIPLIES COUNTS, NOT MULTIPLICITIES.  Every yield at LUMI_FB carries
# it; the per-event multiplicities, the pion-presence fractions and the
# multiplicity compared against the NNLO calculation do NOT, because those
# are properties of the events rather than of the measurement, and folding a
# detector factor into them would break the comparison the NNLO figure makes.
#
# >>> AND THE ERROR DOES NOT SCALE AS 1/sqrt(EFF). <<<  Keeping each hadron
# independently with probability eps DECORRELATES the hadrons within an
# event, so the compound-Poisson inflation shrinks toward one.  With the
# event count Poisson and the kept hadrons Binomial(n, eps),
#
#     y   -> eps * y            (y   = lambda <n>)
#     var -> eps^2 * var + eps * (1 - eps) * y      (var = lambda <n^2>)
#
# which is what hadron_yield() below applies.  Numerically the clustering
# factor falls from 1.24 to 1.08 on pions, so a 60% efficiency costs less
# statistical error than the naive sqrt would say -- and quoting "1.24 times
# the naive sqrt(N)" alongside an efficiency-corrected yield would be wrong.
HADRON_SELECTION_EFF = 0.6


def hadron_yield(y, var=None):
    """A generated hadron count (and its compound-Poisson variance) as the
    detector would select it.  Pass var=None for a count with no error."""
    y = np.asarray(y, dtype=float)
    ye = HADRON_SELECTION_EFF * y
    if var is None:
        return ye
    var = np.asarray(var, dtype=float)
    return ye, (HADRON_SELECTION_EFF ** 2 * var
                + HADRON_SELECTION_EFF * (1.0 - HADRON_SELECTION_EFF) * y)


LUMI_FB = 300.0
REGION = sample_layout.REGION                         # q4w3
# the SIDIS ladder, plus POWHEG-V2's points below 300 GeV (2026-09-21,
# beams.SIDIS_ENERGIES_LOW); a generator without them simply has no
# spectra there and keeps the held-constant extrapolation
ENERGIES = beams.SIDIS_ENERGIES_LOW + beams.SIDIS_ENERGIES

# THE REGIONS: the region and its subsets, under the SIDIS study's names.
SELS = {"sidis_w": selmod.SIDIS_W, "sidis_e": selmod.SIDIS_E_W3,
        "sidis_ex": selmod.SIDIS_E_W3_X, "sidis_paper": selmod.SIDIS_PAPER}
# the selection each one IS (same cuts, selection.py)
REGION_ALIAS = {"sidis_w": "q4w3", "sidis_e": "q4w3_faser_e"}
# THE ENCLOSING REGION each one is normalised through: a subset's rate is
# its parent's times the fraction this pass measures, so a subset can never
# come out larger than the set.  sidis_w is the root: it IS the region
# and has the tracked cross-section ladder.
RATE_PARENT = {"sidis_w": None, "sidis_e": "sidis_w", "sidis_ex": "sidis_e",
               "sidis_paper": "sidis_w"}
# the regions whose efficiency is exactly zero at E_beam <= 200 GeV
TIER_E_REGIONS = ("sidis_e", "sidis_ex")
ROOT = "sidis_w"

HADRONS = {"pi": 211, "K": 321}
SIGNED = {"pip": 211, "pim": -211, "Kp": 321, "Km": -321}
# "h" IS EVERY STABLE CHARGED HADRON (user, 2026-09-30: "FASERnu cannot tell
# apart pions from kaons ... the sum of all charged hadrons"): analyze's
# CHARGED_HADRONS, pi K p Sigma+- Xi- Omega- and their antiparticles -- the
# same set Tier E counts tracks with -- charge summed, since the emulsion does
# not measure the charge either.
SPECIES = ["pi", "K", "pip", "pim", "Kp", "Km", "h"]
TAN_EMULSION = 0.5

# Histogram bins.  z in [0, 1]; E_h logarithmic from 0.3 GeV to the beam;
# p_T linear to 3 GeV; multiplicity 0..30.
Z_EDGES = np.linspace(0.0, 1.0, 21)
PT_EDGES = np.linspace(0.0, 3.0, 31)
N_EDGES = np.arange(0, 32) - 0.5

# The hadron-level z floor the "does this event have a pion?" question is
# asked above (user, 2026-09-07).  Same 0.1 as the yields are quoted for, and
# it is a BIN EDGE of Z_EDGES, which is what makes the count exact.
Z_TAG = 0.1

# generator key -> label, per current; the keys are the result keys
GENERATORS = {
    "nu": [("POWHEG-V2", "powheg_nu"), ("Herwig 7 (POWHEG)", "herwig_nlo_full"),
           ("Sherpa MC@NLO", "sherpa_nlo"), ("GENIE (GRV98LO)", "genie_lo"),
           ("GENIE (HEDIS)", "genie")],
    "mu": [("POWHEG-RES", "powheg"), ("Herwig 7 (POWHEG)", "herwig_nlo_powheg_full"),
           ("Sherpa MC@NLO", "sherpa"), ("GENIE (GRV98LO)", "genie")],
    # THE ANTINEUTRINO ROWS (user, 2026-10-04: "add the nubar contribution
    # to the SIDIS yields").  Same generators, own samples: POWHEG-V2 and
    # Sherpa from the FASER-ladder points (Q2 > 4 at generation, no W/y cut),
    # Herwig from its FASER ladder (tools/herwig_faser_ladder.sh), GENIE from
    # tools/genie_nubar_sidis.sh.  The region (q4w3) is applied here exactly
    # as for the neutrino; the cross-section is the generator's total times
    # the weighted fraction in the region (nubar_region_sigma).
    "nubar": [("POWHEG-V2", "powheg_nubar"), ("Herwig 7 (POWHEG)", "herwig_nlo_nubar"),
              ("Sherpa MC@NLO", "sherpa_nlo_nubar"), ("GENIE (GRV98LO)", "genie_lo_nubar"),
              ("GENIE (HEDIS)", "genie_nubar")],
}
# the neutrino row each antineutrino row is added to (faser_sidis)
NUBAR_OF = {"powheg_nu": "powheg_nubar", "herwig_nlo_full": "herwig_nlo_nubar",
            "sherpa_nlo": "sherpa_nlo_nubar", "genie_lo": "genie_lo_nubar",
            "genie": "genie_nubar"}


def e_edges(energy):
    return np.geomspace(0.3, float(energy), 31)


def resdir(current):
    return "results" if current == "mu" else "results_nu"


def spectra_path(current, key, energy, t):
    """results{,_nu}/pions_<key>_q4w3_<t>[_TAG].json"""
    return (f"{BASE}/{resdir(current)}/pions_"
            f"{beams.at_energy(sample_layout.result_key(key, t, REGION), energy)}.json")


def species_in(spectra):
    """The SPECIES every energy of these spectra carries.  "h" (all charged
    hadrons, 2026-09-30) was extracted for the neutrino current only -- the
    SIDIS study draws that current alone (CONVENTIONS.md 2b carve-out) -- so a
    muon row, whose spectra predate it, goes on without it rather than
    failing; re-extracting the muon samples would add it there too."""
    first = [spectra[e]["selections"] for e in spectra]
    return [h for h in SPECIES
            if all(f"{h}_emul_z" in sel for s_ in first for sel in s_.values())]


def events_sidecar_path(current, key, energy, t):
    """results{,_nu}/pions_events_<key>_q4w3_<t>[_TAG].npz -- the per-event
    record analysis/mhou_hadron_z.py joins to POWHEG's scale weights."""
    return (f"{BASE}/{resdir(current)}/pions_events_"
            f"{beams.at_energy(sample_layout.result_key(key, t, REGION), energy)}.npz")


# THE PER-EVENT SIDECAR (2026-09-30, for the MHOU band on the inclusive-hadron
# yields).  Only where scale weights exist to be joined: POWHEG-V2 at the
# 1 TeV anchor.  For each event in the yields region it keeps the file, the
# Les Houches number and the charged-hadron counts per z bin inside the
# emulsion acceptance; mhou_hadron_z.py does the join.
SIDECAR = {("nu", "powheg_nu", 1000.0)}
SIDECAR_REGION, SIDECAR_SPECIES = "sidis_e", "h_emul"


# ------------------------------------------------------------- the samples
def sample_files(current, key, energy, t):
    """The event files of one generator row on one nucleon (sample_layout).
    Herwig is its POSITIVE half: a per-event spectrum has no use for the
    negative-weight half, which is 0.1-3% of the rate and normalises the
    cross-section through the tracked histos, not through this pass."""
    if current == "nubar":
        return nubar_files(key, energy, t)
    if key in ("powheg", "powheg_nu"):
        files = job_files(f"{BASE}/powheg",
                          beams.at_energy(sample_layout.powheg_job_base(current, t), energy),
                          "events.hepmc")
        files = [f for f in files if os.path.exists(f"{os.path.dirname(f)}/V2_OK")]
    elif key.startswith("genie"):
        files = sample_layout.genie_files(current, key, t, energy, "events.hepmc")
    elif key.startswith("herwig"):
        files = job_files(f"{BASE}/herwig7",
                          beams.at_energy(sample_layout.herwig_job_base(current, t), energy),
                          "events.hepmc")
    elif key.startswith("sherpa"):
        files = job_files(sample_layout.sherpa_rundir(current, t, energy), "job", "evtfull")
    else:
        raise SystemExit(f"no sample rule for {current} {key}")
    return files


def nubar_files(key, energy, t):
    """The antineutrino samples of one row on one nucleon (2026-10-04)."""
    e = int(round(float(energy)))
    if key == "powheg_nubar":
        f = f"{paths.POWHEG_V2}/ladder-faser-v2/nubar_{t}_E{e}/sidis_events.hepmc"
        ok = os.path.exists(os.path.join(os.path.dirname(f), "SIDIS_SHOWER_OK"))
        return [f] if ok and os.path.exists(f) else []
    if key == "herwig_nlo_nubar":
        tag = beams.Beams("nu", float(energy)).tag
        return [f for f in (f"{BASE}/herwig7/flnubarpwg_{t}_job_{tag}_{j}/events.hepmc"
                            for j in (1, 2)) if os.path.exists(f)]
    if key == "sherpa_nlo_nubar":
        d = (f"{paths.SHERPA_RUNS}/NuDIS_NLO_FASER_nubar_p_E{e}" if t == "p"
             else f"{paths.SHERPA_RUNS}/V2_NuDIS_NLO_FASER_nubar_n_E{e}")
        ok = os.path.exists(f"{d}/SIDIS_GEN_OK") and os.path.exists(f"{d}/evtsidis")
        return [f"{d}/evtsidis"] if ok else []
    if key in ("genie_lo_nubar", "genie_nubar"):
        cfg = "nubar_grv" if key == "genie_lo_nubar" else "nubar_hedis"
        # job_files, NOT a glob of "<base>_*": at 1 TeV the base is untagged
        # and "_*" also catches job_400GeV_N (it did, 2026-10-04)
        base = beams.at_energy(f"v2g_{cfg}_{t}_job", float(energy))
        return [f for f in job_files(f"{BASE}/genie", base, "events.hepmc")
                if os.path.exists(os.path.join(os.path.dirname(f), "V2_OK"))]
    raise SystemExit(f"no antineutrino sample rule for {key}")


def nubar_total_sigma(key, energy, t):
    """(sigma_pb, err_pb): the generator's total over its generation region,
    on nucleon t, for the antineutrino rows.  Times the weighted fraction of
    the parsed events in the region, it is the region's cross-section."""
    e = int(round(float(energy)))
    if key == "powheg_nubar":
        import faser_dimuon_cutflow as fdc
        return fdc.powheg_total_pb(f"{paths.POWHEG_V2}/ladder-faser-v2/nubar_{t}_E{e}"), 0.0
    if key == "herwig_nlo_nubar":
        sys.path.insert(0, os.path.join(BASE, "herwig7"))
        import herwig_xsec as hx
        tag = beams.Beams("nu", float(energy)).tag
        pos = hx.read_jobs(f"{BASE}/herwig7/flnubarpwg_{t}_job_{tag}_*/FL-nubar-PWG-{t}-{tag}-S*.out", verbose=False)
        neg = hx.read_jobs(f"{BASE}/herwig7/flnubarpwgneg_{t}_job_{tag}_*/FL-nubar-PWGNEG-{t}-{tag}-S*.out", verbose=False)
        sp, ep, _n, _s = hx.combine(pos)
        sn, en, _n, _s = hx.combine(neg) if neg else (0.0, 0.0, 0, 1.0)
        return sp - sn, math.hypot(ep, en)
    if key == "sherpa_nlo_nubar":
        from faser_emulsion_shapes import sherpa_total_pb
        d = (f"{paths.SHERPA_RUNS}/NuDIS_NLO_FASER_nubar_p_E{e}" if t == "p"
             else f"{paths.SHERPA_RUNS}/V2_NuDIS_NLO_FASER_nubar_n_E{e}")
        return sherpa_total_pb(d), 0.0
    if key in ("genie_lo_nubar", "genie_nubar"):
        # from the job directories' events_xsec.json, which outlive the
        # events: NOT via nubar_files(), which is empty once the extraction
        # has deleted the hepmc (it was, 2026-10-04: sigma 0 at every point)
        cfg = "nubar_grv" if key == "genie_lo_nubar" else "nubar_hedis"
        base = beams.at_energy(f"v2g_{cfg}_{t}_job", float(energy))
        vals = []
        for f in job_files(f"{BASE}/genie", base, "events_xsec.json"):
            with open(f) as fh:
                vals.append(json.load(fh)["sigma_gen_mb"] * 1e9)
        if len(vals) != 4:
            raise SystemExit(f"{key} {energy} {t}: {len(vals)} GENIE cross-sections, want 4")
        return float(np.mean(vals)), 0.0
    raise SystemExit(f"no antineutrino cross-section rule for {key}")


# --------------------------------------------------------------- extraction
def extract(current, key, energy, t):
    if t not in ("p", "n"):
        raise SystemExit("a sample is ONE nucleon: --target p or n")
    files = sample_files(current, key, energy, t)
    if not files:
        raise SystemExit(f"no files for {current} {key} {t} at {energy:g} GeV")
    lep_out = 13
    if current == "nubar":
        lep_out = -13
        analyze_nu.NU_BEAM_PID = -14
        events = lambda fn: analyze_nu.parse_hepmc3(fn)          # noqa: E731
    elif current == "nu":
        analyze_nu.NU_BEAM_PID = 14
        events = lambda fn: analyze_nu.parse_hepmc3(fn)          # noqa: E731
    else:
        events = lambda fn: analyze.parse_hepmc3(fn, beam_pid=13)  # noqa: E731
    ee = e_edges(energy)
    acc = {}
    for sname in SELS:
        acc[sname] = {"sum_w": 0.0, "n": 0, "sum_w2": 0.0}
        for h in SPECIES:
            for reg in ("all", "emul"):
                acc[sname][f"{h}_{reg}_z"] = np.zeros(len(Z_EDGES) - 1)
                # SECOND MOMENT of the per-event count in each z bin, for the
                # statistical error a real measurement would have: pions in a
                # bin are NOT independent counts, several come from the same
                # event, so the variance of the bin total is N_ev <n^2> and
                # not N_ev <n>.  See faser_sidis.stat_error().
                acc[sname][f"{h}_{reg}_z2"] = np.zeros(len(Z_EDGES) - 1)
                acc[sname][f"{h}_{reg}_E"] = np.zeros(len(ee) - 1)
                acc[sname][f"{h}_{reg}_pt"] = np.zeros(len(PT_EDGES) - 1)
                acc[sname][f"{h}_{reg}_n"] = np.zeros(len(N_EDGES) - 1)
                # THE SAME COUNT, but only over hadrons above Z_TAG.  It
                # answers "what fraction of DIS events has a pion in it?"
                # (user, 2026-09-07) for the z range a measurement would
                # actually use: 1 - hist[0] is that fraction.
                acc[sname][f"{h}_{reg}_nz"] = np.zeros(len(N_EDGES) - 1)
                acc[sname][f"{h}_{reg}_sum"] = 0.0
    n_parsed = 0
    # THE GENERATION-CUT TRIPWIRE, kept from the earlier study: a sample generated
    # inside a y window holds nothing below y = 0.2, and any no-y-window
    # number from it is silently TRUNCATED.  Every sample is generated
    # without a window, so this now guards against a wrong sample table.
    n_q2 = n_lowy = 0
    w_lowy = 0.0
    sum_w_all = 0.0          # every parsed event: the antineutrino rows
                             # normalise by it (nubar_total_sigma)
    side = ((current, key, float(energy)) in SIDECAR)
    side_rows = []           # (file index, lhe_index, weight, counts per z bin)
    for ifn, fn in enumerate(files):
        print(f"  parsing {fn}", flush=True)
        for w, k, P, parts, _dmes, _hard in events(fn):
            n_parsed += 1
            w = 1.0 if w is None else w
            sum_w_all += w
            best, best_e = None, -1.0
            for pid, p in parts:
                if pid == lep_out:
                    e = lab_energy(p, P)
                    if e > best_e:
                        best, best_e = p, e
            if best is None:
                continue
            q2, y, xbj, kP = dis_invariants(k, P, best)
            # the pre-filter is the Q2 floor alone; every region applies its
            # own W and x cuts below
            if not selmod.ALL_Y.passes(q2, y):
                continue
            n_q2 += 1
            if y < 0.2:
                n_lowy += 1
                w_lowy += w
            theta = lepton_theta(best_e, best)
            nu_had = kP / M_P - best_e
            # W^2 = m_p^2 + 2 P.q - Q^2, and P.q = m_p nu because the target
            # is at rest in the frame lab_energy() is measured in.
            w2 = M_P * M_P + 2.0 * M_P * nu_had - q2
            phi_mu = math.atan2(best[2], best[1])
            lep_ok = {sname: (sel.passes(q2, y, w2=w2, x=xbj)
                              and sel.passes_lepton(best_e, theta))
                      for sname, sel in SELS.items()}
            n05 = n01 = 0
            sum_px = sum_py = 0.0
            had = []          # (species, signed species, E_lab, pT, inside acceptance)
            for pid, p in parts:
                a = abs(pid)
                if a not in CHARGED_HADRONS:
                    continue
                elab_h = dot(p, P) / M_P
                # signed p_z,lab, and Tier E counts only TRACKS above 1 GeV
                # (2026-09-19, review S1; analyze.track_tan)
                pt_h, tan_h = analyze.track_tan(p, elab_h, P, k)
                if tan_h < 0.5 and elab_h > TRACK_E_MIN:
                    n05 += 1
                    if tan_h < 0.1:
                        n01 += 1
                sum_px += p[1]
                sum_py += p[2]
                if a == HADRONS["pi"]:
                    had.append(("pi", "pip" if pid > 0 else "pim",
                                elab_h, pt_h, tan_h < TAN_EMULSION))
                elif a == HADRONS["K"]:
                    had.append(("K", "Kp" if pid > 0 else "Km",
                                elab_h, pt_h, tan_h < TAN_EMULSION))
                # and every charged hadron, pions and kaons included, once
                had.append(("h", None, elab_h, pt_h, tan_h < TAN_EMULSION))
            dphi_x = (abs(math.remainder(math.atan2(sum_py, sum_px) - phi_mu, 2.0 * math.pi))
                      if (sum_px or sum_py) else None)
            passed = {sname: (ok and SELS[sname].passes_hadrons(n05, n01, dphi_x))
                      for sname, ok in lep_ok.items()}
            if not any(passed.values()):
                continue
            # the hadrons do not depend on the region: binned ONCE per event
            # into sparse per-event counts, added to every passing region
            counts = {}          # pre -> total pions
            countz = {}          # pre -> pions above Z_TAG
            zc = {}              # pre -> {z bin: n in this event}
            ec = {}              # pre -> {E bin: n}
            pc = {}              # pre -> {pT bin: n}
            for sp, sp_q, e_h, pt_h, in_em in had:
                z = e_h / nu_had if nu_had > 0 else 0.0
                iz = np.searchsorted(Z_EDGES, z, side="right") - 1
                ie = np.searchsorted(ee, e_h, side="right") - 1
                ip = np.searchsorted(PT_EDGES, pt_h, side="right") - 1
                for reg in (("all", "emul") if in_em else ("all",)):
                    for pre in ((f"{sp}_{reg}", f"{sp_q}_{reg}") if sp_q
                                else (f"{sp}_{reg}",)):
                        counts[pre] = counts.get(pre, 0) + 1
                        if z >= Z_TAG:
                            countz[pre] = countz.get(pre, 0) + 1
                        if 0 <= iz < len(Z_EDGES) - 1:
                            d = zc.setdefault(pre, {})
                            d[iz] = d.get(iz, 0) + 1
                        if 0 <= ie < len(ee) - 1:
                            d = ec.setdefault(pre, {})
                            d[ie] = d.get(ie, 0) + 1
                        if 0 <= ip < len(PT_EDGES) - 1:
                            d = pc.setdefault(pre, {})
                            d[ip] = d.get(ip, 0) + 1
            if side and passed.get(SIDECAR_REGION):
                li = _hard.get("lhe_index") if isinstance(_hard, dict) else None
                row = np.zeros(len(Z_EDGES) - 1)
                for i, n in zc.get(SIDECAR_SPECIES, {}).items():
                    row[i] = n
                side_rows.append((ifn, -1 if li is None else int(li), w, row))
            for sname, ok in passed.items():
                if not ok:
                    continue
                a_ = acc[sname]
                a_["sum_w"] += w
                a_["sum_w2"] += w * w
                a_["n"] += 1
                for pre, c in counts.items():
                    a_[pre + "_sum"] += w * c
                for pre, d in zc.items():
                    hz, hz2 = a_[pre + "_z"], a_[pre + "_z2"]
                    for i, n in d.items():
                        hz[i] += w * n
                        hz2[i] += w * n * n
                for pre, d in ec.items():
                    h = a_[pre + "_E"]
                    for i, n in d.items():
                        h[i] += w * n
                for pre, d in pc.items():
                    h = a_[pre + "_pt"]
                    for i, n in d.items():
                        h[i] += w * n
                for h_ in SPECIES:
                    for reg in ("all", "emul"):
                        pre = f"{h_}_{reg}"
                        a_[pre + "_n"][min(counts.get(pre, 0),
                                           len(N_EDGES) - 2)] += w
                        a_[pre + "_nz"][min(countz.get(pre, 0),
                                            len(N_EDGES) - 2)] += w
    out = {"current": current, "key": key, "energy_gev": float(energy),
           "target": t, "region": REGION, "region_aliases": REGION_ALIAS,
           "n_parsed": n_parsed, "sum_w_all": sum_w_all,
           "y_coverage": {"n_above_q2_floor": n_q2, "n_below_y02": n_lowy,
                          "sum_w_below_y02": w_lowy,
                          "frac_below_y02": (n_lowy / n_q2) if n_q2 else 0.0,
                          "note": "fraction of the sample below y = 0.2, "
                                  "above the Q2 floor. Zero means the sample "
                                  "was GENERATED inside a y window and "
                                  "cannot serve a region without one."},
           "lumi_note": "per-event spectra, no luminosity",
           "z_edges": Z_EDGES.tolist(), "e_edges": ee.tolist(),
           "pt_edges": PT_EDGES.tolist(), "n_edges": N_EDGES.tolist(),
           "tan_emulsion": TAN_EMULSION, "z_tag": Z_TAG,
           "inputs": analyze.input_manifest(files), "selections": {}}
    for sname, a_ in acc.items():
        sw = a_["sum_w"]
        rec = {"n_selected": a_["n"], "sum_w": sw,
               "n_eff": (sw * sw / a_["sum_w2"]) if a_["sum_w2"] else 0.0}
        for kname, v in a_.items():
            if kname in ("sum_w", "n", "sum_w2"):
                continue
            if isinstance(v, np.ndarray):
                rec[kname] = (v / sw).tolist() if sw else v.tolist()
            else:
                rec[kname] = v / sw if sw else 0.0      # mean multiplicity
        out["selections"][sname] = rec
    p = spectra_path(current, key, energy, t)
    with open(p, "w") as f:
        json.dump(out, f)
    if side:
        sp_ = events_sidecar_path(current, key, energy, t)
        np.savez_compressed(
            sp_, files=np.array([os.path.relpath(f, BASE) for f in files]),
            file_idx=np.array([r[0] for r in side_rows], dtype=np.int32),
            lhe_index=np.array([r[1] for r in side_rows], dtype=np.int64),
            w=np.array([r[2] for r in side_rows], dtype=float),
            **{SIDECAR_SPECIES: np.array([r[3] for r in side_rows],
                                         dtype=np.int32).reshape(-1, len(Z_EDGES) - 1)})
        print(f"  per-event sidecar: {len(side_rows)} events -> {sp_}")
    print(f"  y coverage: {n_lowy}/{n_q2} events below y = 0.2 "
          f"({100.0 * n_lowy / n_q2 if n_q2 else 0.0:.2f}%)")
    for sname in SELS:
        r = out["selections"][sname]
        print(f"  {sname:12s} {r['n_selected']:8d} events, pi/event "
              f"{r['pi_all_sum']:.3f} (emulsion acc. {r['pi_emul_sum']:.3f}), "
              f"pi+ {r['pip_emul_sum']:.3f} pi- {r['pim_emul_sum']:.3f}")
    print(f"{n_parsed} events; written {p}")


# ----------------------------------------------------------- the combination
def region_sigma(current, key, energy, t):
    """(sigma, err) in pb per nucleon of the region on nucleon t: the
    tracked histos_<key>_q4w3_<t>[_TAG].json.  For the antineutrino rows:
    the generator's total times this pass's weighted fraction in the region."""
    if current == "nubar":
        with open(spectra_path(current, key, energy, t)) as f:
            sp = json.load(f)
        tot, err = nubar_total_sigma(key, energy, t)
        if not tot > 0:
            raise SystemExit(f"{key} {energy} {t}: generator total {tot} pb -- refusing")
        frac = sp["selections"][ROOT]["sum_w"] / sp["sum_w_all"] if sp["sum_w_all"] else 0.0
        return tot * frac, err * frac
    p = (f"{BASE}/{resdir(current)}/histos_"
         f"{beams.at_energy(sample_layout.result_key(key, t, REGION), energy)}.json")
    if not os.path.exists(p):
        raise SystemExit(f"no cross-section for {current} {key} {t} at "
                         f"{energy:g} GeV ({os.path.relpath(p, BASE)}): run "
                         f"tools/analyse_production.py --energies {energy:g}")
    with open(p) as f:
        d = json.load(f)
    if d.get("target") != t:
        raise SystemExit(f"{p} is stamped target {d.get('target')!r}, not {t!r}")
    return float(d["sigma_fid_pb"]), float(d.get("sigma_fid_err_pb", 0.0))


def combine(current, key, energy):
    """Tungsten per nucleon from the proton and neutron passes."""
    parts = {}
    for t in ("p", "n"):
        p = spectra_path(current, key, energy, t)
        if not os.path.exists(p):
            raise SystemExit(f"missing {os.path.relpath(p, BASE)} -- run extract")
        with open(p) as f:
            parts[t] = json.load(f)
        if parts[t]["target"] != t:
            raise SystemExit(f"{p} is stamped target {parts[t]['target']!r}")
    for k in ("energy_gev", "z_edges", "e_edges", "pt_edges", "n_edges"):
        if parts["p"][k] != parts["n"][k]:
            raise SystemExit(f"p and n spectra differ in {k}")
    sig = {t: region_sigma(current, key, energy, t) for t in ("p", "n")}
    zn = {"p": float(target.Z_W), "n": float(target.N_W)}
    out = {k: v for k, v in parts["p"].items() if k != "selections"}
    out.update({"target": "W", "n_parsed": {t: parts[t]["n_parsed"] for t in parts},
                "inputs": {t: parts[t]["inputs"] for t in parts},
                "y_coverage": {t: parts[t]["y_coverage"] for t in parts},
                "combination": f"(74 sigma_p + 110 sigma_n)/184 of each region; "
                               f"the region's sigma_t is the tracked q4w3 "
                               f"cross-section of nucleon t times this pass's "
                               f"fraction of the {ROOT} events in it",
                "selections": {}})
    root = {t: parts[t]["selections"][ROOT]["sum_w"] for t in parts}
    for sname in SELS:
        rec = {"n_selected": sum(parts[t]["selections"][sname]["n_selected"] for t in parts),
               "sum_w": {t: parts[t]["selections"][sname]["sum_w"] for t in parts}}
        wgt, sig_r, err_r = {}, {}, {}
        for t in parts:
            s = parts[t]["selections"][sname]
            frac = s["sum_w"] / root[t] if root[t] else 0.0
            # the fraction's own (binomial) error on n_eff events
            neff = parts[t]["selections"][ROOT]["n_eff"]
            dfrac = math.sqrt(max(frac * (1.0 - frac), 0.0) / neff) if neff else 0.0
            sig_r[t] = sig[t][0] * frac
            err_r[t] = math.hypot(sig[t][1] * frac, sig[t][0] * dfrac)
            wgt[t] = zn[t] * sig_r[t]
            rec[f"frac_{t}"] = frac
            rec[f"sigma_{t}_pb"] = sig_r[t]
        tot = wgt["p"] + wgt["n"]
        rec["sigma_pb"] = tot / target.A_W
        rec["sigma_err_pb"] = math.hypot(zn["p"] * err_r["p"], zn["n"] * err_r["n"]) / target.A_W
        rec["weight_p"] = wgt["p"] / tot if tot else 0.5
        sp_, sn_ = parts["p"]["selections"][sname], parts["n"]["selections"][sname]
        for kname, v in sp_.items():
            if kname in ("n_selected", "sum_w", "n_eff"):
                continue
            if isinstance(v, list):
                rec[kname] = ((rec["weight_p"] * np.array(v)
                               + (1.0 - rec["weight_p"]) * np.array(sn_[kname])).tolist())
            else:
                rec[kname] = rec["weight_p"] * v + (1.0 - rec["weight_p"]) * sn_[kname]
        out["selections"][sname] = rec
    p = spectra_path(current, key, energy, "W")
    with open(p, "w") as f:
        json.dump(out, f)
    for sname in SELS:
        r = out["selections"][sname]
        print(f"  {sname:12s} sigma_W {r['sigma_pb']:.5g} +- {r['sigma_err_pb']:.2g} pb "
              f"(p {r['sigma_p_pb']:.5g}, n {r['sigma_n_pb']:.5g}; frac p "
              f"{r['frac_p']:.4f} n {r['frac_n']:.4f})  pi/event {r['pi_emul_sum']:.3f}")
    print(f"written {p}")


# ------------------------------------------------------------------- rates
def load_spectra(current, key, t="W"):
    out = {}
    for e in ENERGIES:
        p = spectra_path(current, key, e, t)
        if os.path.exists(p):
            with open(p) as f:
                out[e] = json.load(f)
    return out


def yadism_shape(current):
    """(E, sigma) of the YADISM NLO (ZM-VFNS, TMC on) tungsten ladder of the
    region, the shape every generator's rate rides across the flux.

    The antineutrino rows ride the NEUTRINO ladder (resdir("nubar") is
    results_nu): it is only the interpolation template, and each row's own
    ratio to it at its six energies (300 GeV-4 TeV) carries the antineutrino
    energy dependence; outside them the ratio is held at the end values."""
    es, ss = [], []
    for e in ENERGIES:
        p = (f"{BASE}/{resdir(current)}/histos_"
             f"{beams.at_energy('yadism_nlo_' + REGION + '_W_tmc', e)}.json")
        if os.path.exists(p):
            with open(p) as f:
                es.append(e)
                ss.append(float(json.load(f)["sigma_fid_pb"]))
    if len(es) < 2:
        raise SystemExit(f"no YADISM {REGION} tungsten ladder for {current}")
    return np.array(es), np.array(ss)


def _loglog(e, xs, ys):
    """log-log interpolation, extrapolated with the end slopes."""
    le, lx, ly = np.log(np.atleast_1d(e)), np.log(xs), np.log(ys)
    out = np.interp(le, lx, ly)
    lo, hi = le < lx[0], le > lx[-1]
    out[lo] = ly[0] + (le[lo] - lx[0]) * (ly[1] - ly[0]) / (lx[1] - lx[0])
    out[hi] = ly[-1] + (le[hi] - lx[-1]) * (ly[-1] - ly[-2]) / (lx[-1] - lx[-2])
    return np.exp(out)


def region_sigma_curve(current, key, region, spectra):
    """sigma^W_region(E) as a callable over the flux, and how it was built."""
    es = np.array(sorted(e for e in spectra if region in spectra[e]["selections"]))
    if len(es) == 0:
        return None, None
    ye, ys = yadism_shape(current)
    root = np.array([spectra[e]["selections"][ROOT]["sigma_pb"] for e in es])
    r_root = root / _loglog(es, ye, ys)
    le = np.log(es)
    if region == ROOT:
        def sig(e):
            e = np.atleast_1d(np.asarray(e, dtype=float))
            r = np.exp(np.interp(np.clip(np.log(e), le[0], le[-1]), le, np.log(r_root)))
            return _loglog(e, ye, ys) * r
        return sig, ("YADISM NLO (ZM-VFNS, TMC) shape of the benchmark region x this "
                     f"generator's ratio at {', '.join(f'{x:g}' for x in es)} GeV")
    frac = np.array([spectra[e]["selections"][region]["sigma_pb"] for e in es]) / root
    root_sig, _how = region_sigma_curve(current, key, ROOT, spectra)
    e_zero = beams.TIER_E_ZERO_ENERGY

    def sig(e):
        e = np.atleast_1d(np.asarray(e, dtype=float))
        eff = np.interp(np.clip(np.log(e), le[0], le[-1]), le, frac)
        if region in TIER_E_REGIONS:
            # exactly zero at E <= 200 GeV, linear in log E up to the first
            # measured point when that point is above the zero (every
            # generator but POWHEG-V2, whose ladder reaches 20 GeV and
            # MEASURES the zero at 200 GeV)
            if es[0] > e_zero:
                below = e < es[0]
                eff[below] = frac[0] * np.clip(np.log(e[below] / e_zero)
                                               / np.log(es[0] / e_zero), 0.0, None)
            eff[e <= e_zero] = 0.0
        return root_sig(e) * eff
    return sig, (f"the {ROOT} rate x this pass's {region} fraction at "
                 f"{', '.join(f'{x:g}' for x in es)} GeV"
                 + (f", zero at {e_zero:g} GeV" if region in TIER_E_REGIONS else ""))


def flux_weights(current, region, key, spectra, pid="14"):
    """(E points, events per E point at LUMI_FB, how) for one generator.

    pid picks the neutrino FLAVOUR of the flux (14 = nu_mu, 12 = nu_e).  The
    cross-section is flavour blind -- massless-lepton CC DIS is the same
    calculation for an electron and for a muon in the final state -- so a
    nu_e rate is the same sigma against a different flux.  It is ignored on
    the muon current, which has one flux."""
    e_pts, phi, norm = fr.flux_weights("nu" if current == "nubar" else current, pid)
    sig, how = region_sigma_curve(current, key, region, spectra)
    if sig is None:
        return None, None, None
    n = phi * sig(e_pts) * norm * (LUMI_FB / fr.LUMI_FB)
    return e_pts, n, how


def shape_energies(spectra, sel):
    """The energies at which a selection has a SHAPE: those with events.

    AN EMPTY SELECTION HAS NO SHAPE, NOT A ZERO ONE.  POWHEG-V2's ladder
    reaches below the Tier E zero (200 GeV, beams.SIDIS_ENERGIES_LOW),
    where every per-event histogram and mean is zero because no event
    passes; interpolating through those points would pull the MULTIPLICITY
    towards zero between 200 and 300 GeV while the RATE already ramps to
    zero there, suppressing that band twice (measured: -0.14% on the
    POWHEG-V2 Tier E pions before this was fixed).  So shapes and means are
    read only where the selection has events, and held constant outside.
    """
    return [e for e in sorted(spectra)
            if spectra[e]["selections"][sel].get("sigma_pb", 1.0) > 0]


def interp_spectrum(spectra, sel, name, e_pts):
    """The per-event histogram at every flux energy: linear in log E between
    the computed energies, held constant outside (the E_h axis, whose bins
    scale with the beam, is interpolated in the bin INDEX, i.e. in E_h/E)."""
    es = shape_energies(spectra, sel)
    arr = np.array([spectra[e]["selections"][sel][name] for e in es])   # (nE, nbins)
    le = np.log(np.array(es))
    x = np.clip(np.log(e_pts), le[0], le[-1])
    out = np.empty((len(e_pts), arr.shape[1]))
    for j in range(arr.shape[1]):
        out[:, j] = np.interp(x, le, arr[:, j])
    return out


def rates():
    out = {"lumi_fb": LUMI_FB, "tan_emulsion": TAN_EMULSION, "region": REGION,
           "hadron_selection_eff": HADRON_SELECTION_EFF,
           "region_aliases": REGION_ALIAS, "energies": list(ENERGIES),
           "target": "tungsten per nucleon, 74 p + 110 n",
           "z_edges": Z_EDGES.tolist(), "pt_edges": PT_EDGES.tolist(),
           "n_edges": N_EDGES.tolist(), "currents": {}}
    for current in ("nu", "mu"):
        cur = {"generators": {}}
        for label, key in GENERATORS[current]:
            sp = load_spectra(current, key)
            if not sp:
                print(f"{current} {key}: no spectra yet")
                continue
            rec = {"label": label, "energies": sorted(sp), "selections": {}}
            for sel in SELS:
                e_pts, n, how = flux_weights(current, sel, key, sp)
                if n is None:
                    print(f"{current} {key} {sel}: no rate")
                    continue
                srec = {"events": float(n.sum()), "normalisation": how,
                        "sigma_pb": {f"{e:g}": sp[e]["selections"][sel]["sigma_pb"]
                                     for e in sorted(sp)}}
                for h in species_in(sp):
                    for reg in ("all", "emul"):
                        pre = f"{h}_{reg}"
                        for obs in ("z", "pt", "n"):
                            spec = interp_spectrum(sp, sel, f"{pre}_{obs}", e_pts)
                            srec[f"{pre}_{obs}"] = hadron_yield(
                                (n[:, None] * spec).sum(axis=0)).tolist()
                        mult = np.array([spectra_mean(sp, sel, pre + "_sum", e)
                                         for e in e_pts])
                        gen_total = float((n * mult).sum())
                        # THE MULTIPLICITY IS THE GENERATED ONE, the yield is
                        # not: per_event is a property of the event, the total
                        # is what the detector would select.
                        srec[f"{pre}_per_event"] = (gen_total / srec["events"]
                                                    if srec["events"] else 0.0)
                        srec[f"{pre}_total"] = float(hadron_yield(gen_total))
                    if h in HADRONS:
                        e_spec, srec["E_edges_common"] = absolute_e_spectrum(
                            sp, sel, f"{h}_emul_E", e_pts, n)
                        srec[f"{h}_emul_E"] = hadron_yield(e_spec).tolist()
                rec["selections"][sel] = srec
                print(f"{current:2s} {key:24s} {sel:12s} events {srec['events']:8.1f}  "
                      f"pi/event {srec['pi_all_per_event']:.3f} (emulsion acc. "
                      f"{srec['pi_emul_per_event']:.3f})  pions {srec['pi_emul_total']:9.0f}  "
                      f"kaons {srec['K_emul_total']:8.0f}")
            cur["generators"][key] = rec
        out["currents"][current] = cur
    p = f"{BASE}/results_nu/faser_pions.json"
    with open(p, "w") as f:
        json.dump(out, f)
    print(f"written {p}")
    return out


def spectra_mean(spectra, sel, name, e):
    es = shape_energies(spectra, sel)
    vals = [spectra[x]["selections"][sel][name] for x in es]
    le = np.log(np.array(es))
    return float(np.interp(np.clip(np.log(e), le[0], le[-1]), le, vals))


COMMON_E_EDGES = np.geomspace(0.3, 6000.0, 44)


def absolute_e_spectrum(spectra, sel, name, e_pts, n):
    """dN/dE_h on a fixed log grid: each computed energy's histogram is
    rebinned onto it, interpolated in log E between energies, and weighted."""
    es = shape_energies(spectra, sel)
    reb = []
    for e in es:
        d = spectra[e]
        src_edges = np.array(d["e_edges"])
        vals = np.array(d["selections"][sel][name])
        tgt = np.zeros(len(COMMON_E_EDGES) - 1)
        for i in range(len(vals)):
            lo, hi = np.log(src_edges[i]), np.log(src_edges[i + 1])
            for j in range(len(tgt)):
                a, b = np.log(COMMON_E_EDGES[j]), np.log(COMMON_E_EDGES[j + 1])
                ov = max(0.0, min(hi, b) - max(lo, a))
                if ov > 0:
                    tgt[j] += vals[i] * ov / (hi - lo)
        reb.append(tgt)
    reb = np.array(reb)
    le = np.log(np.array(es))
    x = np.clip(np.log(e_pts), le[0], le[-1])
    out = np.zeros(reb.shape[1])
    for j in range(reb.shape[1]):
        out[j] = float((n * np.interp(x, le, reb[:, j])).sum())
    return out.tolist(), COMMON_E_EDGES.tolist()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    ex = sub.add_parser("extract")
    ex.add_argument("--current", required=True)
    ex.add_argument("--key", required=True)
    ex.add_argument("--energy", type=float, required=True)
    ex.add_argument("--target", required=True, choices=("p", "n"))
    co = sub.add_parser("combine")
    co.add_argument("--current", required=True)
    co.add_argument("--key", required=True)
    co.add_argument("--energy", type=float, required=True)
    sub.add_parser("rates")
    a = ap.parse_args()
    if a.cmd == "extract":
        extract(a.current, a.key, a.energy, a.target)
    elif a.cmd == "combine":
        combine(a.current, a.key, a.energy)
    else:
        rates()


if __name__ == "__main__":
    main()
