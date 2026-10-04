#!/usr/bin/env python3
"""Analysis for the CC NEUTRINO-DIS generator benchmark (1 TeV nu_mu on p).

Charged-current counterpart of analysis/analyze.py.  Deliberately a SEPARATE
file writing into a SEPARATE directory (results_nu/) so the muon results can
never be clobbered.  It imports BINS and the kinematic helpers from analyze.py
so the bin edges and the fiducial definition are guaranteed identical to the
muon study and the two can be overlaid.

Differences from the muon analysis:
  * the beam lepton is a nu_mu (pid 14), the scattered lepton is a mu- (pid 13)
  * Sherpa neutrino-beam samples carry a factor-2 helicity-average error (see
    SHERPA_NU_SPIN_FACTOR below and Runs/NuDIS_LO/Sherpa.yaml)

Fiducial region (identical): Q2 > 4 GeV2, 0.2 < y < 0.9, all invariants built
from the beam records of each event, so the generation frame is irrelevant.

Usage:
  analyze_nu.py pythia    # pythia8/nu_job_*/events.hepmc + *_xsec.json
  analyze_nu.py sherpa_lo # Sherpa Runs/NuDIS_LO/job_*/evtfull
  analyze_nu.py genie     # genie/nucc_job_*/events.hepmc + *_xsec.json
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402
import beams  # noqa: E402
import selection as selection_mod  # noqa: E402  -- config.sh is the one source of every path
from analyze import (track_tan, TRACK_E_MIN, NEUTRINO_PIDS, subtracts_shower_charm, stamp_shower_charm, QED_ARMS, sherpa_arm_sigma_pb, BINS, ENERGY, SELECTION, at_energy, job_files, sel_suffix, CHARGED_HADRONS, D_0_PID, D_CH_PID, D_HIST_KEYS, n_opposite_sign_muons, fmt_sigma, flag_low_stats,
                     closure_check, delivered_sigma_pb, input_manifest,
                     herwig_lhe_delivered_pb, lhe_events_offered,
                     NO_DPHI, dphi_summary,
                     CC_BOSONS, find_hard_vertex, is_heavy_hadron,
                     pooled_sigma_pb,
                     D_PIDS_SIGNED, M_MU, M_P, Q2_MIN, W2_MIN, Y_MAX, Y_MIN,
                     dedup_dmesons, weakly_decaying_charm, dmeson_summary, dot,
                     dis_invariants, lab_energy,
                     lepton_theta, check_anchor_literals)  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = f"{BASE}/results_nu"
SHERPA_RUNS = paths.SHERPA_RUNS

# --- Sherpa neutrino-beam normalisation ------------------------------------
# TWO independent corrections are needed for the Sherpa CC neutrino samples.
#
# (1) FACTOR 2.  Sherpa 3.0.5 averages over two initial lepton helicities even
#     for a neutrino beam, where only the left-handed one exists, so its ME is
#     exactly half the physical value.  Verified channel by channel against the
#     analytic parton model (0.4997 / 0.5000 / 0.5001 in the flat, the (1-y)^2
#     and the summed channel) and by the fact that the same setup with an e-
#     beam is correct.  Unlike Herwig this CANNOT be fixed in the card:
#     BEAM_POLARIZATIONS is parsed (BEAM/Main/Beam_Parameters.C) but the
#     function that would use it, Beam_Spectra_Handler::Weight(pol_types,dofs),
#     is COMMENTED OUT in BEAM/Main/Weight_Base.C ("The following only if we
#     ever go back to polarised beams"), and Comix hard-codes 2 polarisations
#     for every spin-1/2 flavour (COMIX/Amplitude/Amplitude.C:113,538).  So the
#     factor stays a named constant here.  It is a pure normalisation.
SHERPA_NU_SPIN_FACTOR = 2.0
# (2) The samples are normalised to Sherpa's EVENT-LEVEL cross-section
#     ("Nominal" in the end-of-run table of job_*/sherpa.log), NOT to the
#     integrator value in integ.log.  With RESPECT_MASSIVE_FLAG: true the two
#     agree; before that fix 23% of the events were vetoed at the ME->shower
#     interface and the two differed by 22% (see NEUTRINO_NOTES.md).  Using the
#     event-level number means a repeat of that pathology shows up immediately
#     instead of being hidden by the integrator.
SHERPA_NU_LO_SIGMA_RAW_PB = 1.59530        # event level, mean of 8 jobs
SHERPA_NU_LO_SIGMA_RAW_ERR_PB = 0.00284
SHERPA_NU_LO_INTEG_RAW_PB = 1.59692        # integrator, +- 0.00118
SHERPA_NU_LO_SIGMA_PB = SHERPA_NU_SPIN_FACTOR * SHERPA_NU_LO_SIGMA_RAW_PB
SHERPA_NU_LO_SIGMA_ERR_PB = SHERPA_NU_SPIN_FACTOR * SHERPA_NU_LO_SIGMA_RAW_ERR_PB
# same, ME-level run dir (Runs/NuDIS_LO_ME): unaffected by the massive-PS bug
# (SHOWER_GENERATOR: None -> the interface returns before the veto)
SHERPA_NU_LO_ME_SIGMA_PB = SHERPA_NU_SPIN_FACTOR * 1.59691
SHERPA_NU_LO_ME_SIGMA_ERR_PB = SHERPA_NU_SPIN_FACTOR * 0.00118

# --- Herwig neutrino-beam normalisation ------------------------------------
# Herwig 7.3.0 MEDISCC has the SAME factor-2 neutrino helicity-average problem
# as Sherpa (MEChargedCurrentDIS.cc:212, unconditional `me *= 0.25`), BUT
# unlike Sherpa it can be fixed IN THE CARD: declaring the beam as a
# ThePEG::PolarizedBeamParticleData with LongitudinalPolarization = -1 makes
# the ME take the `menew.average(rho[0],rho[1])` branch instead.  That is what
# herwig7/DIS-nu-POL.in does, and it reproduces exactly 2 x the unpolarised
# result (ME level: 3.18838 vs 2 x 1.59419 pb, identical to 6 digits).
# The samples in herwig7/nu_job_* are generated with DIS-nu-POL.in, so NO
# hand-coded factor is applied any more.  Set this back to 2.0 only for a
# sample made with the old unpolarised DIS-nu.in (archived in
# herwig7/old_nu_unpol/).
HERWIG_NU_SPIN_FACTOR = 1.0

# --- Sherpa CC MC@NLO (Runs/NuDIS_NLO_ckm3) --------------------------------
# Same event-level (not integrator) convention as the LO sample above.
# Integrator for reference: BVI 1.62089 + RS -0.09807 = 1.52282 +- 0.00107 pb;
# the event-level mean is 0.9975 of it, the usual MC@NLO S/H-event difference.
# The spin factor 2 applies here too, and that was verified NON-CIRCULARLY by
# a K-factor test rather than assumed: Sherpa's internal
# K = 1.52370/1.59510(LO, same full CKM) = 0.95524 against YADISM's
# 3.0464/3.1896 = 0.95510, agreeing to 0.1%.  A ratio internal to Sherpa is
# independent of the overall helicity factor, so this shows the factor is
# common to B, V, I and RS -- a B-only run would have proved nothing.
# FULL CKM since 2026-08-25 (Runs/NuDIS_NLO_ckm3).  The predecessor was
# forced to CKM Order 0 because Sherpa's virtual-ME getter refused the Cabibbo
# cross-family channels; patches/sherpa-dy-qcd-virtual-ckm.diff lifts that, so
# this sample is now like-for-like with the LO card and CAN be used for the CC
# charm strand.  Integrator moved 1.52282 -> 1.52370 pb (+0.058%, 0.50 sigma),
# i.e. unchanged, as CKM unitarity requires; the point of the change is the
# charm COMPOSITION, which unit CKM got wrong by ~3 points.
# >>> ONE TABLE, KEYED BY ENERGY -- NOT THREE LOOSE CONSTANTS. <<<
# The sample path and the normalisation MUST become energy-aware together.
# Making one of them energy-aware while the other stays a 1 TeV literal is the
# single failure this benchmark has repeated most often: it normalises one
# energy's events with another energy's cross-section, nothing errors, and the
# result looks entirely plausible.  Keeping the numbers in a dict indexed by
# the same tag the run directory uses means a missing energy raises instead.
#
# raw = before the factor-2 helicity correction, which is what the event
# weights carry.  Event level is the mean over the 8 jobs; the integrator is
# BVI + RS from that energy's own integ.log.
SHERPA_NU_NLO = {
    #  tag       event-level        err        integrator
    "400GeV": (0.61452,          0.00187,     0.61503),
    "1TeV":   (1.51994,          0.00368,     1.52370),
    "4TeV":   (5.46126,          0.00975,     5.48364),
}
# Cross-check at the anchor: the pooled sum(w)/sum(ntrials) over all 200k
# events is 1.51990 pb, agreeing with the job-wise mean to 0.002%.  The two
# estimators weight the jobs differently, so their agreeing says the jobs are
# homogeneous -- which is what the per-job spread is quoted from.


def _sherpa_nu_nlo(which):
    """Per-energy Sherpa CC MC@NLO normalisation, raw (pre-spin-factor)."""
    tag = beams.Beams("nu", ENERGY).tag
    try:
        return SHERPA_NU_NLO[tag][which]
    except KeyError:
        raise SystemExit(
            f"no Sherpa CC MC@NLO normalisation for {tag} -- add it to "
            f"SHERPA_NU_NLO in analyze_nu.py, from that energy's own "
            f"integ.log and job spread") from None


# RESOLVED ONLY AT AN ENERGY THE TABLE HAS.  These used to be evaluated
# unconditionally at import, so ANY analysis at an energy outside the earlier scan
# -- a Herwig or POWHEG sample at 700 GeV, which has nothing to do with
# Sherpa -- died with "no Sherpa CC MC@NLO normalisation".  The sherpa_nlo
# branch calls _sherpa_nu_nlo() itself when they are None, so a Sherpa the earlier production
# analysis at such an energy still refuses exactly as before.
if beams.Beams("nu", ENERGY).tag in SHERPA_NU_NLO:
    SHERPA_NU_NLO_SIGMA_RAW_PB = _sherpa_nu_nlo(0)
    SHERPA_NU_NLO_SIGMA_RAW_ERR_PB = _sherpa_nu_nlo(1)
    SHERPA_NU_NLO_SIGMA_PB = SHERPA_NU_SPIN_FACTOR * SHERPA_NU_NLO_SIGMA_RAW_PB
    SHERPA_NU_NLO_SIGMA_ERR_PB = (SHERPA_NU_SPIN_FACTOR
                                  * SHERPA_NU_NLO_SIGMA_RAW_ERR_PB)
    SHERPA_NU_NLO_INTEG_RAW_PB = _sherpa_nu_nlo(2)
else:
    SHERPA_NU_NLO_SIGMA_RAW_PB = SHERPA_NU_NLO_SIGMA_RAW_ERR_PB = None
    SHERPA_NU_NLO_SIGMA_PB = SHERPA_NU_NLO_SIGMA_ERR_PB = None
    SHERPA_NU_NLO_INTEG_RAW_PB = None

# Integrator value per run, for the delivery-closure gate.  These samples are
# deliberately normalised to the EVENT-LEVEL cross-section rather than to the
# integrator (see note (2) above), so the gate is what actually compares the
# two -- it is the tripwire that would have caught RESPECT_MASSIVE_FLAG on day
# one.  Raw, i.e. before the factor-2 helicity correction, because that is what
# the event weights carry.
SHERPA_NU_INTEG_RAW_PB = {
    "sherpa_lo": SHERPA_NU_LO_INTEG_RAW_PB,
    "sherpa_lo_me": 1.59691,
    "sherpa_nlo": SHERPA_NU_NLO_INTEG_RAW_PB,
}

# The CC shower-model arms.  Mirrors SHERPA_SHOWER_ARMS in analyze.py; the run
# directory for the NLO one is NuDIS_NLO_ckm3, NOT NuDIS_NLO -- the latter is
# the superseded unit-CKM run kept as documentation, and pointing an arm at it
# would compare the shower against a sample with a known CKM defect.
SHERPA_NU_SHOWER_ARMS = {
    "sherpa_nlo_dire": ("NuDIS_NLO_ckm3_Dire", True,
                        "Sherpa 3.0.5 MC@NLO, Dire"),
    "sherpa_lo_dire":  ("NuDIS_LO_Dire", False, "Sherpa 3.0.5 LO, Dire"),
}


# --- POWHEG-V2 nu-DIS normalisation ----------------------------------------
# POWHEG-V2 (software/powheg-cmass/nu-DIS-master), card
# powheg/cards/POWHEG-V2/nu1TeV-prod/powheg.input, showered with Pythia 8.311
# through powheg/main_powheg.  NB the project distinguishes POWHEG-V2 from
# POWHEG-RES (powheg-dis-main/DIS_v) -- see README; the muon NC NLO entry is
# the RES code, this neutrino one is V2.
# Integrator sigma from prod-nu1TeV/pwg-stat.dat, "total (btilde+remnants)".
# The GENERATION region (Q2 > 2.25 on the Born, NO y cut) is WIDER than the
# fiducial region, so sigma_fid = integrator sigma x fiducial weight fraction,
# exactly as for the muon POWHEG-RES entry in analyze.py.  No neutrino spin
# factor here: POWHEG-V2 takes a single-helicity neutrino beam correctly.
# 2026-08-21: REGENERATED after the CKM fix (patches/powheg-v2-nu-dis-ckm.diff).
# fun_ckm2 in Born.f had no branch for an incoming b and returned an
# UNINITIALISED value, and gave an outgoing b a Cabibbo weight (|V_cs|^2 for
# c -> b) instead of |V_cb|^2.  Both only bite with numflav 5, which this
# benchmark sets.  Effect of the fix, measured:
#     incoming-b events        6.25%  ->  0.00%
#     negative-weight fraction 7.41%  ->  2.66%
#     |negative| weights       0.388  ->  0.125 pb
#     total cross section      4.46801 -> 4.45357 pb   (-0.32%)
#     Pythia aborts in the shower  92  ->  11
# The old sample is kept in powheg/old_nu_buggy_ckm/ and is NOT globbed.
POWHEG_V2_NU_SIGMA_PB = 4.45357
POWHEG_V2_NU_SIGMA_ERR_PB = 0.00797

# --- POWHEG-RES nu-DIS normalisation ---------------------------------------
# The CROSS-CHECK partner of the entry above: the SAME neutrino CC process
# through the OTHER code (powheg-dis-main/DIS_v, card
# powheg/cards/POWHEG-RES/nu1TeV-CC), so the neutrino side has two independent
# NLO matchings from the two POWHEG families, as the muon side now does.
# From parallel-nu1TeV-CC/pwg-0001-st3-stat.dat, "grand total (pos.-|neg.|)":
# btilde 4.0178 + remnants 0.4264 = 4.44424 +- 0.00250 pb, i.e. 0.22% from
# POWHEG-V2's 4.45357 on the same wide generation region.
#
# >>> INCLUSIVE SELECTION ONLY.  THIS SAMPLE HAS NO CHARM ROW. <<<
# DIS_v/init_processes.f hardcodes CKM_diag = .true. in three places with no
# input keyword, so the CKM matrix never reaches the matrix elements and
# d -> c is absent altogether.  CKM unitarity protects the INCLUSIVE rate --
# which is exactly why it agrees with POWHEG-V2 to 0.22% -- and does not
# protect the charm composition.  The identical defect in Sherpa cost 3.4
# points on the charm fraction (18.69% against the correct 22.10%), and there
# is no equivalent patch here.  analyse_nu refuses the charm variant rather
# than producing a number that would sit in a table looking comparable.
POWHEG_RES_NU_SIGMA_PB = 4.44424
POWHEG_RES_NU_SIGMA_ERR_PB = 0.00250

NU_BEAM_PID = 14      # nu_mu beam
LEP_OUT_PID = 13      # mu- scattered lepton

# Herwig cross-section bookkeeping is shared with the muon side so that the
# two cannot drift apart again: this file used to hold its own copy of the
# rule while herwig7/make_histos.py used the WRONG line of Herwig's table for
# months.  The rule, and the evidence for it, are in herwig7/herwig_xsec.py.
sys.path.insert(0, os.path.join(BASE, "herwig7"))
from herwig_xsec import combine as hw_combine  # noqa: E402
from herwig_xsec import read_jobs as hw_read_jobs  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp
import sample_layout  # noqa: E402  -- where the paper-plots samples live


def pwg_runname(stem):
    """Herwig run name for the POWHEG samples at the ACTIVE energy.

    The job directory alone cannot scope these: "nupwg_job_*" matches
    nupwg_job_400GeV_* perfectly happily.  The run name can, because Herwig
    names the run after the card's saverun -- "DIS-nu-PWG-S*" cannot match
    "DIS-nu-PWG-400GeV-S*".  Same discriminator the LO branch uses, and the
    same trap that pooled two energies' cross-sections in the genie_lo branch
    until it was caught by the per-job spread warning.
    """
    return (stem if ENERGY == beams.ANCHOR_ENERGY
            else f"{stem}-{beams.Beams('nu', ENERGY).tag}")


def reconstruct_beams(nu3):
    """Beams for events without status-4 records (Sherpa BEAM_REMNANTS: false):
    the beam neutrino is the hardest status-3 nu_mu, and the proton beam is
    back-to-back with it in the generation (c.m.) frame."""
    e_p = math.sqrt(nu3[1]**2 + nu3[2]**2 + nu3[3]**2 + M_P*M_P)
    return nu3, (e_p, -nu3[1], -nu3[2], -nu3[3])


def parse_hepmc3(fname, run_stats=None, hard_graph=False):
    """Yield (weight, nu_beam, p_beam, particles, dmesons, hard) per
    event; momenta (E,px,py,pz).  particles = list of (pid, p4) for
    final-state (status 1); dmesons = de-duplicated (pid, p4) list of D
    mesons counted AT PRODUCTION; hard = {"in", "out", "chad"}: the PDG codes
    of the hard-process incoming/outgoing records (Pythia status 21/23, or the
    graph walk for generators that write neither), plus "chad" -- whether the
    event contains ANY charm hadron after shower and hadronisation, which is
    the benchmark's final-state-charm tag (see analyze.is_charm_hadron).

    TRAP (identical to the muon side, see analyze.py): with the ctau > 10 mm
    stability convention a D meson always decays, so there is not a single
    status-1 D in any sample -- they must be taken from the status-2 records
    and de-duplicated with dedup_dmesons(), which needs the vertex incoming
    lists.  A final-state filter returns exactly zero everywhere."""
    w = None
    nu_beam = p_beam = nu3 = None
    parts = []
    dcand = []          # every D-flavoured record, any status
    chad = []           # every prompt charm-HADRON record (hard_graph only)
    vin = {}            # vertex id -> raw "id,id,..." of its incoming particles
    hard_out, hard_in = [], []
    # Heavy flavour in the final state: fs_heavy holds 4 when the event has
    # PROMPT charm and 5 when it has bottom at all.  The b -> c cascade is
    # excluded by ancestry -- a bottom hadron always decays to charm, so a
    # naive test would make every b event a charm event.  Identical machinery
    # and identical reasoning to analyze.py; see the long note there.  Lists
    # so the parse loop can write to them without a nonlocal.
    fs_heavy = [set()]
    b_tainted = [set()]     # particle ids: a bottom hadron or its descendant
    v_tainted = [set()]     # vertex ids with a tainted incoming particle
    pall, vout, vgraph = {}, {}, {}
    lhe_index = [None]
    have_event = False

    def finish():
        ho, hi = hard_out, hard_in
        if hard_graph and not ho:
            # Sherpa / Herwig write no Pythia status-21/23 records; walk the
            # graph instead.  CC, so the exchanged boson is a W.
            hi, ho = find_hard_vertex(pall, vgraph, vout, NU_BEAM_PID,
                                      CC_BOSONS)
        ho = {"in": hi, "out": ho, "fs_heavy": frozenset(fs_heavy[0]),
              "lhe_index": lhe_index[0]}
        if hard_graph:
            ho["charm_final"] = weakly_decaying_charm(chad, vin)
        if nu_beam is not None and p_beam is not None:
            return w, nu_beam, p_beam, parts, dedup_dmesons(dcand, vin), ho
        if nu_beam is None and p_beam is None and nu3 is not None:
            return ((w,) + reconstruct_beams(nu3)
                    + (parts, dedup_dmesons(dcand, vin), ho))
        return None

    with open(fname) as f:
        for line in f:
            c0 = line[0:2]
            if c0 == "E ":
                if have_event:
                    ev = finish()
                    if ev is not None:
                        yield ev
                w, nu_beam, p_beam, nu3, parts = None, None, None, None, []
                dcand, vin, hard_out, hard_in = [], {}, [], []
                chad = []
                fs_heavy, b_tainted, v_tainted = [set()], [set()], [set()]
                pall, vout, vgraph = {}, {}, {}
                lhe_index = [None]
                have_event = True
            elif c0 == "A " and line.startswith("A 0 lhe_index "):
                # main_powheg's Les Houches event number: the join key for
                # every reweighting (memory: lhe-index-and-wrong-lhe), so it
                # travels with the event rather than being re-derived
                lhe_index[0] = int(line.split()[3])
            elif c0 == "W " and w is None:
                f_ = line.split()
                try:
                    w = float(f_[1])
                except ValueError:
                    pass          # the "W Weight|EXTRA__..." header line
                else:
                    if run_stats is not None:
                        run_stats["sum_w"] = run_stats.get("sum_w", 0.0) + w
                        if len(f_) > 4:      # Sherpa: NTrials is column 4
                            try:
                                nt = float(f_[4])
                            except ValueError:
                                pass
                            else:
                                run_stats["sum_ntrials"] = (
                                    run_stats.get("sum_ntrials", 0.0) + nt)
                                run_stats["n_with_trials"] = (
                                    run_stats.get("n_with_trials", 0) + 1)
            elif c0 == "P ":
                c = line.split()
                if hard_graph:
                    pall[int(c[1])] = (int(c[3]), int(c[9]))
                    par = int(c[2])
                    if par < 0:
                        vout.setdefault(par, []).append(int(c[1]))
                pid, status = int(c[3]), int(c[9])
                if status == 1 or status == 4:
                    p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                    if status == 4:
                        if pid == NU_BEAM_PID:
                            nu_beam = p
                        elif pid == 2212 or pid == 2112:
                            # 2112: the FASER-rate ladders scatter off a free
                            # NEUTRON too (tools/genie_dimuon_ladder.sh); the
                            # kinematics use M_P, 0.14% off, immaterial there
                            p_beam = p
                    else:
                        parts.append((pid, p))
                elif status == 3 and pid == NU_BEAM_PID:
                    p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                    if nu3 is None or p[0] > nu3[0]:
                        nu3 = p
                elif status == 23:
                    hard_out.append(pid)
                elif status == 21:
                    hard_in.append(pid)
                # deliberately outside the status branches: charm hadrons
                # are always decayed (status 2), never final state, under the
                # ctau > 10 mm convention -- a status-1 filter finds none
                if is_heavy_hadron(pid, 5) or (status == 1 and abs(pid) == 5):
                    fs_heavy[0].add(5)
                    b_tainted[0].add(int(c[1]))
                else:
                    tainted = False
                    if b_tainted[0]:
                        par_ = int(c[2])
                        tainted = ((par_ in v_tainted[0]) if par_ < 0
                                   else (par_ in b_tainted[0]))
                        if tainted:
                            b_tainted[0].add(int(c[1]))
                    if not tainted and 4 not in fs_heavy[0] and (
                            is_heavy_hadron(pid, 4) or
                            (status == 1 and abs(pid) == 4)):
                        fs_heavy[0].add(4)
                if pid in D_PIDS_SIGNED:
                    dcand.append(
                        (pid, (float(c[7]), float(c[4]), float(c[5]),
                               float(c[6])), int(c[1]), int(c[2])))
                # EVERY prompt charm hadron with its momentum, for the charm
                # ANATOMY (2026-09-01).  Same block as analyze.py, kept in
                # step by tools/check_parser_drift.py.  Opt-in via hard_graph
                # so nothing on the production path changes.
                if hard_graph and is_heavy_hadron(pid, 4) \
                        and int(c[1]) not in b_tainted[0]:
                    chad.append(
                        (pid, (float(c[7]), float(c[4]), float(c[5]),
                               float(c[6])), int(c[1]), int(c[2])))
            elif c0 == "V ":
                c = line.split()
                raw = c[3][1:-1]
                vin[int(c[1])] = raw
                if b_tainted[0] and any(int(t) in b_tainted[0]
                                        for t in raw.split(",") if t):
                    v_tainted[0].add(int(c[1]))
                if hard_graph:
                    vgraph[int(c[1])] = [int(t) for t in raw.split(",") if t]
        if have_event:
            ev = finish()
            if ev is not None:
                yield ev


def analyze(files, charm_tag=False, tag_stats=None, run_stats=None,
            charm_mode="out", subtract_shower_charm=False):
    """charm_tag: keep only charm-production events.  charm_mode picks what
    that means -- "final" (the benchmark definition since 2026-08-21: at least
    one charm quark in the FINAL state, i.e. a charm hadron after shower and
    hadronisation), or the two hard-process tags "out" (charm outgoing from the
    matrix element, the CC channel W + d/s/b -> c) and "any" (charm on either
    side of it).  All three are counted in tag_stats whichever is selected, so
    a run reports how far apart they are.  n_events still counts ALL parsed events so the caller's
    sigma_generated / N_generated normalisation gives the channel's absolute
    contribution to sigma_fid.  tag_stats: optional dict, filled with the
    tagged-event count and the D-meson capture fraction of the tag."""
    obs = {k: [] for k in BINS}
    obs["nd_s"] = []      # D_s multiplicity: reported as a mean, not binned
    # THE LES HOUCHES EVENT NUMBER OF EVERY SURVIVING EVENT.  Not a histogram
    # axis -- it is the join key that lets a scale- or PDF-reweighted LHE be
    # married to the SHOWERED events, which is the only way to put a band on
    # an observable the LHE cannot compute (a hadron multiplicity, a Delta phi
    # against the hadron system).  -1 where the sample carries none, as GENIE
    # and Sherpa do.  Same precedent as nd_s above: an obs key that no BINS
    # entry consumes.
    obs["lhe_index"] = []
    weights = []
    n_events = 0
    sum_w_all = 0.0
    if tag_stats is not None:
        tag_stats.setdefault("n_tagged", 0)
        tag_stats.setdefault("n_d_all", 0)
        tag_stats.setdefault("n_d_tagged", 0)
    for fn in files:
        print(f"  parsing {fn}", flush=True)
        want_hard = charm_tag or tag_stats is not None or subtract_shower_charm
        for w, k, P, parts, dmes, hard in parse_hepmc3(
                fn, run_stats, hard_graph=want_hard):
            if w is None:
                w = 1.0
            n_events += 1
            sum_w_all += w
            # shower-only charm (Sherpa): see analyze.SHOWER_CHARM_SUBTRACTED_PREFIXES
            shower_only = False
            raw_final_event = False
            if subtract_shower_charm and isinstance(hard, dict):
                _hin, _hout = hard.get("in", []), hard.get("out", [])
                shower_only = (4 in hard.get("fs_heavy", frozenset())
                               and not any(abs(q) == 4 for q in _hin)
                               and not any(abs(q) == 4 for q in _hout))
            if charm_tag or tag_stats is not None:
                # CHARM TAG MODES.
                #  "final" -- a charm HADRON anywhere in the event record,
                #           i.e. charm in the final state after shower and
                #           hadronisation.  THE BENCHMARK DEFINITION.
                #  "out" -- charm in the hard-process FINAL state (s/d -> c).
                #  "any" -- charm ANYWHERE in the hard process, initial or
                #           final.  This is the apples-to-apples match to the
                #           YADISM "all terms proportional to V_cd/V_cs"
                #           definition, which counts nu cbar -> mu- sbar.
                #
                # PHYSICS NOTE (user, 2026-08-20): the proton carries ZERO net
                # charm number, so an incoming cbar from the sea leaves its
                # partner c in the REMNANT.  Such an event therefore does have
                # a charm quark in the final state of the FULL EVENT -- it is
                # only absent from the hard-process final state.  So "out" and
                # "any" differ at hard-process level, while a tag on
                # final-state charm HADRONS would see both.  That third
                # definition also picks up shower g -> c cbar, which is a
                # different process again.
                hin = hard.get("in", []) if isinstance(hard, dict) else []
                hout = hard.get("out", []) if isinstance(hard, dict) else hard
                tagged_out = any(abs(p) == 4 for p in hout)
                tagged_any = tagged_out or any(abs(p) == 4 for p in hin)
                # PROMPT charm: the b -> c cascade was already removed by
                # the ancestry pass in parse_hepmc3().  Diagnosed on the muon
                # POWHEG-RES sample, where the naive test added 1338 events of
                # which 1219 were b-INITIATED, carrying net NEGATIVE weight.
                heavy = (hard.get("fs_heavy", frozenset())
                         if isinstance(hard, dict) else frozenset())
                tagged_final = 4 in heavy
                tagged = {"any": tagged_any, "final": tagged_final}.get(
                    charm_mode, tagged_out)
                if charm_mode == "final" and subtract_shower_charm \
                        and shower_only and tagged_final:
                    # Sherpa: a shower-only charm event is dropped from the
                    # published result but counted in the raw sum -- after
                    # the selection cuts, so the raw number is the raw number
                    # of THIS selection
                    raw_final_event = True
                if tag_stats is not None:
                    tag_stats["n_tagged_out"] = (
                        tag_stats.get("n_tagged_out", 0) + int(tagged_out))
                    tag_stats["n_tagged_any"] = (
                        tag_stats.get("n_tagged_any", 0) + int(tagged_any))
                    tag_stats["n_tagged_final"] = (
                        tag_stats.get("n_tagged_final", 0) + int(tagged_final))
                    tag_stats["n_charm_with_bottom"] = (
                        tag_stats.get("n_charm_with_bottom", 0)
                        + int(4 in heavy and 5 in heavy))
                if tag_stats is not None:
                    tag_stats["n_d_all"] += len(dmes)
                    if tagged:
                        tag_stats["n_tagged"] += 1
                        tag_stats["n_d_tagged"] += len(dmes)
                if charm_tag and not tagged:
                    continue
            # hardest final-state mu- by invariant lab energy
            best, best_elab = None, -1.0
            for pid, p in parts:
                if pid == LEP_OUT_PID:
                    elab = lab_energy(p, P)
                    if elab > best_elab:
                        best, best_elab = p, elab
            if best is None:
                continue
            # through the SAME helpers the muon side uses -- not a second
            # copy of them.  This block WAS a second copy, byte for byte.
            Q2, y, xbj, kP = dis_invariants(k, P, best)
            # The W^2 branch is here for the same reason it is in analyze.py:
            # rule 2b asks that a capability added on one current exist on the
            # other, so that a future charged-current comparison against a
            # (Q, W) region costs a selection and not a code change.  W2_MIN is
            # None for every selection the benchmark publishes.
            if W2_MIN is not None and M_P*M_P + 2.0*y*kP - Q2 < W2_MIN:
                continue
            if Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX:
                continue
            e_lab = best_elab
            theta = lepton_theta(e_lab, best)
            # FASER tiers, kept identical to the muon side (CONVENTIONS.md rule 2b)
            if not SELECTION.passes_lepton(e_lab, theta):
                continue
            # The DIMUON tier, through the SAME helper the muon side uses --
            # not a second copy of it.  Here the primary is the mu- of the CC
            # vertex and the partner is a mu+ from charm decay, which is the
            # historical dimuon charm signal in neutrino DIS.
            if SELECTION.dimuon and not SELECTION.passes_dimuon(
                    n_opposite_sign_muons(parts, P, LEP_OUT_PID, SELECTION)):
                continue
            if shower_only:
                # Sherpa: charm from its massless shower is not counted, nor
                # is a dimuon it decays into (analyze.py says why)
                if SELECTION.dimuon:
                    continue
                dmes = []
            ebeam_lab = kP / M_P

            # lab azimuth of the scattered lepton; a z-boost leaves phi alone,
            # so the generation-frame components serve
            phi_mu = math.atan2(best[2], best[1])

            nch = nch1 = 0
            elead = 0.0
            # Tier E: charged-hadron angular counts and the summed hadron
            # system.  The LAB angle is rebuilt from invariants (p_T is
            # boost-invariant, E_lab = p.P/M_p, m^2 = p.p), because theta is
            # NOT invariant under the boost from the generation frame.
            n05 = n01 = 0
            nch05 = 0
            # THE NUCLEAR OBSERVABLES, identical to the muon side (rule 2b).
            # Counted over EVERY final-state particle, not only the charged
            # hadrons, since a neutron is neither; Ehad sums the lab energy of
            # everything that is not a neutrino and the primary lepton is
            # removed afterwards.  See analyze.py's bins_for() for why
            # antiparticles are counted with particles.
            nprot = nneut = npi = 0
            ehad = 0.0
            # dphi_min: the closest single charged hadron in azimuth.  ADDED
            # 2026-09-08 -- the muon parser has had it since the beginning and
            # this one never did, so every dphi figure was muon-only, which is
            # the rule-2b asymmetry that is invisible in any single number.
            dphi_min = math.inf
            sum_px = sum_py = 0.0
            for pid, p in parts:
                a = abs(pid)
                if a == 2212:
                    nprot += 1
                elif a == 2112:
                    nneut += 1
                elif a == 211:
                    npi += 1
                if a not in NEUTRINO_PIDS:
                    ehad += dot(p, P) / M_P
                if abs(pid) in CHARGED_HADRONS:
                    elab_h = dot(p, P) / M_P
                    nch += 1
                    if elab_h > 1.0:
                        nch1 += 1
                    if elab_h > elead:
                        elead = elab_h
                    # sign of p_z,lab and the 1 GeV track threshold: see
                    # analyze.py's copy of this count (2026-09-19, review S1)
                    _pt, tan_h = track_tan(p, elab_h, P, k)
                    if tan_h < 0.5 and elab_h > TRACK_E_MIN:
                        n05 += 1
                        nch05 += 1
                        if tan_h < 0.1:
                            n01 += 1
                    sum_px += p[1]
                    sum_py += p[2]
                    d = abs(math.remainder(math.atan2(p[2], p[1]) - phi_mu,
                                           2.0 * math.pi))
                    if d < dphi_min:
                        dphi_min = d

            if sum_px or sum_py:
                dphi_x = abs(math.remainder(
                    math.atan2(sum_py, sum_px) - phi_mu, 2.0 * math.pi))
            else:
                dphi_x = None
            if not SELECTION.passes_hadrons(n05, n01, dphi_x):
                continue

            # D mesons at production; lab energies via the invariant p.P/M_p.
            # Ed_* is the LEADING D of the event (one entry per event, so the
            # per-event weight array stays valid).  D_s is counted but NOT
            # folded into the charged-D observable: AHADIC makes ~2x the D_s of
            # the other generators, which would import a hadronisation
            # difference into a channel-rate plot.
            nd_ch = nd_0 = nd_s = 0
            ed_ch = ed_0 = 0.0
            for pdg, p in dmes:
                a = abs(pdg)
                if a == D_CH_PID:
                    nd_ch += 1
                    e_d = dot(p, P) / M_P
                    if e_d > ed_ch:
                        ed_ch = e_d
                elif a == D_0_PID:
                    nd_0 += 1
                    e_d = dot(p, P) / M_P
                    if e_d > ed_0:
                        ed_0 = e_d
                else:
                    nd_s += 1

            vals = {"Q2": Q2, "xbj": xbj, "y": y, "Emu": e_lab,
                    "theta": theta, "nch": nch, "nch1": nch1,
                    "nch05": nch05,
                    "dphi": dphi_min if nch else NO_DPHI,
                    "dphix": dphi_x if dphi_x is not None else NO_DPHI,
                    "nu": ebeam_lab - e_lab, "Elead": elead,
                    "nd_ch": nd_ch, "nd_0": nd_0, "nd_s": nd_s,
                    "Ed_ch": ed_ch, "Ed_0": ed_0,
                    "nprot": nprot, "nneut": nneut, "npi": npi,
                    "Ehad": ehad - e_lab}
            # >>> THE SUBTRACTED EVENT MUST NOT REACH `obs` EITHER. <<<
            # This block used to sit AFTER the loop that appends to obs, so a
            # Sherpa shower-only charm event was recorded in every observable
            # list and then skipped for `weights`.  The arrays then had
            # different lengths AND, worse, were MISALIGNED from the first
            # subtracted event on: obs[i] no longer belonged to weights[i].
            # np.average refuses shapes that differ, which is the only reason
            # this was ever visible -- it made the charm-tagged Sherpa results
            # fail to be written rather than be written wrong, and the failure
            # was swallowed by a driver that greps its output.  Found on
            # 2026-09-08 while asking why sherpa_dire_charmfinal had no Tier S
            # or Tier E result; introduced with the subtraction itself
            # (24407c1, 2026-09-03).  CONVENTIONS.md rule 2, once more.
            if raw_final_event:
                if tag_stats is not None:
                    tag_stats["sum_w_final_raw"] = tag_stats.get("sum_w_final_raw", 0.0) + w
                continue
            for kk, v in vals.items():
                if kk in obs:
                    obs[kk].append(v)
            li = hard.get("lhe_index") if isinstance(hard, dict) else None
            obs["lhe_index"].append(-1 if li is None else int(li))
            weights.append(w)
            if tag_stats is not None and charm_tag and charm_mode == "final":
                tag_stats["sum_w_final_raw"] = tag_stats.get("sum_w_final_raw", 0.0) + w
    return obs, np.asarray(weights), n_events, sum_w_all


def powheg_v2_nu_sigma_pb():
    """POWHEG-V2 integrator cross-section at the configured energy, in pb.

    Read from that run's pwg-stat.dat, "total (btilde+remnants)".  The ANCHOR
    keeps its literal: the committed 4.45357 does not equal what the current
    stat file says (4.46801), so the literal carries history the file no
    longer does, and replacing it would move a published number.

    Without this the delivery-closure gate would compare a 400 GeV sample
    against the 1 TeV integrator, which is how the same trap was caught six
    times already in this scan.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return POWHEG_V2_NU_SIGMA_PB
    tag = beams.Beams("nu", ENERGY).tag
    f = f"{paths.POWHEG_V2}/prod-nu{tag}/pwg-stat.dat"
    if not os.path.exists(f):
        raise SystemExit(f"no POWHEG-V2 statistics at {f} -- run "
                         f"powheg/run_powheg_v2.sh {ENERGY:g} first")
    for line in open(f, errors="replace"):
        m = re.search(r"total \(btilde\+remnants\) cross section in pb\s+"
                      r"([0-9.eEdD+-]+)", line)
        if m:
            v = float(m.group(1).replace("D", "E").replace("d", "e"))
            print(f"  [powheg-v2] integrator {v:.5f} pb (from pwg-stat.dat)")
            return v
    raise SystemExit(f"could not read the total cross section from {f}")


def powheg_v2mc_sigma_pb():
    """POWHEG-V2mc integrator cross-section at the configured energy, in pb.

    A SEPARATE reader, not a flag on the one above, because it reads a
    different run directory and a mix-up between the two would be invisible:
    both are POWHEG-V2 numbers of the same order for the same beam.
    """
    tag = beams.Beams("nu", ENERGY).tag
    f = f"{paths.POWHEG_V2}/prod-numc{tag}/pwg-stat.dat"
    if not os.path.exists(f):
        raise SystemExit(f"no POWHEG-V2mc statistics at {f} -- run "
                         f"powheg/run_powheg_v2.sh --mc {ENERGY:g} first")
    for line in open(f, errors="replace"):
        m = re.search(r"total \(btilde\+remnants\) cross section in pb\s+"
                      r"([0-9.eEdD+-]+)", line)
        if m:
            v = float(m.group(1).replace("D", "E").replace("d", "e"))
            print(f"  [powheg-v2mc] integrator {v:.5f} pb (from pwg-stat.dat)")
            return v
    raise SystemExit(f"could not read the total cross section from {f}")


def powheg_res_nu_sigma_pb():
    """POWHEG-RES neutrino CC integrated cross-section for this energy, in pb.

    Read as `grand total pos. weights - grand total |neg.| weights` from seed
    0001's stage-3 statistics, the same rule and the same file shape as the
    muon POWHEG-RES entry in analyze.py -- each seed integrates the full phase
    space, so one seed is a complete estimate rather than a partial sum.

    The anchor keeps its literal so nothing published can move.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return POWHEG_RES_NU_SIGMA_PB
    tag = beams.Beams("nu", ENERGY).tag
    f = f"{paths.POWHEG_RES}/parallel-nu{tag}-CC/pwg-0001-st3-stat.dat"
    if not os.path.exists(f):
        raise SystemExit(f"no POWHEG-RES statistics at {f} -- run "
                         f"powheg/run_powheg_res.sh --current nu {ENERGY:g} "
                         f"first")
    pos = neg = None
    for line in open(f, errors="replace"):
        m = re.search(r"grand total pos\.\s+weights:\s+([0-9.eE+-]+)", line)
        if m:
            pos = float(m.group(1))
        m = re.search(r"grand total \|neg\.\|\s+weights:\s+([0-9.eE+-]+)",
                      line)
        if m:
            neg = float(m.group(1))
    if pos is None or neg is None:
        raise SystemExit(f"could not read the grand totals from {f}")
    sig = pos - neg
    print(f"  [powheg-res-nu] integrator {sig:.5f} pb "
          f"(pos {pos:.5f} - |neg| {neg:.5f}, from {os.path.basename(f)})")
    return sig


# Where each neutrino POWHEG entry's LHE files live, for the closure gate's
# denominator (see analyze.delivered_sigma_pb).
POWHEG_NU_LHE_DIR = {
    # 1 TeV LIVES IN prod-nu1TeV-ckm, exactly as powheg_variants.sh says.
    # The shell side was repointed on 2026-08-28 and THIS COPY WAS NOT, so the
    # closure gate went on reading the superseded 2026-08-13 directory's
    # .nevents sidecar.  It was invisible while both held 50000 events -- the
    # gate got the right denominator from the wrong file -- and surfaced only
    # when the sample was tripled and the gate reported +194%.
    # A second copy of a path is how this benchmark keeps breaking; the
    # condition is written the same way in both places on purpose.
    "powheg_nu": lambda t: (f"{paths.POWHEG_V2}/prod-nu1TeV-ckm"
                            if t == "1TeV" and os.path.isdir(
                                f"{paths.POWHEG_V2}/prod-nu1TeV-ckm")
                            else f"{paths.POWHEG_V2}/prod-nu{t}"),
    "powheg_nu_mc": lambda t: f"{paths.POWHEG_V2}/prod-numc{t}",
    "powheg_res_nu": lambda t: f"{paths.POWHEG_RES}/parallel-nu{t}-CC",
    # paper plots: one integration directory per nucleon, from sample_layout
    # (the tag argument is ignored: sample_layout derives it from ENERGY itself)
    "v2_powheg_nu": lambda t: sample_layout.powheg_rundir(
        "nu", sample_layout.require()[1], ENERGY),
    "v2_powheg_nu_mc": lambda t: sample_layout.powheg_rundir(
        "nu", sample_layout.require()[1], ENERGY, mc=True),
}


def qed_arm_nu(gen):
    """The arm for a `powheg_res_nu_qed_<arm>` key, or None.

    The CHARGED CURRENT HAS NO `fsrisr` ARM: initial-state radiation off the
    lepton line needs a charged incoming lepton and the beam here is a
    neutrino.  That is refused rather than silently treated as `fsr`, because
    the two would then be indistinguishable in the results and the asymmetry
    between the currents -- which is a finding of this study -- would look
    like a duplicate row.
    """
    pre = "powheg_res_nu_qed_"
    if not gen.startswith(pre):
        return None
    arm = gen[len(pre):].replace("_charm", "")
    if arm == "fsrisr":
        sys.exit("powheg_res_nu_qed_fsrisr does not exist: ISR off the lepton "
                 "line\n  needs a CHARGED incoming lepton, and the charged "
                 "current starts from a\n  neutrino.  The muon side has this "
                 "arm and the neutrino side cannot;\n  that asymmetry is a "
                 "result of the QED study, not a gap in it.")
    if arm not in QED_ARMS:
        sys.exit(f"no such QED arm {arm!r} -- the charged current has "
                 f"{'fsr'!r} and {'full'!r}")
    return arm


def sherpa_nu_lo_sigma_pb():
    """Sherpa CC LO integrator cross-section at the configured energy, in pb.

    The anchor keeps its literal, which carries the measured event-level value
    and its history.  Other energies read that run's own integ.log, exactly as
    the muon side does in analyze.sherpa_lo_sigma_pb().

    THE SPIN FACTOR STILL APPLIES.  Sherpa halves the neutrino CC cross-section
    (it averages over two incoming lepton helicities where only one couples),
    so the integrator number must be multiplied by SHERPA_NU_SPIN_FACTOR here
    as it is in the literal above -- forgetting it would look like a clean
    factor-two disagreement with everyone else.
    """
    if ENERGY == beams.ANCHOR_ENERGY:
        return SHERPA_NU_LO_SIGMA_PB
    log = f"{SHERPA_RUNS}/{at_energy('NuDIS_LO')}/integ.log"
    if not os.path.exists(log):
        raise SystemExit(f"no Sherpa CC integration log at {log}")
    ansi = re.compile(r"\x1b\[[0-9;]*m")
    pat = re.compile(r":\s*([0-9.eE+-]+)\s*pb\s*\+-")
    hit = None
    with open(log, errors="replace") as f:
        for line in f:
            m = pat.search(ansi.sub("", line))
            if m:
                hit = m
    if hit is None:
        raise SystemExit(f"could not find the integrated cross-section in {log}")
    raw = float(hit.group(1))
    print(f"  [sherpa nu] integrator {raw:.5f} pb raw "
          f"x {SHERPA_NU_SPIN_FACTOR} spin = {raw * SHERPA_NU_SPIN_FACTOR:.5f} pb")
    return raw * SHERPA_NU_SPIN_FACTOR



def nuwro_sigma_and_n(files, generator="nuwro"):
    """(sigma_pb, N_generated) for a converted sample, from its own sidecars.

    Shared by NuWro and GiBUU: both are converted into HepMC3 by a script in
    this repository, and both of those write the same sidecar beside the file
    they produce.

    >>> THE COUNT MATTERS AS MUCH AS THE CROSS-SECTION, because the delivery
    closure gate cannot help here. <<<  That gate compares the delivered sum
    of weights against the integrator; nuwro2hepmc DEFINES the weight as
    sigma/N, so the two agree by construction whatever happened to the file.
    A truncated HepMC would sail through it.  What catches that instead is
    this count, checked against the events actually parsed.

    Jobs are pooled by ADDING cross-section-weighted counts, not by averaging
    the cross-sections: two jobs of a NuWro sample differ only in seed, so
    their sigmas agree to the integration error, and what the fiducial
    fraction needs is the total generated.
    """
    sig, n = None, 0
    for f in files:
        j = f.replace(".hepmc", "_xsec.json")
        if not os.path.exists(j):
            sys.exit(f"no {j} -- rerun nuwro/nuwro2hepmc for that job")
        with open(j) as fh:
            d = json.load(fh)
        if d.get("generator") != generator:
            sys.exit(f"{j} says generator {d.get('generator')!r}, not "
                     f"{generator!r} -- that sample is not what was asked for")
        s_j = float(d["sigma_pb"])
        if sig is None:
            sig = s_j
        elif abs(s_j / sig - 1.0) > 0.02:
            sys.exit(f"{j} quotes sigma = {s_j:g} pb against {sig:g} pb for "
                     f"the first job -- these are not the same run condition")
        n += int(d["n_events"])
    if sig is None:
        return None, 0
    return sig, n


def main():
    # Before anything is normalised: every anchor literal must still
    # agree with the integration file it was taken from.  Cheap (four
    # small files), and it is the check whose absence let a
    # re-integrated anchor be normalised against the superseded
    # number on 2026-08-27.
    check_anchor_literals()
    gen = sys.argv[1] if len(sys.argv) > 1 else "pythia"
    # <gen>_charm      = tag charm OUT of the hard process (the original tag)
    # <gen>_charmany   = tag charm ANYWHERE in the hard process
    # <gen>_charmfinal = tag charm in the FINAL STATE, after shower and
    #                    hadronisation -- the benchmark definition
    charm_mode = "out"
    for _suf, _mode in (("_charmany", "any"), ("_charmfinal", "final")):
        if gen.endswith(_suf):
            gen, charm_mode = gen[:-len(_suf)] + "_charm", _mode
            break
    # not every branch reads per-job sidecars; the manifest tolerates empty
    sig_pb, job_counts = [], []
    charm_tag, tag_stats = False, None
    # UNWEIGHTED sample: sigma_fid = sigma_generation_region x n_fid/n_generated.
    # Set explicitly per branch.  It used to be inferred from a list of
    # generator names, which silently stopped matching as soon as a charm-tag
    # suffix was appended to `gen` a few lines below.
    unweighted = False
    if gen in ("pythia", "pythia_charm"):
        files = job_files(f"{BASE}/pythia8", at_energy("nu_job"), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/pythia8", at_energy("nu_job"), "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = "Pythia 8.311 LO"
        if gen == "pythia_charm":
            # CC charm production W + d/s/b -> c, tagged on the outgoing
            # hard-process charm (status 23) before shower/hadronisation
            charm_tag, tag_stats = True, {}
            label = "Pythia 8.311 LO (charm)"
    elif gen in ("mg5_ps", "mg5_ps_charm"):
        # MG5_aMC LO matrix element, showered and hadronised by Pythia 8
        # (mg5/run_mg5.sh).  Unweighted like the Pythia sample above, and the
        # cross-section is MG5's own integrated weight, which Pythia reports
        # back unchanged through the Les Houches file.
        #
        # THE SAMPLE IS GENERATED INSIDE THE FIDUCIAL REGION: Q2 and y are cut
        # at generation by mg5/dis_hooks.f, since MG5's run card has no DIS
        # cuts at all.  So n_fiducial/n_parsed is 1 by construction, and a
        # number materially below it would mean the shower is moving the
        # leptonic invariants -- which, with QED radiation off and the dipole
        # recoil on, it must not.
        files = job_files(f"{BASE}/mg5", at_energy("ps_nu_job"), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/mg5", at_energy("ps_nu_job"),
                            "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            # THE MEAN DELIVERED WEIGHT, NOT Pythia's sigma_gen.  The events
            # are unweighted and every Les Houches weight IS MG5's integrated
            # cross-section, so their mean is that number exactly.  Pythia's
            # own sigmaGen divides the same sum by the events it TRIED, which
            # includes the handful whose hard process it could not construct
            # (30 in 200000 here, "setting mass failed"), while the fiducial
            # fraction below is taken over the events DELIVERED -- mixing the
            # two denominators would lose that fraction twice.
            sig_pb.append(d["sum_weights_pb"] / d["n_accepted"])
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = "MG5_aMC 3.7.2 LO + Pythia8"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("herwig", "herwig_charm"):
        # unweighted; sigma = the ME integral over the generation cuts, which
        # equal the fiducial cuts (see read_herwig_sigma_pb for why the other
        # number in Herwig's table must NOT be used).  Samples come from
        # DIS-nu-POL.in, so HERWIG_NU_SPIN_FACTOR = 1.
        files = job_files(f"{BASE}/herwig7", at_energy("nu_job"), "events.hepmc")
        # SCOPE THE .out GLOB TO THIS ENERGY.  "nu_job_*/DIS-nu*-S*.out"
        # matches every energy's directories AND run names at once: at
        # 400 GeV it found 8 events.hepmc against 19 .out files, and the
        # mismatch guard below refused the sample.  The run name is the
        # reliable discriminator, since Herwig names the run after the card's
        # saverun -- "DIS-nu-S*" cannot match "DIS-nu-400GeV-S*".
        # ...and the run name is the POLARIZED card's, at EVERY energy.
        # run_herwig_nu.sh uses DIS-nu-POL.in at the anchor and
        # DIS-nu-POL-<tag>.in elsewhere, so Herwig's saverun writes
        # DIS-nu-POL[-<tag>]-S*.out.  The first version of this scoping said
        # "DIS-nu" for the anchor, which matches NOTHING -- the old glob
        # DIS-nu*-S*.out had covered the POL suffix by accident.  It fails
        # loudly on the MISMATCH below rather than silently, but it fails at
        # the anchor too, so keep this in step with the run script's CARD=.
        runname = ("DIS-nu-POL" if ENERGY == beams.ANCHOR_ENERGY
                   else f"DIS-nu-POL-{beams.Beams('nu', ENERGY).tag}")
        jobs = hw_read_jobs(
            f"{BASE}/herwig7/{at_energy('nu_job')}_*/{runname}-S*.out")
        if len(jobs) != len(files):
            sys.exit(f"MISMATCH: {len(files)} events.hepmc but {len(jobs)} "
                     f".out files -- refusing to normalise a mixed sample.")
        sig, err, n_gen_tot, survival = hw_combine(jobs)
        print(f"  herwig nu: survival {survival:.4f} "
              f"({100*(1-survival):.2f}% vetoed and regenerated)")
        sigma_sample_pb = HERWIG_NU_SPIN_FACTOR * sig
        sigma_err_pb = HERWIG_NU_SPIN_FACTOR * err
        unweighted = True
        label = "Herwig 7.3.0 LO"
        if gen.endswith("_charm"):
            # CC: the W changes flavour, so charm is a FINAL-state particle
            # only (s -> c).  The NC self-check "c in AND c out" would reject
            # every event here.
            charm_tag, tag_stats = True, {}
            label = "Herwig 7.3.0 LO (charm)"
    elif gen in ("herwig_nlo_neg", "herwig_nlo_neg_charm"):
        # THE NEGATIVE HALF of Herwig's NLO cross section (DIS-nu-PWGNEG.in,
        # PowhegMEDISCC:Contribution 2).  DISBase splits the NLO weight into
        # max(0,w) and max(0,-w) and generates them as two separate
        # positive-weight samples; Herwig's defaults ship Contribution 1, so
        # one run is only half the answer.  Combined with the positive half by
        # herwig7/combine_nlo.py -- see that script's header.
        files = job_files(f"{BASE}/herwig7", at_energy("nupwgneg_job"), "events.hepmc")
        jobs = hw_read_jobs(
            f"{BASE}/herwig7/{at_energy('nupwgneg_job')}_*/"
            f"{pwg_runname('DIS-nu-PWGNEG')}-S*.out")
        if len(jobs) != len(files):
            sys.exit(f"MISMATCH: {len(files)} events.hepmc but {len(jobs)} "
                     f".out files -- refusing to normalise a mixed sample.")
        sig, err, n_gen_tot, survival = hw_combine(jobs)
        print(f"  herwig nu NLO (negative): survival {survival:.4f}")
        sigma_sample_pb = HERWIG_NU_SPIN_FACTOR * sig
        sigma_err_pb = HERWIG_NU_SPIN_FACTOR * err
        unweighted = True
        label = "Herwig 7.3.0 NLO (POWHEG, negative)"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "Herwig 7.3.0 NLO (POWHEG, negative, charm)"
    elif gen in ("herwig_nlo", "herwig_nlo_charm"):
        # DIS-nu-POWHEG.in: PowhegMEDISCC, which Herwig's defaults create as
        # Herwig::MEChargedCurrentDIS with Contribution 1 -- the SAME class as
        # the LO matrix element, hence the same code path, hence the polarized
        # beam IS honoured and HERWIG_NU_SPIN_FACTOR = 1 as at LO.  Matchbox
        # does NOT honour it: its CC Born came out at 0.4935 of the analytic
        # value and declaring a left-handed beam changed nothing.
        #
        # Bookkeeping is identical to LO because POWHEG's hardest-emission
        # generation delivers UNIT WEIGHTS -- there is no weight tail here, in
        # contrast to the muon NC Matchbox NLO sample (N_eff 6-10%,
        # max|w| ~ 1131).  Generation cuts equal the fiducial cuts.
        files = job_files(f"{BASE}/herwig7", at_energy("nupwg_job"), "events.hepmc")
        jobs = hw_read_jobs(
            f"{BASE}/herwig7/{at_energy('nupwg_job')}_*/"
            f"{pwg_runname('DIS-nu-PWG')}-S*.out")
        if len(jobs) != len(files):
            sys.exit(f"MISMATCH: {len(files)} events.hepmc but {len(jobs)} "
                     f".out files -- refusing to normalise a mixed sample.")
        sig, err, n_gen_tot, survival = hw_combine(jobs)
        print(f"  herwig nu NLO: survival {survival:.4f} "
              f"({100*(1-survival):.2f}% vetoed and regenerated)")
        sigma_sample_pb = HERWIG_NU_SPIN_FACTOR * sig
        sigma_err_pb = HERWIG_NU_SPIN_FACTOR * err
        unweighted = True
        label = "Herwig 7.3.0 NLO (POWHEG)"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "Herwig 7.3.0 NLO (POWHEG, charm)"
    elif gen in ("genie_lo", "genie_lo_charm"):
        # GENIE's CLASSIC LO path: tune G18_02a_00_000, Bodek-Yang DIS built on
        # GRV98LO, AGKY/Pythia8 hadronisation through the benchmark config
        # overlay (genie/run_genie_nucc_lo.sh).  This is the LO GENIE entry;
        # the "genie" branch below is HEDIS/BGR18, which is NLO.
        #
        # The event-generator list is CCDISCHARM, i.e. DIS-CC *and*
        # DIS-CC-CHARM: QPMDISPXSec subtracts the charm piece from the
        # inclusive DIS cross-section internally, so running DIS-CC alone
        # would deliver a charm-free sample with a correspondingly short rate.
        # sigma_gen is the sum over both channels' splines at 1 TeV
        # (genie/spline_to_json.py), and the generation region is the FULL
        # phase space -- much wider than fiducial, hence the n_fid/n_gen
        # scaling below.
        files = job_files(f"{BASE}/genie", at_energy("nucclo_job"), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        # job_files(), NOT a raw glob of "nucclo_job_*".  The events.hepmc line
        # above already used it; this loop did not, and the moment
        # nucclo_job_400GeV_* existed the two disagreed -- the EVENTS came from
        # one energy while the CROSS-SECTION was pooled over both.  The
        # anchor was equally exposed, since "nucclo_job_*" matches the tagged
        # directories too.  This is the job_*/job_400GeV_* trap that job_files()
        # was written for, surviving one line below a correct call to it.
        # Caught by the per-job sigma_gen spread warning (84%), not by reading.
        for fn in job_files(f"{BASE}/genie", at_energy("nucclo_job"),
                            "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = "GENIE LO (GRV98LO)"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "GENIE LO (GRV98LO, charm)"
    elif gen in ("genie_nnpdf", "genie_nnpdf_charm"):
        # THE SAME GENIE LO PATH AS genie_lo, WITH ONE THING CHANGED: the
        # parton densities.  G18_02a builds Bodek-Yang on GRV98 LO, so the
        # GENIE row differs from every other LO generator here in TWO ways at
        # once -- the PDF, and everything else GENIE does differently
        # (Bodek-Yang low-Q2 modelling, AGKY hadronisation, its own
        # structure-function treatment).  Holding the PDF common to the rest
        # of the benchmark leaves the residual, which is the quantity of
        # interest (user request, 2026-08-27).
        #
        # genie/config/nnpdf swaps Uncorr-PDF-Set to LHAPDF6/NNPDF40_nnlo_0
        # inside BYPDF; genie/run_genie_nnpdf_nu.sh composes it with the
        # nu_cc_lo overlay that supplies the CCDISCHARM list.  Measured effect
        # on the total CC DIS cross-section at 1 TeV: 4.464 -> 4.723 pb,
        # i.e. +5.8%.
        #
        # The muon side has had this since the energy scan (nnpdf_job_*,
        # analyze.py's genie_nnpdf branch); this is the CC counterpart that
        # CONVENTIONS.md rule 2b asks for.
        files = job_files(f"{BASE}/genie", at_energy("nnpdfnu_job"),
                          "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        # job_files(), NOT a raw glob -- see the note in the genie_lo branch,
        # where a raw glob pooled two energies' cross-sections.
        for fn in job_files(f"{BASE}/genie", at_energy("nnpdfnu_job"),
                            "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = "GENIE LO (NNPDF4.0)"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "GENIE LO (NNPDF4.0, charm)"
    elif gen in ("v2_genie", "v2_genie_charm", "v2_genie_lo",
                 "v2_genie_lo_charm", "v2_genie_nnpdf", "v2_genie_nnpdf_charm"):
        # >>> PAPER PLOTS V2: the three earlier GENIE neutrino rows on ONE nucleon
        # (user, 2026-09-14, for pp03/pp04 ; mirror of analyze.py's
        # v2_genie, rule 2b). <<<  genie/production/genie_job.sh nu_hedis / nu_grv /
        # nu_nnpdf -- the earlier production's tunes, lists and overlays with the target as the only
        # change.  As in the earlier production, "genie" is HEDIS/BGR18 and "genie_lo" the classic
        # G18_02a.  Unweighted, full generation phase space, sigma_gen from
        # this nucleon's own spline.
        _sel, _t = sample_layout.require()
        _k = gen[3:].removesuffix("_charm")
        files = sample_layout.genie_files("nu", _k, _t, ENERGY, "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in sample_layout.genie_files("nu", _k, _t, ENERGY,
                                       "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"] * 1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        if not files:
            sys.exit(f"no complete GENIE {_k} jobs for {_t} at {ENERGY:g} GeV")
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = {"genie": "GENIE HEDIS (BGR18 NLO)",
                 "genie_lo": "GENIE LO (GRV98LO)",
                 "genie_nnpdf": "GENIE LO (NNPDF4.0)"}[_k]
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = label[:-1] + ", charm)"
    elif re.fullmatch(r"nugen_genie(_lo)?_(p|W)", gen):
        # >>> THE NEUTRINO-GENERATOR APPENDIX (user, 2026-10-01). <<<
        # genie/run_genie_nugen.sh: G18_02a (nugen_genie_lo_*, FASER's tune)
        # or HEDIS (nugen_genie_*), on a free proton or on the tungsten
        # NUCLEUS with GENIE's own nuclear model -- beside NuWro and GiBUU.
        # The beam record of a nuclear event is the struck nucleon
        # (genie/gtohepmc3.cc), so the invariants are per nucleon, and the
        # driver writes sigma_gen PER NUCLEON (checked here, not assumed).
        _cfg = "nu_grv" if "_lo_" in gen else "nu_hedis"
        _t = gen.rsplit("_", 1)[1]
        _base = at_energy(f"nugen_{_cfg}_{_t}_job")
        files, sig_pb, job_counts, n_gen_tot = [], [], [], 0
        for fn in job_files(f"{BASE}/genie", _base, "events.hepmc"):
            d_ = os.path.dirname(fn)
            if not os.path.exists(f"{d_}/NUGEN_OK"):
                continue                      # short job, kept for inspection
            with open(f"{d_}/events_xsec.json") as f:
                d = json.load(f)
            if d.get("per_nucleon_of") != (184 if _t == "W" else 1):
                sys.exit(f"{d_}/events_xsec.json is not per nucleon")
            files.append(fn)
            sig_pb.append(d["sigma_gen_mb"] * 1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        if not files:
            sys.exit(f"no complete GENIE {_cfg} jobs on {_t} at {ENERGY:g} GeV")
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = ("GENIE G18_02a (GRV98)" if _cfg == "nu_grv"
                 else "GENIE HEDIS (BGR18 NLO)") + (
                     " (free proton)" if _t == "p" else " (tungsten)")
    elif gen in ("genie", "genie_charm"):
        # unweighted, same per-job layout and normalisation as pythia: the
        # generation region (all Q2, y at 1 TeV) is WIDER than fiducial and
        # sigma_gen comes from the gmkspl spline summed at 1 TeV
        # (genie/spline_to_json.py).  HEDIS with BGR18 structure functions is
        # NLO in QCD, so the natural analytic reference is YADISM NLO.
        # genie_charm: outgoing hard charm tag, from the status-23 record the
        # converter writes out of GENIE's interaction summary (FinalQuarkPdg).
        files = job_files(f"{BASE}/genie", at_energy("nucc_job"), "events.hepmc")
        sig_pb, job_counts, n_gen_tot = [], [], 0
        for fn in job_files(f"{BASE}/genie", at_energy("nucc_job"), "events_xsec.json"):
            with open(fn) as f:
                d = json.load(f)
            sig_pb.append(d["sigma_gen_mb"]*1e9)
            job_counts.append(d["n_accepted"])
            n_gen_tot += d["n_accepted"]
        sigma_sample_pb = pooled_sigma_pb(sig_pb, job_counts, f" [{gen}]")
        sigma_err_pb = None
        unweighted = True
        label = "GENIE HEDIS (BGR18 NLO)"
        if gen == "genie_charm":
            charm_tag, tag_stats = True, {}
            label = "GENIE HEDIS (BGR18 NLO, charm)"
    elif gen in ("sherpa_lo", "sherpa_lo_charm"):
        files = job_files(f"{SHERPA_RUNS}/{at_energy('NuDIS_LO')}", "job", "evtfull")
        sigma_sample_pb = sherpa_nu_lo_sigma_pb()
        sigma_err_pb = SHERPA_NU_LO_SIGMA_ERR_PB
        n_gen_tot = None
        label = "Sherpa 3.0.5 LO"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "Sherpa 3.0.5 LO (charm)"
    elif gen == "sherpa_lo_me":
        # same process/cuts as NuDIS_LO, shower + hadronisation + remnants off
        files = job_files(f"{SHERPA_RUNS}/{at_energy('NuDIS_LO_ME')}", "job", "evtfull")
        sigma_sample_pb = SHERPA_NU_LO_ME_SIGMA_PB
        sigma_err_pb = SHERPA_NU_LO_ME_SIGMA_ERR_PB
        n_gen_tot = None
        label = "Sherpa 3.0.5 LO (ME)"
    elif gen in ("v2_powheg_nu", "v2_powheg_nu_charm"):
        # >>> PAPER PLOTS V2 (user, 2026-09-13): POWHEG-V2 + Pythia 8 in the
        # final region, one nucleon per run (PRODUCTION.md).  Mirror of
        # analyze.py's v2_powheg (rule 2b). <<<  Born q2cut 2.25, no y cut,
        # so WIDER than q4w3: integrator x fiducial weight fraction, the
        # integrator from this nucleon's own pwg-stat.dat (its error too).
        # The closure denominator counts the added -b<j> batches, which
        # lhe_events_offered() already globs.  Neutron runs have a genuine
        # 2112 beam: POWHEG-V2 (ih2 2) and Pythia swap the proton set.
        _sel, _t = sample_layout.require()
        files = job_files(f"{BASE}/powheg",
                          at_energy(sample_layout.powheg_job_base("nu", _t)),
                          "events.hepmc")
        # THE CHARM FIX IS REQUIRED, not assumed (2026-09-14): a job showered
        # with Pythia's default LesHouches:matchInOut = on lost ~2% of its
        # events, nearly all charm (powheg/powheg_nu_v2.cmnd).  Pythia prints
        # every changed setting in its init table, so the log says which.
        for _fn in files:
            _log = os.path.join(os.path.dirname(_fn), "shower.log")
            if not re.search(r"LesHouches:matchInOut\s*\|\s*off",
                             open(_log, errors="replace").read()
                             if os.path.exists(_log) else ""):
                sys.exit(f"{_log}: not showered with LesHouches:matchInOut "
                         f"= off -- run powheg/production/reshower_nu_matchinout.sh")
        _f = f"{sample_layout.powheg_rundir('nu', _t, ENERGY)}/pwg-stat.dat"
        sigma_sample_pb = sigma_err_pb = None
        for _line in (open(_f, errors="replace") if os.path.exists(_f)
                      else []):
            _m = re.search(r"total \(btilde\+remnants\) cross section in pb"
                           r"\s+([0-9.eEdD+-]+)\s+\+-\s+([0-9.eEdD+-]+)",
                           _line)
            if _m:
                sigma_sample_pb, sigma_err_pb = (
                    float(g.replace("D", "E").replace("d", "e"))
                    for g in _m.groups())
        if sigma_sample_pb is None:
            sys.exit(f"no POWHEG-V2 total cross section in {_f}")
        print(f"  [powheg-v2 {_t}] integrator {sigma_sample_pb:.5f} "
              f"+- {sigma_err_pb:.5f} pb (from {_f})")
        n_gen_tot = None
        label = "POWHEG-V2 NLO + Pythia8"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("v2_powheg_nu_mc", "v2_powheg_nu_mc_charm"):
        # >>> PAPER PLOTS V2: POWHEG-V2mc, MASSIVE charm (qmass 1.51, numflav
        # 4), one nucleon per run (powheg/production/run_v2_mc.sh, user 2026-09-14).
        # THE SAMPLE IS CHARM PRODUCTION ONLY -- every LHE event carries the
        # massive charm (see powheg_nu_mc below) -- so its sigma sits beside
        # charm rows, never inclusive; `_charmfinal` asks whether a charm
        # HADRON survived, the benchmark definition.  Otherwise the massless
        # branch below it, line for line (the copy is the text that follows):
        # POWHEG-V2 + Pythia 8 in the
        # final region, one nucleon per run (PRODUCTION.md).  Mirror of
        # analyze.py's v2_powheg (rule 2b). <<<  Born q2cut 2.25, no y cut,
        # so WIDER than q4w3: integrator x fiducial weight fraction, the
        # integrator from this nucleon's own pwg-stat.dat (its error too).
        # The closure denominator counts the added -b<j> batches, which
        # lhe_events_offered() already globs.  Neutron runs have a genuine
        # 2112 beam: POWHEG-V2 (ih2 2) and Pythia swap the proton set.
        _sel, _t = sample_layout.require()
        files = job_files(f"{BASE}/powheg",
                          at_energy(sample_layout.powheg_job_base("nu", _t, mc=True)),
                          "events.hepmc")
        # THE CHARM FIX IS REQUIRED, not assumed (2026-09-14): a job showered
        # with Pythia's default LesHouches:matchInOut = on lost ~2% of its
        # events, nearly all charm (powheg/powheg_nu_v2.cmnd).  Pythia prints
        # every changed setting in its init table, so the log says which.
        for _fn in files:
            _log = os.path.join(os.path.dirname(_fn), "shower.log")
            if not re.search(r"LesHouches:matchInOut\s*\|\s*off",
                             open(_log, errors="replace").read()
                             if os.path.exists(_log) else ""):
                sys.exit(f"{_log}: not showered with LesHouches:matchInOut "
                         f"= off -- run powheg/production/reshower_nu_matchinout.sh")
        _f = f"{sample_layout.powheg_rundir('nu', _t, ENERGY, mc=True)}/pwg-stat.dat"
        sigma_sample_pb = sigma_err_pb = None
        for _line in (open(_f, errors="replace") if os.path.exists(_f)
                      else []):
            _m = re.search(r"total \(btilde\+remnants\) cross section in pb"
                           r"\s+([0-9.eEdD+-]+)\s+\+-\s+([0-9.eEdD+-]+)",
                           _line)
            if _m:
                sigma_sample_pb, sigma_err_pb = (
                    float(g.replace("D", "E").replace("d", "e"))
                    for g in _m.groups())
        if sigma_sample_pb is None:
            sys.exit(f"no POWHEG-V2mc total cross section in {_f}")
        print(f"  [powheg-v2mc {_t}] integrator {sigma_sample_pb:.5f} "
              f"+- {sigma_err_pb:.5f} pb (from {_f})")
        n_gen_tot = None
        label = "POWHEG-V2mc NLO + Pythia8 (massive charm)"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif (gen.removesuffix("_charm").startswith("v2_powheg_nu_")
          and gen.removesuffix("_charm")[len("v2_powheg_nu_"):] in sample_layout.ARMS
          and gen.removesuffix("_charm") != "v2_powheg_nu_mc"):
        # >>> PAPER PLOTS 9-11: THE SHOWER ARMS, neutrino side (rule 2b).
        # Mirror of analyze.py's v2_powheg_<arm>. <<<  The SAME 1 TeV
        # POWHEG-V2 Les Houches batches as v2_powheg_nu, on the arm subset of
        # its jobs: `sub` the Pythia baseline, `qedfsr/qedfull` the QED
        # re-showers (no fsrisr: a neutrino beam has no ISR photon), `hw`
        # Herwig 7.  Integrator x fiducial weight fraction over the subset
        # (each batch is an independent estimate of the integral).
        _sel, _t = sample_layout.require()
        _arm = gen.removesuffix("_charm")[len("v2_powheg_nu_"):]
        files = sample_layout.arm_files("nu", _arm, _t)
        if _arm != "hw":
            for _fn in files:
                _log = os.path.join(os.path.dirname(_fn), "shower.log")
                if not re.search(r"LesHouches:matchInOut\s*\|\s*off",
                                 open(_log, errors="replace").read()
                                 if os.path.exists(_log) else ""):
                    sys.exit(f"{_log}: not showered with LesHouches:matchInOut = off")
        _f = f"{sample_layout.powheg_rundir('nu', _t, ENERGY)}/pwg-stat.dat"
        sigma_sample_pb = sigma_err_pb = None
        for _line in (open(_f, errors="replace") if os.path.exists(_f) else []):
            _m = re.search(r"total \(btilde\+remnants\) cross section in pb"
                           r"\s+([0-9.eEdD+-]+)\s+\+-\s+([0-9.eEdD+-]+)", _line)
            if _m:
                sigma_sample_pb, sigma_err_pb = (
                    float(g.replace("D", "E").replace("d", "e")) for g in _m.groups())
        if sigma_sample_pb is None:
            sys.exit(f"no POWHEG-V2 total cross section in {_f}")
        n_gen_tot = None
        label = {"sub": "POWHEG-V2 NLO + Pythia8",
                 "qedfsr": "POWHEG-V2 NLO + Pythia8, QED FSR (lepton line)",
                 "qedfull": "POWHEG-V2 NLO + Pythia8, QED FSR + quark line",
                 "hw": "POWHEG-V2 NLO + Herwig7"}[_arm]
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("v2_sherpa_nlo_dire", "v2_sherpa_nlo_dire_charm"):
        # >>> PAPER PLOT 9: Sherpa CC MC@NLO with the DIRE shower, one
        # nucleon, 1 TeV.  Mirror of analyze.py's v2_sherpa_dire. <<<
        _sel, _t = sample_layout.require()
        _rd = sample_layout.sherpa_rundir("nu", _t, ENERGY) + "_Dire"
        files = job_files(_rd, "job", "evtfull")
        _s, _e = sherpa_arm_sigma_pb(os.path.relpath(_rd, SHERPA_RUNS), True,
                                     with_err=True)
        sigma_sample_pb = SHERPA_NU_SPIN_FACTOR * _s
        sigma_err_pb = SHERPA_NU_SPIN_FACTOR * _e
        n_gen_tot = None
        label = "Sherpa 3.0.5 MC@NLO, Dire"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("v2_sherpa_nlo", "v2_sherpa_nlo_charm"):
        # >>> PAPER PLOTS V2: Sherpa CC MC@NLO in the final region, one
        # nucleon.  Mirror of analyze.py's v2_sherpa (rule 2b). <<<
        # NORMALISED TO THE INTEGRATOR (BVI + RS from this run's own
        # integ.log) x the helicity factor 2, as the muon side is -- NOT to a
        # hand-copied event-level table like the earlier production's SHERPA_NU_NLO, which is the
        # 1 TeV-literal shape this benchmark has hit six times.  The
        # event-level number is still compared, by the closure gate below.
        # q4w3 is a subset of the generation region (INEL at the W = 3 edge
        # of the Q2 floor), so the weighted branch takes the weight fraction.
        _sel, _t = sample_layout.require()
        _rd = sample_layout.sherpa_rundir("nu", _t, ENERGY)
        files = job_files(_rd, "job", "evtfull")
        _s, _e = sherpa_arm_sigma_pb(os.path.relpath(_rd, SHERPA_RUNS), True,
                                     with_err=True)
        sigma_sample_pb = SHERPA_NU_SPIN_FACTOR * _s
        sigma_err_pb = SHERPA_NU_SPIN_FACTOR * _e
        n_gen_tot = None
        label = "Sherpa 3.0.5 MC@NLO"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("v2_herwig_nlo", "v2_herwig_nlo_neg",
                 "v2_herwig_nlo_charm", "v2_herwig_nlo_neg_charm"):
        # >>> PAPER PLOTS V2: Herwig 7 POWHEG NLO, one nucleon, one half. <<<
        # Mirror of herwig7/make_histos_nlo.py's modes (rule 2b).  Unit
        # weights, generated at MinQ2 3.5 / MinW2 5 -- looser than q4w3, so
        # the unweighted branch's sigma_gen x N_fid / N_gen is the right
        # recipe.  Polarised beam in the card, so no spin factor.
        _sel, _t = sample_layout.require()
        _neg = gen.removesuffix("_charm").endswith("_neg")
        _jb = at_energy(sample_layout.herwig_job_base("nu", _t, _neg))
        files = job_files(f"{BASE}/herwig7", _jb, "events.hepmc")
        jobs = hw_read_jobs(
            f"{BASE}/herwig7/{_jb}_*/"
            f"{sample_layout.herwig_runname('nu', _t, ENERGY, _neg)}-S*.out")
        if len(jobs) != len(files):
            sys.exit(f"MISMATCH: {len(files)} events.hepmc but {len(jobs)} "
                     f".out files -- refusing to normalise a mixed sample.")
        sig, err, n_gen_tot, survival = hw_combine(jobs)
        print(f"  herwig nu NLO {_t}{' (negative)' if _neg else ''}: "
              f"survival {survival:.4f}")
        sigma_sample_pb = HERWIG_NU_SPIN_FACTOR * sig
        sigma_err_pb = HERWIG_NU_SPIN_FACTOR * err
        unweighted = True
        label = ("Herwig 7.3.0 NLO (POWHEG" + (", negative" if _neg else "")
                 + ")")
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen in ("sherpa_nlo", "sherpa_nlo_charm"):
        # Sherpa 3.0.5 CC MC@NLO (Runs/NuDIS_NLO); generation cuts == fiducial
        files = job_files(f"{SHERPA_RUNS}/{at_energy('NuDIS_NLO_ckm3')}",
                          "job", "evtfull")
        if SHERPA_NU_NLO_SIGMA_PB is None:
            _sherpa_nu_nlo(0)                 # raises: no entry for this energy
        sigma_sample_pb = SHERPA_NU_NLO_SIGMA_PB
        sigma_err_pb = SHERPA_NU_NLO_SIGMA_ERR_PB
        n_gen_tot = None
        label = "Sherpa 3.0.5 MC@NLO"
        if gen.endswith("_charm"):
            # Full CKM as of 2026-08-25, so this row is finally like-for-like
            # with the LO card: s -> c carries |V_cs|^2 and d -> c off the
            # valence d is present, rather than being dropped.
            charm_tag, tag_stats = True, {}
            label = "Sherpa 3.0.5 MC@NLO (charm)"
    elif gen.replace("_charm", "") in SHERPA_NU_SHOWER_ARMS:
        # The CC half of the NLO shower study (CONVENTIONS.md rule 2b: whatever is
        # done on one current is owed to the other).  Same card, same matrix
        # element, SHOWER_GENERATOR overridden on the command line; the arm
        # reads its OWN integ.log because the MC@NLO subtraction is built from
        # the shower's splitting kernels.
        base = gen.replace("_charm", "")
        rundir, is_nlo, lab = SHERPA_NU_SHOWER_ARMS[base]
        files = job_files(f"{SHERPA_RUNS}/{rundir}", "job", "evtfull")
        # THE SPIN FACTOR APPLIES HERE TOO.  Sherpa's neutrino beam carries a
        # factor-2 helicity-average error, and every other CC Sherpa
        # normalisation on this page is SHERPA_NU_SPIN_FACTOR times the raw
        # integrator value.  sherpa_arm_sigma_pb reads the raw integ.log, so
        # the factor has to be applied here -- omitting it put Dire/CSS at
        # EXACTLY 0.5000, which is the only reason it was noticed.
        _raw, _rawerr = sherpa_arm_sigma_pb(rundir, is_nlo, with_err=True)
        sigma_sample_pb = SHERPA_NU_SPIN_FACTOR * _raw
        sigma_err_pb = SHERPA_NU_SPIN_FACTOR * _rawerr
        n_gen_tot = None
        label = lab
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label += " (charm)"
    elif gen.startswith("gibuu"):
        # >>> GiBUU: A TRANSPORT MODEL, RUN WELL ABOVE ITS DOCUMENTED RANGE.
        # <<<  Its neutrino module is documented for 1-50 GeV
        # (arXiv:1205.1061) and this benchmark runs at 400 GeV - 4 TeV, so
        # every GiBUU number here is an extrapolation and every figure says
        # so.  It does run, and the energy dependence is right; see
        # gibuu/README.md, which also records the two patches it needed.
        #
        # WEIGHTED, WITH THE GENERATION REGION WIDER THAN FIDUCIAL.  The
        # converter writes each event's own cross-section contribution in pb,
        # so the weights of the whole file sum to sigma_tot and the weights of
        # the selected events sum to sigma_fid.  Keys as for NuWro: `gibuu`
        # free proton, `gibuu_n` free neutron, `gibuu_W` tungsten with the
        # hadrons transported out of the nucleus.
        tgt = gen[len("gibuu"):].lstrip("_") or "p"
        base = at_energy(f"gibuu_job_nu_{tgt}")
        files = job_files(f"{BASE}/gibuu", base, "events.hepmc")
        sigma_sample_pb, n_gen_tot = nuwro_sigma_and_n(files, "gibuu")
        sigma_err_pb = None
        label = {"p": "GiBUU (free proton)", "n": "GiBUU (free neutron)",
                 "W": "GiBUU (tungsten, transport)"}.get(tgt, f"GiBUU ({tgt})")
    elif gen.startswith("nuwro"):
        # >>> NuWro: A DEDICATED NEUTRINO GENERATOR, CONVERTED TO HepMC3. <<<
        # nuwro/run_nuwro.sh generates and nuwro/nuwro2hepmc converts, so the
        # events reach this function through exactly the same parser as every
        # other sample (see nuwro/README.md for why a converter and not a
        # third reader).  The key names the target:
        #     nuwro      free proton      -- like for like with the rest of
        #                                   this benchmark
        #     nuwro_n    free neutron
        #     nuwro_W    tungsten, Fermi gas and intranuclear cascade ON
        # and the difference between the first and the last is the nuclear
        # physics this benchmark otherwise omits, which is what NuWro is in
        # the comparison for.
        #
        # EQUALLY WEIGHTED, SO IT IS AN "unweighted" SAMPLE HERE.  The
        # converter writes sigma/N on every event and records both numbers in
        # events_xsec.json; the generation region is the whole DIS-CC phase
        # space, which is WIDER than the fiducial cuts, so sigma_fid is
        # sigma_tot x N_selected/N_generated -- the same recipe as GENIE.
        tgt = gen[len("nuwro"):].lstrip("_") or "p"
        base = at_energy(f"nuwro_job_{'nubar' if 'bar' in tgt else 'nu'}"
                         f"_{tgt.replace('bar', '') or 'p'}")
        files = job_files(f"{BASE}/nuwro", base, "events.hepmc")
        sigma_sample_pb, n_gen_tot = nuwro_sigma_and_n(files)
        sigma_err_pb = None
        unweighted = True
        label = {"p": "NuWro (free proton)", "n": "NuWro (free neutron)",
                 "W": "NuWro (tungsten, FSI)"}.get(tgt, f"NuWro ({tgt})")
    elif gen in ("powheg_nu", "powheg_nu_charm"):
        # POWHEG-V2 nu-DIS NLO + Pythia 8.311 shower (powheg/main_powheg with
        # powheg_nu1TeV.cmnd).  Weighted +-sigma_abs events; generation region
        # (Born Q2 > 2.25, full y) is WIDER than the fiducial cuts.
        files = job_files(f"{BASE}/powheg", at_energy("nu_job"), "events.hepmc")
        sigma_sample_pb = powheg_v2_nu_sigma_pb()
        sigma_err_pb = POWHEG_V2_NU_SIGMA_ERR_PB
        n_gen_tot = None
        label = "POWHEG-V2 NLO + Pythia8"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "POWHEG-V2 NLO + Pythia8 (charm)"
    elif gen in ("powheg_nu_hw", "powheg_nu_hw_charm"):
        # >>> THE SAME POWHEG-V2 LHE EVENTS, SHOWERED BY HERWIG 7. <<<
        # powheg/run_powheg_herwig.sh --current nu, card herwig7/LHE-nu.in.
        # The charged-current half of the hadronisation study (user,
        # 2026-09-08); CONVENTIONS.md rule 2b requires it beside the muon
        # `powheg_hw` and it carries the same argument -- one matrix element,
        # one matching, one parton set, two showers and two hadronisation
        # models, so the difference between this and `powheg_nu` is the
        # fragmentation and nothing else.
        #
        # NO BOOST WAS NEEDED HERE.  POWHEG-V2 writes this side already in
        # the lepton-proton c.m. frame; POWHEG-RES writes the muon side
        # fixed-target, with the proton at rest and therefore with no
        # momentum, which ThePEG cannot build parton bins from.  See
        # tools/lhe_boost_cm.py, which decides between the two and says which
        # it did.
        files = job_files(f"{BASE}/powheg", at_energy("hw_nu_job"),
                          "events.hepmc")
        sigma_sample_pb = powheg_v2_nu_sigma_pb()
        sigma_err_pb = POWHEG_V2_NU_SIGMA_ERR_PB
        n_gen_tot = None
        label = "POWHEG-V2 NLO + Herwig7"
        if gen.endswith("_charm"):
            charm_tag, tag_stats = True, {}
            label = "POWHEG-V2 NLO + Herwig7 (charm)"
    elif gen in ("powheg_nu_mc", "powheg_nu_mc_charm"):
        # POWHEG-V2mc: POWHEG-V2 with MASSIVE charm (qmass = 1.51, numflav 4),
        # charged current only.  Card cards/POWHEG-V2mc/nu<TAG>-prod.
        #
        # THIS SAMPLE IS CHARM PRODUCTION, NOT AN INCLUSIVE CROSS-SECTION, and
        # that is the single most important thing to know about it.  With a
        # non-zero qmass the massive final-state quark IS the process: measured
        # on the LHE, 20000 of 20000 events carry an outgoing charm quark,
        # against 20.2% for the massless POWHEG-V2 sample.  So its
        # cross-section belongs beside the CHARM rows and beside YADISM FONLL
        # charm, never beside an inclusive number -- comparing it with one
        # would look like a factor-five discrepancy and be a category error.
        #
        # No `_charm` variant is registered for the same reason: the whole
        # sample is the charm sample, so a charm tag on top of it would be a
        # second, redundant cut and the two keys would differ only by the
        # tagging efficiency.
        files = job_files(f"{BASE}/powheg", at_energy("numc_job"),
                          "events.hepmc")
        sigma_sample_pb = powheg_v2mc_sigma_pb()
        sigma_err_pb = None
        n_gen_tot = None
        label = "POWHEG-V2mc NLO + Pythia8 (massive charm)"
        if gen.endswith("_charm"):
            # THE HADRON-LEVEL TAG ON TOP OF A PARTON-LEVEL CHARM SAMPLE, and
            # it is not redundant.  The sample is charm by construction at
            # PARTON level; `_charmfinal` asks whether a charm HADRON survived
            # the shower and hadronisation, which is the benchmark's own charm
            # definition and the one `powheg_nu_charmfinal` uses.  The two
            # differ by about 7% in this project's measurement, so comparing
            # the massive sample's parton-level rate against the massless
            # sample's hadron-level tag mixes definitions -- which is exactly
            # what made the massive sample look FURTHER from FONLL than the
            # massless one.  This key is what makes that comparison honest.
            charm_tag, tag_stats = True, {}
            label = "POWHEG-V2mc NLO + Pythia8 (massive charm, charm)"
    elif gen in ("powheg_res_nu", "powheg_res_nu_charm"):
        # POWHEG-RES (powheg-dis-main/DIS_v) on the neutrino CC -- the second
        # independent NLO matching on this side, the cross-check partner of
        # the POWHEG-V2 entry above.  Card cards/POWHEG-RES/nu1TeV-CC,
        # showered through powheg/main_powheg with powheg_nu1TeV.cmnd.
        # Generation region (Born Q2 > 2.25 via Qmin 1.5, full y) is WIDER
        # than the fiducial cuts, so the same recipe as powheg_nu applies.
        if gen.endswith("_charm"):
            sys.exit(
                "powheg_res_nu has NO CHARM ROW, and this is not an oversight."
                "\n  DIS_v hardcodes CKM_diag = .true.: the CKM matrix never "
                "reaches the\n  matrix elements and d -> c is absent, so the "
                "charm composition is wrong\n  even though CKM unitarity keeps "
                "the inclusive rate right to 0.22%.\n  The identical defect in "
                "Sherpa cost 3.4 points on the charm fraction.\n  Use "
                "powheg_nu_charm (POWHEG-V2), which has the full CKM.")
        files = job_files(f"{BASE}/powheg", at_energy("resnu_job"),
                          "events.hepmc")
        sigma_sample_pb = powheg_res_nu_sigma_pb()
        sigma_err_pb = POWHEG_RES_NU_SIGMA_ERR_PB
        n_gen_tot = None
        label = "POWHEG-RES NLO + Pythia8 (diagonal CKM, inclusive only)"
    elif qed_arm_nu(gen):
        # A QED arm on the CC: the same LHE files as powheg_res_nu, re-showered
        # with photon radiation on.  Normalisation, generation region and
        # closure are powheg_res_nu's, unchanged.
        arm = qed_arm_nu(gen)
        if gen.endswith("_charm"):
            sys.exit(
                "the powheg_res_nu QED arms have NO CHARM ROW, for the same "
                "reason\n  powheg_res_nu itself does not: DIS_v hardcodes "
                "CKM_diag = .true., so\n  d -> c is absent from the matrix "
                "elements.  Re-showering cannot add a\n  channel the matrix "
                "element never had.  Use the inclusive arm.")
        files = job_files(f"{BASE}/powheg", at_energy(f"qed{arm}_resnu_job"),
                          "events.hepmc")
        sigma_sample_pb = powheg_res_nu_sigma_pb()
        sigma_err_pb = POWHEG_RES_NU_SIGMA_ERR_PB
        n_gen_tot = None
        label = (f"POWHEG-RES NLO + Pythia8, {QED_ARMS[arm]} "
                 f"(diagonal CKM, inclusive only)")
    else:
        sys.exit(f"unknown generator {gen}")
    if not files:
        sys.exit(f"no input files found for {gen}")

    run_stats = {}
    subtract = subtracts_shower_charm(gen)
    if subtract and tag_stats is None:
        tag_stats = {}
    obs, w, n_events, sum_w_all = analyze(files, charm_tag, tag_stats,
                                          run_stats=run_stats,
                                          charm_mode=charm_mode,
                                          subtract_shower_charm=subtract)
    print(f"{gen}: parsed {n_events} events, {len(w)} in fiducial region")
    # >>> NuWro's ONE REAL INTEGRITY CHECK. <<<  Its per-event weight is
    # sigma/N by construction, so the delivery closure gate is vacuous for it
    # (see nuwro_sigma_and_n).  A truncated or partially converted HepMC would
    # therefore pass every other test and simply report a smaller
    # cross-section.  The generated count recorded at conversion time is what
    # catches it, and a mismatch is an error rather than a warning.
    if (gen.startswith(("nuwro", "gibuu")) and n_gen_tot
            and n_events != n_gen_tot):
        sys.exit(f"{gen}: parsed {n_events} events but the conversion "
                 f"recorded {n_gen_tot} -- the HepMC file is short, and the "
                 f"fiducial fraction would be computed against the wrong "
                 f"denominator")
    # >>> NuWro's WEIGHT MUST NOT BE APPLIED TWICE (found 2026-10-01). <<<
    # The branch above treats NuWro as unweighted, normalising by
    # sigma/N x N_bin -- but the parser returns the per-event weight the
    # converter wrote, which IS sigma/N.  sigma_fid (len(w)) was right, every
    # histogram was sigma/N = 2e-5 too small, and nothing noticed because the
    # only figure reading them normalised each curve to unity.  So the weights
    # are checked to be the uniform sigma/N and set to one.
    if unweighted and gen.startswith("nuwro") and len(w):
        _wn = sigma_sample_pb / n_gen_tot
        if not np.allclose(w, _wn, rtol=1e-4):
            sys.exit(f"{gen}: event weights are not the uniform sigma/N = "
                     f"{_wn:g} pb the converter writes")
        w = np.ones_like(w)
    if tag_stats is not None:
        cap = (tag_stats["n_d_tagged"] / tag_stats["n_d_all"]
               if tag_stats["n_d_all"] else float("nan"))
        print(f"{gen}: {tag_stats['n_tagged']} charm-tagged events "
              f"({100*tag_stats['n_tagged']/n_events:.2f}%); tagged events "
              f"carry {100*cap:.1f}% of all D mesons")
        print(f"{gen}: tag comparison on the same sample -- "
              f"hard-out {tag_stats.get('n_tagged_out', 0)}, "
              f"hard-any {tag_stats.get('n_tagged_any', 0)}, "
              f"final-state {tag_stats.get('n_tagged_final', 0)} "
              f"of {n_events} parsed; "
              f"{tag_stats.get('n_charm_with_bottom', 0)} vetoed as "
              f"charm-from-bottom")

    if charm_mode == "any":
        gen = gen + "any"
        label = label.replace("(charm)", "(charm, hard-process any)")
    elif charm_mode == "final":
        # the label is left alone: "final" is the benchmark definition, so
        # every charm row carries it and decorating one label but not another
        # (some already read "(POWHEG, charm)" or "(charm, unit CKM)") would
        # make the report table look like a mix of definitions.  The mode is
        # recorded in the JSON as charm_tag_mode.
        gen = gen + "final"
    # THE DISPATCH KEY IS THE RUN, NOT THE TAG (pipeline review 2026-09-19,
    # finding C1): with the charm suffix on, v2_powheg_nu_charmfinal and
    # v2_powheg_nu_mc_charmfinal missed the POWHEG branch below and fell into
    # the generic one -- same sigma formula, but "closure NOT checked", so
    # POWHEG-V2mc had never been through the gate.  analyze.py strips the tag
    # the same way (pbase).
    _gbase = re.sub(r"_charm(any|final)?$", "", gen)
    out = {"generator": gen, "label": label, "n_parsed": n_events,
           # see analyze.py: one field, one name, with both legacy keys kept
           "tag_mode": charm_mode,
           "charm_tag_mode": charm_mode,
           "hard_tag_mode": charm_mode,
           "n_fiducial": len(w)}
    out["inputs"] = input_manifest(files, sig_pb if sig_pb else None,
                                   job_counts if job_counts else None,
                                   f" [{gen}]")
    if tag_stats is not None:
        out["tag_stats"] = tag_stats
    if unweighted:
        # unweighted: dsigma = (N_bin / N_generated) * sigma_generated / width.
        # Pythia's and GENIE's generation regions (all y, GENIE also all Q2)
        # are WIDER than fiducial; Herwig's generation cuts already equal the
        # fiducial ones.  A charm tag needs no special case: it removes events
        # from N_bin and from len(w) but not from N_generated, so what comes
        # out is the channel's absolute contribution to sigma_fid.
        norm = sigma_sample_pb / n_gen_tot
        out["sigma_fid_pb"] = norm * len(w)
        # MC statistical error, plus the integrator error where it is known
        # (Herwig quotes one; Pythia's is negligible next to the count error)
        stat = norm * math.sqrt(len(w))
        if sigma_err_pb is not None:
            stat = math.hypot(stat, sigma_err_pb * len(w) / n_gen_tot)
        out["sigma_fid_err_pb"] = stat
    elif gen.startswith("gibuu"):
        # weighted, generation region wider than fiducial, and the weights
        # already carry the absolute normalisation -- so sigma_fid is simply
        # their sum over the selected events.  The closure between that and
        # the sidecar's sigma_tot is exact by construction (both are sums of
        # the same weights), which is why the integrity check that matters is
        # the event count, as for NuWro.
        # `norm` is 1: the weights ALREADY carry the absolute normalisation,
        # so the histogram code below multiplies by unity rather than by a
        # ratio.  It is set explicitly because every other branch sets it and
        # the histogram loop reads it unconditionally.
        norm = 1.0
        out["sigma_fid_pb"] = float(w.sum())
        out["sigma_fid_err_pb"] = float(np.sqrt((w * w).sum()))
        out["sigma_total_pb"] = sigma_sample_pb
        out["fiducial_weight_fraction"] = (float(w.sum()) / sum_w_all
                                           if sum_w_all else None)
        print(f"{gen}: fiducial weight fraction = "
              f"{out['fiducial_weight_fraction']:.4f}")
    elif (_gbase in ("powheg_nu", "powheg_nu_mc", "powheg_res_nu",
                     "powheg_nu_hw", "v2_powheg_nu", "v2_powheg_nu_mc")
          or qed_arm_nu(_gbase)
          or (_gbase.startswith("v2_powheg_nu_")
              and _gbase[len("v2_powheg_nu_"):] in sample_layout.ARMS)):
        # weighted, generation region WIDER than fiducial: sigma_fid is the
        # integrator sigma times the fiducial weight fraction (same recipe as
        # the muon POWHEG-RES branch of analyze.py)
        norm = sigma_sample_pb / sum_w_all
        out["sigma_fid_pb"] = norm * float(w.sum())
        out["sigma_fid_err_pb"] = norm * float(np.sqrt((w * w).sum()))
        out["mean_weight_pb"] = sum_w_all / n_events
        out["fiducial_weight_fraction"] = float(w.sum()) / sum_w_all
        print(f"{gen}: fiducial weight fraction = "
              f"{out['fiducial_weight_fraction']:.4f}, negative-weight events "
              f"= {int((w < 0).sum())}/{len(w)}")
        if _gbase.startswith("powheg_nu_hw") or _gbase == "v2_powheg_nu_hw":
            # THE HERWIG ARM CANNOT USE THE MEAN-WEIGHT CLOSURE: ThePEG
            # normalises the HepMC weights to the largest weight in the file,
            # so they carry no absolute scale and their mean is a pure
            # number.  Herwig's own end-of-run table quotes the cross-section
            # it delivered, which is the independent measurement the pinned
            # normalisation has to be checked against.  Identical to the
            # muon side (analyze.py), which is the point.
            dlv, n_w, n_a = herwig_lhe_delivered_pb(files)
            out["n_lhe_offered"] = n_a
            out["herwig_shower_survival"] = (n_w / n_a) if n_a else None
            closure_check(gen, dlv, sigma_sample_pb, w, out,
                          "Herwig's own delivered sigma over the Les Houches "
                          "events consumed")
        else:
            # an arm reads powheg_res_nu's LHE directory, by construction
            lhe_key = "powheg_res_nu" if qed_arm_nu(_gbase) else _gbase
            if (_gbase.startswith("v2_powheg_nu_")
                    and _gbase[len("v2_powheg_nu_"):] in sample_layout.ARMS):
                # a arm or its baseline: the subset's own offered events
                n_read = sample_layout.powheg_offered("nu", sample_layout.require()[1],
                                                 ENERGY, files)
            else:
                n_read = lhe_events_offered(
                    POWHEG_NU_LHE_DIR[lhe_key](beams.Beams("nu", ENERGY).tag))
            out["n_lhe_offered"] = n_read
            dlv, how = delivered_sigma_pb(run_stats, n_events, n_read)
            closure_check(gen, dlv, sigma_sample_pb, w, out, how)
    else:
        # weighted, generation cuts == fiducial cuts: shape from the weights,
        # normalisation pinned to the (factor-2-corrected) integrator sigma
        norm = sigma_sample_pb / sum_w_all
        # A SELECTION NARROWER THAN THE GENERATION CUTS breaks the pinning
        # exactly as a channel tag does; "generation cuts == fiducial cuts"
        # holds only for the INCLUSIVE selection these samples were generated
        # for.  Kept identical to analyze.py (CONVENTIONS.md rule 2b): the muon side
        # was fixed first and the neutrino side still reported Sherpa's Tier S
        # as 3.1906 pb, its inclusive value to the digit.
        # A RAISED Q2 FLOOR IS A SUBSET TOO -- see the note in analyze.py:
        # the selection's NAME does not change, so the name test alone pinned
        # a Q2 > 11 result to the Q2 > 4 integrator cross-section.
        subset = (charm_tag or SELECTION.name != selection_mod.DEFAULT
                  or SELECTION.q2_min != selection_mod.Q2_FLOOR_DEFAULT)
        if subset:
            # sum_w_all still counts every parsed event, so the absolute
            # contribution of the surviving subset is the pinned sigma times
            # the surviving weight fraction
            out["sigma_fid_pb"] = norm * float(w.sum())
            out["sigma_fid_err_pb"] = norm * float(np.sqrt((w * w).sum()))
            out["selected_weight_fraction"] = float(w.sum()) / sum_w_all
            if charm_tag:
                out["tagged_weight_fraction"] = float(w.sum()) / sum_w_all
        else:
            out["sigma_fid_pb"] = sigma_sample_pb
            out["sigma_fid_err_pb"] = sigma_err_pb
        out["mean_weight_pb"] = sum_w_all / n_events
        # the events carry the RAW weights, so compare against the raw
        # integrator and let the spin factor cancel out of the ratio
        # a channel tag does not change the run, so strip it for the lookup
        # (including the charm-mode suffix appended a few lines above)
        base = re.sub(r"_charm(any|final)?$", "", gen)
        integ_raw = SHERPA_NU_INTEG_RAW_PB.get(base)
        # PER-ENERGY RAW INTEGRATOR.  The table above holds 1 TeV literals, so
        # away from the anchor the gate compared 400 GeV events (0.6451 pb)
        # against the 1 TeV integrator (1.59692) and reported -59.6%.  That is
        # the SIXTH instance of the 1 TeV-literal trap in this scan, and the
        # closure gate caught it exactly as it caught the others.  The raw
        # value is the spin-factor-corrected one divided back out, so the
        # comparison stays raw-against-raw.
        if integ_raw is not None and ENERGY != beams.ANCHOR_ENERGY:
            if base == "sherpa_lo":
                integ_raw = sherpa_nu_lo_sigma_pb() / SHERPA_NU_SPIN_FACTOR
            elif base == "sherpa_nlo":
                # per-energy since the CC MC@NLO scan was run; SHERPA_NU_NLO
                # is keyed by beam tag, so this cannot silently fall back to
                # the anchor's number the way a bare literal would
                integ_raw = SHERPA_NU_NLO_INTEG_RAW_PB
            else:
                integ_raw = None      # no per-energy number for this run yet
        if base.startswith("v2_sherpa_nlo"):
            # (_dire included: its own integ.log, same helicity factor --
            # without this the Dire twins were never closure-checked;
            # pipeline review 2026-09-19, finding C1)
            # read from this run's own integ.log above; raw = before the
            # helicity factor, which is what the event weights carry
            integ_raw = sigma_sample_pb / SHERPA_NU_SPIN_FACTOR
        if integ_raw is not None:
            dlv, how = delivered_sigma_pb(run_stats, n_events)
            # MC@NLO delivers S and H events whose weights sum to slightly
            # less than the BVI+RS integral; ~0.6% here is expected and was
            # documented before this gate existed.
            note = ("expected for MC@NLO: the S/H event split does not "
                    "reproduce BVI+RS exactly"
                    if base in ("sherpa_nlo", "v2_sherpa_nlo",
                                "v2_sherpa_nlo_dire") else "")
            closure_check(gen, dlv, integ_raw, w, out, how, note=note)
        else:
            print(f"{gen}: no integrator reference recorded -- closure "
                  f"NOT checked")
    # same autoscaling as the muon side (rule 2b) -- see analyze.fmt_sigma
    print(f"{gen}: sigma_fid = {fmt_sigma(out['sigma_fid_pb'])}"
          f"  ({out.get('n_fiducial', 0)} events selected)")
    flag_low_stats(out, gen)   # same guard as the muon side (rule 2b)

    out["means"] = {kk: float(np.average(np.asarray(obs[kk], dtype=float),
                                         weights=w))
                    for kk in ("nch", "nch1", "nch05", "Q2", "xbj", "y")}
    # The DeltaPhi bookkeeping the muon side has always written: how many
    # fiducial events had no charged hadron at all, and the mean of the rest.
    # Added on the neutrino side 2026-09-08 with the observable itself.
    if "dphi" in obs and len(obs["dphi"]):
        dsum = dphi_summary(obs, w)
        out.update(dsum)
        out["means"]["dphi"] = dsum["mean_dphi"]

    msum = dmeson_summary(obs, w)
    out.update(msum)
    out["means"].update({"nd_ch": msum.get("mean_n_d_ch"),
                         "nd_0": msum.get("mean_n_d_0"),
                         "nd_s": msum.get("mean_n_d_s")})
    print(f"{gen}: <N(D+-)> = {msum.get('mean_n_d_ch', 0):.4f}, "
          f"<N(D0)> = {msum.get('mean_n_d_0', 0):.4f}, "
          f"<N(Ds)> = {msum.get('mean_n_d_s', 0):.4f}; "
          f"sigma fraction with >=1: D+- {100*msum.get('frac_ge1_d_ch', 0):.2f}%"
          f", D0 {100*msum.get('frac_ge1_d_0', 0):.2f}%")

    hists = {}
    for kk, edges in BINS.items():
        if not obs[kk]:
            continue
        v = np.asarray(obs[kk], dtype=float)
        cnt, _ = np.histogram(v, bins=edges, weights=w)
        err2, _ = np.histogram(v, bins=edges, weights=w*w)
        widths = np.diff(edges)
        hists[kk] = {"edges": np.asarray(edges).tolist(),
                     "dsig": (cnt*norm/widths).tolist(),
                     "err": (np.sqrt(err2)*norm/widths).tolist()}
    if not msum.get("d_supported"):
        for kk in D_HIST_KEYS:
            hists.pop(kk, None)
        print(f"{gen}: NO D mesons in the sample -- D histograms omitted")
    out["hists"] = hists

    os.makedirs(OUTDIR, exist_ok=True)
    # PER-ENERGY RESULT KEY.  Without this a 400 GeV run overwrites the
    # 1 TeV file: it happened, and three anchor results had to be restored
    # from git.  The energy_gev stamp is what made it visible -- 3.191 pb at
    # 1000 GeV replaced by 0.645 pb at 400 GeV, under the same filename.
    ofn = (f"{OUTDIR}/histos_{at_energy(gen + sel_suffix())}"
           f"{selection_mod.q2_suffix()}.json")
    if gen.startswith("v2_"):
        # histos_<gen>_q4w3_<p|n>[_tag].json -- analysis/sample_layout.py
        ofn = (f"{OUTDIR}/histos_"
               f"{at_energy(sample_layout.result_key(gen[3:], _t, sample_layout.require()[0].name))}.json")
        out["target"] = _t
    with open(ofn, "w") as f:
        stamp(out, current="nu")
        json.dump(stamp_shower_charm(out, tag_stats, w, subtract), f)
    print(f"wrote {ofn}")


if __name__ == "__main__":
    main()
