#!/usr/bin/env python3
"""Locate the hard 2 -> 2 scattering in a HepMC3 event, for any of the three
generators in this benchmark.

The channel tags in analyze.py / analyze_nu.py originally keyed on PYTHIA's
status codes -- 21 = incoming to the hard process, 23 = outgoing from it.
Neither Sherpa nor Herwig writes those, so the tags silently returned "no hard
record" on their samples.  Checked directly in the event files:

    Pythia   status 21 / 23 records                     -> use them
    Sherpa   hard partons at status 3; the hard vertex
             has incoming (beam-flavour lepton, parton), both status 3
    Herwig   everything at status 11; the hard vertex is
             the one whose incoming list holds the EXCHANGED BOSON
             (gamma/Z for NC, W for CC) alongside the struck parton

Example, Sherpa muon NC (Runs/MuonDIS_LO/job_1/evtfull, 2nd event):

    P  9 -6   13 ... 3        <- incoming lepton, beam energy
    P 10 -6    4 ... 3        <- incoming CHARM: a charm-initiated event
    V -5  1 [9,10]            <- the hard vertex
    P 11 -5   13 ... 3        <- outgoing lepton
    P 12 -5    4 ... 3        <- outgoing charm

Example, Herwig muon NC (herwig7/job_1/events.hepmc, 2nd event):

    P  2 -1   22 ... 11       <- exchanged photon
    P  6  4   -2 ... 11       <- struck ubar out of the proton
    P  7  6   -2 ... 11
    V -3    [2,7]             <- the hard vertex
    P  8 -3   -2 ... 11       <- outgoing ubar

WHY THE TAG DIFFERS BETWEEN NC AND CC.  In NC DIS at LO the t-channel
exchange conserves flavour, so a charm event has charm BOTH incoming and
outgoing, and "c in and c out" is a self-check on the tag.  In CC the W
changes flavour: charm is a FINAL-state particle only (s -> c), and the NC
self-check would reject every event.  See charm-production-benchmark and
neutrino-cc-benchmark in the project notes.
"""

LEPTON_PIDS = {11, 12, 13, 14, 15, 16}
NC_BOSONS = {22, 23}
CC_BOSONS = {24, -24}
PARTON_MAX = 6          # |pdg| <= 6 is a quark; 21 is a gluon


def is_parton(pdg):
    return abs(pdg) <= PARTON_MAX or abs(pdg) == 21


def hard_partons(pids):
    """Quark/gluon PDG codes among a list of hard-process records."""
    return [p for p in pids if abs(p) not in LEPTON_PIDS and abs(p) != 22
            and abs(p) != 23 and abs(p) != 24]


def find_hard_vertex(pall, vin, vout, beam_pid, bosons):
    """Return (incoming pdgs, outgoing pdgs) of the hard vertex, or ([], []).

    pall: particle id -> (pdg, status)
    vin:  vertex id -> [incoming particle ids]
    vout: vertex id -> [outgoing particle ids]
    """
    # --- Herwig: the vertex that has the exchanged boson coming in ----------
    for vid, ins in vin.items():
        pdgs = [pall[i][0] for i in ins if i in pall]
        if any(p in bosons for p in pdgs) and any(is_parton(p) for p in pdgs):
            outs = [pall[i][0] for i in vout.get(vid, ()) if i in pall]
            return pdgs, outs

    # --- Sherpa: the vertex fed by a status-3 lepton of the beam flavour ----
    # Several vertices carry a lepton; only the hard one has it at status 3
    # together with a status-3 parton (the shower/hadronisation vertices use
    # status 11).  Prefer the candidate with exactly two incoming legs.
    best = None
    for vid, ins in vin.items():
        recs = [pall[i] for i in ins if i in pall]
        if len(recs) < 2:
            continue
        has_lep = any(abs(p) == abs(beam_pid) and s == 3 for p, s in recs)
        parts = [p for p, s in recs if is_parton(p) and s == 3]
        if has_lep and parts:
            cand = ([p for p, _ in recs],
                    [pall[i][0] for i in vout.get(vid, ()) if i in pall])
            if len(recs) == 2:
                return cand
            if best is None:
                best = cand
    if best is not None:
        return best
    return [], []
