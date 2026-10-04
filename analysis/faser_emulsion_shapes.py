#!/usr/bin/env python3
"""FASER's own emulsion selection, on our events, at FASER's own settings.

WHAT THIS IS FOR (user, 2026-09-08): "when comparing with FASER emulsion data
we should use exactly their settings, cuts etc, and not those of the benchmark
setting.  IF this requires to regenerate some numbers, so be it."  And, on the
four questions that raised: the 681.1 kg Moriond note rather than the 128.6 kg
paper, "shapes normalised to observed count, so that efficiency cancels",
"only POWHEG-V2 and GENIE are sufficient", and "yes, we want both nu and nubar
samples, to compare as correctly as possible with the data".

>>> THIS IS NOT THE BENCHMARK'S FIDUCIAL REGION, AND THAT IS THE POINT. <<<
Every other selection in this repository carries Q2 > 4 GeV2 and
0.2 < y < 0.9.  FASER imposes neither.  What they do impose, from
arXiv:2403.12520 and CERN-FASER-CONF-2026-002, is a VERTEX TOPOLOGY:

  * tracks are reconstructed with tan(theta) <= 0.5 and momentum above 1 GeV;
  * a vertex needs MORE THAN FOUR such tracks;
  * MORE THAN THREE of them must have tan(theta) <= 0.1;
  * one track is the lepton: momentum above 200 GeV with tan(theta) > 0.005;
  * and the lepton must be back to back with the rest, Delta phi > pi/2
    against the sum of the other tracks in the vertex.

THE TRACK COUNTS INCLUDE THE LEPTON, which is where this differs from the
benchmark's Tier E by one track.  Tier E was built from the same paper but
counts CHARGED HADRONS -- n05 >= 5 and n01 >= 4 -- so it is one track
stricter than FASER on both counts.  That was the right call for a benchmark
selection, whose job is to be stated once and never move; it is the wrong one
for a comparison with the data, so this module counts tracks the way the
measurement does and leaves Tier E alone.

WHY THE SHAPES AND NOT THE YIELD.  Their selection efficiency is 20-40% and
comes from their detector simulation; we cannot reproduce it and will not
guess it.  A distribution normalised to the observed number of events divides
it out to the extent that it is flat across the bin, which is the honest
comparison available -- and it is the comparison their own Figs. 9 and 10
make, where the simulation is "normalised to the number of observed events".

WHY A FLUX FOLD AND NOT ONE BEAM ENERGY.  The 33 muon-neutrino candidates
span 200 GeV to about 3 TeV.  A fixed-energy sample would get the lepton
momentum spectrum wrong by construction, so each generator is run at a ladder
of beam energies on a proton and on a neutron, for neutrinos and
antineutrinos, and the per-event tables are combined with weight
    flux(E) x sigma_selected(E) x (74 p + 110 n),
the same recipe the dimuon cut flow uses.

Usage:
  analysis/faser_emulsion_shapes.py extract --beam=14 --target=p --energy=1000
      --generator=powheg_v2 <events.hepmc> <out.npz>
  analysis/faser_emulsion_shapes.py combine        -> results_nu/faser_emulsion_shapes.json

The ladders (paper plot 12, user 2026-09-14): the POWHEG-V2 ladder is the one
re-showered with LesHouches:matchInOut = off, with its NEUTRON points
regenerated on a genuine 2112 beam ($POWHEG_V2/ladder-faser-v2,
tools/faser_emulsion_ladder.sh); Sherpa's neutron points are the 2112-beam
re-runs ($SHERPA_RUNS/V2_NuDIS_NLO_FASER_*_n_*); GENIE's ladder always used
the neutron target code and is unchanged.  The fixed-energy scale band comes
from the 1 TeV Tier E join (mhou_hadron_q4w3_faser_e_W.json); the energy
remix piece still reads faser_powheg_rates_rwgt.json, whose sigma(E) scale ratios
do not depend on the remnant or the shower.
"""
import json
import math
import os
import re
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import analyze_nu                                            # noqa: E402
import beams                                                 # noqa: E402
import paths                                                 # noqa: E402
import selection                                             # noqa: E402
import faser_rates as fr                                     # noqa: E402
from analyze import (track_tan, CHARGED_HADRONS, M_P, dis_invariants,   # noqa: E402
                     dot, lab_energy, lepton_theta)

# FASER's vertex selection, as published.  Named constants rather than
# literals in the code, so the one place they are written is the one place
# they would have to change.
TRACK_E_MIN = selection.TRACK_E_MIN   # GeV, the one track threshold
TRACK_TAN_MAX = 0.5          # track angular acceptance
TRACK_TAN_NARROW = 0.1       # the narrow subset the background cut uses
N_TRACKS_MIN = 5             # "more than four tracks", lepton included
N_NARROW_MIN = 4             # "more than three", lepton included
LEP_P_MIN = 200.0            # GeV
LEP_TAN_MIN = 0.005
DPHI_MIN = math.pi / 2.0


def extract(hepmc, out, beam, target, energy, generator):
    """One HepMC file -> one per-event table of FASER's vertex observables.

    Kept per event: the neutrino energy, the weight, the lepton momentum and
    tan(theta), the two track counts INCLUDING the lepton, the azimuth
    between the lepton and the summed other tracks, and whether the event
    passes FASER's selection.  Everything their Figs. 9 and 10 plot, plus
    what is needed to re-cut.
    """
    beam = int(beam)
    lep_out = 13 if beam > 0 else -13
    analyze_nu.NU_BEAM_PID = beam
    keys = ("enu", "w", "q2", "y", "plep", "tanlep", "ntr", "ntrn",
            "dphi", "pass")
    cols = {k: [] for k in keys}
    n_parsed = n_no_primary = 0
    for w, k, P, parts, _dmes, _hard in analyze_nu.parse_hepmc3(hepmc):
        n_parsed += 1
        w = 1.0 if w is None else w
        best, best_e = None, -1.0
        for pid, p in parts:
            if pid == lep_out:
                e = lab_energy(p, P)
                if e > best_e:
                    best, best_e = p, e
        if best is None:
            n_no_primary += 1
            continue
        q2, y, _x, kP = dis_invariants(k, P, best)
        tan_lep = math.tan(lepton_theta(best_e, best))
        phi_lep = math.atan2(best[2], best[1])
        # THE LAB ANGLE IS REBUILT FROM INVARIANTS, never read off a
        # component: theta is not invariant under the boost from the
        # generation frame, and Sherpa and Herwig generate in the c.m.  p_T
        # is unchanged by a boost along z, the lab energy is p.P/M_p and the
        # mass is p.p, so |p|_lab and p_z,lab follow.  Same construction as
        # analyze.py's Tier E counts, deliberately.
        n_tr = n_narrow = 0
        sum_px = sum_py = 0.0
        for pid, p in parts:
            if abs(pid) not in CHARGED_HADRONS:
                continue
            e_h = dot(p, P) / M_P
            if e_h <= TRACK_E_MIN:
                continue
            _pt, tan_h = track_tan(p, e_h, P, k)   # signed p_z,lab (2026-09-19)
            if tan_h < TRACK_TAN_MAX:
                n_tr += 1
                if tan_h < TRACK_TAN_NARROW:
                    n_narrow += 1
                sum_px += p[1]
                sum_py += p[2]
        # THE LEPTON IS A TRACK.  It is reconstructed like any other and it
        # is inside both angular windows whenever it passes its own cuts, so
        # it counts once in each.  This is the one place this module differs
        # from Tier E; see the module docstring.
        if tan_lep < TRACK_TAN_MAX:
            n_tr += 1
            if tan_lep < TRACK_TAN_NARROW:
                n_narrow += 1
        dphi = (abs(math.remainder(math.atan2(sum_py, sum_px) - phi_lep,
                                   2.0 * math.pi))
                if (sum_px or sum_py) else -1.0)
        ok = (best_e > LEP_P_MIN and tan_lep > LEP_TAN_MIN
              and n_tr >= N_TRACKS_MIN and n_narrow >= N_NARROW_MIN
              and dphi >= DPHI_MIN)
        cols["enu"].append(kP / M_P)
        cols["w"].append(w)
        cols["q2"].append(q2)
        cols["y"].append(y)
        cols["plep"].append(best_e)
        cols["tanlep"].append(tan_lep)
        cols["ntr"].append(n_tr)
        cols["ntrn"].append(n_narrow)
        cols["dphi"].append(dphi)
        cols["pass"].append(ok)
    meta = {"beam": beam, "target": target, "energy": float(energy),
            "generator": generator, "source": os.path.abspath(hepmc),
            "n_parsed": n_parsed, "n_no_primary": n_no_primary,
            "n_kept": len(cols["w"]),
            "cuts": {"track_e_min": TRACK_E_MIN,
                     "track_tan_max": TRACK_TAN_MAX,
                     "track_tan_narrow": TRACK_TAN_NARROW,
                     "n_tracks_min": N_TRACKS_MIN,
                     "n_narrow_min": N_NARROW_MIN,
                     "lep_p_min": LEP_P_MIN, "lep_tan_min": LEP_TAN_MIN,
                     "dphi_min": DPHI_MIN,
                     "note": "FASER's vertex selection; track counts include "
                             "the lepton; no Q2 cut and no y window"}}
    arrays = {k: np.asarray(v) for k, v in cols.items()}
    np.savez_compressed(out, meta=json.dumps(meta), **arrays)
    n = len(arrays["w"])
    eff = float(arrays["pass"].mean()) if n else 0.0
    print(f"{n_parsed} parsed, {n_no_primary} without a primary lepton, "
          f"{n} kept, {100*eff:.2f}% pass FASER's selection")


# --------------------------------------------------------------- combining
LADDER_ES = [30, 60, 100, 200, 400, 700, 1000, 2000, 4000, 6800]
# FASER's own numbers, from CERN-FASER-CONF-2026-002 (681.1 kg, 9.5 fb^-1).
OBSERVED = {"numu": 33, "nue": 7}
# Table 2: observed muon-neutrino events per RECONSTRUCTED energy bin, and
# their own truth-level expectation.  Tabulated, not read off a figure --
# which is why this is the one shape comparison made against numbers.
ENERGY_BINS = [(10.0, 200.0, 0, 0.00), (200.0, 1000.0, 19, 22.3),
               (1000.0, 1e9, 14, 17.9)]

# The observables, their binning, and the axis each is drawn on.  The binning
# follows FASER's Figs. 9 and 10 as closely as their axes can be read: the
# multiplicity in unit steps, the lepton angle to tan(theta) = 0.1, the lepton
# momentum to 2 TeV, and the azimuth in degrees above the cut.
# >>> THE BINNING IS FASER'S OWN, READ FROM THE DIGITISED FILE. <<<  A shape
# comparison in bins the measurement does not use is a comparison of two
# different objects; taking the edges from data/faser_emulsion/ means every
# panel is bin-for-bin against the published one, and it means the binning
# cannot drift away from the data by an edit here.  Delta phi is stored in
# DEGREES in that file and in radians in the event tables, so the conversion
# is done once, here.
DATA_FILE = "data/faser_emulsion/faser_conf2026_shapes.json"
OBS_KEY = {"ntracks": "ntr", "tanlep": "tanlep", "plep": "plep",
           "dphi": "dphi"}
DEG = {"dphi": math.pi / 180.0}


def _data():
    with open(f"{BASE}/{DATA_FILE}") as f:
        return json.load(f)


def _edges(data, flavour, obs):
    """FASER's own bin edges for one observable, in the event tables' units."""
    e = np.array(data[flavour][obs]["edges"], dtype=float)
    return e * DEG.get(obs, 1.0)


def _point(path):
    """(meta, arrays) for one ladder point, or None."""
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=False)
    return json.loads(str(z["meta"])), {k: z[k] for k in z.files
                                        if k != "meta"}


def _point_shapes(d, bins):
    """(selected weight fraction, {obs: normalised histogram}) at one energy."""
    w = d["w"].astype(float)
    ok = d["pass"].astype(bool)
    tot = float(w.sum())
    if tot <= 0 or not ok.any():
        return 0.0, {k: np.zeros(len(e) - 1) for k, e in bins.items()}
    frac = float(w[ok].sum()) / tot
    out = {}
    for k, edges in bins.items():
        h, _ = np.histogram(d[OBS_KEY[k]][ok], bins=edges, weights=w[ok])
        s = h.sum()
        out[k] = h / s if s > 0 else h
    return frac, out


# THE SHERPA LADDER IS ITS OWN RUN DIRECTORY PER POINT, and its cross-section
# has to be read out of the run log: the event file is deleted after the
# per-event table is extracted (rule 1), so there is no header left to read.
# The line wanted is the nominal entry of the end-of-run summary table, not
# the running "XS = ..." estimates above it, which are the integrator's
# progress and stop at whatever the last event left.
_SHERPA_XS = re.compile(r"Nominal\s.*?\s([0-9.eE+-]+)\s")


def sherpa_total_pb(rundir):
    """Sherpa's own delivered cross-section for one ladder point, in pb.

    >>> WITH THE NEUTRINO SPIN FACTOR APPLIED. <<<  Sherpa 3.0.5 averages the
    matrix element over two initial-state lepton helicities even for a
    neutrino beam, where only one exists, so its integrator sigma is exactly
    half the physical one; analyze_nu carries the same factor for every other
    Sherpa neutrino number here.  It cancels out of a normalised shape, but a
    number that is half the cross-section should not be called one.
    """
    p = os.path.join(rundir, "gen.log")
    if not os.path.exists(p):
        return None
    txt = open(p, errors="replace").read()
    # strip the ANSI colouring the summary table is drawn with
    txt = re.sub(r"\x1b\[[0-9;]*m", "", txt)
    m = _SHERPA_XS.search(txt)
    if not m:
        return None
    v = float(m.group(1)) * analyze_nu.SHERPA_NU_SPIN_FACTOR
    # A REGEX THAT MATCHES THE WRONG COLUMN RETURNS A PERFECTLY PLAUSIBLE
    # NUMBER.  The summary row also carries a relative deviation of "0 %" and
    # an absolute error; picking either would give a cross-section of zero or
    # of a few per cent of the right one, and the fold would quietly reweight
    # this energy point to nothing.  A non-positive value is impossible here,
    # so it is an error rather than a missing point.
    if v <= 0.0:
        raise SystemExit(f"{p}: read a cross-section of {v} -- the summary "
                         f"line was parsed wrongly")
    return v


def herwig_ladder_sigma(cur, tgt, e):
    """(sigma, err) in pb of one Herwig FASER-ladder point: the positive
    minus the negative half (herwig7/herwig_xsec.py, inverse variance over
    jobs), or None if the point is incomplete."""
    sys.path.insert(0, os.path.join(BASE, "herwig7"))
    import herwig_xsec as hx
    tag = beams.Beams("nu", float(e)).tag
    pos = hx.read_jobs(f"{BASE}/herwig7/fl{cur}pwg_{tgt}_job_{tag}_*/"
                       f"FL-{cur}-PWG-{tgt}-{tag}-S*.out", verbose=False)
    neg = hx.read_jobs(f"{BASE}/herwig7/fl{cur}pwgneg_{tgt}_job_{tag}_*/"
                       f"FL-{cur}-PWGNEG-{tgt}-{tag}-S*.out", verbose=False)
    if len(pos) < 2 or not neg:
        return None
    sp, ep, _n, _s = hx.combine(pos)
    sn, en, _n, _s = hx.combine(neg)
    return sp - sn, math.hypot(ep, en)


def _herwig_point(cur, tgt, e):
    """(meta, merged per-event table of the positive jobs), sigma -- the
    Herwig FASER ladder (tools/herwig_faser_ladder.sh, 2026-10-04).  The
    negative half (~0.1% of the CC rate) enters the cross-section only."""
    tag = beams.Beams("nu", float(e)).tag
    parts = [_point(f"{BASE}/herwig7/fl{cur}pwg_{tgt}_job_{tag}_{j}/faserdata_events.npz")
             for j in (1, 2)]
    sig = herwig_ladder_sigma(cur, tgt, e)
    if not all(parts) or sig is None:
        return None, None
    merged = {k: np.concatenate([p[1][k] for p in parts]) for k in parts[0][1]}
    return (parts[0][0], merged), sig[0]


def _ladder(generator, cur, tgt, bins):
    """[(E, sigma_tot_pb, selected fraction, {obs: shape})] for one arm."""
    import faser_dimuon_cutflow as fdc
    pid = 14 if cur == "nu" else -14
    pts = []
    for e in LADDER_ES:
        if generator == "powheg_v2":
            r = f"{paths.POWHEG_V2}/ladder-faser-v2/{cur}_{tgt}_E{e}"
            got = _point(os.path.join(r, "faserdata_events.npz"))
            sig = fdc.powheg_total_pb(r) if got else None
        elif generator == "sherpa_nlo":
            # neutrons: the genuine-2112 re-runs; protons: 2212 beams
            r = (f"{paths.SHERPA_RUNS}/V2_NuDIS_NLO_FASER_{cur}_{tgt}_E{e}"
                 if tgt == "n"
                 else f"{paths.SHERPA_RUNS}/NuDIS_NLO_FASER_{cur}_{tgt}_E{e}")
            got = _point(os.path.join(r, "faserdata_events.npz"))
            sig = sherpa_total_pb(r) if got else None
        elif generator == "herwig_nlo":
            got, sig = _herwig_point(cur, tgt, e)
        else:
            got = _point(f"{BASE}/genie/fdladder_job_{cur}_{tgt}_E{e}_1/"
                         f"faserdata_events.npz")
            sig = (float(fdc.genie_sigma_tot(pid, tgt, np.array([float(e)]))[0])
                   if got else None)
        if not got or not sig:
            continue
        frac, shapes = _point_shapes(got[1], bins)
        pts.append((float(e), float(sig), frac, shapes))
    return pts


def _fold(generator, flavour, kscale=None):
    """The flux-folded shapes, and the events per energy bin, for one flavour.

    `kscale`, when given, is {arm: (E grid, factor)} from `_scale_factors()`
    for ONE of the seven scale choices: the fiducial cross-section at that
    scale over the nominal one, per energy.  It reweights the fold and
    nothing else, which is the point -- see the MHOU block below.

    THE FOLD IS OVER THE FLUX AND OVER BOTH CHARGES AND BOTH NUCLEONS.  At
    each flux point the number of selected interactions is
    flux(E) x sigma_tot(E) x (selected fraction)(E), summed over nu and nubar
    with their own fluxes and over 74 protons and 110 neutrons; the shape at
    that point is the ladder shape interpolated in log E.  Normalisation
    cancels at the end -- the result is scaled to the observed count -- so
    what this has to get right is the RELATIVE weighting, which is exactly
    where the charge mixture and the energy spectrum enter.
    """
    import faser_dimuon_cutflow as fdc
    data = _data()
    bins = {o: _edges(data, flavour, o) for o in OBS_KEY}
    pid_nu, pid_nb = (14, -14) if flavour == "numu" else (12, -12)
    grid = np.geomspace(20.0, 8000.0, 300)
    shapes = {k: np.zeros(len(e) - 1) for k, e in bins.items()}
    ebins = np.zeros(len(ENERGY_BINS))
    total = 0.0
    for cur, pid in (("nu", pid_nu), ("nubar", pid_nb)):
        # THE FLUX IS THE MUON-NEUTRINO ONE FOR BOTH FLAVOURS ONLY IF THE
        # FILE EXISTS; the electron-neutrino flux is a different file and a
        # different spectrum, so it is read by its own pid.
        e_f, phi = fr.flux(str(pid))
        w_flux = np.interp(grid, e_f, phi, left=0.0, right=0.0)
        for tgt, n_nuc in (("p", fdc.Z_W), ("n", fdc.N_W)):
            # THE LADDER IS ALWAYS THE MUON-NEUTRINO ONE.  At these energies
            # the charged-current cross-section and the hadronic final state
            # are flavour independent to far better than the comparison's
            # precision -- the lepton mass is the only difference and it is
            # 0.1 GeV against 200 -- so the electron-neutrino prediction uses
            # the same events with the electron-neutrino FLUX.  Said here
            # because it is an approximation, not an identity.
            lad = _ladder(generator, cur, tgt, bins)
            if not lad:
                continue
            es = np.array([p[0] for p in lad])
            sig = np.array([p[1] for p in lad])
            frac = np.array([p[2] for p in lad])
            # sigma_tot log-log, zero below the ladder; the selected fraction
            # linear in log E and held at the ends -- the convention
            # faser_dimuon_cutflow's convolve() uses, kept identical here.
            with np.errstate(divide="ignore"):
                ls = np.interp(np.log(grid), np.log(es), np.log(sig))
            sig_g = np.where(grid < es[0], 0.0, np.exp(ls))
            if kscale is not None:
                ke, kv = kscale[f"{cur}_{tgt}"]
                sig_g = sig_g * np.interp(grid, ke, kv)
            frac_g = np.interp(np.log(grid), np.log(es), frac)
            wk = w_flux * sig_g * frac_g * n_nuc
            total += float(wk.sum())
            for i, (lo, hi, _o, _t) in enumerate(ENERGY_BINS):
                ebins[i] += float(wk[(grid >= lo) & (grid < hi)].sum())
            # the shape at each grid point: nearest two ladder points in log E
            idx = np.clip(np.searchsorted(es, grid), 1, len(es) - 1)
            lo_e, hi_e = es[idx - 1], es[idx]
            t = np.clip((np.log(grid) - np.log(lo_e))
                        / (np.log(hi_e) - np.log(lo_e)), 0.0, 1.0)
            for k in bins:
                m = np.array([p[3][k] for p in lad])
                interp = (1.0 - t)[:, None] * m[idx - 1] + t[:, None] * m[idx]
                shapes[k] += (wk[:, None] * interp).sum(axis=0)
    if total <= 0:
        return None
    for k in shapes:
        s = shapes[k].sum()
        if s > 0:
            shapes[k] = shapes[k] / s
    return {"shapes": shapes, "bins": {k: v.tolist() for k, v in bins.items()},
            "energy_bins": (ebins / ebins.sum()).tolist(),
            "n_selected_arbitrary": total}



# --------------------------------------------------------------- the MHOU
# >>> WHY THE BAND ON THIS FIGURE IS BUILT FROM TWO PIECES (user, 2026-09-08:
# "I don't see MHOUs in the POWHEG predictions?"). <<<
#
# The shapes here are NORMALISED to the observed candidate count, so an
# overall scale shift cancels and only a shape distortion survives.  Two
# things distort the shape:
#
#   (a) THE ENERGY REMIX.  Each energy enters the fold with weight
#       flux(E) x sigma(E) x (selected fraction)(E), and sigma(E) moves with
#       the scale by a per-cent or two that is ENERGY DEPENDENT -- so the
#       mixture of energies changes and with it the folded shape.  This piece
#       is computed EXACTLY, from POWHEG-V2's own reweighting of every ladder
#       point (results_nu/faser_powheg_rates_rwgt.json).
#
#   (b) THE SHAPE AT FIXED ENERGY.  At one energy the distributions themselves
#       move with the scale.  Computing this exactly here would need the
#       Les Houches scale weights joined to the SHOWERED events by lhe_index,
#       and the ladder's HepMC files were pruned once their per-event tables
#       existed -- so it would mean reweighting, re-showering and re-extracting
#       all twenty ladder points.  Instead it is taken from
#       results_nu/mhou_hadron_faser_e.json, which is that join done properly
#       at 1 TeV under the benchmark's own emulsion-like tier.
#
# >>> (b) IS AN APPROXIMATION AND THE FIGURE SAYS SO. <<<  It is measured at
# ONE energy and under a selection that differs from FASER's in three ways:
# the Q2 floor, the y window, and counting the lepton among the tracks.  It is
# also the LARGER of the two pieces by a wide margin -- (a) is below a per
# cent -- so the band is dominated by a number imported from a neighbouring
# region.  What it is NOT is absent, and what it shows is the right size.
MHOU_MAP = {"ntracks": ("nch05", 1), "plep": ("Emu", 0),
            "dphi": ("dphix", 0), "tanlep": ("theta", 0)}


def _scale_factors():
    """k[arm][j](E): sigma at scale j over sigma at the nominal, per energy."""
    p = f"{BASE}/results_nu/faser_powheg_rates_rwgt.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    out = {}
    for arm, pts in d["ladder"].items():
        if not all(q.get("sigma_fid_by_scale_pb") for q in pts):
            return None
        e = np.array([q["E_GeV"] for q in pts], float)
        v = np.array([q["sigma_fid_by_scale_pb"]["ally"] for q in pts], float)
        out[arm] = (e, v / v[:, :1])
    return out


def _fixed_energy_band(obs, edges):
    """(rel_lo, rel_hi) per bin of `edges`, from the 1 TeV lhe_index join.

    >>> IT GOES THROUGH mhou_band, NOT STRAIGHT AT THE STORED HISTOGRAM. <<<
    That module merges adjacent bins until the envelope beats its own Monte
    Carlo error, and here that matters: read raw, one Delta phi bin came out
    at 29% -- which is one event carrying a reweighting factor of 156, not a
    theory uncertainty.  Merged, the same region is 5%.

    Two mappings are needed.  The emulsion's track count INCLUDES the lepton
    and the benchmark's nch05 does not, so the offset shifts the edges.  And
    the two grids do not span the same range -- FASER's lepton momentum runs
    to 2 TeV where the benchmark's stops at 800 GeV -- so a bin off the end
    takes the band of the nearest bin that exists rather than zero, which
    would read as "no uncertainty here".
    """
    import mhou_band
    if obs not in MHOU_MAP:
        return None
    key, off = MHOU_MAP[obs]
    b = mhou_band.rel("results_nu", key, "q4w3_faser_e", esuffix="_W")
    if b is None:
        return None
    be, blo, bhi = b
    mid = 0.5 * (np.asarray(edges[:-1]) + np.asarray(edges[1:])) - off
    idx = np.clip(np.searchsorted(be, mid) - 1, 0, len(blo) - 1)
    return blo[idx], bhi[idx]


def _mhou(generator, flavour, bins, nominal):
    """{obs: (rel_lo, rel_hi)} for one arm, the two pieces in quadrature.

    Piece (a) is seven folds -- one per scale choice -- of the SAME ladder
    shapes with the SAME flux, differing only in sigma(E); the envelope of
    the resulting normalised shapes about the nominal is the energy remix,
    exactly.  Piece (b) is `_fixed_energy_band`.  They are independent in the
    sense that matters here (one is a reweighting between energies, the other
    a distortion within one) so they are added in quadrature; in practice (b)
    dominates by an order of magnitude and the choice barely shows.
    """
    if generator != "powheg_v2":
        return None
    k = _scale_factors()
    if k is None:
        return None
    nsc = len(next(iter(k.values()))[1][0])
    lo = {o: np.zeros(len(bins[o]) - 1) for o in bins}
    hi = {o: np.zeros(len(bins[o]) - 1) for o in bins}
    for j in range(1, nsc):
        r = _fold(generator, flavour,
                  kscale={a: (e, v[:, j]) for a, (e, v) in k.items()})
        if r is None:
            return None
        for o in bins:
            with np.errstate(divide="ignore", invalid="ignore"):
                d = np.where(nominal[o] > 0,
                             r["shapes"][o] / np.where(nominal[o] > 0,
                                                       nominal[o], 1.0) - 1.0,
                             0.0)
            hi[o] = np.maximum(hi[o], d)
            lo[o] = np.minimum(lo[o], d)
    out = {}
    for o in bins:
        b = _fixed_energy_band(o, bins[o])
        if b is None:
            out[o] = (lo[o].tolist(), hi[o].tolist())
            continue
        blo, bhi = b
        out[o] = ((-np.hypot(lo[o], blo)).tolist(),
                  np.hypot(hi[o], bhi).tolist())
    return out


def combine():
    """Write results_nu/faser_emulsion_shapes.json."""
    out = {"what": "FASER's own emulsion selection, flux-folded over nu and "
                   "nubar on tungsten, shapes normalised to the observed "
                   "candidate counts",
           "source": "CERN-FASER-CONF-2026-002 (681.1 kg, 9.5 fb^-1)",
           "observed": OBSERVED,
           "energy_bins": [{"lo": lo, "hi": hi, "observed": o,
                            "faser_truth_expected": t}
                           for lo, hi, o, t in ENERGY_BINS],
           "data": _data(),
           "cuts": {"track_e_min": TRACK_E_MIN,
                    "track_tan_max": TRACK_TAN_MAX,
                    "n_tracks_min": N_TRACKS_MIN,
                    "n_narrow_min": N_NARROW_MIN,
                    "lep_p_min": LEP_P_MIN, "lep_tan_min": LEP_TAN_MIN,
                    "dphi_min": DPHI_MIN},
           "generators": {}}
    for gen in ("powheg_v2", "genie", "sherpa_nlo", "herwig_nlo"):
        out["generators"][gen] = {}
        for fl in ("numu", "nue"):
            r = _fold(gen, fl)
            if r is None:
                continue
            out["generators"][gen][fl] = {
                "shapes": {k: (v * OBSERVED[fl]).tolist()
                           for k, v in r["shapes"].items()},
                "bins": r["bins"],
                "energy_bin_fractions": r["energy_bins"]}
            bins = {o: np.array(v) for o, v in r["bins"].items()}
            m = _mhou(gen, fl, bins, r["shapes"])
            if m is not None:
                out["generators"][gen][fl]["mhou_rel"] = m
    p = f"{BASE}/results_nu/faser_emulsion_shapes.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {p}")
    for gen, d in out["generators"].items():
        for fl, v in d.items():
            print("  %-10s %-5s energy split %s (FASER truth %s, observed %s)"
                  % (gen, fl,
                     ", ".join("%.3f" % x for x in v["energy_bin_fractions"]),
                     ", ".join("%.3f" % (t / 40.2)
                               for _l, _h, _o, t in ENERGY_BINS),
                     ", ".join("%.3f" % (o / OBSERVED["numu"])
                               for _l, _h, o, _t in ENERGY_BINS)
                     if fl == "numu" else "n/a, single bin"))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "extract":
        kw = {}
        pos = []
        for a in sys.argv[2:]:
            if a.startswith("--"):
                k, _, v = a[2:].partition("=")
                kw[k] = v
            else:
                pos.append(a)
        if len(pos) != 2:
            sys.exit(__doc__)
        extract(pos[0], pos[1], kw.get("beam", 14), kw.get("target", "p"),
                kw.get("energy", 1000.0), kw.get("generator", "unknown"))
    elif len(sys.argv) > 1 and sys.argv[1] == "combine":
        combine()
    else:
        sys.exit(__doc__)
