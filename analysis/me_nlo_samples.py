#!/usr/bin/env python3
"""NLO-side matrix-element-level samples for the muon-DIS benchmark.

Two samples, both parton level (no shower, no hadronisation), written in the
schema of analysis/analyze.py (BINS imported read-only, so the binning can
never drift):

  sherpa_nlo_me  fixed-order NLO QCD from Sherpa 3.0.5
                 (Runs/MuonDIS_NLO_ME: NLO_Mode Fixed_Order, NLO_Part BVIRS).
                 This is a genuine fixed-order calculation, NOT MC@NLO with
                 the shower switched off -- MC@NLO subtraction terms are
                 defined w.r.t. the parton shower, so an unshowered MC@NLO
                 sample is not a well-defined fixed-order prediction.
                 -> results/histos_sherpa_nlo_me.json

  powheg_lhe     the POWHEG-DIS LHE files (Born + POWHEG hardest emission),
                 read directly at parton level.  NLO + first emission, i.e.
                 NOT a pure fixed-order NLO matrix element: it is NLO
                 accurate for inclusive quantities but Sudakov-resummed in
                 the emission variable, so it need not (and does not) sit on
                 top of the analytic NLO curve everywhere.
                 -> results/histos_powheg_lhe.json

Only the six muon-side DIS observables (Q2, xbj, y, Emu, theta, nu) are
filled: there are no hadrons at ME level, so N_ch / E_lead / DeltaPhi have no
meaning (same convention as analysis/yadism_calc.py).

Usage:  me_nlo_samples.py [sherpa_nlo_me|powheg_lhe|both]
"""
import glob
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- config.sh is the one source of every path
from analyze import (BINS, ENERGY, at_energy, job_files, M_MU,  # noqa: E402
                     M_P, Q2_MIN, Y_MAX, Y_MIN, BEAMS, check_anchor_literals,
                     closure_check, dot, powheg_res_sigma_pb, result_path)
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHERPA_RUNS = paths.SHERPA_RUNS
# ENERGY-AWARE, and it was not until 2026-08-27.  This module calls
# at_energy() for the Sherpa event files while pinning BOTH its POWHEG inputs
# and its Sherpa integration log to 1 TeV, which is the "input path made
# energy-aware while something else was left fixed" trap that has now fired
# six times here.  It was latent only because the ME-level samples exist at
# the anchor alone; it would have fired the moment the scan was extended.
POWHEG_LHE = (f"{paths.POWHEG_RES}/parallel-mu{BEAMS.tag}-wide/"
              f"pwgevents-*.lhe")

# Sherpa fixed-order NLO integrator result, summed over the two process
# groups (BVI + RS).  integ.log is the run that read back the converged grids
# from Results.zip and quoted both groups; init.log is the original (killed)
# integration and only ever quoted the BVI group.
SHERPA_NLO_ME_LOG = f"{SHERPA_RUNS}/{at_energy('MuonDIS_NLO_ME')}/integ.log"

# The beam follows BENCH_ENERGY like every other analyser here.  It was a
# 1000.0 literal, so at any other energy this module reconstructed x_Bj and
# the lepton angles against the wrong 2k.P while reading that energy's events.
EBEAM = ENERGY                       # GeV, lab
TWO_KP = BEAMS.two_kP                # exact for a target at rest
ME_KEYS = ("Q2", "xbj", "y", "Emu", "theta", "nu")


# --------------------------------------------------------------- observables
def lab_observables(Q2, y):
    """(Emu, nu) in the lab from (Q2, y) -- exact for a target at rest:
    P.k' = (1-y) P.k with P.k = M_P E_beam."""
    return EBEAM * (1.0 - y), EBEAM * y


def theta_from_pt(pt, e_lab):
    """analyze.py's lab scattering angle: asin(pT / |p_lab|).  pT is invariant
    under the boost along z that connects the generation frame to the lab."""
    pmag = math.sqrt(max(e_lab * e_lab - M_MU * M_MU, 1e-12))
    return math.asin(min(pt / pmag, 1.0))


def theta_from_q2y(Q2, y):
    """Same angle from the invariants (yadism_calc.py convention); used as a
    consistency check on theta_from_pt."""
    ep = EBEAM * (1.0 - y)
    k = math.sqrt(EBEAM**2 - M_MU**2)
    kp = math.sqrt(max(ep * ep - M_MU**2, 1e-12))
    c = (2.0 * EBEAM * ep - 2.0 * M_MU**2 - Q2) / (2.0 * k * kp)
    return math.acos(min(max(c, -1.0), 1.0))


# ------------------------------------------------------- POWHEG LHE reading
def parse_lhe(fname):
    """Yield (weight_pb, k_in, k_out) per event; momenta as (E, px, py, pz).

    The POWHEG LHE is written with massless collinear beams along +-z (the
    incoming muon is (E,0,0,E), the incoming parton (E',0,0,-E')), which is a
    boost of the fixed-target lab along z.  Q2 is a pure invariant of the
    muon line; y follows from the light-cone ratio (see analyse_lhe).
    """
    with open(fname) as f:
        in_ev = False
        head = False
        npart = 0
        w = None
        k_in = k_out = None
        for line in f:
            s = line.strip()
            if s.startswith("<event"):
                in_ev, head, k_in, k_out = True, True, None, None
                continue
            if not in_ev:
                continue
            if s.startswith("</event>"):
                if w is not None and k_in is not None and k_out is not None:
                    yield w, k_in, k_out
                in_ev = False
                continue
            if s.startswith("#") or s.startswith("<"):
                continue
            c = s.split()
            if head:
                if len(c) < 6:
                    continue
                npart, w = int(c[0]), float(c[2])
                head = False
                continue
            if npart <= 0 or len(c) < 10:
                continue
            npart -= 1
            pid, istup = int(c[0]), int(c[1])
            if pid != 13:
                continue
            p = (float(c[9]), float(c[6]), float(c[7]), float(c[8]))
            if istup == -1:
                k_in = p
            elif istup == 1 and (k_out is None or p[0] > k_out[0]):
                k_out = p


def analyse_lhe(files):
    """Histogram the POWHEG-RES LHE parton-level events."""
    obs = {k: [] for k in ME_KEYS}
    wts = []
    n_events = 0
    sum_w_all = 0.0
    dth_max = 0.0
    for fn in files:
        print(f"  parsing {fn}", flush=True)
        for w, k, kp in parse_lhe(fn):
            n_events += 1
            sum_w_all += w
            q = (k[0]-kp[0], k[1]-kp[1], k[2]-kp[2], k[3]-kp[3])
            Q2 = -dot(q, q)
            # P is exactly massless along -z in this frame, P ~ (1,0,0,-1):
            # y = P.q/P.k = (q0 + qz)/(k0 + kz), independent of |P|
            y = (q[0] + q[3]) / (k[0] + k[3])
            if Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX:
                continue
            xbj = Q2 / (y * TWO_KP)          # benchmark convention 2k.P
            e_lab, nu = lab_observables(Q2, y)
            th = theta_from_pt(math.hypot(kp[1], kp[2]), e_lab)
            dth_max = max(dth_max, abs(th - theta_from_q2y(Q2, y)) / th)
            for kk, v in (("Q2", Q2), ("xbj", xbj), ("y", y), ("Emu", e_lab),
                          ("theta", th), ("nu", nu)):
                obs[kk].append(v)
            wts.append(w)
    print(f"  max relative theta(pT) vs theta(Q2,y) mismatch: {dth_max:.2e}")
    return obs, np.asarray(wts), n_events, sum_w_all


# -------------------------------------------------- Sherpa fixed-order NLO
def parse_sherpa_fo(fname):
    """Yield (subs, ntrials, k_beam, P_beam) per *event*, where subs is the
    list of (summed weight, muon momentum) of the correlated NLO sub-events,
    merged over sub-events that share the same muon momentum.

    Sherpa writes each sub-event of an RS event as a separate HepMC3 GenEvent
    carrying the same event number (SHERPA/Tools/HepMC3_Interface.C,
    SubEvtList2ShortHepMC).  In DIS all Catani-Seymour dipoles are
    initial-final / final-initial, so the mapped Born momenta leave the muon
    line untouched: every sub-event of an event has *identical* muon momenta
    and therefore identical Q2, y, ... .  Summing their weights before
    filling is then both the manual's prescription (Sherpa manual,
    "Computing (differential) cross sections of real correction events")
    and exact -- no cancellation at the fiducial cut boundary.
    """
    ev_no = None
    subs = []            # [[weight, muon], ...] merged by muon momentum
    ntr = 0.0
    beams = {}
    mu3 = None
    cur_w = None
    cur_mu = None

    def add_sub():
        """Merge the sub-event just read into subs (by muon momentum)."""
        if cur_w is None or cur_mu is None:
            return
        for s in subs:
            if all(abs(a - b) <= 1e-9 * (abs(a) + abs(b) + 1e-12)
                   for a, b in zip(s[1], cur_mu)):
                s[0] += cur_w
                return
        subs.append([cur_w, cur_mu])

    def flush():
        if ev_no is None or not subs:
            return None
        kb = beams.get(13, mu3)
        pb = beams.get(2212)
        if kb is None:
            return None
        if pb is None:                       # c.m. frame, back-to-back beams
            e_p = math.sqrt(kb[1]**2 + kb[2]**2 + kb[3]**2 + M_P*M_P)
            pb = (e_p, -kb[1], -kb[2], -kb[3])
        return list(subs), ntr, kb, pb

    with open(fname) as f:
        for line in f:
            c0 = line[0:2]
            if c0 == "E ":
                add_sub()
                c = line.split()
                new_no = int(c[1])
                if new_no != ev_no:          # a new event: emit the old one
                    out = flush()
                    if out is not None:
                        yield out
                    ev_no, subs, ntr, beams, mu3 = new_no, [], 0.0, {}, None
                cur_w, cur_mu = None, None
            elif c0 == "W ":
                c = line.split()
                # the HepMC3 header carries a "W <name>\|<name>\|..." line
                # naming the weights; skip it (it precedes any "E " line)
                if ev_no is None or "\\|" in c[1]:
                    continue
                cur_w = float(c[1])
                # EXTRA__NTrials is the 4th named weight, identical for all
                # sub-events of an event -> count it once per event
                if not subs:
                    ntr = float(c[4]) if len(c) > 4 else 1.0
            elif c0 == "P ":
                c = line.split()
                pid, status = int(c[3]), int(c[9])
                p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                if status == 4 or status == 11:
                    if pid == 13 and (13 not in beams or p[0] > beams[13][0]):
                        beams[13] = p
                    elif pid == 2212:
                        beams[2212] = p
                elif status == 3 and pid == 13:
                    if mu3 is None or p[0] > mu3[0]:
                        mu3 = p
                elif status == 1 and pid == 13:
                    if cur_mu is None or p[0] > cur_mu[0]:
                        cur_mu = p
        add_sub()
        out = flush()
        if out is not None:
            yield out


def analyse_sherpa_fo(files):
    """Fill the histograms from the fixed-order NLO event files.

    Book-keeping follows the Sherpa manual: an "event" is a full real
    correction event with all its sub-events; sub-event weights are summed
    before filling.  In DIS every sub-event carries the same muon (see
    parse_sherpa_fo), so the sum is over one kinematic point and the
    statistical error is the plain weighted one; n_multi counts the events
    where that is not the case (expected: zero).
    """
    obs = {k: [] for k in ME_KEYS}
    wts = []
    n_events = 0
    n_sub = 0
    n_multi = 0
    sum_w_all = 0.0
    sum_ntrials = 0.0
    for fn in files:
        print(f"  parsing {fn}", flush=True)
        for subs, ntr, k, P in parse_sherpa_fo(fn):
            n_events += 1
            n_sub += len(subs)
            if len(subs) > 1:
                n_multi += 1
            sum_ntrials += ntr
            for w, kp in subs:
                sum_w_all += w
                q = (k[0]-kp[0], k[1]-kp[1], k[2]-kp[2], k[3]-kp[3])
                Q2 = -dot(q, q)
                kP = dot(k, P)
                Pq = dot(P, q)
                y = Pq / kP
                if Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX:
                    continue
                xbj = Q2 / (2.0 * Pq)
                e_lab = dot(kp, P) / M_P
                th = theta_from_pt(math.hypot(kp[1], kp[2]), e_lab)
                nu = kP / M_P - e_lab
                for kk, v in (("Q2", Q2), ("xbj", xbj), ("y", y),
                              ("Emu", e_lab), ("theta", th), ("nu", nu)):
                    obs[kk].append(v)
                wts.append(w)
    print(f"  {n_events} events, {n_sub} distinct kinematic points "
          f"({n_multi} events with more than one -> {100.0*n_multi/max(n_events,1):.3f}%)")
    return obs, np.asarray(wts), n_events, sum_w_all, sum_ntrials


# ------------------------------------------------------------------- output
def write_json(gen, label, obs, w, norm, sigma_fid, sigma_err, n_parsed,
               extra=None):
    out = {"generator": gen, "label": label, "n_parsed": n_parsed,
           "n_fiducial": len(w), "sigma_fid_pb": sigma_fid,
           "sigma_fid_err_pb": sigma_err}
    if extra:
        out.update(extra)
    out["means"] = {kk: float(np.average(np.asarray(obs[kk], dtype=float),
                                         weights=w))
                    for kk in ("Q2", "xbj", "y")}
    hists = {}
    for kk in ME_KEYS:
        edges = np.asarray(BINS[kk], dtype=float)
        v = np.asarray(obs[kk], dtype=float)
        cnt, _ = np.histogram(v, bins=edges, weights=w)
        err2, _ = np.histogram(v, bins=edges, weights=w*w)
        widths = np.diff(edges)
        hists[kk] = {"edges": edges.tolist(),
                     "dsig": (cnt*norm/widths).tolist(),
                     "err": (np.sqrt(err2)*norm/widths).tolist()}
    out["hists"] = hists
    os.makedirs(f"{BASE}/results", exist_ok=True)
    # result_path(), NOT a bare histos_<genkey>.json.  Making the INPUT job
    # base energy-aware while leaving the OUTPUT name fixed is precisely how a
    # 4 TeV run overwrites the 1 TeV anchor -- which is exactly what happened
    # here on 2026-08-25, caught by the energy_gev stamp rather than by any
    # error.  The two must move together, always.
    ofn = result_path(gen)
    with open(ofn, "w") as f:
        stamp(out, current="mu")
        json.dump(out, f)
    print(f"wrote {ofn}")
    return out


def run_powheg_lhe():
    files = sorted(glob.glob(POWHEG_LHE))
    if not files:
        sys.exit("no POWHEG-RES LHE files found")
    obs, w, n_events, sum_w_all = analyse_lhe(files)
    print(f"powheg_lhe: parsed {n_events} events, {len(w)} in fiducial region")
    # weighted, generation region (Q2 > 2.25, full y) WIDER than fiducial:
    # sigma_fid = integrator sigma x fiducial weight fraction -- exactly the
    # convention of analyze.py's powheg branch (POWHEG_SIGMA_PB)
    sigma_gen = powheg_res_sigma_pb()
    norm = sigma_gen / sum_w_all
    sigma_fid = norm * float(w.sum())
    sigma_err = norm * float(np.sqrt((w * w).sum()))
    print(f"powheg_lhe: <w> = {sum_w_all/n_events:.1f} pb, "
          f"sigma_fid = {sigma_fid/1e3:.3f} +- {sigma_err/1e3:.3f} nb "
          f"({100*w.sum()/sum_w_all:.2f}% of the generated sigma)")
    # THE GATE.  This module wrote two published entries for a year without
    # one, printing the comparison and continuing whatever it said -- which is
    # precisely how a sample normalised against the wrong integrator gets onto
    # the page looking reasonable.  An LHE is read directly, so there is no
    # shower dropout here and the mean weight IS the delivered cross-section.
    extra = {"mean_weight_pb": sum_w_all / n_events, "sigma_gen_pb": sigma_gen}
    closure_check("powheg_lhe", sum_w_all / n_events, sigma_gen, w, extra,
                  how="mean weight over the LHE (no shower, nothing dropped)")
    write_json("powheg_lhe", "POWHEG-RES NLO (LHE, +1st emission)", obs, w, norm,
               sigma_fid, sigma_err, n_events, extra)


def sherpa_sigma_from_log(fname):
    """Total fixed-order cross section from a Sherpa init/run log:
    sum of the per-process-group results (BVI + RS) with errors."""
    tot, err2, parts = 0.0, 0.0, []
    with open(fname, errors="ignore") as f:
        for line in f:
            if "exp. eff:" not in line or " pb " not in line:
                continue
            txt = line.replace("\x1b[34m", "").replace("\x1b[0m", "")
            for tag in ("\x1b[1m", "\x1b[31m", "\x1b[32m"):
                txt = txt.replace(tag, "")
            try:
                name = txt.split(":")[0].strip()
                v = float(txt.split(":")[1].split("pb")[0].strip())
                e = float(txt.split("+-")[1].split("pb")[0].strip("( "))
            except (IndexError, ValueError):
                continue
            parts.append((name, v, e))
    # keep the last quote for each process group
    seen = {}
    for name, v, e in parts:
        seen[name] = (v, e)
    for name, (v, e) in seen.items():
        print(f"  {name}: {v:.2f} +- {e:.2f} pb")
        tot += v
        err2 += e * e
    return tot, math.sqrt(err2)


def run_sherpa_nlo_me(sigma_pb=None, sigma_err_pb=None):
    files = job_files(f"{SHERPA_RUNS}/{at_energy('MuonDIS_NLO_ME')}", "job", "evtme")
    if not files:
        sys.exit("no Sherpa fixed-order NLO event files found")
    if sigma_pb is None:
        sigma_pb, sigma_err_pb = sherpa_sigma_from_log(SHERPA_NLO_ME_LOG)
    obs, w, n_events, sum_w_all, sum_ntrials = analyse_sherpa_fo(files)
    print(f"sherpa_nlo_me: parsed {n_events} events, "
          f"{len(w)} in fiducial region")
    sigma_evt = sum_w_all / sum_ntrials if sum_ntrials else float("nan")
    print(f"sherpa_nlo_me: sum(w)/sum(ntrials) = {sigma_evt/1e3:.4f} nb "
          f"vs integrator {sigma_pb/1e3:.4f} nb")
    # weighted, generation cuts == fiducial cuts: shape from the weights,
    # normalisation pinned to the integrator (as for the other Sherpa samples)
    norm = sigma_pb / sum_w_all
    sigma_fid = sigma_pb
    # sigma_fid is pinned to the integrator, so its uncertainty is the
    # integrator's; the (weighted) event-sample error is kept alongside as an
    # independent cross-check and is what the histogram errors reflect
    stat = norm * math.sqrt(float((w * w).sum()))
    sigma_err = sigma_err_pb or stat
    print(f"sherpa_nlo_me: sigma_fid = {sigma_fid/1e3:.4f} +- "
          f"{sigma_err/1e3:.4f} nb (integrator); event-sample stat error "
          f"{stat/1e3:.4f} nb")
    extra = {"mean_weight_pb": sum_w_all / n_events,
             "sigma_event_estimate_pb": sigma_evt,
             "sigma_event_stat_err_pb": stat,
             "sigma_integrator_pb": sigma_pb,
             "sigma_integrator_err_pb": sigma_err_pb}
    # sigma_fid here is PINNED to the integrator, so the gate is the only
    # thing comparing the events against it -- exactly the configuration
    # CONVENTIONS.md rule 2 was written about, where the cross-section cannot look
    # wrong because it is not measured from the events at all.
    closure_check("sherpa_nlo_me", sigma_evt, sigma_pb, w, extra,
                  how="sum(w)/sum(ntrials)")
    write_json("sherpa_nlo_me", "Sherpa 3.0.5 NLO (ME)", obs, w, norm,
               sigma_fid, sigma_err, n_events, extra)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("powheg_lhe", "both"):
        run_powheg_lhe()
    if which in ("sherpa_nlo_me", "both"):
        run_sherpa_nlo_me()


if __name__ == "__main__":
    main()
