#!/usr/bin/env python3
"""Inspect POWHEG LHE files and showered HepMC samples for the production.

A delivery check, not the analysis: it answers the questions the drivers must
not take on trust (CONVENTIONS.md rule 2) -- which beams are in the file, how many
events, whether the target remnant is the nucleon the run was meant to have --
and it measures the signed-weight fraction of events in the final region,
Q2 > 4 GeV2 and W > 3 GeV, which is what the generation-cut validation
compares.  The analysis proper lives in analysis/ and is not touched here.

    inspect_sample.py lhe   FILE.lhe [...]            [--e-lab E]
    inspect_sample.py hepmc DIR_OR_FILE [...]         [--max N] [--lhe FILE.lhe]
    inspect_sample.py beam  DIR_OR_FILE               -> prints "<lepton> <hadron>"

KINEMATICS.  Q2 = -(k - k')^2 is invariant.  y = q.P / k.P uses the TARGET
record, so it is frame-independent: POWHEG-RES writes a lab-like frame with a
light-like target direction, POWHEG-V2 the c.m. frame, and Pythia keeps
whichever frame the LHE had (memory: "apply every cut in the target rest
frame").  W2 = M_P^2 + y 2kP - Q2 with 2kP = 2 E_lab M_P, as in
analysis/selection.py; for a HepMC file 2 k.P is taken from the beam records.

With --lhe, a HepMC sample is joined to its LHE by the lhe_index attribute
main_powheg writes, and the per-event change of Q2 across the shower is
reported (the Pythia heavy-quark re-massing moves the scattered lepton).

--max-lhe caps how many LHE events are read (default: all).
"""
import glob
import math
import os
import sys

M_P = 0.938272
LEPTONS = {11, 12, 13, 14, 15, 16}
Q2_MIN, W_MIN = 4.0, 3.0


# ------------------------------------------------------------------ charges
_QCHARGE = {1: -1, 2: 2, 3: -1, 4: 2, 5: -1, 6: 2}     # in units of e/3


def charge3(pdg):
    """Three times the electric charge of a PDG code (quarks, leptons,
    gauge bosons, mesons, baryons, diquarks).  Nuclei are not needed here."""
    a = abs(pdg)
    s = 1 if pdg > 0 else -1
    if a <= 6:
        return s * _QCHARGE[a]
    if a in (11, 13, 15):
        return -3 * s
    if a in (12, 14, 16, 21, 22, 23, 25, 111, 130, 310):
        return 0
    if a == 24:
        return 3 * s
    nq3, nq2, nq1 = (a // 1000) % 10, (a // 100) % 10, (a // 10) % 10
    if nq3 == 0:                                   # meson q qbar
        if nq2 == 0 or nq1 == 0:
            return 0
        # the heavier quark (nq2) is the quark for even-up-type convention
        c = _QCHARGE[nq2] - _QCHARGE[nq1]
        if nq2 % 2 == 1:                           # down-type heavier: flip
            c = -c
        return s * c
    if nq1 == 0:                                   # diquark
        return s * (_QCHARGE[nq3] + _QCHARGE[nq2])
    return s * (_QCHARGE[nq3] + _QCHARGE[nq2] + _QCHARGE[nq1])


# ---------------------------------------------------------------------- LHE
def read_lhe(path, nmax=None):
    """Yield (init_ids, init_energies) once, then per event
    (weight, k, k', flavours_out)."""
    ids = energies = None
    with open(path) as f:
        in_init = in_ev = False
        first = False
        for line in f:
            if line.startswith("<init"):
                in_init = True
                first = True
                continue
            if in_init:
                if first:
                    p = line.split()
                    ids = (int(p[0]), int(p[1]))
                    energies = (float(p[2]), float(p[3]))
                    yield ("init", ids, energies)
                    first = False
                if line.startswith("</init"):
                    in_init = False
                continue
            if line.startswith("<event"):
                in_ev, head, k, kp, fl = True, True, None, None, []
                continue
            if line.startswith("</event"):
                in_ev = False
                yield ("event", w, k, kp, fl)
                if nmax is not None:
                    nmax -= 1
                    if nmax <= 0:
                        return
                continue
            if in_ev:
                if head:
                    w = float(line.split()[2])
                    head = False
                    continue
                if line.lstrip().startswith(("#", "<")):
                    continue
                p = line.split()
                if len(p) < 13:
                    continue
                pid, st = int(p[0]), int(p[1])
                mom = tuple(map(float, p[6:10]))
                if abs(pid) in LEPTONS and st == -1:
                    k = mom
                elif abs(pid) in LEPTONS and st == 1:
                    kp = mom
                elif st == 1:
                    fl.append(pid)
        if in_ev:
            yield ("truncated",)


def q2_y(k, kp, pdir=None, P=None):
    """Q2 and y.  Either a target four-momentum P, or (LHE) the target taken
    light-like along -z: q.P / k.P = (q0 + qz) / (k0 + kz)."""
    q = [k[i] - kp[i] for i in range(4)]
    Q2 = -(q[3] ** 2 - q[0] ** 2 - q[1] ** 2 - q[2] ** 2)
    if P is None:
        y = (q[3] + q[2]) / (k[3] + k[2])
        kP = None
    else:
        qP = q[3] * P[3] - q[0] * P[0] - q[1] * P[1] - q[2] * P[2]
        kP = k[3] * P[3] - k[0] * P[0] - k[1] * P[1] - k[2] * P[2]
        y = qP / kP
    return Q2, y, kP


def lhe_summary(path, e_lab=None, nmax=None):
    ids = en = None
    n = 0
    sw = sw_fid = 0.0
    q2s = []
    complete = True
    for rec in read_lhe(path, nmax):
        if rec[0] == "init":
            ids, en = rec[1], rec[2]
            continue
        if rec[0] == "truncated":
            complete = False
            break
        _, w, k, kp, fl = rec
        if k is None or kp is None:
            continue
        Q2, y, _ = q2_y(k, kp)
        n += 1
        sw += w
        q2s.append(Q2)
        if e_lab:
            two_kP = 2.0 * e_lab * M_P
            W2 = M_P ** 2 + y * two_kP - Q2
            if Q2 > Q2_MIN and W2 > W_MIN ** 2:
                sw_fid += w
    if nmax is None:
        with open(path, "rb") as f:
            f.seek(max(0, os.path.getsize(path) - 200))
            complete = complete and b"</LesHouchesEvents>" in f.read()
    q2s.sort()
    return dict(file=path, init_ids=ids, init_energies=en, n_events=n,
                complete=complete, min_q2=q2s[0] if q2s else None,
                frac_q2_below_4=(sum(1 for x in q2s if x < 4) / n) if n else None,
                frac_fid=(sw_fid / sw) if (e_lab and sw) else None)


# -------------------------------------------------------------------- HepMC
def hepmc_events(path, nmax=None):
    """Yield per event dict(weight, lhe_index, particles={id: (pdg, st, p4,
    prodvtx)}, vertices={vid: [incoming ids]})."""
    ev = None
    n = 0
    with open(path) as f:
        for line in f:
            c = line[:2]
            if c == "E ":
                if ev is not None:
                    yield ev
                    n += 1
                    if nmax and n >= nmax:
                        return
                ev = dict(weight=None, lhe_index=None, particles={},
                          vertices={})
            elif ev is None:
                continue
            elif c == "W ":
                ev["weight"] = float(line.split()[1])
            elif c == "P ":
                p = line.split()
                ev["particles"][int(p[1])] = (int(p[3]), int(p[9]),
                                              tuple(map(float, p[4:8])),
                                              int(p[2]))
            elif c == "V ":
                p = line.split()
                ins = p[3].strip("[]").split(",")
                ev["vertices"][int(p[1])] = [int(x) for x in ins if x]
            elif line.startswith("A 0 lhe_index"):
                ev["lhe_index"] = int(line.split()[3])
    if ev is not None:
        yield ev


def event_kinematics(ev):
    """(Q2, W2, beams, charge_sum3, remnant ids) of one showered event."""
    parts = ev["particles"]
    beams = {pid: v for pid, v in parts.items() if v[1] == 4}
    blep = [(pid, v) for pid, v in beams.items() if abs(v[0]) in LEPTONS]
    bhad = [(pid, v) for pid, v in beams.items() if abs(v[0]) not in LEPTONS]
    if len(blep) != 1 or len(bhad) != 1:
        return None
    k, P = blep[0][1][2], bhad[0][1][2]
    # scattered lepton: the final-state lepton whose lepton ancestry reaches
    # the outgoing hard-process lepton (status 23).  Pythia can insert
    # recoil copies (status 42/44) between the two, so "produced AT the hard
    # vertex" is not a usable test.  In HepMC3 ascii the parent field is a
    # particle id when positive and a vertex id when negative.
    def from_hard(pid):
        seen = 0
        while pid and seen < 50:
            seen += 1
            pdg, st, _, par = parts[pid]
            if abs(pdg) not in LEPTONS:
                return False
            if st == 23:
                return True
            if par > 0:
                pid = par
            elif par < 0:
                vin = ev["vertices"].get(par, [])
                # the hard vertex itself: its incoming lepton is status 21
                if any(abs(parts[i][0]) in LEPTONS and parts[i][1] == 21
                       for i in vin):
                    return True
                ins = [i for i in vin if abs(parts[i][0]) in LEPTONS]
                pid = ins[0] if len(ins) == 1 else None
            else:
                return False
        return False
    cand = [v[2] for pid, v in parts.items()
            if abs(v[0]) in LEPTONS and v[1] == 1 and from_hard(pid)]
    if len(cand) != 1:
        return None
    Q2, y, kP = q2_y(k, cand[0], P=P)
    W2 = M_P ** 2 + 2.0 * y * kP - Q2
    ch3 = sum(charge3(v[0]) for v in parts.values() if v[1] == 1)
    remn = [v[0] for v in parts.values() if v[1] == 63]
    return Q2, W2, (blep[0][1][0], bhad[0][1][0]), ch3, remn


def hepmc_files(arg):
    if os.path.isdir(arg):
        return [os.path.join(arg, "events.hepmc")]
    return sorted(glob.glob(arg))


def hepmc_summary(path, nmax=None, lhe=None):
    from collections import Counter
    beams = Counter()
    charge = Counter()
    remn = Counter()
    n = bad = 0
    sw = sw_fid = 0.0
    lheq2 = None
    if lhe:
        lheq2 = []
        for rec in read_lhe(lhe):
            if rec[0] == "event":
                _, w, k, kp, fl = rec
                lheq2.append(q2_y(k, kp)[0] if k and kp else None)
    moved = out_ = in_ = matched = 0
    maxup = 0.0
    for ev in hepmc_events(path, nmax):
        r = event_kinematics(ev)
        n += 1
        if r is None:
            bad += 1
            continue
        Q2, W2, b, ch3, rm = r
        beams[b] += 1
        charge[ch3] += 1
        for x in rm:
            remn[x] += 1
        w = ev["weight"] or 0.0
        sw += w
        if Q2 > Q2_MIN and W2 > W_MIN ** 2:
            sw_fid += w
        if lheq2 is not None and ev["lhe_index"] is not None \
                and ev["lhe_index"] < len(lheq2) and lheq2[ev["lhe_index"]]:
            ql = lheq2[ev["lhe_index"]]
            matched += 1
            rel = (Q2 - ql) / ql
            maxup = max(maxup, rel)
            moved += abs(rel) > 0.01
            out_ += (ql > Q2_MIN and Q2 <= Q2_MIN)
            in_ += (ql <= Q2_MIN and Q2 > Q2_MIN)
    out = dict(file=path, n_events=n, unparsed=bad, beams=dict(beams),
               charge_sum_times3=dict(charge),
               remnants=dict(remn.most_common(6)),
               frac_fid=(sw_fid / sw) if sw else None)
    if lheq2 is not None:
        out.update(lhe_matched=matched, q2_moved_gt_1pct=moved,
                   crossed_out_of_q2_4=out_, crossed_into_q2_4=in_,
                   max_upward_rel_shift=maxup)
    return out


def main():
    a = sys.argv[1:]
    if not a:
        raise SystemExit(__doc__)
    opt = {}
    for key in ("--e-lab", "--max", "--lhe", "--max-lhe"):
        if key in a:
            i = a.index(key)
            opt[key] = a[i + 1]
            del a[i:i + 2]
    mode, args = a[0], a[1:]
    if mode == "lhe":
        e = float(opt["--e-lab"]) if "--e-lab" in opt else None
        nmax = int(opt["--max-lhe"]) if "--max-lhe" in opt else None
        for f in args:
            print(lhe_summary(f, e, nmax))
    elif mode == "hepmc":
        nmax = int(opt["--max"]) if "--max" in opt else None
        for arg in args:
            for f in hepmc_files(arg):
                print(hepmc_summary(f, nmax, opt.get("--lhe")))
    elif mode == "beam":
        for ev in hepmc_events(hepmc_files(args[0])[0], 1):
            r = event_kinematics(ev)
            if r is None:
                raise SystemExit("could not read the beams of the first event")
            print(r[2][0], r[2][1])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
