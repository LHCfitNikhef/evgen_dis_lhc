#!/usr/bin/env python3
"""Per-job event check for the paper-plots Sherpa samples.

    tools/sherpa_check_events.py <evtfull> --target p|n [--beam 13|14]

Reads a Sherpa HepMC3 evtfull and FAILS (exit 1) unless every event has

  * the status-4 beams (lepton, target) the run was configured for:
    2212 on a proton run, 2112 on a neutron run;
  * a target breakup vertex whose outgoing partons and diquarks carry the
    target's VALENCE content: (net u, net d) = (2, 1) for the proton,
    (1, 2) for the neutron;
  * no diquark the target cannot have: 1103 (dd) on a proton, 2203 (uu) on
    a neutron.

WHY IT EXISTS.  The neutron runs use a patched COPY of Sherpa
(tools/build_sherpa_npatch.sh), selected through the environment.  Had the
copy not been the one loaded, a 2112 beam would segfault -- but a 2212 beam
with the neutron PDF, which is what the neutron dirs held before
(user decision 2026-09-13), produces perfectly plausible events with a PROTON
remnant, and nothing in the cross-section shows it.  The remnant content is
the one thing that does.  It is checked on the events themselves, event by
event, rather than inferred from the configuration.

It also prints sum(W)/sum(NTrials), the event-level cross-section the
analysis's closure gate uses, the negative-weight fraction, and the diquark
census.  Prints one summary line starting OK or FAIL.
"""
import argparse
import sys
from collections import Counter

DQ = {1103: (0, 2), 2101: (1, 1), 2103: (1, 1), 2203: (2, 0),
      3101: (0, 1), 3103: (0, 1), 3201: (1, 0), 3203: (1, 0), 3303: (0, 0)}
TARGET = {"p": (2212, (2, 1), 1103), "n": (2112, (1, 2), 2203)}


def ud(pid):
    a, s = abs(pid), (1 if pid > 0 else -1)
    if a == 2:
        return (s, 0)
    if a == 1:
        return (0, s)
    if a in DQ:
        return (s * DQ[a][0], s * DQ[a][1])
    return (0, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("evtfull")
    ap.add_argument("--target", required=True, choices=("p", "n"))
    ap.add_argument("--beam", type=int, choices=(13, 14, -14))
    args = ap.parse_args()
    tpid, valence, forbidden = TARGET[args.target]

    nev = nneg = 0
    sw = st = 0.0
    beams, content, diq = Counter(), Counter(), Counter()
    bad = []

    def finish(ev):
        nonlocal nev, nneg, sw, st
        if ev is None:
            return
        nev += 1
        sw += ev["w"]
        st += ev["nt"]
        nneg += ev["w"] < 0
        b4 = tuple(sorted(pid for (_, _, pid, s) in ev["P"] if s == 4))
        beams[b4] += 1
        want = tuple(sorted((args.beam, tpid))) if args.beam else None
        if (want and b4 != want) or tpid not in b4:
            bad.append((ev["n"], f"beams {b4}"))
            return
        ids = {i for (i, _, pid, s) in ev["P"] if s == 11 and pid == tpid}
        vtx = [v for v, ins in ev["V"].items() if ins & ids]
        if not ids or not vtx:
            bad.append((ev["n"], "no target breakup vertex"))
            return
        out = [pid for (i, par, pid, s) in ev["P"] if par == vtx[0] or par in ids]
        c = (sum(ud(x)[0] for x in out), sum(ud(x)[1] for x in out))
        content[c] += 1
        for x in out:
            if abs(x) in DQ:
                diq[x] += 1
        if c != valence:
            bad.append((ev["n"], f"valence {c}"))
        if any(abs(x) == forbidden for x in out):
            bad.append((ev["n"], f"diquark {forbidden}"))

    ev = None
    with open(args.evtfull) as f:
        for line in f:
            c = line[:2]
            if c == "E ":
                finish(ev)
                ev = {"n": int(line.split()[1]), "P": [], "V": {}, "w": 0.0, "nt": 0.0}
            elif ev is None:
                continue
            elif c == "W " and "Weight" not in line:
                x = line.split()
                ev["w"], ev["nt"] = float(x[1]), float(x[4])
            elif c == "P ":
                x = line.split()
                ev["P"].append((int(x[1]), int(x[2]), int(x[3]), int(x[9])))
            elif c == "V ":
                x = line.split()
                if len(x) > 3 and x[3].startswith("["):
                    ev["V"][int(x[1])] = {int(t) for t in x[3].strip("[]").split(",") if t}
    finish(ev)

    ok = nev > 0 and not bad
    print(f"{'OK' if ok else 'FAIL'} {args.evtfull}: target {args.target} ({tpid}), "
          f"{nev} events, bad {len(bad)}")
    print(f"  beams {dict(beams)}")
    print(f"  target-vertex valence (u, d): {dict(content)}")
    print(f"  diquarks: {dict(diq)}")
    if st > 0:
        print(f"  sum(W)/sum(NTrials) = {sw / st:.6g} pb, negative weights "
              f"{nneg / max(nev, 1):.1%}")
    for n, why in bad[:5]:
        print(f"  event {n}: {why}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
