#!/usr/bin/env python3
"""Inspect finished Herwig jobs: is this sample what it claims to be?

    herwig7/production/inspect_sample.py <jobdir> [<jobdir> ...] [--max-events N]

A DIAGNOSTIC, not the analysis: it writes nothing, and the analysis proper
(selection q4w3) is not touched.  Per job it prints

  * the end-of-run bookkeeping through herwig_xsec.py (the one place that
    rule lives): integrated sigma and error, events written / attempted;
  * the event file: event count, closing tag, distinct event weights;
  * the BEAMS as written (status-4 PDG ids) -- a neutron run must say 2112;
  * the REMNANT DIQUARK census.  Extracting a valence quark from a proton
    leaves uu (2203) or ud (2101/2103); from a neutron dd (1103) or ud.  A
    neutron beam with a proton's remnant, or the reverse, shows up here;
  * the PDF ACTUALLY USED, event by event: GenPdfInfo's xf for the struck
    parton against LHAPDF, evaluated for the target the job CLAIMS (proton
    set isospin-swapped for 2112).  A double swap -- neutron beam plus the
    neutron set -- fails this check by the size of u/d;
  * the region: Q2, y and W2 from the scattered lepton in the TARGET REST
    FRAME (W2 = M_P^2 + 2 P.q - Q2, the analysis convention, P from the
    event's own target record), the fractions outside Q2 > 4 / W2 > 9, the
    smallest W2 generated, and the generator-cut margin it implies;
  * the error census from the .log: discarded events by reason, warnings.

Usage of --pdfset: the set the PDF check evaluates, default the proton set;
the check applies ThePEG's own isospin rule for a 2112 beam.
"""
import glob
import math
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "herwig7"))
sys.path.insert(0, os.path.join(REPO, "analysis"))
import target  # noqa: E402,F401  -- puts data/pdfs/lhapdf on LHAPDF_DATA_PATH
from herwig_xsec import read_job  # noqa: E402

M_P = 0.938272
DIQUARKS = {1103: "dd1", 2101: "ud0", 2103: "ud1", 2203: "uu1"}


def mdot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]


def events(path, maxev):
    """Yield (weights, pdfinfo, particles) with particles = [(id, st, p4)]."""
    ev = None
    n = 0
    with open(path) as f:
        for line in f:
            c = line[:2]
            if c == "E ":
                if ev is not None:
                    yield ev
                    n += 1
                    if maxev and n >= maxev:
                        return
                ev = {"w": None, "pdf": None, "p": []}
            elif ev is None:
                continue
            elif c == "W ":
                ev["w"] = tuple(float(x) for x in line.split()[1:])
            elif line.startswith("A 0 GenPdfInfo"):
                ev["pdf"] = line.split()[3:]
            elif c == "P ":
                t = line.split()
                ev["p"].append((int(t[3]), int(t[9]),
                                (float(t[7]), float(t[4]), float(t[5]),
                                 float(t[6]))))
    if ev is not None and (not maxev or n < maxev):
        yield ev


def inspect(d, maxev, pdfset):
    outs = sorted(glob.glob(os.path.join(d, "*-S*.out")))
    print(f"\n=== {os.path.relpath(d, REPO)}")
    if len(outs) != 1:
        print(f"  expected ONE .out, found {len(outs)}")
        if not outs:
            return
    j = read_job(outs[0])
    unit, sc = ("nb", 1e3) if j.sigma_pb > 1e3 else ("pb", 1.0)
    print(f"  sigma (integrated) = {j.sigma_pb/sc:.4f} +- {j.err_pb/sc:.4f} "
          f"{unit}   written {j.n_written}  attempted {j.n_attempted}  "
          f"survival {j.survival:.4f}")

    hep = os.path.join(d, "events.hepmc")
    with open(hep, "rb") as f:
        f.seek(max(0, os.path.getsize(hep) - 200))
        closed = b"END_EVENT_LISTING" in f.read()
    print(f"  events.hepmc {os.path.getsize(hep)/1e6:.1f} MB, closing tag "
          f"{'present' if closed else 'MISSING'}, "
          f"{os.path.getsize(hep)/1e3/max(j.n_written,1):.2f} kB/event")

    import lhapdf
    lhapdf.setVerbosity(0)
    pdf = lhapdf.mkPDF(pdfset, 0)
    swap = {1: 2, 2: 1, -1: -2, -2: -1}

    beams, wts, diq = Counter(), Counter(), Counter()
    nev = nq2 = nw2 = ngap = 0
    w2min = 1e9
    ymax = 0.0
    pdfdev, npdf = 0.0, 0
    for ev in events(hep, maxev):
        nev += 1
        wts[ev["w"]] += 1
        st4 = [(pid, p) for pid, st, p in ev["p"] if st == 4]
        beams[tuple(sorted(pid for pid, _ in st4))] += 1
        lep_in = [p for pid, p in st4 if abs(pid) in (13, 14)]
        tgt = [(pid, p) for pid, p in st4 if abs(pid) in (2212, 2112)]
        for pid, st, p in ev["p"]:
            if abs(pid) in DIQUARKS:
                diq[abs(pid)] += 1
        # scattered lepton: the hardest final-state mu
        fs = [p for pid, st, p in ev["p"] if st == 1 and abs(pid) == 13]
        if lep_in and tgt and fs:
            k = lep_in[0]
            kp = max(fs, key=lambda v: v[0])
            P = tgt[0][1]
            q = tuple(a - b for a, b in zip(k, kp))
            Q2 = -mdot(q, q)
            y = mdot(P, q) / mdot(P, k)
            W2 = M_P**2 + 2.0*mdot(P, q) - Q2
            ymax = max(ymax, y)
            w2min = min(w2min, W2)
            nq2 += Q2 <= 4.0
            nw2 += W2 <= 9.0
            ngap += W2 <= 9.0 and Q2 > 4.0
        # PDF actually used: GenPdfInfo id1 id2 x1 x2 Q xf1 xf2
        if ev["pdf"] and tgt:
            i1, i2 = int(ev["pdf"][0]), int(ev["pdf"][1])
            x2, Q, xf2 = float(ev["pdf"][3]), float(ev["pdf"][4]), \
                float(ev["pdf"][6])
            if abs(i2) <= 5 and xf2 > 0:
                fl = swap.get(i2, i2) if tgt[0][0] == 2112 else i2
                ref = pdf.xfxQ(fl, x2, Q)
                if ref > 0:
                    pdfdev = max(pdfdev, abs(xf2/ref - 1.0))
                    npdf += 1

    print(f"  parsed {nev} events; beams {dict(beams)}")
    print(f"  distinct weights {len(wts)}: {list(wts.items())[:3]}")
    tot = sum(diq.values()) or 1
    print("  remnant diquarks " + ", ".join(
        f"{DIQUARKS[k]}({k}) {v/tot:.1%}" for k, v in sorted(diq.items())))
    print(f"  PDF check vs {pdfset}"
          f"{' (isospin-swapped for 2112)' if beams and 2112 in next(iter(beams)) else ''}: "
          f"max |xf/LHAPDF - 1| = {pdfdev:.2e} over {npdf} events")
    if nev:
        print(f"  region (target rest frame): Q2<=4 {nq2/nev:.4%}, "
              f"W2<=9 {nw2/nev:.4%} (of which Q2>4: {ngap/nev:.4%}), "
              f"min W2 {w2min:.3f}, max y {ymax:.5f}")

    log = outs[0][:-4] + ".log"
    txt = open(log).read() if os.path.exists(log) else ""
    disc = Counter()
    for m in re.finditer(r"\*\* An event exception[^\n]*\n([^\n]*)\n"
                         r"The event will be discarded", txt):
        disc[m.group(1).strip()[:70]] += 1
    warn = Counter(m.group(1).strip()[:70] for m in re.finditer(
        r"\* A warning exception[^\n]*\n([^\n]*)", txt))
    print(f"  discarded (as logged, first reports only): {dict(disc)}")
    print(f"  warnings (first reports only): {dict(warn)}")
    tail = re.search(r"The following exception classes were reported in this "
                     r"run:\n((?:.*\n)*)", txt)
    if tail:
        print("  exception classes: " + " | ".join(
            l.strip() for l in tail.group(1).splitlines() if l.strip()))
    other = [l for l in open(os.path.join(d, "herwig_run.log"))
             if re.search(r"error|abort|segmentation|fatal", l, re.I)] \
        if os.path.exists(os.path.join(d, "herwig_run.log")) else []
    print(f"  herwig_run.log error-like lines: {len(other)}")


def main():
    args = sys.argv[1:]
    maxev = 0
    pdfset = target.PDFSET_P
    if "--max-events" in args:
        i = args.index("--max-events")
        maxev = int(args[i + 1])
        del args[i:i + 2]
    if "--pdfset" in args:
        i = args.index("--pdfset")
        pdfset = args[i + 1]
        del args[i:i + 2]
    for d in args:
        inspect(os.path.abspath(d), maxev, pdfset)


if __name__ == "__main__":
    main()
