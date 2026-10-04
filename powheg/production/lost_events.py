#!/usr/bin/env python3
"""Which LHE events did the Pythia shower NOT deliver, and does it bias sigma_fid?

    powheg/production/lost_events.py <sample_jobbase> [...]     e.g. v2nu_p_job_400GeV
    powheg/production/lost_events.py --jobs powheg/v2nu_p_job_400GeV_1 [...] [--dump-lost out.lhe]

--dump-lost writes the lost LHE events of the listed jobs to one file, so
they can be re-showered with a candidate fix.

The normalisation is

    sigma_fid = sigma_int x sum(w, delivered & fiducial) / sum(w, delivered)

which is unbiased only if the events Pythia throws away ("ProcessContainer::
constructProcess: setting mass failed", "parton+hadronLevel failed; giving
up") have the same fiducial fraction as the ones it keeps.  This measures it
directly, at LHE level so that the shower's own migration does not enter:

  * lost events = LHE indices absent from the HepMC `lhe_index` attributes
    (main_powheg stamps info.nSelected() - 1, the 0-based position of the LHE
    event just read -- the join analysis/mhou_hadron.py relies on);
  * the JOIN IS CHECKED, not assumed: for the first --check-n delivered
    events the Q2 of the hard-process lepton in the HepMC must equal the LHE
    Q2 for the overwhelming majority (a permuted join would disagree for
    essentially all of them);
  * for delivered, lost and offered (= all) events: the signed-weight
    fraction inside q4w3 (Q2 > 4, W2 = M_P^2 + y 2kP - Q2 > 9, y = q.P/k.P with
    the target light-like along -z, as in inspect_sample.py), and the
    heavy-flavour content (any c or b parton, incoming or outgoing);
  * the bias of the delivered-event estimator relative to the offered-event
    one: f_delivered / f_offered - 1, with its binomial error.

--w2min X moves the W2 edge of the fiducial test (default 9, i.e. q4w3), to ask
whether a harder W cut removes the heavy-flavour loss (W > 5 GeV: --w2min 25).

Only jobs with a V2_OK marker are used.  Nothing is written.
"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_sample as I  # noqa: E402

M_P = I.M_P
POWHEG = os.path.dirname(HERE)
W2MIN = 9.0                   # fiducial W2 edge; --w2min overrides


def lhe_path(jobdir):
    for line in open(os.path.join(jobdir, "shower.log")):
        if "Beams:LHEF" in line:
            return line.split("|")[2].strip()
    raise SystemExit(f"no Beams:LHEF in {jobdir}/shower.log")


def delivered_indices(hepmc):
    idx = []
    with open(hepmc) as f:
        for line in f:
            if line.startswith("A 0 lhe_index"):
                idx.append(int(line.split()[3]))
    return idx


def lhe_events(path):
    """Yield (two_kP, weight, Q2, y, W2, heavy) per event, in file order."""
    two_kP = None
    with open(path) as f:
        in_init = in_ev = False
        for line in f:
            if line.startswith("<init"):
                in_init, first = True, True
                continue
            if in_init:
                if first:
                    p = line.split()
                    e1, e2 = float(p[2]), float(p[3])
                    if e2 < 1.0:                       # POWHEG-RES fixed target
                        two_kP = 2.0 * e1 * M_P
                    else:                              # POWHEG-V2, sqrt(s)/2 each
                        two_kP = 4.0 * e1 * e2 - M_P ** 2
                    first = False
                if line.startswith("</init"):
                    in_init = False
                continue
            if line.startswith("<event"):
                in_ev, head, k, kp, heavy = True, True, None, None, ""
                nout = 0
                continue
            if line.startswith("</event"):
                in_ev = False
                Q2, y, _ = I.q2_y(k, kp)
                W2 = M_P ** 2 + y * two_kP - Q2
                yield two_kP, w, Q2, y, W2, (heavy, nout)
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
                if abs(pid) in I.LEPTONS:
                    if st == -1:
                        k = mom
                    elif st == 1:
                        kp = mom
                else:
                    if st == 1:
                        nout += 1
                    if abs(pid) in (4, 5):
                        # classify: an OUTGOING c or b wins over an incoming one
                        tag = ("out" if st == 1 else "in") + ("c" if abs(pid) == 4 else "b")
                        if not heavy or tag.startswith("out"):
                            heavy = tag if not heavy.startswith("out") else heavy


def join_check(hepmc, lhe_q2, n):
    """(agree, agree_if_shifted, total): |dQ2/Q2| < 1% between the HepMC
    scattered lepton and the LHE event at lhe_index, and at lhe_index + 1.

    1%, not tighter: the shower legitimately moves Q2 for the ~5-10% of
    events with a re-massed heavy quark, and POWHEG-RES's LHE lepton energy
    (1000.469 at 1 TeV) is not Pythia's beam record's.  The shifted join is
    the control -- a permuted join is what this must be able to exclude."""
    good = bad = tot = 0
    for ev in I.hepmc_events(hepmc, n):
        r = I.event_kinematics(ev)
        i = ev["lhe_index"]
        if r is None or i is None or i + 1 >= len(lhe_q2):
            continue
        tot += 1
        good += abs(r[0] - lhe_q2[i]) <= 1e-2 * lhe_q2[i]
        bad += abs(r[0] - lhe_q2[i + 1]) <= 1e-2 * lhe_q2[i + 1]
    return good, bad, tot


class Acc:
    def __init__(self):
        from collections import Counter
        self.cat = Counter()
        self.n = 0
        self.sw = self.swf = 0.0
        self.swf_c = 0.0          # fiducial AND an outgoing charm quark
        self.swf_b = 0.0          # fiducial AND an outgoing bottom quark
        self.n_heavy = 0
        self.lowW = self.lowQ = 0
        self.q2 = []
        self.w2 = []

    def add(self, w, Q2, W2, hv, keep=False):
        heavy, nout = hv
        self.cat[(heavy or "light") + f"/{nout}p"] += 1
        self.n += 1
        self.sw += w
        fid = Q2 > 4.0 and W2 > W2MIN
        if fid:
            self.swf += w
            if heavy == "outc":
                self.swf_c += w
            elif heavy == "outb":
                self.swf_b += w
        if heavy:
            self.n_heavy += 1
        self.lowW += W2 <= 9.0
        self.lowQ += Q2 <= 4.0
        if keep:
            self.q2.append(Q2)
            self.w2.append(W2)

    @property
    def f(self):
        return self.swf / self.sw if self.sw else float("nan")


def med(a):
    a = sorted(a)
    return a[len(a) // 2] if a else float("nan")


def analyse(jobs, check_n, dump=None):
    dl, lo = Acc(), Acc()
    joins = []
    for j in jobs:
        lhe = lhe_path(j)
        hep = os.path.join(j, "events.hepmc")
        idx = delivered_indices(hep)
        have = set(idx)
        if len(have) != len(idx):
            raise SystemExit(f"{j}: duplicate lhe_index values -- join unusable")
        q2s = []
        n_off = 0
        for i, (_, w, Q2, y, W2, heavy) in enumerate(lhe_events(lhe)):
            n_off += 1
            q2s.append(Q2)
            if i in have:
                dl.add(w, Q2, W2, heavy)
            else:
                lo.add(w, Q2, W2, heavy, keep=True)
        if max(idx) >= n_off:
            raise SystemExit(f"{j}: lhe_index {max(idx)} beyond {n_off} LHE events")
        nside = open(lhe + ".nevents").read().strip() if os.path.exists(lhe + ".nevents") else "?"
        g, b, t = join_check(hep, q2s, check_n)
        joins.append((os.path.basename(j), n_off, nside, len(idx), n_off - len(idx), g, b, t))
        if dump is not None:
            dump_lost(lhe, have, dump)
    return dl, lo, joins


def report(label, jobs, check_n, dump=None):
    dl, lo, joins = analyse(jobs, check_n, dump)
    off_sw = dl.sw + lo.sw
    f_off = (dl.swf + lo.swf) / off_sw
    f_del, f_lost = dl.f, lo.f
    bias = f_del / f_off - 1.0
    # binomial error on the bias: the lost set is a subsample of size n_lost
    p = lo.n / (dl.n + lo.n)
    # signed weights can put f_lost outside [0, 1] once the cut removes most of
    # the lost set; the binomial error is then only indicative
    err = p * math.sqrt(abs(f_lost * (1 - f_lost)) / max(lo.n, 1)) / f_off if lo.n else 0.0
    print(f"== {label}: {len(jobs)} job(s)")
    for name, n_off, nside, nd, nl, g, b, t in joins:
        print(f"   {name}: offered {n_off} (sidecar {nside}), delivered {nd}, lost {nl}; "
              f"join check (|dQ2|<1%): {g}/{t}, shifted-join control {b}/{t}")
    print(f"   lost: {lo.n} events = {100*lo.n/(dl.n+lo.n):.3f}% of offered, "
          f"{100*lo.sw/off_sw:.3f}% of the signed weight")
    print(f"   q4w3 fraction (LHE level, signed weights): delivered {f_del:.5f}, "
          f"lost {f_lost:.5f}, offered {f_off:.5f}")
    print(f"   BIAS of the delivered-event estimator: {100*bias:+.3f}% +- {100*err:.3f}%")
    print(f"   heavy flavour (c/b anywhere): delivered {100*dl.n_heavy/dl.n:.1f}%, "
          f"lost {100*lo.n_heavy/max(lo.n,1):.1f}%")
    for name, acc in (("delivered", dl), ("lost", lo)):
        tops = ", ".join(f"{k} {100*v/max(acc.n,1):.1f}%" for k, v in acc.cat.most_common(7))
        print(f"   {name} by flavour/outgoing partons: {tops}")
    print(f"   Q2 <= 4: delivered {100*dl.lowQ/dl.n:.1f}%, lost {100*lo.lowQ/max(lo.n,1):.1f}%;  "
          f"W2 <= 9: delivered {100*dl.lowW/dl.n:.1f}%, lost {100*lo.lowW/max(lo.n,1):.1f}%;  "
          f"lost medians Q2 {med(lo.q2):.2f}, W2 {med(lo.w2):.1f}")
    # the same estimator restricted to an outgoing c (b) inside q4w3 -- the
    # parton-level stand-in for a charm-tagged fiducial cross-section
    for q, a_d, a_l in (("c", dl.swf_c, lo.swf_c), ("b", dl.swf_b, lo.swf_b)):
        if a_d + a_l == 0:
            continue
        bq = (a_d / dl.sw) / ((a_d + a_l) / off_sw) - 1.0
        print(f"   q4w3 with an outgoing {q} quark: {100*a_l/(a_d+a_l):.2f}% of its weight "
              f"lost; delivered-estimator bias {100*bq:+.2f}%")
    return bias, err


def dump_lost(lhe, have, out):
    """Append the LOST events of one LHE (header of the first file) to out."""
    first = not os.path.exists(out) or os.path.getsize(out) == 0
    with open(lhe) as f, open(out, "a") as o:
        i = -1
        buf = None
        for line in f:
            if line.startswith("<event"):
                i += 1
                buf = [line]
                continue
            if buf is not None:
                buf.append(line)
                if line.startswith("</event"):
                    if i not in have:
                        o.writelines(buf)
                    buf = None
                continue
            if first and not line.startswith("</LesHouchesEvents"):
                o.write(line)


def main():
    a = sys.argv[1:]
    check_n = 5000
    if "--check-n" in a:
        i = a.index("--check-n")
        check_n = int(a[i + 1])
        del a[i:i + 2]
    global W2MIN
    if "--w2min" in a:
        i = a.index("--w2min")
        W2MIN = float(a[i + 1])
        del a[i:i + 2]
        print(f"fiducial edge W2 > {W2MIN:g} GeV2 (the 'q4w3' labels below mean Q2 > 4 and this)")
    dump = None
    if "--dump-lost" in a:
        i = a.index("--dump-lost")
        dump = a[i + 1]
        del a[i:i + 2]
    if a and a[0] == "--jobs":
        report("explicit jobs", a[1:], check_n, dump)
        if dump:
            open(dump, "a").write("</LesHouchesEvents>\n")
        return
    for base in a:
        jobs = sorted(d for d in glob.glob(os.path.join(POWHEG, base + "_*"))
                      if re.fullmatch(r"\d+", d[len(os.path.join(POWHEG, base)) + 1:])
                      and os.path.exists(os.path.join(d, "V2_OK")))
        if not jobs:
            print(f"== {base}: no complete jobs")
            continue
        report(base, jobs, check_n)


if __name__ == "__main__":
    main()
