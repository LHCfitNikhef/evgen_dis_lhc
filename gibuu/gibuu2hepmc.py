#!/usr/bin/env python3
"""GiBUU's Les Houches output -> HepMC3, for analysis/analyze_nu.py.

GiBUU 2017 writes events in Les Houches, OSCAR 2013 or Shanghai format; there
is no HepMC.  Converting is the same choice made for NuWro and for the same
reason: every number in this benchmark goes through `parse_hepmc3`, and a
third event reader would put a whole generator behind code no gate watches.

>>> TWO THINGS ABOUT GiBUU'S LES HOUCHES FILE THAT A NAIVE READER GETS WRONG.
<<<

  1. THE OUTGOING LEPTON IS NOT IN THE PARTICLE BLOCK.  Only hadrons are
     listed.  The neutrino, the scattered lepton and the struck nucleon are on
     the `# 5` comment line that closes each event, written by
     LH_write_additionalInfo as

         '# 5 ', evtType, Weight, momLepIn, momLepOut, momNuc

     with each four-momentum in (E, px, py, pz) order.  An analysis that took
     the particle block as the final state would find no lepton at all and
     silently select nothing.

  2. THE WEIGHT IS IN UNITS OF 1e-38 cm^2, not pb and not mb
     (code/init/neutrino/neutrinoXsection.f90).  1 pb = 1e-36 cm^2, so the
     conversion is a factor 0.01 -- read out of the source rather than
     inferred from the numbers looking about right.

WHAT THE STRUCK NUCLEON'S CHARGE IS, AND WHY IT HAS TO BE SUPPLIED.
`NeutrinoProdInfo_Get` returns `chrg_nuc`, and GiBUU does NOT write it to the
file.  On a free-nucleon run we know it from the jobcard, so `--nucleon`
carries it through.  On a nucleus it is genuinely unavailable, and `--nucleon
mixed` writes a proton record while saying so on stderr: the four-momentum is
right either way -- which is what the DIS invariants are built from -- and
only the 0.14% mass difference is affected, the same approximation analyze.py
already documents for the GENIE ladders.

>>> THE OUTPUT IS VERIFIED BY THE CONSUMER, NOT BY EYE. <<<  The file is
written by hand (there is no HepMC3 binding in this Python), so it is then
read back with analyze_nu.parse_hepmc3 and the event count, the beams and the
summed weight are checked against the Les Houches file.  Hand-written formats
that "look right" are how this project has been bitten before.

>>> ONE LES HOUCHES FILE PER RUN, AND THE WEIGHTS ARE PER RUN. <<<  GiBUU
writes EventOutput.Pert.<run>.lhe for each of `num_runs_SameEnergy`, and the
weights in EACH file sum to the cross-section on their own.  Summing all the
files therefore gives R times sigma, not sigma -- a factor that would sail
straight through every downstream check, because a cross-section five times
too large is still a perfectly plausible number.  So every input file is
named on the command line, the weights are divided by how many there are, and
their individual sums are compared with one another as a cross-check.

Usage: gibuu2hepmc.py [--current nu|nubar] [--nucleon p|n|mixed]
                      --out <out.hepmc> <in1.lhe> [in2.lhe ...]
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

# 1 pb = 1e-36 cm^2 and GiBUU quotes 1e-38 cm^2.
W_TO_PB = 1.0e-2
NUCLEON_PID = {"p": 2212, "n": 2112, "mixed": 2212}


def events(path):
    """Yield (weight_pb, beam_nu, lep_out, nucleon, hadrons, evt_type).

    Four-momenta come out as (E, px, py, pz); the particle block is stored
    (px, py, pz, E, m) and is reordered here, once.
    """
    with open(path, errors="replace") as f:
        in_ev = False
        hadrons = []
        n_want = 0
        for line in f:
            s = line.strip()
            if s == "<event>":
                in_ev, hadrons, n_want = True, [], -1
                continue
            if not in_ev:
                continue
            if s.startswith("# 5"):
                c = s.split()
                # '# 5 <evtType> <Weight> <momLepIn x4> <momLepOut x4>
                #                                       <momNuc x4>'
                evt_type = int(c[2])
                w = float(c[3])
                lin = tuple(float(x) for x in c[4:8])
                lout = tuple(float(x) for x in c[8:12])
                nuc = tuple(float(x) for x in c[12:16])
                yield w * W_TO_PB, lin, lout, nuc, hadrons, evt_type
                in_ev = False
                continue
            if s.startswith("</event>"):
                in_ev = False
                continue
            c = s.split()
            if n_want < 0:
                n_want = int(c[0])          # NUP, the header line
                continue
            if len(c) >= 11:
                px, py, pz, e = (float(c[6]), float(c[7]),
                                 float(c[8]), float(c[9]))
                hadrons.append((int(c[0]), (e, px, py, pz), float(c[10])))


def main():
    argv = sys.argv[1:]
    nucleon, current = "p", "nu"
    for flag, default in (("--nucleon", "p"), ("--current", "nu")):
        if flag in argv:
            i = argv.index(flag)
            val = argv[i + 1]
            del argv[i:i + 2]
            if flag == "--nucleon":
                nucleon = val
            else:
                current = val
    if "--out" not in argv:
        sys.exit("usage: gibuu2hepmc.py [--current nu|nubar] "
                 "[--nucleon p|n|mixed] --out <out.hepmc> <in1.lhe> ...")
    i = argv.index("--out")
    dst = argv[i + 1]
    del argv[i:i + 2]
    srcs = argv
    if (not srcs or nucleon not in NUCLEON_PID
            or current not in ("nu", "nubar")):
        sys.exit("usage: gibuu2hepmc.py [--current nu|nubar] "
                 "[--nucleon p|n|mixed] --out <out.hepmc> <in1.lhe> ...")
    n_runs = len(srcs)
    # >>> THE CURRENT IS A FLAG AND NOT A GUESS. <<<  Nothing in the Les
    # Houches file says whether the beam was a neutrino or an antineutrino:
    # the outgoing lepton is written by four-momentum only, with no charge,
    # and GiBUU's process code 34 is "DIS" for both.  It has to come from the
    # jobcard, so the driver passes it -- and the beam and lepton PDG codes
    # below are the only place it enters.
    beam_pid, lep_pid = (14, 13) if current == "nu" else (-14, -13)
    if nucleon == "mixed":
        print("!! the struck nucleon's charge is not in the Les Houches file; "
              "writing a proton record.  The four-momentum is right, so the "
              "DIS invariants are; only the 0.14% mass difference is not.",
              file=sys.stderr)

    n, sum_w, types = 0, 0.0, {}
    per_file = []
    with open(dst, "w") as out:
        out.write("HepMC::Version 3.02.05\n")
        out.write("HepMC::Asciiv3-START_EVENT_LISTING\n")
        for src in srcs:
          s_file = 0.0
          for w_raw, lin, lout, nuc, hadrons, et in events(src):
            # >>> THE DIVISION BY THE NUMBER OF RUNS, and it is the whole
            # reason the files are named individually. <<<
            w = w_raw / n_runs
            s_file += w_raw
            parts = [(beam_pid, lin, 0.0, 4),
                     (NUCLEON_PID[nucleon], nuc, abs(nuc[0]), 4),
                     (lep_pid, lout, 0.10566, 1)]
            parts += [(pid, p, m, 1) for pid, p, m in hadrons]
            out.write(f"E {n} 1 {len(parts)}\n")
            out.write("U GEV MM\n")
            out.write(f"W {w:.17e}\n")
            pid_i = 0
            for pdg, p, m, st in parts:
                pid_i += 1
                par = 0 if st == 4 else -1
                out.write(f"P {pid_i} {par} {pdg} {p[1]:.10e} {p[2]:.10e} "
                          f"{p[3]:.10e} {p[0]:.10e} {m:.10e} {st}\n")
                if pid_i == 2:
                    out.write("V -1 0 [1,2]\n")
            n += 1
            sum_w += w
            types[et] = types.get(et, 0) + 1
          per_file.append(s_file)
        out.write("HepMC::Asciiv3-END_EVENT_LISTING\n")

    # EACH RUN IS AN INDEPENDENT ESTIMATE OF THE SAME CROSS-SECTION, so they
    # must agree; if they do not, the runs were not the same run condition and
    # dividing by their number is meaningless.
    if len(per_file) > 1:
        lo, hi = min(per_file), max(per_file)
        if lo <= 0 or hi / lo - 1.0 > 0.25:
            sys.exit(f"the per-run cross-sections span {lo:.4g} to {hi:.4g} "
                     f"(x{hi / max(lo, 1e-300):.2f}) -- these runs are not "
                     f"the same condition and must not be pooled")

    # ------------------------------------------------ read it back and check
    import analyze_nu as an
    an.NU_BEAM_PID = beam_pid
    n_back, w_back, seen_lep = 0, 0.0, 0
    for ww, k, P, ps, _d, _h in an.parse_hepmc3(dst):
        n_back += 1
        w_back += ww
        if k is None or P is None:
            sys.exit(f"event {n_back} lost a beam record in conversion")
        if any(pid == lep_pid for pid, _p in ps):
            seen_lep += 1
    if n_back != n:
        sys.exit(f"wrote {n} events, read back {n_back}")
    if abs(w_back / sum_w - 1.0) > 1e-9:
        sys.exit(f"weights do not survive the round trip: {sum_w} -> {w_back}")
    if seen_lep != n:
        sys.exit(f"{n - seen_lep} of {n} events have no outgoing lepton in "
                 f"the converted file -- the analysis would select nothing")

    js = (dst[:-len(".hepmc")] if dst.endswith(".hepmc") else dst) \
        + "_xsec.json"
    with open(js, "w") as f:
        json.dump({"generator": "gibuu", "sigma_pb": sum_w, "n_events": n,
                   "gibuu_event_types": types, "nucleon": nucleon,
                   "current": current}, f)
    print(f"{n} events, sigma = {sum_w:.6g} pb, GiBUU event types {types} "
          f"-> {dst}")


if __name__ == "__main__":
    main()
