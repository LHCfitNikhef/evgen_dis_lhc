#!/usr/bin/env python3
"""Copy the first N events of a Les Houches file, header and footer intact.

WHY THIS EXISTS.  The 400 GeV POWHEG-V2 sample was raised to 1,000,000 events
so the DIMUON tier would have the statistics to be plotted -- the dimuon
selection keeps a few events in ten thousand.  The PDF reweighting does not
need that sample size and cannot afford it: measured here, pwhg_main
reweights about 5,900 events a minute, so one pass over 1M events takes ~170
minutes and the 54 passes come to SIX AND A HALF DAYS, against about seven
hours at the 50,000 events the other two energies use.

REWEIGHTING A SUBSET IS SOUND, and it is worth being explicit about why.  The
PDF band is a spread ACROSS MEMBERS evaluated on the SAME events, so the
per-event Monte Carlo error very largely cancels in it: what is wanted is the
ratio member/central, not either cross-section on its own.  50,000 events is
also exactly what the 400 GeV band is being compared against at the other two
energies, so capping here makes the three energies like-for-like rather than
introducing an asymmetry.

WHAT IT DOES NOT DO.  It does not touch the FULL sample, which stays as it is
for the cross-section and the dimuon tier; only the copy handed to the
reweighter is capped.  Any analysis joining these weights back to showered
events must therefore drop events whose LHE index is beyond the cap -- that is
a real restriction, not a rounding, and the cap is written to a sidecar JSON
so the join can read it rather than assume it.

Usage: truncate_lhe.py <in.lhe> <out.lhe> <n_events>
"""
import json
import os
import sys


def main():
    src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    kept = 0
    with open(src, errors="replace") as f, open(dst, "w") as g:
        # header: everything up to the first <event>
        for line in f:
            if line.lstrip().startswith("<event"):
                break
            g.write(line)
        else:
            sys.exit(f"{src}: no <event> found")
        # the first <event> line has already been consumed
        g.write("<event>\n" if line.strip() == "<event>" else line)
        in_ev = True
        for line in f:
            s = line.strip()
            if s.startswith("<event"):
                if kept >= n:
                    break
                in_ev = True
            if in_ev:
                g.write(line)
            if s.startswith("</event>"):
                kept += 1
                if kept >= n:
                    break
                in_ev = False
        g.write("</LesHouchesEvents>\n")
    side = os.path.splitext(dst)[0] + ".cap.json"
    with open(side, "w") as f:
        json.dump({"source": os.path.basename(src), "events_kept": kept,
                   "cap_requested": n,
                   "note": "weights from this file cover LHE indices "
                           "0..events_kept-1 only; a join must drop the rest"},
                  f, indent=1)
    print(f"  wrote {dst}: {kept} events (cap {n}); manifest {side}")
    if kept < n:
        print(f"  note: source held only {kept} events, fewer than the cap")


if __name__ == "__main__":
    main()
