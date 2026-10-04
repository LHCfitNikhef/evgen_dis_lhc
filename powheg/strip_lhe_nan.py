#!/usr/bin/env python3
"""Drop LHE events carrying NaN or Infinity, loudly.

WHY THIS EXISTS (2026-08-26).  POWHEG-RES writes a NaN four-momentum into a
small fraction of neutrino CC events -- 3 events in 600000 on the first
production run.  That on its own would be harmless at 5e-6 of the sample.
What is not harmless is how Pythia reacts: its LHEF reader cannot parse the
token, and reports

    Abort from Pythia::next: reached end of Les Houches Events File

So a bad event does not skip, it ENDS THE JOB, and the shower stops there
claiming a clean end of file.  On that run three of ten jobs stopped early --
at events 10630, 31179 and 51557 of 60000 -- and 86819 events, 14.5% of the
sample, vanished with no error message anywhere in the chain.  The per-job
event counts in the shower log were the only visible symptom, and only
because they were printed side by side.

This is the CONVENTIONS.md rule 2 failure mode exactly: the run completed, every
exit status was zero, and the cross-section printed in the LHE header was
untouched, because the header is the integrator's number and does not know
how many events were read.

The published muon NC POWHEG-RES sample and both POWHEG-V2 samples are clean,
so this is specific to the RES code in its CC configuration -- but the guard
runs on every sample, since the cost of checking is one pass over a file that
is about to be read anyway.

    powheg/strip_lhe_nan.py [--check] <file.lhe> [...]

Rewrites each file in place with the offending events removed, and reports
what it dropped.  --check reports without writing.  A clean file is never
rewritten, so mtimes stay meaningful for analyze.py's input manifest.

Dropping is the right response rather than repairing: the events are
unphysical, there are a handful of them, and the analysis normalises by the
delivered sum of weights over the delivered event count, so removing an event
before the shower is exactly equivalent to never having generated it.

IT ALSO WRITES `<file>.nevents`, THE COUNT THE SHOWER WILL BE OFFERED.  That
number is not recoverable afterwards and the closure gate needs it.  Pythia
answers a `partonLevel failed; try again` on LHE input by DISCARDING that
event and reading the next one, so events vanish between the file and the
HepMC with nothing counting them -- 1187 of 50000 on the POWHEG-V2 muon NC
sample.  Divide the delivered sum of weights by the number of events that
came OUT and the loss cancels out of its own check; divide by the number that
went IN and it shows.  On the published POWHEG-V2 neutrino sample the right
denominator turns a reassuring 1.0008 into 0.9820 -- 1.8% of the weight is
genuinely missing and the gate never said so.  See analyze.delivered_sigma_pb.
"""
import re
import sys

BAD = re.compile(r"\b(NaN|nan|[+-]?Infinity|[+-]?inf)\b")


def strip(path, check=False):
    """Return (n_events, n_dropped, first_dropped_index)."""
    out, block, in_event = [], [], False
    n = dropped = 0
    first = None
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("<event"):
                in_event, block = True, [line]
                continue
            if in_event:
                block.append(line)
                if line.lstrip().startswith("</event>"):
                    n += 1
                    if any(BAD.search(b) for b in block):
                        dropped += 1
                        if first is None:
                            first = n
                    else:
                        out.extend(block)
                    in_event = False
                continue
            out.append(line)
    if in_event:                      # truncated final block: never keep it
        print(f"  {path}: file ends inside an <event> block -- tail discarded",
              file=sys.stderr)
    if dropped and not check:
        with open(path, "w") as f:
            f.writelines(out)
    if not check:
        # The count the shower will be offered.  Written even when nothing was
        # dropped -- a clean file still needs its denominator recorded.
        with open(path + ".nevents", "w") as f:
            f.write(f"{n - dropped}\n")
    return n, dropped, first


def main(argv):
    check = "--check" in argv
    paths = [a for a in argv if not a.startswith("-")]
    if not paths:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        print("usage: strip_lhe_nan.py [--check] <file.lhe> [...]",
              file=sys.stderr)
        return 2
    total_n = total_dropped = 0
    for p in paths:
        n, dropped, first = strip(p, check)
        total_n += n
        total_dropped += dropped
        if dropped:
            verb = "would drop" if check else "dropped"
            print(f"  {p}: {verb} {dropped}/{n} event(s), first at #{first}"
                  f" -- Pythia would have STOPPED THERE, losing {n - first + 1}")
    if total_dropped:
        print(f"strip_lhe_nan: {total_dropped}/{total_n} event(s) carried "
              f"NaN/Inf across {len(paths)} file(s)"
              f"{' (--check: nothing written)' if check else ''}")
    else:
        print(f"strip_lhe_nan: clean -- {total_n} events across "
              f"{len(paths)} file(s), nothing to drop")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
