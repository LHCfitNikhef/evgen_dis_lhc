#!/usr/bin/env python3
"""Digest a gevgen log into rejection counts, then shrink it.

WHY.  gevgen prints every event record at its default message level, so a
125000-event job writes a 1 GB log that nothing reads.  What the non-DIS
study DOES need from it is the count of events GENIE generated, found
unphysical and redrew -- because that is the mechanism by which the
delivered channel mix departs from the spline cross-sections (the resonance
kinematics generator gives up on a (W,Q2) pair after 1001 tries at these
energies, about 10% of the time at 400 GeV, and the replacement is drawn
from the whole list).  So each job's log is reduced to `rejections.json`
(counts by kinematics generator) plus its last 200 lines, and the full log
is deleted (CONVENTIONS.md rule 1: intermediates do not outlive their results).

Usage: nondis_log_digest.py <jobdir> [...]
Idempotent: a job already digested is skipped.
"""
import json
import os
import re
import sys

PAT = {
    "rejected": re.compile(r"unphysical event is rejected"),
    "res_kinematics_failed": re.compile(r"RESKinematicsGenerator.*Could not select"),
    "qel_kinematics_failed": re.compile(r"QELKinematicsGenerator.*Couldn't"),
    "dis_kinematics_failed": re.compile(r"DISKinematicsGenerator.*Couldn't|DISKinematicsGenerator.*Could not"),
    "dfr_kinematics_failed": re.compile(r"DFRKinematicsGenerator.*Couldn't|DFRKinematicsGenerator.*Could not"),
    "other_iteration_warnings": re.compile(r"iterations"),
}


def digest(d):
    log = os.path.join(d, "genie.log")
    out = os.path.join(d, "rejections.json")
    if os.path.exists(out) and not os.path.exists(log):
        return
    if not os.path.exists(log):
        print(f"{d}: no genie.log")
        return
    n = {k: 0 for k in PAT}
    tail = []
    with open(log, errors="replace") as f:
        for line in f:
            for k, p in PAT.items():
                if p.search(line):
                    n[k] += 1
            tail.append(line)
            if len(tail) > 200:
                tail.pop(0)
    m = re.search(r"Number of events requested: (\d+)", "".join(tail))
    n["events_requested"] = int(m.group(1)) if m else None
    with open(out, "w") as f:
        json.dump(n, f, indent=1)
    with open(log + ".tail", "w") as f:
        f.writelines(tail)
    size = os.path.getsize(log)
    os.remove(log)
    print(f"{d}: {n['rejected']} rejected (RES {n['res_kinematics_failed']}, "
          f"QEL {n['qel_kinematics_failed']}); log {size/1e6:.0f} MB -> tail")


if __name__ == "__main__":
    for d in sys.argv[1:]:
        digest(d)
