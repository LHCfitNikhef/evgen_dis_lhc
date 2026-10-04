#!/usr/bin/env python3
"""Pull the <rwgt> weights out of a reweighted LHE into a compact array.

WHY THE WEIGHTS ARE HARVESTED AND THE LHE THROWN AWAY.  The PDF members have to
be reweighted in independent passes -- pwhg_main aborts once the file it reads
carries about twenty weights -- so a full set of 309 members means about thirty
output files of 60 MB each, per beam energy, all holding the SAME events and
differing only in which weights they carry.  Keeping them would be two
gigabytes of duplicated event records to store a few megabytes of numbers.

Each pass is therefore reduced to an (n_events, n_weights) array as soon as it
finishes, and merged at the end.

THE MERGE JOINS ON POSITION, and that is safe here for a reason worth stating:
every pass reads the same pristine pwgevents.lhe and POWHEG's reweighting is a
per-event operation that neither drops nor reorders events.  The merge CHECKS
that every pass reports the same event count and refuses otherwise, because a
positional join that has silently slipped produces perfectly plausible numbers
attached to the wrong events -- which is exactly how the first attempt at the
LHE-to-HepMC join failed.

Usage:
  harvest_weights.py <in.lhe> <out.npz>
  harvest_weights.py --merge <dir_of_npz> <out.npz>
"""
import os
import re
import sys

import numpy as np

WGT = re.compile(r"<wgt id='([^']+)'>\s*([-+0-9.eEdD]+)\s*</wgt>")


def harvest(lhe, out):
    ids, rows, cur = None, [], {}
    in_ev = False
    with open(lhe, errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("<event"):
                in_ev, cur = True, {}
            elif s.startswith("</event>"):
                if in_ev:
                    if ids is None:
                        ids = sorted(cur)
                    rows.append([cur.get(i, np.nan) for i in ids])
                in_ev = False
            elif in_ev:
                m = WGT.search(s)
                if m:
                    cur[m.group(1)] = float(
                        m.group(2).replace("D", "E").replace("d", "e"))
    if not rows:
        sys.exit(f"no weights found in {lhe}")
    arr = np.array(rows, dtype=np.float64)
    np.savez_compressed(out, ids=np.array(ids), weights=arr)
    print(f"    {os.path.basename(out)}: {arr.shape[0]} events x "
          f"{arr.shape[1]} weights")


def merge(d, out):
    files = sorted(f for f in os.listdir(d) if f.endswith(".npz"))
    if not files:
        sys.exit(f"no .npz files in {d}")
    ids, blocks, n = [], [], None
    for fn in files:
        z = np.load(os.path.join(d, fn), allow_pickle=False)
        w, i = z["weights"], [str(x) for x in z["ids"]]
        if n is None:
            n = w.shape[0]
        elif w.shape[0] != n:
            sys.exit(f"{fn} has {w.shape[0]} events, expected {n} -- the "
                     f"positional join would be misaligned; refusing")
        ids.extend(i)
        blocks.append(w)
    allw = np.concatenate(blocks, axis=1)
    # a duplicate id would mean two passes computed the same weight; keep the
    # first and say so rather than silently averaging or overwriting
    seen, keep = set(), []
    for k, i in enumerate(ids):
        if i in seen:
            continue
        seen.add(i)
        keep.append(k)
    if len(keep) != len(ids):
        print(f"  note: {len(ids)-len(keep)} duplicate weight id(s) dropped")
    np.savez_compressed(out, ids=np.array([ids[k] for k in keep]),
                        weights=allw[:, keep])
    print(f"  merged {len(files)} pass(es): {n} events x {len(keep)} weights")


def main():
    if sys.argv[1:2] == ["--merge"]:
        merge(sys.argv[2], sys.argv[3])
    else:
        harvest(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
