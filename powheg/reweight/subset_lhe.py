#!/usr/bin/env python3
"""Extract a chosen set of events from a Les Houches file, by index.

WHY THIS IS EXACT AND NOT AN APPROXIMATION.  A PDF band on a selection is

    band over members of   S(member) = sum_{selected events} w_member(event)

divided by a normalisation that does NOT depend on the member (the nominal
sum, fixed once).  So every event that fails the selection contributes exactly
nothing to the member dependence, and reweighting it is wasted work.  Handing
the reweighter only the selected events therefore returns the SAME numbers,
not an estimate of them.

WHAT IT IS FOR.  The 400 GeV sample carries 1,000,000 events because the
dimuon tier keeps roughly one event in 12,000, and reweighting all of them
over 319 members was measured at 32 HOURS.  The events that actually enter the
dimuon band number about eighty.  Reweighting those takes seconds.

THE MANIFEST IS THE POINT.  A subset file loses the original indices, so the
row -> lhe_index map is written beside it as JSON and every consumer reads it
rather than assuming the order.  A positional join that has silently slipped
is the single failure mode this whole area keeps producing, so the mapping is
recorded rather than reconstructed.

Usage:
  subset_lhe.py <in.lhe> <out.lhe> <indices.json|comma list>
"""
import json
import os
import sys


def wanted_from(arg):
    if os.path.exists(arg):
        with open(arg) as f:
            d = json.load(f)
        idx = d["lhe_index"] if isinstance(d, dict) else d
    else:
        idx = [int(x) for x in arg.split(",") if x.strip()]
    return sorted(set(int(i) for i in idx))


def main():
    src, dst, spec = sys.argv[1], sys.argv[2], sys.argv[3]
    want = wanted_from(spec)
    want_set = set(want)
    if not want:
        sys.exit("no indices requested")
    kept = []
    with open(src, errors="replace") as f, open(dst, "w") as g:
        for line in f:                      # header, up to the first <event>
            if line.lstrip().startswith("<event"):
                break
            g.write(line)
        else:
            sys.exit(f"{src}: no <event> found")
        n = 0                               # index of the event being read
        buf = [line]
        in_ev = True
        for line in f:
            s = line.strip()
            if in_ev:
                buf.append(line)
                if s.startswith("</event>"):
                    if n in want_set:
                        g.writelines(buf)
                        kept.append(n)
                    n += 1
                    in_ev = False
                    buf = []
            elif s.startswith("<event"):
                in_ev, buf = True, [line]
        g.write("</LesHouchesEvents>\n")
    manifest = os.path.splitext(dst)[0] + ".index.json"
    with open(manifest, "w") as f:
        json.dump({"source": os.path.basename(src),
                   "rows": len(kept), "lhe_index": kept,
                   "note": "row k of any weight array harvested from this "
                           "file belongs to lhe_index[k] of the source"},
                  f)
    print(f"  wrote {dst}: {len(kept)} of {len(want)} requested events "
          f"(source held {n})")
    if len(kept) != len(want):
        miss = sorted(want_set - set(kept))[:5]
        print(f"  !! {len(want)-len(kept)} requested index/indices not found, "
              f"e.g. {miss} -- the source may be shorter than expected")
        sys.exit(1)
    print(f"  manifest {manifest}")


if __name__ == "__main__":
    main()
