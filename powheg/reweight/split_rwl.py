#!/usr/bin/env python3
"""Split a POWHEG rwl file into batches small enough for pwhg_main to survive.

pwhg_main aborts with a bare SIGABRT after roughly twenty PDF members in one
reweight file, so the members are added a few at a time and `rwl_add 1` makes
each pass append to what is already there.

THE FIRST BATCH CARRIES THE SCALE POINTS AND THE FOUR CENTRAL IDS, so a run
interrupted after one pass still produces exactly what the old central-only
file produced.  That ordering is deliberate: a partial result should be a
subset of the old one, not something new and half-formed.

Prints the number of batches on stdout, which is what the driver reads.

Usage: split_rwl.py <rwl_members.xml> <dest_dir> <weights_per_batch>
"""
import os
import re
import sys


def main():
    src, dest, per = sys.argv[1], sys.argv[2], int(sys.argv[3])
    text = open(src).read()
    weights = re.findall(r"<weight id='(\d+)'>(.*?)</weight>", text)
    if not weights:
        sys.exit(f"no <weight> entries in {src}")

    # scale points and the four central members first, then everything else
    def rank(wid):
        return 0 if wid.startswith(("1", "2")) else 1
    weights.sort(key=lambda w: (rank(w[0]), int(w[0])))

    batches = [weights[i:i + per] for i in range(0, len(weights), per)]
    for b, chunk in enumerate(batches, 1):
        lines = ["<initrwgt>",
                 f"<weightgroup name='batch{b}' combine='none'>"]
        for wid, body in chunk:
            lines.append(f"<weight id='{wid}'>{body}</weight>")
        lines += ["</weightgroup>", "</initrwgt>"]
        with open(os.path.join(dest, f"rwl_b{b:03d}.xml"), "w") as f:
            f.write("\n".join(lines) + "\n")
    print(len(batches))


if __name__ == "__main__":
    main()
