#!/usr/bin/env python3
"""Does every result still read the WHOLE sample it was analysed from?

WHY (2026-09-04).  Every result JSON carries an input manifest -- the exact
files, sizes and mtimes it was built from -- and `analyze.input_manifest`
warns when those files span more than a day, which catches a result built
from a MIX of productions.  It cannot catch a sample that GREW: when a
production is extended with new job directories, none of the files the old
result listed changes, so its manifest stays perfectly self-consistent while
the result quietly describes a third of the sample.  That is what a 2026-09-04
audit found: 24 GENIE tier results reading 4 of 12 job files (the sample was
tripled on 08-29, after the tier pass) and 6 POWHEG-RES 4 TeV tier results
reading 10 of 22.  Rule 1b's "more GENIE events" were on disk, unanalysed.

WHAT IT DOES.  For each result with a manifest, take the job directories its
files came from (<parent>/<base>_<N>/<file>, N an integer, the same rule as
`analyze.job_files`), count the matching directories on disk that hold the
same file, and complain when there are MORE on disk than in the manifest.
Fewer is not an error: samples are pruned once analysed (rule 1), and a
result whose inputs have vanished is the only record of them.

Exit status 1 if any result is stale, so a driver can gate on it.
tools/run_checks.sh runs it.

Usage: tools/check_manifests.py [-v]        (-v also lists the clean ones)
"""
import glob
import json
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def jobdirs_on_disk(parent, base, fname):
    """Directories <parent>/<base>_<N>/ holding <fname>, N a bare integer --
    the exact rule of analyze.job_files, so job_400GeV_N never counts as a
    job of the anchor."""
    pat = re.compile(rf"^{re.escape(base)}_\d+$")
    n = 0
    for d in glob.glob(f"{parent}/{base}_*"):
        if os.path.isdir(d) and pat.match(os.path.basename(d)) \
                and os.path.exists(os.path.join(d, fname)):
            n += 1
    return n


def check_one(path):
    """Return (status, message).  status: 'stale' | 'ok' | 'pruned' | 'skip'."""
    with open(path) as f:
        d = json.load(f)
    if not isinstance(d, dict):
        return "skip", "not a result record"
    man = d.get("inputs")
    if not isinstance(man, dict) or not man.get("files"):
        return "skip", "no manifest"
    # group the manifest's files by (parent, base, fname)
    groups = {}
    for e in man["files"]:
        rel = e["path"]
        full = os.path.normpath(os.path.join(BASE, rel))
        jobdir = os.path.dirname(full)
        fname = os.path.basename(full)
        m = re.match(r"^(.*)_(\d+)$", os.path.basename(jobdir))
        if not m:
            continue                     # not a job directory: single file
        key = (os.path.dirname(jobdir), m.group(1), fname)
        groups[key] = groups.get(key, 0) + 1
    if not groups:
        return "skip", "no job-directory inputs"
    stale, pruned, msgs = False, False, []
    # >>> A DELIBERATE SUBSET IS NOT A STALE RESULT, BUT IT IS STILL CHECKED.
    # The shower-arm baselines ("powheg[_nu]_sub", paper plots 9-11) read
    # the production jobs the QED and Herwig arms re-showered and no others,
    # so that each ratio is taken on the same events.  They must read EXACTLY
    # that subset (analysis/sample_layout.arm_jobnums): fewer or more is stale.
    sub = re.match(r"^histos_powheg(_nu)?_sub_", os.path.basename(path))
    if sub:
        sys.path.insert(0, os.path.join(BASE, "analysis"))
        import sample_layout
        cur = "nu" if sub.group(1) else "mu"
        t = d.get("target")
        # No arm jobs on disk (pruned, or a fresh clone without the samples)
        # is not a stale result: fall through to the generic comparison,
        # which reports it as pruned.
        try:
            want = len(sample_layout.arm_jobnums(cur, t)) if t in ("p", "n") else None
        except SystemExit:
            want = None
    for (parent, base, fname), n_man in sorted(groups.items()):
        if sub and want is not None:
            ok = n_man == want
            stale = stale or not ok
            msgs.append(f"{os.path.relpath(parent, BASE)}/{base}_N/{fname}: "
                        f"{n_man} = the arm subset of {want}" if ok else
                        f"{os.path.relpath(parent, BASE)}/{base}_N/{fname}: "
                        f"manifest {n_man}, arm subset {want}")
            continue
        n_disk = jobdirs_on_disk(parent, base, fname)
        where = os.path.relpath(parent, BASE)
        if n_disk > n_man:
            stale = True
            msgs.append(f"{where}/{base}_N/{fname}: manifest {n_man}, "
                        f"on disk {n_disk}")
        elif n_disk < n_man:
            pruned = True
            msgs.append(f"{where}/{base}_N/{fname}: manifest {n_man}, "
                        f"on disk {n_disk} (pruned)")
        else:
            msgs.append(f"{where}/{base}_N/{fname}: {n_man}")
    status = "stale" if stale else ("pruned" if pruned else "ok")
    return status, "; ".join(msgs)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    verbose = "-v" in argv
    files = sorted(glob.glob(f"{BASE}/results/*.json")
                   + glob.glob(f"{BASE}/results_nu/*.json"))
    counts = {"stale": 0, "ok": 0, "pruned": 0, "skip": 0}
    for path in files:
        try:
            status, msg = check_one(path)
        except (ValueError, OSError) as exc:
            status, msg = "skip", f"unreadable: {exc}"
        counts[status] += 1
        rel = os.path.relpath(path, BASE)
        if status == "stale":
            print(f"STALE   {rel}: {msg}")
        elif verbose and status != "skip":
            print(f"{status:7s} {rel}: {msg}")
    print(f"manifests: {counts['ok']} whole, {counts['pruned']} pruned, "
          f"{counts['stale']} STALE, {counts['skip']} without job inputs")
    if counts["stale"]:
        print("*** a result above reads FEWER job files than its sample has: "
              "re-run that analysis (tools/resample.sh) ***", file=sys.stderr)
    return 1 if counts["stale"] else 0


if __name__ == "__main__":
    sys.exit(main())
