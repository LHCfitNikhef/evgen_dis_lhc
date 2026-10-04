#!/usr/bin/env python3
"""Analyse a freshly generated sample, then delete the event files.

THE POINT.  122 GB of HepMC currently backs about 20 MB of JSON -- 6000:1.
Nothing downstream reads an event file once analyze.py has run, so on a
cluster the sample is an intermediate that should not outlive its own
analysis.  Keeping the pipeline generate -> analyse -> delete is what makes
"more energies, higher statistics" a CPU question rather than a disk one.

    tools/analyse_and_prune.py --gen pythia --energy 400
    tools/analyse_and_prune.py --gen pythia --energy 400 --dry-run
    tools/analyse_and_prune.py --gen pythia --energy 400 --keep-one

IT WILL NOT DELETE ANYTHING UNLESS THE ANALYSIS ACTUALLY SUCCEEDED.  Every
variant must have written a result whose n_parsed is non-zero and whose mtime
is newer than the sample it came from.  A deletion that races ahead of a
failed or stale analysis is unrecoverable, so the checks are the point of the
script rather than a nicety -- the same reasoning as the delivery-closure gate
in analyze.py.

--keep-one leaves the first job directory intact.  Worth it while a
configuration is still being debugged: one job is a few GB and can regenerate
the argument, where the whole sample cannot.

NOT for the 1 TeV anchor samples while the FASER selection is still outstanding:
a second fiducial region means re-parsing every sample, so those are still
live inputs.  The script refuses the anchor unless --allow-anchor is given.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "analysis"))

import beams  # noqa: E402

# generator -> (analysis script, sample dir, job-dir base, event filename,
#               the result variants to produce before deleting)
GENERATORS = {
    "pythia": ("analysis/analyze.py", "pythia8", "job", "events.hepmc",
               ["pythia", "pythia_charmfinal"]),
    "pythia_nu": ("analysis/analyze_nu.py", "pythia8", "nu_job", "events.hepmc",
                  ["pythia", "pythia_charmfinal"]),
    "genie": ("analysis/analyze.py", "genie", "job", "events.hepmc",
              ["genie", "genie_charmfinal"]),
    "powheg": ("analysis/analyze.py", "powheg", "job", "events.hepmc",
               ["powheg", "powheg_charmfinal"]),
    # Sherpa lives in the install tree and its event file is "evtfull"; the
    # run directory, not a job base, carries the energy.
    "sherpa_lo": ("analysis/analyze.py", None, "job", "evtfull",
                  ["sherpa_lo", "sherpa_lo_charmfinal"]),
}


def sample_dirs(gen, energy):
    """The job directories for this generator at this energy."""
    script, subdir, base, evt, _ = GENERATORS[gen]
    tagged = beams.at_energy(base, energy)
    if subdir is None:                      # Sherpa: energy is in the run dir
        import paths
        run = beams.at_energy("MuonDIS_LO", energy)
        parent = f"{paths.SHERPA_RUNS}/{run}"
        tagged = "job"
    else:
        parent = os.path.join(REPO, subdir)
    sys.path.insert(0, os.path.join(REPO, "analysis"))
    from analyze import job_files
    return [os.path.dirname(f) for f in job_files(parent, tagged, evt)], evt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", required=True, choices=sorted(GENERATORS))
    ap.add_argument("--energy", type=float, required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--keep-one", action="store_true",
                    help="leave the first job directory intact for debugging")
    ap.add_argument("--allow-anchor", action="store_true",
                    help="permit pruning the 1 TeV anchor (see the header)")
    a = ap.parse_args()

    if a.energy == beams.ANCHOR_ENERGY and not a.allow_anchor:
        sys.exit(f"refusing to prune the {beams.ANCHOR_ENERGY:g} GeV anchor: the "
                 f"FASER selection still has to re-parse those samples.\n"
                 f"Pass --allow-anchor once that is done.")

    script, _, _, _, variants = GENERATORS[a.gen]
    dirs, evt = sample_dirs(a.gen, a.energy)
    if not dirs:
        sys.exit(f"no {a.gen} sample found at {a.energy:g} GeV")
    total = sum(os.path.getsize(os.path.join(d, evt))
                for d in dirs if os.path.exists(os.path.join(d, evt)))
    print(f"{a.gen} at {a.energy:g} GeV: {len(dirs)} job dirs, "
          f"{total/1e9:.1f} GB of {evt}")

    env = dict(os.environ, BENCH_ENERGY=str(a.energy))
    oldest_sample = min(os.path.getmtime(os.path.join(d, evt)) for d in dirs
                        if os.path.exists(os.path.join(d, evt)))

    # --- analyse every variant, and require each to have really worked -----
    for v in variants:
        name = beams.at_energy(v, a.energy)
        out = os.path.join(REPO, "results_nu" if a.gen.endswith("_nu")
                           else "results", f"histos_{name}.json")
        print(f"  analysing {v} -> {os.path.basename(out)} ... ", end="", flush=True)
        if a.dry_run:
            print("(dry run)")
            continue
        r = subprocess.run([sys.executable, os.path.join(REPO, script), v],
                           cwd=REPO, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            print("FAILED")
            sys.exit("analysis failed; NOTHING deleted:\n" +
                     (r.stderr or r.stdout)[-800:])
        if not os.path.exists(out):
            sys.exit(f"analysis wrote no {out}; NOTHING deleted")
        d = json.load(open(out))
        if not d.get("n_parsed"):
            sys.exit(f"{out} has n_parsed={d.get('n_parsed')!r}; NOTHING deleted")
        if os.path.getmtime(out) < oldest_sample:
            sys.exit(f"{out} is OLDER than the sample it should describe; "
                     f"NOTHING deleted")
        print(f"ok ({d['n_parsed']} events)")

    # --- only now, delete -------------------------------------------------
    keep = {dirs[0]} if a.keep_one else set()
    freed = 0
    for d in dirs:
        f = os.path.join(d, evt)
        if d in keep or not os.path.exists(f):
            continue
        freed += os.path.getsize(f)
        print(f"  {'would delete' if a.dry_run else 'deleting'} {f}")
        if not a.dry_run:
            os.remove(f)
    print(f"\n{'would free' if a.dry_run else 'freed'} {freed/1e9:.1f} GB"
          + ("  (kept the first job dir)" if a.keep_one else ""))


if __name__ == "__main__":
    main()
