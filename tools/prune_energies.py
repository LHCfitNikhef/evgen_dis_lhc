#!/usr/bin/env python3
"""Delete the event files away from the 1 TeV anchor, once every result
that reads them is written (CONVENTIONS.md rule 1: intermediates do not outlive
their results; user, 2026-09-18: "make sure to afterwards delete all files
which are not needed").

The off-anchor samples exist for INTEGRATED results (the tracked
histos_*_<tag>.json) and, since 2026-09-18, for the per-event pion spectra
of the SIDIS study (results{,_nu}/pions_<key>_q4w3_{p,n}_<tag>.json).  Both
are tracked, so once the pion pass has read a sample its events are
intermediates again.  This script:

  1. lists every event file of every sample at an energy other than the
     anchor (POWHEG, Herwig, GENIE, Sherpa; both nucleons; the SIDIS ladder);
  2. REFUSES to delete a sample whose pion spectra are missing or OLDER
     than any of its event files -- the tripwire tools/analyse_and_prune.py
     uses, so a pass that never ran cannot be mistaken for one that did;
  3. deletes, and appends the record (file, events, size) to
     REGENERATE_ENERGIES.md, the recipe to bring them back
     (tools/regenerate_energies.sh).

The 1 TeV anchors are never touched (they are the differential samples and
the user decides about them).  Only event files go: LHE, grids, sidecars,
V2_OK markers, Herwig .out tables and Sherpa Results.zip stay, which is what
makes the regeneration cheap.

Usage: tools/prune_energies.py [--dry-run] [--only powheg,herwig,genie,sherpa]
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))
import beams                                                # noqa: E402
import faser_pions as fp                                    # noqa: E402
import paths                                                # noqa: E402
import sample_layout                                             # noqa: E402

# the SIDIS ladder, plus POWHEG-V2's points below it (2026-09-21,
# beams.SIDIS_ENERGIES_LOW; no other sample exists there)
ENERGIES = [e for e in beams.SIDIS_ENERGIES_LOW + beams.SIDIS_ENERGIES
            if e != beams.ANCHOR_ENERGY]
GENIE_ROWS = {"nu": ["genie_lo", "genie"], "mu": ["genie"]}     # the pion rows
PION_KEYS = {"nu": ["powheg_nu", "herwig_nlo_full", "sherpa_nlo", "genie_lo", "genie"],
             "mu": ["powheg", "herwig_nlo_powheg_full", "sherpa", "genie"]}


def job_dirs(parent, base):
    pat = re.compile(rf"^{re.escape(base)}_\d+$")
    return sorted(d for d in glob.glob(f"{parent}/{base}_*")
                  if os.path.isdir(d) and pat.match(os.path.basename(d)))


def samples():
    """(row, current, energy, target, key-for-the-pion-pass, [event files])"""
    out = []
    for cur in ("mu", "nu"):
        for e in ENERGIES:
            for t in ("p", "n"):
                at = lambda b: beams.at_energy(b, e)                # noqa: E731
                pk = "powheg_nu" if cur == "nu" else "powheg"
                fs = [f"{d}/events.hepmc" for d in
                      job_dirs(f"{BASE}/powheg", at(sample_layout.powheg_job_base(cur, t)))]
                out.append(("powheg", cur, e, t, pk, fs))
                hk = "herwig_nlo_full" if cur == "nu" else "herwig_nlo_powheg_full"
                fs = [f"{d}/events.hepmc" for neg in (False, True) for d in
                      job_dirs(f"{BASE}/herwig7", at(sample_layout.herwig_job_base(cur, t, neg)))]
                out.append(("herwig", cur, e, t, hk, fs))
                sk = "sherpa_nlo" if cur == "nu" else "sherpa"
                fs = [f"{d}/evtfull" for d in
                      job_dirs(sample_layout.sherpa_rundir(cur, t, e), "job")]
                out.append(("sherpa", cur, e, t, sk, fs))
                for gk in GENIE_ROWS[cur]:
                    fs = [f"{d}/events.hepmc" for d in
                          job_dirs(f"{BASE}/genie", at(sample_layout.genie_job_base(cur, gk, t)))]
                    out.append(("genie", cur, e, t, gk, fs))
                # the rows the pion study does not read (GENIE NNPDF, the
                # POWHEG-V2mc charm-only samples) have no per-event reader at
                # all: their integrated results are tracked, so they go too
                for gk in ("genie_nnpdf",):
                    fs = [f"{d}/events.hepmc" for d in
                          job_dirs(f"{BASE}/genie", at(sample_layout.genie_job_base(cur, gk, t)))]
                    out.append(("genie", cur, e, t, None, fs))
                if cur == "nu":
                    fs = [f"{d}/events.hepmc" for d in
                          job_dirs(f"{BASE}/powheg", at(sample_layout.powheg_job_base(cur, t, mc=True)))]
                    out.append(("powheg", cur, e, t, None, fs))
    return out


def n_events(f):
    if f.endswith("evtfull"):
        return None
    try:
        return int(subprocess.run(["grep", "-c", "^E ", f], capture_output=True,
                                  text=True).stdout.strip() or 0)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", default="powheg,herwig,genie,sherpa")
    a = ap.parse_args()
    rows = set(a.only.split(","))
    todo, refused, total = [], [], 0
    for row, cur, e, t, key, fs in samples():
        fs = [f for f in fs if os.path.exists(f) and os.path.getsize(f) > 0]
        if row not in rows or not fs:
            continue
        if key is not None:
            res = fp.spectra_path(cur, key, e, t)
            if (not os.path.exists(res)
                    or os.path.getmtime(res) < max(os.path.getmtime(f) for f in fs)):
                refused.append((row, cur, e, t, os.path.relpath(res, BASE)))
                continue
        for f in fs:
            todo.append((row, cur, e, t, f, os.path.getsize(f)))
            total += os.path.getsize(f)
    for r in refused:
        print("REFUSED (spectra missing or older than the events):", *r)
    print(f"{len(todo)} event file(s), {total / 1e9:.2f} GB"
          + (" would be deleted" if a.dry_run else " to delete"))
    if a.dry_run or not todo:
        return 1 if refused else 0
    rec = [f"\n## Deleted {time.strftime('%Y-%m-%d %H:%M')} by tools/prune_energies.py\n",
           "After the SIDIS pion passes (results{,_nu}/pions_*_q4w3_*.json) had read them. "
           "Regenerate: tools/regenerate_energies.sh (LHE, grids, sidecars kept).\n\n",
           "| file | events | GB |\n|---|---|---|\n"]
    for row, cur, e, t, f, sz in todo:
        n = n_events(f)
        rec.append(f"| `{os.path.relpath(f, BASE)}` | {n if n is not None else '-'} | {sz / 1e9:.2f} |\n")
        os.remove(f)
    with open(f"{BASE}/REGENERATE_ENERGIES.md", "a") as fh:
        fh.writelines(rec)
    print(f"deleted {len(todo)} file(s), {total / 1e9:.2f} GB; recorded in REGENERATE_ENERGIES.md")
    return 1 if refused else 0


if __name__ == "__main__":
    sys.exit(main())
