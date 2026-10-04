#!/usr/bin/env python3
"""The antineutrino SIDIS pass (paper plot 15, user 2026-10-04): extract each
antineutrino sample's hadron spectra as soon as it is complete, verify the
result, DELETE the event file at once (CONVENTIONS.md rule 1, "don't wait until
the very end"), and combine p and n into tungsten when both are done.

    tools/nubar_sidis.py watch     # loop until every point is extracted
    tools/nubar_sidis.py status

A sample counts as complete when its producer marked it so:
  powheg_nubar      <ladder point>/SIDIS_SHOWER_OK   (tools/powheg_nubar_sidis.sh)
  herwig_nlo_nubar  both positive jobs extracted by tools/herwig_faser_ladder.sh
  sherpa_nlo_nubar  <run dir>/SIDIS_GEN_OK            (tools/sherpa_nubar_sidis.sh)
  genie_*_nubar     4 jobs with V2_OK                 (tools/genie_nubar_sidis.sh)
The event file is deleted only if the spectra JSON exists, parsed every
event, and is newer than the file.  Herwig keeps its faserdata_events.npz and
.out files (the pp12/12b ladder); POWHEG keeps its Les Houches file.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "analysis"))
import beams          # noqa: E402
import faser_pions as fp   # noqa: E402
import paths          # noqa: E402

ES = [300.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0]
KEYS = ["powheg_nubar", "herwig_nlo_nubar", "sherpa_nlo_nubar",
        "genie_lo_nubar", "genie_nubar"]
PY = os.environ.get("BENCH_PYTHON", sys.executable)


def points():
    for k in KEYS:
        for e in ES:
            if k == "genie_lo_nubar" and e > 1000.0:
                continue          # G18_02a's declared validity
            for t in ("p", "n"):
                yield k, e, t


def ready(k, e, t):
    if k == "herwig_nlo_nubar":
        tag = beams.Beams("nu", e).tag
        ds = [f"{BASE}/herwig7/flnubarpwg_{t}_job_{tag}_{j}" for j in (1, 2)]
        return all(os.path.exists(f"{d}/faserdata_events.npz")
                   and os.path.exists(f"{d}/events.hepmc") for d in ds)
    if k.startswith("genie"):
        return len(fp.nubar_files(k, e, t)) == 4
    return bool(fp.nubar_files(k, e, t))


def done(k, e, t):
    p = fp.spectra_path("nubar", k, e, t)
    return os.path.exists(p) and os.path.getsize(p) > 0


def extract_and_prune(k, e, t):
    files = fp.nubar_files(k, e, t)
    log = f"{BASE}/logs/pions_nubar_{k}_{e:g}_{t}.log"
    r = subprocess.run([PY, f"{BASE}/analysis/faser_pions.py", "extract",
                        "--current", "nubar", "--key", k, "--energy", f"{e:g}",
                        "--target", t], stdout=open(log, "w"),
                       stderr=subprocess.STDOUT)
    p = fp.spectra_path("nubar", k, e, t)
    if r.returncode != 0 or not os.path.exists(p):
        print(f"  EXTRACT FAILED {k} {e:g} {t}  (see {log})", flush=True)
        return False
    d = json.load(open(p))
    newest = max(os.path.getmtime(f) for f in files)
    if d.get("n_parsed", 0) <= 0 or os.path.getmtime(p) < newest:
        print(f"  NOT PRUNED {k} {e:g} {t}: suspicious result", flush=True)
        return False
    gb = sum(os.path.getsize(f) for f in files) / 1e9
    for f in files:
        os.remove(f)
    print(f"  done {k:18s} {e:6g} {t}: {d['n_parsed']} events, "
          f"{gb:.1f} GB of events deleted", flush=True)
    return True


def combine_ready():
    for k in KEYS:
        for e in ES:
            if k == "genie_lo_nubar" and e > 1000.0:
                continue
            w = fp.spectra_path("nubar", k, e, "W")
            if os.path.exists(w):
                continue
            if all(done(k, e, t) for t in ("p", "n")):
                r = subprocess.run([PY, f"{BASE}/analysis/faser_pions.py", "combine",
                                    "--current", "nubar", "--key", k,
                                    "--energy", f"{e:g}"], capture_output=True, text=True)
                print(f"  combined {k} {e:g}: "
                      f"{'ok' if r.returncode == 0 else 'FAILED ' + r.stdout[-300:] + r.stderr[-300:]}",
                      flush=True)


def status():
    pts = list(points())
    n_done = sum(done(*p) for p in pts)
    n_ready = sum((not done(*p)) and ready(*p) for p in pts)
    print(f"{n_done}/{len(pts)} extracted, {n_ready} ready and waiting")
    return n_done == len(pts)


def watch():
    while True:
        for k, e, t in points():
            if not done(k, e, t) and ready(k, e, t):
                extract_and_prune(k, e, t)
        combine_ready()
        if status():
            print("all antineutrino SIDIS points extracted and combined", flush=True)
            return
        time.sleep(120)


if __name__ == "__main__":
    {"watch": watch, "status": status}[sys.argv[1] if len(sys.argv) > 1 else "status"]()
