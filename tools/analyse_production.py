#!/usr/bin/env python3
"""Analyse the paper-plots samples, in dependency order, into results.

The whole chain from the six generator rows to the per-nucleon tungsten
results pp01 reads, written down ONCE (the same reason tools/resample.sh
exists for the earlier production):

  1. per nucleon (p, n), per energy:
       mu  analysis/analyze.py      v2_powheg, v2_sherpa
           herwig7/make_histos_nlo.py v2powheg, v2powhegneg
       nu  analysis/analyze_nu.py   v2_powheg_nu, v2_sherpa_nlo,
                                    v2_herwig_nlo, v2_herwig_nlo_neg
  2. Herwig's full NLO = positive - negative, per nucleon
       (herwig7/combine_nlo.py -> herwig_nlo[_powheg]_full_q4w3_<t>)
  3. tungsten per nucleon = (74 p + 110 n)/184
       (analysis/combine_target.py -> <gen>_q4w3_W)

A STEP IS SKIPPED ONLY IF ITS RESULT IS NEWER THAN EVERY INPUT FILE it reads
(event files for step 1, result JSONs for 2 and 3), so a sample that grows or
is regenerated is re-analysed, and an interrupted run resumes.  A sample that
is not complete yet (its generator's own OK marker or file count missing) is
reported and left alone rather than analysed half-written.

Usage: tools/analyse_production.py [--energies 400,1000] [--targets p,n]
                           [--only powheg,sherpa,herwig,genie] [-j 3] [--force]
                           [--region q4w3] [--tag charmfinal]

--region runs the chain under another selection nested in q4w3 (sample_layout.
BENCH_SELECTIONS): q4w5 for the charm figures, q4w3_faser_{s,e} and q4w5_faser_dimuon for the
FASER tiers (the dimuon tier is q4w5_faser_dimuon, a charm tag).  --tag appends a heavy-flavour tag mode to every sample key
(analyze.py's suffix: charmfinal is the benchmark's charm definition), and to
every result name.
"""
import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))
import beams                                                # noqa: E402
import sample_layout                                             # noqa: E402

PY = sys.executable


def at(key, e):
    return beams.at_energy(key, e)


REGION = sample_layout.REGION     # set from --region in main()


def res(current, key, t, e):
    d = "results" if current == "mu" else "results_nu"
    return f"{BASE}/{d}/histos_{at(sample_layout.result_key(key, t, REGION), e)}.json"


def job_dirs(parent, base):
    import re
    pat = re.compile(rf"^{re.escape(base)}_\d+$")
    return sorted(d for d in glob.glob(f"{parent}/{base}_*")
                  if os.path.isdir(d) and pat.match(os.path.basename(d)))


def sample_inputs(row, current, t, e):
    """(event files, complete?) for one sample."""
    if row in ("powheg", "powhegmc"):
        dirs = job_dirs(f"{BASE}/powheg",
                        at(sample_layout.powheg_job_base(current, t,
                                                    mc=row == "powhegmc"), e))
        files = [f"{d}/events.hepmc" for d in dirs]
        ok = bool(dirs) and all(os.path.exists(f"{d}/V2_OK") for d in dirs)
    elif row.startswith("genie"):
        # genie/production/genie_job.sh; complete = every job directory carries V2_OK
        # (the driver leaves a short job in place for inspection)
        key = GENIE_ROW_KEY[row]
        base = at(sample_layout.genie_job_base(current, key, t), e)
        dirs = job_dirs(f"{BASE}/genie", base)
        files = [f"{d}/events.hepmc" for d in dirs]
        ok = bool(dirs) and all(os.path.exists(f"{d}/V2_OK") for d in dirs)
    elif row.startswith("herwig"):
        neg = row.endswith("neg")
        dirs = job_dirs(f"{BASE}/herwig7",
                        at(sample_layout.herwig_job_base(current, t, neg), e))
        files = [f"{d}/events.hepmc" for d in dirs]
        name = sample_layout.herwig_runname(current, t, e, neg)
        ok = bool(dirs) and all(glob.glob(f"{d}/{name}-S*.out") for d in dirs)
    else:
        rd = sample_layout.sherpa_rundir(current, t, e)
        dirs = job_dirs(rd, "job")
        files = [f"{d}/evtfull" for d in dirs]
        # the Sherpa driver writes the event file as it goes; a job is done
        # when its log carries the end-of-run table
        ok = bool(dirs) and all(
            any("Nominal" in ln for ln in open(f"{d}/sherpa.log",
                                               errors="replace"))
            if os.path.exists(f"{d}/sherpa.log") else False for d in dirs)
    return [f for f in files if os.path.exists(f)], ok


STEP1 = [
    # (current, row, command, result key)
    ("mu", "powheg",    ["analysis/analyze.py", "v2_powheg"], "powheg"),
    ("mu", "sherpa",    ["analysis/analyze.py", "v2_sherpa"], "sherpa"),
    ("mu", "herwig",    ["herwig7/make_histos_nlo.py", "v2powheg"],
     "herwig_nlo_powheg"),
    ("mu", "herwigneg", ["herwig7/make_histos_nlo.py", "v2powhegneg"],
     "herwig_nlo_powheg_neg"),
    ("nu", "powheg",    ["analysis/analyze_nu.py", "v2_powheg_nu"],
     "powheg_nu"),
    ("nu", "sherpa",    ["analysis/analyze_nu.py", "v2_sherpa_nlo"],
     "sherpa_nlo"),
    ("nu", "herwig",    ["analysis/analyze_nu.py", "v2_herwig_nlo"],
     "herwig_nlo"),
    ("nu", "herwigneg", ["analysis/analyze_nu.py", "v2_herwig_nlo_neg"],
     "herwig_nlo_neg"),
    # POWHEG-V2mc, massive charm: a CHARM-ONLY sample, so it is run only with
    # --tag charmfinal (an untagged result would be its parton-level charm
    # rate under an inclusive-looking name)
    ("nu", "powhegmc",  ["analysis/analyze_nu.py", "v2_powheg_nu_mc"],
     "powheg_nu_mc"),
]
# GENIE (pp03/pp04): row name -> earlier result key, per current.  G18_02a rows
# exist only up to 1 TeV (its declared validity), so higher energies report
# "not complete" for them, which is correct.
GENIE_ROW_KEY = {"genie_grv": "genie", "genie_nnpdf": "genie_nnpdf",
                 "genie_lo": "genie_lo", "genie_hedis": "genie"}
STEP1 += [
    ("mu", "genie_grv",   ["analysis/analyze.py", "v2_genie"], "genie"),
    ("mu", "genie_nnpdf", ["analysis/analyze.py", "v2_genie_nnpdf"],
     "genie_nnpdf"),
    ("nu", "genie_lo",    ["analysis/analyze_nu.py", "v2_genie_lo"], "genie_lo"),
    ("nu", "genie_nnpdf", ["analysis/analyze_nu.py", "v2_genie_nnpdf"],
     "genie_nnpdf"),
    ("nu", "genie_hedis", ["analysis/analyze_nu.py", "v2_genie"], "genie"),
]
HERWIG_FULL = {"mu": ("herwig_nlo_powheg", "herwig_nlo_powheg_neg",
                      "herwig_nlo_powheg_full"),
               "nu": ("herwig_nlo", "herwig_nlo_neg", "herwig_nlo_full")}
TUNGSTEN = {"mu": ("powheg", "sherpa", "herwig_nlo_powheg_full",
                   "genie", "genie_nnpdf"),
            "nu": ("powheg_nu", "sherpa_nlo", "herwig_nlo_full",
                   "genie_lo", "genie_nnpdf", "genie", "powheg_nu_mc")}


def fresh(out, inputs, force):
    if force or not os.path.exists(out) or not inputs:
        return False
    return os.path.getmtime(out) > max(os.path.getmtime(f) for f in inputs)


def run(cmd, env, log):
    with open(log, "w") as fh:
        rc = subprocess.call([PY] + cmd, cwd=BASE, env=env, stdout=fh,
                             stderr=subprocess.STDOUT)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--energies", default=",".join(
        f"{e:g}" for e in beams.BENCH_ENERGIES))
    ap.add_argument("--targets", default="p,n")
    ap.add_argument("--only", default="powheg,sherpa,herwig,genie")
    ap.add_argument("-j", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--region", default=sample_layout.REGION,
                    choices=sample_layout.BENCH_SELECTIONS)
    ap.add_argument("--tag", default="", choices=("", "charmfinal"))
    a = ap.parse_args()
    global REGION
    REGION = a.region
    sfx = f"_{a.tag}" if a.tag else ""
    es = [float(x) for x in a.energies.split(",")]
    ts = a.targets.split(",")
    only = a.only.split(",")
    os.makedirs(f"{BASE}/logs/analyse_production", exist_ok=True)

    tasks, status = [], []
    for e in es:
        for t in ts:
            for cur, row, cmd, key in STEP1:
                if not any(row.startswith(o) for o in only):
                    continue
                if row == "powhegmc" and not sfx:
                    continue
                if sfx:
                    # analyze.py / analyze_nu.py read the tag as a key suffix;
                    # make_histos_nlo.py as a second argument
                    if cmd[0] == "herwig7/make_histos_nlo.py":
                        cmd = cmd + [a.tag]
                    else:
                        cmd = cmd[:-1] + [cmd[-1] + sfx]
                key = key + sfx
                files, ok = sample_inputs(row, cur, t, e)
                out = res(cur, key, t, e)
                tag = f"{key} {t} {e:g}"
                if not ok:
                    status.append(f"  not complete: {tag}")
                    continue
                if fresh(out, files, a.force):
                    status.append(f"  up to date:   {tag}")
                    continue
                env = dict(os.environ, BENCH_SELECTION=REGION,
                           BENCH_TARGET=t, BENCH_ENERGY=f"{e:g}")
                env.pop("BENCH_Q2MIN", None)
                log = (f"{BASE}/logs/analyse_production/{key}_{REGION}_{t}_"
                       f"{beams.Beams('mu', e).tag}.log")
                tasks.append((tag, cmd, env, log))
    print("\n".join(status))
    fails = []
    with ThreadPoolExecutor(max_workers=a.j) as ex:
        futs = {ex.submit(run, c, env, log): (tag, log)
                for tag, c, env, log in tasks}
        for f, (tag, log) in futs.items():
            rc = f.result()
            print(f"  {'ok  ' if rc == 0 else 'FAIL'} {tag}  ({os.path.relpath(log, BASE)})")
            if rc:
                fails.append(tag)

    # step 2: Herwig positive - negative
    if any(o.startswith("herwig") for o in only):
        for e in es:
            for t in ts:
                for cur, (pk, nk, fk) in HERWIG_FULL.items():
                    pk, nk, fk = pk + sfx, nk + sfx, fk + sfx
                    p, n, f = (res(cur, k, t, e) for k in (pk, nk, fk))
                    if not (os.path.exists(p) and os.path.exists(n)):
                        continue
                    if fresh(f, [p, n], a.force):
                        continue
                    rc = subprocess.call([PY, "herwig7/combine_nlo.py",
                                          p, n, f], cwd=BASE)
                    if rc:
                        fails.append(f"combine_nlo {fk} {t} {e:g}")
    # step 3: tungsten per nucleon
    for e in es:
        for cur, keys in TUNGSTEN.items():
            for k in keys:
                if not any(k.startswith(o) for o in only):
                    continue
                k = k + sfx
                p, n, w = (res(cur, k, t, e) for t in ("p", "n", "W"))
                if not (os.path.exists(p) and os.path.exists(n)):
                    continue
                if fresh(w, [p, n], a.force):
                    continue
                rc = subprocess.call([PY, "analysis/combine_target.py",
                                      k, cur, f"{e:g}", "--region", REGION],
                                     cwd=BASE)
                if rc:
                    fails.append(f"combine_target {k} {e:g}")
    if fails:
        print("FAILED:\n  " + "\n  ".join(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
