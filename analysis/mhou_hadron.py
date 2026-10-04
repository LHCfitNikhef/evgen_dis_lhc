#!/usr/bin/env python3
"""The 7-point scale band on a HADRON-LEVEL distribution, bin by bin.

WHAT THIS ADDS THAT NOTHING ELSE HERE HAS (user, 2026-09-08: "in pp_09, the
MHOUs from POWHEG+Pythia8 are missing, can you please add them?" and, of the
FASER comparison, "an important aspect of NLO calculations is a proper
estimate of MHOUs, so we should not forget about this").

`analysis/powheg_nlo_uncertainty.py` gives POWHEG scale bands on sigma_fid and
REFUSES the FASERnu selection, correctly: it reads the Les Houches file, an
LHE has no hadrons, and that selection cuts on the charged-hadron multiplicity
and on Delta phi against the summed hadron system.  `analysis/mhou_diff.py`
gives per-bin bands on the ANALYTIC reference, which reaches the leptonic
kinematics and nothing else.  Neither can put a band on a track multiplicity.

>>> SO THE WEIGHTS ARE MARRIED TO THE SHOWERED EVENTS BY `lhe_index`. <<<
POWHEG's `storeinfo_rwgt 1` lets the matrix element be re-evaluated at another
scale after the fact, and powheg/reweight/ appends the seven points to the LHE
as extra weights; main_powheg writes the Les Houches event number into the
HepMC as `A 0 lhe_index N`, and analyze.analyze() now carries it through the
selection.  Joining the two gives every SELECTED event its seven weights, and
the band is then the envelope of seven histograms of the same events.

THE JOIN IS CHECKED, NOT ASSUMED.  A positional join that has slipped produces
perfectly plausible numbers attached to the wrong events, which is how the
first attempt at this join failed once already (memory:
lhe-index-and-wrong-lhe).  Here the nominal weight of the reweighted LHE must
reproduce the HepMC weight of the event it is joined to, up to one overall
constant, for every event -- and if it does not, this refuses rather than
returning a band.

WHY ONLY THE BASE SAMPLE IS REWEIGHTED, AND WHY THAT IS ENOUGH.  The neutrino
statistics were raised by adding independent BATCHES in their own run
directories (powheg/add_v2_batches.sh), and only the base directory's LHE is
reweighted here.  That costs nothing that matters: the band is a RATIO of the
same events at two scales, so the statistical fluctuation cancels between
numerator and denominator, and a relative band measured on 150k events is
applied to a central value measured on 980k.  The number of events behind the
band is recorded in the output.

Usage:
  BENCH_SELECTION=faser_e analysis/mhou_hadron.py --current nu
      -> results_nu/mhou_hadron_faser_e.json
"""
import glob
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import beams                                                 # noqa: E402
import selection as selection_mod                            # noqa: E402

CURRENT = (sys.argv[sys.argv.index("--current") + 1]
           if "--current" in sys.argv else "nu")
if CURRENT not in ("mu", "nu"):
    sys.exit("--current takes mu or nu")

# The nominal point of rwl_scale.xml: xiR = xiF = 1.  Everything is measured
# relative to it, so it has to be named rather than assumed to be first.
NOMINAL_ID = "1001"

def join_scale_weights(w, idx, W, ids, where):
    """Join showered events (weights w, Les Houches numbers idx) to the seven
    scale weights W of their own LHE, and return (F, worst): F[k] is each
    event's weight ratio to the nominal at scale point k.

    THE JOIN IS CHECKED, NOT ASSUMED (see the module docstring).  The nominal
    weight must reproduce the shower's up to one constant to 1e-4 -- the LHE
    writes six significant digits, so a correct join closes to ~3e-5 -- and
    the SAME test on the join shifted by one event must FAIL, since
    POWHEG-RES weights are nearly all equal up to sign and the tolerance
    alone would not expose a slipped join.  Either failure refuses.
    """
    inom = ids.index(NOMINAL_ID)
    if (len(idx) != len(w) or (idx < 0).any()
            or (len(idx) and idx.max() >= W.shape[0])):
        sys.exit(f"{where}: lhe_index cannot be joined to its weights")
    wn = W[idx, inom]
    good = wn != 0
    ratio = w[good] / wn[good]
    c = float(np.median(ratio)) if len(ratio) else 1.0
    worst = float(np.max(np.abs(ratio / c - 1.0))) if len(ratio) else 0.0
    if worst > 1e-4:
        sys.exit(f"{where}: the lhe_index join does not close ({worst:.3g}); refusing")
    js = np.clip(idx + 1, 0, W.shape[0] - 1)
    ws = W[js, inom]
    gs = ws != 0
    dsh = np.abs((w[gs] / ws[gs]) / c - 1.0)
    if len(dsh) and float((dsh > 1e-4).mean()) < 0.005:
        sys.exit(f"{where}: a join shifted by one event also closes -- "
                 f"the closure test cannot tell a slipped join; refusing")
    F = np.stack([np.where(wn != 0, W[idx, k] / np.where(wn != 0, wn, 1.0), 1.0)
                  if ids[k] != NOMINAL_ID else np.ones(len(w))
                  for k in range(len(ids))])
    return F, worst


def main():
    """PAPER PLOTS: the band on the  1 TeV samples, tungsten per nucleon.

        BENCH_SELECTION=q4w3_faser_e analysis/mhou_hadron.py --current nu
            -> results_nu/mhou_hadron_q4w3_faser_e_W.json

    The band is measured on EACH ROW'S OWN SAMPLE -- POWHEG-RES on the muon
    current (reweighted seed by seed, powheg/production/reweight_scale.sh, the first
    ten seeds), POWHEG-V2 on the neutrino current (its base batch) -- so the earlier production's
    "measured on the POWHEG-V2 cross-variant" caveat is gone.

    PER JOB, THEN PER NUCLEON, THEN TUNGSTEN.  Each job's showered events are
    joined to that job's own weights by lhe_index, with the earlier production's closure test (the
    nominal weight reproduces the showered weight up to one constant, to
    1e-6) applied job by job.  A nucleon's seven histograms are the sums over
    its jobs, then scaled so the nominal integral equals that nucleon's
    published sigma_fid in this selection (histos_<key>_<sel>_<t>.json) --
    a scale common to all seven curves, so every ratio is untouched and only
    the relative weight of p and n in the sum is set.  Tungsten is then
    (74 p + 110 n)/184 of the scaled curves, the sums of squares with the
    squared factors, so analysis/mhou_band.py reads it unchanged.
    """
    import re
    sel = selection_mod.get()
    if CURRENT == "nu":
        import analyze_nu as an
        key, resdir = "powheg_nu", "results_nu"
    else:
        import analyze as an
        key, resdir = "powheg", "results"
    Z, A = 74, 184
    per = {}
    for t in ("p", "n"):
        wdir = os.path.join(BASE, resdir, "mhou_weights")
        npzs = sorted(glob.glob(f"{wdir}/scale_{CURRENT}_{t}_job*.npz"),
                      key=lambda x: int(re.search(r"job(\d+)", x).group(1)))
        if not npzs:
            sys.exit(f"no scale weights for {CURRENT} {t} in {wdir} -- run "
                     f"powheg/production/reweight_scale.sh {CURRENT} {t}")
        acc, ids0, n_ev, worst_all, jobs = None, None, 0, 0.0, []
        for wf in npzs:
            n = int(re.search(r"job(\d+)", wf).group(1))
            job = f"{BASE}/powheg/v2{CURRENT}_{t}_job_{n}"
            z = np.load(wf)
            ids = [str(x) for x in z["ids"]]
            if ids0 is None:
                ids0 = ids
            elif ids != ids0:
                sys.exit(f"{wf}: weight ids differ from the first job's")
            W = z["weights"]
            inom = ids.index(NOMINAL_ID)
            obs, w, _np, _sw = an.analyze([f"{job}/events.hepmc"])
            w = np.asarray(w, dtype=float)
            idx = np.asarray(obs["lhe_index"], dtype=int)
            if len(idx) != len(w) or (idx < 0).any() or (len(idx) and idx.max() >= W.shape[0]):
                sys.exit(f"{job}: lhe_index cannot be joined to {wf}")
            F, worst = join_scale_weights(w, idx, W, ids, job)
            worst_all = max(worst_all, worst)
            h = {}
            for obsk, edges in an.BINS.items():
                v = np.asarray(obs[obsk], dtype=float)
                if len(v) != len(w):
                    continue
                h[obsk] = {
                    "curves": np.stack([np.histogram(v, edges, weights=w * F[k])[0] for k in range(len(ids))]),
                    "sw2": np.histogram(v, edges, weights=w * w)[0],
                    "sw2f": np.stack([np.histogram(v, edges, weights=w * w * F[k])[0] for k in range(len(ids))]),
                    "sw2f2": np.stack([np.histogram(v, edges, weights=w * w * F[k] ** 2)[0] for k in range(len(ids))]),
                    "edges": np.asarray(edges)}
            tot = np.array([(w * F[k]).sum() for k in range(len(ids))])
            if acc is None:
                acc = {"h": h, "tot": tot}
            else:
                for obsk in acc["h"]:
                    for q in ("curves", "sw2", "sw2f", "sw2f2"):
                        acc["h"][obsk][q] = acc["h"][obsk][q] + h[obsk][q]
                acc["tot"] = acc["tot"] + tot
            n_ev += len(w)
            jobs.append(os.path.relpath(job, BASE))
        rp = f"{BASE}/{resdir}/histos_{sample_layout_key(key, t, sel.name)}.json"
        with open(rp) as f:
            sig = json.load(f)["sigma_fid_pb"]
        scale = sig / acc["tot"][ids0.index(NOMINAL_ID)]
        per[t] = {"acc": acc, "scale": scale, "ids": ids0, "n": n_ev,
                  "worst": worst_all, "jobs": jobs, "sigma_fid_pb": sig}
        print(f"  {CURRENT} {t}: {len(jobs)} jobs, {n_ev} selected events, "
              f"join closes to {worst_all:.2g}, sigma_fid {sig:.6g} pb")
    ids = per["p"]["ids"]
    if per["n"]["ids"] != ids:
        sys.exit("p and n carry different weight ids")
    out = {"what": "7-point scale band on hadron-level distributions: "
                   "POWHEG weights joined per job to the showered events by "
                   "lhe_index, tungsten per nucleon (74 p + 110 n)/184",
           "current": CURRENT, "generator": {"mu": "POWHEG-RES", "nu": "POWHEG-V2"}[CURRENT],
           "sample_key": key, "selection": sel.as_dict(), "target": "W",
           "energy_gev": 1000.0, "nominal_id": NOMINAL_ID, "ids": ids,
           "n_events_band": {t: per[t]["n"] for t in per},
           "jobs": {t: per[t]["jobs"] for t in per},
           "join_closure": max(per[t]["worst"] for t in per), "hists": {}}
    fz, fn_ = Z / A, (A - Z) / A
    wt = {"p": fz * per["p"]["scale"], "n": fn_ * per["n"]["scale"]}
    inom = ids.index(NOMINAL_ID)
    for obsk in per["p"]["acc"]["h"]:
        hp, hn = per["p"]["acc"]["h"][obsk], per["n"]["acc"]["h"][obsk]
        curves = wt["p"] * hp["curves"] + wt["n"] * hn["curves"]
        sw2 = wt["p"] ** 2 * hp["sw2"] + wt["n"] ** 2 * hn["sw2"]
        sw2f = wt["p"] ** 2 * hp["sw2f"] + wt["n"] ** 2 * hn["sw2f"]
        sw2f2 = wt["p"] ** 2 * hp["sw2f2"] + wt["n"] ** 2 * hn["sw2f2"]
        base = curves[inom]
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(base != 0, curves / np.where(base != 0, base, 1.0), 1.0) - 1.0
        others = [k for k in range(len(ids)) if k != inom]
        out["hists"][obsk] = {"edges": hp["edges"].tolist(),
                              "rel_hi": np.max(np.vstack([r[others], np.zeros(len(base))]), axis=0).tolist(),
                              "rel_lo": np.min(np.vstack([r[others], np.zeros(len(base))]), axis=0).tolist(),
                              "curves": curves.tolist(), "sw2": sw2.tolist(),
                              "sw2f": sw2f.tolist(), "sw2f2": sw2f2.tolist()}
    tot = wt["p"] * per["p"]["acc"]["tot"] + wt["n"] * per["n"]["acc"]["tot"]
    envs = [tot[k] / tot[inom] - 1.0 for k in range(len(ids)) if k != inom]
    out["sigma_rel_hi"], out["sigma_rel_lo"] = max(envs), min(envs)
    p = os.path.join(BASE, resdir, f"mhou_hadron_{sel.name}_W.json")
    with open(p, "w") as f:
        json.dump(out, f)
    print(f"{out['generator']}, {sel.name}, tungsten: integrated band "
          f"+{100*out['sigma_rel_hi']:.2f}/{100*out['sigma_rel_lo']:.2f}%")
    print(f"wrote {p}")


def sample_layout_key(key, t, region):
    import sample_layout
    return sample_layout.result_key(key, t, region)


if __name__ == "__main__":
    main()
