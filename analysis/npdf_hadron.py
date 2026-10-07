#!/usr/bin/env python3
"""Nuclear-PDF uncertainties on HADRON-LEVEL distributions under the FASERnu
selection, by POWHEG reweighting joined to the showered events.

>>> NEUTRINO CURRENT ONLY (user decision, 2026-10-07). <<<  The study is an
illustration of how nPDF uncertainties propagate to hadron-level observables,
and one current suffices for that; the user restricted it to nu_mu CC on
tungsten.  This is a deliberate carve-out from CONVENTIONS.md rule 2b (muon and
neutrino DIS kept in sync), like the single-inclusive pion study's, and the
figure (ppA2b) says so.  The muon side's leptonic nPDF impact is in Fig. C.1.

WHAT THIS ADDS.  Appendix C (analysis/npdf_impact.py, paper plot A2) gives the
nPDF impact on the inclusive leptonic distributions from YADISM, which has no
hadrons.  The FASERnu selection cuts on the charged-track multiplicity and on
Delta phi against the hadron system, so a band on N_ch or E_lead needs the
matrix-element weight of each SHOWERED event under every nPDF member.  This is
the construction of analysis/mhou_hadron.py with PDF members in place of scale
points: POWHEG re-evaluates each Les Houches event with another PDF
(`storeinfo_rwgt 1`), and the weights are married to the showered events by
`lhe_index`, with the join checked rather than assumed (mhou_hadron's
join_scale_weights, reused here).

>>> ONLY THE PROTON SAMPLE IS REWEIGHTED, AND THAT IS EXACT FOR THE AVERAGE
    NUCLEON. <<<  The nPDF grids (nNNPDF3.0 A184_Z74, EPPS21 W184,
nCTEQ15HQ 184_74) and their free baselines (*_W184free) are all the AVERAGE
nucleon of tungsten, (74 p + 110 n)/184 (memory nuclear-pdf-impact: the valence
sum rules give 1.40 and 1.60, checked).  The production's neutron samples are
run with `ih2 2`, which tells POWHEG to swap u <-> d in the grid it is given
-- so handing an average-nucleon grid to the NEUTRON sample would produce the
mirror nucleus (Z = 110), not tungsten.  The PROTON sample, ih2 = 1, used with
an average-nucleon grid IS a sample of the average nucleon: the cross-section
is linear in the PDF and every partonic channel is present in it.  The one
thing that stays proton-like is the beam remnant handed to the shower, which
enters the RATIOS drawn here (each nPDF over its own free baseline, both on
the same events) at no visible level.  The check is built in: the free-nucleon
NNPDF4.0 average (NNPDF40_nnlo_as_01180_W184free member 0) over the nominal
proton must reproduce the tungsten/proton ratio of the published
sigma_fid -- recorded in the output as `isospin_check`.

THE SUBSET.  Only selected events are reweighted (powheg/reweight/subset_lhe.py
explains why this is exact), and only a few jobs: every curve drawn is a ratio
of the SAME events under two PDFs, so the statistical fluctuation cancels
between numerator and denominator and the residual error is what the second
moments recorded here give.

Steps:
  analysis/npdf_hadron.py dump <job>     -> selected events of
        powheg/v2nu_p_job_<job>: lhe_index, weight and the four observables
        (results_nu/npdf_weights/events_nu_p_job<N>.npz) and the index list
        for subset_lhe.py (split into chunks indices_nu_p_job<N>c<k>.json to
        reweight in parallel; the production used the 29580 selected events
        with lhe_index < 66000 in six chunks)
  powheg/reweight/run_reweight_npdf.sh   (the POWHEG passes)
  analysis/npdf_hadron.py combine        -> results_nu/npdf_hadron_q4w3_faser_e.json
"""
import glob
import json
import os
import re
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

SELECTION = "q4w3_faser_e"
OBS = ["nch05", "Elead", "dphix", "Emu"]
NOMINAL_PDF = "NNPDF40_nnlo_as_01180"      # what the samples were generated with
BASELINE = "NNPDF40_nnlo_as_01180_W184free"
NUCLEAR = ["nNNPDF30_nlo_as_0118_A184_Z74", "EPPS21nlo_CT18Anlo_W184",
           "nCTEQ15HQ_FullNuc_184_74"]
OWN_BASELINE = {"nNNPDF30_nlo_as_0118_A184_Z74": "nNNPDF30_nlo_as_0118_p_W184free",
                "EPPS21nlo_CT18Anlo_W184": "CT18ANLO_W184free",
                "nCTEQ15HQ_FullNuc_184_74": "nCTEQ15HQ_1_1_W184free"}
RWL_JSON = os.path.join(BASE, "powheg", "reweight", "rwl_npdf_hadron.json")


CUR = "nu"
RESDIR = os.path.join(BASE, "results_nu")
WDIR = os.path.join(RESDIR, "npdf_weights")


def dump(job):
    os.environ["BENCH_SELECTION"] = SELECTION
    import analyze_nu as an
    d = f"{BASE}/powheg/v2{CUR}_p_job_{job}"
    obs, w, _np, _sw = an.analyze([f"{d}/events.hepmc"])
    w = np.asarray(w, dtype=float)
    idx = np.asarray(obs["lhe_index"], dtype=int)
    if len(idx) != len(w) or (idx < 0).any():
        sys.exit(f"{d}: lhe_index missing or short -- cannot join")
    os.makedirs(WDIR, exist_ok=True)
    arrs = {k: np.asarray(obs[k], dtype=float) for k in OBS}
    for k, v in arrs.items():
        if len(v) != len(w):
            sys.exit(f"{d}: observable {k} has {len(v)} entries for {len(w)} events")
    np.savez(f"{WDIR}/events_{CUR}_p_job{job}.npz", w=w, lhe_index=idx,
             edges=json.dumps({k: list(map(float, an.BINS[k])) for k in OBS}),
             **arrs)
    with open(f"{WDIR}/indices_{CUR}_p_job{job}.json", "w") as f:
        json.dump({"lhe_index": sorted(set(idx.tolist())), "selection": SELECTION,
                   "sample": os.path.relpath(d, BASE)}, f)
    print(f"{d}: {len(w)} selected events, lhe_index {idx.min()}..{idx.max()}")


def _band(name, vals):
    """(central, err_plus, err_minus) at 68% CL, from LHAPDF's own rule for the
    set's error type (replicas, or Hessian rescaled from 90% CL)."""
    import lhapdf
    lhapdf.setVerbosity(0)
    u = lhapdf.getPDFSet(name).uncertainty(list(vals), 68.268949)
    return u.central, u.errplus, u.errminus


def load(rwl):
    """[(job, observables, w, F, ids)] for every reweighted job, joined.

    Weight files are weights_nu_p_job<J>.npz or, for a job split into
    chunks, weights_nu_p_job<J>c<k>.npz, each with its subset manifest
    (.index.json: row k -> lhe_index).  Chunks are stacked, and only the
    showered events whose lhe_index was reweighted are kept -- the neutrino
    job is reweighted on its first 66k Les Houches events only."""
    from mhou_hadron import join_scale_weights
    files = glob.glob(f"{WDIR}/weights_{CUR}_p_job*.npz")
    by_job = {}
    for p in files:
        m = re.search(r"job(\d+)(c\d+)?\.npz$", p)
        by_job.setdefault(int(m.group(1)), []).append(p)
    if not by_job:
        sys.exit(f"no reweighted jobs for {CUR} in {WDIR}")
    out = []
    for j in sorted(by_job):
        Ws, idxs, ids0 = [], [], None
        for p in sorted(by_job[j]):
            z = np.load(p)
            man = json.load(open(p[:-4] + ".index.json"))
            ids = [str(x) for x in z["ids"]]
            if ids0 is None:
                ids0 = ids
            elif ids != ids0:
                sys.exit(f"{p}: weight ids differ from the job's first chunk")
            if z["weights"].shape[0] != man["rows"]:
                sys.exit(f"{p}: {z['weights'].shape[0]} rows, manifest says {man['rows']}")
            Ws.append(z["weights"]); idxs += man["lhe_index"]
        W = np.concatenate(Ws)
        row_of = {li: r for r, li in enumerate(idxs)}
        if len(row_of) != len(idxs):
            sys.exit(f"job {j}: an lhe_index is reweighted twice")
        ev = np.load(f"{WDIR}/events_{CUR}_p_job{j}.npz")
        keep = np.array([int(i) in row_of for i in ev["lhe_index"]])
        rows = np.array([row_of[int(i)] for i in ev["lhe_index"][keep]], dtype=int)
        w = ev["w"][keep]
        F, wst = join_scale_weights(w, rows, W, ids0, f"{CUR} job {j}",
                                    nominal_id=rwl["closure_id"])
        obs = {k: ev[k][keep] for k in OBS}
        edges = json.loads(str(ev["edges"]))
        out.append((j, obs, w, F, ids0, edges, wst))
    return out


def combine():
    rwl = json.load(open(RWL_JSON))
    data = load(rwl)
    ids = data[0][4]
    if any(d[4] != ids for d in data):
        sys.exit("jobs carry different weight ids")
    col = {i: k for k, i in enumerate(ids)}

    def members(name):
        return [col[i] for i in rwl["sets"][name]["weight_ids"]]

    ib = members(BASELINE)[0]
    ic = col[rwl["closure_id"]]
    edges = data[0][5]
    # the `lhapdf` keyword trap: PDF weights identical to the closure
    for d in data:
        if np.max(np.abs(d[3][ib] - 1.0)) < 1e-6:
            sys.exit(f"job {d[0]}: the PDF weights equal the closure -- the "
                     "`lhapdf` keyword trap (powheg/reweight/README.md)")
    # pooled arrays: every job is the same generator at the same energy with
    # flat-magnitude weights, so pooling the events IS accepted-count weighting
    v = {k: np.concatenate([d[1][k] for d in data]) for k in OBS}
    w = np.concatenate([d[2] for d in data])
    F = np.concatenate([d[3] for d in data], axis=1)
    wf = w[None, :] * F                         # weights x events
    n_ev = len(w)
    out = {"what": "nPDF bands on hadron-level distributions, POWHEG weights "
                   "joined to the showered events by lhe_index; proton sample "
                   "reweighted with average-nucleon grids = tungsten per nucleon",
           "current": CUR, "generator": "POWHEG-V2",
           "selection": SELECTION, "energy_gev": 1000.0,
           "jobs": [d[0] for d in data], "n_events": n_ev,
           "join_closure": max(d[6] for d in data), "obs": {}, "integrated": {}}
    tot = wf.sum(axis=1)
    # NORMALISATION: the free-nucleon NNPDF4.0 curve is scaled to the
    # published tungsten sigma_fid of this selection (p and n samples run
    # separately, histos_<key>_<sel>_W.json), a factor common to every curve,
    # so no ratio moves.  The ISOSPIN CHECK compares the reweighted
    # W184free/proton with the published W/p: the two agree at parton level by
    # construction, and differ at hadron level only through the neutron
    # remnant's effect on the track-multiplicity selection, which reweighting
    # the proton sample cannot see (it is common to every nPDF and its
    # baseline, so the drawn ratios are not affected).
    key = "powheg_nu"
    pub = {t: json.load(open(os.path.join(RESDIR, f"histos_{key}_{SELECTION}_{t}.json")))["sigma_fid_pb"]
           for t in ("p", "W")}
    norm = pub["W"] / tot[ib]
    wf = wf * norm
    tot = tot * norm
    out["normalisation"] = {"published_sigma_fid_pb": pub, "factor": float(norm)}
    out["isospin_check"] = {"W184free_over_proton": float(tot[ib] / tot[ic]),
                            "published_W_over_p": pub["W"] / pub["p"]}

    def ratio_err(k, a, b):
        """MC error of sum(a)/sum(b) per bin, a and b the SAME events under
        two weights: sqrt(sum (a - R b)^2) / sum b."""
        e = np.asarray(edges[k])
        A = np.histogram(v[k], e, weights=a)[0]
        B = np.histogram(v[k], e, weights=b)[0]
        with np.errstate(divide="ignore", invalid="ignore"):
            R = np.where(B != 0, A / B, 1.0)
            bi = np.clip(np.digitize(v[k], e) - 1, 0, len(e) - 2)
            d = a - R[bi] * b
            var = np.histogram(v[k], e, weights=d * d)[0]
            return np.where(B != 0, np.sqrt(var) / np.abs(B), np.nan)

    for k in OBS:
        e = np.asarray(edges[k])
        s = np.stack([np.histogram(v[k], e, weights=wf[i])[0] for i in range(len(ids))])
        s2 = np.histogram(v[k], e, weights=(w * norm) ** 2)[0]
        o = {"edges": edges[k], "baseline": s[ib].tolist(),
             "proton_nominal": s[ic].tolist(), "sumw2_nominal": s2.tolist(),
             "n_events": np.histogram(v[k], e)[0].tolist()}
        for name in NUCLEAR:
            own = members(OWN_BASELINE[name])[0]
            m = members(name)
            bands = [_band(name, s[m][:, b]) for b in range(s.shape[1])]
            cen, ep, em = (np.array([x[i] for x in bands]) for i in range(3))
            ob = s[own]
            with np.errstate(divide="ignore", invalid="ignore"):
                o[name] = {"central": cen.tolist(), "errplus": ep.tolist(),
                           "errminus": em.tolist(), "own_baseline": ob.tolist(),
                           "ratio": np.where(ob > 0, cen / ob, np.nan).tolist(),
                           "ratio_hi": np.where(ob > 0, (cen + ep) / ob, np.nan).tolist(),
                           "ratio_lo": np.where(ob > 0, (cen - em) / ob, np.nan).tolist(),
                           "to_nnpdf40": np.where(s[ib] > 0, cen / s[ib], np.nan).tolist(),
                           "ratio_mc_error": ratio_err(k, wf[m[0]], wf[own]).tolist()}
        out["obs"][k] = o
    for name in NUCLEAR:
        c, p, mm = _band(name, tot[members(name)])
        own = tot[members(OWN_BASELINE[name])[0]]
        out["integrated"][name] = {"to_own_baseline": c / own, "errplus": p / own,
                                   "errminus": mm / own,
                                   "to_nnpdf40_W184free": c / tot[ib]}
    fn = os.path.join(RESDIR, f"npdf_hadron_{SELECTION}.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print(f"{CUR}: {n_ev} events in {len(data)} jobs, join closes to "
          f"{out['join_closure']:.2g}, W184free/proton = "
          f"{out['isospin_check']['W184free_over_proton']:.4f} (published W/p "
          f"{out['isospin_check']['published_W_over_p']:.4f})")
    for name, x in out["integrated"].items():
        print(f"  {name:34s} / own baseline {x['to_own_baseline']:.4f} "
              f"+{x['errplus']:.4f} -{x['errminus']:.4f}   / NNPDF4.0 "
              f"{x['to_nnpdf40_W184free']:.4f}")
    print(f"wrote {fn}")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("dump", "combine"):
        sys.exit(__doc__)
    if sys.argv[1] == "dump":
        dump(int(sys.argv[2]))
    else:
        combine()
