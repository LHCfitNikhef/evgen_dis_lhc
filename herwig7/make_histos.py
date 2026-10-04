#!/usr/bin/env python3
"""Build results/histos_herwig.json from the Herwig job_*/events.hepmc
samples, in exactly the schema written by analysis/analyze.py.

Normalisation (unweighted events, like the pythia branch): the cross-section
of the generation region is Herwig's INTEGRATED total (the sampler's
"Total integrated xsec", which includes events later vetoed during
showering), combined over jobs by inverse variance; then
    dsigma = (N_bin / N_written) * sigma_int / width.
The generation cuts equal the fiducial cuts (Q2 > 4 GeV2, 0.2 < y < 0.9),
which analyze.py re-applies on the reconstructed muon.

SELECTIONS.  analyze() applies the active BENCH_SELECTION itself, so obs/w
already carry only the surviving events; this script has to do two things
with that.  The result name goes through analyze.result_path(), so a FASER
tier cannot overwrite the inclusive result -- the overwrite trap that energy
tagging had to fix twice already.  And the normalisation needs NO special
case: fiducial_sigma() scales the integrator total by n_fid/n_written, so
unlike the WEIGHTED branches of analyze.py it was never pinned to the
integrator and a narrower selection is handled by construction.  (That
pinning bug is why Sherpa's Tier S first came out at its inclusive value to
the digit -- see analyze.py and CONVENTIONS.md rule 2.)

CORRECTED 2026-08-20: this script previously used the "Total (from generated
events)" line, which is the same integral scaled DOWN by the shower/remnant
veto survival rate (~1.2%).  Because the generation cuts are the fiducial
cuts, the vetoed events were fiducial events, so that convention reported
sigma_fid low by the veto rate -- it is what made the hadron-level sample
(0.9861 x YADISM LO) disagree with the veto-free ME-level sample (0.9975).
The rule and the evidence live in herwig_xsec.py; analyze_nu.py had already
applied it on the neutrino side.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "analysis"))
from analyze import (BINS, D_HIST_KEYS, SELECTION, analyze,  # noqa: E402
                     dmeson_summary, dphi_summary, input_manifest,
                     result_path, fmt_sigma, flag_low_stats)
import selection  # noqa: E402
from hardproc import hard_partons  # noqa: E402
from herwig_xsec import combine, fiducial_sigma, read_jobs  # noqa: E402
import sys as _sys, os as _os  # noqa: E402
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "analysis"))
import beams  # noqa: E402
from analyze import ENERGY, at_energy, job_files  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp


HARD_TAGS = {"charm": 4, "bottom": 5}


def main():
    # optional channel tag: `make_histos.py charm` selects the heavy-quark-
    # INITIATED channel (gamma*/Z + Q -> Q off the heavy-quark PDF), tagged on
    # the hard process before shower and hadronisation.  Herwig writes no
    # Pythia 21/23 status codes, so analysis/hardproc.py locates the hard
    # vertex from the graph instead -- the one with the exchanged boson coming
    # in.  NOT photon-gluon fusion, which is O(alpha_s) and absent at LO.
    #
    # Since 2026-08-21 the tag may carry a MODE suffix, matching
    # analysis/analyze.py: `charmfinal` is the benchmark definition (a charm
    # quark in the FINAL state after shower and hadronisation, i.e. a charm
    # hadron), `charmin` is the charm-INITIATED one, `charmany` either side,
    # and a bare `charm` the original out-of-the-hard-process tag.
    tag = sys.argv[1].lower() if len(sys.argv) > 1 else None
    mode = "out"
    if tag is not None and tag not in HARD_TAGS:
        for suf in ("final", "in", "any"):
            if tag.endswith(suf) and tag[:-len(suf)] in HARD_TAGS:
                tag, mode = tag[:-len(suf)], suf
                break
    if tag is not None and tag not in HARD_TAGS:
        sys.exit(f"unknown tag '{sys.argv[1]}' (use "
                 f"{'|'.join(HARD_TAGS)}, optionally + final|in|any)")
    hard_flavour = HARD_TAGS[tag] if tag else None
    tag_stats = {} if tag else None
    genkey = "herwig" + (f"_{tag}{'' if mode == 'out' else mode}"
                         if tag else "")

    # Per-energy job dirs and run name.  The anchor keeps "job" and
    # "DIS-mu"; other energies are job_<tag>_N holding DIS-mu-<tag>-S*.out,
    # because Herwig's saverun names the run after the card.
    jbase = at_energy("job")
    runname = "DIS-mu" if ENERGY == beams.ANCHOR_ENERGY else \
        f"DIS-mu-{beams.Beams('mu', ENERGY).tag}"
    files = job_files(f"{HERE}", jbase, "events.hepmc")
    if not files:
        sys.exit(f"no herwig7/{jbase}_*/events.hepmc files found")
    jobs = read_jobs(f"{HERE}/{jbase}_*/{runname}-S*.out")
    if len(jobs) != len(files):
        sys.exit(f"MISMATCH: {len(files)} events.hepmc files but {len(jobs)} "
                 f".out files -- a job produced events without a cross-section "
                 f"or vice versa.  Refusing to normalise a mixed sample.")
    sigma_sample_pb, sigma_err_pb, n_gen_tot, survival = combine(jobs)
    print(f"herwig: sigma_gen_region = {sigma_sample_pb/1e3:.4f} "
          f"+- {sigma_err_pb/1e3:.4f} nb over {len(jobs)} jobs, "
          f"survival {survival:.4f} ({100*(1-survival):.2f}% vetoed in "
          f"shower/hadronisation and regenerated)")

    obs, w, n_events, sum_w_all = analyze(
        files, hard_flavour=hard_flavour, tag_stats=tag_stats,
        hard_tag_mode=mode)
    print(f"herwig: parsed {n_events} events, {len(w)} in fiducial region")
    if n_events != n_gen_tot:
        print(f"WARNING: parsed {n_events} != generated {n_gen_tot}")

    label = "Herwig 7.3.0 LO" + (f" ({tag})" if tag else "")
    out = {"generator": genkey, "label": label,
           "hard_tag_mode": mode,
           "n_parsed": n_events, "n_fiducial": len(w),
           "sigma_gen_region_pb": sigma_sample_pb,
           "sigma_gen_region_err_pb": sigma_err_pb,
           "veto_survival": survival}
    out["inputs"] = input_manifest(
        files, [j.sigma_pb for j in jobs],
        [j.n_written for j in jobs], " [herwig]")
    # unweighted: dsigma = (N_bin / N_written) * sigma_int / width
    norm = sigma_sample_pb / n_gen_tot
    out["sigma_fid_pb"], out["sigma_fid_err_pb"] = fiducial_sigma(
        sigma_sample_pb, sigma_err_pb, len(w), n_gen_tot)
    if tag_stats is not None:
        tot = sum(tag_stats["in_flavours"].values())
        print("herwig: hard-process incoming flavours: " + ", ".join(
            f"|{k}|: {100*v/tot:.2f}%"
            for k, v in sorted(tag_stats["in_flavours"].items(),
                               key=lambda kv: int(kv[0]))))
        print(f"herwig: {tag_stats['n_no_hard_record']} events with no hard "
              f"record; {tag_stats['n_tagged']} tagged, of which "
              f"{tag_stats['n_in_out_mismatch']} WITHOUT the same flavour "
              f"outgoing (at LO the t-channel conserves flavour, so under the "
              f"hard-process tag this must be 0; under the final-state tag it "
              f"is not checked)")
        print(f"herwig: tag comparison on the same sample -- hard-out "
              f"{tag_stats.get('n_tagged_out', 0)}, final-state "
              f"{tag_stats.get('n_tagged_final', 0)} of {n_events} parsed")
        out["tag_stats"] = tag_stats
        out["hard_flavour"] = hard_flavour
    if SELECTION.name != selection.DEFAULT:
        # parity with analyze.py, which records this for every narrower-than-
        # generation selection; it is the efficiency the report divides by.
        out["selected_weight_fraction"] = len(w) / n_events
    print(f"herwig: sigma_fid = {fmt_sigma(out['sigma_fid_pb'])} "
          f"+- {fmt_sigma(out['sigma_fid_err_pb'])} "
          f"({len(w)} events selected)")
    flag_low_stats(out, "herwig")

    if not len(w):
        # An EMPTY selection is a legitimate result, not a crash: the means
        # and the D-meson summaries below all divide by sum(w).
        out["hists"] = {}
        os.makedirs(f"{BASE}/results", exist_ok=True)
        ofn = result_path(genkey)
        with open(ofn, "w") as f:
            stamp(out, current="mu")
            json.dump(out, f)
        print(f"herwig: EMPTY SELECTION -- 0 events pass; wrote {ofn}")
        return

    out["means"] = {kk: float(np.average(np.asarray(obs[kk], dtype=float),
                                         weights=w))
                    for kk in ("nch", "nch1", "Q2", "xbj", "y")}
    dsum = dphi_summary(obs, w)
    out.update(dsum)
    out["means"]["dphi"] = dsum["mean_dphi"]
    print(f"herwig: <dphi_min> = {dsum['mean_dphi']:.4f} rad, "
          f"{dsum['n_no_hadron']} events ({100*dsum['frac_no_hadron']:.3f}% "
          f"of sigma) with no charged hadron -> DeltaPhi unfilled")
    msum = dmeson_summary(obs, w)
    out.update(msum)
    out["means"].update({"nd_ch": msum.get("mean_n_d_ch"),
                         "nd_0": msum.get("mean_n_d_0"),
                         "nd_s": msum.get("mean_n_d_s")})
    print(f"herwig: <N(D+-)> = {msum.get('mean_n_d_ch', 0):.4f}, "
          f"<N(D0)> = {msum.get('mean_n_d_0', 0):.4f}, "
          f"<N(Ds)> = {msum.get('mean_n_d_s', 0):.4f}; "
          f"sigma fraction with >=1: D+- {100*msum.get('frac_ge1_d_ch', 0):.2f}%"
          f", D0 {100*msum.get('frac_ge1_d_0', 0):.2f}%")

    hists = {}
    for kk, edges in BINS.items():
        v = np.asarray(obs[kk], dtype=float)
        cnt, _ = np.histogram(v, bins=edges, weights=w)
        err2, _ = np.histogram(v, bins=edges, weights=w * w)
        widths = np.diff(edges)
        hists[kk] = {"edges": edges.tolist(),
                     "dsig": (cnt * norm / widths).tolist(),
                     "err": (np.sqrt(err2) * norm / widths).tolist()}
    if not msum.get("d_supported"):
        for kk in D_HIST_KEYS:
            hists.pop(kk, None)
        print("herwig: NO D mesons in the sample -- D histograms omitted")
    out["hists"] = hists

    os.makedirs(f"{BASE}/results", exist_ok=True)
    ofn = result_path(genkey)
    with open(ofn, "w") as f:
        stamp(out, current="mu")
        json.dump(out, f)
    print(f"wrote {ofn}")


if __name__ == "__main__":
    main()
