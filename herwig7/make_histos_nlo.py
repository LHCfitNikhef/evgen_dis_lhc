#!/usr/bin/env python3
"""Build results/histos_herwig_nlo_powheg*.json from the Herwig POWHEG-matched
samples in herwig7/{mu,nu}pwg[neg]_job_*/events.hepmc, in exactly the schema written by
analysis/analyze.py (BINS imported from there, read-only).

Normalisation (WEIGHTED events, generation region WIDER than fiducial --
Q2 > 4 GeV2 with NO y cut, see the TRAP 2 comment in DIS-mu-NLO.in):

    norm       = sigma_gen / sum_of_all_weights
    sigma_fid  = norm * sum_of_fiducial_weights
    sigma_err  = norm * sqrt(sum_of_fiducial_weights^2)

sigma_gen is the INTEGRATED cross-section of the generation region,
combined over jobs by inverse variance (see herwig_xsec.py).  It excludes
the ~2.4% of events vetoed during showering, which are regenerated rather
than lost from the physics.  This is the same recipe as the "powheg" branch
of analyze.py.

WARNING: this sample remains UNUSABLE and is wired into no plot.  Its
weights have a pathological tail at y -> 1 (0.37% of events carry 81% of
the variance; one event at w = -494), so both the central value and the
error are artefacts of where those rare events happened to fall.  The
normalisation rule below is correct; the sample is not.
"""
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "analysis"))
from analyze import BINS, analyze, input_manifest  # noqa: E402
from analyze import result_path  # noqa: E402 -- the one name rule
from analyze import dmeson_summary  # noqa: E402
from herwig_xsec import combine, read_jobs  # noqa: E402
import sys as _sys, os as _os  # noqa: E402
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "analysis"))
from analyze import job_files  # noqa: E402
from analyze import ENERGY, at_energy  # noqa: E402 -- the one energy rule
from analyze import SELECTION  # noqa: E402 -- the one fiducial region
import selection  # noqa: E402  -- the one Q2 floor
import beams  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp
try:
    from analyze import dphi_summary
except ImportError:                                   # older analyze.py
    dphi_summary = None

# `make_histos_nlo.py ycut` reads the LeptonicDISCut production instead of the
# original no-y-cut sample.  The two differ in generation region (0.15 < y <
# 0.93 vs no y cut at all), so the fiducial weight fraction differs, but the
# recipe -- sigma_gen x fiducial weight fraction -- is unchanged.
# An optional SECOND argument is a heavy-flavour tag, with the same mode
# suffixes as analysis/analyze.py: `charmfinal` is the benchmark definition
# (a charm quark in the final state after shower and hadronisation).  It
# selects a subset of the sample, so sigma_fid becomes sigma_gen times the
# tagged weight fraction even when the generation cuts equal the fiducial
# ones -- see the normalisation branch below.
HARD_TAGS = {"charm": 4, "bottom": 5}


def pwg_runname(stem):
    """Herwig run name for the POWHEG samples at the ACTIVE energy.

    The anchor keeps its original card and run names so re-running it cannot
    move a published number; other energies carry the beam tag, which is what
    stops "mupwg_job_*" quietly folding 400 GeV into the 1 TeV result.
    """
    return (stem if ENERGY == beams.ANCHOR_ENERGY
            else f"{stem}-{beams.Beams('mu', ENERGY).tag}")
_tagarg = sys.argv[2].lower() if len(sys.argv) > 2 else None
TAG, TAGMODE = None, "out"
if _tagarg is not None:
    TAG = _tagarg
    if TAG not in HARD_TAGS:
        for _suf in ("final", "in", "any"):
            if TAG.endswith(_suf) and TAG[:-len(_suf)] in HARD_TAGS:
                TAG, TAGMODE = TAG[:-len(_suf)], _suf
                break
    if TAG not in HARD_TAGS:
        sys.exit(f"unknown tag '{_tagarg}' (use "
                 f"{'|'.join(HARD_TAGS)}, optionally + final|in|any)")

# MATCHBOX MC@NLO WAS REMOVED (user decision, 2026-08-25): a failed road.
# Its sample had N_eff of 6-10% and max|w| ~ 1131, so differential
# distributions were unusable, and the wide-generation variant drove the
# fiducial weight fraction to 1.2582 -- above one, which is unphysical.  The
# POWHEG-matched path below carries unit weights and is the NLO entry.
# NOTE the neutrino side's genkey "herwig_nlo" is the POWHEG sample, NOT this
# one; only the MUON "herwig_nlo" was Matchbox.  See analyze_nu.py.
if len(sys.argv) > 1 and sys.argv[1] == "powheg":
    # DIS-mu-POWHEG.in: PowhegMEDISNC = Herwig::MENeutralCurrentDIS with
    # Contribution 1, i.e. the CLASSIC matrix-element path rather than
    # Matchbox.  Two things change.  (i) The events carry UNIT WEIGHTS, so the
    # variance pathology of the Matchbox sample (N_eff 6-10%, max|w| ~ 1131)
    # is simply absent and differential distributions are usable.  (ii)
    # Generation cuts equal the fiducial cuts, so again no ratio estimator.
    JOBBASE, RUNNAME = at_energy("mupwg_job"), pwg_runname("DIS-mu-PWG")
    LABEL, GENKEY = "Herwig 7.3.0 NLO (POWHEG)", "herwig_nlo_powheg"
elif len(sys.argv) > 1 and sys.argv[1] == "powhegneg":
    # THE NEGATIVE HALF of the NLO cross section.  Herwig's DISBase splits the
    # NLO weight into its positive and negative parts and generates them as
    # two separate POSITIVE-weight samples, selected by the Contribution
    # switch (DISBase.cc::NLOWeight returns max(0,w) for 1 and max(0,-w) for
    # 2).  Herwig's own defaults put Contribution 1 on PowhegMEDISNC, so a
    # single run out of the box is only half the answer.  The physical result
    # is pos - neg, assembled by herwig7/combine_nlo.py.
    JOBBASE, RUNNAME = at_energy("mupwgneg_job"), pwg_runname("DIS-mu-PWGNEG")
    LABEL, GENKEY = "Herwig 7.3.0 NLO (POWHEG, negative)", "herwig_nlo_powheg_neg"
elif len(sys.argv) > 1 and sys.argv[1] in ("v2powheg", "v2powhegneg"):
    # PAPER PLOTS V2 (user, 2026-09-13): one nucleon per run, final region
    # (PRODUCTION.md).  Generated at MinQ2 3.5 and MinW2 5, looser than
    # q4w3 on purpose (herwig7/production/), so this is the projection branch below:
    # sigma_fid = sigma_gen x fiducial fraction.  result_path() turns the
    # "v2_" key into histos_herwig_nlo_powheg[_neg]_q4w3_<t>[_tag].json.
    import sample_layout  # noqa: E402
    _neg = sys.argv[1].endswith("neg")
    _t = sample_layout.require()[1]
    JOBBASE = at_energy(sample_layout.herwig_job_base("mu", _t, _neg))
    RUNNAME = sample_layout.herwig_runname("mu", _t, ENERGY, _neg)
    LABEL = ("Herwig 7.3.0 NLO (POWHEG" + (", negative" if _neg else "")
             + ")")
    GENKEY = "v2_herwig_nlo_powheg" + ("_neg" if _neg else "")
elif len(sys.argv) > 1 and sys.argv[1] == "powheglo":
    # DIAGNOSTIC, for TODO P3b.  DIS-mu-POWHEG-LO.in is the SAME card as the
    # NLO one with `PowhegMEDISNC:Contribution 0` -- a LEADING-ORDER matrix
    # element -- while `ShowerHandler:HardEmission POWHEG` stays ON.
    #
    # It exists to separate the two remaining explanations of Herwig NLO's
    # 3.7% <Q2> excess.  If the excess is the POWHEG hardest emission and its
    # recoil under the angular-ordered shower, this sample shows it, because
    # the hard emission is still there.  If the excess is the Contribution 1/2
    # weight split, this sample is CLEAN, because there is no split at
    # Contribution 0.  Herwig's ordinary LO sample is already known clean
    # (+0.70% at 1 TeV, like Pythia LO), so the two bracket the question.
    JOBBASE, RUNNAME = "pwgvar_LO", pwg_runname("DIS-mu-PWGLO")
    LABEL = "Herwig 7.3.0 LO ME + POWHEG hard emission (diagnostic)"
    GENKEY = "herwig_pwg_lo"
elif len(sys.argv) > 1 and sys.argv[1].startswith("knob"):
    # A VARIANT FROM run_herwig_knob.sh: the same card with a few `set` lines
    # inserted, run into knob<label>_<LO|NLO>/.  Used for the scan over the
    # settings that might stop Herwig discarding 5.9% of its events; the label
    # travels into the result key so several variants coexist on disk.
    # knob<label>_<lo|nlo>
    _spec = sys.argv[1][4:]
    _lab, _, _ord = _spec.rpartition("_")
    _ord = (_ord or "lo").upper()
    if not _lab.isalnum() or _ord not in ("LO", "NLO", "NEG"):
        sys.exit("usage: make_histos_nlo.py knob<label>_<lo|nlo|neg> "
                 "[charm tag]")
    JOBBASE = f"knob{_lab}_{_ord}"
    _stem = {"LO": "DIS-mu-PWGLO", "NLO": "DIS-mu-PWG",
             "NEG": "DIS-mu-PWGNEG"}[_ord]
    RUNNAME = pwg_runname(f"{_stem}-KNOB{_lab}")
    LABEL = f"Herwig 7.3.0 {_ord} + POWHEG, variant {_lab}"
    GENKEY = f"herwig_knob{_lab}_{_ord.lower()}"
elif len(sys.argv) > 1 and sys.argv[1].startswith("maxtry"):
    # THE Q2-TILT EXPERIMENT (user, 2026-09-01).  Same card, same matrix
    # element, same hard emission -- only ShowerHandler:MaxTry raised above
    # the default 100, so that events Herwig would have DISCARDED after 100
    # failed shower attempts are kept instead.  The tilt tracks the discard
    # rate across every sample here (see run_herwig_maxtry.sh's header), so
    # this is the direct test of whether it is bookkeeping or physics.
    # maxtry<N>_<lo|nlo>
    _spec = sys.argv[1][6:]
    _n, _, _ord = _spec.rpartition("_")
    _ord = (_ord or "lo").upper()
    if not _n.isdigit() or _ord not in ("LO", "NLO"):
        sys.exit("usage: make_histos_nlo.py maxtry<N>_<lo|nlo> [charm tag]")
    JOBBASE = f"maxtry{_n}_{_ord}"
    RUNNAME = pwg_runname(f"DIS-mu-PWG{'LO' if _ord == 'LO' else ''}"
                          f"-MAXTRY{_n}")
    LABEL = f"Herwig 7.3.0 {_ord} + POWHEG, ShowerHandler:MaxTry {_n}"
    GENKEY = f"herwig_maxtry{_n}_{_ord.lower()}"
elif len(sys.argv) > 1 and sys.argv[1].startswith("evo"):
    # P3b DIAGNOSTIC: the same POWHEG-matched run with a different shower
    # EVOLUTION SCHEME, which in Herwig 7.3 is the DIS recoil choice
    # (arXiv:1904.11866; the older ReconstructionOption switch is deprecated
    # in favour of it).  `evoPT_lo`, `evoQ2_nlo` and so on -- the scheme and
    # the matrix-element order, matching run_herwig_evoscheme.sh's directory
    # naming.
    # evo<Scheme>[HE<HardEmission>]_<lo|nlo>
    _spec = sys.argv[1][3:]
    _sch, _, _ord = _spec.rpartition("_")
    _ord = (_ord or "lo").upper()
    _he = ""
    if "HE" in _sch:
        _sch, _, _he = _sch.partition("HE")
    if _sch not in ("DotProduct", "Q2", "pT") or _ord not in ("LO", "NLO"):
        sys.exit("usage: make_histos_nlo.py "
                 "evo<DotProduct|Q2|pT>[HE<MECorrection|None|POWHEG>]_<lo|nlo>")
    JOBBASE = f"evo{_sch}{'HE' + _he if _he else ''}_{_ord}"
    RUNNAME = pwg_runname(f"DIS-mu-PWG{'LO' if _ord == 'LO' else ''}"
                          f"-EVO{_sch}" + (f"-HE{_he}" if _he else ""))
    LABEL = (f"Herwig 7.3.0 {_ord} + POWHEG, EvolutionScheme {_sch}"
             + (f", HardEmission {_he}" if _he else "") + " (diagnostic)")
    GENKEY = (f"herwig_evo_{_sch.lower()}"
              + (f"_he{_he.lower()}" if _he else "") + f"_{_ord.lower()}")
else:
    sys.exit("usage: make_histos_nlo.py "
             "<powheg|powhegneg|powheglo|evo...|maxtry<N>_<lo|nlo>> "
             "[charm tag]\n"
             "  (the Matchbox MC@NLO modes were removed -- see above)")


def main():
    # exact job-number match, not a prefix glob -- see
    # analyze.job_files(): "<base>_*" also matches the energy scan's
    # own <base>_400GeV_N directories.
    files = job_files(f"{HERE}", JOBBASE, "events.hepmc")
    if not files:
        sys.exit(f"no herwig7/{JOBBASE}_*/events.hepmc files found")
    jobs = read_jobs(f"{HERE}/{JOBBASE}_*/{RUNNAME}-S*.out")
    if len(jobs) != len(files):
        sys.exit(f"MISMATCH: {len(files)} events.hepmc files but {len(jobs)} "
                 f".out files -- refusing to normalise a mixed sample.")
    sigma_sample_pb, sigma_err_pb, n_gen_tot, survival = combine(jobs)
    print(f"  sigma_gen_region = {sigma_sample_pb/1e3:.4f} "
          f"+- {sigma_err_pb/1e3:.4f} nb over {len(jobs)} jobs, "
          f"{n_gen_tot} events, survival {survival:.4f}")

    hard_flavour = HARD_TAGS[TAG] if TAG else None
    tag_stats = {} if TAG else None
    genkey = GENKEY + (f"_{TAG}{'' if TAGMODE == 'out' else TAGMODE}"
                       if TAG else "")
    label = LABEL + (f" ({TAG})" if TAG else "")
    obs, w, n_events, sum_w_all = analyze(files, hard_flavour=hard_flavour,
                                          tag_stats=tag_stats,
                                          hard_tag_mode=TAGMODE)
    print(f"{genkey}: parsed {n_events} events, {len(w)} in fiducial region")
    if tag_stats is not None:
        print(f"{genkey}: tag comparison -- hard-out "
              f"{tag_stats.get('n_tagged_out', 0)}, final-state "
              f"{tag_stats.get('n_tagged_final', 0)} of {n_events} parsed")
    if n_events != n_gen_tot:
        print(f"WARNING: parsed {n_events} != generated {n_gen_tot}")

    out = {"generator": genkey, "label": label,
           "hard_tag_mode": TAGMODE,
           "n_parsed": n_events, "n_fiducial": len(w),
           "sigma_gen_region_pb": sigma_sample_pb,
           "sigma_gen_region_err_pb": sigma_err_pb,
           "veto_survival": survival}
    if GENKEY.startswith("v2_"):
        out["target"] = sample_layout.require()[1]
    out["inputs"] = input_manifest(
        files, [j.sigma_pb for j in jobs],
        [j.n_written for j in jobs], f" [{genkey}]")
    # weighted, generation region wider than fiducial
    norm = sigma_sample_pb / sum_w_all
    # "CUTS COINCIDE" IS A STATEMENT ABOUT THE FIDUCIAL REGION, not about the
    # job directory.  The cards generate at Q2 > 4, 0.2 < y < 0.9 under the
    # inclusive selection, so pinning sigma_fid to the integrator is right
    # THERE and nowhere else.  Testing only the job base pinned a Q2 > 11 run
    # to the Q2 > 4 integrated cross-section, reporting 36221 pb with the
    # event count correctly down to a third -- the same silent failure as the
    # Tier S pinning bug and as analyze.py's `subset` test, found on the first
    # raised-floor run (2026-09-01).  Three copies of one rule, so all three
    # now say the same thing.
    coincide = (SELECTION.name == selection.DEFAULT
                and SELECTION.q2_min == selection.Q2_FLOOR_DEFAULT)
    if JOBBASE.startswith(("mupwg_job", "mupwgneg_job")) \
            and hard_flavour is None and coincide:
        # cuts coincide: no projection, and the quoted error is the job
        # spread rather than the (over-optimistic) inverse-variance one
        out["sigma_fid_pb"] = sigma_sample_pb
    else:
        # A CHANNEL TAG always needs the projection, even when the generation
        # cuts equal the fiducial ones: the tag selects a subset, so the
        # channel's contribution is sigma_gen times the TAGGED weight
        # fraction.  Handing back the full sigma_gen here would have reported
        # the inclusive rate as the charm rate.
        out["sigma_fid_pb"] = norm * float(w.sum())
    # MC error on the fiducial weight sum, plus the integration error on the
    # overall normalisation (independent, so in quadrature)
    stat = float(np.sqrt((w * w).sum())) / float(w.sum()) if w.sum() else 0.0
    rel_int = sigma_err_pb / sigma_sample_pb if sigma_sample_pb else 0.0
    out["sigma_fid_err_pb"] = out["sigma_fid_pb"] * math.hypot(stat, rel_int)
    out["mean_weight_pb"] = sum_w_all / n_events
    # N_eff = (sum w)^2 / sum w^2.  This sample's weight distribution has a
    # pathological tail at y -> 1 (see the project notes): N_eff is a small
    # fraction of N, which is why the quoted error is large and why the
    # central value is NOT trustworthy regardless of the normalisation rule.
    out["n_eff"] = (float(w.sum())**2 / float((w * w).sum())
                    if (w * w).sum() else 0.0)
    print(f"{genkey}: N_eff = {out['n_eff']:.0f} of {len(w)} fiducial events "
          f"({100*out['n_eff']/max(len(w), 1):.1f}%)")
    print(f"{genkey}: fiducial weight fraction = {w.sum()/sum_w_all:.4f}, "
          f"negative-weight events = {int((w < 0).sum())}/{len(w)}")
    print(f"{genkey}: sigma_fid = {out['sigma_fid_pb']/1e3:.3f} +- "
          f"{out['sigma_fid_err_pb']/1e3:.3f} nb")

    out["means"] = {kk: float(np.average(np.asarray(obs[kk], dtype=float),
                                         weights=w))
                    for kk in ("nch", "nch1", "nch05", "Q2", "xbj", "y",
                               "nd_ch", "nd_0", "nd_s")}
    # D-species summary (frac_ge1_d_ch / frac_ge1_d_0 / ...), the same block
    # every other sample carries.  It was missing here, which left the Herwig
    # NLO row out of the D-species table and out of the charm-fraction plot's
    # D panels.
    out.update(dmeson_summary(obs, w))
    if dphi_summary is not None and "dphi" in BINS:
        dsum = dphi_summary(obs, w)
        out.update(dsum)
        out["means"]["dphi"] = dsum["mean_dphi"]

    hists = {}
    for kk, edges in BINS.items():
        v = np.asarray(obs[kk], dtype=float)
        cnt, _ = np.histogram(v, bins=edges, weights=w)
        err2, _ = np.histogram(v, bins=edges, weights=w * w)
        widths = np.diff(edges)
        hists[kk] = {"edges": edges.tolist(),
                     "dsig": (cnt * norm / widths).tolist(),
                     "err": (np.sqrt(err2) * norm / widths).tolist()}
    out["hists"] = hists

    os.makedirs(f"{BASE}/results", exist_ok=True)
    # result_path(), NOT a bare histos_<genkey>.json.  Making the INPUT job
    # base energy-aware while leaving the OUTPUT name fixed is precisely how a
    # 4 TeV run overwrites the 1 TeV anchor -- which is exactly what happened
    # here on 2026-08-25, caught by the energy_gev stamp rather than by any
    # error.  The two must move together, always.
    ofn = result_path(genkey)
    with open(ofn, "w") as f:
        stamp(out, current="mu")
        json.dump(out, f)
    print(f"wrote {ofn}")


if __name__ == "__main__":
    main()
