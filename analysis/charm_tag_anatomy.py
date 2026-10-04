#!/usr/bin/env python3
"""Take the charm tag apart: WHICH charm each definition counts, and how soft.

WHY (user, 2026-09-01).  Two questions were asked together, and they turn out
to be one question:

  * "Investigate why Sherpa NLO predicts more charm production, also at the
    total charm cross-section level ... it is odd that Sherpa gives the
    correct charm cross-sections in muon DIS but not neutrino DIS."
  * "You write that an initial-state anti-charm leaves its partner charm in
    the beam remnant, but maybe this partner charm is very soft and would not
    be reconstructed in a real experiment?  We should for sure correct for
    this gap."

THE BENCHMARK'S CHARM TAG IS HADRON-LEVEL: an event counts as charm when a
prompt charm hadron exists after shower and hadronisation, with the b -> c
cascade removed.  That definition is right for a comparison with an
experiment, and it is the one every published number here uses -- but it says
nothing about WHERE the charm came from or whether anything could see it.
This script answers both, per event, on the same samples the results are
built from:

  hard_out    the hard process has an outgoing charm quark
  hard_any    the hard process has charm anywhere, INCLUDING an incoming
              c or cbar -- the beam-remnant partner lives here
  final       a prompt charm hadron exists (the benchmark tag)
  E, theta    of the LEADING prompt charm hadron, in the lab

and then reports what an energy threshold on that hadron would do to each
generator's charm cross-section.  A threshold is the closest thing to "would
a detector see it" that this benchmark can express without a detector
simulation, and it is applied identically to every generator.

Usage:
  analysis/charm_tag_anatomy.py [--current mu|nu] [--max N] [key ...]

With no keys it runs the NLO matchings of the chosen current.  Writes
results{,_nu}/charm_anatomy_<key>.json and prints the table.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                        # noqa: E402

import analyze                                            # noqa: E402
import beams                                              # noqa: E402
import paths                                              # noqa: E402
import selection                                          # noqa: E402

BASE = paths.REPO

# The energy thresholds scanned on the leading charm hadron.  0 is the
# benchmark's own tag and must come first, so every other row is read as a
# correction to it.
E_CUTS = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0)

# key -> beam pid.  THE FILE LIST IS NOT A GLOB (2026-09-01, after a glob
# cost a wrong answer): `powheg/job_*` matches job_400GeV_* and job_4TeV_*
# as well as the anchor, and `herwig7/nupwg_job_*` matched eighteen
# directories across three beam energies.  Mixing energies in a charm
# fraction is exactly the failure this benchmark keeps meeting -- it
# completes, it looks plausible, and it is wrong.
#
# So the sample is read from the PUBLISHED result's own input manifest:
# whatever files produced `histos_<key>_charmfinal.json` are the files the
# anatomy explains, by construction, and a moved sample raises rather than
# quietly changing the answer.
SAMPLES = {
    "mu": [("powheg", 13), ("sherpa", 13), ("herwig_nlo_powheg", 13),
           # the LO samples (2026-09-03): the same split at leading order,
           # where there is no matching to move charm between the pieces
           ("sherpa_lo", 13), ("pythia", 13), ("herwig", 13)],
    "nu": [("powheg_nu", 14), ("sherpa_nlo", 14), ("herwig_nlo", 14),
           ("sherpa_lo", 14), ("pythia", 14), ("herwig", 14),
           ("genie_lo", 14), ("genie", 14)],
}
# The beam energy, through the same knob as every other analysis: the keys
# above name the 1 TeV anchor and BENCH_ENERGY=400 or 4000 moves them to
# the scan's other points, with results written under the energy-suffixed
# name (charm_anatomy_<key>_400GeV.json).
ENERGY = float(os.environ.get("BENCH_ENERGY", beams.ANCHOR_ENERGY))


def at_e(name):
    return beams.at_energy(name, ENERGY)


def sample_files(resdir, key):
    """The exact files the published charm result for `key` was built from."""
    p = f"{resdir}/histos_{at_e(key + '_charmfinal')}.json"
    if not os.path.exists(p):
        raise SystemExit(f"no published charm result at {p} -- the anatomy "
                         f"explains a number, so that number must exist")
    with open(p) as f:
        d = json.load(f)
    files = [os.path.join(BASE, r["path"])
             for r in (d.get("inputs") or {}).get("files", [])]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        raise SystemExit(f"{key}: {len(missing)} file(s) in the manifest are "
                         f"gone, first {missing[0]}")
    return files, d


def lab_kinematics(p4, P):
    """(E, theta) of a particle IN THE PROTON REST FRAME.

    >>> THE EVENT RECORD IS NOT IN THE LAB. <<<  POWHEG writes a symmetric
    pair at sqrt(s)/2 -- 21.66 GeV a side here -- so reading p4[0] as an
    energy gives charm hadrons of a few GeV where the lab spectrum runs to
    hundreds, and every conclusion about whether the charm is "soft" would
    be drawn in the wrong frame.  The analysis has always known this:
    `lab_energy` is dot(p, P) / M_P, covariant and frame-independent, and
    the same construction is used here.

    theta follows `lepton_theta`: the boost between the generation frame and
    the lab is along z, so p_T is invariant, and |p| in the lab comes from
    the covariant energy and the particle's own invariant mass.  Forward and
    backward are not distinguished, which is immaterial for the forward charm
    this script is about.
    """
    e_lab = analyze.dot(p4, P) / analyze.M_P
    m2 = max(analyze.dot(p4, p4), 0.0)
    pmag = np.sqrt(max(e_lab * e_lab - m2, 1e-12))
    pt = np.hypot(p4[1], p4[2])
    return float(e_lab), float(np.arcsin(min(pt / pmag, 1.0)))


def anatomy(files, beam_pid, max_events=None):
    """One pass over a sample; returns the per-event arrays it needs.

    The fiducial region is the benchmark's own, read from `selection`, so this
    script and the cross-sections it is compared against cut identically.
    """
    w, ho, ha, fin, ce, cth = [], [], [], [], [], []
    n_seen = 0
    for fn in files:
        for ev in analyze.parse_hepmc3(fn, beam_pid=beam_pid, hard_graph=True):
            weight, k, P, parts, _dm, hard = ev
            n_seen += 1
            if max_events and n_seen > max_events:
                break
            # THE FIDUCIAL CUT IS THE ANALYSIS' OWN, imported rather than
            # rewritten: `lab_energy`, `dis_invariants`, `lepton_theta` and
            # `SELECTION` are the same objects analyze.py cuts with, so this
            # script cannot drift away from the cross-sections it explains.
            # The outgoing lepton is the beam species in both currents here
            # (mu- NC, and mu- from nu_mu CC).
            best, best_elab = None, -1.0
            for pid, p in parts:
                if pid == 13:
                    el = analyze.lab_energy(p, P)
                    if el > best_elab:
                        best, best_elab = p, el
            if best is None:
                continue
            q2, y, _x, _kP = analyze.dis_invariants(k, P, best)
            sel = analyze.SELECTION
            if not sel.passes(q2, y):
                continue
            if not sel.passes_lepton(best_elab,
                                     analyze.lepton_theta(best_elab, best)):
                continue
            hin, hout = hard.get("in", []), hard.get("out", [])
            charm = [(pdg, p4) + lab_kinematics(p4, P)
                     for pdg, p4 in (hard.get("charm_final") or [])]
            lead = max(charm, key=lambda c: c[2]) if charm else None
            w.append(weight)
            ho.append(any(abs(p) == 4 for p in hout))
            ha.append(any(abs(p) == 4 for p in hin + hout))
            fin.append(4 in hard.get("fs_heavy", frozenset()))
            ce.append(lead[2] if lead else 0.0)
            cth.append(lead[3] if lead else -1.0)
        if max_events and n_seen > max_events:
            break
    return (np.array(w), np.array(ho), np.array(ha), np.array(fin),
            np.array(ce), np.array(cth))


def summarise(w, ho, ha, fin, ce, cth):
    """Everything the report and the paper need from one sample."""
    sw = float(w.sum())
    out = {"n_events": int(w.size), "sum_w": sw,
           "frac_hard_out": float(w[ho].sum() / sw) if sw else 0.0,
           "frac_hard_any": float(w[ha].sum() / sw) if sw else 0.0,
           "frac_final": float(w[fin].sum() / sw) if sw else 0.0}
    # THE THREE DISJOINT PIECES of the hadron-level tag.  This is the split
    # the investigation turns on: `remnant` is the beam-remnant partner charm
    # the standing explanation blames, `shower` is charm the hard record does
    # not contain at all.
    out["split"] = {
        "hard_out": float(w[fin & ho].sum() / sw) if sw else 0.0,
        "remnant": float(w[fin & ha & ~ho].sum() / sw) if sw else 0.0,
        "shower": float(w[fin & ~ha].sum() / sw) if sw else 0.0,
        "hard_only": float(w[ho & ~fin].sum() / sw) if sw else 0.0,
    }
    # what an energy threshold on the leading charm hadron would remove, and
    # from WHICH piece -- applied identically to every generator
    out["e_cuts"] = {}
    for c in E_CUTS:
        keep = fin & (ce > c)
        out["e_cuts"][f"{c:g}"] = {
            "frac": float(w[keep].sum() / sw) if sw else 0.0,
            "hard_out": float(w[keep & ho].sum() / sw) if sw else 0.0,
            "remnant": float(w[keep & ha & ~ho].sum() / sw) if sw else 0.0,
            "shower": float(w[keep & ~ha].sum() / sw) if sw else 0.0,
        }
    # the spectra themselves, so the added charm can be SHOWN to be soft
    # rather than asserted
    edges = np.geomspace(0.5, 2000.0, 41)
    for nm, m in (("all", fin), ("hard_out", fin & ho),
                  ("remnant", fin & ha & ~ho), ("shower", fin & ~ha)):
        h, _ = np.histogram(ce[m], bins=edges, weights=w[m])
        out[f"E_lead_{nm}"] = [float(v) for v in h]
        out[f"E_lead_{nm}_median"] = (
            float(np.median(ce[m])) if m.sum() else None)
    out["E_lead_edges"] = [float(v) for v in edges]
    tedges = np.geomspace(1e-4, 1.0, 31)
    h, _ = np.histogram(cth[fin & (cth > 0)], bins=tedges,
                        weights=w[fin & (cth > 0)])
    out["theta_lead"] = [float(v) for v in h]
    out["theta_lead_edges"] = [float(v) for v in tedges]
    return out


def main():
    argv = sys.argv[1:]
    current = "nu"
    max_events = None
    keys = []
    i = 0
    while i < len(argv):
        if argv[i] == "--current":
            current = argv[i + 1]
            i += 2
        elif argv[i] == "--max":
            max_events = int(argv[i + 1])
            i += 2
        else:
            keys.append(argv[i])
            i += 1
    table = SAMPLES[current]
    if keys:
        table = [row for row in table if row[0] in keys]
    resdir = f"{BASE}/results" if current == "mu" else f"{BASE}/results_nu"
    for key, pid in table:
        files, published = sample_files(resdir, key)
        print(f"{key}: {len(files)} file(s) from the published manifest")
        res = summarise(*anatomy(files, pid, max_events))
        res["key"] = key
        res["files"] = [os.path.relpath(f, BASE) for f in files]
        # the number this anatomy explains, carried alongside it so a reader
        # (and tools/check_charm_anatomy.py) can see they are the same sample
        res["published_charm_pb"] = published.get("sigma_fid_pb")
        incl = f"{resdir}/histos_{at_e(key)}.json"
        if os.path.exists(incl):
            with open(incl) as f:
                res["published_inclusive_pb"] = json.load(f).get(
                    "sigma_fid_pb")
        res["selection"] = dict(analyze.SELECTION.__dict__)
        res["energy_gev"] = ENERGY
        with open(f"{resdir}/charm_anatomy_{at_e(key)}.json", "w") as f:
            json.dump(res, f, indent=1)
        s = res["split"]
        print(f"  tagged final {100*res['frac_final']:.2f}%  "
              f"= hard-out {100*s['hard_out']:.2f}"
              f" + remnant {100*s['remnant']:.2f}"
              f" + shower {100*s['shower']:.2f}"
              f"   (hard-out but not final: {100*s['hard_only']:.2f})")
        print(f"  median E of the leading charm hadron: "
              f"all {res['E_lead_all_median']}, "
              f"remnant {res['E_lead_remnant_median']}, "
              f"shower {res['E_lead_shower_median']}")
        for c in E_CUTS:
            e = res["e_cuts"][f"{c:g}"]
            print(f"    E > {c:5g} GeV: {100*e['frac']:6.2f}%  "
                  f"(hard-out {100*e['hard_out']:5.2f}, "
                  f"remnant {100*e['remnant']:5.2f}, "
                  f"shower {100*e['shower']:5.2f})")
        print(f"  wrote {resdir}/charm_anatomy_{key}.json")


if __name__ == "__main__":
    main()
