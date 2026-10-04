#!/usr/bin/env python3
"""The charged-particle multiplicity DONUT measured, predicted by our samples.

WHAT IS BEING COMPARED (user, 2026-09-10).  DONUT (arXiv:0711.0728) published
in its FIG. 9 the number of charged particles at the primary vertex of all 578
neutrino interactions it located in emulsion, against its own LEPTO-based
simulation.  This predicts the same distribution from the benchmark's
generators and puts them beside it.

>>> THE MULTIPLICITY AT THE PRIMARY VERTEX IS FLAVOUR-BLIND, WHICH IS WHAT
MAKES THIS POSSIBLE AT ALL. <<<  DONUT's beam is 57% nu_mu, 38% nu_e and 5%
nu_tau (its FIG. 2), and this benchmark has no electron- or tau-neutrino
samples.  It does not need them: every charged-current event puts exactly ONE
lepton track at the vertex, whether that lepton is an electron, a muon or a
tau -- the tau decays microns downstream and DONUT counts the parent.  The
hadronic side does not know which lepton was made.  So a nu_mu sample predicts
the multiplicity of the whole beam, and the user's own argument for the
cross-sections ("neutrino cross-sections are flavour independent") carries over
to the multiplicity for a reason of its own.

WHAT IS NOT CARRIED OVER, said here rather than left implicit: the tau mass
suppresses the nu_tau cross-section at low y, which tilts the hadronic energy
of the 5% of the sample that is nu_tau.  And the electron of a nu_e CC event
showers in the steel, which is why those events are the ones DONUT most often
FAILED to locate -- an efficiency, not a multiplicity.

THE DEFINITION OF n_ch, AND WHY IT NEEDS A CUT.  DONUT counts tracks
RECONSTRUCTED IN EMULSION.  A track has to leave the 1 mm steel plate it was
born in and be seen in the next emulsion layer, and the automated scanning
accepts a limited angular range, so a generator's bare charged multiplicity is
not the measured object: it counts particles too slow to escape the plate and
too steep to be scanned.  The default here is

    charged, stable (the benchmark's ctau > 10 mm), p > 0.3 GeV/c,
    tan(theta) < 0.5 to the beam, and the primary lepton counted

and three variants are computed beside it so the sensitivity to that choice is
a number rather than an assumption.  The paper states no acceptance, so this
is the honest way round: quote the spread.

>>> AND THE PREDICTION HAS NO DETECTOR IN IT. <<<  DONUT's own Monte Carlo
carries the trigger, the scanning and the LOCATION EFFICIENCY, which is 0.31
for shower-like events against 0.77 for the rest (their Section VII D) and so
depends on the topology -- and therefore on the multiplicity itself.  Ours
carries none of that.  The comparison is therefore against DONUT's MC as much
as against its data: where our curve and their MC differ, that is physics;
where their MC and their data differ, that is either physics or their
detector, and we cannot tell which.  Everything is normalised to the same
number of events, so only shapes are compared.

Usage:
  analysis/donut_nch.py genie    genie/donut_job_*/events.hepmc
  analysis/donut_nch.py powheg   --ladder <ladder-donut dir>
  analysis/donut_nch.py --merge
Writes results_nu/donut_nch_<generator>.json, and --merge builds
results_nu/donut_nch.json with the data beside every generator.
"""
import glob
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NEUTRINOS = {12, 14, 16}

DATA = f"{BASE}/data/donut/nch.json"
RESULTS = f"{BASE}/results_nu"

# Charged, stable under the benchmark's ctau > 10 mm convention.  The leptons
# are here too: at the primary vertex a muon and a pion are one track each,
# and DONUT could not tell them apart without the spectrometer.
CHARGED = {211, 321, 2212, 3112, 3222, 3312, 3334, 11, 13, 15}

# (name, p_min in GeV/c, tan(theta) max).  The first is the default; the rest
# exist so that "how much does the acceptance matter" is answered with a
# number.  `raw` is the generator's own charged multiplicity with nothing
# applied, which is what a comparison that ignored the emulsion would show.
VARIANTS = [("default", 0.3, 0.5),
            ("raw", 0.0, float("inf")),
            ("soft", 0.1, 0.5),
            ("wide", 0.3, 1.0)]
NBIN = 13                     # FIG. 9 draws n_ch = 1 .. 13


def die(msg):
    sys.exit(f"donut_nch: {msg}")


def read_events(fname):
    """Yield (target four-momentum, [(pid, p4), ...]) per event, from HepMC3.

    >>> A DEDICATED READER, AND THE REASON IS THE FRAME. <<<  The two samples
    this study compares are written in DIFFERENT frames and nothing in the
    file says so: GENIE generates with the nucleus AT REST, so its records are
    already the laboratory; POWHEG-V2 writes a symmetric massless pair at
    sqrt(s)/2 and Pythia showers it there, so its records are the centre of
    mass.  An angular cut applied to both without asking would compare a lab
    angle with a centre-of-mass one -- and it does not error, it just returns
    a multiplicity four times too small, which is how this was found.

    So the TARGET RECORD is read from each event and everything is boosted to
    its rest frame.  When the target is already at rest the boost is the
    identity, so one code path serves both.

    analyze.parse_hepmc3 is not used because it identifies the target beam by
    `pid == 2212 or pid == 2112` and DONUT's target is a NUCLEUS -- iron in
    GENIE's records.  Widening that function would mean widening the copy
    tools/check_parser_drift.py holds identical to it, for one study.
    """
    beams, parts = [], []
    started = False
    with open(fname, errors="replace") as fh:
        for line in fh:
            if line.startswith("E "):
                if started:
                    yield beams, parts
                beams, parts, started = [], [], True
            elif line.startswith("P "):
                # P <id> <parent> <pid> <px> <py> <pz> <E> <m> <status>
                c = line.split()
                if len(c) < 10:
                    continue
                pid, status = int(c[3]), int(c[9])
                if status not in (1, 4):
                    continue
                p = (float(c[7]), float(c[4]), float(c[5]), float(c[6]))
                (beams if status == 4 else parts).append((pid, p))
    if started:
        yield beams, parts


def target_of(beams):
    """The struck object's four-momentum: the status-4 record that is not the
    neutrino.  Returns None when the file carries no beam records at all."""
    cand = [p for pid, p in beams if abs(pid) not in NEUTRINOS]
    if not cand:
        return None
    # the heaviest, so a nucleus wins over anything else the writer put there
    return max(cand, key=lambda p: p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3])


def to_rest(p, beta):
    """Boost one four-momentum by -beta (beta = target velocity)."""
    bx, by, bz = beta
    b2 = bx*bx + by*by + bz*bz
    if b2 < 1e-16:
        return p
    g = 1.0 / math.sqrt(1.0 - b2)
    e, px, py, pz = p
    bp = bx*px + by*py + bz*pz
    k = (g - 1.0) / b2
    return (g * (e - bp),
            px + (k * bp - g * e) * bx,
            py + (k * bp - g * e) * by,
            pz + (k * bp - g * e) * bz)


def count(parts, p_min, tan_max):
    """n_ch of one event under one acceptance."""
    n = 0
    for pid, p in parts:
        if abs(pid) not in CHARGED:
            continue
        e, px, py, pz = p
        pt = math.hypot(px, py)
        mom = math.sqrt(px*px + py*py + pz*pz)
        if mom < p_min:
            continue
        # tan(theta) to the beam axis.  A backward track (pz < 0) has no
        # finite tan(theta) in the forward sense and is never scanned, so it
        # is dropped rather than folded to a small angle.
        if pz <= 0.0:
            continue
        if pt / pz > tan_max:
            continue
        n += 1
    return n


def histogram(files):
    """The raw (NBIN+2, ) count arrays per variant, plus the event totals."""
    hist = {name: np.zeros(NBIN + 2) for name, _p, _t in VARIANTS}
    n_ev = n_cc = n_noframe = 0
    for f in files:
        for beams, parts in read_events(f):
            P = target_of(beams)
            if P is None or P[0] <= 0.0:
                n_noframe += 1
                continue
            n_ev += 1
            beta = (P[1] / P[0], P[2] / P[0], P[3] / P[0])
            lab = [(pid, to_rest(p, beta)) for pid, p in parts]
            # CC or NC, read off the final state rather than from a process
            # code: the HepMC conversion does not carry GENIE's channel, and
            # a charged lepton at the vertex IS the charged current.
            if any(abs(pid) in (11, 13, 15) for pid, _p in lab):
                n_cc += 1
            for name, p_min, tan_max in VARIANTS:
                n = count(lab, p_min, tan_max)
                # bin 0 holds n_ch = 0, bins 1..13 the drawn range, and the
                # last bin the overflow, kept rather than dropped so the
                # normalisation below can say what fraction it discards
                hist[name][min(n, NBIN + 1)] += 1.0
    if n_noframe:
        # NOT skipped quietly: an event with no target record cannot be
        # boosted, and a sample of them would silently become a small sample
        die(f"{n_noframe} event(s) carry no status-4 target record, so the "
            f"frame is unknown; refusing to mix them with the rest")
    return hist, n_ev, n_cc


def fold_ladder(ladder):
    """Fold a ladder of fixed-energy samples with DONUT's own FIG. 2 spectrum.

    THE SHAPE IS INTERPOLATED IN log E, not assigned to the nearest rung.  The
    multiplicity moves fastest at the bottom of the spectrum, where the rungs
    are furthest apart in energy and closest in log; a nearest-rung assignment
    would put a step in the folded distribution that is the ladder's spacing
    rather than the physics.

    FIG. 2 IS ALREADY sigma(E) x flux(E) x acceptance -- it is the spectrum of
    neutrinos that INTERACTED -- so each rung enters with the spectrum's own
    weight and POWHEG's cross-section never appears.  Multiplying by it as
    well would count the energy dependence of the cross-section twice.
    """
    flux = json.load(open(f"{BASE}/data/donut/flux.json"))
    edges = np.array(flux["e_edges_gev"])
    centres = 0.5 * (edges[:-1] + edges[1:])
    weight = np.zeros(len(centres))
    for k in flux["spectra"]:
        weight += np.array(flux["spectra"][k])

    rungs = sorted(ladder)
    logs = np.log(rungs)
    out = {name: np.zeros(NBIN + 2) for name, _p, _t in VARIANTS}
    n_ev = n_cc = 0
    used = 0.0
    for e_c, w in zip(centres, weight):
        if w <= 0.0:
            continue
        used += w
        le = math.log(max(e_c, 1e-3))
        for name in out:
            shapes = np.array([ladder[r]["shape"][name] for r in rungs])
            # numpy interpolates each bin independently, clamping outside the
            # ladder -- the lowest rung stands in for everything below it and
            # says so in the result
            out[name] += w * np.array(
                [np.interp(le, logs, shapes[:, j]) for j in range(NBIN + 2)])
    for r in rungs:
        n_ev += ladder[r]["n_ev"]
        n_cc += ladder[r]["n_cc"]
    return out, n_ev, n_cc, {"rungs_gev": rungs,
                             "spectrum_weight": float(used),
                             "below_lowest_rung_frac": float(
                                 weight[centres < min(rungs)].sum()
                                 / weight.sum()),
                             "above_highest_rung_frac": float(
                                 weight[centres > max(rungs)].sum()
                                 / weight.sum())}


def run_ladder(gen, d):
    """Every E<GeV> subdirectory of a ladder, folded into one prediction."""
    import re
    dirs = sorted(glob.glob(os.path.join(d, "E*")))
    ladder = {}
    for sub in dirs:
        m = re.match(r"E([0-9.]+)$", os.path.basename(sub))
        f = os.path.join(sub, "events.hepmc")
        if not m or not os.path.exists(f):
            continue
        h, n_ev, n_cc = histogram([f])
        ladder[float(m.group(1))] = {
            "shape": {k: (v / v.sum()) for k, v in h.items()},
            "n_ev": n_ev, "n_cc": n_cc}
        print(f"    E = {m.group(1):>5s} GeV: {n_ev} events, "
              f"<n_ch> = {sum((i)*h['default'][i] for i in range(NBIN+2)) / n_ev:.2f}")
    if len(ladder) < 2:
        die(f"{d}: found {len(ladder)} ladder point(s); at least two are "
            f"needed to interpolate")
    hist, n_ev, n_cc, meta = fold_ladder(ladder)
    write(gen, hist, n_ev, n_cc, [os.path.relpath(d, BASE)], meta)


def write(gen, hist, n_ev, n_cc, files, meta=None):
    print(f"donut_nch: {gen}, {n_ev} events, "
          f"{100.0*n_cc/max(n_ev,1):.1f}% charged current")
    out = {"generator": gen, "n_events": n_ev, "cc_fraction": n_cc / max(n_ev, 1),
           "nch": list(range(1, NBIN + 1)), "variants": {},
           "files": files, "ladder": meta,
           "definition": {"charged_pids": sorted(CHARGED),
                          "variants": {n: {"p_min_gev": p, "tan_theta_max": t}
                                       for n, p, t in VARIANTS}}}
    for name, _p, _t in VARIANTS:
        h = hist[name]
        drawn = h[1:NBIN + 1]
        out["variants"][name] = {
            "counts": drawn.tolist(),
            "shape": (drawn / drawn.sum()).tolist(),
            "mean": float(sum((i + 1) * drawn[i] for i in range(NBIN))
                          / drawn.sum()),
            "frac_zero": float(h[0] / h.sum()),
            "frac_overflow": float(h[NBIN + 1] / h.sum()),
        }
        v = out["variants"][name]
        print(f"    {name:8s} <n_ch> = {v['mean']:.2f}   "
              f"n_ch = 0: {100*v['frac_zero']:.1f}%,  "
              f"> {NBIN}: {100*v['frac_overflow']:.1f}%")
    p = f"{RESULTS}/donut_nch_{gen}.json"
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"  wrote {os.path.relpath(p, BASE)}")


def run(gen, files):
    if not files:
        die(f"no event files given for {gen}")
    hist, n_ev, n_cc = histogram(files)
    if not n_ev:
        die(f"{gen}: the event files yielded no events")
    write(gen, hist, n_ev, n_cc,
          [os.path.relpath(f, BASE) for f in files])


def merge():
    """The data and every generator in one file, normalised to the data."""
    if not os.path.exists(DATA):
        die(f"missing {DATA} -- run tools/digitise_donut.py")
    data = json.load(open(DATA))
    n_data = sum(data["data"])
    out = {"what": ("DONUT FIG. 9 charged-particle multiplicity against the "
                    "benchmark's generators, every curve normalised to the "
                    "data's own total"),
           "source": data["source"], "nch": data["nch"],
           "data": data["data"], "data_err": data["data_err"],
           "donut_mc": data["mc"], "n_data": n_data, "generators": {}}
    got = sorted(glob.glob(f"{RESULTS}/donut_nch_*.json"))
    for f in got:
        if f.endswith("donut_nch.json"):
            continue
        d = json.load(open(f))
        g = {"label": d["generator"], "n_events": d["n_events"],
             "cc_fraction": d["cc_fraction"], "variants": {}}
        for name, v in d["variants"].items():
            g["variants"][name] = {
                "events": [n_data * s for s in v["shape"]],
                "mean": v["mean"], "frac_zero": v["frac_zero"],
                "frac_overflow": v["frac_overflow"]}
        out["generators"][d["generator"]] = g
    if not out["generators"]:
        die("no per-generator results to merge -- run the generators first")
    # the data's own mean, over the same drawn range, so the comparison of
    # means is like for like
    tot = sum(data["data"])
    out["data_mean"] = sum(n * v for n, v in zip(data["nch"], data["data"])) / tot
    out["donut_mc_mean"] = (sum(n * v for n, v in zip(data["nch"], data["mc"]))
                            / sum(data["mc"]))
    # CHI-SQUARE AGAINST THE DATA'S OWN POISSON ERROR, on the drawn range.
    # It is a SHAPE test: every curve was normalised to the data's total, so
    # one degree of freedom is spent and the count is N - 1.  No theory error
    # enters, which makes it a floor rather than a fit quality -- the point is
    # to order the curves, not to claim a p-value.
    def chi2(pred):
        terms = [((p - v) / e) ** 2
                 for p, v, e in zip(pred, data["data"], data["data_err"])
                 if e > 0]
        return sum(terms), len(terms) - 1
    c, nd = chi2(data["mc"])
    out["donut_mc_chi2"], out["ndf"] = c, nd
    print(f"donut_nch: data <n_ch> = {out['data_mean']:.2f}, "
          f"DONUT MC {out['donut_mc_mean']:.2f} "
          f"(chi2/N = {c/nd:.1f})")
    for name, g in out["generators"].items():
        for vname, v in g["variants"].items():
            v["chi2"] = chi2(v["events"])[0]
    for name, g in out["generators"].items():
        d0 = g["variants"]["default"]
        print(f"  {name:12s} <n_ch> = {d0['mean']:.2f} "
              f"(default acceptance), {g['variants']['raw']['mean']:.2f} raw, "
              f"chi2/N = {d0['chi2']/nd:.1f}")
    p = f"{RESULTS}/donut_nch.json"
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"  wrote {os.path.relpath(p, BASE)}")


def main():
    if "--merge" in sys.argv:
        merge()
        return
    if len(sys.argv) < 3:
        die(__doc__)
    if sys.argv[2] == "--ladder":
        run_ladder(sys.argv[1], sys.argv[3])
        return
    run(sys.argv[1], sys.argv[2:])


if __name__ == "__main__":
    main()
