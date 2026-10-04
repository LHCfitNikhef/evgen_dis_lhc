#!/usr/bin/env python3
"""Turn an MG5_aMC LHE file into the benchmark's histos_*.json schema.

MG5 is a parton-level (matrix-element) generator here -- no shower, no
hadronisation -- so it belongs with the ME-LEVEL samples (Sherpa LO ME,
Herwig LO ME, Pythia LO ME) and only the six DIS-kinematics observables are
meaningful.  nch/nch1/Elead/dphi/D-mesons are left out entirely rather than
written as zeros, which would read as "no hadrons produced" instead of "this
sample cannot express that observable".

The observables use analyze.py's definitions verbatim, with one wrinkle: an
LHE records the incoming PARTON, not the beam proton, so the proton momentum
is rebuilt from the run card's beam energy along -z.  Every quantity below is
then the same invariant the HepMC analyses compute.

usage: mg5_lhe_histos.py <run_dir> <out_tag> [--nu]
       run_dir = an MG5 process directory (its Events/*/unweighted_events.lhe.gz
       and the cross-section from the run are read automatically)
"""
import glob
import gzip
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from analyze import BINS, M_MU, M_P, Q2_MIN, Y_MAX, Y_MIN, dot  # noqa: E402
from runmeta import stamp  # noqa: E402  -- the one (selection, energy) stamp

OBS = ("Q2", "xbj", "y", "Emu", "theta", "nu")


def read_beam_energies(run_dir):
    """(E1, E2) from the run card, in GeV."""
    card = f"{run_dir}/Cards/run_card.dat"
    e = {}
    for line in open(card):
        m = re.match(r"\s*([0-9.eEdD+-]+)\s*=\s*(ebeam[12])\b", line)
        if m:
            e[m.group(2)] = float(m.group(1))
    return e["ebeam1"], e["ebeam2"]


def read_xsec_pb(evt_dir):
    """Cross-section in pb, from THIS run's banner.

    Must come from the same run directory as the LHE: a process dir
    accumulates run_01, run_02, ... and globbing across all of them silently
    pairs new events with an old cross-section (which is exactly what
    happened here -- the pre-scale-fix 19026 pb against post-fix events).
    """
    for fn in sorted(glob.glob(f"{evt_dir}/*banner*.txt")):
        txt = open(fn, errors="ignore").read()
        m = re.search(r"Integrated weight \(pb\)\s*:\s*([0-9.eE+-]+)", txt)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"no banner cross-section in {evt_dir}")


def events(lhe, lep_in, lep_out, e2):
    """Yield (Q2, xbj, y, Emu, theta, nu) per generated event, or None.

    EVERY generated event is yielded, cut or not, and None marks one whose
    lepton legs could not be read.  The fiducial cut is applied by the caller
    because the GENERATION REGION IS WIDER than the fiducial one on purpose
    (mg5/install_hooks.py: Pythia moves the reconstructed Q2 slightly when it
    puts the Les Houches particles on their mass shells, and a sample
    generated exactly on the boundary can only lose events across it).  So
    the fiducial fraction here is a number below one that MUST divide the
    generated cross-section -- it is not the survival check it used to be.
    """
    # proton beam along -z, rebuilt from the run card (the LHE records the
    # incoming PARTON, not the beam hadron)
    pz = math.sqrt(max(e2 * e2 - M_P * M_P, 0.0))
    P = (e2, 0.0, 0.0, -pz)
    op = gzip.open if lhe.endswith(".gz") else open
    inev, parts = False, []
    for line in op(lhe, "rt"):
        t = line.strip()
        if t.startswith("<event"):
            inev, parts = True, []
            continue
        if t.startswith("</event"):
            inev = False
            k = next((p for p in parts if p[1] == -1 and p[0] == lep_in), None)
            kp = [p for p in parts if p[1] == 1 and p[0] == lep_out]
            if k is None or not kp:
                yield None
                continue
            # hardest outgoing lepton, by invariant lab energy
            best = max(kp, key=lambda p: dot(p[2:], P) / M_P)
            kv, bv = k[2:], best[2:]
            q = tuple(kv[i] - bv[i] for i in range(4))
            Q2 = -dot(q, q)
            kP, Pq = dot(kv, P), dot(P, q)
            if kP <= 0 or Pq <= 0:
                yield None
                continue
            y = Pq / kP
            e_lab = dot(bv, P) / M_P
            pt = math.hypot(bv[1], bv[2])
            pmag = math.sqrt(max(e_lab * e_lab - M_MU * M_MU, 1e-12))
            yield (Q2, Q2 / (2.0 * Pq), y, e_lab,
                   math.asin(min(pt / pmag, 1.0)), kP / M_P - e_lab)
            continue
        if inev and t and not t.startswith("<"):
            f = t.split()
            if len(f) >= 13:
                try:
                    parts.append((int(f[0]), int(f[1]), float(f[9]),
                                  float(f[6]), float(f[7]), float(f[8])))
                except ValueError:
                    pass


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    run_dir, tag = sys.argv[1], sys.argv[2]
    is_nu = "--nu" in sys.argv
    lep_in, lep_out = (14, 13) if is_nu else (13, 13)
    outdir = f"{BASE}/{'results_nu' if is_nu else 'results'}"

    _e1, e2 = read_beam_energies(run_dir)
    lhe = sorted(glob.glob(f"{run_dir}/Events/*/unweighted_events.lhe.gz"))
    if not lhe:
        sys.exit(f"no unweighted_events.lhe.gz under {run_dir}/Events")
    lhe = lhe[-1]                      # newest run
    sigma_pb = read_xsec_pb(os.path.dirname(lhe))
    print(f"  {os.path.basename(os.path.dirname(lhe))}: "
          f"sigma_gen = {sigma_pb:.4f} pb")

    cols = {k: [] for k in OBS}
    n_all = n_reco = 0
    for row in events(lhe, lep_in, lep_out, e2):
        n_all += 1
        if row is None:
            continue
        n_reco += 1
        Q2, _x, y = row[0], row[1], row[2]
        if Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX:
            continue
        for k, v in zip(OBS, row):
            cols[k].append(v)
    n_fid = len(cols["Q2"])
    if not n_fid:
        sys.exit("no events survived the fiducial cut")

    # TWO FRACTIONS, AND THEY MEAN DIFFERENT THINGS.
    #
    # RECONSTRUCTION (n_reco / n_all) must be 1: every generated event has an
    # incoming and an outgoing lepton by construction, so anything missing is
    # a parsing failure, not physics.  It is checked because nothing else
    # would notice -- the same shape as the Sherpa RESPECT_MASSIVE_FLAG
    # discard, which threw away 23% of a sample in silence (CONVENTIONS.md rule 2).
    #
    # FIDUCIAL (n_fid / n_all) is BELOW one on purpose: the generation region
    # is deliberately wider than the fiducial one, so the generated
    # cross-section has to be scaled by this fraction.  It used to be pinned
    # instead, which was right only while the two regions coincided.
    reco = n_reco / n_all if n_all else 0.0
    if reco < 0.99:
        msg = (f"{tag}: only {n_reco} of {n_all} generated events could be "
               f"reconstructed ({100*reco:.2f}%).  Every MG5 event has both "
               f"lepton legs, so this is a parsing failure.")
        if os.environ.get("ALLOW_BAD_CLOSURE") != "1":
            sys.exit(msg + "\n  Set ALLOW_BAD_CLOSURE=1 to analyse anyway.")
        print("WARNING: " + msg)
    fid = n_fid / n_all if n_all else 0.0
    print(f"  {tag}: {100*reco:.2f}% reconstructed, "
          f"{100*fid:.2f}% of them inside the fiducial region")
    sigma_fid_pb = sigma_pb * fid
    out = {"generator": tag, "label": f"MG5_aMC 3.7.2 LO (ME{', nu CC' if is_nu else ''})",
           "n_parsed": n_all, "n_fiducial": n_fid,
           "sigma_gen_pb": sigma_pb,
           "sigma_fid_pb": sigma_fid_pb,
           "sigma_fid_err_pb": sigma_fid_pb / math.sqrt(n_fid),
           "means": {k: float(np.mean(cols[k])) for k in ("Q2", "xbj", "y")},
           "hists": {}}
    norm = sigma_fid_pb / n_fid
    for key in OBS:
        edges = np.asarray(BINS[key], dtype=float)
        cnt, _ = np.histogram(np.asarray(cols[key]), bins=edges)
        widths = np.diff(edges)
        out["hists"][key] = {
            "edges": edges.tolist(),
            "dsig": (cnt * norm / widths).tolist(),
            "err": (np.sqrt(cnt) * norm / widths).tolist()}

    os.makedirs(outdir, exist_ok=True)
    ofn = f"{outdir}/histos_{tag}.json"
    with open(ofn, "w") as f:
        stamp(out, current="mu")
        json.dump(out, f)
    print(f"{tag}: sigma_fid = {sigma_fid_pb:.4f} pb "
          f"(sigma_gen {sigma_pb:.4f} pb) from {n_fid} events -> {ofn}")


if __name__ == "__main__":
    main()
