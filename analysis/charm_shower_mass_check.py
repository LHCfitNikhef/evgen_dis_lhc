#!/usr/bin/env python3
"""N4: does Sherpa's excess charm come from a massless charm in its shower?

The charm-tag anatomy splits every generator's hadron-level charm into the
hard process, the beam remnant and shower g -> c cbar.  Sherpa's shower piece
is an order of magnitude above Pythia's and Herwig's on both currents, and it
is the piece that grows with beam energy.  Its logs say why: "Massive PS
flavours: (none)" -- the benchmark's RESPECT_MASSIVE_FLAG: true, which makes
the ME massless as the convention requires, makes the SHOWER massless too,
so g -> c cbar has no charm-mass threshold there.

This script compares the same LO card with MASSIVE_PS: [4, 5]
(sherpa/run_sherpa_massive_ps.sh) against the published LO sample, on both
currents, through the anatomy's own event pass, and reports the shower
piece of each.  It also reports the price of the diagnostic: with a massive
shower charm the ME -> shower interface discards events whose incoming
charm it cannot put on shell (the trap of 2026-08-13), so the delivered
cross-section falls below the integrator's and the sample is depleted in
charm-INITIATED events.  The shower piece, measured on events whose hard
process has no charm at all, is not touched by that.

Usage: analysis/charm_shower_mass_check.py [--max N]
Writes results_nu/charm_shower_mass_check.json and prints the table.
"""
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import charm_tag_anatomy as cta  # noqa: E402
import paths  # noqa: E402

RUNS = paths.SHERPA_RUNS
ARMS = {"mu": [("MuonDIS_LO", "LO, massless shower (benchmark)"),
               ("MuonDIS_LO_massivePS", "LO, MASSIVE_PS: [4, 5]")],
        "nu": [("NuDIS_LO", "LO, massless shower (benchmark)"),
               ("NuDIS_LO_massivePS", "LO, MASSIVE_PS: [4, 5]")],
        # the same at MC@NLO: tests whether the massless first emission of
        # the S-events is what lifts Sherpa's HARD-process charm at NLO
        "mu_nlo": [("MuonDIS_NLO", "MC@NLO, massless shower (benchmark)"),
                   ("MuonDIS_NLO_massivePS", "MC@NLO, MASSIVE_PS: [4, 5]")],
        "nu_nlo": [("NuDIS_NLO_ckm3", "MC@NLO, massless shower (benchmark)"),
                   ("NuDIS_NLO_ckm3_massivePS", "MC@NLO, MASSIVE_PS: [4, 5]")]}
PID = {"mu": 13, "nu": 14, "mu_nlo": 13, "nu_nlo": 14}


def discards(rundir):
    """Events the ME -> shower interface threw away, from the job logs."""
    n_fail = n_tot = 0
    for log in glob.glob(f"{rundir}/job_*/sherpa.log"):
        for line in open(log, errors="replace"):
            m = re.search(r'From "Jet_Evolution:CSS": ~?(\d+) \(~?(\d+)\)', line)
            if m:
                n_fail += int(m.group(1))
                n_tot += int(m.group(2))
    return n_fail, n_tot


def main():
    max_events = None
    if "--max" in sys.argv:
        max_events = int(sys.argv[sys.argv.index("--max") + 1])
    out = {}
    which = [a for a in sys.argv[1:] if a in ARMS] or ["mu", "nu"]
    prev = {}
    if os.path.exists(f"{BASE}/results_nu/charm_shower_mass_check.json"):
        with open(f"{BASE}/results_nu/charm_shower_mass_check.json") as f:
            prev = json.load(f)
    out.update(prev)
    for cur in which:
        out[cur] = {}
        for run, label in ARMS[cur]:
            files = sorted(glob.glob(f"{RUNS}/{run}/job_*/evtfull"))
            if not files:
                print(f"{cur} {run}: no events yet")
                continue
            res = cta.summarise(*cta.anatomy(files, PID[cur], max_events))
            nf, nt = discards(f"{RUNS}/{run}")
            s = res["split"]
            out[cur][run] = {"label": label, "files": len(files),
                             "n_events": res["n_events"],
                             "frac_final": res["frac_final"],
                             "hard_out": s["hard_out"], "remnant": s["remnant"],
                             "shower": s["shower"],
                             "shower_over_nonhard": s["shower"] / (1 - res["frac_hard_any"]),
                             "discarded": nf, "discarded_of": nt,
                             "E_lead_shower_median": res.get("E_lead_shower_median")}
            print(f"{cur} {run:22s} [{label}] {res['n_events']} events: "
                  f"final {100*res['frac_final']:.2f}% = hard {100*s['hard_out']:.2f} "
                  f"+ remnant {100*s['remnant']:.2f} + shower {100*s['shower']:.2f}; "
                  f"shower per non-charm hard event {100*s['shower']/(1-res['frac_hard_any']):.3f}%; "
                  f"interface discards {nf} of {nt}")
        a = out[cur]
        if len(a) == 2:
            k0, k1 = [r for r, _ in ARMS[cur]]
            if k0 in a and k1 in a:
                r = a[k1]["shower_over_nonhard"] / a[k0]["shower_over_nonhard"]
                out[cur]["shower_ratio_massive_over_massless"] = r
                print(f"{cur}: shower charm per non-charm hard event, massive / massless = {r:.3f}")
    with open(f"{BASE}/results_nu/charm_shower_mass_check.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote results_nu/charm_shower_mass_check.json")


if __name__ == "__main__":
    main()
