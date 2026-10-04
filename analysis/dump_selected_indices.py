#!/usr/bin/env python3
"""The lhe_index values that pass a selection, for exact subset reweighting.

Pairs with `powheg/reweight/subset_lhe.py`.  The PDF band on a selection
depends only on the events that pass it -- every other event contributes
identically to numerator and denominator and cancels -- so the reweighter
only ever needs those.  This writes the list.

The selection rule is NOT reimplemented here: it comes from `selection.py`
and the same helpers `analyze.py` uses, because a second copy of a selection
rule is how this benchmark has broken before.

Usage:
  BENCH_ENERGY=400 BENCH_SELECTION=faser_dimuon \\
      analysis/dump_selected_indices.py [out.json]
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import beams                                                # noqa: E402
import selection as selection_mod                           # noqa: E402
from analyze import (parse_hepmc3, dis_invariants, lab_energy,  # noqa: E402
                     lepton_theta, n_opposite_sign_muons)
from powheg_selection_pdf import hepmc_indices, LEP_OUT      # noqa: E402
from powheg_selection_members import sample_for              # noqa: E402


def main():
    energy = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
    tag = beams.Beams("nu", energy).tag
    sel = selection_mod.get()
    sample = sample_for(tag)
    if not os.path.exists(sample):
        sys.exit(f"no indexed twin at {sample} -- run "
                 f"powheg/reshower_v2_indexed.sh {energy:g} first")
    out = sys.argv[1] if len(sys.argv) > 1 else \
        f"{BASE}/results_nu/selected_indices_{sel.name}_{tag}.json"

    idx = hepmc_indices(sample)
    keep, n_ev = [], 0
    for ievt, (w, k, P, parts, _dm, _hard) in enumerate(
            parse_hepmc3(sample, beam_pid=14)):
        n_ev += 1
        if ievt >= len(idx) or idx[ievt] is None:
            continue
        best, be = None, -1.0
        for pid, p in parts:
            if pid == LEP_OUT:
                e = lab_energy(p, P)
                if e > be:
                    best, be = p, e
        if best is None:
            continue
        Q2, y, _x, _kP = dis_invariants(k, P, best)
        if Q2 < sel.q2_min or y < sel.y_min or y > sel.y_max:
            continue
        if not sel.passes_lepton(be, lepton_theta(be, best)):
            continue
        if sel.dimuon and not sel.passes_dimuon(
                n_opposite_sign_muons(parts, P, LEP_OUT, sel)):
            continue
        keep.append(int(idx[ievt]))
    if n_ev != len(idx):
        sys.exit(f"parser disagreement: {n_ev} events vs {len(idx)} indices")
    with open(out, "w") as f:
        json.dump({"current": "nu", "energy_gev": energy, "beam_tag": tag,
                   "selection": sel.name, "n_showered": n_ev,
                   "n_selected": len(keep), "lhe_index": sorted(keep)}, f)
    print(f"  {len(keep)} of {n_ev} showered events pass {sel.name}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
