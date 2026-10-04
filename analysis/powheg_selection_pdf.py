#!/usr/bin/env python3
"""PDF and scale dependence of a HADRON-LEVEL selection, from POWHEG-V2.

WHY THIS EXISTS.  `powheg_nlo_uncertainty.py` reweights the Les Houches file
directly and is therefore limited to selections that cut on the scattered
lepton alone -- it REFUSES Tier E and the dimuon tier by name, because an LHE
has no shower and no hadrons.  But the dimuon tier is exactly where a PDF
study is most interesting: it is a charm tag built on a decay muon, so it
probes the strange and charm densities far more directly than the inclusive
rate does.

HOW THE TWO HALVES ARE JOINED.  The weights live in the LHE; the selection
needs the showered event.  `main_powheg` writes an `lhe_index` attribute on
every HepMC event, exact because every `pythia.next()` consumes exactly one
LHE event whether it succeeds or aborts.  This script reads the weights from
the LHE by index, reads the selection from the HepMC, and joins them.

WHY NOT USE PYTHIA'S PROPAGATED WEIGHTS, which are also in the HepMC.  Because
they are wrong for this file, measured:

  * Pythia exposes only THREE of the four PDF weights -- id 2004 is dropped
    silently, and 2004 is ATLASpdf21, the set furthest from the nominal.
  * Its weight NAMES follow its own canonical (muR, muF) enumeration while the
    VALUES come back in LHE order, so the two are MISMATCHED: it calls the
    second weight MUR1.0_MUF2.0 where the header says 1002 is renscfact=2d0
    facscfact=1d0, i.e. (2, 1).

Reading the LHE with this repository's own parser gets all eleven, correctly
labelled.  The HepMC weights are ignored here on purpose.

THE NORMALISATION IS PINNED TO THE INTEGRATOR, exactly as every other POWHEG
result on this page is:

    sigma_X = sigma_integrator * sum_selected(w_X) / sum_all(w_nominal)

with both sums over the SHOWERED events.  The denominator is the NOMINAL sum
for every variation, not that variation's own sum: dividing each variation by
its own total would renormalise it away and leave only its shape.

The obvious alternative, sum(w)/N_offered, was tried first and is wrong here
by 5%: this sample's delivered cross-section sits 1.8% under the integrator
(its stored `integrator_closure` is 0.9820) and the showered file is missing
the 940 events that aborted.  Pinning to the integrator is what analyze_nu.py
does for `powheg_nu`, and the inclusive selection is checked against that
published result below -- it must reproduce it, since it is the same sample,
the same selection and the same events.

Usage:
  BENCH_SELECTION=faser_dimuon powheg_selection_pdf.py
  powheg_selection_pdf.py --all          # every selection in turn
Writes results_nu/powheg_seldep_<selection>.json
"""
import json
import math
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import beams  # noqa: E402
import paths  # noqa: E402
import selection as selection_mod  # noqa: E402
from analyze import (parse_hepmc3, dis_invariants, lepton_theta, lab_energy,
                     n_opposite_sign_muons, M_P)  # noqa: E402
from powheg_nlo_uncertainty import (SCALE_IDS, SCALE_LABEL, PDF_IDS,
                                    PDF_LABEL, NOMINAL_ID, PDF_NOMINAL_ID,
                                    read_init, parse_events)  # noqa: E402
from analyze_nu import POWHEG_V2_NU_SIGMA_PB  # noqa: E402 -- one source

RESULTS_NU = f"{BASE}/results_nu"
# THE INDEXED TWIN, and only ever the one reshower_v2_indexed.sh writes.
# `powheg/rwgtnu_job_1` -- the original 1 TeV twin -- must NOT be used: it was
# showered two minutes before the lhe_index fix landed, so its indices count
# Pythia's retries, AND it came from the superseded pre-CKM-fix LHE.  With the
# scrambled join only 28% of its dimuon-selected events had an outgoing charm
# quark in the matching LHE event, against 98% for a correct sample -- and a
# dimuon in CC DIS IS a charm decay.
ENERGY = float(os.environ.get("BENCH_ENERGY")
                or beams.ANCHOR_ENERGY)
_TAG = beams.Beams("nu", ENERGY).tag
SAMPLE = f"{BASE}/powheg/rwgtnu_job_{_TAG}_1/events.hepmc"
LEP_OUT = 13          # mu- in the charged current


def npz_weights(path, want):
    """{lhe_index: {weight id: weight}} from the harvested member array.

    WHY NOT THE LHE ANY MORE.  This used to stream the inline <wgt> tags out
    of pwgevents-rwgt.lhe.  `run_reweight_members.sh` deletes that file after
    each of its 54 passes -- it is 60 MB of duplicated event records to carry
    a few numbers -- so once a member run has finished there is no reweighted
    LHE left to read.  weights.npz holds the same numbers for all 319 weights
    rather than the 11 the old inline file carried.
    """
    import numpy as np
    z = np.load(path, allow_pickle=False)
    ids = [str(x) for x in z["ids"]]
    W = z["weights"]
    col = {i: k for k, i in enumerate(ids)}
    missing = [i for i in want if i not in col]
    if missing:
        sys.exit(f"{path} carries no weight id(s) {missing}")
    table = {n: {i: float(W[n, col[i]]) for i in want}
             for n in range(W.shape[0])}
    return table, W.shape[0]


def hepmc_indices(path):
    """The lhe_index of each event in the HepMC, in file order.

    A separate pass rather than a new parameter on parse_hepmc3: that function
    exists in two copies which tools/check_parser_drift.py holds identical, so
    widening it means widening both.  Reading the attribute here keeps the
    shared parser untouched, and the two passes see the same file in the same
    order so a positional zip is exact -- which is asserted below, not assumed.
    """
    out, started, cur = [], False, None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("E "):
                if started:
                    out.append(cur)
                started, cur = True, None
            elif line.startswith("A 0 lhe_index "):
                cur = int(line.split()[3])
    if started:
        out.append(cur)
    return out


# Tier E cuts on the CHARGED HADRONIC FINAL STATE -- multiplicity above two
# thresholds and the azimuthal separation to the nearest charged hadron --
# and that logic lives inside analyze.py's event loop.  Reproducing it here
# would be a SECOND COPY of a selection rule, which is the single thing this
# repository has most reliably got wrong: two copies of one rule drift, and a
# drifted selection is invisible because both halves still produce plausible
# numbers.  So it is refused rather than approximated.  The closure gate found
# this on its first run -- Tier E came out +27.9% high on 26324 events against
# the published 20568, because the hadronic cuts simply were not applied.
UNSUPPORTED = {"faser_e"}


def run(sel):
    if sel.name in UNSUPPORTED:
        sys.exit(f"{sel.name} is not supported here.\n"
                 f"  It cuts on the charged hadronic final state, and that "
                 f"rule lives in analyze.py's\n"
                 f"  event loop.  Copying it into this script would create a "
                 f"second copy of a\n"
                 f"  selection rule, which is exactly how this benchmark has "
                 f"broken before.\n"
                 f"  Supported: inclusive, faser_s, faser_dimuon -- all of "
                 f"which cut on the\n"
                 f"  scattered lepton and the muons alone.")
    if not os.path.exists(SAMPLE):
        sys.exit(f"no showered reweighted sample at {SAMPLE}\n"
                 f"  shower $POWHEG_V2/rwgt-nu1TeV/pwgevents-rwgt.lhe with "
                 f"powheg/main_powheg first")
    npz = f"{paths.POWHEG_V2}/rwgt-nu{_TAG}/weights.npz"
    if not os.path.exists(npz):
        sys.exit(f"no harvested weights at {npz}\n"
                 f"  run powheg/reweight/run_reweight_members.sh first")
    print(f"powheg_selection_pdf: {_TAG}, selection {sel.name}")
    wtab, n_offered = npz_weights(npz, SCALE_IDS + PDF_IDS)
    print(f"  weights: {n_offered} events, "
          f"{len(wtab[0])} of them used here")

    idx = hepmc_indices(SAMPLE)
    all_ids = SCALE_IDS + PDF_IDS
    sums = {i: 0.0 for i in all_ids}
    sum_all_nom = 0.0          # denominator: the NOMINAL sum over all showered
    n_sel = n_ev = 0
    missing_idx = 0

    for ievt, (w, k, P, parts, _dm, _hard) in enumerate(
            parse_hepmc3(SAMPLE, beam_pid=14)):
        n_ev += 1
        # the positional zip is exact only if both passes agree on the event
        # count; a mismatch means one of them skipped an event and every
        # weight after it would belong to the wrong event
        if ievt >= len(idx) or idx[ievt] is None:
            missing_idx += 1
            continue
        wg = wtab.get(idx[ievt])
        if wg is None or any(i not in wg for i in all_ids):
            missing_idx += 1
            continue
        sum_all_nom += wg[NOMINAL_ID]
        best, best_elab = None, -1.0
        for pid, p in parts:
            if pid == LEP_OUT:
                e = lab_energy(p, P)
                if e > best_elab:
                    best, best_elab = p, e
        if best is None:
            continue
        Q2, y, _x, _kP = dis_invariants(k, P, best)
        if Q2 < sel.q2_min or y < sel.y_min or y > sel.y_max:
            continue
        theta = lepton_theta(best_elab, best)
        if not sel.passes_lepton(best_elab, theta):
            continue
        if sel.dimuon and not sel.passes_dimuon(
                n_opposite_sign_muons(parts, P, LEP_OUT, sel)):
            continue
        n_sel += 1
        for i in all_ids:
            sums[i] += wg[i]

    if n_ev != len(idx):
        sys.exit(f"parser disagreement: parse_hepmc3 yielded {n_ev} events "
                 f"but {len(idx)} lhe_index attributes were found.\n"
                 f"  The positional join would be misaligned; refusing.")
    if missing_idx:
        print(f"  !! {missing_idx} event(s) had no usable index/weights")
    if not n_sel:
        sys.exit(f"no events pass {sel.name}")

    sig = {i: POWHEG_V2_NU_SIGMA_PB * sums[i] / sum_all_nom
           for i in all_ids}
    nom = sig[NOMINAL_ID]
    self_check = sig[PDF_NOMINAL_ID] / nom - 1.0
    if abs(self_check) > 1e-6:
        sys.exit(f"weight {PDF_NOMINAL_ID} should reproduce the nominal and "
                 f"is off by {self_check:+.3e} -- the reweighting did not do "
                 f"what the XML asked.")

    print(f"  {n_sel} of {n_ev} showered events pass; "
          f"sigma_fid = {nom:.5g} pb")
    # CLOSURE against the published result for the same sample and selection.
    # Not decoration: it is the only thing that shows the LHE->HepMC join is
    # aligned.  A misaligned join would still produce a plausible number.
    # THE ENERGY SUFFIX.  Without it the 400 GeV and 4 TeV runs compared
    # themselves against the 1 TeV published numbers and the gate reported
    # -42% and +300% -- the same missing-suffix trap this scan has produced
    # half a dozen times.  beams gives the suffix; it is not spelled by hand.
    _esuf = "" if ENERGY == beams.ANCHOR_ENERGY else f"_{_TAG}"
    _ssuf = "" if sel.name == "inclusive" else f"_{sel.name}"
    pub = f"{RESULTS_NU}/histos_powheg_nu{_ssuf}{_esuf}.json"
    if os.path.exists(pub):
        with open(pub) as f:
            ref = json.load(f)
        r = nom / ref["sigma_fid_pb"]
        print(f"  closure vs published {os.path.basename(pub)}: {r:.4f} "
              f"({100*(r-1):+.2f}%; {nom:.4f} vs {ref['sigma_fid_pb']:.4f} pb, "
              f"{n_sel} vs {ref.get('n_fiducial')} events)")
        if abs(r - 1.0) > 0.01:
            sys.exit(f"  closure FAILED: the join or the normalisation is "
                     f"wrong, and a band on top of that would look fine.")
    print("  PDF sets:")
    for i in PDF_IDS:
        print(f"      {PDF_LABEL[i]:>22}  {sig[i]:9.5g} pb  "
              f"{100*(sig[i]/nom - 1):+7.2f}%")
    print("  scale points:")
    for i in SCALE_IDS:
        print(f"      {SCALE_LABEL[i]:>12}  {sig[i]:9.5g} pb  "
              f"{100*(sig[i]/nom - 1):+7.2f}%")

    out = {
        "what": ("PDF and scale dependence of a hadron-level selection, "
                 "POWHEG-V2 nu CC, by joining LHE weights to showered events "
                 "on lhe_index"),
        "generator": "powheg_nu",
        "selection": sel.as_dict(),
        "energy_gev": beams.ANCHOR_ENERGY,
        "n_offered": n_offered, "n_showered": n_ev, "n_selected": n_sel,
        "sigma_fid_pb": nom,
        "pdf": {"ids": PDF_IDS, "labels": [PDF_LABEL[i] for i in PDF_IDS],
                "sigma_fid_pb": [sig[i] for i in PDF_IDS],
                "rel": [sig[i]/nom - 1.0 for i in PDF_IDS]},
        "scale": {"ids": SCALE_IDS,
                  "labels": [SCALE_LABEL[i] for i in SCALE_IDS],
                  "sigma_fid_pb": [sig[i] for i in SCALE_IDS],
                  "rel": [sig[i]/nom - 1.0 for i in SCALE_IDS]},
    }
    os.makedirs(RESULTS_NU, exist_ok=True)
    p = f"{RESULTS_NU}/powheg_seldep_{sel.name}.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"  wrote {os.path.relpath(p, BASE)}")


def main():
    if "--all" in sys.argv:
        for name in ("inclusive", "faser_s", "faser_dimuon"):
            os.environ["BENCH_SELECTION"] = name
            import importlib
            importlib.reload(selection_mod)
            run(selection_mod.get())
            print()
        return
    run(selection_mod.get())


if __name__ == "__main__":
    main()
