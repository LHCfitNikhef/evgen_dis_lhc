#!/usr/bin/env python3
"""PDF bands on a HADRON-LEVEL selection, over every PDF member.

WHAT THIS ADDS to `powheg_selection_pdf.py`, which it otherwise follows
exactly.  That module joins the ELEVEN weights carried inline in the
reweighted LHE -- seven scale points and four central members -- so it gives a
set-to-set spread but no BAND.  This one joins the 319 weights harvested into
`weights.npz`, so every member of every set is available and a real PDF
uncertainty can be formed on the dimuon tier.

WHY THE DIMUON TIER NEEDS THE SHOWER AT ALL.  It asks for an opposite-sign
second muon from a charm decay, which does not exist in a Les Houches file.
So the weights (in the LHE) must be joined to the selection (in the HepMC),
and the join is by `lhe_index` -- NOT by position, because Pythia aborts some
events and the showered file is shorter than the LHE.  After the first abort a
positional join would attach every weight to the wrong event and still produce
a perfectly plausible band.

THE PUBLISHED PRODUCTION SAMPLES CARRY NO lhe_index; they predate it.  The
indexed twins live in `powheg/rwgtnu_job_<TAG>_1`, re-showered from the SAME
LHE by `powheg/reshower_v2_indexed.sh`.  The closure below is what shows the
twin reproduces the published sample rather than merely resembling it.

AT 400 GeV THE WEIGHTS COVER A PREFIX of a longer sample -- the full 1M events
cannot be reweighted member-by-member in reasonable time -- so events beyond
the cap are DROPPED, not silently mismatched, and the count that survives is
reported.  The dimuon tier keeps roughly one event in 12000 there, so this is
the tier where the cap actually bites.

Usage:
  BENCH_ENERGY=1000 BENCH_SELECTION=faser_dimuon \\
      analysis/powheg_selection_members.py
Writes results_nu/powheg_selmembers_<selection>_<tag>.json
"""
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import beams                                                # noqa: E402
import paths                                                # noqa: E402
import selection as selection_mod                           # noqa: E402
from analyze import (dis_invariants, lab_energy, lepton_theta,  # noqa: E402
                     parse_hepmc3, n_opposite_sign_muons)
from powheg_selection_pdf import hepmc_indices, LEP_OUT      # noqa: E402

CL68 = 100.0 * math.erf(1.0 / math.sqrt(2.0))
NOMINAL_ID = "1001"
SCALE_IDS = ["1001", "1002", "1003", "1004", "1005", "1006", "1007"]


def sample_for(tag):
    """The INDEXED showered twin for one beam tag.

    `rwgtnu_job_1` is the ORIGINAL 1 TeV twin and must not be used: it was
    showered two minutes before the lhe_index fix landed, so its indices count
    Pythia's retries and the join is scrambled.  Measured, 28% of its
    dimuon-selected events had an outgoing charm quark in the matching LHE
    against 98% for a correctly indexed sample -- and a dimuon in CC DIS IS a
    charm decay, so the right answer is ~100%.  It was also showered from the
    superseded pre-CKM-fix LHE.  The twin to use at every energy is the one
    `reshower_v2_indexed.sh` writes, `rwgtnu_job_<TAG>_1`.
    """
    return f"{BASE}/powheg/rwgtnu_job_{tag}_1/events.hepmc"


def _inclusive_denominator(energy, tag, sel):
    """(per-member inclusive weight sums, nominal tier fraction).

    The sums come from powheg_pdf_members.py, which is the only pass that
    sees every fiducial event with every member weight.  The nominal fraction
    is the ratio of two MEASURED cross-sections -- the tier's and the
    inclusive one, as analyze_nu.py wrote them -- so the absolute scale of
    this figure is a benchmark result rather than a reweighting artefact.

    Returns ({}, None) if either input is missing, and the caller then simply
    does not write a fraction: a figure with a missing panel is better than
    one with a panel built from whatever happened to be on disk.
    """
    suf = "" if tag == beams.Beams("nu", beams.ANCHOR_ENERGY).tag \
        else f"_{tag}"
    p = f"{BASE}/results_nu/powheg_pdf_members{suf}.json"
    if not os.path.exists(p):
        print(f"  no {os.path.basename(p)} -- no inclusive denominator, "
              f"so no fraction is written")
        return {}, None
    with open(p) as f:
        d = json.load(f)
    sums = d.get("fid_sums_by_id") or {}
    if not sums:
        print(f"  {os.path.basename(p)} predates fid_sums_by_id -- rerun "
              f"analysis/powheg_pdf_members.py")
        return {}, None
    a = f"{BASE}/results_nu/histos_powheg_nu_{sel.name}{suf}.json"
    b = f"{BASE}/results_nu/histos_powheg_nu{suf}.json"
    if not (os.path.exists(a) and os.path.exists(b)):
        return sums, None
    with open(a) as f:
        num = json.load(f)
    with open(b) as f:
        den = json.load(f)
    if num.get("stat_insufficient") or not den.get("sigma_fid_pb"):
        return sums, None
    frac = num["sigma_fid_pb"] / den["sigma_fid_pb"]
    print(f"  nominal {sel.name} fraction = {100*frac:.4f}% of the "
          f"inclusive rate")
    return sums, frac


def main():
    energy = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
    tag = beams.Beams("nu", energy).tag
    sel = selection_mod.get()
    sample = sample_for(tag)
    # PREFER THE SELECTED-SUBSET HARVEST when one exists for this selection.
    # It is the same calculation -- events failing the selection contribute
    # nothing to the member spread, proven to 3e-13 by
    # tools/check_subset_reweight.py -- and at 400 GeV it is the only
    # affordable route: the full 1M-event reweighting was measured at 32 h,
    # against seconds for the ~80 events that actually enter the band.
    sub = f"{paths.POWHEG_V2}/rwgtsel-nu{tag}-{sel.name}"
    row_of = None
    if os.path.exists(f"{sub}/weights.npz"):
        npz = f"{sub}/weights.npz"
        with open(f"{sub}/pwgevents-subset.index.json") as f:
            man = json.load(f)
        row_of = {int(li): r for r, li in enumerate(man["lhe_index"])}
        print(f"  using the {sel.name} subset harvest ({man['rows']} events)")
    else:
        npz = f"{paths.POWHEG_V2}/rwgt-nu{tag}/weights.npz"
    for p in (sample, npz):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    z = np.load(npz, allow_pickle=False)
    ids = [str(x) for x in z["ids"]]
    W = z["weights"]
    col = {i: k for k, i in enumerate(ids)}
    # ROW -> lhe_index.  For a full-sample harvest the two coincide, because
    # every pass reads the same pristine LHE in order.  For a SELECTED-SUBSET
    # harvest they do not, and the mapping is read from the manifest rather
    # than reconstructed -- a positional join that has silently slipped is the
    # single failure mode this area keeps producing.
    if row_of is None:
        row_of = {j: j for j in range(W.shape[0])}
    with open(f"{BASE}/powheg/reweight/rwl_members.json") as f:
        rwl = json.load(f)
    print(f"powheg_selection_members: nu {energy:g} GeV [{tag}], "
          f"selection {sel.name}")
    print(f"  weights {W.shape[0]} events x {W.shape[1]}; sample "
          f"{os.path.relpath(sample, BASE)}")

    idx = hepmc_indices(sample)
    nw = W.shape[1]
    sums = np.zeros(nw)
    sum_all_nom = 0.0
    n_ev = n_sel = beyond = 0
    ncol = col[NOMINAL_ID]
    for ievt, (w, k, P, parts, _dm, _hard) in enumerate(
            parse_hepmc3(sample, beam_pid=14)):
        n_ev += 1
        if ievt >= len(idx) or idx[ievt] is None:
            continue
        j = row_of.get(int(idx[ievt]))
        if j is None or j >= W.shape[0]:
            beyond += 1     # outside the reweighted prefix or subset
            continue
        row = W[j]
        # NOT summed over a subset: with a subset harvest this denominator
        # would run over the selected events only.  It is member-independent
        # either way, so the BAND is unaffected -- and the band is all this
        # module reports.  Recorded so nobody later reads it as a cross-section
        # normalisation, which it is not.
        sum_all_nom += row[ncol]
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
        if not sel.passes_lepton(best_elab, lepton_theta(best_elab, best)):
            continue
        if sel.dimuon and not sel.passes_dimuon(
                n_opposite_sign_muons(parts, P, LEP_OUT, sel)):
            continue
        n_sel += 1
        sums += row
    if n_ev != len(idx):
        sys.exit(f"parser disagreement: {n_ev} events vs {len(idx)} indices")
    if beyond:
        # With a SUBSET harvest this is simply every event the selection does
        # not keep, and is expected -- it is not lost data.  With a CAPPED
        # full harvest it is events outside the reweighted prefix, which IS a
        # restriction worth reading.  The two are distinguished so a routine
        # count is not mistaken for a warning.
        why = ("not in the selected subset (expected)" if row_of and
               len(row_of) < n_ev else
               f"beyond the reweighted prefix of {W.shape[0]}")
        print(f"  {beyond} showered event(s) {why}")
    print(f"  {n_sel} of {n_ev} showered events pass {sel.name}")
    if n_sel < 10:
        sys.exit(f"  only {n_sel} selected event(s) -- below the statistics "
                 f"floor (analyze.MIN_SELECTED_EVENTS = 10). A band on this "
                 f"would be the choice of events, not the PDF.")

    def sig(i):
        c = col.get(i)
        return None if c is None or not sum_all_nom else sums[c]/sum_all_nom

    nom = sig(NOMINAL_ID)

    # ---- the INCLUSIVE denominator, member by member -----------------------
    # THE TIER'S RATE AS A FRACTION OF THE INCLUSIVE ONE (user, 2026-09-01),
    # so this figure reads like "PDF dependence of the CC charm fraction":
    # a fraction on top, the ratio to NNPDF4.0 beneath, both with bands.
    #
    # The denominator CANNOT be formed here.  This module reads a SELECTED
    # SUBSET of the events, so its own "sum over all events" is a sum over
    # the subset.  powheg_pdf_members.py is the pass that sees every fiducial
    # event with every member weight, and it now records that sum per member
    # (fid_sums_by_id), normalised to the nominal.
    #
    # Both sums are divided by their own nominal before the ratio is taken,
    # so every member-independent constant cancels -- the sample sizes, the
    # integrator normalisation, the fact that one pass reads showered events
    # and the other the Les Houches file.  What survives is exactly the
    # member-to-member variation of sigma_tier/sigma_inclusive.  The ABSOLUTE
    # scale then comes from the benchmark's own two cross-sections, which are
    # measured rather than reweighted.
    incl_by_id, frac_nominal = _inclusive_denominator(energy, tag, sel)
    out = {"current": "nu", "energy_gev": energy, "beam_tag": tag,
           "selection": sel.name, "n_selected": int(n_sel),
           "n_showered": int(n_ev), "dropped_beyond_prefix": int(beyond),
           "confidence_level_percent": CL68,
           "note": ("relative to the nominal weight; an absolute "
                    "cross-section is not formed here because only the BAND "
                    "is used, and the band is a ratio on a fixed event set"),
           "sets": {}}
    fs = [sig(i) for i in SCALE_IDS if sig(i)]
    if nom and fs:
        out["scale"] = {"rel_up": max(fs)/nom - 1.0,
                        "rel_down": 1.0 - min(fs)/nom,
                        "rel_half_width": 0.5*(max(fs)-min(fs))/nom}
    import lhapdf
    lhapdf.setVerbosity(0)
    for name, rec in rwl["sets"].items():
        vals = [sig(i) for i in rec["weight_ids"] if i in col]
        vals = [v for v in vals if v]
        if len(vals) < 2:
            continue
        pset = lhapdf.getPDFSet(name)
        try:
            u = pset.uncertainty(list(vals), CL68)
            c, ep, em = u.central, u.errplus, u.errminus
        except Exception as exc:                              # noqa: BLE001
            print(f"    !! LHAPDF uncertainty() failed ({exc})")
            c = vals[0]
            ep = em = math.sqrt(sum((v-c)**2 for v in vals[1:]))
        rec_out = {"rel_to_nominal": c/nom if nom else None,
                   "rel_plus": ep/c if c else None,
                   "rel_minus": em/c if c else None,
                   "members_used": len(vals)}
        # the same members again, as a FRACTION of the inclusive rate, with
        # the correlation between numerator and denominator kept because the
        # ratio is formed member by member before the envelope is taken
        if incl_by_id and frac_nominal:
            fr = []
            for i in rec["weight_ids"]:
                d = incl_by_id.get(str(i))
                v = sig(i)
                if d and v and nom:
                    fr.append((v / nom) / d * frac_nominal)
            if len(fr) >= 2:
                try:
                    uf = pset.uncertainty(list(fr), CL68)
                    fc, fp, fm = uf.central, uf.errplus, uf.errminus
                except Exception:                             # noqa: BLE001
                    fc = fr[0]
                    fp = fm = math.sqrt(sum((v-fc)**2 for v in fr[1:]))
                rec_out.update({"fraction": fc, "fraction_err_plus": fp,
                                "fraction_err_minus": fm,
                                "fraction_members_used": len(fr)})
        out["sets"][name] = rec_out
        print(f"  {name:24s} {c/nom:7.4f} x nominal   "
              f"band +{100*ep/c:5.2f}% -{100*em/c:5.2f}%  "
              f"({len(vals)} members)")
    fn = f"{BASE}/results_nu/powheg_selmembers_{sel.name}_{tag}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {fn}")


if __name__ == "__main__":
    main()
