#!/usr/bin/env python3
"""NLO scale and PDF uncertainties on sigma_fid, from POWHEG's LHE weights.

WHAT THIS ADDS THAT NOTHING ELSE HERE GIVES.  The benchmark's uncertainty
machinery is otherwise analytic: YADISM supplies PDF bands
(`analysis/pdf_dependence.py`) and a 7-point MHOU (`analysis/mhou.py`), both
nearly free because the YADISM operator is independent of the PDF and of the
scale.  What YADISM cannot give is either of those for a MATCHED NLO+PS
sample under a fiducial selection.  POWHEG can, because `storeinfo_rwgt 1`
keeps enough per event to re-evaluate the matrix element at a different scale
or with a different PDF, and `powheg/reweight/` has already appended the
results as extra weights.  The hard process is never rerun.

WHY THIS READS THE LHE AND NOT THE HepMC.  The weights are appended to the
LHE.  The showered HepMC does not carry them -- Pythia keeps one weight per
event -- so consuming them means either re-showering every variation (11x the
generation, and 11 samples to store) or reading the LHE directly.  This reads
the LHE.

WHICH SELECTIONS THAT ALLOWS, AND WHICH IT DOES NOT.  An LHE has no shower and
no hadrons, so only cuts on the SCATTERED LEPTON are defined:

    inclusive   Q2 > 4, 0.2 < y < 0.9                        -> supported
    faser_s     + E' > 100 GeV, theta < 25 mrad              -> supported
    faser_e     + charged-hadron multiplicity, Delta phi     -> REFUSED
    faser_dimuon+ an opposite-sign muon from charm decay     -> REFUSED

and the two supported ones are not an approximation.  The benchmark
convention switches QED radiation OFF entirely, so the shower does not touch
the scattered lepton: its Q2, y, energy and angle are the same before and
after showering.  The two refused tiers are refused in place, with the reason,
rather than silently returning a number from a partial final state --
`sigma_fid` under Tier E depends on the hadron multiplicity, which is the
whole point of that tier.

THE FRAME.  POWHEG-V2 generates a symmetric MASSLESS pair at sqrt(s)/2, so the
LHE is not in the lab frame and its records carry no proton at all -- only the
struck parton.  The beams are rebuilt here from the `<init>` energies with the
particles put ON SHELL, which is what Pythia does when it reads the same file,
and it reproduces 2 k.P = 2 E m_p to a part in 10^6.  Every kinematic quantity
then comes from analyze.py's own helpers, so this file cannot drift away from
the fiducial region the rest of the benchmark uses.

AND IT CLOSES AGAINST THE SHOWERED RESULT.  The nominal weight here must
reproduce the published sigma_fid, which was measured on the HepMC by a
different code path over a different file.  That comparison is computed and
stored on every run, and a disagreement beyond CLOSURE_TOL is fatal.  Without
it the frame reconstruction above would be an untested assumption, and a
plausible-looking band is exactly the kind of result CONVENTIONS.md rule 2 exists
to distrust.

Usage:
  analysis/powheg_nlo_uncertainty.py [--current mu|nu]
  BENCH_ENERGY=400 BENCH_SELECTION=faser_s analysis/powheg_nlo_uncertainty.py
"""
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import beams                                              # noqa: E402
import paths                                              # noqa: E402
import selection as selection_mod                         # noqa: E402
from analyze import (M_P, dis_invariants, lab_energy,      # noqa: E402
                     lepton_theta)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The nominal must come back equal to itself through a different code path, a
# different file and a different frame reconstruction.  1% is loose enough to
# absorb the shower dropout (the HepMC is missing the events Pythia discarded,
# which are weight-neutral but not exactly so at finite statistics) and tight
# enough that a wrong beam axis or a mis-signed y could not hide under it.
CLOSURE_TOL = 0.01

# The weight ids written by powheg/reweight/rwl_scale_pdf.xml.  Stated here as
# well so a silently truncated <rwgt> block is an error rather than a smaller
# envelope: POWHEG writes the weights it managed to compute, and a band built
# from five of seven points is narrower and looks fine.
SCALE_IDS = ["1001", "1002", "1003", "1004", "1005", "1006", "1007"]
SCALE_LABEL = {"1001": "(1, 1)", "1002": "(2, 1)", "1003": "(1/2, 1)",
               "1004": "(1, 2)", "1005": "(1, 1/2)", "1006": "(2, 2)",
               "1007": "(1/2, 1/2)"}
PDF_IDS = ["2001", "2002", "2003", "2004"]
PDF_LABEL = {"2001": "NNPDF40_nnlo_as_01180", "2002": "CT18NNLO",
             "2003": "MSHT20nnlo_as118", "2004": "ATLASpdf21_T1"}
NOMINAL_ID = "1001"
# 2001 re-evaluates the event with the SAME PDF it was generated with, so it
# must return the nominal.  It is the check that the PDF key was honoured at
# all: `rwl_setup_params_weights` tests only `lhapdf`, and the production
# card's own spelling `lhans1` is accepted by the XML parser and then ignored,
# giving four weights that all equal the nominal while the scale weights vary
# correctly beside them.  See powheg/reweight/README.md.
PDF_NOMINAL_ID = "2001"
PDF_SELFCHECK_TOL = 1e-6

SUPPORTED_SELECTIONS = ("inclusive", "faser_s")

EVENT_RE = re.compile(r"<event", re.I)
WGT_RE = re.compile(r"<wgt\s+id=['\"]([^'\"]+)['\"]\s*>\s*([0-9eEdD.+-]+)\s*<")


def die(msg):
    sys.exit(f"powheg_nlo_uncertainty: {msg}")


def lhe_path(current, energy):
    """The reweighted LHE for one (current, energy), or a reason it is absent."""
    tag = beams.Beams(current, energy).tag
    root = paths.POWHEG_V2
    d = f"{root}/rwgt-{current}{tag}"
    f = f"{d}/pwgevents-rwgt.lhe"
    if not os.path.exists(f):
        alt = f"{d}/pwgevents.lhe"
        if os.path.exists(alt):
            f = alt
        else:
            die(f"no reweighted LHE at {d}\n"
                f"  run: powheg/reweight/run_reweight.sh --current {current} "
                f"{energy:g}")
    return f


def read_init(fh, current):
    """Beam PDG ids and energies from <init>, with the beams put on shell.

    POWHEG writes a symmetric massless pair, so taking the record at face
    value would put the proton at pz = -E with no mass.  Pythia restores the
    masses when it reads the same file (the HepMC beam records carry
    m = 0.93827 and 0.10566), and the fiducial cuts are specified in the
    proton rest frame, so the masses have to come back here too.
    """
    ids = energies = None
    for line in fh:
        if "<init>" in line:
            first = next(fh).split()
            ids = (int(first[0]), int(first[1]))
            energies = (float(first[2]), float(first[3]))
            break
    if ids is None:
        die("no <init> block -- is this an LHE?")
    if ids[1] != 2212:
        die(f"second beam is {ids[1]}, expected a proton (2212)")
    m_lep = beams.LEPTONS[current][1]
    e_lep, e_p = energies
    k = (e_lep, 0.0, 0.0, math.sqrt(max(e_lep*e_lep - m_lep*m_lep, 0.0)))
    P = (e_p, 0.0, 0.0, -math.sqrt(max(e_p*e_p - M_P*M_P, 0.0)))
    return ids, k, P


def parse_events(fh, lep_out):
    """Yield (lepton four-momentum, {weight id: weight}, nominal, has_charm).

    Streams: the reweighted files are ~60 MB and there is no reason to hold
    one in memory.

    THE CHARM TAG IS PARTON LEVEL, and that is the right level here rather
    than a compromise (user, 2026-08-27).  An LHE has no shower and no
    hadrons, so the hadron-level `_charmfinal` definition used elsewhere on the
    page cannot be evaluated on it.  But the quantity this study is compared
    against -- YADISM's F2_charm -- is itself a parton-level object, so an
    outgoing charm QUARK is the like-for-like match to it, and arguably a
    closer one than the hadron-level tag would be.

    "Outgoing charm" means status +1 with |pid| = 4, which at this order picks
    up both the Born s -> c and the real-emission configurations.
    """
    in_event = False
    npart = 0
    read = 0
    lep = None
    lep_e3 = -1.0
    wgts = {}
    nominal = None
    has_charm = False
    for line in fh:
        s = line.strip()
        if not in_event:
            if EVENT_RE.match(s):
                in_event, npart, read = True, 0, 0
                lep, lep_e3, wgts, nominal = None, -1.0, {}, None
                has_charm = False
            continue
        if s.startswith("</event>"):
            in_event = False
            if lep is not None:
                yield lep, wgts, nominal, has_charm
            continue
        if npart == 0 and not s.startswith("#") and not s.startswith("<"):
            f = s.split()
            if len(f) >= 6:
                npart = int(f[0])
                nominal = float(f[2])
            continue
        if read < npart and not s.startswith("#") and not s.startswith("<"):
            f = s.split()
            if len(f) >= 13:
                read += 1
                pid, status = int(f[0]), int(f[1])
                if status == 1 and abs(pid) == 4:
                    has_charm = True
                if status == 1 and pid == lep_out:
                    # px, py, pz, E -> the (E, px, py, pz) order dot() wants
                    px, py, pz, e = (float(f[6]), float(f[7]),
                                     float(f[8]), float(f[9]))
                    if e > lep_e3:
                        lep, lep_e3 = (e, px, py, pz), e
            continue
        m = WGT_RE.search(s)
        if m:
            wgts[m.group(1)] = float(m.group(2).replace("D", "E")
                                     .replace("d", "e"))


def hepmc_reference(current, energy, sel_name):
    """The published sigma_fid for the same sample, measured on the HepMC.

    Returns (value_pb, path) or (None, reason).  The result files are the
    tracked half of this repository, so this is reading a committed number,
    not recomputing one.
    """
    gen = "powheg_nu" if current == "nu" else "powheg_v2"
    out = "results_nu" if current == "nu" else "results"
    parts = [gen]
    if sel_name != "inclusive":
        parts.append(sel_name)
    tag = beams.Beams(current, energy).tag
    if energy != beams.ANCHOR_ENERGY:
        parts.append(tag)
    f = f"{BASE}/{out}/histos_{'_'.join(parts)}.json"
    if not os.path.exists(f):
        return None, f"no showered result at {os.path.relpath(f, BASE)}"
    return json.load(open(f))["sigma_fid_pb"], f


def main():
    argv = sys.argv[1:]
    current = "nu"
    if argv and argv[0] == "--current":
        if len(argv) < 2 or argv[1] not in ("mu", "nu"):
            die("--current takes mu or nu")
        current, argv = argv[1], argv[2:]
    energy = float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)
    sel = selection_mod.get()

    if sel.name not in SUPPORTED_SELECTIONS:
        die(f"selection {sel.name!r} is not defined on an LHE.\n"
            f"  An LHE has no shower and no hadrons, and {sel.name} cuts on "
            f"the hadronic final state\n"
            f"  ({'a second muon from charm decay' if sel.dimuon else 'charged-hadron multiplicity and Delta phi'}), "
            f"which does not exist before showering.\n"
            f"  Supported: {', '.join(SUPPORTED_SELECTIONS)} -- both cut on "
            f"the scattered lepton only,\n"
            f"  and the benchmark's no-QED-radiation convention makes those "
            f"shower-invariant.")

    lep_out = 13          # mu- in both currents
    f = lhe_path(current, energy)
    print(f"powheg_nlo_uncertainty: {current} {energy:g} GeV, "
          f"selection {sel.name}")
    print(f"  reading {f} ({os.path.getsize(f)/1e6:.1f} MB)")

    with open(f, errors="replace") as fh:
        ids, k, P = read_init(fh, current)
        e_lab_beam = (k[0]*P[0] - k[3]*P[3]) / M_P
        print(f"  beams {ids[0]} + {ids[1]}, rebuilt on shell -> "
              f"E_lab = {e_lab_beam:.4f} GeV")
        if abs(e_lab_beam/energy - 1.0) > 1e-3:
            die(f"rebuilt beam energy {e_lab_beam:.4f} GeV disagrees with the "
                f"requested {energy:g} GeV")

        all_ids = SCALE_IDS + PDF_IDS
        sum_fid = {i: 0.0 for i in all_ids}
        sum_fid2 = {i: 0.0 for i in all_ids}
        sum_dif2 = {i: 0.0 for i in all_ids}
        sum_all = {i: 0.0 for i in all_ids}
        # the charm-tagged subset, summed with the SAME weights so the
        # fraction can be formed per weight id -- the fraction's variation is
        # not obtainable from the two variations separately, since numerator
        # and denominator move together under a PDF or scale change
        sum_charm = {i: 0.0 for i in all_ids}
        sum_charm_dif2 = {i: 0.0 for i in all_ids}
        n_events = n_fid = n_charm = 0
        n_missing = 0
        for lep, wgts, nominal, has_charm in parse_events(fh, lep_out):
            n_events += 1
            missing = [i for i in all_ids if i not in wgts]
            if missing:
                n_missing += 1
                if n_missing <= 3:
                    print(f"  !! event {n_events} is missing weights "
                          f"{missing[:5]}")
                continue
            for i in all_ids:
                sum_all[i] += wgts[i]
            Q2, y, xbj, kP = dis_invariants(k, P, lep)
            if Q2 < sel.q2_min or y < sel.y_min or y > sel.y_max:
                continue
            e_lab = lab_energy(lep, P)
            if not sel.passes_lepton(e_lab, lepton_theta(e_lab, lep)):
                continue
            n_fid += 1
            if has_charm:
                n_charm += 1
            wnom = wgts[NOMINAL_ID]
            for i in all_ids:
                w = wgts[i]
                sum_fid[i] += w
                sum_fid2[i] += w*w
                # The SHIFT against the nominal, event by event.  Every
                # variation is evaluated on the SAME events as the nominal, so
                # sigma^X and sigma^nom are almost perfectly correlated and
                # the error on their DIFFERENCE is far smaller than either
                # one's own.  Quoting the raw MC error instead would say a
                # 50000-event sample cannot resolve a 1% envelope (0.56% per
                # point, so a 2% band would be "2 sigma"), which is wrong: the
                # band is a ratio on a fixed event set.
                d = w - wnom
                sum_dif2[i] += d*d
            if has_charm:
                for i in all_ids:
                    sum_charm[i] += wgts[i]
                    dc = wgts[i] - wnom
                    sum_charm_dif2[i] += dc*dc

    if n_missing:
        die(f"{n_missing} of {n_events} events carry an incomplete <rwgt> "
            f"block.\n  A partial block would silently NARROW the envelope, "
            f"so this is fatal rather than skipped.")
    if not n_events:
        die("no events parsed")

    # POWHEG's XWGTUP is absolute: sigma = <w> over the file, so no ratio
    # estimator is needed here and each variation normalises itself.
    def sigma(sums, i):
        return sums[i] / n_events

    sig_fid = {i: sigma(sum_fid, i) for i in all_ids}
    err_fid = {i: math.sqrt(sum_fid2[i]) / n_events for i in all_ids}
    # MC error on the SHIFT against the nominal, in pb.  Divided by the
    # nominal below to give a relative figure; keeping it in pb here means the
    # units cannot drift between the print and the stored value.
    err_shift = {i: math.sqrt(sum_dif2[i]) / n_events for i in all_ids}
    sig_tot = {i: sigma(sum_all, i) for i in all_ids}
    nom = sig_fid[NOMINAL_ID]
    # The charm-tagged piece and, per weight id, the FRACTION.  Formed id by
    # id rather than from two separate bands: numerator and denominator share
    # the same quark densities and move together under a PDF or scale change,
    # so a fraction built from independent envelopes would be roughly twice
    # too wide.  Same argument as pdf_dependence.py's member-by-member ratio.
    sig_charm = {i: sigma(sum_charm, i) for i in all_ids}
    frac = {i: (sig_charm[i] / sig_fid[i]) if sig_fid[i] else float("nan")
            for i in all_ids}
    nom_charm, nom_frac = sig_charm[NOMINAL_ID], frac[NOMINAL_ID]

    # The PDF self-check: 2001 is the generation PDF re-evaluated, so it must
    # return the nominal exactly.  If it does not, the lhapdf key was ignored
    # and every PDF weight is the nominal in disguise.
    self_check = sig_fid[PDF_NOMINAL_ID] / nom - 1.0
    if abs(self_check) > PDF_SELFCHECK_TOL:
        die(f"weight {PDF_NOMINAL_ID} should reproduce the nominal and is off "
            f"by {self_check:+.3e}.\n  That weight re-evaluates the event with "
            f"the PDF it was generated with, so a difference means the "
            f"reweighting did not do what the XML asked.")
    spread = [sig_fid[i]/nom - 1.0 for i in PDF_IDS if i != PDF_NOMINAL_ID]
    if max(abs(x) for x in spread) < PDF_SELFCHECK_TOL:
        die("every PDF weight equals the nominal.  The `lhapdf` key was "
            "accepted and ignored -- see powheg/reweight/README.md, which "
            "documents exactly this silent failure.")

    scale_vals = [sig_fid[i] for i in SCALE_IDS]
    lo, hi = min(scale_vals), max(scale_vals)
    pdf_vals = [sig_fid[i] for i in PDF_IDS]

    ref, ref_src = hepmc_reference(current, energy, sel.name)
    closure = None
    if ref:
        closure = nom / ref
        print(f"  closure vs the showered result: {closure:.4f} "
              f"({100*(closure-1):+.2f}%; LHE {nom:.4f} pb vs HepMC "
              f"{ref:.4f} pb)")
        if abs(closure - 1.0) > CLOSURE_TOL:
            die(f"the nominal LHE weight does not reproduce the published "
                f"sigma_fid.\n  {os.path.relpath(ref_src, BASE)} says "
                f"{ref:.4f} pb, this says {nom:.4f} pb "
                f"({100*(closure-1):+.2f}%).\n  The frame reconstruction or "
                f"the selection has drifted; a band on top of that would look "
                f"perfectly reasonable.")
    else:
        print(f"  !! no closure check: {ref_src}")

    print(f"  {n_events} events, {n_fid} in the fiducial region")
    print(f"  sigma_fid (nominal)  = {nom:.5g} pb  "
          f"+- {err_fid[NOMINAL_ID]:.3g} (MC)")
    print(f"  sigma_tot (nominal)  = {sig_tot[NOMINAL_ID]:.5g} pb")
    print("  scale variation (xiR, xiF):")
    for i in SCALE_IDS:
        print(f"      {SCALE_LABEL[i]:>12}  {sig_fid[i]:9.5g} pb  "
              f"{100*(sig_fid[i]/nom - 1):+7.2f}%"
              f"  +- {100*err_shift[i]/nom:.3f}% (correlated)")
    print(f"      envelope      {100*(lo/nom-1):+.2f}% / "
          f"{100*(hi/nom-1):+.2f}%")
    print("  PDF sets:")
    for i in PDF_IDS:
        print(f"      {PDF_LABEL[i]:>22}  {sig_fid[i]:9.5g} pb  "
              f"{100*(sig_fid[i]/nom - 1):+7.2f}%"
              f"  +- {100*err_shift[i]/nom:.3f}% (correlated)")
    print(f"      spread        {100*(min(pdf_vals)/nom-1):+.2f}% / "
          f"{100*(max(pdf_vals)/nom-1):+.2f}%")
    print(f"  charm (parton-level tag): {n_charm} events, "
          f"sigma = {nom_charm:.5g} pb, fraction = {100*nom_frac:.2f}%")
    print("  charm fraction under the PDF sets:")
    for i in PDF_IDS:
        print(f"      {PDF_LABEL[i]:>22}  {100*frac[i]:7.3f}%  "
              f"{100*(frac[i]/nom_frac - 1):+7.2f}%")
    print("  charm fraction under the scale points:")
    for i in SCALE_IDS:
        print(f"      {SCALE_LABEL[i]:>12}  {100*frac[i]:7.3f}%  "
              f"{100*(frac[i]/nom_frac - 1):+7.2f}%")

    out = {
        "generator": "powheg_nu" if current == "nu" else "powheg_v2",
        "current": current,
        "energy_gev": energy,
        "beam_tag": beams.Beams(current, energy).tag,
        "selection": sel.as_dict(),
        "source_lhe": {
            "path": os.path.relpath(f, os.path.dirname(BASE)),
            "bytes": os.path.getsize(f),
            "mtime": os.path.getmtime(f)},
        "n_events": n_events,
        "n_fiducial": n_fid,
        "sigma_fid_pb": nom,
        "sigma_fid_err_pb": err_fid[NOMINAL_ID],
        "sigma_tot_pb": sig_tot[NOMINAL_ID],
        "closure_vs_hepmc": closure,
        "hepmc_sigma_fid_pb": ref,
        "n_charm": n_charm,
        "sigma_charm_pb": nom_charm,
        "charm_fraction": nom_frac,
        "charm_tag": "parton level: >=1 outgoing charm quark in the LHE",
        "charm_fraction_by_scale": [frac[i] for i in SCALE_IDS],
        "charm_fraction_by_pdf": [frac[i] for i in PDF_IDS],
        "sigma_charm_by_pdf_pb": [sig_charm[i] for i in PDF_IDS],
        "sigma_charm_by_scale_pb": [sig_charm[i] for i in SCALE_IDS],
        "scale": {
            "ids": SCALE_IDS,
            "labels": [SCALE_LABEL[i] for i in SCALE_IDS],
            "sigma_fid_pb": [sig_fid[i] for i in SCALE_IDS],
            "rel_lo": lo/nom - 1.0,
            "rel_hi": hi/nom - 1.0,
            "rel_err_correlated": [err_shift[i]/nom for i in SCALE_IDS]},
        "pdf": {
            "ids": PDF_IDS,
            "labels": [PDF_LABEL[i] for i in PDF_IDS],
            "sigma_fid_pb": [sig_fid[i] for i in PDF_IDS],
            "rel_lo": min(pdf_vals)/nom - 1.0,
            "rel_hi": max(pdf_vals)/nom - 1.0,
            "rel_err_correlated": [err_shift[i]/nom for i in PDF_IDS],
            "self_check": self_check},
    }
    outdir = f"{BASE}/{'results_nu' if current == 'nu' else 'results'}"
    parts = ["nlo_unc", out["generator"]]
    if sel.name != "inclusive":
        parts.append(sel.name)
    if energy != beams.ANCHOR_ENERGY:
        parts.append(out["beam_tag"])
    path = f"{outdir}/{'_'.join(parts)}.json"
    with open(path, "w") as fh:
        json.dump(out, fh)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
