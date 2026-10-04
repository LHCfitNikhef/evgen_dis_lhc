#!/usr/bin/env python3
"""The FASER dimuon cut flow, reproduced with our own machinery.

WHAT THIS IS FOR.  The opposite-sign dimuon signal of charm production in
nu_mu charged-current DIS at FASER, for the full Run 3 (300 fb^-1), as a cut
flow computed with this benchmark's machinery, with our GENIE and with
POWHEG-V2, under the FASER spectrometer fiducial cuts.  It was built to be
compared with a FASER-internal GENIE estimate; that reference is not public
and is read, only when present, from data/faser_internal/ (REF below).

HOW.  The same convolution as the rest of that tab (analysis/faser_rates.py,
faser_genie_rates.py, faser_powheg_rates.py): the vendored FASERnu flux, the
1.1 t of tungsten as 74 protons and 110 neutrons behind the aperture the flux
was counted through, and OUR cross-section -- but here every step after the
first needs the hadronic final state, so the cross-section is not a spline
value but a full event sample at each energy of a ladder:

    GENIE     gevgen, default tune G18_02a_00_000, the CC event-generator list
              (every CC process, as the reference's own GENIE would have run),
              Pythia8 decays through the benchmark overlay, at 12 energies from
              10 GeV to 6.8 TeV on p and on n, nu_mu and nubar_mu
              (tools/genie_dimuon_ladder.sh);
    POWHEG-V2 NLO, showered with Pythia 8 under the production steering card,
              at the 10 ladder energies of tools/powheg_v2_faser_ladder.sh,
              same four beams (tools/powheg_v2_dimuon_ladder.sh).

At each point the events are reduced to a compact table (extract mode: the
event kinematics, the primary muon, every opposite-sign muon), the HepMC is
deleted, and the rates mode convolves each step's cross-section
sigma_tot(E) x eff_step(E) with the flux -- sigma_tot log-log across the
ladder as everywhere on the tab, the efficiency linear in log E and held
constant beyond the ends.

>>> THE GENIE NUMBER IS ALSO COMPUTED WITHOUT ANY LADDER. <<<  gevgen can be
handed the flux histogram itself (-f) and then generates events distributed
as Phi(E) sigma(E): the flux-averaged efficiency of every step is then a
plain count on those events, with no interpolation anywhere.  Both routes are
run and the report shows their agreement; the ladder is what POWHEG-V2 can
do, so the ladder is the common method and the flux mode is its check.

THE CUT FLOW, AND WHAT EACH STEP MEANS HERE.
  1. CC interactions.  For GENIE, every CC interaction (QEL, RES, DIS, ...),
     which is what the reference counts; for POWHEG-V2, DIS with Q2 > 4 GeV2 and
     no y window -- the region the tab compares to arXiv:2402.13318 on,
     because below Q2 = 4 a perturbative generator is not to be trusted
     (user, 2026-09-03).  GENIE is ALSO shown restricted to Q2 > 4, so that
     the GENIE-versus-POWHEG comparison is like for like.
  2. With a charm hadron: the benchmark's final-state charm tag (a prompt
     charm hadron after hadronisation, the b -> c cascade excluded).
  3. Charm hadron -> mu: at least one muon of the sign OPPOSITE to the
     primary, any momentum.  Under the benchmark's ctau > 10 mm stability
     convention pi and K never decay, so an opposite-sign muon comes only
     from charm (or, above 1 TeV, at the per-cent level from bottom): the
     step measures the semileptonic branching ratio averaged over the
     generator's charm-species mix, and nothing else.
  4-6. Both muons with |p| > 20, 50, 100 GeV AND both inside theta < 25 mrad:
     the primary and at least one opposite-sign partner, in the FASER
     spectrometer acceptance.  The reference does not write the angle, but it is
     part of the spectrometer tier its rows live in (user, 2026-09-07), so
     these are the like-for-like rows.  The reference's "fiducial volume" is a
     cut on the vertex position in the emulsion, which no generator-level
     number can reproduce; its size is unknown, so the first row's ratio to
     the reference is reported as a measurement of that volume fraction and every
     later row is compared through its cumulative EFFICIENCY, which the
     volume cancels out of.
  4'-6'. The same momentum cuts WITHOUT the angular cut, kept beside them as
     the measure of what the acceptance costs.
  7. The benchmark's own dimuon tier (Q2 > 4, 0.2 < y < 0.9, both muons
     E > 100 GeV and theta < 25 mrad), so that this table joins onto the
     rest of the page.

WHAT IS NOT REPRODUCED, and is said in place: the reference's target is FASERnu
tungsten in an unspecified fiducial volume and, for the rate, an older flux;
nuclear effects are neglected here as on the rest of the tab; the reference's
row 1 may or may not include antineutrinos, so both nu_mu alone and
nu_mu + nubar_mu are given.

Usage:
  faser_dimuon_cutflow.py extract --beam 14|-14 --target p|n --energy E|flux \
                                  --generator genie|powheg_v2 in.hepmc out.npz
  faser_dimuon_cutflow.py rates          # -> results_nu/faser_dimuon_cutflow.json
  faser_dimuon_cutflow.py subset in.npz out.json   # the events to reweight
  faser_dimuon_cutflow.py pdf            # -> results_nu/faser_dimuon_cutflow_pdf.json
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import analyze_nu                                            # noqa: E402
import faser_rates as fr                                     # noqa: E402
import paths                                                 # noqa: E402
from analyze import (M_P, M_MU, dis_invariants, lab_energy,  # noqa: E402
                     lepton_theta)
from faser_genie_rates import NUCLEI_PER_CM2, Z_W, N_W, spline_total  # noqa: E402

# The reference's luminosity: the full Run 3.  Every number here is at THIS
# luminosity, not the 250 fb^-1 of the rest of the tab, and says so.
LUMI_FB = 300.0

# The FASER-internal reference cut flow (not distributed): compared with, and
# written into the result, only where the file exists.
_REF_FILE = os.path.join(BASE, "data", "faser_internal", "dimuon_cutflow_reference.json")
REF = None
if os.path.exists(_REF_FILE):
    with open(_REF_FILE) as _f:
        REF = json.load(_f)

# The cut flow.  (key, label, parent-step for the step efficiency)
#
# THE MOMENTUM ROWS CARRY THE SPECTROMETER ACCEPTANCE, theta < 25 mrad ON
# BOTH MUONS (user, 2026-09-07): the reference does not write it, but it is
# part of the FASER spectrometer tier the reference's rows are defined in, so
# the like-for-like rows are p20/p50/p100 WITH the angle, and the rows
# without it are kept beside them as the measure of what the acceptance
# costs.
#
# THE PARENT OF EVERY MOMENTUM ROW IS THE charm -> mu ROW, as on the reference:
# its "step efficiency" is 52.8% / 29.1% / 12.8% = 27.7 / 15.3 / 6.7 over
# 52.5, i.e. each momentum cut is quoted against the dimuon sample and not
# against the looser momentum cut above it.  Matching that convention is
# what makes the step columns comparable.
STEPS = [
    ("cc",       "CC interactions",                          None),
    ("charm",    "with a charm hadron",                      "cc"),
    ("charm_mu", "charm hadron -> mu (opposite sign)",       "charm"),
    ("p20",      "both mu |p| > 20 GeV, theta < 25 mrad",    "charm_mu"),
    ("p50",      "both mu |p| > 50 GeV, theta < 25 mrad",    "charm_mu"),
    ("p100",     "both mu |p| > 100 GeV, theta < 25 mrad",   "charm_mu"),
    ("p20_noth", "both mu |p| > 20 GeV, no angular cut",     "charm_mu"),
    ("p50_noth", "both mu |p| > 50 GeV, no angular cut",     "charm_mu"),
    ("p100_noth", "both mu |p| > 100 GeV, no angular cut",   "charm_mu"),
    ("bench",    "benchmark dimuon tier (Q2 > 4, 0.2 < y < 0.9, "
                 "both E > 100 GeV and theta < 25 mrad)",     "cc"),
]
THETA_MAX = 0.025
Q2_MIN = 4.0

# The three columns.  "base" says what row 1 is.
COLUMNS = [
    ("genie_total", "GENIE, all CC",            "genie",     "all"),
    ("genie_q2",    "GENIE, Q2 > 4 GeV2",       "genie",     "q2"),
    ("powheg_q2",   "POWHEG-V2, Q2 > 4 GeV2",   "powheg_v2", "q2"),
]

GENIE_LADDER_ES = [10, 20, 30, 60, 100, 200, 400, 700, 1000, 2000, 4000, 6800]
POWHEG_LADDER_ES = [30, 60, 100, 200, 400, 700, 1000, 2000, 4000, 6800]


# ------------------------------------------------------------------ extract
def extract(hepmc, out, beam, target, energy, generator):
    """One HepMC file -> one compact per-event table.

    Everything a later cut on the muons can ask for is kept, so the HepMC
    can go: E_nu (lab), Q2, y, W2, weight, the charm and bottom tags, the
    primary muon's lab energy and angle, up to two opposite-sign muons (lab
    energy and angle, hardest first) and how many there were, and how many
    EXTRA same-sign muons the event carried (a diagnostic: with pi and K
    stable those come only from b -> c -> mu chains).
    """
    beam = int(beam)
    lep_out = 13 if beam > 0 else -13
    # the parser reads these module globals at run time; the antineutrino
    # beam and the neutron target are the two things the production
    # analysis never had to know about
    analyze_nu.NU_BEAM_PID = beam
    cols = {k: [] for k in ("idx", "enu", "w", "q2", "y", "w2", "emu", "thmu",
                            "charm", "bottom", "nos", "eos1", "thos1",
                            "eos2", "thos2", "nss")}
    n_parsed = n_no_primary = 0
    for w, k, P, parts, _dmes, hard in analyze_nu.parse_hepmc3(hepmc):
        n_parsed += 1
        w = 1.0 if w is None else w
        best, best_e = None, -1.0
        os_mu, n_ss = [], 0
        for pid, p in parts:
            if pid == lep_out:
                e = lab_energy(p, P)
                if e > best_e:
                    if best is not None:
                        n_ss += 1
                    best, best_e = p, e
                else:
                    n_ss += 1
            elif pid == -lep_out:
                e = lab_energy(p, P)
                os_mu.append((e, lepton_theta(e, p)))
        if best is None:
            n_no_primary += 1
            continue
        q2, y, _x, kP = dis_invariants(k, P, best)
        os_mu.sort(reverse=True)
        heavy = hard.get("fs_heavy", frozenset()) if isinstance(hard, dict) \
            else frozenset()
        # the Les Houches event number where the sample has one (POWHEG), so
        # a PDF-reweighted weight column can be joined onto this table later
        # without re-parsing; -1 for GENIE, whose events carry none
        li = hard.get("lhe_index") if isinstance(hard, dict) else None
        cols["idx"].append(-1 if li is None else li)
        cols["enu"].append(kP / M_P)
        cols["w"].append(w)
        cols["q2"].append(q2)
        cols["y"].append(y)
        cols["w2"].append(M_P * M_P + 2.0 * y * kP - q2)
        cols["emu"].append(best_e)
        cols["thmu"].append(lepton_theta(best_e, best))
        cols["charm"].append(4 in heavy)
        cols["bottom"].append(5 in heavy)
        cols["nos"].append(len(os_mu))
        cols["eos1"].append(os_mu[0][0] if os_mu else -1.0)
        cols["thos1"].append(os_mu[0][1] if os_mu else -1.0)
        cols["eos2"].append(os_mu[1][0] if len(os_mu) > 1 else -1.0)
        cols["thos2"].append(os_mu[1][1] if len(os_mu) > 1 else -1.0)
        cols["nss"].append(n_ss)
    meta = {"beam": beam, "target": target, "energy": energy,
            "generator": generator, "source": os.path.abspath(hepmc),
            "n_parsed": n_parsed, "n_no_primary": n_no_primary,
            "n_kept": len(cols["w"])}
    arrays = {k: np.asarray(v) for k, v in cols.items()}
    np.savez_compressed(out, meta=json.dumps(meta), **arrays)
    ch = float(arrays["charm"].mean()) if len(arrays["w"]) else 0.0
    mu = float(((arrays["nos"] > 0) & arrays["charm"]).mean()) if len(arrays["w"]) else 0.0
    print(f"{n_parsed} events parsed, {n_no_primary} without a primary muon, "
          f"charm {100*ch:.2f}%, charm with an opposite-sign muon {100*mu:.3f}%")


# ------------------------------------------------- the reweighting subset
# WHICH EVENTS NEED A PDF WEIGHT.  Every row from "charm hadron -> mu" down
# is a subset of the charm_mu events, so reweighting those (about 1% of a
# point) gives every dimuon row EXACTLY (powheg/reweight/run_reweight_selected.sh
# says why: an event outside the row contributes nothing to the row's member
# spread).  The CC row is every event and the charm row every charm event;
# each gets a fixed random sample of RANDOM_N events, which gives its band to
# the per-mille since a member weight is a smooth function of (x, Q2).  The
# three parts go into one index list, listed separately so the analysis
# normalises each to its own denominator.  Measured before this was done:
# 77 s per pass of six weights on 21k events, i.e. 12 hours for the ladder;
# 7k events is a third of that.
RANDOM_N = 1500


def subset(npz, out, seed=1):
    d, meta = load_tables([npz])
    idx = d["idx"]
    if (idx < 0).any():
        raise SystemExit(f"{npz}: events without an lhe_index cannot be joined")
    rng = np.random.default_rng(seed)
    charm_mu = sorted(int(i) for i in idx[d["charm"] & (d["nos"] >= 1)])
    where_charm = np.flatnonzero(d["charm"])
    pick = rng.choice(where_charm, size=min(RANDOM_N, len(where_charm)), replace=False)
    rand_charm = sorted(int(i) for i in idx[pick])
    pick = rng.choice(len(idx), size=min(RANDOM_N, len(idx)), replace=False)
    rand_cc = sorted(int(i) for i in idx[pick])
    allidx = sorted(set(charm_mu) | set(rand_charm) | set(rand_cc))
    with open(out, "w") as f:
        json.dump({"source": meta["source"], "n_events": int(len(idx)),
                   "n_charm": int(d["charm"].sum()),
                   "lhe_index": allidx, "charm_mu": charm_mu,
                   "random_charm": rand_charm, "random_cc": rand_cc,
                   "random_seed": seed}, f)
    print(f"{len(allidx)} of {len(idx)} events to reweight "
          f"({len(charm_mu)} charm->mu, {len(rand_charm)} random charm, "
          f"{len(rand_cc)} random for the CC row)")


# ---------------------------------------------------------------- the masks
def step_masks(d, base):
    """Cumulative boolean masks for every step, on one table.

    base: "all" (every event) or "q2" (Q2 > 4 GeV2).  Momenta from the lab
    energies; a muon's |p| is sqrt(E^2 - m^2), which at 20 GeV is E to 1e-5.
    """
    pmu = np.sqrt(np.maximum(d["emu"] ** 2 - M_MU ** 2, 0.0))
    p1 = np.sqrt(np.maximum(d["eos1"] ** 2 - M_MU ** 2, 0.0))
    p2 = np.sqrt(np.maximum(d["eos2"] ** 2 - M_MU ** 2, 0.0))
    has1, has2 = d["nos"] >= 1, d["nos"] >= 2
    th_mu = d["thmu"] < THETA_MAX
    th1 = has1 & (d["thos1"] < THETA_MAX)
    th2 = has2 & (d["thos2"] < THETA_MAX)
    m = {}
    m["cc"] = np.ones(len(d["w"]), bool) if base == "all" else d["q2"] > Q2_MIN
    m["charm"] = m["cc"] & d["charm"]
    m["charm_mu"] = m["charm"] & has1
    for pcut in (20.0, 50.0, 100.0):
        partner = (has1 & (p1 > pcut)) | (has2 & (p2 > pcut))
        m[f"p{pcut:.0f}_noth"] = m["charm_mu"] & (pmu > pcut) & partner
        partner_th = (th1 & (p1 > pcut)) | (th2 & (p2 > pcut))
        m[f"p{pcut:.0f}"] = m["charm_mu"] & (pmu > pcut) & th_mu & partner_th
    partner_b = (th1 & (d["eos1"] > 100.0)) | (th2 & (d["eos2"] > 100.0))
    m["bench"] = ((d["q2"] > Q2_MIN) & (d["y"] > 0.2) & (d["y"] < 0.9)
                  & (d["emu"] > 100.0) & th_mu & partner_b)
    return m


def load_tables(files):
    """Concatenate the tables of one ladder point (several job files)."""
    parts, meta = [], None
    for f in sorted(files):
        z = np.load(f)
        meta = json.loads(str(z["meta"]))
        parts.append({k: z[k] for k in z.files if k != "meta"})
    if not parts:
        return None, None
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    return d, meta


def fractions(d, base):
    """Weighted fraction of ALL events in each step, its MC error, and the
    raw selected count.  The error is binomial on the effective number of
    events (sum w)^2 / sum w^2, which is exact for GENIE's unit weights."""
    w = d["w"]
    sw, sw2 = float(w.sum()), float((w * w).sum())
    neff = sw * sw / sw2 if sw2 > 0 else 0.0
    out = {}
    for key, mask in step_masks(d, base).items():
        f = float(w[mask].sum()) / sw if sw else 0.0
        err = math.sqrt(max(f * (1.0 - f), 0.0) / neff) if neff else 0.0
        out[key] = (f, err, int(mask.sum()))
    return out


# --------------------------------------------------------------- the ladders
def genie_sigma_tot(pid, tgt, es):
    """GENIE default-tune total CC per nucleon [pb] at the ladder energies,
    from the same splines the events were generated with."""
    e_knots, sig, _names = spline_total(str(pid), tgt)
    ok = sig > 0
    return np.exp(np.interp(np.log(es), np.log(e_knots[ok]), np.log(sig[ok])))


def powheg_total_pb(rundir):
    for line in open(os.path.join(rundir, "pwg-stat.dat"), errors="replace"):
        if "total (btilde+remnants) cross section in pb" in line:
            v = line.split("pb")[1].split()[0]
            return float(v.replace("D", "E").replace("d", "e"))
    raise SystemExit(f"no total in {rundir}/pwg-stat.dat")


def hadroniser_failures(files):
    """How many charm events GENIE's AGCharm hadroniser gave up on ("Got an
    empty particle list. Hadronizer failed!" -- the event is redrawn, so the
    delivered charm fraction is short by that many; per mille at 10 GeV and
    zero above).  Counted from the reduced logs beside the tables; the
    digest keeps at most 300 complaint lines, so a saturated count is a
    lower bound and is flagged as such."""
    n, sat = 0, False
    for f in files:
        dg = os.path.join(os.path.dirname(f), "genie_log_digest.txt")
        if not os.path.exists(dg):
            continue
        lines = [ln for ln in open(dg, errors="replace")
                 if "Hadronizer failed" in ln]
        n += len(lines)
        if sum(1 for ln in open(dg, errors="replace")
               if ln.split(":", 1)[0].isdigit()) >= 300:
            sat = True
    return n, sat


DIAG_KEYS = ("frac_extra_same_sign_mu", "frac_bottom")


def diagnostics(d):
    """Per-point diagnostics: the fraction of events with an EXTRA same-sign
    muon (a light-meson mu+ mu- decay makes one of each, so this bounds the
    fake opposite-sign tags) and the fraction with a bottom hadron."""
    return {"frac_extra_same_sign_mu": float((d["nss"] > 0).mean()),
            "frac_bottom": float(d["bottom"].mean())}


def genie_ladder(cur, tgt):
    """[(E, sigma_tot_pb, fractions{all}, fractions{q2}, n_events, n_files,
    extras)] sorted."""
    pid = 14 if cur == "nu" else -14
    pts = []
    for e in GENIE_LADDER_ES:
        files = glob.glob(f"{BASE}/genie/dmuladder_job_{cur}_{tgt}_E{e}_[0-9]*/dimuon_events.npz")
        d, meta = load_tables(files)
        if d is None:
            continue
        nh, sat = hadroniser_failures(files)
        extra = dict(diagnostics(d), n_hadroniser_failed=nh,
                     n_hadroniser_failed_is_lower_bound=sat)
        pts.append((float(e), d, meta, len(files), extra))
    if not pts:
        return []
    es = np.array([p[0] for p in pts])
    sig = genie_sigma_tot(pid, tgt, es)
    return [(e, float(s), fractions(d, "all"), fractions(d, "q2"), len(d["w"]), nf, extra)
            for (e, d, meta, nf, extra), s in zip(pts, sig)]


def powheg_ladder(cur, tgt):
    pts = []
    for e in POWHEG_LADDER_ES:
        r = f"{paths.POWHEG_V2}/ladder-dimuon/{cur}_{tgt}_E{e}"
        f = os.path.join(r, "dimuon_events.npz")
        if not os.path.exists(f):
            continue
        d, meta = load_tables([f])
        pts.append((float(e), powheg_total_pb(r), None, fractions(d, "q2"),
                    len(d["w"]), 1, diagnostics(d)))
    return pts


def genie_flux_mode(cur, tgt):
    files = glob.glob(f"{BASE}/genie/dmuflux_job_{cur}_{tgt}_[0-9]*/dimuon_events.npz")
    d, meta = load_tables(files)
    if d is None:
        return None
    nh, sat = hadroniser_failures(files)
    return dict({"all": fractions(d, "all"), "q2": fractions(d, "q2"),
                 "n_events": len(d["w"]), "n_files": len(files),
                 "mean_enu": float(d["enu"].mean()),
                 "n_hadroniser_failed": nh,
                 "n_hadroniser_failed_is_lower_bound": sat}, **diagnostics(d))


# ------------------------------------------------------------ the convolution
def convolve(pid, es, sig_tot, eff, lumi_fb=LUMI_FB, per_bin=False):
    """Events at lumi_fb on the tungsten target: the flux (150 fb^-1 files)
    times sigma_tot(E) x eff(E) per NUCLEON, times the column density.
    sigma_tot log-log across the ladder and ZERO below its first point, as
    faser_powheg_rates.rate does; eff linear in log E, held at the ends.
    per_bin=True returns the contribution of every flux bin instead."""
    e, phi = fr.flux(str(pid))
    phi = phi * lumi_fb / 150.0
    s = np.exp(np.interp(np.log(e), np.log(es), np.log(sig_tot),
                         left=-np.inf, right=np.log(sig_tot[-1])))
    f = np.interp(np.log(e), np.log(es), eff)
    v = phi * s * f * 1e-36 * NUCLEI_PER_CM2
    return v if per_bin else float(np.sum(v))


def column_rates(lad_p, lad_n, pid, which):
    """Every step's events for one beam from a p and an n ladder.

    which: index 2 (fractions on all events) or 3 (fractions with Q2 > 4).
    Tungsten = Z sigma_p + N sigma_n per nucleus.  The error propagates each
    ladder point's binomial error through the convolution one at a time and
    sums in quadrature; points are independent samples.
    """
    out = {}
    es = np.array([p[0] for p in lad_p])
    assert np.allclose(es, [p[0] for p in lad_n]), "p and n ladders differ"
    sp = np.array([p[1] for p in lad_p])
    sn = np.array([p[1] for p in lad_n])
    for key, _lab, _par in STEPS:
        fp = np.array([p[which][key][0] for p in lad_p])
        fn = np.array([p[which][key][0] for p in lad_n])
        dfp = np.array([p[which][key][1] for p in lad_p])
        dfn = np.array([p[which][key][1] for p in lad_n])
        n_p = Z_W * convolve(pid, es, sp, fp)
        n_n = N_W * convolve(pid, es, sn, fn)
        var = 0.0
        for i in range(len(es)):
            fpi, fni = fp.copy(), fn.copy()
            fpi[i] += dfp[i]
            fni[i] += dfn[i]
            var += (Z_W * convolve(pid, es, sp, fpi) - n_p) ** 2
            var += (N_W * convolve(pid, es, sn, fni) - n_n) ** 2
        # where in the flux the events come from, bin by bin (for the figure)
        spec = (Z_W * convolve(pid, es, sp, fp, per_bin=True)
                + N_W * convolve(pid, es, sn, fn, per_bin=True))
        out[key] = {"events": n_p + n_n, "err": math.sqrt(var),
                    "events_p": n_p, "events_n": n_n,
                    "n_selected": int(sum(p[which][key][2] for p in lad_p)
                                      + sum(p[which][key][2] for p in lad_n)),
                    "spectrum": [float(x) for x in spec]}
    out["spectrum_E_GeV"] = [float(x) for x in fr.flux(str(pid))[0]]
    return out


def flux_mode_rates(fm_p, fm_n, pid, base):
    """GENIE flux mode: the flux-averaged fraction is a plain count, and the
    total is the spline convolution of faser_genie_rates on each nucleon."""
    es_p, sig_p, _ = spline_total(str(pid), "p")
    es_n, sig_n, _ = spline_total(str(pid), "n")
    e, phi = fr.flux(str(pid))
    phi = phi * LUMI_FB / 150.0

    def tot(es, sig):
        ok = sig > 0
        s = np.exp(np.interp(np.log(e), np.log(es[ok]), np.log(sig[ok])))
        s[e < es[ok][0]] = 0.0
        return float(np.sum(phi * s * 1e-36 * NUCLEI_PER_CM2))
    tp, tn = Z_W * tot(es_p, sig_p), N_W * tot(es_n, sig_n)
    out = {"total_p": tp, "total_n": tn}
    for key, _lab, _par in STEPS:
        fp, dfp, np_ = fm_p[base][key]
        fn, dfn, nn_ = fm_n[base][key]
        out[key] = {"events": tp * fp + tn * fn,
                    "err": math.hypot(tp * dfp, tn * dfn),
                    "n_selected": np_ + nn_}
    return out


def add_efficiencies(rows):
    """Step and cumulative efficiencies in per cent, in place."""
    for key, _lab, par in STEPS:
        r = rows[key]
        n0 = rows["cc"]["events"]
        r["cum_eff_pct"] = 100.0 * r["events"] / n0 if n0 else 0.0
        if par is None:
            r["step_eff_pct"] = 100.0
        else:
            np_ = rows[par]["events"]
            r["step_eff_pct"] = 100.0 * r["events"] / np_ if np_ else 0.0
    return rows


def combine(a, b):
    """nu_mu + nubar_mu."""
    out = {}
    for key, _lab, _par in STEPS:
        out[key] = {"events": a[key]["events"] + b[key]["events"],
                    "err": math.hypot(a[key]["err"], b[key]["err"]),
                    "n_selected": a[key]["n_selected"] + b[key]["n_selected"]}
        if "spectrum" in a[key] and "spectrum" in b[key]:
            out[key]["spectrum"] = [x + y for x, y in
                                    zip(a[key]["spectrum"], b[key]["spectrum"])]
    if "spectrum_E_GeV" in a:
        out["spectrum_E_GeV"] = a["spectrum_E_GeV"]
    return add_efficiencies(out)


# --------------------------------------------------------------------- rates
def rates():
    out = {"reference": REF, "lumi_fb": LUMI_FB,
           "target": {"mass_g": fr.FASERNU_MASS_G, "Z": Z_W, "N": N_W,
                      "nuclei_per_cm2": NUCLEI_PER_CM2,
                      "note": "1.1 t of tungsten behind the aperture the flux "
                              "was counted through (faser_rates.py)"},
           "steps": [{"key": k, "label": l, "parent": p} for k, l, p in STEPS],
           "theta_max_rad": THETA_MAX, "q2_min": Q2_MIN,
           "columns": {}, "ladders": {}, "flux_mode": {}}
    lad = {}
    for gen, fn in (("genie", genie_ladder), ("powheg_v2", powheg_ladder)):
        for cur in ("nu", "nubar"):
            for tgt in ("p", "n"):
                pts = fn(cur, tgt)
                lad[(gen, cur, tgt)] = pts
                out["ladders"][f"{gen}_{cur}_{tgt}"] = [
                    dict({"E_GeV": e, "sigma_tot_pb": s, "n_events": n, "n_files": nf,
                          "frac_all": ({k: list(v) for k, v in fa.items()} if fa else None),
                          "frac_q2": {k: list(v) for k, v in fq.items()}}, **extra)
                    for e, s, fa, fq, n, nf, extra in pts]
                if pts:
                    print(f"{gen} {cur} on {tgt}: {len(pts)} ladder points, "
                          f"{sum(p[4] for p in pts)} events")
    # the three columns, per beam and summed
    for col, label, gen, base in COLUMNS:
        which = 2 if base == "all" else 3
        rec = {"label": label, "generator": gen, "base": base}
        ok = True
        for cur, pid in (("nu", 14), ("nubar", -14)):
            lp, ln = lad[(gen, cur, "p")], lad[(gen, cur, "n")]
            if len(lp) < 3 or len(ln) < 3 or len(lp) != len(ln):
                print(f"  {col}: ladder for {cur} incomplete ({len(lp)} p, {len(ln)} n)")
                ok = False
                break
            rec[str(pid)] = add_efficiencies(column_rates(lp, ln, pid, which))
        if not ok:
            continue
        rec["sum"] = combine(rec["14"], rec["-14"])
        out["columns"][col] = rec
    # GENIE flux mode, both bases
    fm = {(cur, tgt): genie_flux_mode(cur, tgt)
          for cur in ("nu", "nubar") for tgt in ("p", "n")}
    if all(v is not None for v in fm.values()):
        for base, col in (("all", "genie_total"), ("q2", "genie_q2")):
            rec = {"column": col, "base": base,
                   "n_events": {f"{c}_{t}": v["n_events"] for (c, t), v in fm.items()},
                   "mean_enu": {f"{c}_{t}": v["mean_enu"] for (c, t), v in fm.items()}}
            for cur, pid in (("nu", 14), ("nubar", -14)):
                rec[str(pid)] = add_efficiencies(
                    flux_mode_rates(fm[(cur, "p")], fm[(cur, "n")], pid, base))
            rec["sum"] = combine(rec["14"], rec["-14"])
            # closure of the ladder against the flux mode, per step
            if col in out["columns"]:
                rec["ladder_over_flux"] = {
                    k: (out["columns"][col]["sum"][k]["events"] / rec["sum"][k]["events"]
                        if rec["sum"][k]["events"] else None)
                    for k, _l, _p in STEPS}
            out["flux_mode"][col] = rec
    else:
        print("  GENIE flux mode: not all four (beam, nucleon) runs present")
    # against the reference
    cmp = {}
    for col in (out["columns"] if REF else []):
        for who in ("14", "sum"):
            rows = out["columns"][col][who]
            c = {}
            for key, (nev, step, cum) in REF["rows"].items():
                r = rows[key]
                c[key] = {"slide_events": nev, "ours_events": r["events"],
                          "ours_over_slide": r["events"] / nev,
                          "slide_cum_pct": cum, "ours_cum_pct": r["cum_eff_pct"],
                          "cum_ratio": r["cum_eff_pct"] / cum,
                          "slide_step_pct": step, "ours_step_pct": r["step_eff_pct"],
                          "step_ratio": r["step_eff_pct"] / step}
            cmp[f"{col}_{who}"] = c
    out["vs_slide"] = cmp
    # the per-flavour split of row 1 from the flux, for the reader
    p = f"{BASE}/results_nu/faser_dimuon_cutflow.json"
    # >>> NEVER OVERWRITE A COMPLETE RESULT WITH AN INCOMPLETE ONE. <<<  The
    # ladders are event samples outside the repository; on a machine (or a
    # fresh clone) without them every column above is skipped, and writing
    # would silently replace the tracked result with an empty one -- found
    # 2026-09-30 running tools/make_all_plots.sh in a clean clone.
    missing = [c for c, *_x in COLUMNS if c not in out["columns"]]
    if missing and os.path.exists(p):
        sys.exit(f"columns {missing} could not be built (ladder samples not on "
                 f"this machine); keeping the existing {os.path.relpath(p, BASE)}")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print_table(out)
    print(f"\nwritten {p}")
    return out


def print_table(out):
    cols = [c for c, *_ in COLUMNS if c in out["columns"]]
    print(f"\nFASER dimuon cut flow at {LUMI_FB:.0f} fb^-1, 1.1 t of tungsten "
          f"(events; cumulative efficiency %)")
    hdr = f"{'step':42s} {'reference':>14s}"
    for c in cols:
        hdr += f" {c + ' nu':>18s} {c + ' nu+nubar':>20s}"
    print(hdr)
    for key, lab, _par in STEPS:
        line = f"{lab[:42]:42s} "
        if REF and key in REF["rows"]:
            nev, _s, cum = REF["rows"][key]
            line += f"{nev:7.1f} {cum:6.3f}"
        else:
            line += f"{'':14s}"
        for c in cols:
            for who in ("14", "sum"):
                r = out["columns"][c][who][key]
                line += f" {r['events']:9.1f} {r['cum_eff_pct']:8.3f}"
        print(line)
    for col, rec in out.get("flux_mode", {}).items():
        print(f"\nGENIE flux mode, {col}: nu / nu+nubar events and ladder/flux ratio")
        for key, lab, _par in STEPS:
            lf = rec.get("ladder_over_flux", {}).get(key)
            print(f"  {lab[:42]:42s} {rec['14'][key]['events']:9.1f} "
                  f"{rec['sum'][key]['events']:9.1f}   "
                  + (f"{lf:.3f}" if lf else "--"))


# ------------------------------------------------------- the PDF dependence
# POWHEG-V2's cut flow under six PDF sets with their bands (user, 2026-09-07),
# from the member weights tools/powheg_v2_dimuon_ladder_pdf.sh appends to
# the reweighting subset of every ladder point.  For every weight id m and
# every row:
#
#     sigma_row(m) = sigma_tot,nominal x S_row(m) / sum_all w_nominal
#     S_row(m)     = sum_{subset & row} w_tab x (w_m / w_nominal) x scale_row
#
# where w_tab is the event weight in the table, w_m / w_nominal the member's
# factor from the reweighted LHE (a RATIO, so the LHE and HepMC weight
# conventions cancel), and scale_row = sum_{all & row} w_tab / sum_{subset &
# row} w_tab -- exactly 1 for the dimuon rows, whose events are all in the
# subset, and the sampling correction for the CC and charm rows, whose
# subsets are random samples.  The per-member cross-sections then go through
# the same tungsten convolution as the nominal ones, and each set's band is
# LHAPDF's own uncertainty() over its members at 68% CL, as every PDF band
# in this benchmark is.
NOMINAL_ID = "3000"        # NNPDF4.0 member 0, the set the ladder was made with


def load_point_weights(rundir):
    """(table, ids, per-row member factors w_m / w_nominal, subset groups)."""
    d, _meta = load_tables([os.path.join(rundir, "dimuon_events.npz")])
    z = np.load(os.path.join(rundir, "pdf", "weights.npz"))
    ids, w = list(z["ids"]), z["weights"]
    with open(os.path.join(rundir, "pdf", "pwgevents-subset.index.json")) as f:
        lhe_rows = json.load(f)["lhe_index"]
    with open(os.path.join(rundir, "pdf", "subset.json")) as f:
        groups = json.load(f)
    if len(lhe_rows) != w.shape[0]:
        raise SystemExit(f"{rundir}: {len(lhe_rows)} indices for {w.shape[0]} weight rows")
    pos = {int(i): k for k, i in enumerate(d["idx"])}
    rows = np.array([pos[i] for i in lhe_rows])          # table row per weight row
    nom = w[:, ids.index(NOMINAL_ID)]
    fac = w / nom[:, None]
    return d, ids, rows, fac, groups


def pdf_point(rundir):
    """S_row(m) / sum_all w_tab for every row and every weight id, one point."""
    d, ids, rows, fac, groups = load_point_weights(rundir)
    masks = step_masks(d, "q2")
    w_tab = d["w"]
    tot = float(w_tab.sum())
    in_sub = {g: np.zeros(len(w_tab), bool) for g in ("charm_mu", "random_charm", "random_cc")}
    pos = {int(i): k for k, i in enumerate(d["idx"])}
    for g in in_sub:
        for i in groups[g]:
            in_sub[g][pos[i]] = True
    sub_row = np.zeros(len(w_tab), int) - 1
    sub_row[rows] = np.arange(len(rows))
    out = {}
    for key, _lab, _par in STEPS:
        m = masks[key]
        if key == "cc":
            g = in_sub["random_cc"]
        elif key == "charm":
            g = in_sub["random_charm"]
        else:
            g = in_sub["charm_mu"]
        sel = m & g
        denom = float(w_tab[sel].sum())
        if denom == 0.0:
            out[key] = np.zeros(fac.shape[1])
            continue
        scale = float(w_tab[m].sum()) / denom
        r = sub_row[sel]
        assert (r >= 0).all(), f"{rundir}: a {key} event is outside the reweighted subset"
        out[key] = (w_tab[sel][:, None] * fac[r]).sum(axis=0) * scale / tot
    return ids, out


def pdf_rates():
    import lhapdf
    lhapdf.setVerbosity(0)
    with open(f"{BASE}/powheg/reweight/rwl_dimuon_p.json") as f:
        rwl = json.load(f)
    lad = {}
    ids_ref = None
    for cur in ("nu", "nubar"):
        for tgt in ("p", "n"):
            pts = []
            for e in POWHEG_LADDER_ES:
                r = f"{paths.POWHEG_V2}/ladder-dimuon/{cur}_{tgt}_E{e}"
                if not os.path.exists(os.path.join(r, "pdf", "weights.npz")):
                    continue
                ids, fr = pdf_point(r)
                if ids_ref is None:
                    ids_ref = ids
                assert ids == ids_ref, f"{r}: weight ids differ"
                pts.append((float(e), powheg_total_pb(r), fr))
            lad[(cur, tgt)] = pts
            print(f"pdf: {cur} on {tgt}: {len(pts)} reweighted points")
    if any(len(v) < 3 for v in lad.values()):
        raise SystemExit("pdf: fewer than 3 reweighted points on a ladder -- "
                         "run tools/powheg_v2_dimuon_ladder_pdf.sh")
    nw = len(ids_ref)
    # events per row per weight id, per beam and summed
    ev = {}
    for cur, pid in (("nu", 14), ("nubar", -14)):
        lp, ln = lad[(cur, "p")], lad[(cur, "n")]
        es = np.array([p[0] for p in lp])
        assert np.allclose(es, [p[0] for p in ln])
        sp, sn = np.array([p[1] for p in lp]), np.array([p[1] for p in ln])
        ev[str(pid)] = {}
        for key, _lab, _par in STEPS:
            vals = np.zeros(nw)
            for j in range(nw):
                fp = np.array([p[2][key][j] for p in lp])
                fn = np.array([p[2][key][j] for p in ln])
                vals[j] = Z_W * convolve(pid, es, sp, fp) + N_W * convolve(pid, es, sn, fn)
            ev[str(pid)][key] = vals
    ev["sum"] = {k: ev["14"][k] + ev["-14"][k] for k, _l, _p in STEPS}
    # per set: central and 68% CL band from LHAPDF's own prescription
    out = {"lumi_fb": LUMI_FB, "nominal_set": "NNPDF40_nnlo_as_01180",
           "confidence_level_percent": 100.0 * math.erf(1.0 / math.sqrt(2.0)),
           "steps": [k for k, _l, _p in STEPS], "sets": {}, "scale": {},
           "n_points": {f"{c}_{t}": len(v) for (c, t), v in lad.items()}}
    cl = out["confidence_level_percent"]
    nomv = {who: {k: float(ev[who][k][ids_ref.index(NOMINAL_ID)]) for k, _l, _p in STEPS}
            for who in ("14", "-14", "sum")}
    for name, rec in rwl["sets"].items():
        cols = [ids_ref.index(i) for i in rec["weight_ids"]]
        pset = lhapdf.getPDFSet(name)
        srec = {"error_type": rec["error_type"], "members": len(cols), "rows": {}}
        for who in ("14", "-14", "sum"):
            srec["rows"][who] = {}
            for key, _l, _p in STEPS:
                vals = [float(ev[who][key][c]) for c in cols]
                if len(vals) >= 2:
                    try:
                        u = pset.uncertainty(vals, cl)
                        c0, ep, em = u.central, u.errplus, u.errminus
                    except Exception as exc:                  # noqa: BLE001
                        print(f"  !! {name}: uncertainty() failed ({exc})")
                        c0 = vals[0]
                        ep = em = float(np.std(vals[1:]))
                else:
                    c0, ep, em = vals[0], 0.0, 0.0
                srec["rows"][who][key] = {
                    "events": c0, "err_plus": ep, "err_minus": em,
                    "rel_to_nominal": c0 / nomv[who][key] if nomv[who][key] else None}
        out["sets"][name] = srec
    # the 7-point scale envelope on the nominal set
    sc = [ids_ref.index(i) for i in rwl["scales"]]
    for who in ("14", "-14", "sum"):
        out["scale"][who] = {}
        for key, _l, _p in STEPS:
            vals = [float(ev[who][key][c]) for c in sc]
            out["scale"][who][key] = {"central": vals[0], "lo": min(vals), "hi": max(vals)}
    # closure: the nominal member must reproduce the unweighted cut flow
    p_main = f"{BASE}/results_nu/faser_dimuon_cutflow.json"
    if os.path.exists(p_main):
        with open(p_main) as f:
            main_ = json.load(f)["columns"]["powheg_q2"]
        out["closure_vs_cutflow"] = {
            k: nomv["sum"][k] / main_["sum"][k]["events"] if main_["sum"][k]["events"] else None
            for k, _l, _p in STEPS}
    p = f"{BASE}/results_nu/faser_dimuon_cutflow_pdf.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nPOWHEG-V2 cut flow under the PDF sets, nu_mu + nubar_mu, {LUMI_FB:.0f} fb^-1 "
          f"(events, +/- 68% CL PDF band, ratio to NNPDF4.0)")
    for key, lab, _p in STEPS:
        line = f"{lab[:38]:38s}"
        for name, srec in out["sets"].items():
            r = srec["rows"]["sum"][key]
            line += f"  {name[:8]:8s} {r['events']:7.1f} +{r['err_plus']:4.1f} -{r['err_minus']:4.1f} ({r['rel_to_nominal']:.3f})"
        print(line)
    if "closure_vs_cutflow" in out:
        print("closure of the nominal member against the unweighted cut flow:",
              " ".join(f"{k} {v:.4f}" for k, v in out["closure_vs_cutflow"].items() if v))
    print(f"written {p}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    ex = sub.add_parser("extract")
    ex.add_argument("--beam", type=int, required=True, help="14 or -14")
    ex.add_argument("--target", required=True, help="p or n")
    ex.add_argument("--energy", required=True, help="beam energy in GeV, or 'flux'")
    ex.add_argument("--generator", required=True, help="genie or powheg_v2")
    ex.add_argument("hepmc")
    ex.add_argument("out")
    sub.add_parser("rates")
    sb = sub.add_parser("subset")
    sb.add_argument("npz")
    sb.add_argument("out")
    sub.add_parser("pdf")
    a = ap.parse_args()
    if a.cmd == "extract":
        extract(a.hepmc, a.out, a.beam, a.target, a.energy, a.generator)
    elif a.cmd == "subset":
        subset(a.npz, a.out)
    elif a.cmd == "pdf":
        pdf_rates()
    else:
        rates()


if __name__ == "__main__":
    main()
