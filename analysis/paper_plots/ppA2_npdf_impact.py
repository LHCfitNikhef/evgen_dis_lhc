#!/usr/bin/env python3
"""Appendix plot 2: the impact of nuclear PDFs on the inclusive DIS
distributions, in the final region.

The successor of the earlier production's nuclear-PDF study (cmp_npdf_dis_{mu,nu}.png, deleted
with the earlier results in ea61871).

  * THE CALCULATION: YADISM at NLO in ZM-VFNS with target-mass corrections --
    the reference itself -- evaluated member by member from its PineAPPL
    grids (analysis/npdf_impact.py).  earlier reweighted POWHEG-V2 events; the
    muon generator, POWHEG-RES, cannot be reweighted, and an analytic
    route treats both currents alike.  Every observable is leptonic.
  * THE REGION: Q2 > 4 GeV2 and W > 3 GeV, no y cut, 1 TeV, per nucleon of
    tungsten.
  * THE BASELINE is NNPDF4.0 averaged over the free nucleons of tungsten,
    (74 p + 110 n)/184 -- nuclear effects assumed to vanish, as in every 
    reference -- with its replica band.  Against it: nNNPDF3.0 (replicas),
    EPPS21 and nCTEQ15 (Hessian, rescaled from 90% to 68% CL by LHAPDF), all
    for tungsten and all already the average bound nucleon, so no isospin is
    applied on top (checked by the valence sum rules, a CLAIM below).
  * THE LOWER PANELS are ratios to EACH SET'S OWN free-nucleon baseline
    (user, 2026-09-30), each set with its own 68% band: the nuclear
    modification proper, without the difference between the proton fits.
    NNPDF4.0 on free nucleons is its own baseline and sits at one.  The
    CLAIMS still quote shifts against the common NNPDF4.0 baseline, as the
    text does, and say which.
  * BINS HOLDING LESS THAN 0.3% OF THE BASELINE RATE ARE NOT DRAWN (MIN_SHARE):
    EPPS21's ratio in the last x bin (x > 0.79, 0.07% of the neutrino rate)
    is 4.7, the large-x rise of its bound-nucleon valence, and on a common
    axis it would flatten every per-cent effect the figure is for.  The same
    rise lifts EPPS21 to 1.17 in the neutrino bin Q2 > 1000 GeV2 (0.26% of
    the rate), where only large x is reachable: physical, but a one-bin spike
    that reads as a fluctuation in a paper figure (user, 2026-09-21: no
    visible fluctuations), so the floor was raised from 0.1% to 0.3%.  It
    also trims the muon Q2 tail above 100 GeV2 and the two or three lowest-x
    bins of each current, each below 0.26% of the rate.  The calculation is
    analytic, so there is no Monte Carlo noise; SMOOTHNESS is a CLAIM.
  * EACH SET'S OWN FREE-NUCLEON BASELINE (nNNPDF3.0's proton fit, CT18A,
    nCTEQ15's proton) is computed and quoted in MESSAGE, to separate the
    nuclear modification proper from the proton-fit difference, but not
    drawn: the figure is about what a FASER prediction would change by.

GENIE IS ABSENT, under rule 1b's paper-plot carve-out: a nuclear-PDF effect is
a property of the parton distributions, not of a generator, and GENIE cannot
be re-run per PDF member (438 members per current).  The analytic reference
is the calculation every generator is compared against, so the ratios carry
over to them.  GENIE's own comparison with that reference is pp03/pp04.

Inputs: results{,_nu}/npdf_impact_q4w3_W_tmc.json (analysis/npdf_impact.py),
        results{,_nu}/histos_yadism_nlo_q4w3_W_tmc.json (the closure).

Usage: analysis/paper_plots/ppA2_npdf_impact.py
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE, "analysis"))

PAPER_SECTION = "appendix"
SLUG = "ppA2_npdf_impact"
V1_NUMBER = "A2"
IN_PAPER = True
TITLE = "Impact of nuclear PDFs"
OUTPUT = "pp_npdf_impact.png"
RESULTS = "results_nu"

SURFACE = "#ffffff"
FILES = {"mu": f"{BASE}/results/npdf_impact_q4w3_W_tmc.json",
         "nu": f"{BASE}/results_nu/npdf_impact_q4w3_W_tmc.json"}
REFS = {"mu": f"{BASE}/results/histos_yadism_nlo_q4w3_W_tmc.json",
        "nu": f"{BASE}/results_nu/histos_yadism_nlo_q4w3_W_tmc.json"}

BASELINE = "NNPDF40_nnlo_as_01180_W184free"
NUCLEAR = ["nNNPDF30_nlo_as_0118_A184_Z74", "EPPS21nlo_CT18Anlo_W184",
           "nCTEQ15HQ_FullNuc_184_74"]
LABEL = {BASELINE: "NNPDF4.0, free nucleons",
         "nNNPDF30_nlo_as_0118_A184_Z74": "nNNPDF3.0",
         "EPPS21nlo_CT18Anlo_W184": "EPPS21",
         "nCTEQ15HQ_FullNuc_184_74": "nCTEQ15HQ"}
# Okabe-Ito, with the baseline in near-black like every reference
COLOUR = {BASELINE: "#111111",
          "nNNPDF30_nlo_as_0118_A184_Z74": "#0072b2",
          "EPPS21nlo_CT18Anlo_W184": "#d55e00",
          "nCTEQ15HQ_FullNuc_184_74": "#009e73"}
LS = {BASELINE: "-", "nNNPDF30_nlo_as_0118_A184_Z74": "--",
      "EPPS21nlo_CT18Anlo_W184": "-.",
      "nCTEQ15HQ_FullNuc_184_74": (0, (5.0, 1.3, 1.0, 1.3, 1.0, 1.3))}
MIN_SHARE = 3e-3

OBS = [("xbj", r"$x_{\rm Bj}$", r"d\sigma/dx_{\rm Bj}"),
       ("Q2", r"$Q^2$  [GeV$^2$]", r"d\sigma/dQ^2"),
       ("y", r"$y$", r"d\sigma/dy")]
ROWS = [("mu", "Muon DIS", 1000.0, "nb"), ("nu", "Neutrino DIS", 1.0, "pb")]
UNIT_X = {"xbj": "", "Q2": "/GeV$^2$", "y": ""}

CAPTION = (
    "Nuclear parton distributions and the inclusive DIS distributions at "
    "1 TeV, per nucleon of tungsten, for Q&sup2; &gt; 4 GeV&sup2; and "
    "W &gt; 3 GeV with no y cut: Bjorken <i>x</i>, momentum transfer and "
    "inelasticity for muon neutral-current (top) and neutrino charged-current "
    "(bottom) scattering. YADISM (ZM-VFNS) at NLO with target-mass "
    "corrections, with NNPDF4.0 averaged over the free protons and neutrons "
    "of tungsten (no nuclear effects) and with the tungsten sets of "
    "nNNPDF3.0, EPPS21 and nCTEQ15HQ. Beneath each distribution, the ratio "
    "of each set to its own free-nucleon baseline (nNNPDF3.0's proton fit, "
    "CT18A, nCTEQ15HQ's free proton), i.e. the nuclear modification proper; "
    "the bands are each set's 68% confidence-level PDF uncertainty. Bins "
    "holding less than 0.3% of the rate are not shown.")

MESSAGE = """<b>Nuclear PDFs lower the fiducial rate on both currents, by more for muons.</b> Against NNPDF4.0 on free nucleons, whose own PDF uncertainty is 0.8% on either current, the muon neutral-current rate falls by 5.4% (nNNPDF3.0), 7.0% (EPPS21) and 9.1% (nCTEQ15HQ), and the neutrino charged-current rate by 2.3%, 3.7% and 0.6%. The nuclear-PDF bands are 2.3-4.3 times the free-nucleon one on the muon side and 2.8-5.9 times on the neutrino side, nNNPDF3.0's the widest (3.5% and 4.4%), so its neutrino shift is well within its own band.

<b>The two currents sample different nuclear effects.</b> Half of the muon rate lies below x = 0.04, in the shadowing region: there all three sets fall with x, to 0.66-0.86 of the free-nucleon NNPDF4.0 value in the lowest x bin, and the muon suppression grows with y, where x is smallest at fixed Q&sup2;. The neutrino rate has its median near x = 0.17, between antishadowing -- the sets rise by 2-5% at x of about 0.1 -- and the EMC suppression above x = 0.3, which largely cancel. In nNNPDF3.0 and EPPS21 the EMC suppression wins, leaving a deficit almost flat in y (3.5-3.8% in EPPS21 above y = 0.1); nCTEQ15HQ's is much weaker (0.98 of its baseline from x = 0.3), and its neutrino rate is within 1% of free nucleons.

<b>Part of each shift is the proton fit, not the nucleus.</b> Against their own free-nucleon baselines -- what the figure draws -- the nuclear modification proper is -4.8% (nNNPDF3.0), -6.0% (EPPS21) and -4.3% (nCTEQ15HQ) on the muon side and -2.1%, -1.8% and +0.7% on the neutrino side: 2-3 times larger for muons in the first two, while nCTEQ15HQ has none for neutrinos. nCTEQ15HQ's free-proton baseline alone sits 5.0% below NNPDF4.0 on the muon side. Against NNPDF4.0 the sets agree on the sign everywhere and disagree most where the data are weakest: at the lowest x the spread between them, 0.20 in ratio on the muon side, is more than twice any one set's band. The one step in the common-baseline ratio, EPPS21 turning up to 1.04 in the neutrino bin 0.62 < x < 0.79, is its proton fit and not the nucleus: its CT18A baseline is itself 1.13 of NNPDF4.0 there, and EPPS21 is 0.92 of CT18A, an ordinary EMC suppression."""


def _load(cur):
    with open(FILES[cur]) as f:
        return json.load(f)


def _ref(cur):
    with open(REFS[cur]) as f:
        return json.load(f)


def _integ(cur, name):
    return _load(cur)["sets"][name]["integrated"]


def _shift(cur, name):
    """sigma_fid(set) / sigma_fid(baseline) - 1."""
    return (_integ(cur, name)["sigma_fid_pb"]
            / _integ(cur, BASELINE)["sigma_fid_pb"] - 1.0)


def _band(cur, name):
    """Mean of the + and - 68% relative half-widths of sigma_fid."""
    v = _integ(cur, name)
    return 0.5 * (v["err_plus_pb"] + v["err_minus_pb"]) / v["sigma_fid_pb"]


def _own(cur, name):
    """The nuclear modification proper: set / its own free-nucleon baseline."""
    d = _load(cur)
    ob = d["own_baseline"][name]
    return (d["sets"][name]["integrated"]["sigma_fid_pb"]
            / d["sets"][ob]["integrated"]["sigma_fid_pb"] - 1.0)


def _obs(cur, name, obs):
    d = _load(cur)
    s = d["sets"][name]["observables"][obs]
    e = np.asarray(d["binning"][obs]["edges"])
    return e, np.asarray(s["dsig"]), np.asarray(s["err_plus"]), \
        np.asarray(s["err_minus"])


def _share(cur, obs):
    e, b, _u, _d = _obs(cur, BASELINE, obs)
    r = b * np.diff(e)
    return r / r.sum()


def _vis(cur, obs):
    return _share(cur, obs) >= MIN_SHARE


def _ratio(cur, name, obs):
    """(edges, ratio, ratio of +band, ratio of -band) to the baseline central."""
    e, b, _u, _d = _obs(cur, BASELINE, obs)
    _e, y, up, dn = _obs(cur, name, obs)
    return e, y / b, up / b, dn / b


def _own_ratio(cur, name, obs):
    """(edges, ratio, +band, -band) to the CENTRAL value of the set's OWN
    free-nucleon baseline -- what the lower panels draw (user, 2026-09-30:
    "switch to each fit's own baseline").  The free-nucleon NNPDF4.0 is its
    own baseline, so it sits at one with its band."""
    d = _load(cur)
    ob = d["own_baseline"].get(name, name)
    e = np.asarray(d["binning"][obs]["edges"])
    b = np.asarray(d["sets"][ob]["observables"][obs]["dsig"])
    _e, y, up, dn = _obs(cur, name, obs)
    return e, y / b, up / b, dn / b


def _median(cur, obs):
    """The baseline's median in one observable (log-interpolated for x, Q2)."""
    e = np.asarray(_load(cur)["binning"][obs]["edges"])
    c = np.concatenate([[0.0], np.cumsum(_share(cur, obs))])
    le = np.log(e) if obs in ("xbj", "Q2") else e
    v = float(np.interp(0.5, c, le))
    return float(np.exp(v)) if obs in ("xbj", "Q2") else v


def _closure(cur):
    """Largest |baseline member 0 / tracked reference - 1| over sigma_fid and
    every x, Q2, y bin -- re-derived from both files, not read from a flag."""
    d, r = _load(cur), _ref(cur)
    b = d["sets"][BASELINE]
    worst = abs(b["integrated"]["member0_pb"] / r["sigma_fid_pb"] - 1.0)
    for o, *_x in OBS:
        a = np.asarray(b["observables"][o]["dsig_member0"])
        t = np.asarray(r["hists"][o]["dsig"])
        m = t != 0
        worst = max(worst, float(np.max(np.abs(a[m] / t[m] - 1.0))))
    return worst


def _valence_ok():
    """Every set's valence integrals match the average nucleon of tungsten,
    ((2Z+N)/A, (Z+2N)/A), within 3% -- and not the proton's (2, 1)."""
    out = []
    for cur in ("mu", "nu"):
        d = _load(cur)
        z, n, a = d["target"]["Z"], d["target"]["N"], d["target"]["A"]
        want = ((2 * z + n) / a, (z + 2 * n) / a)
        for name in [BASELINE] + NUCLEAR + list(d["own_baseline"].values()):
            v = d["target_composition"]["sets"][name]
            out.append(max(abs(v["u_v"] - want[0]), abs(v["d_v"] - want[1])))
    return max(out)


def _lowest_x(cur):
    """The nuclear sets' ratios in the lowest drawn x bin."""
    i = int(np.nonzero(_vis(cur, "xbj"))[0][0])
    return [float(_ratio(cur, n, "xbj")[1][i]) for n in NUCLEAR]


def _peak(cur, name, lo, hi):
    """Largest ratio over drawn x bins with lower edge in [lo, hi)."""
    e, r, _u, _d = _ratio(cur, name, "xbj")
    m = _vis(cur, "xbj") & (e[:-1] >= lo) & (e[:-1] < hi)
    return float(r[m].max())


def _emc(name):
    """Neutrino ratios in the x bin starting at 0.30 and the one after."""
    e, r, _u, _d = _ratio("nu", name, "xbj")
    i = int(np.argmin(np.abs(e[:-1] - 0.30)))
    return float(r[i]), float(r[i + 1])


def _near(got, quoted, tol=0.06):
    """Each derived number rounds to the one quoted in MESSAGE (within tol,
    so that -8.95 counts as the -9.0 it is quoted as)."""
    return len(got) == len(quoted) and all(abs(g - q) <= tol
                                           for g, q in zip(got, quoted))


def _kink(cur, obs, name):
    """max |second difference| of one ratio curve over the drawn bins."""
    r = _ratio(cur, name, obs)[1][_vis(cur, obs)]
    return float(np.max(np.abs(np.diff(r, 2)))) if len(r) > 2 else 0.0


def _worst_kink(obs_list, skip_last_x=False):
    out = 0.0
    for cur in ("mu", "nu"):
        for o in obs_list:
            for name in _load(cur)["sets"]:
                if name == BASELINE:
                    continue
                r = _ratio(cur, name, o)[1][_vis(cur, o)]
                if skip_last_x and o == "xbj":
                    r = r[:-1]
                if len(r) > 2:
                    out = max(out, float(np.max(np.abs(np.diff(r, 2)))))
    return out


def _zigzag(v):
    """Largest amplitude of a zig-zag: three consecutive first differences
    alternating in sign, measured by the smallest of the three.  A turnover
    (one sign change) is shape; two in a row is a fluctuation."""
    d = np.diff(v)
    z = [min(abs(d[i]), abs(d[i + 1]), abs(d[i + 2]))
         for i in range(len(d) - 2)
         if np.sign(d[i]) != np.sign(d[i + 1])
         and np.sign(d[i + 1]) != np.sign(d[i + 2])]
    return max(z) if z else 0.0


def _worst_zigzag():
    """(amplitude, where) of the worst zig-zag over every drawn ratio curve
    AND band edge of the figure."""
    out = (0.0, "none")
    for cur in ("mu", "nu"):
        for o, *_x in OBS:
            vis = _vis(cur, o)
            for n in [BASELINE] + NUCLEAR:
                _e, r, u, d = _ratio(cur, n, o)
                for lab, v in (("central", r), ("upper", r + u),
                               ("lower", r - d)):
                    z = _zigzag(v[vis])
                    if z > out[0]:
                        out = (z, f"{cur} {o} {LABEL[n]} {lab}")
    return out


def _quad(cur):
    """The driver's refined-quadrature test (npdf_impact.quadrature_check)."""
    return _load(cur)["quadrature_check"]


def _last_x(name):
    """(ratio to the NNPDF4.0 baseline) in the last drawn neutrino x bin."""
    return float(_ratio("nu", name, "xbj")[1][_vis("nu", "xbj")][-1])


CLAIMS = [
    {"what": "the EPPS21 step in the last drawn neutrino x bin (1.04) is its "
             "CT18A proton baseline (1.13 of NNPDF4.0 there); EPPS21 over "
             "CT18A is 0.92",
     "check": lambda: (1.03 < _last_x("EPPS21nlo_CT18Anlo_W184") < 1.045
                       and 1.12 < _last_x("CT18ANLO_W184free") < 1.14
                       and 0.91 < _last_x("EPPS21nlo_CT18Anlo_W184")
                       / _last_x("CT18ANLO_W184free") < 0.93),
     "detail": lambda: "EPPS21 %.3f, CT18A %.3f, ratio %.3f" % (
         _last_x("EPPS21nlo_CT18Anlo_W184"), _last_x("CT18ANLO_W184free"),
         _last_x("EPPS21nlo_CT18Anlo_W184") / _last_x("CT18ANLO_W184free"))},
    {"what": "the fiducial integration is converged: a refined quadrature "
             "(+4 Gauss points, panels and y lattice halved) moves no bin "
             "holding 0.1% of the rate by 1e-4 and no ratio by 1e-5, for "
             "every set on both currents",
     "check": lambda: all(max(_quad(c)["bin"].values()) < 1e-4
                          and max(_quad(c)["ratio"].values()) < 1e-5
                          and set(_quad(c)["bin"]) == {BASELINE, *NUCLEAR}
                          for c in ("mu", "nu")),
     "detail": lambda: "; ".join(
         f"{c}: bin {max(_quad(c)['bin'].values()):.1e}, ratio "
         f"{max(_quad(c)['ratio'].values()):.1e}" for c in ("mu", "nu"))},
    {"what": "no bin-to-bin fluctuation: no drawn ratio curve or band edge "
             "zig-zags (two sign flips in a row) by more than 1e-3",
     "check": lambda: _worst_zigzag()[0] < 1e-3,
     "detail": lambda: "worst %.1e (%s)" % _worst_zigzag()},
    {"what": "no visible fluctuation: every drawn ratio curve is smooth in Q2 "
             "and y (bin-to-bin second difference below 6%), and in x apart "
             "from the last drawn bin, where the bound-nucleon valence rises "
             "at large x (Fermi motion)",
     "check": lambda: (_worst_kink(("Q2", "y")) < 0.06
                       and _worst_kink(("xbj",), skip_last_x=True) < 0.06),
     "detail": lambda: (f"Q2,y {_worst_kink(('Q2', 'y')):.4f}; x without the "
                        f"last bin {_worst_kink(('xbj',), True):.4f}; "
                        f"x with it {_worst_kink(('xbj',)):.4f}")},
    {"what": "both inputs are YADISM NLO in ZM-VFNS with TMC 3, the q4w3 "
             "region, tungsten per nucleon, 1 TeV",
     "check": lambda: all(
         _load(c)["calculation"]["scheme"] == "ZM-VFNS"
         and _load(c)["calculation"]["pto"] == 1
         and _load(c)["calculation"]["tmc"] == 3
         and _load(c)["selection"]["name"] == "q4w3"
         and _load(c)["target"]["name"] == "W"
         and _load(c)["energy_gev"] == 1000.0 for c in ("mu", "nu")),
     "detail": lambda: ", ".join(
         f"{c}: {_load(c)['selection']['label']}, {_load(c)['calculation']}"
         for c in ("mu", "nu"))},
    {"what": "closure: the baseline's member 0 reproduces the tracked "
             "YADISM (ZM-VFNS) NLO reference, sigma_fid and every x, Q2, y "
             "bin, to better than 1e-6",
     "check": lambda: max(_closure("mu"), _closure("nu")) < 1e-6,
     "detail": lambda: f"mu {_closure('mu'):.1e}, nu {_closure('nu'):.1e}"},
    {"what": "every set is the average nucleon of tungsten (valence sum "
             "rules (2Z+N)/A = 1.402, (Z+2N)/A = 1.598 within 0.03): no "
             "isospin applied twice",
     "check": lambda: _valence_ok() < 0.03,
     "detail": lambda: f"largest deviation {_valence_ok():.4f}"},
    {"what": "the free-nucleon NNPDF4.0 band is 0.8% on either current",
     "check": lambda: all(round(100 * _band(c, BASELINE), 1) == 0.8
                          for c in ("mu", "nu")),
     "detail": lambda: ", ".join(f"{c} {100*_band(c, BASELINE):.2f}%"
                                 for c in ("mu", "nu"))},
    {"what": "muon rate: -5.4% (nNNPDF3.0), -7.0% (EPPS21), -9.1% (nCTEQ15HQ)",
     "check": lambda: _near([100 * _shift("mu", n) for n in NUCLEAR],
                            [-5.4, -7.0, -9.1]),
     "detail": lambda: ", ".join(f"{100*_shift('mu', n):+.2f}%"
                                 for n in NUCLEAR)},
    {"what": "neutrino rate: -2.3%, -3.7%, -0.6%",
     "check": lambda: _near([100 * _shift("nu", n) for n in NUCLEAR],
                            [-2.3, -3.7, -0.6]),
     "detail": lambda: ", ".join(f"{100*_shift('nu', n):+.2f}%"
                                 for n in NUCLEAR)},
    {"what": "nuclear bands are 2.3-4.3x the free-nucleon one (muon) and "
             "2.8-5.9x (neutrino); nNNPDF3.0's the widest, 3.5% and 4.4%",
     "check": lambda: (
         _near([f(_band(c, n) / _band(c, BASELINE) for n in NUCLEAR)
                for c in ("mu", "nu") for f in (min, max)],
               [2.3, 4.3, 2.8, 5.9])
         and all(max(NUCLEAR, key=lambda n: _band(c, n)) == NUCLEAR[0]
                 for c in ("mu", "nu"))
         and round(100 * _band("mu", NUCLEAR[0]), 1) == 3.5
         and round(100 * _band("nu", NUCLEAR[0]), 1) == 4.4),
     "detail": lambda: "; ".join(
         f"{c}: " + ", ".join(f"{_band(c, n)/_band(c, BASELINE):.2f}x"
                              for n in NUCLEAR) for c in ("mu", "nu"))},
    {"what": "nNNPDF3.0's neutrino shift is within its own band",
     "check": lambda: abs(_shift("nu", NUCLEAR[0])) < _band("nu", NUCLEAR[0]),
     "detail": lambda: f"{100*_shift('nu', NUCLEAR[0]):+.2f}% vs "
                       f"{100*_band('nu', NUCLEAR[0]):.2f}%"},
    {"what": "half the muon rate lies below x = 0.04; the neutrino median "
             "is near x = 0.17",
     "check": lambda: 0.035 < _median("mu", "xbj") < 0.0405
     and 0.165 < _median("nu", "xbj") < 0.18,
     "detail": lambda: f"mu {_median('mu', 'xbj'):.4f}, "
                       f"nu {_median('nu', 'xbj'):.4f}"},
    {"what": "muon, lowest drawn x bin: all three sets at 0.66-0.86 of the "
             "free-nucleon NNPDF4.0 value, and falling with decreasing x",
     "check": lambda: (round(min(_lowest_x("mu")), 2) == 0.66
                       and round(max(_lowest_x("mu")), 2) == 0.86
                       and all(_ratio("mu", n, "xbj")[1][:6].tolist()
                               == sorted(_ratio("mu", n, "xbj")[1][:6])
                               for n in NUCLEAR)),
     "detail": lambda: ", ".join(f"{v:.3f}" for v in _lowest_x("mu"))},
    {"what": "muon suppression grows with y (each set lower in the last y "
             "bin than at y = 0.05-0.075)",
     "check": lambda: all(_ratio("mu", n, "y")[1][-1]
                          < _ratio("mu", n, "y")[1][2] - 0.05 for n in NUCLEAR),
     "detail": lambda: ", ".join(
         f"{_ratio('mu', n, 'y')[1][2]:.3f} -> {_ratio('mu', n, 'y')[1][-1]:.3f}"
         for n in NUCLEAR)},
    {"what": "antishadowing: the sets peak 2-5% up near x = 0.1 (the two "
             "Hessian sets on both currents, all three on the neutrino)",
     "check": lambda: all(0.015 < _peak(c, n, 0.05, 0.2) - 1 < 0.05
                          for c in ("mu", "nu") for n in NUCLEAR[1:])
     and 0.015 < _peak("nu", NUCLEAR[0], 0.05, 0.2) - 1 < 0.05,
     "detail": lambda: "; ".join(
         f"{c}: " + ", ".join(f"{_peak(c, n, 0.05, 0.2):.3f}" for n in NUCLEAR)
         for c in ("mu", "nu"))},
    {"what": "EMC suppression above x = 0.3 on the neutrino side: nNNPDF3.0 "
             "and EPPS21 below 0.95 in the bin from x = 0.30 and below 0.91 "
             "in the next; nCTEQ15HQ much weaker, 0.98 and 0.96",
     "check": lambda: (all(_emc(n)[0] < 0.95 and _emc(n)[1] < 0.91
                           for n in NUCLEAR[:2])
                       and _near(list(_emc(NUCLEAR[2])), [0.98, 0.96], 0.006)),
     "detail": lambda: ", ".join(f"{_emc(n)[0]:.3f}/{_emc(n)[1]:.3f}"
                                 for n in NUCLEAR)},
    {"what": "the neutrino deficit is almost flat in y above y = 0.1: "
             "3.5-3.8% in EPPS21, at most 1% in nCTEQ15HQ",
     "check": lambda: (
         3.5 <= round(100 * (1 - _ratio("nu", NUCLEAR[1], "y")[1][4:].max()), 1)
         and round(100 * (1 - _ratio("nu", NUCLEAR[1], "y")[1][4:].min()), 1) <= 3.8
         and round(100 * (1 - _ratio("nu", NUCLEAR[2], "y")[1][4:].min()), 1) <= 1.0),
     "detail": lambda: ", ".join(
         f"{min(_ratio('nu', n, 'y')[1][4:]):.3f}-"
         f"{max(_ratio('nu', n, 'y')[1][4:]):.3f}" for n in NUCLEAR)},
    {"what": "nuclear modification proper: muon -4.8/-6.0/-4.3%, neutrino "
             "-2.1/-1.8/+0.7%; 2-3x larger for muons in nNNPDF3.0 and "
             "EPPS21, while nCTEQ15HQ has none on the neutrino side",
     "check": lambda: _near([100 * _own("mu", n) for n in NUCLEAR],
                            [-4.8, -6.0, -4.3])
     and _near([100 * _own("nu", n) for n in NUCLEAR], [-2.1, -1.8, 0.7])
     and all(1.9 < _own("mu", n) / _own("nu", n) < 3.5 for n in NUCLEAR[:2]),
     "detail": lambda: "; ".join(
         f"{c}: " + ", ".join(f"{100*_own(c, n):+.2f}%" for n in NUCLEAR)
         for c in ("mu", "nu"))},
    {"what": "nCTEQ15HQ's free-proton baseline alone is 5.0% below NNPDF4.0 "
             "on the muon side",
     "check": lambda: round(100 * _shift(
         "mu", _load("mu")["own_baseline"][NUCLEAR[2]]), 1) == -5.0,
     "detail": lambda: f"{100*_shift('mu', _load('mu')['own_baseline'][NUCLEAR[2]]):+.2f}%"},
    {"what": "against NNPDF4.0 the sets agree on the sign of the integrated "
             "shift; at the lowest muon x their spread (0.20) is over twice "
             "any band",
     "check": lambda: all(_shift(c, n) < 0 for c in ("mu", "nu")
                          for n in NUCLEAR)
     and _near([max(_lowest_x("mu")) - min(_lowest_x("mu"))], [0.20], 0.006)
     and (max(_lowest_x("mu")) - min(_lowest_x("mu")))
     > 2.0 * max(0.5 * (_ratio("mu", n, "xbj")[2][_vis("mu", "xbj")][0]
                        + _ratio("mu", n, "xbj")[3][_vis("mu", "xbj")][0])
                 for n in NUCLEAR),
     "detail": lambda: f"spread {max(_lowest_x('mu')) - min(_lowest_x('mu')):.3f}"},
]


def main():
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=(15.0, 10.4))
    outer = fig.add_gridspec(2, 3, hspace=0.26, wspace=0.26, left=0.065,
                             right=0.99, top=0.865, bottom=0.07)
    for r, (cur, title, div, unit) in enumerate(ROWS):
        for c, (obs, xlab, ylab) in enumerate(OBS):
            inner = outer[r, c].subgridspec(2, 1, height_ratios=[1.75, 1.0],
                                            hspace=0.0)
            ax = fig.add_subplot(inner[0])
            ar = fig.add_subplot(inner[1], sharex=ax)
            vis = _vis(cur, obs)
            lo, hi = np.inf, -np.inf
            for name in [BASELINE] + NUCLEAR:
                e, y, up, dn = _obs(cur, name, obs)
                col = COLOUR[name]
                yv = np.where(vis, y / div, np.nan)
                ax.stairs(yv, e, color=col, lw=1.6, ls=LS[name])
                # each bin's band is its own span: fill_between(step="post")
                # stops at the last LEFT edge and would drop the last bin
                _e, rr, ru, rd = _own_ratio(cur, name, obs)
                rv = np.where(vis, rr, np.nan)
                ar.stairs(rv, e, color=col, lw=1.6, ls=LS[name])
                for i in np.nonzero(vis)[0]:
                    ar.fill_between([e[i], e[i + 1]], rr[i] - rd[i],
                                    rr[i] + ru[i], color=col, alpha=0.18,
                                    lw=0)
                    ax.fill_between([e[i], e[i + 1]], (y[i] - dn[i]) / div,
                                    (y[i] + up[i]) / div, color=col,
                                    alpha=0.18, lw=0)
                lo = min(lo, np.nanmin((rr - rd)[vis]))
                hi = max(hi, np.nanmax((rr + ru)[vis]))
            ar.axhline(1.0, color="#666666", lw=0.8, ls=":")
            e = np.asarray(_load(cur)["binning"][obs]["edges"])
            iv = np.nonzero(vis)[0]
            ax.set_xlim(e[iv[0]], e[iv[-1] + 1])
            if obs in ("xbj", "Q2"):
                ax.set_xscale("log")
            b = _obs(cur, BASELINE, obs)[1][vis] / div
            if b.max() / b.min() > 20.0:
                ax.set_yscale("log")
                ax.set_ylim(b.min() * 0.4, b.max() * 2.5)
            else:
                ax.set_ylim(0.0, b.max() * 1.25)
            pad = 0.08 * (hi - lo)
            ar.set_ylim(lo - pad, hi + pad)
            ax.set_ylabel(tex(f"${ylab}$ [{unit}{UNIT_X[obs]}]"),
                          fontsize=plotstyle.FS_YLABEL - 2)
            if c == 0:
                ar.set_ylabel(tex("ratio to own") + "\n" + tex("free nucleons"),
                              fontsize=plotstyle.FS_YLABEL - 4)
            ar.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL - 1)
            if c == 1:
                # larger (user, 2026-09-30: "increase font size in 'muon DIS'
                # and 'neutrino DIS'")
                ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE + 5)
            plotstyle.ticks(ax, labelbottom=False)
            plotstyle.ticks(ar)

    hs = [(Patch(color=COLOUR[n], alpha=0.18, lw=0),
           Line2D([], [], color=COLOUR[n], lw=1.6, ls=LS[n]))
          for n in [BASELINE] + NUCLEAR]
    labs = [tex(LABEL[n]) for n in [BASELINE] + NUCLEAR]
    # THE BEAM ENERGY IS IN THE TITLE (user, 2026-09-21), as on pp02/pp07:
    # both beams are at E_l = 1 TeV.  Legend enlarged on the same request.
    # (user, 2026-09-30)
    fig.suptitle(tex(r"Nuclear PDF corrections in Tungsten target, $E_{\ell} = 1$ TeV"),
                 y=0.995, fontsize=plotstyle.FS_SUPTITLE + 2)
    fig.legend(hs, labs, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=4, frameon=False, handlelength=3.0,
               fontsize=plotstyle.FS_LEGEND + 6)
    out = f"{BASE}/{RESULTS}/{OUTPUT}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
