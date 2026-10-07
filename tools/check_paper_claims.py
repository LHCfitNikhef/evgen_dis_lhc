#!/usr/bin/env python3
"""Verify the paper's quantitative statements against the results on disk.

The report has had such a checker since 2026-08-27 and it has caught prose
going stale half a dozen times -- most recently six numbers that moved when
the samples were tripled.  The PAPER had no equivalent, and it is the document
that gets published, so the same numbers were being maintained by hand in two
places with only one of them checked.

Each entry names the claim, the literal as it appears in the LaTeX, and a
callable that re-derives it from the result JSONs.  A tolerance is given
because the paper quotes rounded values.

Usage: check_paper_claims.py; exits non-zero if a claim no longer holds.
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))


def _ngr():
    with open(f"{BASE}/results_nu/faser_genie_rates.json") as f:
        return json.load(f)


def _npr():
    with open(f"{BASE}/results_nu/faser_powheg_rates_rwgt.json") as f:
        return json.load(f)


def _dmp():
    with open(f"{BASE}/results_nu/faser_dimuon_cutflow_pdf.json") as f:
        return json.load(f)


def _dmp_row(name, key, field="rel_to_nominal"):
    return _dmp()["sets"][name]["rows"]["sum"][key][field]


CLAIMS = [
    # ONLY THE FINAL REGION SURVIVES (user, 2026-09-19): every claim on the
    # 0.2 < y < 0.9 results, and every claim read through the deleted the earlier production
    # figure scripts in analysis/paper_plots/, was removed with them.
    # 2026-09-23: the "Expected event rates" section (flux reproduction,
    # published-prediction comparison, Run-3/4 table) was removed from the
    # paper when Secs. 7 and 8 merged into "Applications to FASER"; its
    # claims went with it.  The figures' own numbers quoted in Secs. 5-7 are
    # the MESSAGE literals of analysis/paper_plots/, checked by
    # tools/check_paper_plots.py.
]


# --- Appendix: impact of non-DIS processes (2026-09-04) --------------------
def _nd(cur):
    with open(f"{BASE}/{'results_nu' if cur != 'mu' else 'results'}/genie_nondis.json") as f:
        d = json.load(f)
    return d["antineutrino"] if cur == "nubar" else d


def _nd_share(cur, e, g="non-DIS", t="p"):
    b = _nd(cur)["benchmark"][t][f"{e:g}"]
    return 100.0 * (b["nondis_frac"] if g == "non-DIS" else b["frac"][g])


_NDV2 = {"nu400": "results_nu/genie_nondis_diff_W_400GeV.json",
         "nu": "results_nu/genie_nondis_diff_W.json",
         "mu400": "results/genie_nondis_diff_W_400GeV.json",
         "mu": "results/genie_nondis_diff_W.json"}


def _ndv2(k):
    with open(f"{BASE}/{_NDV2[k]}") as f:
        return json.load(f)


def _ndv2_share(k, reg):
    """Non-DIS percentage of the rate in region `reg` ("total" = no cut),
    on the tungsten combination of paper plot A1."""
    d = _ndv2(k)
    groups = ["DIS", "DIS charm", "RES", "QEL", "DFR"]
    a = {g: d["sigma_pb"][g] * (1.0 if reg == "total" else d["regions"][g][reg])
         for g in groups}
    return 100.0 * (a["RES"] + a["QEL"] + a["DFR"]) / sum(a.values())


_ND_TABLE = {   # the LaTeX table, (beam, energy) -> {channel: quoted %}
    ("mu", 400.0): {"DIS": 99.20, "RES": 0.75, "QEL": 0.05, "non-DIS": 0.80},
    ("mu", 1000.0): {"DIS": 99.42, "RES": 0.54, "QEL": 0.03, "non-DIS": 0.58},
    ("mu", 4000.0): {"DIS": 99.62, "RES": 0.36, "QEL": 0.02, "non-DIS": 0.38},
    ("nu", 400.0): {"DIS": 88.99, "DIS charm": 10.26, "RES": 0.58, "QEL": 0.01, "DFR": 0.16, "non-DIS": 0.75},
    ("nu", 1000.0): {"DIS": 88.00, "DIS charm": 11.68, "RES": 0.24, "QEL": 0.005, "DFR": 0.08, "non-DIS": 0.32},
    ("nu", 4000.0): {"DIS": 86.13, "DIS charm": 13.77, "RES": 0.07, "QEL": 0.001, "DFR": 0.03, "non-DIS": 0.09},
    ("nubar", 400.0): {"DIS": 90.20, "DIS charm": 8.28, "RES": 0.78, "QEL": 0.56, "DFR": 0.18, "non-DIS": 1.52},
    ("nubar", 1000.0): {"DIS": 89.56, "DIS charm": 9.81, "RES": 0.32, "QEL": 0.23, "DFR": 0.09, "non-DIS": 0.63},
    ("nubar", 4000.0): {"DIS": 87.95, "DIS charm": 11.88, "RES": 0.09, "QEL": 0.06, "DFR": 0.03, "non-DIS": 0.17},
}


def _nd_table_bad():
    bad = []
    for (cur, e), row in _ND_TABLE.items():
        for g, v in row.items():
            got = _nd_share(cur, e, g)
            tol = 0.0051
            if abs(got - v) > tol:
                bad.append(f"{cur} {e:g} {g}: table {v} vs {got:.4f}")
    return bad


CLAIMS += [
    ("the non-DIS table (all entries) matches the splines to the quoted precision",
     lambda: not _nd_table_bad(),
     lambda: "; ".join(_nd_table_bad()) or "all entries agree"),
    ("non-DIS share of the nu_mu total: 0.75 / 0.32 / 0.09% at 400 GeV / 1 TeV / 4 TeV",
     lambda: all(abs(_nd_share("nu", e) - v) < 0.006 for e, v in
                 ((400.0, 0.75), (1000.0, 0.32), (4000.0, 0.09))),
     lambda: ", ".join(f"{_nd_share('nu', e):.3f}" for e in (400.0, 1000.0, 4000.0))),
    ("the antineutrino non-DIS share is about twice the neutrino one",
     lambda: all(1.5 < _nd_share("nubar", e) / _nd_share("nu", e) < 2.5
                 for e in (400.0, 1000.0)),
     lambda: ", ".join(f"{_nd_share('nubar', e)/_nd_share('nu', e):.2f}" for e in (400.0, 1000.0))),
    ("resonances are 21% of the nu_mu proton rate at 10 GeV",
     lambda: abs(100.0 * _nd("nu")["targets"]["p"]["sigma_pb"]["RES"][0]
                 / _nd("nu")["targets"]["p"]["total_pb"][0] - 21.0) < 0.6
     and abs(_nd("nu")["energies_gev"][0] - 10.0) < 1e-6,
     lambda: f"{100.0 * _nd('nu')['targets']['p']['sigma_pb']['RES'][0] / _nd('nu')['targets']['p']['total_pb'][0]:.2f}%"),
    # moved with the 2026-10-07 flux switch (EPOS-LHC light + POWHEG charm,
    # was the 2021 average: 0.9% / 2.0%); not quoted in the paper text
    ("flux-weighted non-DIS share of the FASERnu rate: 1.07% nu_mu, 2.19% nubar_mu",
     lambda: abs(100 * _nd("nu")["flux_weighted_W"]["nondis"] - 1.07) < 0.06
     and abs(100 * _nd("nubar")["flux_weighted_W"]["nondis"] - 2.19) < 0.06,
     lambda: f"{100*_nd('nu')['flux_weighted_W']['nondis']:.3f}, {100*_nd('nubar')['flux_weighted_W']['nondis']:.3f}"),
    ("differential paragraph (tungsten): inside Q2 > 4 AND W > 3 no RES or "
     "QEL event on either current, no non-DIS muon event; the nu_mu non-DIS share "
     "of the fiducial rate 0.005% / 0.004% at 400 GeV / 1 TeV against 0.69% / "
     "0.29% of the total, and 0.08% / 0.03% above the Q2 floor alone",
     lambda: all(_ndv2(k)["regions"][g]["fiducial"] == 0.0
                 for k in _NDV2 for g in ("RES", "QEL"))
     and all(_ndv2(k)["regions"]["share_of_fiducial"]["non-DIS"] == 0.0
             for k in ("mu400", "mu"))
     and all(abs(_ndv2_share(k, reg) - v) < 0.0006 if reg == "fiducial"
             else abs(_ndv2_share(k, reg) - v) < 0.006
             for k, reg, v in (("nu400", "fiducial", 0.005), ("nu", "fiducial", 0.004),
                               ("nu400", "total", 0.69), ("nu", "total", 0.29),
                               ("nu400", "q2_above_floor", 0.08),
                               ("nu", "q2_above_floor", 0.03))),
     lambda: ", ".join(f"{k} {reg} {_ndv2_share(k, reg):.4f}%"
                       for k in ("nu400", "nu")
                       for reg in ("fiducial", "q2_above_floor", "total"))),
    ("the EM Q2 floor of this build is 4 GeV2, recorded in the result",
     lambda: _nd("mu")["em_q2_floor_gev2"] == 4.0,
     lambda: str(_nd("mu")["em_q2_floor_gev2"])),
]


def _dmc():
    with open(f"{BASE}/results_nu/faser_dimuon_cutflow.json") as f:
        return json.load(f)


CLAIMS += [
    ("FASER dimuon cut flow (Sec. faser-dimuon): POWHEG-V2 charm fraction 12.9% vs "
     "GENIE 8.6% at Q2>4 (8.9% with the 2021 flux, 9.8% before the 2026-10-01 charm fix), CC rates agree to "
     "0.4%, charm->mu step within 10%; p>100 GeV in 25 mrad 7.0 (GENIE) / "
     "11.5 (POWHEG-V2) events with the antineutrino (EPOS-LHC + POWHEG flux since "
     "2026-10-07; 10.5 / 17 with the 2021 average)",
     lambda: (abs(_dmc()["columns"]["powheg_q2"]["14"]["charm"]["cum_eff_pct"] - 12.9) < 0.05
              and abs(_dmc()["columns"]["genie_q2"]["14"]["charm"]["cum_eff_pct"] - 8.6) < 0.05
              and abs(_dmc()["columns"]["powheg_q2"]["14"]["cc"]["events"]
                      / _dmc()["columns"]["genie_q2"]["14"]["cc"]["events"] - 1) < 0.005
              and abs(_dmc()["columns"]["powheg_q2"]["14"]["charm_mu"]["step_eff_pct"]
                      / _dmc()["columns"]["genie_q2"]["14"]["charm_mu"]["step_eff_pct"] - 1) < 0.10
              and abs(_dmc()["columns"]["genie_q2"]["sum"]["p100"]["events"] - 7.0) < 0.1
              and abs(_dmc()["columns"]["powheg_q2"]["sum"]["p100"]["events"] - 11.5) < 0.5),
     lambda: "charm %.1f vs %.1f; CC ratio %.4f; step ratio %.3f; p100 %.1f / %.1f" % (
         _dmc()["columns"]["powheg_q2"]["14"]["charm"]["cum_eff_pct"],
         _dmc()["columns"]["genie_q2"]["14"]["charm"]["cum_eff_pct"],
         _dmc()["columns"]["powheg_q2"]["14"]["cc"]["events"]
         / _dmc()["columns"]["genie_q2"]["14"]["cc"]["events"],
         _dmc()["columns"]["powheg_q2"]["14"]["charm_mu"]["step_eff_pct"]
         / _dmc()["columns"]["genie_q2"]["14"]["charm_mu"]["step_eff_pct"],
         _dmc()["columns"]["genie_q2"]["sum"]["p100"]["events"],
         _dmc()["columns"]["powheg_q2"]["sum"]["p100"]["events"])),
]


def _pio(cur, key, field, sel="sidis_e"):
    """sidis_e: FASER Tier E nested in the benchmark region (== q4w3_faser_e)."""
    with open(f"{BASE}/results_nu/faser_pions.json") as f:
        return json.load(f)["currents"][cur]["generators"][key]["selections"][sel][field]


def _pio_zratio(cur, key, ref, lo):
    import numpy as np
    a, b = np.array(_pio(cur, key, "pi_emul_z")), np.array(_pio(cur, ref, "pi_emul_z"))
    sl = slice(0, 2) if lo else slice(10, None)
    return float(a[sl].sum() / b[sl].sum())


def _sidis_ladder():
    """Every SIDIS energy some neutrino generator other than POWHEG-V2 ran at."""
    gens = _json("results_nu/faser_pions.json")["currents"]["nu"]["generators"]
    return sorted({float(e) for k, g in gens.items() if k != "powheg_nu"
                   for e in g["energies"]})


def _json(rel):
    with open(f"{BASE}/{rel}") as f:
        return json.load(f)


def _pio_eff(cur, key, energy):
    """The Tier E efficiency of the tungsten record at one ladder energy."""
    import faser_pions as _fp
    with open(_fp.spectra_path(cur, key, energy, "W")) as f:
        s = json.load(f)["selections"]
    return s["sidis_e"]["sigma_pb"] / s["sidis_w"]["sigma_pb"]


CLAIMS += [
    ("single-inclusive pions (Sec. faser-sidis, NEUTRINO ONLY by user decision "
     "2026-09-07; benchmark region 2026-09-18): six-energy ladder from 300 GeV to 4 TeV; "
     "Tier E efficiency 0.48 / 0.14 / 0.03 at 1 TeV / 400 / 300 GeV",
     lambda: (len(_sidis_ladder()) == 6
              and min(_sidis_ladder()) == 300.0
              and max(_sidis_ladder()) == 4000.0
              and abs(_pio_eff("nu", "powheg_nu", 1000.0) - 0.48) < 0.005
              and abs(_pio_eff("nu", "powheg_nu", 400.0) - 0.14) < 0.005
              and abs(_pio_eff("nu", "powheg_nu", 300.0) - 0.03) < 0.005),
     lambda: "ladder %s; eff %.3f/%.3f/%.3f" % (
         _sidis_ladder(),
         _pio_eff("nu", "powheg_nu", 1000.0), _pio_eff("nu", "powheg_nu", 400.0),
         _pio_eff("nu", "powheg_nu", 300.0))),
]


# ---- GENIE differentially (Sec. inclusive, pp04 ; 2026-09-18) --------------
def _pp04_range(resdir, key, obs):
    """(min, max) of GENIE/FONLL over the bins pp04 DISPLAYS."""
    import importlib.util
    import numpy as np
    spec = importlib.util.spec_from_file_location(
        "pp04", f"{BASE}/analysis/paper_plots/pp04_genie_dis_distributions.py")
    pp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pp)
    e, ya, _ = pp.hist(resdir, key, obs)
    _e, yb, _ = pp.hist(resdir, pp.REF, obs)
    x0, x1 = pp.x_lo(obs, resdir, e, yb), pp.x_hi(obs, resdir, e)
    m = (e[:-1] >= x0 - 1e-12) & (e[1:] <= x1 + 1e-12) & (yb > 0)
    r = (ya / yb)[m]
    return float(r.min()), float(r.max())


CLAIMS += [
    ("GENIE differentially (Sec. inclusive): the default tune spans -16% to +22% "
     "(muon) and -15% to +8% (neutrino) of the reference in x_Bj, and reaches "
     "-6% (muon) and -10% (neutrino) in theta, over the bins the figure shows",
     lambda: (abs(_pp04_range("results", "genie_q4w3_W", "xbj")[0] - 0.84) < 0.006
              and abs(_pp04_range("results", "genie_q4w3_W", "xbj")[1] - 1.22) < 0.006
              and abs(_pp04_range("results_nu", "genie_lo_q4w3_W", "xbj")[0] - 0.85) < 0.006
              and abs(_pp04_range("results_nu", "genie_lo_q4w3_W", "xbj")[1] - 1.08) < 0.006
              and abs(_pp04_range("results", "genie_q4w3_W", "theta")[0] - 0.94) < 0.006
              and abs(_pp04_range("results_nu", "genie_lo_q4w3_W", "theta")[0] - 0.90) < 0.006),
     lambda: "x mu %.3f-%.3f nu %.3f-%.3f; theta mu %.3f nu %.3f" % (
         *_pp04_range("results", "genie_q4w3_W", "xbj"),
         *_pp04_range("results_nu", "genie_lo_q4w3_W", "xbj"),
         _pp04_range("results", "genie_q4w3_W", "theta")[0],
         _pp04_range("results_nu", "genie_lo_q4w3_W", "theta")[0])),
]


# ---- the NNLO SIDIS comparison (Sec. rates-sidis, 2026-09-07;  2026-09-18) --
_SIDIS = {}


def _sid():
    if not _SIDIS:
        with open(f"{BASE}/results_nu/faser_sidis.json") as f:
            _SIDIS.update(json.load(f))
    return _SIDIS


def _sd(cur, fl, key, region, field):
    return _sid()["currents"][cur]["flavours"][fl]["generators"][key][
        "regions"][region][field]


def _sd_gens():
    return list(_sid()["currents"]["nu"]["flavours"]["nu_mu"]["generators"])


def _sd_eff(cur, fl, key, region, parent="sidis_w"):
    return (_sd(cur, fl, key, region, "events")
            / _sd(cur, fl, key, parent, "events"))


def _sd_cmp(key, field):
    return _sid()["comparison"]["generators"][key][field]


def _sd_pres(cur, fl, key, region, kind):
    return _sd(cur, fl, key, region, "pion_presence")[kind]


def _sd_cluster():
    """How much the intra-event correlation inflates the statistical error in
    the first bin above the z cut."""
    import numpy as np
    r = _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_z_err")
    p = _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_z_err_poisson")
    edges = np.array(_sid()["z_edges"])
    i = int(np.argmax(edges[:-1] >= _sid()["z_min"] - 1e-12))
    return r[i] / p[i]


def _ghp():
    """The GENIE-HEDIS-PDF cross-check (analysis/genie_hedis_pdf.py)."""
    p = f"{BASE}/results_nu/genie_hedis_pdf.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def _sd_kpi(key, cur="nu", fl="nu_mu"):
    """Charged kaons as a percentage of the charged pions above the z cut."""
    return (100 * _sd(cur, fl, key, "sidis_e", "K_total_zcut")
            / _sd(cur, fl, key, "sidis_e", "pi_total_zcut"))


def _sd_chr(key, plus, minus, cur="nu", fl="nu_mu"):
    """A charge ratio above the z cut: h+/h-."""
    return (_sd(cur, fl, key, "sidis_e", f"{plus}_total_zcut")
            / _sd(cur, fl, key, "sidis_e", f"{minus}_total_zcut"))


def _sd_reach(cur, fl):
    """The upper edge of the last z bin whose statistical error is under 10%."""
    import numpy as np
    key = _sid()["reference"][cur]
    reg = _sid()["yield_region"]
    y = np.array(_sd(cur, fl, key, reg, "pi_z"))
    e = np.array(_sd(cur, fl, key, reg, "pi_z_err"))
    edges = np.array(_sid()["z_edges"])
    ok = [i for i in range(len(y))
          if edges[i] >= _sid()["z_min"] - 1e-12 and y[i] > 0
          and 100 * e[i] / y[i] < 10.0]
    return float(edges[max(ok) + 1]) if ok else float("nan")


CLAIMS += [
    # 2026-10-04: the pion and kaon yield entries that stood here were
    # removed -- the paper no longer quotes those numbers (its pion and kaon
    # paragraphs went in the 09-23/09-30 rewrites; the charged-hadron yield
    # it does quote is checked below), and with the antineutrino added the
    # yields moved.  The figures' own numbers stay checked in pp15/pp16.
    # Why GENIE HEDIS sits 5% low (Sec. inclusive, user question 2026-09-21).
    # The explanation is a whole calculation, so it is checked against the
    # JSON analysis/genie_hedis_pdf.py writes rather than retyped.
    ("GENIE HEDIS is its PDF (Sec. inclusive and Fig. genie-vs-energy): the "
     "HEDIS set's light quarks 4-8% below NNPDF4.0 over 0.02 < x < 0.5; the "
     "same YADISM (ZM-VFNS) calculation with it gives 0.947-0.951 from 400 GeV "
     "to 4 TeV against the 0.946-0.951 observed; the drawn FONLL curve carries "
     "the HEDIS points to 0.2% at every energy and 0.04% from 1 TeV up",
     lambda: (_ghp() is not None
              and _ghp()["hedis_pdf"] == "NNPDF31sx_nlo_as_0118_LHCb_nf_6"
              and all(0.92 <= v <= 0.96 for k, v in _ghp()["quark_ratio_x"].items()
                      if 0.02 <= float(k) <= 0.5 and float(k) != 0.1)
              and abs(min(r["pdf_only"] for r in _ghp()["rows"]) - 0.947) < 0.001
              and abs(max(r["pdf_only"] for r in _ghp()["rows"]) - 0.951) < 0.001
              and abs(min(r["ratio_published"] for r in _ghp()["rows"]) - 0.946) < 0.001
              and abs(max(r["ratio_published"] for r in _ghp()["rows"]) - 0.951) < 0.001
              and all(abs(r["genie_over_like_for_like"] - 1.0) < 0.001
                      for r in _ghp()["rows"])
              and all(abs(r["genie_over_fonll_same_pdf"] - 1.0) < 0.002
                      for r in _ghp()["rows"])
              and all(abs(r["genie_over_fonll_same_pdf"] - 1.0) < 0.0005
                      for r in _ghp()["rows"] if r["energy_gev"] >= 1000)
              and abs(_ghp()["fonll_factorisation_check"]
                      ["genie_over_direct_fonll"] - 1.0) < 0.0006),
     lambda: ", ".join(f"{r['energy_gev']:g} {r['pdf_only']:.4f}/"
                       f"{r['ratio_published']:.4f}->"
                       f"{r['genie_over_fonll_same_pdf']:.4f}"
                       for r in _ghp()["rows"])),
    ("dimuon PDF dependence (Sec. faser-dimuon): CC 6%, charm 22% (CT18 0.78, ATLAS 0.79, "
     "MSHT 0.90, GRV98 0.56); p>100 in 25 mrad 11.5+-0.5 / 9.6 / 11.5 / 10.6 / 10.9 / 6.6, "
     "spread 17%, largest band 16%, scale 4% (EPOS-LHC + POWHEG flux, 2026-10-07)",
     lambda: (abs(100 * (max(_dmp_row(n, "cc") for n in _dmp()["sets"])
                         / min(_dmp_row(n, "cc") for n in _dmp()["sets"]) - 1) - 6) < 0.5
              and abs(100 * (1 - min(_dmp_row(n, "charm") for n in ("CT18NNLO", "MSHT20nnlo_as118", "ATLASpdf21_T1", "ABMP16als118_5_nnlo"))) - 22) < 0.5
              and abs(_dmp_row("CT18NNLO", "charm") - 0.78) < 0.005
              and abs(_dmp_row("ATLASpdf21_T1", "charm") - 0.79) < 0.005
              and abs(_dmp_row("MSHT20nnlo_as118", "charm") - 0.90) < 0.005
              and abs(_dmp_row("GRV98lo", "charm") - 0.56) < 0.005
              and all(abs(_dmp_row(n, "p100", "events") - v) < 0.05 for n, v in (
                  ("NNPDF40_nnlo_as_01180", 11.5), ("CT18NNLO", 9.6), ("MSHT20nnlo_as118", 11.5),
                  ("ATLASpdf21_T1", 10.6), ("ABMP16als118_5_nnlo", 10.9), ("GRV98lo", 6.6)))
              and abs(_dmp_row("NNPDF40_nnlo_as_01180", "p100", "err_plus") - 0.46) < 0.05
              and abs(100 * (1 - min(_dmp_row(n, "p100") for n in ("CT18NNLO", "MSHT20nnlo_as118", "ATLASpdf21_T1", "ABMP16als118_5_nnlo"))) - 17) < 0.5
              and abs(100 * _dmp_row("CT18NNLO", "p100", "err_plus") / _dmp_row("CT18NNLO", "p100", "events") - 16) < 0.5
              and abs(100 * (_dmp()["scale"]["sum"]["p100"]["hi"] / _dmp()["scale"]["sum"]["p100"]["central"] - 1) - 4) < 0.5),
     lambda: "cc spread %.1f%%, charm min %.3f, p100 %s" % (
         100 * (max(_dmp_row(n, "cc") for n in _dmp()["sets"]) / min(_dmp_row(n, "cc") for n in _dmp()["sets"]) - 1),
         min(_dmp_row(n, "charm") for n in ("CT18NNLO", "MSHT20nnlo_as118", "ATLASpdf21_T1", "ABMP16als118_5_nnlo")),
         ", ".join("%.1f" % _dmp_row(n, "p100", "events") for n in _dmp()["sets"]))),
]


# --- Appendix: comparison with dedicated neutrino generators (2026-10-01) --
def _nug():
    """The figure script itself, so the text and the figure read one source."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_ppA3", os.path.join(BASE, "analysis", "paper_plots",
                              "ppA3_nu_generators.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CLAIMS += [
    ("App. nugen: free proton, sigma_fid below G18_02a by 7.4% (HEDIS), "
     "9.0% (NuWro), 9.9% (GiBUU); tungsten 5.9% (HEDIS), 11.8% (NuWro)",
     lambda: all(abs(-100 * _nug().rel(s, t) - v) < 0.06 for s, t, v in (
         ("nugen_genie_p", "p", 7.4), ("nuwro", "p", 9.0),
         ("gibuu", "p", 9.9), ("nugen_genie_W", "W", 5.9),
         ("nuwro_W", "W", 11.8))),
     lambda: ", ".join(f"{s} {-100 * _nug().rel(s, t):.2f}" for s, t in (
         ("nugen_genie_p", "p"), ("nuwro", "p"), ("gibuu", "p"),
         ("nugen_genie_W", "W"), ("nuwro_W", "W")))),
    ("App. nugen: free-proton mean N_ch 5.85, 5.70, 6.00, 5.87",
     lambda: all(abs(_nug().mean(s, "p", "nch") - v) < 0.005 for s, v in (
         ("nugen_genie_lo_p", 5.85), ("nugen_genie_p", 5.70),
         ("nuwro", 6.00), ("gibuu", 5.87))),
     lambda: ", ".join(f"{_nug().mean(s, 'p', 'nch'):.3f}" for s in (
         "nugen_genie_lo_p", "nugen_genie_p", "nuwro", "gibuu"))),
    ("App. nugen: tungsten mean N_ch 12.5 / 10.6 / 5.23 and N_p 5.2 / 4.0 / "
     "0.65 (G18_02a / NuWro / HEDIS)",
     lambda: all(abs(_nug().mean(s, "W", o) - v) < tol for s, o, v, tol in (
         ("nugen_genie_lo_W", "nch", 12.5, 0.05),
         ("nuwro_W", "nch", 10.6, 0.05),
         ("nugen_genie_W", "nch", 5.23, 0.005),
         ("nugen_genie_lo_W", "nprot", 5.2, 0.05),
         ("nuwro_W", "nprot", 4.0, 0.05),
         ("nugen_genie_W", "nprot", 0.65, 0.005))),
     lambda: "; ".join(f"{s} {_nug().mean(s, 'W', 'nch'):.3f} / "
                       f"{_nug().mean(s, 'W', 'nprot'):.3f}" for s in (
         "nugen_genie_lo_W", "nuwro_W", "nugen_genie_W"))),
]


class NeedsSamples(Exception):
    """A claim that re-derives its number from EVENT samples outside the
    repository (not from a tracked result) cannot run in a clone that lacks
    them: it is SKIPPED there, not failed.  On a machine with the samples it
    is checked as before."""


def _fe_remix():
    """The energy-remix piece of the emulsion band alone, in per cent."""
    sys.path.insert(0, os.path.join(BASE, "analysis"))
    import faser_emulsion_shapes as fes
    keep = fes._fixed_energy_band
    fes._fixed_energy_band = lambda obs, edges: None
    try:
        try:
            data = fes._data()
        except FileNotFoundError:
            raise NeedsSamples("the digitised FASER data (binning) are not "
                               "part of the public release")
        bins = {o: fes._edges(data, "numu", o) for o in fes.OBS_KEY}
        r = fes._fold("powheg_v2", "numu")
        if r is None:
            raise NeedsSamples("the POWHEG-V2 FASER ladder's event files "
                               "($POWHEG_V2/ladder-faser-v2) are not on this machine")
        m = fes._mhou("powheg_v2", "numu", bins, r["shapes"])
    finally:
        fes._fixed_energy_band = keep
    return 100 * max(max(abs(min(lo)), abs(max(hi))) for lo, hi in m.values())


CLAIMS += [
    # AND THE CAPTION'S CLAIM ABOUT THE BAND IT DRAWS.  The figure's band is
    # two pieces and the caption says which is exact; the exact one is the
    # energy remix, so that is the number checked here.  It is recomputed
    # rather than read back, so what is verified is the piece the caption
    # names and not the quadrature sum the figure shows.
    ("the emulsion figure's band (Sec. faser-emulsion): the energy-remix "
     "piece is computed exactly and stays below 0.7% in every bin "
     "(0.6% with the 2021 flux)",
     lambda: _fe_remix() < 0.7,
     lambda: "energy remix at most %.2f%%" % _fe_remix()),
]



# --- Sec. faser, SIDIS yields (2026-10-04: nu + nubar; 2026-10-06: 80% eff) -
def _sidis_h(fl):
    with open(f"{BASE}/results_nu/faser_sidis.json") as f:
        d = json.load(f)
    return d["currents"]["nu"]["flavours"][fl]["generators"]["powheg_nu"]["regions"][d["yield_region"]]


CLAIMS += [
    ("SIDIS yields (Sec. faser, Fig. sidis-yields): POWHEG-V2 expects around 900 (3700) "
     "charged hadrons above z = 0.1 for nu_e + nubar_e (nu_mu + nubar_mu) at 300 fb^-1, "
     "the antineutrino included (EPOS-LHC + POWHEG flux, 2026-10-07; 1900 / 4700 "
     "with the 2021 average)",
     lambda: (abs(_sidis_h("nu_e")["h_total_zcut"] - 900) < 30
              and abs(_sidis_h("nu_mu")["h_total_zcut"] - 3700) < 50
              and _sidis_h("nu_e")["events_nubar"] > 0 and _sidis_h("nu_mu")["events_nubar"] > 0),
     lambda: "nu_e %.0f, nu_mu %.0f; nubar share of events %.3f / %.3f" % (
         _sidis_h("nu_e")["h_total_zcut"], _sidis_h("nu_mu")["h_total_zcut"],
         _sidis_h("nu_e")["events_nubar"] / _sidis_h("nu_e")["events"],
         _sidis_h("nu_mu")["events_nubar"] / _sidis_h("nu_mu")["events"])),
]


# --- App. nondis, Table tab:kin-regions (Felix's comment, 2026-10-07) --------
def _kr(cur, where):
    with open(f"{BASE}/results_nu/genie_kinematic_regions.json") as f:
        d = json.load(f)[cur]
    return d["flux_weighted"] if where == "flux" else d["by_energy"][where]


_KR_TABLE = {   # region: (nu 100, nu 1000, nu flux, nubar 100, nubar 1000, nubar flux), in %
    "dis_benchmark": (79.5, 96.5, 90.8, 62.9, 93.4, 82.3),
    "dis_w2to3": (1.3, 0.14, 0.56, 2.3, 0.28, 0.97),
    "soft_dis": (15.5, 3.0, 7.2, 27.5, 5.7, 13.5),
    "nondis": (3.6, 0.31, 1.5, 7.3, 0.62, 3.2),
}
_KR_COLS = [("nu", "100"), ("nu", "1000"), ("nu", "flux"),
            ("nubar", "100"), ("nubar", "1000"), ("nubar", "flux")]


def _kr_ok():
    for reg, vals in _KR_TABLE.items():
        for (cur, w), v in zip(_KR_COLS, vals):
            got = 100 * _kr(cur, w)[reg]
            # each entry as printed: half a unit in its last digit
            nd = len(str(v).split(".")[1]) if "." in str(v) else 0
            if abs(got - v) > 0.5 * 10 ** -nd + 1e-9:
                return False
    return True


CLAIMS += [
    ("App. nondis, Table kin-regions: GENIE nu_mu / nubar_mu CC on tungsten in the "
     "benchmark region, 2 < W < 3, SIS and the resonance region W < 2, at 100 GeV, 1 TeV and flux-weighted "
     "(EPOS-LHC + POWHEG); text: 91% (82%) inside, SIS 7% (13%), resonance region 1.5% (3%), "
     "outside 3% (7%) at 1 TeV, 20% (37%) at 100 GeV, about half (three quarters) at "
     "30 GeV, the 2-3 GeV strip at most 5%; conclusions: 9% (18%) of the rate outside",
     lambda: (_kr_ok()
              and round(100 * _kr("nu", "flux")["dis_benchmark"]) == 91
              and round(100 * _kr("nubar", "flux")["dis_benchmark"]) == 82
              and round(100 * _kr("nu", "flux")["soft_dis"]) == 7
              and round(100 * _kr("nubar", "flux")["soft_dis"]) == 13
              and round(100 * (1 - _kr("nu", "1000")["dis_benchmark"])) == 3
              and round(100 * (1 - _kr("nubar", "1000")["dis_benchmark"])) == 7
              and round(100 * (1 - _kr("nu", "100")["dis_benchmark"])) == 20
              and round(100 * (1 - _kr("nubar", "100")["dis_benchmark"])) == 37
              and 0.4 < 1 - _kr("nu", "30")["dis_benchmark"] < 0.55
              and 0.7 < 1 - _kr("nubar", "30")["dis_benchmark"] < 0.8
              and max(_kr(c, e)["dis_w2to3"] for c in ("nu", "nubar")
                      for e in ("30", "100", "1000")) < 0.0525
              and round(100 * (1 - _kr("nu", "flux")["dis_benchmark"])) == 9
              and round(100 * (1 - _kr("nubar", "flux")["dis_benchmark"])) == 18),
     lambda: "flux nu %s; nubar %s; outside @30 GeV %.2f / %.2f" % (
         {k: round(100 * v, 2) for k, v in _kr("nu", "flux").items() if k != "mean_energy_gev"},
         {k: round(100 * v, 2) for k, v in _kr("nubar", "flux").items() if k != "mean_energy_gev"},
         1 - _kr("nu", "30")["dis_benchmark"], 1 - _kr("nubar", "30")["dis_benchmark"])),
]


# --- App. npdf, Fig. npdf-hadron (Felix's comment, 2026-10-07; nu only) ------
def _nph():
    """The figure script itself, so the text and the figure read one source."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_ppA2b", os.path.join(BASE, "analysis", "paper_plots",
                               "ppA2b_npdf_hadron.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CLAIMS += [
    ("App. npdf, Fig. npdf-hadron: nuclear modifications of the hadron-level "
     "observables at most 4% in Delta phi and below 3% elsewhere, inside nPDF "
     "bands of 2-5%",
     lambda: (round(_nph()._maxdev("dphix")) <= 4
              and max(_nph()._maxdev(o) for o in ("nch05", "Elead", "Emu")) < 3.0
              and all(2.0 <= x <= 5.0 for n in _nph().NUCLEAR for x in _nph()._half(n))),
     lambda: "max dev: " + ", ".join(f"{o} {_nph()._maxdev(o):.2f}%" for o in
                                     ("nch05", "Elead", "dphix", "Emu"))
     + "; half-bands " + ", ".join("+%.2f/-%.2f" % _nph()._half(n) for n in _nph().NUCLEAR)),
]


def main():
    nbad = nskip = 0
    for what, check, detail in CLAIMS:
        try:
            ok = bool(check())
            d = detail()
        except NeedsSamples as exc:
            print(f"  SKIP   {what}\n         {exc}")
            nskip += 1
            continue
        except Exception as exc:                              # noqa: BLE001
            print(f"  ERROR  {what}\n         {type(exc).__name__}: {exc}")
            nbad += 1
            continue
        print(f"  {'ok  ' if ok else 'FAIL'}   {what}   [{d}]")
        if not ok:
            nbad += 1
    print()
    if nbad:
        print(f"{nbad} paper claim(s) no longer hold. Fix the paper, or this "
              f"table if the sentence was rewritten.")
        return 1
    print(f"all {len(CLAIMS) - nskip} paper claim(s) hold"
          + (f" ({nskip} skipped: samples not on this machine)" if nskip else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

