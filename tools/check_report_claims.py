#!/usr/bin/env python3
"""Verify every number quoted in the report's PROSE against the results.

WHY THIS EXISTS.  The report's blurbs carry roughly a hundred hand-written
numbers -- cross-sections, ratios, K-factors, charm fractions.  Every one of
them was true when it was written, and there is nothing to keep them true
afterwards: regenerate a sample, fix a matrix-element bug, subtract a
negative NLO half, and the tables update themselves while the prose beside
them does not.  A second reading of the page found seven such claims, one of
them flatly wrong -- the neutrino LO blurb still said Herwig sat "2% low"
long after a polarised-beam card had moved it to 1.0005.  Nothing errored,
and the sentence sat directly above a table contradicting it.

This is the prose counterpart of analyze.py's closure gate: the gate stops a
SAMPLE being trusted because its cross-section looks right, and this stops a
SENTENCE being trusted for the same reason.

Each claim below names what it asserts, the number in the prose, and how to
recompute it from results/.  Add a claim whenever a number is written into a
blurb; a claim with no matching prose is reported too, so a deleted sentence
does not leave a check behind pretending to guard it.

Usage: tools/check_report_claims.py [--verbose]
Exit status 1 if any claim fails.
"""
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")    # faser_pions, for _pio_w

RESULTS = f"{BASE}/results"
RESULTS_NU = f"{BASE}/results_nu"
REPORT_SRC = f"{BASE}/analysis/make_report.py"


def _dmc():
    with open(f"{RESULTS_NU}/faser_dimuon_cutflow.json") as f:
        return json.load(f)


def _dmp():
    with open(f"{RESULTS_NU}/faser_dimuon_cutflow_pdf.json") as f:
        return json.load(f)


def _dmp_row(name, key, field="rel_to_nominal", who="sum"):
    return _dmp()["sets"][name]["rows"][who][key][field]


def _dmc_row(col, who, key, field):
    return _dmc()["columns"][col][who][key][field]


def _pio(cur, key, field, sel="sidis_e"):
    """sidis_e is FASER Tier E nested in the benchmark region (== q4w3_faser_e)."""
    with open(f"{RESULTS_NU}/faser_pions.json") as f:
        return json.load(f)["currents"][cur]["generators"][key]["selections"][sel][field]


def _pio_w(cur, key, energy):
    """The tungsten per-event record of one ladder energy."""
    import faser_pions as _fp
    with open(_fp.spectra_path(cur, key, energy, "W")) as f:
        return json.load(f)["selections"]


def _pio_eff(cur, key, energy):
    s = _pio_w(cur, key, energy)
    return s["sidis_e"]["sigma_pb"] / s["sidis_w"]["sigma_pb"]


def _pio_sigma(cur, key, energy):
    return _pio_w(cur, key, energy)["sidis_e"]["sigma_pb"]


def _pio_zratio(cur, key, ref, lo):
    """Pions below z = 0.1 (lo=True) or above 0.5 (lo=False), over the reference."""
    import numpy as np
    a, b = np.array(_pio(cur, key, "pi_emul_z")), np.array(_pio(cur, ref, "pi_emul_z"))
    sl = slice(0, 2) if lo else slice(10, None)
    return float(a[sl].sum() / b[sl].sum())


def _json_field(path, field):
    with open(path) as f:
        return json.load(f)[field]


_SID = {}


def _sid():
    if not _SID:
        with open(f"{RESULTS_NU}/faser_sidis.json") as f:
            _SID.update(json.load(f))
    return _SID


def _sd(cur, fl, key, region, field):
    """One field of the SIDIS yield record."""
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
    """How much the intra-event correlation inflates the statistical error,
    in the first bin above the z cut. It is the same on every beam."""
    import numpy as np
    r = _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_z_err")
    p = _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_z_err_poisson")
    zc = _sid()["z_min"]
    edges = np.array(_sid()["z_edges"])
    i = int(np.argmax(edges[:-1] >= zc - 1e-12))
    return r[i] / p[i]


def _sd_reach(cur, fl):
    """The upper edge of the last z bin whose statistical error is under 10%."""
    import numpy as np
    key = _sid()["reference"][cur]
    reg = _sid()["yield_region"]
    y = np.array(_sd(cur, fl, key, reg, "pi_z"))
    e = np.array(_sd(cur, fl, key, reg, "pi_z_err"))
    edges = np.array(_sid()["z_edges"])
    zc = _sid()["z_min"]
    ok = [i for i in range(len(y))
          if edges[i] >= zc - 1e-12 and y[i] > 0 and 100 * e[i] / y[i] < 10.0]
    return float(edges[max(ok) + 1]) if ok else float("nan")


# --- the NOMAD charm-scheme scan (2026-09-08) ------------------------------
def _load(rel):
    with open(f"{BASE}/{rel}") as f:
        return json.load(f)


def _nom():
    return _load("results_nu/nomad_dimuon_fonll.json")


def _nom_chi2(key):
    d = _nom()
    return (d["chi2_per_point"] if key == "r_th"
            else d[f"chi2_{key}"] / d["ndat"])


def _nom_ratio(key, which):
    p = _nom()["points"]
    r = [x[key] / x["r"] for x in p]
    return r[0] if which == "first" else r[-1]


def _donut_flux_frac(flavour):
    """One flavour's share of DONUT's interacting-neutrino spectrum, in %."""
    with open(f"{BASE}/data/donut/flux.json") as f:
        d = json.load(f)["spectra"]
    tot = sum(sum(v) for v in d.values())
    return 100.0 * sum(d[flavour]) / tot


def _sidis_common_ladder():
    """The SIDIS ladder of the generator comparison: every energy some
    neutrino generator OTHER than POWHEG-V2 was run at (GENIE's default tune
    stops at 1 TeV, so this is a union, not an intersection)."""
    with open(f"{RESULTS_NU}/faser_pions.json") as f:
        gens = json.load(f)["currents"]["nu"]["generators"]
    return sorted({float(e) for k, g in gens.items() if k != "powheg_nu"
                   for e in g["energies"]})


def _eff():
    """The Tier E efficiency factors (analysis/faser_sidis_efficiency.py)."""
    p = f"{RESULTS_NU}/faser_sidis_efficiency.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def _eff_z(fl, h, field, z):
    """One efficiency factor in the bin containing z, the mean of the five
    generators -- the value the figure and the table quote (user,
    2026-09-21)."""
    d = _eff()
    if d is None:
        return None
    import bisect
    i = bisect.bisect_right(d["z_edges"], z + 1e-12) - 1
    return d["flux_averaged_mean"][fl][h][field][i]


def _eff_event(fl, gen=None):
    """The flux-averaged Tier E event efficiency: the five-generator mean,
    or one generator's (the POWHEG-V2 yield prose quotes its own)."""
    d = _eff()
    if d is None:
        return None
    if gen:
        return d["flux_averaged"][fl][gen]["tier_e_event_efficiency"]
    return d["flux_averaged_mean"][fl]["tier_e_event_efficiency"]


def _eff_at_energy(e):
    d = _eff()
    return None if d is None else \
        d["per_energy"]["powheg_nu"][f"{e:g}"]["tier_e_event_efficiency"]


def _flux_mean(pid):
    import numpy as np
    sys.path.insert(0, f"{BASE}/analysis")
    import faser_rates as fr
    e, phi, _n = fr.flux_weights("nu", pid)
    e, phi = np.asarray(e), np.asarray(phi)
    return float((e * phi).sum() / phi.sum())


def claims():
    M, N = RESULTS, RESULTS_NU
    return [
        # ONLY THE FINAL REGION SURVIVES (user, 2026-09-19): every claim on
        # the 0.2 < y < 0.9 results and on the deleted no-y-cut proton studies
        # went with the report sections that quoted them.
        # --- DONUT (user, 2026-09-10) ------------------------------------
        ("DONUT located events", "578",
         json.load(open(f"{BASE}/data/donut/nch.json"))["n_located"], 0.5),
        ("DONUT beam: nu_mu share [%]", "57", _donut_flux_frac("nu_mu"), 1.0),
        ("DONUT beam: nu_e share [%]", "38", _donut_flux_frac("nu_e"), 1.0),
        ("DONUT beam: nu_tau share [%]", "5", _donut_flux_frac("nu_tau"), 1.0),
        # --- the kaon share of the pion study (2026-09-04) ---------------
        ("pions: kaon over pion fraction, POWHEG-V2", "12",
         100 * _pio("nu", "powheg_nu", "K_emul_total") / _pio("nu", "powheg_nu", "pi_emul_total"), 0.5),
        ("pions: kaon over pion fraction, GENIE", "15",
         100 * _pio("nu", "genie_lo", "K_emul_total") / _pio("nu", "genie_lo", "pi_emul_total"), 0.5),
        # --- the PDF dependence of the dimuon rates (2026-09-07) ----------------
        ("dimuon PDF: nominal closure", "1.0000",
         max(abs(v - 1) for v in _dmp()["closure_vs_cutflow"].values()) + 1.0, 5e-5),
        ("dimuon PDF: CC row spread", "6.2",
         100 * (max(_dmp_row(n, "cc") for n in _dmp()["sets"])
                / min(_dmp_row(n, "cc") for n in _dmp()["sets"]) - 1), 0.1),
        ("dimuon PDF: CC row min", "0.956", min(_dmp_row(n, "cc") for n in _dmp()["sets"]), 0.001),
        ("dimuon PDF: CC row max", "1.015", max(_dmp_row(n, "cc") for n in _dmp()["sets"]), 0.001),
        ("dimuon PDF: NNPDF band on CC", "0.5",
         100 * _dmp_row("NNPDF40_nnlo_as_01180", "cc", "err_plus") / _dmp_row("NNPDF40_nnlo_as_01180", "cc", "events"), 0.1),
        ("dimuon PDF: charm row spread among the NNLO sets", "22",
         100 * (1 - min(_dmp_row(n, "charm") for n in ("CT18NNLO", "MSHT20nnlo_as118", "ATLASpdf21_T1", "ABMP16_5_nnlo"))), 0.5),
        ("dimuon PDF: CT18 charm", "0.78", _dmp_row("CT18NNLO", "charm"), 0.005),
        ("dimuon PDF: ATLAS charm", "0.80", _dmp_row("ATLASpdf21_T1", "charm"), 0.005),
        ("dimuon PDF: ABMP16 charm", "0.86", _dmp_row("ABMP16_5_nnlo", "charm"), 0.005),
        ("dimuon PDF: MSHT20 charm", "0.91", _dmp_row("MSHT20nnlo_as118", "charm"), 0.005),
        ("dimuon PDF: GRV98 charm", "0.56", _dmp_row("GRV98lo", "charm"), 0.005),
        ("dimuon PDF: NNPDF p100", "17.3", _dmp_row("NNPDF40_nnlo_as_01180", "p100", "events"), 0.05),
        ("dimuon PDF: NNPDF p100 band", "0.7", _dmp_row("NNPDF40_nnlo_as_01180", "p100", "err_plus"), 0.05),
        ("dimuon PDF: CT18 p100", "14.3", _dmp_row("CT18NNLO", "p100", "events"), 0.05),
        ("dimuon PDF: MSHT20 p100", "17.1", _dmp_row("MSHT20nnlo_as118", "p100", "events"), 0.05),
        ("dimuon PDF: ATLAS p100", "15.7", _dmp_row("ATLASpdf21_T1", "p100", "events"), 0.05),
        ("dimuon PDF: ABMP16 p100", "16.0", _dmp_row("ABMP16_5_nnlo", "p100", "events"), 0.05),
        ("dimuon PDF: GRV98 p100", "9.9", _dmp_row("GRV98lo", "p100", "events"), 0.05),
        ("dimuon PDF: p100 spread among the NNLO sets", "17",
         100 * (1 - min(_dmp_row(n, "p100") for n in ("CT18NNLO", "MSHT20nnlo_as118", "ATLASpdf21_T1", "ABMP16_5_nnlo"))), 0.5),
        ("dimuon PDF: largest band on p100 (CT18)", "15",
         100 * _dmp_row("CT18NNLO", "p100", "err_plus") / _dmp_row("CT18NNLO", "p100", "events"), 0.5),
        ("dimuon PDF: NNPDF band on p100", "4",
         100 * _dmp_row("NNPDF40_nnlo_as_01180", "p100", "err_plus") / _dmp_row("NNPDF40_nnlo_as_01180", "p100", "events"), 0.5),
        ("dimuon PDF: scale envelope on p100, down", "3.6",
         100 * (1 - _dmp()["scale"]["sum"]["p100"]["lo"] / _dmp()["scale"]["sum"]["p100"]["central"]), 0.1),
        ("dimuon PDF: scale envelope on p100, up", "4.0",
         100 * (_dmp()["scale"]["sum"]["p100"]["hi"] / _dmp()["scale"]["sum"]["p100"]["central"] - 1), 0.1),
        # --- NOMAD: the charm-scheme scan (2026-09-08) --------------------
        ("NOMAD: chi2/N with massless charm", "216", _nom_chi2("r_zm"), 2.0),
        ("NOMAD: chi2/N with FONLL", "126", _nom_chi2("r_th"), 1.5),
        ("NOMAD: chi2/N with massive nf=3 charm", "53",
         _nom_chi2("r_ffns"), 1.5),
        ("NOMAD: FONLL theory/data at the lowest energy", "2.16",
         _nom_ratio("r_th", "first"), 0.03),
        ("NOMAD: FONLL theory/data at the highest energy", "1.37",
         _nom_ratio("r_th", "last"), 0.03),
        ("NOMAD: massless theory/data, lowest", "2.54",
         _nom_ratio("r_zm", "first"), 0.03),
        ("NOMAD: massless theory/data, highest", "1.45",
         _nom_ratio("r_zm", "last"), 0.03),
        ("NOMAD: massive theory/data, lowest", "1.36",
         _nom_ratio("r_ffns", "first"), 0.03),
        ("NOMAD: massive theory/data, highest", "1.32",
         _nom_ratio("r_ffns", "last"), 0.03),
        # --- single-inclusive pions at FASER, region (2026-09-18) -----------
        ("pions: POWHEG-V2 Tier E events", "2313", _pio("nu", "powheg_nu", "events"), 1.0),
        ("pions: POWHEG-V2 pions per event, emulsion", "6.88",
         _pio("nu", "powheg_nu", "pi_emul_per_event"), 0.006),
        ("pions: POWHEG-V2 pions, thousand (60% selection efficiency)", "9.5",
         _pio("nu", "powheg_nu", "pi_emul_total") / 1e3, 0.05),
        ("pions: Herwig pions per event", "7.37",
         _pio("nu", "herwig_nlo_full", "pi_emul_per_event"), 0.006),
        ("pions: Sherpa pions per event", "7.11",
         _pio("nu", "sherpa_nlo", "pi_emul_per_event"), 0.006),
        ("pions: GENIE default pions per event", "5.80",
         _pio("nu", "genie_lo", "pi_emul_per_event"), 0.006),
        ("pions: GENIE HEDIS pions per event", "5.85",
         _pio("nu", "genie", "pi_emul_per_event"), 0.006),
        ("pions: spread of the multiplicity", "1.27",
         _pio("nu", "herwig_nlo_full", "pi_emul_per_event")
         / _pio("nu", "genie_lo", "pi_emul_per_event"), 0.01),
        ("pions: GENIE pions, thousand (60% selection efficiency)", "6.7",
         _pio("nu", "genie_lo", "pi_emul_total") / 1e3, 0.05),
        ("pions: GENIE/POWHEG-V2 below z=0.1", "0.67",
         _pio_zratio("nu", "genie_lo", "powheg_nu", True), 0.005),
        ("pions: GENIE/POWHEG-V2 above z=0.5", "1.47",
         _pio_zratio("nu", "genie_lo", "powheg_nu", False), 0.005),
        ("pions: Sherpa/POWHEG-V2 below z=0.1", "1.16",
         _pio_zratio("nu", "sherpa_nlo", "powheg_nu", True), 0.005),
        ("pions: Sherpa/POWHEG-V2 above z=0.5", "0.72",
         _pio_zratio("nu", "sherpa_nlo", "powheg_nu", False), 0.005),
        ("pions: Herwig/POWHEG-V2 below z=0.1", "1.10",
         _pio_zratio("nu", "herwig_nlo_full", "powheg_nu", True), 0.005),
        ("pions: Herwig/POWHEG-V2 above z=0.5", "0.94",
         _pio_zratio("nu", "herwig_nlo_full", "powheg_nu", False), 0.005),
        ("pions: kaons, thousand, low", "1.0",
         min(_pio("nu", k, "K_emul_total") for k in
             ("powheg_nu", "herwig_nlo_full", "sherpa_nlo", "genie_lo", "genie")) / 1e3, 0.05),
        ("pions: kaons, thousand, high", "1.4",
         max(_pio("nu", k, "K_emul_total") for k in
             ("powheg_nu", "herwig_nlo_full", "sherpa_nlo", "genie_lo", "genie")) / 1e3, 0.05),
        # the ladder the note describes: six energies, 300 GeV the lowest
        # the ladder EVERY generator has; POWHEG-V2 alone also runs
        # down to 20 GeV for the region-only yield (beams.SIDIS_ENERGIES_LOW)
        ("pions: SIDIS ladder points", "6", len(_sidis_common_ladder()), 1e-9),
        ("pions: lowest ladder point, GeV", "300", min(_sidis_common_ladder()), 1e-9),
        # the tier efficiencies the note quotes (POWHEG-V2, tungsten)
        ("pions: Tier E efficiency at 1 TeV", "0.48", _pio_eff("nu", "powheg_nu", 1000.0), 0.005),
        ("pions: Tier E efficiency at 400 GeV", "0.14", _pio_eff("nu", "powheg_nu", 400.0), 0.005),
        ("pions: Tier E efficiency at 300 GeV", "0.03", _pio_eff("nu", "powheg_nu", 300.0), 0.005),
        # --- the luminosity of the "Predictions for FASER" tab -------------
        # The tab intro once said "everything is quoted at 250 fb^-1" while
        # three of the five sub-tabs were at 300; these pin each sub-tab's
        # JSON to the number the intro names, so it cannot drift again.
        ("FASER tab: dimuon cut-flow luminosity", "300",
         _json_field(f"{RESULTS_NU}/faser_dimuon_cutflow.json", "lumi_fb"), 1e-9),
        ("FASER tab: pion study luminosity", "300",
         _json_field(f"{RESULTS_NU}/faser_pions.json", "lumi_fb"), 1e-9),
        ("FASER tab: SIDIS study luminosity", "300",
         _json_field(f"{RESULTS_NU}/faser_sidis.json", "lumi_fb"), 1e-9),
        # A CLOSURE ACROSS ANALYSES, not a restatement: the pion pass's Tier E
        # tungsten cross-section at 1 TeV must equal the tracked  Tier E
        # result of the same sample (histos_*_q4w3_faser_e_W.json), which was
        # analysed by a different script.
        ("pions: Tier E closes against the Tier E result, nu", "1.0000",
         _pio_sigma("nu", "powheg_nu", 1000.0)
         / _json_field(f"{RESULTS_NU}/histos_powheg_nu_q4w3_faser_e_W.json", "sigma_fid_pb"), 1e-4),
        ("pions: Tier E closes against the Tier E result, mu", "1.0000",
         _pio_sigma("mu", "powheg", 1000.0)
         / _json_field(f"{RESULTS}/histos_powheg_q4w3_faser_e_W.json", "sigma_fid_pb"), 1e-4),
        # --- charged pions in z against arXiv:2504.05376, region (2026-09-18)
        # No y window and no x cut on any yield (user, 2026-09-07).
        ("sidis: x>0.1 would keep, nu Tier E", "0.72",
         _sd_eff("nu", "nu_mu", "powheg_nu", "sidis_ex", "sidis_e"), 0.005),
        ("sidis: nu_mu Tier E events", "2813",
         _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "events"), 1.0),
        ("sidis: nu_mu pi+ above z=0.1", "1513",
         _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pip_total_zcut"), 1.0),
        ("sidis: nu_mu pi- above z=0.1", "1167",
         _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pim_total_zcut"), 1.0),
        ("sidis: nu_mu pions above z=0.1", "2680",
         _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_total_zcut"), 1.0),
        ("sidis: nu_mu statistical error", "2.1",
         100 * _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_total_zcut_err")
         / _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pi_total_zcut"), 0.05),
        ("sidis: nu_e events", "1154",
         _sd("nu", "nu_e", "powheg_nu", "sidis_e", "events"), 1.0),
        ("sidis: nu_e pions above z=0.1", "1074",
         _sd("nu", "nu_e", "powheg_nu", "sidis_e", "pi_total_zcut"), 1.0),
        ("sidis: nu_mu events across the generators, low", "2292",
         min(_sd("nu", "nu_mu", k, "sidis_e", "events") for k in _sd_gens()), 1.0),
        ("sidis: nu_mu events across the generators, high", "3156",
         max(_sd("nu", "nu_mu", k, "sidis_e", "events") for k in _sd_gens()), 1.0),
        ("sidis: nu_mu pions across the generators, low", "2091",
         min(_sd("nu", "nu_mu", k, "sidis_e", "pi_total_zcut") for k in _sd_gens()), 1.0),
        ("sidis: nu_mu pions across the generators, high", "2925",
         max(_sd("nu", "nu_mu", k, "sidis_e", "pi_total_zcut") for k in _sd_gens()), 1.0),
        ("sidis: nu charge ratio above z=0.1", "1.30",
         _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pip_total_zcut")
         / _sd("nu", "nu_mu", "powheg_nu", "sidis_e", "pim_total_zcut"), 0.005),
        ("sidis: charge ratio across the generators, low", "1.25",
         min(_sd("nu", "nu_mu", k, "sidis_e", "pip_total_zcut")
             / _sd("nu", "nu_mu", k, "sidis_e", "pim_total_zcut") for k in _sd_gens()), 0.005),
        ("sidis: charge ratio across the generators, high", "1.31",
         max(_sd("nu", "nu_mu", k, "sidis_e", "pip_total_zcut")
             / _sd("nu", "nu_mu", k, "sidis_e", "pim_total_zcut") for k in _sd_gens()), 0.005),
        # --- the Tier E efficiency factors (user, 2026-09-21) ------------
        # The bridge between a calculation in the region alone and a FASER
        # yield: every number the SIDIS tab's blurb quotes for it.
        ("eff: eps_tier at z -> 0, nu_mu", "0.40",
         _eff_z("nu_mu", "pi", "eps_tier", 0.02), 0.005),
        ("eff: eps_tier at z = 0.8, nu_mu", "0.18",
         _eff_z("nu_mu", "pi", "eps_tier", 0.78), 0.005),
        ("eff: emulsion acceptance above z = 0.05", "1",
         _eff_z("nu_mu", "pi", "eps_acc", 0.30), 1e-6),
        ("eff: emulsion acceptance, lowest z bin", "0.88",
         _eff_z("nu_mu", "pi", "eps_acc", 0.02), 0.005),
        ("eff: Tier E event efficiency at 300 GeV", "0.03",
         _eff_at_energy(300), 0.005),
        ("eff: Tier E event efficiency at 1 TeV", "0.48",
         _eff_at_energy(1000), 0.005),
        ("eff: Tier E event efficiency at 2 TeV", "0.57",
         _eff_at_energy(2000), 0.005),
        ("eff: flux-averaged Tier E efficiency, nu_mu", "0.30",
         _eff_event("nu_mu"), 0.005),
        ("eff: flux-averaged Tier E efficiency, nu_e", "0.40",
         _eff_event("nu_e"), 0.005),
        ("eff: flux-averaged Tier E efficiency, nu_mu, POWHEG-V2", "0.31",
         _eff_event("nu_mu", "powheg_nu"), 0.005),
        ("eff: flux-averaged Tier E efficiency, nu_e, POWHEG-V2", "0.42",
         _eff_event("nu_e", "powheg_nu"), 0.005),
        ("eff: mean flux energy, nu_e [GeV]", "733", _flux_mean("12"), 1.0),
        ("eff: mean flux energy, nu_mu [GeV]", "312", _flux_mean("14"), 1.0),
        # the statistical reach, and the clustering factor that sets it
        ("sidis: clustering factor", "1.15", _sd_cluster(), 0.005),
        ("sidis: 10% reach, nu_mu", "0.40", _sd_reach("nu", "nu_mu"), 1e-9),
        ("sidis: 10% reach, nu_e", "0.25", _sd_reach("nu", "nu_e"), 1e-9),
        # how often a DIS event has a pion at all (user, 2026-09-07)
        ("sidis: nu events with a charged pion", "98.6",
         100 * _sd_pres("nu", "nu_mu", "powheg_nu", "sidis_w", "any_z"), 0.05),
        ("sidis: nu events with a charged pion, GENIE", "98.9",
         100 * _sd_pres("nu", "nu_mu", "genie_lo", "sidis_w", "any_z"), 0.05),
        ("sidis: nu events with a pion above the z cut", "78.5",
         100 * _sd_pres("nu", "nu_mu", "powheg_nu", "sidis_w", "zcut_emulsion"), 0.05),
        ("sidis: nu events with a pion above the z cut, GENIE", "75.3",
         100 * _sd_pres("nu", "nu_mu", "genie_lo", "sidis_w", "zcut_emulsion"), 0.05),
        # --- the FASER dimuon cut flow (2026-09-07) ---------------------------
        ("dimuon: POWHEG-V2 charm fraction, Q2>4", "13.3",
         _dmc_row("powheg_q2", "14", "charm", "cum_eff_pct"), 0.05),
        ("dimuon: GENIE charm fraction, Q2>4", "8.9",
         _dmc_row("genie_q2", "14", "charm", "cum_eff_pct"), 0.05),
        ("dimuon: POWHEG-V2 charm->mu step", "9.31",
         _dmc_row("powheg_q2", "14", "charm_mu", "step_eff_pct"), 0.01),
        ("dimuon: GENIE p>100 in 25 mrad, nu_mu", "7.9",
         _dmc_row("genie_q2", "14", "p100", "events"), 0.05),
        ("dimuon: POWHEG-V2 p>100 in 25 mrad, nu_mu", "13.1",
         _dmc_row("powheg_q2", "14", "p100", "events"), 0.05),
        ("dimuon: GENIE p>100 in 25 mrad, nu+nubar", "10.5",
         _dmc_row("genie_q2", "sum", "p100", "events"), 0.05),
        ("dimuon: POWHEG-V2 p>100 in 25 mrad, nu+nubar", "17.3",
         _dmc_row("powheg_q2", "sum", "p100", "events"), 0.05),
        ("dimuon: POWHEG-V2 over GENIE, p>100 in 25 mrad", "1.64",
         _dmc_row("powheg_q2", "sum", "p100", "events")
         / _dmc_row("genie_q2", "sum", "p100", "events"), 0.005),
        ("dimuon: ladder over flux mode, worst row", "1.02",
         max(abs(v) for v in _dmc()["flux_mode"]["genie_total"]["ladder_over_flux"].values()), 0.005),
    ]


# Claims on prose that compares with material NOT distributed with the public
# release (data/faser_internal/): a FASER-internal dimuon estimate and the NNLO
# SIDIS curves of arXiv:2504.05376 as provided by its authors.  Checked only
# where that material, and the notes quoting it, exist.
PRIVATE_NOTES = [f"{BASE}/data/faser_internal/dimuon_report_note.html",
                 f"{BASE}/data/faser_internal/sidis_nnlo_report_note.html"]


def private_claims():
    if not all(os.path.exists(p) for p in PRIVATE_NOTES):
        return []
    return [
        # the comparison, at their beam energy, on the proton sample
        ("sidis: published NNLO integral above z=0.1", "0.883",
         _sid()["comparison"]["reference_total_zcut"], 5e-4),
        ("sidis: comparison beam energy", "300",
         _sid()["comparison"]["energy_gev"], 1e-9),
        ("sidis: every generator compared at 300 GeV", "5",
         sum(1 for g in _sid()["comparison"]["generators"].values() if g["matched_energy"]), 1e-9),
        ("sidis: POWHEG-V2 integral at their energy", "0.860",
         _sd_cmp("powheg_nu", "total_zcut"), 5e-4),
        ("sidis: POWHEG-V2 over NNLO integral", "0.97",
         _sd_cmp("powheg_nu", "ratio_total_zcut"), 0.005),
        ("sidis: Herwig integral", "0.888",
         _sd_cmp("herwig_nlo_full", "total_zcut"), 5e-4),
        ("sidis: Sherpa integral", "0.913",
         _sd_cmp("sherpa_nlo", "total_zcut"), 5e-4),
        ("sidis: GENIE integral", "0.804",
         _sd_cmp("genie_lo", "total_zcut"), 5e-4),
        ("sidis: GENIE HEDIS integral", "0.787",
         _sd_cmp("genie", "total_zcut"), 5e-4),
        ("sidis: all five within, per cent", "11",
         100 * max(abs(g["ratio_total_zcut"] - 1.0)
                   for g in _sid()["comparison"]["generators"].values()), 1.0),
        ("sidis: POWHEG-V2 mean ratio, z>0.2", "1.33",
         _sd_cmp("powheg_nu", "mean_ratio_z_gt_0p2"), 0.005),
        ("sidis: GENIE mean ratio, z>0.2", "1.84",
         _sd_cmp("genie_lo", "mean_ratio_z_gt_0p2"), 0.005),
        ("sidis: GENIE HEDIS mean ratio, z>0.2", "1.90",
         _sd_cmp("genie", "mean_ratio_z_gt_0p2"), 0.005),
        ("sidis: Herwig mean ratio, z>0.2", "1.01",
         _sd_cmp("herwig_nlo_full", "mean_ratio_z_gt_0p2"), 0.005),
        ("sidis: Sherpa mean ratio, z>0.2", "0.67",
         _sd_cmp("sherpa_nlo", "mean_ratio_z_gt_0p2"), 0.005),
        ("sidis: their x cut moves the multiplicity by at most", "1.8",
         100 * max(abs(g["total_zcut_no_x"] / g["total_zcut"] - 1.0)
                   for g in _sid()["comparison"]["generators"].values()), 0.05),
        ("dimuon: GENIE nu_mu CC at 300 fb^-1", "8043",
         _dmc_row("genie_total", "14", "cc", "events"), 1.0),
        ("dimuon: slide row 1 over ours", "0.69",
         _dmc()["vs_slide"]["genie_total_14"]["cc"]["slide_events"]
         / _dmc_row("genie_total", "14", "cc", "events"), 0.005),
        ("dimuon: GENIE charm fraction", "8.74",
         _dmc_row("genie_total", "14", "charm", "cum_eff_pct"), 0.01),
        ("dimuon: GENIE charm->mu step", "9.04",
         _dmc_row("genie_total", "14", "charm_mu", "step_eff_pct"), 0.01),
        ("dimuon: GENIE cumulative p>20 in 25 mrad", "0.221",
         _dmc_row("genie_total", "14", "p20", "cum_eff_pct"), 0.001),
        ("dimuon: GENIE cumulative p>50 in 25 mrad", "0.160",
         _dmc_row("genie_total", "14", "p50", "cum_eff_pct"), 0.001),
        ("dimuon: GENIE cumulative p>100 in 25 mrad", "0.100",
         _dmc_row("genie_total", "14", "p100", "cum_eff_pct"), 0.001),
        ("dimuon: p>20 ours over slide", "0.44",
         _dmc()["vs_slide"]["genie_total_14"]["p20"]["cum_ratio"], 0.005),
        ("dimuon: p>50 ours over slide", "0.58",
         _dmc()["vs_slide"]["genie_total_14"]["p50"]["cum_ratio"], 0.005),
        ("dimuon: p>100 ours over slide", "0.83",
         _dmc()["vs_slide"]["genie_total_14"]["p100"]["cum_ratio"], 0.005),
        ("dimuon: GENIE cumulative p>20, no angle", "0.448",
         _dmc_row("genie_total", "14", "p20_noth", "cum_eff_pct"), 0.001),
        ("dimuon: GENIE cumulative p>50, no angle", "0.264",
         _dmc_row("genie_total", "14", "p50_noth", "cum_eff_pct"), 0.001),
        ("dimuon: GENIE cumulative p>100, no angle", "0.142",
         _dmc_row("genie_total", "14", "p100_noth", "cum_eff_pct"), 0.001),
        ("dimuon: p>100 no angle, ours over slide", "1.18",
         _dmc_row("genie_total", "14", "p100_noth", "cum_eff_pct") / 0.120, 0.01),
    ]


def main():
    verbose = "--verbose" in sys.argv
    with open(REPORT_SRC) as f:
        prose = f.read()
    for p in PRIVATE_NOTES:
        if os.path.exists(p):
            with open(p) as f:
                prose += f.read()
    bad, orphan, ok = [], [], 0
    for what, quoted, computed, tol in claims() + private_claims():
        if computed is None:
            bad.append(f"{what}: cannot recompute (result missing)")
            continue
        # the prose must actually still contain the string being checked
        if quoted not in prose:
            orphan.append(f"{what}: prose no longer contains {quoted!r}")
            continue
        # a quoted count may carry the thin-space separator the page uses
        numeric = quoted.replace("&thinsp;", "").replace(",", "")
        if abs(float(numeric) - computed) > tol:
            bad.append(f"{what}: prose says {quoted}, results give "
                       f"{computed:.6g}")
        else:
            ok += 1
            if verbose:
                print(f"  ok  {what:44s} {quoted} == {computed:.6g}")
    print(f"\n{ok} claim(s) verified against results/")
    for o in orphan:
        print(f"ORPHAN  {o}")
    for b in bad:
        print(f"WRONG   {b}")
    if bad:
        print(f"\n{len(bad)} claim(s) in the report prose disagree with the "
              f"results. Fix the prose, or the claim table here if the "
              f"sentence was rewritten.")
        return 1
    if orphan:
        print(f"\n{len(orphan)} claim(s) guard prose that no longer exists; "
              f"drop them from the table.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
