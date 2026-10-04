#!/usr/bin/env python3
"""Impact of non-DIS processes, layer 2: the differential distributions.

WHAT THIS ANSWERS.  Layer 1 (genie_nondis.py) split GENIE's INTEGRATED
cross-section channel by channel and found the explicit non-DIS channels to
be per-mille at FASER energies.  This layer asks WHERE in the kinematic
plane those channels sit, and how much of the benchmark's own fiducial
region they contaminate -- and, the larger question, how much of what GENIE
labels DIS lies outside the benchmark's region (the shallow-inelastic
corner: low Q2, low W).  Six observables, Q2, x, y, W, E' and theta of the
scattered lepton, at 400 GeV and 1 TeV, fully inclusive: no selection and no
acceptance cut at generation (user, 2026-09-04).

INPUTS.  Two GENIE runs per (current, energy) from genie/run_genie_nondis.sh,
as gntpc summary ntuples (gst), on a free proton:

  ndfull_<cur>_job[_tag]_N   the stock default list (CC / EM): every channel
                             at its own rate.  Provides the DIS shapes and
                             the CLOSURE: its per-channel event fractions
                             must reproduce the spline fractions.
  ndonly_<cur>_job[_tag]_N   the same list minus DIS: the non-DIS shapes with
                             full statistics.  Each channel is normalised to
                             its own spline cross-section (layer 1's JSON),
                             which is what gevgen does internally anyway.

The kinematics are the lepton-based ones the benchmark analysis also uses
(gst branches Q2, x, y, W, El, cthl -- computed from the outgoing lepton,
not the generator's "selected" Q2s/xs/ys/Ws, which coincide for QEL/RES and
differ only by the hadronic-system bookkeeping for DIS).

>>> THE MUON SIDE IS Q2 > 4 GeV2 ONLY, by this build's electromagnetic
    floor (patches/genie-em-q2-floor-and-pythia8-teardown.diff): "fully
    inclusive" for muons means every EM channel above 4 GeV2, so the muon
    fiducial fractions are y-window fractions of a Q2 > 4 sample. <<<

OUTPUT, per (current, energy):
  results[_nu]/genie_nondis_diff[_400GeV].json   the histograms (per channel,
        dsigma/dX in pb per unit), the region fractions, the closure
  results[_nu]/cmp_genie_nondis_diff_<cur>[_400GeV].png   six observables,
        each with the channel share per bin beneath

Usage: analysis/genie_nondis_diff.py [--no-plot] [mu|nu] [100|400|1000]

THE TARGET IS TUNGSTEN (user, 2026-09-14: "regenerate GENIE's non-DIS sample
on p+n, so corresponding to a tungsten nucleus").  The inputs
are v2ndfull/v2ndonly_<cur>_<p|n>_job (genie/run_genie_nondis.sh ... <p|n>),
each nucleon normalised to its OWN layer-1 spline cross-sections, and the
per-nucleon tungsten result is (74 p + 110 n)/184: cross-sections and
histograms add with those weights, region fractions are re-formed from the
combined channel cross-sections.  The "fiducial" region is the benchmark
region, Q2 > 4 GeV2 and W > 3 GeV with no y window, and the files
are genie_nondis_diff_<p|n|W>[_tag].json.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import beams  # noqa: E402
import selection  # noqa: E402

GENIE = f"{BASE}/genie"
GROUPS = ["DIS", "DIS charm", "RES", "QEL", "DFR"]
NONDIS = ["RES", "QEL", "DFR"]
COLOURS = {"DIS": "#2b6cb0", "DIS charm": "#c2317b", "RES": "#3fa66a",
           "QEL": "#eb6834", "DFR": "#8e5bd0"}
TOTAL_COLOUR = "#1f1f1f"
ENERGIES = (400.0, 1000.0)
# the benchmark region, from the one place that defines it
Q2_FLOOR = selection.Q2_FLOOR
Y_LO, Y_HI = 0.2, 0.9
W_DIS = 2.0          # the conventional resonance/DIS boundary; NOT a benchmark cut
M_P = beams.M_P


def _resdir(cur):
    return f"{BASE}/results_nu" if cur == "nu" else f"{BASE}/results"


def _tagged(stem, e):
    return beams.at_energy(stem, e)


W_MIN = 3.0          # the benchmark region's W floor (selection q4w3)
Z_W, A_W = 74, 184


def _jobs(kind, cur, e, t=None):
    base = beams.at_energy(f"v2{kind}_{cur}_{t}_job" if t
                           else f"{kind}_{cur}_job", e)
    # the anchor's directories carry no tag, and a tagged glob must not fold
    # the other energies in (the recurring *job_* trap): match _<digits> only
    out = []
    for d in sorted(glob.glob(f"{GENIE}/{base}_[0-9]*")):
        tail = d[len(f"{GENIE}/{base}_"):]
        if tail.isdigit() and os.path.exists(f"{d}/events.gst.root"):
            out.append(f"{d}/events.gst.root")
    return out


def read(files, ebeam):
    """Concatenated arrays from the gst trees: kinematics and channel."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    cols = ["Q2", "x", "y", "W", "El", "cthl", "Ev", "qel", "res", "dis",
            "dfr", "charm", "coh", "mec"]
    parts = {c: [] for c in cols}
    for fn in files:
        df = ROOT.RDataFrame("gst", fn)
        arr = df.AsNumpy(cols)
        for c in cols:
            parts[c].append(np.asarray(arr[c]))
    a = {c: np.concatenate(parts[c]) for c in cols}
    if not np.allclose(a["Ev"], ebeam):
        raise SystemExit(f"beam energy in the ntuples is not {ebeam} GeV: "
                         f"{a['Ev'].min()}..{a['Ev'].max()}")
    if a["coh"].any() or a["mec"].any():
        raise SystemExit("coherent or MEC events in a free-nucleon sample")
    g = np.full(len(a["Q2"]), "", dtype=object)
    g[a["qel"].astype(bool)] = "QEL"
    g[a["res"].astype(bool)] = "RES"
    g[a["dfr"].astype(bool)] = "DFR"
    d = a["dis"].astype(bool)
    g[d & ~a["charm"].astype(bool)] = "DIS"
    g[d & a["charm"].astype(bool)] = "DIS charm"
    if (g == "").any():
        raise SystemExit(f"{(g == '').sum()} events with no channel flag")
    a["group"] = g
    a["theta"] = np.arccos(np.clip(a["cthl"], -1.0, 1.0))
    return a


def spline_sigmas(cur, e, t="p"):
    """{group: sigma_pb} on nucleon t at this energy, from layer 1."""
    with open(f"{_resdir(cur)}/genie_nondis.json") as f:
        d = json.load(f)
    b = d["benchmark"][t][f"{e:g}"]
    return {g: b["sigma_pb"][g] for g in GROUPS}, b["total_pb"]


def observables(cur, e):
    """(name, label, edges, log?) for the six distributions."""
    s = M_P ** 2 + 2.0 * M_P * e
    # the muon sample starts at Q2 = 4 (this build's EM floor), which puts
    # x above 4 / 2ME and the lepton angle above ~1 mrad: no empty decades
    q2lo = Q2_FLOOR if cur == "mu" else 1e-2
    xlo = -3 if cur == "mu" else -5
    thlo = -3.5 if cur == "mu" else -5
    return [
        ("Q2", r"$Q^2$ [GeV$^2$]", np.logspace(np.log10(q2lo), np.log10(s), 36), True),
        ("x", r"$x_{\rm Bj}$", np.logspace(xlo, 0, 31), True),
        ("y", r"$y$", np.linspace(0.0, 1.0, 26), False),
        ("W", r"$W$ [GeV]", np.logspace(np.log10(0.9), np.log10(np.sqrt(s)), 31), True),
        ("El", r"$E'_\ell$ [GeV]", np.linspace(0.0, e, 26), False),
        ("theta", r"$\theta_\ell$ [rad]", np.logspace(thlo, 0, 31), True),
    ]


def region_observables(cur, e):
    """Q2, W and x INSIDE the benchmark region, Q2 > 4 GeV2 AND W > 3 GeV applied
    together to every one of them (user, 2026-09-21: "Both cuts should be
    applied at once ... I want to quantify impact of non-DIS processes for
    the fiducial cuts of the benchmark").  Q2 and W are booked from their
    cut upwards so that no bin straddles it.  Keys "<name>_region", read by
    ppA1; hist() applies the region to any key with that suffix."""
    s = M_P ** 2 + 2.0 * M_P * e
    xlo = np.log10(Q2_FLOOR / (s - M_P ** 2))    # x >= Q2 / (2 M E) at y = 1
    return [
        ("Q2_region", r"$Q^2$ [GeV$^2$]",
         np.logspace(np.log10(Q2_FLOOR), np.log10(s), 31), True),
        ("W_region", r"$W$ [GeV]",
         np.logspace(np.log10(W_MIN), np.log10(np.sqrt(s)), 26), True),
        ("x_region", r"$x_{\rm Bj}$", np.logspace(np.floor(xlo), 0, 31), True),
    ]


def hist(a, sel, w, edges, name):
    """dsigma/dX in pb per unit of X, and its statistical error."""
    if name.endswith("_region"):
        sel = sel & (a["Q2"] > Q2_FLOOR) & (a["W"] > W_MIN)
    v = a[name.replace("_region", "")][sel]
    n, _ = np.histogram(v, bins=edges)
    widths = np.diff(edges)
    return n * w / widths, np.sqrt(n) * w / widths, n


def region_fractions(a, sel, tungsten=False):
    """Fractions of the events in `sel` landing in each region.  On tungsten the
    "fiducial" region is Q2 > 4 and W > 3 GeV, no y window."""
    q2 = a["Q2"][sel]
    y = a["y"][sel]
    W = a["W"][sel]
    n = max(sel.sum(), 1)
    fid = ((q2 > Q2_FLOOR) & (W > W_MIN) if tungsten
           else (q2 > Q2_FLOOR) & (y > Y_LO) & (y < Y_HI))
    return {
        "n": int(sel.sum()),
        "q2_above_floor": float((q2 > Q2_FLOOR).sum() / n),
        "fiducial": float(fid.sum() / n),
        "y_window": float(((y > Y_LO) & (y < Y_HI)).sum() / n),
        "w_above_2": float((W > W_DIS).sum() / n),
        "q2_above_floor_and_w_above_2": float(((q2 > Q2_FLOOR) & (W > W_DIS)).sum() / n),
        "mean_q2": float(q2.mean()) if n else None,
        "mean_w": float(W.mean()) if n else None,
        "mean_y": float(y.mean()) if n else None,
    }


def delivery(full, only, closure):
    """What GENIE generated but did not deliver: events it found unphysical
    and redrew from the whole list (genie/nondis_log_digest.py).

    The resonance kinematics generator gives up on a (W,Q2) pair after 1001
    tries, and the replacement event is drawn afresh from every channel, so
    the delivered resonance rate falls below the spline one by the failure
    rate.  Both the log counts and the closure deficit are recorded; they
    are two measurements of the same number.
    """
    def counts(files):
        n = {}
        for fn in files:
            p = os.path.join(os.path.dirname(fn), "rejections.json")
            if not os.path.exists(p):
                return None
            with open(p) as f:
                for k, v in json.load(f).items():
                    if isinstance(v, (int, float)):
                        n[k] = n.get(k, 0) + v
        return n
    cf, co = counts(full), counts(only)
    d = {"full_run": cf, "nondis_run": co}
    if cf:
        n_res = closure["RES"]["events"]
        n_qel = closure["QEL"]["events"]
        d["res_failure_rate_full"] = cf["res_kinematics_failed"] / (
            n_res + cf["res_kinematics_failed"]) if n_res else None
        d["qel_failure_rate_full"] = cf["qel_kinematics_failed"] / (
            n_qel + cf["qel_kinematics_failed"]) if n_qel else None
    d["res_delivered_over_spline"] = (closure["RES"]["fraction_events"]
                                      / closure["RES"]["fraction_spline"])
    d["qel_delivered_over_spline"] = (closure["QEL"]["fraction_events"]
                                      / closure["QEL"]["fraction_spline"]
                                      if closure["QEL"]["fraction_spline"] else None)
    return d


def analyse(cur, e, t=None):
    tungsten = t is not None
    full = _jobs("ndfull", cur, e, t)
    only = _jobs("ndonly", cur, e, t)
    if not full or not only:
        raise SystemExit(f"no ntuples for {cur} {e:g} GeV{' on ' + t if tungsten else ''} "
                         f"-- run genie/run_genie_nondis.sh {cur} {e:g}"
                         f"{' 1000000 200000 8 ' + t if tungsten else ''}")
    F = read(full, e)
    O = read(only, e)
    sig, sig_tot = spline_sigmas(cur, e, t or "p")
    nF = len(F["Q2"])
    wF = sig_tot / nF                      # pb per event, inclusive run
    # per-channel weights in the non-DIS run: each channel to its own spline
    # cross-section, which is the normalisation gevgen would have used
    nO = {g: int((O["group"] == g).sum()) for g in NONDIS}
    wO = {g: (sig[g] / nO[g] if nO[g] else 0.0) for g in NONDIS}
    # THE CLOSURE: the inclusive run's channel fractions against the splines
    closure = {}
    for g in GROUPS:
        n_g = int((F["group"] == g).sum())
        closure[g] = {"events": n_g, "fraction_events": n_g / nF,
                      "fraction_spline": sig[g] / sig_tot,
                      "fraction_spline_err_events": (np.sqrt(n_g) / nF if n_g else None)}
    out = {"current": cur, "energy_gev": e, "target": t or "p",
           "region": ("Q2 > 4 GeV2, W > 3 GeV, no y cut" if tungsten
                      else "proton: Q2 > 4 GeV2, 0.2 < y < 0.9"),
           "tune": "G18_02a_00_000",
           "list_full": "CC" if cur == "nu" else "EM",
           "list_nondis": "CCNONDIS" if cur == "nu" else "EMNONDIS",
           "em_q2_floor_gev2": Q2_FLOOR if cur == "mu" else None,
           "n_full": nF, "n_nondis_run": {g: nO[g] for g in NONDIS},
           "sigma_pb": dict(sig, total=sig_tot),
           "q2_floor": Q2_FLOOR, "w_dis": W_DIS,
           **({"w_min": W_MIN} if tungsten else {"y_window": [Y_LO, Y_HI]}),
           "closure": closure, "regions": {}, "hist": {},
           "inputs": {"full": [os.path.relpath(f, BASE) for f in full],
                      "nondis": [os.path.relpath(f, BASE) for f in only]}}
    # region fractions per channel (DIS from the inclusive run, the rest
    # from the dedicated one), then the channel shares of each region
    src = {g: (F, F["group"] == g) for g in ("DIS", "DIS charm")}
    src.update({g: (O, O["group"] == g) for g in NONDIS})
    for g in GROUPS:
        a, sel = src[g]
        out["regions"][g] = region_fractions(a, sel, tungsten)
    # the DIS channel outside the benchmark region: the shallow-inelastic
    # corner, the piece layer 1 could not see
    dis_all = sig["DIS"] + sig["DIS charm"]
    dis_out = sum(sig[g] * (1.0 - out["regions"][g]["q2_above_floor"])
                  for g in ("DIS", "DIS charm"))
    dis_lowW = sum(sig[g] * (1.0 - out["regions"][g]["w_above_2"])
                   for g in ("DIS", "DIS charm"))
    out["delivery"] = delivery(full, only, closure)
    _shares(out)
    out["dis_channel"] = {
        "sigma_pb": dis_all,
        "fraction_below_q2_floor": dis_out / dis_all,
        "fraction_below_w2": dis_lowW / dis_all,
        "fraction_outside_fiducial": 1.0 - sum(
            sig[g] * out["regions"][g]["fiducial"] for g in ("DIS", "DIS charm")) / dis_all}
    # histograms; the per-nucleon run also books Q2, W and x INSIDE the region, both
    # cuts at once, for the appendix figure (user, 2026-09-21)
    for name, label, edges, islog in observables(cur, e) + (
            region_observables(cur, e) if t is not None else []):
        h = {"label": label, "edges": [float(v) for v in edges], "log": islog}
        for g in GROUPS:
            a, sel = src[g]
            w = wF if g.startswith("DIS") else wO[g]
            d, err, n = hist(a, sel, w, edges, name)
            h[g] = {"dsig": [float(v) for v in d], "err": [float(v) for v in err],
                    "n": [int(v) for v in n]}
        out["hist"][name] = h
    fn = f"{_resdir(cur)}/{_tagged('genie_nondis_diff' + (f'_{t}' if tungsten else '_proton'), e)}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {fn}")
    return out


def _shares(out):
    """share_of_<region> and sigma_<region>_pb from the channel cross-sections
    and the per-channel region fractions (also used for the W combination)."""
    sig = out["sigma_pb"]
    for reg in ("fiducial", "q2_above_floor", "w_above_2"):
        parts = {g: sig[g] * (out["regions"][g][reg] or 0.0) for g in GROUPS}
        tot = sum(parts.values())
        out["regions"][f"share_of_{reg}"] = {
            g: (parts[g] / tot if tot else 0.0) for g in GROUPS}
        out["regions"][f"share_of_{reg}"]["non-DIS"] = (
            sum(parts[g] for g in NONDIS) / tot if tot else 0.0)
        out["regions"][f"sigma_{reg}_pb"] = tot


def combine_w(p, n, cur, e):
    """Tungsten per nucleon, (74 p + 110 n)/184."""
    fz, fn_ = Z_W / A_W, (A_W - Z_W) / A_W
    for k in ("current", "energy_gev", "region", "tune"):
        if p[k] != n[k]:
            raise SystemExit(f"p and n disagree on {k}: {p[k]} vs {n[k]}")
    if (p["target"], n["target"]) != ("p", "n"):
        raise SystemExit("combine_w needs a p and an n result")
    out = {k: p[k] for k in ("current", "energy_gev", "region", "tune",
                             "list_full", "list_nondis", "em_q2_floor_gev2",
                             "q2_floor", "w_dis", "w_min")}
    out.update(target="W", per_nucleon=True, Z=Z_W, A=A_W,
               combination="(74 p + 110 n)/184; analysis/genie_nondis_diff.py",
               n_full={"p": p["n_full"], "n": n["n_full"]},
               closure={"p": p["closure"], "n": n["closure"]},
               delivery={"p": p["delivery"], "n": n["delivery"]},
               inputs={"p": p["inputs"], "n": n["inputs"]})
    sig = {g: fz * p["sigma_pb"][g] + fn_ * n["sigma_pb"][g] for g in GROUPS}
    sig["total"] = fz * p["sigma_pb"]["total"] + fn_ * n["sigma_pb"]["total"]
    out["sigma_pb"] = sig
    out["regions"] = {}
    for g in GROUPS:
        rp, rn = p["regions"][g], n["regions"][g]
        sp, sn = fz * p["sigma_pb"][g], fn_ * n["sigma_pb"][g]
        r = {"n": rp["n"] + rn["n"]}
        for k in rp:
            if k == "n":
                continue
            a, b = rp[k], rn[k]
            if not (sp + sn):
                # a channel with no cross-section on either nucleon (muon DIS
                # charm: GENIE's EM DIS has none) -- its fractions weigh zero
                r[k] = 0.0 if k not in ("mean_q2", "mean_w", "mean_y") else None
            else:
                r[k] = (None if a is None or b is None
                        else (sp * a + sn * b) / (sp + sn))
        out["regions"][g] = r
    _shares(out)
    dis_all = sig["DIS"] + sig["DIS charm"]
    out["dis_channel"] = {
        "sigma_pb": dis_all,
        "fraction_below_q2_floor": sum(sig[g] * (1 - out["regions"][g]["q2_above_floor"])
                                       for g in ("DIS", "DIS charm")) / dis_all,
        "fraction_below_w2": sum(sig[g] * (1 - out["regions"][g]["w_above_2"])
                                 for g in ("DIS", "DIS charm")) / dis_all,
        "fraction_outside_fiducial": 1.0 - sum(
            sig[g] * out["regions"][g]["fiducial"] for g in ("DIS", "DIS charm")) / dis_all}
    out["hist"] = {}
    for name, hp in p["hist"].items():
        hn = n["hist"][name]
        if hp["edges"] != hn["edges"]:
            raise SystemExit(f"{name}: p and n histograms have different edges")
        h = {"label": hp["label"], "edges": hp["edges"], "log": hp["log"]}
        for g in GROUPS:
            h[g] = {"dsig": [fz * a + fn_ * b for a, b in zip(hp[g]["dsig"], hn[g]["dsig"])],
                    "err": [float(np.hypot(fz * a, fn_ * b))
                            for a, b in zip(hp[g]["err"], hn[g]["err"])],
                    "n": [a + b for a, b in zip(hp[g]["n"], hn[g]["n"])]}
        out["hist"][name] = h
    fn = f"{_resdir(cur)}/{_tagged('genie_nondis_diff_W', e)}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {fn}")
    return out


def figure(out, cur, e, tungsten=False):
    import plotstyle
    plotstyle.apply()
    from plotstyle import tex
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.lines import Line2D

    obs = observables(cur, e)
    fig = plt.figure(figsize=(15.0, 10.8))
    # two rows of three observables; each observable is a FUSED pair (the
    # distribution and the share beneath, no gap), with room between pairs
    outer = gridspec.GridSpec(2, 3, figure=fig, hspace=0.28, wspace=0.26,
                              left=0.06, right=0.99, top=0.90, bottom=0.06)
    handles = None
    for k, (name, label, edges, islog) in enumerate(obs):
        r, c = divmod(k, 3)
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[r, c], height_ratios=[2.0, 1.0], hspace=0.0)
        ax = fig.add_subplot(inner[0])
        ar = fig.add_subplot(inner[1], sharex=ax)
        h = out["hist"][name]
        ed = np.array(h["edges"])
        tot = sum(np.array(h[g]["dsig"]) for g in GROUPS)
        ax.step(ed, np.r_[tot, tot[-1]], where="post", color=TOTAL_COLOUR,
                lw=1.8, label=tex("total"))
        for g in GROUPS:
            d = np.array(h[g]["dsig"])
            if not d.any():
                continue
            ax.step(ed, np.r_[d, d[-1]], where="post", color=COLOURS[g],
                    lw=1.4, label=tex(g))
            if g in NONDIS:
                share = np.where((tot > 0) & (np.array(h[g]["n"]) >= 5),
                                 100.0 * d / np.where(tot > 0, tot, 1), np.nan)
                ar.step(ed, np.r_[share, share[-1]], where="post",
                        color=COLOURS[g], lw=1.2)
        nd = sum(np.array(h[g]["dsig"]) for g in NONDIS)
        # a share built on fewer than five non-DIS events is noise, not a
        # measurement: the lone high-Q2 resonance event in an empty bin
        # would otherwise draw as a 3% spike
        nnd = sum(np.array(h[g]["n"]) for g in NONDIS)
        share = np.where((tot > 0) & (nnd >= 5),
                         100.0 * nd / np.where(tot > 0, tot, 1), np.nan)
        ar.step(ed, np.r_[share, share[-1]], where="post", color=TOTAL_COLOUR,
                lw=1.8)
        # the benchmark's cuts where the observable carries them (dashed);
        # the conventional W = 2 GeV boundary, which is NOT a cut (dotted)
        cuts = ({"Q2": [(Q2_FLOOR, "--")], "W": [(W_DIS, ":"), (W_MIN, "--")]}
                if tungsten else
                {"Q2": [(Q2_FLOOR, "--")], "y": [(Y_LO, "--"), (Y_HI, "--")],
                 "W": [(W_DIS, ":")]}).get(name, [])
        for v, ls_ in cuts:
            for a in (ax, ar):
                a.axvline(v, color="#666666", lw=1.0, ls=ls_)
        if islog:
            ax.set_xscale("log")
        ax.set_yscale("log")
        ar.set_yscale("log")
        ar.set_ylim(1e-2, 150.0)
        ar.axhline(100.0, color="#bbbbbb", lw=0.6)
        ax.set_xlim(ed[0], ed[-1])
        pos = tot[tot > 0]
        if len(pos):
            ax.set_ylim(pos.min() * 0.3, pos.max() * 3.0)
        ax.set_ylabel(tex(r"d$\sigma$/d" + label.split(" [")[0] + " [pb/unit]"))
        ar.set_ylabel(tex("non-DIS share [%]"))
        ar.set_xlabel(tex(label))
        plotstyle.ticks(ax, labelbottom=False)
        plotstyle.ticks(ar)
        if handles is None:
            handles = ax.get_legend_handles_labels()
    # one legend for the figure, above the panels, so it covers no data
    hs, ls = handles
    hs = hs + [Line2D([], [], color="#666666", lw=1.0, ls="--"),
               Line2D([], [], color="#666666", lw=1.0, ls=":")]
    ls = ls + [tex("benchmark cut"), tex("W = 2 GeV (not a cut)")]
    fig.legend(hs, ls, loc="upper center", ncol=len(hs), frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, 0.955))
    beam = {"nu": r"$\nu_\mu$ CC", "mu": r"$\mu^-$ EM, $Q^2 > 4$ GeV$^2$ only"}[cur]
    etxt = f"{e:g} GeV" if e < 1000 else f"{e/1000:g} TeV"
    tgt = "tungsten (per nucleon)" if tungsten else "a proton"
    fig.suptitle(tex(f"GENIE (GRV98LO), {beam} on {tgt} at {etxt}, fully "
                     "inclusive: channel by channel"), y=0.985)
    fn = f"{_resdir(cur)}/{_tagged('cmp_genie_nondis_diff_' + cur + ('_W' if tungsten else '_proton'), e)}.png"
    fig.savefig(fn, dpi=150)
    plt.close(fig)
    print(f"wrote {fn}")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    plot = "--no-plot" not in sys.argv
    tungsten = True   # tungsten per nucleon: the proton-only mode's samples are gone
    curs = [a for a in args if a in ("mu", "nu")] or ["nu", "mu"]
    es = [float(a) for a in args if a.replace(".", "").isdigit()] or list(ENERGIES)
    for cur in curs:
        for e in es:
            if tungsten:
                out = combine_w(analyse(cur, e, "p"), analyse(cur, e, "n"), cur, e)
            else:
                out = analyse(cur, e)
            sh = out["regions"]["share_of_fiducial"]
            dc = out["dis_channel"]
            print(f"  {cur} {e:g} GeV: non-DIS share of the fiducial rate "
                  f"{100*sh['non-DIS']:.3f}%  (RES {100*sh['RES']:.3f}, QEL "
                  f"{100*sh['QEL']:.3f}, DFR {100*sh['DFR']:.3f}); DIS channel "
                  f"below Q2 floor {100*dc['fraction_below_q2_floor']:.2f}%, "
                  f"below W = 2 {100*dc['fraction_below_w2']:.2f}%, outside "
                  f"fiducial {100*dc['fraction_outside_fiducial']:.2f}%")
            for t_, cl in (out["closure"].items() if tungsten else [("p", out["closure"])]):
                for g, c in cl.items():
                    print(f"     closure {t_} {g:10s} events {100*c['fraction_events']:.3f}%"
                          f"  spline {100*c['fraction_spline']:.3f}%")
            if plot:
                figure(out, cur, e, tungsten)


if __name__ == "__main__":
    main()
