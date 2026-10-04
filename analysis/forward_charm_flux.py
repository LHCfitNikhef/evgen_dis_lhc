#!/usr/bin/env python3
"""The FASERnu neutrino flux from CHARM (and bottom) hadron decay, from the
POWHEG event files of arXiv:2309.12793 -- and the interaction rate it gives
with this benchmark's own cross-section.

>>> WHY (user, 2026-09-04): "the nu_e charm flux comes from POWHEG, this is
    the one that we should be using." <<<

arXiv:2402.13318 Table I builds its combination row out of two pieces:
EPOS-LHC for neutrinos from light-hadron decay, and POWHEG + Pythia 8.3 for
neutrinos from charm-hadron decay.  The second piece is where its whole
uncertainty lives (+911/-372 on nu_e, which is 30% charm) and it is the piece
the vendored histograms in `data/faser_flux/` get wrong: those are the
combination row of arXiv:2105.08270, whose charm comes from an average over
SIBYLL, DPMJET and two Pythia setups, and DPMJET alone predicts 3457 nu_e
against SIBYLL's 901.  Averaging them puts our nu_e 55% above Table I.

arXiv:2309.12793 -- Ref. [46] of that table, i.e. THE calculation it uses --
releases the neutrinos themselves rather than a histogram: energy, parent
hadron, and the transverse position at z = 480 m, for a 1 m x 1 m square
centred on the beam collision axis, with the 7-point scale variation attached
per primary event.  So this module can do what the vendored files cannot:

  * apply the FASERnu Run-3 aperture EXACTLY, from the positions in the file,
    instead of inheriting whatever aperture a histogram was made through;
  * separate charm from bottom;
  * carry the scale envelope that Table I quotes as its uncertainty.

>>> THE ONE THING TO GET RIGHT IS THE ORIGIN, AND IT IS MEASURED, NOT
    ASSUMED. <<<  The file gives (x0, y0) at z = 480 m, and the question is
what they are measured from.  They are NOMINAL coordinates with the crossing
angle already in the neutrino directions: the weighted transverse profile of
the sample peaks at y = -8 +- 1 cm and falls monotonically above it, which is
the 7.7 cm the Run 3 half-crossing angle of 160 urad downwards displaces the
true line of sight by (arXiv:2402.13318 Sec. II).  So the FASERnu aperture is
centred at (1.0, -3.3) cm -- the detector position that paper quotes, used
directly.  `profile()` prints the measurement every run; `ORIGINS` also
carries the two wrong choices, so that "nominal is the one that matches" is
something the output shows rather than something a comment asserts.

THE PROFILE IS FLAT, WHICH IS WHY THE APERTURE HARDLY MATTERS AND THE MASS
MATTERS COMPLETELY.  Over +-25 cm the flux density varies by under 20%, so
the rate is very nearly (flux density) x (target mass): the aperture cancels
between the number of neutrinos collected and the column density they see.
That is the same cancellation `faser_rates.column_density` turns on, and here
it means a 25 x 25 cm and a 25 x 30 cm detector of the SAME mass give the same
answer to a per cent.  Both masses are therefore evaluated: the Run 3
detector's 1.1 t of arXiv:2402.13318, and the 1.2 t of arXiv:2309.12793,
which describes FASERnu as "25 cm x 25 cm x 1 m ... with roughly 1.2 tons".

WHAT THE CROSS-SECTION IS.  Ours: the GENIE G18_02a total charged-current
splines of `tools/genie_faser_splines.sh`, on 74 p + 110 n, exactly as
`faser_genie_rates.py` uses them.  The authors ship their own GENIE tungsten
cross-sections beside the events; those are read too and printed as a ratio,
which is a second, independent version of the closure that module already
does against the flux files' interaction counts.

THE TARGET NEEDS NO ARCHAEOLOGY HERE, unlike the vendored histograms: the
aperture is applied to the events, so the column density is simply the Run-3
detector's, 1.1 t over 25 x 30 cm = 8.83e26 nucleons/cm^2.

Usage:
  tools/fetch_forward_charm.sh          # once, ~300 MB, gitignored
  analysis/forward_charm_flux.py [--bins N]
Writes results_nu/forward_charm_flux.json and data/faser_flux/POWHEG_<...>.txt
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import faser_rates as fr          # noqa: E402
import faser_genie_rates as fg    # noqa: E402
import target as _tgt             # noqa: E402

DATA = f"{BASE}/data/forward_charm"
TAG = "13.6TeV_POWHEG-smallxLHCb_P8-Monash"

# FASERnu in Run 3, arXiv:2402.13318 Sec. II: 25 cm wide, 30 cm high, filled
# with 1.1 t of tungsten, centred at (1.0, -3.3) cm in NOMINAL coordinates
# with the true line of sight 7.7 cm below the nominal one.
W_CM, H_CM = 25.0, 30.0
CENTRE_NOMINAL_CM = (1.0, -3.3)
LOS_SHIFT_CM = (0.0, -7.7)         # true LOS relative to nominal, Run 3
MASS_G = fr.FASERNU_MASS_G         # 1.1 tonnes
def column(mass_g, area_cm2=W_CM * H_CM):
    """nucleons/cm^2 -- the aperture is applied to the events, so this is
    simply the mass over the face it was applied through."""
    return mass_g / _tgt.A_W_GMOL * fr.N_A * _tgt.A_W / area_cm2


COLUMN = column(MASS_G)

# Table I of arXiv:2402.13318, charm-hadron rows, 250 fb^-1
REF_CHARM = {"nue": {"max": 1405.0, "central": 527.0, "min": 294.0},
             "numu": {"max": 1373.0, "central": 511.0, "min": 284.0},
             "nutau": {"max": 76.0, "central": 28.0, "min": 16.0}}
SCALE_COLS = ["11", "12", "21", "22", "H1", "1H", "HH"]

# The three origins the aperture is tried at.  "nominal" is the physical one
# -- the file's coordinates are nominal, which profile() measures -- and the
# other two are kept so the choice is visibly a measurement.
ORIGINS = {
    "nominal": CENTRE_NOMINAL_CM,
    "true_los": (CENTRE_NOMINAL_CM[0] - LOS_SHIFT_CM[0],
                 CENTRE_NOMINAL_CM[1] - LOS_SHIFT_CM[1]),
    "on_axis": (0.0, 0.0),
}
DEFAULT_ORIGIN = "nominal"

# The two target masses in play, and they are not the same detector: the Run 3
# FASERnu of arXiv:2402.13318 is 1.1 t, while arXiv:2309.12793 -- the source of
# these events, and the reference Table I's charm rows are taken from -- says
# 1.2 t.  Since the aperture cancels, the mass is the whole normalisation.
MASSES = {"run3_1.1t": fr.FASERNU_MASS_G, "source_1.2t": 1.2e6}

FLAV = {"nue": (12, -12), "numu": (14, -14), "nutau": (16, -16)}


def load(label):
    """Every neutrino of one parent species, and the scale-variation table."""
    files = sorted(glob.glob(f"{DATA}/events_{label}_{TAG}_*.csv.gzip"))
    if not files:
        raise SystemExit(f"no {label} event files in {DATA} -- run "
                         "tools/fetch_forward_charm.sh")
    ev = pd.concat([pd.read_csv(f, compression="gzip") for f in files],
                   ignore_index=True)
    wf = f"{DATA}/weights_{label}_{TAG}.csv.gzip"
    rw = pd.read_csv(wf, compression="gzip")
    return ev, rw


def in_aperture(ev, centre_cm):
    """The FASERnu face about one origin, in the file's own coordinates."""
    cx, cy = centre_cm[0] / 100.0, centre_cm[1] / 100.0
    hx, hy = W_CM / 200.0, H_CM / 200.0
    return ((ev["x0"] > cx - hx) & (ev["x0"] < cx + hx)
            & (ev["y0"] > cy - hy) & (ev["y0"] < cy + hy))


def sigma_w(pid):
    """Our GENIE total CC per tungsten NUCLEUS [pb], on a fine energy grid.

    None if the splines for that beam do not exist.  A MISSING flavour is
    reported and dropped rather than silently given another one's
    cross-section: the tau's is suppressed by its own mass at these energies
    and borrowing the muon's would overstate the nu_tau rate by a large,
    plausible-looking factor.
    """
    try:
        e, sw, _sp, _sn, _n = fg.sigma_w_per_nucleus(pid)
    except SystemExit as exc:
        print(f"  [no cross-section for pid {pid}] {exc}")
        return None
    return e, sw


def sigma_authors(pid):
    """The authors' own GENIE tungsten cross-section, in cm^2 per NUCLEON.

    Their column is used by their notebook as sigma with rho L / m_p, i.e.
    per nucleon; returned here only for the ratio against ours.
    """
    p = f"{DATA}/xs_GENIE_W_{pid}.txt"
    if not os.path.exists(p):
        return None
    d = np.loadtxt(p)
    return d[:, 0], d[:, 1]


def profile(ev, axis="y0", lo=-25.0, hi=25.0, nb=25):
    """The weighted transverse profile, in cm -- the evidence for the origin.

    Returns (centres, sum of weights per bin).  The sample carries NEGATIVE
    weights (it is an NLO calculation), so the estimator is the sum of
    weights, never a count and never a weighted mean, which the negative tail
    makes meaningless.
    """
    v = ev[axis].to_numpy() * 100.0
    b = np.linspace(lo, hi, nb + 1)
    h, _ = np.histogram(v, bins=b, weights=ev["w"].to_numpy())
    return 0.5 * (b[1:] + b[:-1]), h


def rate_table(ev, rw, lumi_fb, centre_cm, per_nucleus, col=None):
    """Events at `lumi_fb`, per flavour, with the scale envelope.

    `per_nucleus` maps a pid to (E knots, sigma [pb] per W nucleus) -- ours.
    The scale variation multiplies the PRIMARY event weight, so it is applied
    through iEvent, which is what the weights file is keyed on.
    """
    sel = ev[in_aperture(ev, centre_cm)].copy()
    # events at this luminosity, before the scale reweighting: w is in pb
    sel["n"] = sel["w"] * lumi_fb * 1000.0
    out = {}
    for name, pids in FLAV.items():
        if not any(p in per_nucleus for p in pids):
            continue
        part = sel[sel["vpid"].isin(pids)]
        if part.empty:
            out[name] = {"through": 0.0, "interacting": 0.0, "scales": {}}
            continue
        n_int = np.zeros(len(part))
        for pid in pids:
            m = (part["vpid"] == pid).to_numpy()
            if not m.any() or pid not in per_nucleus:
                continue
            e_k, s_k = per_nucleus[pid]
            ok = s_k > 0
            s = np.exp(np.interp(np.log(part["en"].to_numpy()[m]),
                                 np.log(e_k[ok]), np.log(s_k[ok])))
            # sigma [pb] -> cm^2, times nuclei/cm^2
            n_int[m] = (part["n"].to_numpy()[m] * s * fr.PB_TO_CM2
                        * (COLUMN if col is None else col) / _tgt.A_W)
        part = part.assign(n_int=n_int)
        # the scale envelope: one ratio per PRIMARY event, so join on iEvent
        joined = part.merge(rw, on="iEvent", how="left")
        scales = {}
        # NOT `col`: that is this function's column-density argument, and
        # shadowing it here silently fed a two-character string into the rate
        # on the second flavour
        for sc in SCALE_COLS:
            if sc not in joined:
                continue
            r = joined[sc].fillna(1.0).to_numpy()
            scales[sc] = float((joined["n_int"].to_numpy() * r).sum())
        tot = float(part["n_int"].sum())
        vals = list(scales.values()) or [tot]
        out[name] = {"through": float(part["n"].sum()),
                     "interacting": tot,
                     "scales": scales,
                     # the 7-point envelope, stored on EVERY record rather
                     # than only the one main() prints -- a figure that reads
                     # a second configuration must not have to recompute it
                     "envelope": [min(vals), max(vals)],
                     "envelope_over_central": [min(vals) / tot, max(vals) / tot],
                     "n_mc": int(len(part))}
    return out, sel


def main():
    lumi = fr.LUMI_FB
    per_nucleus = {pid: sigma_w(pid) for pid in (12, -12, 14, -14, 16, -16)}
    per_nucleus = {k: v for k, v in per_nucleus.items() if v is not None}
    out = {"source": {"arxiv": "2309.12793",
                      "what": "POWHEG + Pythia 8.3 forward heavy hadrons, "
                              "neutrinos through 1 m^2 at z = 480 m",
                      "tag": TAG},
           "reference": {"arxiv": "2402.13318",
                         "table": "Table I, charm-hadron rows, 250 fb^-1",
                         "charm": REF_CHARM},
           "detector": {"width_cm": W_CM, "height_cm": H_CM,
                        "masses_g": MASSES, "origin": DEFAULT_ORIGIN,
                        "centre_cm": ORIGINS[DEFAULT_ORIGIN],
                        "centre_nominal_cm": CENTRE_NOMINAL_CM,
                        "los_shift_cm": LOS_SHIFT_CM,
                        "nucleons_per_cm2": {k: column(v)
                                             for k, v in MASSES.items()}},
           "lumi_fb": lumi, "cross_section_closure": {}, "species": {}}

    # >>> OUR CROSS-SECTION AGAINST THE FLUX AUTHORS' OWN, DIRECTLY. <<<
    # They ship their GENIE tungsten cross-section beside the events, so this
    # is not the indirect closure faser_genie_rates.py does against their
    # interaction counts -- it is the same quantity, file against file.
    print("our GENIE tungsten cross-section against the authors' own "
          "(data/forward_charm/xs_GENIE_W_*.txt):")
    for pid in sorted(per_nucleus, key=abs):
        a = sigma_authors(pid)
        if a is None:
            continue
        e_k, s_k = per_nucleus[pid]
        ok = s_k > 0
        rows = []
        for E in (30.0, 100.0, 300.0, 1000.0, 3000.0):
            ours = float(np.exp(np.interp(np.log(E), np.log(e_k[ok]),
                                          np.log(s_k[ok] / _tgt.A_W))))
            theirs = float(np.interp(E, a[0], a[1])) / fr.PB_TO_CM2
            rows.append((E, theirs / ours))
        out["cross_section_closure"][str(pid)] = {
            "E_GeV": [r[0] for r in rows], "theirs_over_ours": [r[1] for r in rows]}
        print(f"  pid {pid:3d}: " + "  ".join(f"{int(E)} GeV {r:.4f}"
                                              for E, r in rows))

    print(f"\nFASERnu {W_CM:g} x {H_CM:g} cm at "
          f"({ORIGINS[DEFAULT_ORIGIN][0]:+.1f}, {ORIGINS[DEFAULT_ORIGIN][1]:+.1f}) cm, "
          f"{lumi:g} fb^-1; column "
          + ", ".join(f"{k} {column(v):.3e}" for k, v in MASSES.items()))

    for label in ("charm", "bottom"):
        try:
            ev, rw = load(label)
        except SystemExit as exc:
            print(f"[{label}] {exc}")
            continue
        print(f"\n{label}: {len(ev):,} neutrinos through the 1 m^2 square, "
              f"{len(rw):,} primary events with scale weights")
        if label == "charm":
            c, h = profile(ev)
            k = int(np.argmax(h))
            print("  transverse profile in y [cm], sum of weights normalised "
                  "to its peak:")
            print("   " + " ".join(f"{cc:+.0f}:{hh/h.max():.2f}"
                                   for cc, hh in zip(c, h)))
            print(f"  peak at y = {c[k]:+.1f} cm -- the crossing angle puts "
                  f"the true line of sight at {LOS_SHIFT_CM[1]:+.1f} cm, so "
                  f"the file's coordinates are NOMINAL and the detector "
                  f"centre is used as quoted")
            out["profile_y_cm"] = {"centres": c.tolist(),
                                   "weights": h.tolist(), "peak_cm": float(c[k])}
        rec = {}
        for mname, mass in MASSES.items():
            col = column(mass)
            for oname, centre in ORIGINS.items():
                tab, _sel = rate_table(ev, rw, lumi, centre, per_nucleus, col)
                rec[f"{mname}/{oname}"] = tab
                if oname != DEFAULT_ORIGIN and mname != "run3_1.1t":
                    continue
                line = f"  {mname:11s} {oname:9s}:"
                for name in tab:
                    v = tab[name]["interacting"]
                    ref = REF_CHARM[name]["central"] if label == "charm" else None
                    line += (f"  {name} {v:8.1f}"
                             + (f" ({v/ref:.3f})" if ref else ""))
                print(line)
        out["species"][label] = rec

    # >>> THE ANSWER, on the measured origin and the Run 3 mass <<<
    if "charm" in out["species"]:
        key = f"run3_1.1t/{DEFAULT_ORIGIN}"
        alt = f"source_1.2t/{DEFAULT_ORIGIN}"
        got, got12 = out["species"]["charm"][key], out["species"]["charm"][alt]
        print()
        print(f"{'flavour':8s} {'ours (1.1 t)':>13s} {'ours (1.2 t)':>13s} "
              f"{'scale envelope':>21s} {'Table I charm':>21s} {'1.1t/ref':>9s}")
        for name in got:
            v, v12 = got[name], got12[name]
            lo, hi = v["envelope"]
            r = REF_CHARM[name]
            print(f"{name:8s} {v['interacting']:13.1f} {v12['interacting']:13.1f} "
                  f"{lo:9.1f} -{hi:10.1f} "
                  f"{r['min']:6.0f} {r['central']:6.0f} {r['max']:6.0f}  "
                  f"{v['interacting']/r['central']:9.3f}")
            # the SHAPE of the scale band is a sharper test than its
            # normalisation: it is the calculation, not the detector
            v["over_reference"] = v["interacting"] / r["central"]
            v["reference_envelope_over_central"] = [r["min"] / r["central"],
                                                    r["max"] / r["central"]]
        print("\nthe scale envelope as a RATIO to the central value -- the "
              "test of the calculation rather than of the detector:")
        for name in got:
            a = got[name]["envelope_over_central"]
            b = got[name]["reference_envelope_over_central"]
            print(f"  {name:8s} ours {a[0]:.3f} - {a[1]:.3f}   "
                  f"Table I {b[0]:.3f} - {b[1]:.3f}")

    p = f"{BASE}/results_nu/forward_charm_flux.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
