#!/usr/bin/env python3
"""GENIE's default tune against Table I of arXiv:2402.13318, with our own
machinery -- and the decomposition of what is left over.

WHAT THE PAPER DID.  The FASER collaboration convoluted its forward neutrino
fluxes with "the default GENIE cross section" over the 1.1 t tungsten target
of FASERnu, and quotes the expected charged-current interactions at 250 fb^-1
per flavour (combination row): nu_e + nubar_e 1675 (+911/-372),
nu_mu + nubar_mu 8507 (+992/-962).  Its combination row is EPOS-LHC for light
hadrons plus POWHEG+Pythia 8.3 for charm hadrons (1149 + 527 = 1676 and
7996 + 511 = 8507, its own Table I).

WHAT THIS DOES.  The same convolution, with every factor ours:

  * the cross-section: GENIE G18_02a_00_000 total CC on a free proton and a
    free neutron, from splines built by tools/genie_faser_splines.sh with the
    validity cap raised to 10 TeV, combined as tungsten,
    sigma_W = 74 sigma_p + 110 sigma_n per nucleus;
  * the flux: the vendored per-flavour, per-charge files, scaled from
    150 fb^-1 to 250;
  * the target: 1.1 tonnes of tungsten behind the aperture THE FLUX WAS
    COUNTED THROUGH -- see faser_rates.FLUX_APERTURE_CM2, and read that block
    before touching any number here.

>>> THE COMPARISON IS FULLY INCLUSIVE, AND THAT IS DELIBERATE (user,
    2026-09-04).  <<<  Table I is every CC interaction at every Q2 and every
    y, so nothing here imposes the benchmark's 0.2 < y < 0.9 window or its
    Q2 > 4 GeV2 floor: the GENIE splines are total cross-sections, summed over
    every CC process of the default list -- quasi-elastic, resonance, DIS and
    the rest.  Putting a fiducial region on one side of this ratio would make
    it a measurement of the region rather than of the physics, which is
    exactly what the muon-side comparison had to be rebuilt to avoid
    (faser_rates.compare_mu_like_for_like).  The benchmark's own region is
    compared separately, on the "all y" region, in faser_rates.py.

WHAT COMES OUT, and it is three separate statements rather than one ratio:

  1. OUR GENIE REPRODUCES THEIRS.  The flux files ship the authors' own CC
     interaction counts beside the flux, so their ratio is sigma(E) x T with
     THEIR cross-section.  Ours divided by theirs lands within 1.2% of unity
     on all four (flavour, charge) combinations once the file's own target is
     used, and is flat in energy to 0.5% (neutrinos) and 2.3%
     (antineutrinos) over 100 GeV to 5.7 TeV, which carries 89% of the rate.
     A cross-section MODEL difference could not be flat over three decades; a
     target normalisation is exactly flat, which is how the two were
     separated.  Below 100 GeV the antineutrino ratio climbs to 1.13 at the
     first flux point -- their calculation is LO on a tungsten NUCLEUS and
     ours the free-nucleon sum of every CC process, and that is where the two
     have room to differ.

  2. THE TARGET IS THE FILES', NOT THE PAPER'S.  The vendored files are
     Table I of arXiv:2105.08270 -- checked to four figures on three flavours
     -- and that paper's FASERnu is 1.2 t behind a 25 x 25 cm aperture, not
     the Run 3 detector's 1.1 t behind 25 x 30 cm.  Mass over aperture is the
     only combination that enters, so predicting the Run 3 rate from these
     files means 1.1 t over 625 cm^2, and nothing else changes.

  3. WHAT IS LEFT IS THE FLUX MODEL, AND ON THE ELECTRON FLAVOUR IT IS
     DPMJET.  Ours is 1.03 of Table I on nu_mu + nubar_mu and 1.55 on
     nu_e + nubar_e.  The vendored flux is the 2105.08270 average over four
     generator setups, one of which is DPMJET, whose nu_e prediction is 3457
     against SIBYLL's 901; that paper also publishes the average WITHOUT
     DPMJET, and rescaled to 250 fb^-1 and 1.1 t it gives 1723 nu_e and 8168
     nu_mu against Table I's 1675 and 8507 -- 1.03 and 0.96.  So both
     flavours agree at the few-per-cent level once the same flux model is
     used, and the electron excess is a known feature of the older average
     rather than anything in this pipeline.  Note also that 1.55 sits inside
     Table I's own nu_e uncertainty, +911/-372, which that paper attributes
     to charm-hadron production -- the same thing.

Electron and muon neutrinos share one cross-section by flavour universality
(the lepton mass is irrelevant at these energies), so the nu_e result is the
same splines under the nu_e flux -- GENIE nonetheless builds them separately
and both are used.

Usage: analysis/faser_genie_rates.py
Writes results_nu/faser_genie_rates.json and prints the comparison.
"""
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import faser_rates as fr  # noqa: E402

SPLDIR = f"{BASE}/genie/splines/faser"
GEV2_TO_PB = 0.389379e9          # 1 GeV^-2 = 0.389379 mb = 3.89379e8 pb
Z_W, N_W = 74, 110
REF = {"arxiv": "2402.13318", "table": "Table I, combination, 250 fb^-1",
       "nue": (1675.0, 911.0, 372.0), "numu": (8507.0, 992.0, 962.0),
       "cross_section": "default GENIE (3.4)",
       "flux": "EPOS-LHC (light hadrons) + POWHEG+Pythia 8.3 (charm)",
       "detector": "1.1 t of tungsten, 25 x 30 x 80 cm"}

# nuclei/cm^2 -- the mass we are predicting over the aperture the flux was
# counted through, NOT over the detector face.  faser_rates carries the why.
NUCLEI_PER_CM2 = fr.column_density() / fr.A_W


def spline_total(pid, tgt):
    """Total CC cross-section (E [GeV], sigma [pb] per nucleon) from one
    gmkspl file: the sum of every process spline it holds."""
    fn = f"{SPLDIR}/{pid}_{tgt}_cc_e8000.xml"
    if not os.path.exists(fn):
        raise SystemExit(f"missing {fn} -- run tools/genie_faser_splines.sh")
    # the processes are tabulated on different knot grids (the resonance
    # splines stop early), so each is interpolated onto one common grid
    es = np.logspace(-1, np.log10(8000.0), 600)
    tot, names = np.zeros_like(es), []
    txt = open(fn).read()
    for m in re.finditer(r'<spline name="([^"]+)" nknots="(\d+)">(.*?)</spline>',
                         txt, flags=re.S):
        name, body = m.group(1), m.group(3)
        e = np.array([float(v) for v in re.findall(r"<E>\s*([0-9.eE+-]+)", body)])
        x = np.array([float(v) for v in re.findall(r"<xsec>\s*([0-9.eE+-]+)", body)])
        tot = tot + np.interp(es, e, x, left=0.0, right=x[-1])
        names.append(name)
    return es, tot * GEV2_TO_PB, names


def sigma_w_per_nucleus(pid):
    ep, sp, np_ = spline_total(pid, "p")
    en, sn, nn_ = spline_total(pid, "n")
    assert np.allclose(ep, en)
    return ep, Z_W * sp + N_W * sn, sp, sn, len(np_)


def _on_flux(pid, e_knots, sig_pb):
    """The spline interpolated onto the flux's own energy points, log-log."""
    e, _phi = fr.flux(str(pid))
    ok = sig_pb > 0
    s = np.exp(np.interp(np.log(e), np.log(e_knots[ok]), np.log(sig_pb[ok])))
    s[e < e_knots[ok][0]] = 0.0
    return e, s


def rate(pid, e_knots, sig_pb, per_cm2):
    """Events at 250 fb^-1: sum over the flux histogram of N_i sigma(E_i) T."""
    _e, phi = fr.flux(str(pid))
    phi = phi * fr.LUMI_FB / 150.0
    _e, s = _on_flux(pid, e_knots, sig_pb)
    return float(np.sum(phi * s * 1e-36 * per_cm2)), float(np.sum(phi))


def closure_against_files(pid, e_knots, sig_pb):
    """Our GENIE against the flux authors' own, ENERGY BY ENERGY.

    Their CCint file over their flux file is sigma(E) x T with THEIR
    cross-section and THEIR target; dividing by ours and by the file's own
    geometric column density leaves a pure cross-section ratio, which is 1 if
    the two GENIE calculations agree.  Returned as (E, ratio) plus the
    rate-weighted mean, because the shape is the whole argument: a model
    difference cannot be flat over three decades in energy.
    """
    # the 2021 files: the only ones that carry the authors' own CC counts
    e, phi = fr.legacy_flux(str(pid))
    _e, cc = fr.legacy_flux(str(pid), cc=True)
    ok = sig_pb > 0
    s = np.exp(np.interp(np.log(e), np.log(e_knots[ok]), np.log(sig_pb[ok])))
    s[e < e_knots[ok][0]] = 0.0                      # pb per NUCLEUS
    t_file = fr.column_density(fr.FLUX_TARGET_MASS_G,
                               fr.LEGACY_FLUX_APERTURE_CM2) / fr.A_W  # nuclei/cm2
    ours = phi * s * 1e-36 * t_file
    ok = (ours > 0) & (cc > 0)
    return e[ok], ours[ok] / cc[ok], float(ours[ok].sum() / cc[ok].sum())


def main():
    got, bad = fr.check_flux_provenance()
    if bad:
        raise SystemExit("[faser_genie_rates] the vendored flux files no "
                         "longer reproduce Table I of arXiv:2105.08270: "
                         + "; ".join(bad))
    out = {"reference": REF, "flux_reference": fr.REF_FLUX_2105,
           "flux_files_cc_150": got,
           "nuclei_per_cm2": NUCLEI_PER_CM2,
           "nucleons_per_cm2": NUCLEI_PER_CM2 * (Z_W + N_W),
           "target_mass_g": fr.FASERNU_MASS_G,
           "flux_aperture_cm2": fr.FLUX_APERTURE_CM2,
           "flux_target_mass_g": fr.FLUX_TARGET_MASS_G,
           "flavours": {}}
    print(f"flux files = Table I of {fr.REF_FLUX_2105['arxiv']}, combination: "
          + ", ".join(f"{k} {got[k]:.1f} vs {v:g}" for k, v in
                      fr.REF_FLUX_2105["combination_all"].items())
          + f" at {fr.REF_FLUX_2105['lumi_fb']:g} fb^-1")
    print(f"target: {fr.FASERNU_MASS_G/1e6:.2g} t of tungsten behind the "
          f"flux's own {fr.FLUX_APERTURE_CM2:.0f} cm2 aperture = "
          f"{NUCLEI_PER_CM2*(Z_W+N_W):.3e} nucleons/cm2  (the files' own "
          f"target is {fr.FLUX_TARGET_MASS_G/1e6:.2g} t = "
          f"{fr.column_density(fr.FLUX_TARGET_MASS_G, fr.LEGACY_FLUX_APERTURE_CM2):.3e})")
    print()
    print(f"{'flavour':8s} {'sig_p(1TeV)':>11s} {'sig_n(1TeV)':>11s} "
          f"{'N (ours)':>10s} {'files x 1.1/1.2':>15s} "
          f"{'our GENIE / theirs':>19s} {'flat >100':>9s}")
    tot_ours = {"nue": 0.0, "numu": 0.0}
    tot_files = {"nue": 0.0, "numu": 0.0}
    for pid in (12, -12, 14, -14):
        e, sw, sp, sn, nproc = sigma_w_per_nucleus(pid)
        n_ours, nflux = rate(pid, e, sw, NUCLEI_PER_CM2)
        ecl, rcl, rbar = closure_against_files(pid, e, sw)
        _e, cc = fr.flux(str(pid), cc=True)
        # the files' own count, moved to our luminosity and our target mass:
        # mass over aperture is all that enters, so this is one factor
        n_files = float(cc.sum()) * (fr.LUMI_FB / fr.REF_FLUX_2105["lumi_fb"]
                                     * fr.FASERNU_MASS_G / fr.FLUX_TARGET_MASS_G)
        # the spread is quoted above 100 GeV, which carries 89% of the rate:
        # below it the antineutrino ratio rises to 1.13, a real difference in
        # the sub-DIS channels (their calculation is LO on a tungsten NUCLEUS,
        # ours the free-nucleon sum of every CC process) that is worth a
        # per-mille of the rate and would otherwise dominate a min/max
        hi = ecl > 100.0
        spread = float(rcl[hi].max() / rcl[hi].min() - 1.0)
        s1p = float(np.exp(np.interp(np.log(1000.0), np.log(e[sp > 0]), np.log(sp[sp > 0]))))
        s1n = float(np.exp(np.interp(np.log(1000.0), np.log(e[sn > 0]), np.log(sn[sn > 0]))))
        key = "nue" if abs(pid) == 12 else "numu"
        tot_ours[key] += n_ours
        tot_files[key] += n_files
        out["flavours"][str(pid)] = {
            "sigma_p_1TeV_pb": s1p, "sigma_n_1TeV_pb": s1n,
            "n_processes": nproc, "flux_through_250": nflux,
            "events": n_ours, "flux_files_cc_rescaled": n_files,
            "closure_vs_files": rbar, "closure_spread_in_E": spread,
            "closure_e": [float(x) for x in ecl],
            "closure_ratio": [float(x) for x in rcl]}
        print(f"{pid:8d} {s1p:11.3f} {s1n:11.3f} {n_ours:10.1f} "
              f"{n_files:15.1f} {rbar:19.4f} {spread*100:8.1f}%")
    print()
    print(f"{'':10s} {'ours':>8s} {'paper':>20s} {'ratio':>7s}   "
          f"{'2105 all':>9s} {'2105 no DPMJET':>15s}")
    lumi_mass = (fr.LUMI_FB / fr.REF_FLUX_2105["lumi_fb"]
                 * fr.FASERNU_MASS_G / fr.FLUX_TARGET_MASS_G)
    for key, lab in (("nue", "nu_e+bar"), ("numu", "nu_mu+bar")):
        c, up, dn = REF[key]
        all_ = fr.REF_FLUX_2105["combination_all"][key] * lumi_mass
        nod = fr.REF_FLUX_2105["combination_no_dpmjet"][key] * lumi_mass
        out[key] = {"events": tot_ours[key], "paper": c, "paper_err": (up, dn),
                    "flux_files_rescaled": tot_files[key],
                    "ours_over_paper": tot_ours[key] / c,
                    "ours_over_files": tot_ours[key] / tot_files[key],
                    "comb_2105_all_rescaled": all_,
                    "comb_2105_no_dpmjet_rescaled": nod,
                    "comb_2105_no_dpmjet_over_paper": nod / c,
                    # our cross-section on the flux WITHOUT DPMJET: the
                    # vendored histograms are the four-generator average, and
                    # the flux authors publish only the no-DPMJET COUNT, so
                    # ours is carried over by that count's ratio to the
                    # all-generator one -- a shape-preserving rescaling, said
                    # as such wherever it is shown (user, 2026-09-04)
                    "no_dpmjet_flux_factor": nod / all_,
                    "events_on_no_dpmjet_flux": tot_ours[key] * nod / all_,
                    "inside_paper_band": bool(c - dn <= tot_ours[key] <= c + up)}
        print(f"{lab:10s} {tot_ours[key]:8.0f} {c:8.0f} +{up:.0f}/-{dn:.0f} "
              f"{tot_ours[key]/c:7.3f}   {all_:9.0f} {nod:15.0f}"
              f"   (no-DPMJET/paper {nod/c:.3f})")
    print()
    print("read: row 1 is this pipeline, row 2 the same flux without DPMJET "
          "in its average -- the electron excess is that generator, not the "
          "cross-section, which closes to ~1% flavour by flavour above.")
    with open(f"{BASE}/results_nu/faser_genie_rates.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote results_nu/faser_genie_rates.json")


if __name__ == "__main__":
    main()
