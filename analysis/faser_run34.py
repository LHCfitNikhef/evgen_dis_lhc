#!/usr/bin/env python3
"""Run-3 and Run-3+4 event-rate projections for FASER$\\nu$.

WHAT THIS IS FOR (user, 2026-09-08): "I want a new section which is titled
'Applications to FASER', where we give predictions for FASER Run 3 and
Run 3 + Run 4, using the latest version of the fluxes"; and, on which fluxes:
"the latest version of the fluxes are those of arXiv:2402.13318 rescaled ...
Total luminosity there is still unclear, you can assume 680 fb^-1 for
Run 3 + Run 4, with some uncertainty."

>>> THE PROJECTION IS A SCALING, AND IT IS ONE ON PURPOSE. <<<  Eq. (rate) is
linear in the integrated luminosity: the flux files ARE a count of neutrinos
through their aperture at a stated luminosity, so a longer run multiplies
every number by L / L_flux and changes nothing else.  Nothing here re-derives
a cross section; it reads the rates analysis/faser_rates.py already produced
at 250 fb^-1 and states them again at 680.  That is why the luminosity being
provisional costs nothing: every number below scales as L / 680, and a
reader who learns a different figure rescales the table in their head.

WHAT DOES NOT SCALE, and is therefore the honest uncertainty on these
numbers, is the flux itself.  Ref. arXiv:2402.13318 quotes +992/-962 on
8507 charged-current nu_mu + nubar_mu interactions -- about 11.5% -- from the
spread of the forward-hadron-production models, and that is far larger than
anything the generator comparison in this benchmark has found.  It is
reported with every projected rate.

>>> AND THE DETECTOR IS THE PRESENT ONE. <<<  arXiv:2503.19775 proposes an
UPGRADED FASER neutrino detector for Run 4 and the HL-LHC era, with a
different target mass and different technology.  These numbers are the
present 1.1 t tungsten FASERnu run for longer, which is the conservative
baseline and the only one this benchmark's cross sections can be turned into
without inventing a detector.

Usage: analysis/faser_run34.py   ->  results_nu/faser_run34.json
"""
import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "analysis"))

import faser_rates                                           # noqa: E402

# The two integrated luminosities.  RUN3 is what every rate on the page is
# already quoted at, and is faser_rates.LUMI_FB rather than a second copy of
# the same number -- the two must never be able to drift apart.
LUMI_RUN3 = faser_rates.LUMI_FB
# Run 3 + Run 4, by user instruction and explicitly provisional.
LUMI_RUN34 = 680.0

# The flux uncertainty, from the reference the flux files come from.
#
# >>> IT IS A NEUTRINO NUMBER AND IT IS APPLIED TO NEUTRINOS ONLY. <<<  The
# neutrino flux is Ref. arXiv:2402.13318, which quotes a spread over
# forward-hadron-production models; the MUON flux is a different calculation
# in a different paper (arXiv:2506.13889) and carries no uncertainty in what
# it releases.  Putting the neutrino spread on the muon rates would attach a
# hadron-production uncertainty to a flux that does not come from hadron
# production -- a plausible-looking band belonging to another calculation.
# The muon rows therefore carry no flux band, and say so.
_R = faser_rates.REF_NU
FLUX_REL = (_R["numu_plus_numubar_err"][0] / _R["numu_plus_numubar"],
            _R["numu_plus_numubar_err"][1] / _R["numu_plus_numubar"])

CURRENTS = ("nu", "mu")
SELECTIONS = ("inclusive", "faser_s", "faser_e")


def _load(cur, sel):
    p = f"{BASE}/results_nu/faser_rates_{cur}_{sel}.json"
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def main():
    out = {"what": "FASERnu event-rate projections for Run 3 and Run 3+4",
           "lumi_run3_fb": LUMI_RUN3, "lumi_run34_fb": LUMI_RUN34,
           "lumi_run34_provisional": True,
           "flux_reference": _R["arxiv"],
           "flux_rel_err": {"hi": FLUX_REL[0], "lo": FLUX_REL[1],
                            "applies_to": "the neutrino current only"},
           "detector": "the present FASERnu, 1.1 t of tungsten",
           "upgrade_reference": "2503.19775",
           "rates": {}}
    for cur in CURRENTS:
        for sel in SELECTIONS:
            d = _load(cur, sel)
            if d is None:
                continue
            if float(d["lumi_fb"]) != LUMI_RUN3:
                sys.exit(f"faser_rates_{cur}_{sel}.json is quoted at "
                         f"{d['lumi_fb']} fb^-1, not {LUMI_RUN3} -- the "
                         f"projection would scale from the wrong baseline")
            k = f"{cur}_{sel}"
            out["rates"][k] = {"label": d.get("reference_label"),
                               "generators": {}}
            for g, v in d["generators"].items():
                n3 = v.get("total")
                if n3 is None:
                    continue
                f = LUMI_RUN34 / LUMI_RUN3
                rec = {"label": v.get("label"),
                       "run3": float(n3),
                       "run34": float(n3) * f}
                if cur == "nu":
                    rec["run34_hi"] = float(n3) * f * (1.0 + FLUX_REL[0])
                    rec["run34_lo"] = float(n3) * f * (1.0 - FLUX_REL[1])
                out["rates"][k]["generators"][g] = rec
    p = f"{BASE}/results_nu/faser_run34.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    for k, v in out["rates"].items():
        g = v["generators"]
        if "powheg_nu" in g or "powheg" in g:
            r = g.get("powheg_nu") or g.get("powheg")
            band = (f"(+{100*FLUX_REL[0]:.1f}/-{100*FLUX_REL[1]:.1f}% flux)"
                    if "run34_hi" in r else "(no flux band: muon flux)")
            print(f"  {k:22s} {r['label']:12s} "
                  f"{r['run3']:9.0f} -> {r['run34']:9.0f} {band}")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
