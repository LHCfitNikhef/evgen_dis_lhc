# Forward heavy-hadron neutrinos (POWHEG-BOX heavy-quark + Pythia 8)

Fetched by `tools/fetch_forward_charm.sh` from
**KlingFelix/forward_heavy_hadrons_NLONLLx**, `NeutrinoFluxMC/`, the data
release of

> L. Buonocore, F. Kling, L. Rottoli, J. Sominka, *Predictions for neutrinos
> and new physics from forward heavy hadron production at the LHC*,
> [arXiv:2309.12793](https://arxiv.org/abs/2309.12793), Eur. Phys. J. C **84**
> (2024) 363.

## Why they are here

`data/faser_flux/` averages four hadronic-interaction models for the
charm-hadron neutrinos, and DPMJET's ν_e prediction is far above the others.
The FASER prediction of [arXiv:2402.13318](https://arxiv.org/abs/2402.13318)
(Table I) instead uses EPOS-LHC for light hadrons plus this POWHEG-BOX
heavy-quark + Pythia 8.3 calculation for charm hadrons, so these files are the matching charm
component.  Used by `analysis/forward_charm_flux.py`.

## Files

    events_<charm|bottom>_13.6TeV_POWHEG-smallxLHCb_P8-Monash_<i>.csv.gzip
    weights_<charm|bottom>_13.6TeV_POWHEG-smallxLHCb_P8-Monash.csv.gzip
    xs_GENIE_W_<pid>.txt

One row per neutrino crossing a 1 m × 1 m square centred on the beam axis at
z = 480 m, √s = 13.6 TeV: `vpid`, `hpid` (parent hadron), `x0`, `y0` [m],
`thx`, `thy` [mrad], `en` [GeV], `w` [pb], `iEvent` (key into the weights
file, which holds the seven-point scale variation as ratios).
`xs_GENIE_W_<pid>.txt` are the authors' GENIE tungsten CC cross-sections, cm²
per nucleon against E in GeV.

## Things to know

- **Negative weights**: it is an NLO calculation, so every estimator is a sum
  of weights, never a count or a weighted mean.
- **Coordinates are nominal**, with the crossing angle already in the
  neutrino directions (the profile peaks at y ≈ −8 cm), so the FASERν centre
  (1.0, −3.3) cm of arXiv:2402.13318 is used as it stands.
- **The aperture cancels, the target mass does not**: the flux density is
  flat to under 20% over ±25 cm.  The two papers quote 1.1 t and 1.2 t, a 9%
  difference; both are computed.

## Not in git

~300 MB of generator output, excluded by `.gitignore`.  This README, the
fetch script and the derived `results_nu/forward_charm_flux.json` are tracked.
Cite [2309.12793](https://arxiv.org/abs/2309.12793) for these files.
