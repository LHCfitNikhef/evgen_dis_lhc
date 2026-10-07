# FASERν neutrino flux, Run 3: EPOS-LHC light + POWHEG charm (Kling 2025)

The neutrino flux used for every FASER rate in Sect. 7 of the paper (since
2026-10-07).  It is the FASER default of the Run-3 rate predictions,
arXiv:2402.13318: **EPOS-LHC** for light hadrons (pions, kaons) and the
NLO+NLL_x **POWHEG + Pythia 8** calculation of arXiv:2309.12793 for charm
hadrons, propagated to 480 m with the fast neutrino flux simulation of
F. Kling and L. J. Nevay, arXiv:2105.08270.

## Files

Energy spectra through the **FASERν tungsten face** (25 × 30 cm,
x in [-115, 135] mm, y in [-171, 129] mm around the line of sight), **minus
160 µrad crossing angle** (the configuration FASER's own GENIE production
uses), 13.6 TeV:

    EPOSLHC_light_fasernu25x30.txt                 light-hadron parents
    Powheg_pythia_charm_central_fasernu25x30.txt   charm parents, central scale
    Powheg_pythia_charm_{min,max}_fasernu25x30.txt charm, scale-variation envelope

Columns: `E_lo_GeV E_hi_GeV` then `N err` for nu_e, nubar_e, nu_mu, nubar_mu,
nu_tau, nubar_tau.  N is **neutrinos per fb^-1 in the bin** (not per GeV),
err its MC statistical error.  10 log bins per decade, 10 GeV - 10 TeV.
A complete flux is one light file PLUS one charm file.

## Provenance

Copied verbatim from `$FASER_DATA/fluxes/2025/text/spectra/minus/`, written by
`faser_format/flux_to_text.py` from F. Kling's ROOT ntuples, whose download
links are those of FASER's GENIE fork (branch `faser-R-3_04_00`,
`faser/getFluxNtp_minus_160_urad_crossing_angle.sh`).  See the README there.

## Normalisation

The spectra are COUNTS through the 25 × 30 cm face, so the rate is

    N = L x Phi_file(E) x sigma(E) x M / (25 x 30 cm^2),

with M the target mass (1.1 t of tungsten, arXiv:2402.13318), as
`analysis/faser_rates.py` does (`FLUX_APERTURE_CM2 = 750`).  This replaces
`data/faser_flux/` (the four-generator average of arXiv:2105.08270, counted
through 25 × 25 cm), which is kept only for the cross-checks against that
paper's Table I.

Cite arXiv:2105.08270, arXiv:2309.12793, arXiv:2402.13318 and EPOS-LHC
(arXiv:1306.0121).
