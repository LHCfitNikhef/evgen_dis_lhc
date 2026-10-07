# LHC neutrino fluxes at FASERν -- LEGACY (2021 average)

> **Superseded on 2026-10-07.** Every rate now uses `data/faser_flux_2025/`
> (EPOS-LHC light + POWHEG charm, the arXiv:2402.13318 default; user rule).
> These files remain only for the cross-checks against Table I of
> arXiv:2105.08270 (`faser_rates.legacy_flux`).

Vendored from **KlingFelix/FastNeutrinoFluxSimulation**, `Fluxes/FASERv/`.
Kept in the repository (data, not code) for those cross-checks.

## Files

Energy histograms, two columns, `E [GeV]` and number of neutrinos, for
**LHC Run 3, 150 fb⁻¹**, at the FASERν location 480 m downstream of the ATLAS
interaction point; 28 logarithmic points from 11.3 to 5698 GeV.

    FASER_<pid>.txt         neutrinos passing through the detector
    FASER_CCint_<pid>.txt   their CC interactions, as computed by the authors

PDG codes ±12, ±14, ±16.  The passing-through files are used; the CCint files
are only a cross-check.

## Normalisation trap

The flux is a **count through the authors' own 25 × 25 cm aperture**
(arXiv:2105.08270), not a density.  It must be paired with the target column
density *mass / that aperture* (1.2 t behind 25 × 25 cm), never with the
25 × 30 cm detector face, which silently gives a 17% error.  Fitting the
column density to the CCint files instead comes out 6% high.

## Provenance and citation

Fast neutrino flux simulation of F. Kling and L. J. Nevay, arXiv:2105.08270.
Cite [2105.08270](https://arxiv.org/abs/2105.08270) and
[2402.13318](https://arxiv.org/abs/2402.13318).

The authors' caveat: "These files have been produced some time ago.  They
therefore do not account for recent changes, for example regarding the
location, crossing angle, or luminosity (LHC Run 3 is longer than initially
planned).  For more precise numbers, please contact the authors."  Rates
derived from them are estimates at the tens-of-per-cent level in
normalisation; the benchmark's contribution is the cross-section factor.
