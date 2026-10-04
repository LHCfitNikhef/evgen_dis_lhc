# FASER's electronic-detector muon-neutrino measurement

`faser_electronic.json` holds FASER's measurement of muon-neutrino
charged-current interactions with its active detector components (the muon is
reconstructed in the spectrometer; the FASERν tungsten is only the target),
arXiv:2412.03186, (65.6 ± 1.4) fb⁻¹: 338.1 ± 19.0 ± 8.8 observed,
background-subtracted.  It is read by `analysis/faser_electronic.py` and the
paper's Fig. 7.2 (`analysis/paper_plots/pp12b_faser_electronic.py`).

A second block, `conf`, with the unfolded yields of a later FASER conference
note (Fig. 7.2, right panel), is digitised from that note's figures and is
not distributed with the public release; without it the figure shows the
predictions alone and the claims that need it are skipped.

## The observable

Both bin in **−L/E_ν**, with L the lepton number (+1 for ν, −1 for ν̄): the
sign is the charge and the magnitude is 1/E, so a one-dimensional per-charge
energy flux predicts it.  Six bins: ν_μ 100–300, 300–600, 600–1000 GeV;
both charges above 1 TeV (the spectrometer cannot sign so stiff a track);
ν̄_μ 300–1000 and 100–300 GeV.

## Transcribed from arXiv:2412.03186

- **Table I**: simulated signal, backgrounds and observed counts per bin
  (50, 97, 71, 69, 48, 27; total 362).
- **Table II**: acceptance and reconstruction efficiency per bin; the muon
  cuts p > 100 GeV and θ < 25 mrad are the benchmark's FASER Tier S.
- **Table III**: the simulated flux and the flux-weighted mean CC
  cross-section per nucleon per bin, which the generators predict directly.
