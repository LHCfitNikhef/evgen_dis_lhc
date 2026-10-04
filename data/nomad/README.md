# NOMAD dimuon data

`nomad_table6_enu.txt` is Table 6 of the NOMAD dimuon measurement (Samoylov
et al., *Nucl. Phys.* B876 (2013) 339, arXiv:1308.4750): the ratio of the
dimuon to the inclusive charged-current cross-section, binned in
reconstructed neutrino energy, for ν_μ on an **iron** target.  Read by
`analysis/nomad_dimuon.py`.

## Columns

Whitespace-separated, 19 rows, no header:

| # | quantity |
|---|---|
| 1 | bin lower edge, E_ν [GeV] |
| 2 | bin upper edge, E_ν [GeV] |
| 3 | bin centre, E_ν [GeV] — the energy the theory is evaluated at |
| 4 | R = σ_μμ / σ_CC, in units of 10⁻³ |
| 5 | statistical uncertainty on R, 10⁻³ |
| 6 | systematic uncertainty on R, 10⁻³ |
| 7 | statistical uncertainty, % |
| 8 | systematic uncertainty, % |

## Things that are easy to get wrong

- **The observable is a ratio, in units of 10⁻³** (2.807 to 8.859).
- **No covariance matrix is published** with this table; a χ² from this file
  adds statistical and systematic in quadrature and treats points as
  uncorrelated.  That approximation is recorded in the result file.
- **The tables in x and in √ŝ cannot be combined with this one**: their
  mutual correlations were never published, so only one may be used.  This is
  the one binned in a variable the experiment measures directly.
- **The target is iron.**  `nomad_dimuon.py` builds iron from free-nucleon
  PDFs (no nuclear corrections), as the measurement's own theory did.  A free
  proton roughly doubles the charm fraction, because the CC denominator on a
  proton is driven by d_v.
- **"Dimuon" is defined by an effective branching ratio, not by a cut**:
  B(E) = a / (1 + b/E) with a = 0.097, b = 6.7, from the NOMAD paper.  Its
  energy dependence is the second-muon acceptance, so a Monte Carlo that
  decays charm and applies a muon cut must not also multiply by it.
