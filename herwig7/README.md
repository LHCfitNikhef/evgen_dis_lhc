# Herwig 7.3.0 in the benchmark

Herwig supplies POWHEG-matched NLO samples (angular-ordered shower, cluster
hadronisation) on both currents, on a proton and a genuine neutron, at the
five ladder energies, plus the FASER ladders of Sect. 7.  It also showers the
POWHEG-RES / POWHEG-V2 Les Houches files for the shower comparison.

## The production route: `herwig7/production/`

- `make_cards.py` writes every card into `production/cards/` (`--check`
  reports drift): `V2-<cur>-PWG[NEG]-<p|n>[-<E>].in` for the benchmark,
  `FL-<nu|nubar>-PWG[NEG]-<p|n>-<E>.in` for the FASER ladder, and
  `V2VAL-*` neutron-validation arms.  Each card is a template
  (`DIS-{mu,nu}-POWHEG.in` / `-PWGNEG.in`) with a fixed list of changes made
  by asserted substitutions — a pattern that is not found aborts the script.
- `production.sh` runs the benchmark jobs (`--list` prints them);
  `tools/herwig_faser_ladder.sh` runs the FASER ladder.  Both source `lib.sh`,
  which builds one `.run` per card in `production/build/` (untracked) and
  fingerprints card + plugin: a finished job built from a different card is
  refused, never mixed.  Jobs land in
  `herwig7/v2<cur>pwg[neg]_<p|n>_job[_<E>]_N/` (FASER ladder: `fl<cur>pwg...`).
- `validate_neutron.sh`, `inspect_sample.py`: neutron-beam validation and a
  per-job sanity check.
- Analysis: per half, `make_histos_nlo.py` (muon) or `analysis/analyze_nu.py`
  (neutrino); `combine_nlo.py` forms positive − negative; the chain is
  `tools/analyse_production.py`.
- `herwig_xsec.py` is **the one place Herwig's cross-section bookkeeping
  lives**; every reader imports it.  `plugin/` is the `LeptonicDISCut` ThePEG
  plugin (own Makefile, reads `config.mk`).
- `LHE-{mu,nu}.in`: the card for showering POWHEG-RES / POWHEG-V2 Les Houches files
  (`powheg/run_powheg_herwig.sh`, `powheg/production/herwig_arm.sh`).

**Legacy single-run drivers.**  `run_herwig.sh`, `run_herwig_nu.sh`,
`run_herwig_pwg.sh` and the `DIS-*.in` cards produced the earlier LO,
matrix-element-level and NLO studies listed in the caption of the paper's
generator table (Table 2.1).  They remain as provenance and as templates.

Run Herwig in the scrubbed environment the drivers build from `config.sh`,
not from a conda shell.

## Settings

- Beams in the lepton–nucleon c.m. frame with k·P = E·M_P for either nucleon.
  The neutron is `/Herwig/Particles/n0` **given the proton set**:
  ThePEG::LHAPDF isospin-conjugates the set itself for 2112, so a neutron set
  would swap twice and silently generate a proton.
- NNPDF40_nnlo_as_01180 in the hard process and backward evolution;
  μ_F = μ_R = Q; `EW/Scheme Independent` with α_em(M_Z) = 1/132.119 and
  sin²θ_W = 0.22291; massless quarks; QED radiation and MPI off;
  `MaxLifeTime 10 mm`; HepMC3 output (`GenEventHepMC3`).
- Generation cuts are deliberately looser than the region (MinQ2 3.5,
  MinW2 5, y in (0, 1)); the exact Q² > 4, W > 3 GeV region is applied in the
  analysis.  Herwig moves the scattered lepton in the shower, so a cut
  coincident with the region leaves a silent hole at its edge.
- ThePEG event-error retry (`patches/thepeg-2.3.0-event-error-retry.diff`,
  `MaxEventErrorRetries 100`) and `DISRemnantOption NoLepton`: without them
  Herwig discards events whose remnant cannot be put on shell (mostly
  high-x charm and sea) or moves the lepton to absorb the recoil, depleting
  and distorting the charm channel while the cross-section looks fine.  A
  card with the retry switch cannot be read by an unpatched ThePEG.

## Known traps

1. **SimpleDISCut defaults to `MaxQ2 = 100` GeV²**, removing the high-Q² tail
   at an almost unchanged cross-section; every card sets 1e6.  Its `Maxy`
   must also be set to 1, or a y window is silently re-imposed.
2. **SimpleDISCut builds y from x_parton, not x_Bj.**  At NLO a generation
   y cut accepts a real emission and rejects its own subtraction term (the
   total once came out negative).  Never cut y at generation.
3. **The neutrino factor 2.**  `MEChargedCurrentDIS` averages over two
   neutrino helicities (`me *= 0.25`) unless the beam is declared polarised
   (`PolarizedBeamParticleData`, helicity −1), so an unpolarised CC beam gives
   half the physical cross-section.  The branch exists only in Herwig's native
   POWHEG-matching path, not in Matchbox.
4. **The end-of-run cross-section table has more than one line**; reading the
   one scaled by the shower veto gave a spurious 1.4% deficit.  Hence
   `herwig_xsec.py`.
5. **`Contribution 1` is only the positive half of the NLO cross-section.**
   The negative-weight half (`PWGNEG` cards) is generated and subtracted.
6. **`IntrinsicPtGaussian` defaults to 2.2 GeV** (an LHC value); at these
   energies it makes the remnant handler discard low-Q² events, which looks
   like a Q² tilt.  All cards set 0.5 GeV.
7. **Build traps.**  `Herwig run` was killed with "Code Signature Invalid"
   until `ThePEG/HepMCAnalysis.30.so` was re-signed (`codesign -f -s -`).
   `Format GenEvent` writes HepMC2 even against HepMC3.  The NC POWHEG-matching setup
   needs `Herwig-cache/` copied alongside the `.run`.
8. **A card built by string substitution can silently stay LO**, which is why
   `make_cards.py` asserts every substitution.
