# Sherpa 3.0.5 in the benchmark

Sherpa supplies the MC@NLO samples (Catani–Seymour shower, AHADIC cluster
hadronisation) on both currents, muon NC and neutrino CC, on a proton and a
genuine neutron, at the five ladder energies (400, 700, 1000, 2000, 4000 GeV).
The LO, fixed-order and matrix-element-level runs and the Dire-shower twin are
the additional studies listed in the paper's generator table.

## Where things are

- **Cards**: `sherpa/Runs/<name>/Sherpa.yaml`.  `tools/link_cards.sh` symlinks
  each into the install-side run directory `$SHERPA_RUNS/<name>/`, which holds
  what Sherpa writes (integration grids, `Process/`, `Results.zip`, event
  samples).  Those are not tracked.
- **Production cards** `V2_{MuonDIS,NuDIS}_NLO_<p|n>[_<E>]` are written by
  `tools/make_sherpa_cards.py` (`--check` reports drift) from the templates
  `MuonDIS_NLO` and `NuDIS_NLO_ckm3`, changing only beams, PDF and selectors.
  `tools/sherpa_production.sh run` generates them; `SHOWER=Dire` runs the
  Dire twin (each twin integrates itself: the MC@NLO subtraction depends on
  the shower kernels).
- **FASER ladders** (Sect. 7): `tools/sherpa_faser_ladder.sh`
  (`NuDIS_NLO_FASER_*` run directories) and `tools/sherpa_nubar_sidis.sh`.
- **Earlier studies**: `MuonDIS_*`, `NuDIS_*` (LO, `_ME` matrix-element level,
  `_400GeV`/`_4TeV`), `run_sherpa_shower.sh` (CSS vs Dire) and
  `run_sherpa_massive_ps.sh` (massive-shower diagnostic).  Their event files
  are deleted; `REGENERATE_LO_SAMPLES.md` has the recipe and seeds.
- Analysis: `analysis/analyze.py` (muon) and `analysis/analyze_nu.py`
  (neutrino), chained by `tools/analyse_production.py`.

## Settings that matter

- Generated in the lepton–nucleon c.m. frame (Sherpa's `Fixed_Target` beam
  mode is broken in 3.0.5); the analysis is frame-invariant.  Beam energies
  are chosen so that k·P = E·M_P exactly for either nucleon.
- Generation cuts: Q² > 4 GeV² and an `INEL` y floor at the W = 3 GeV edge
  (Sherpa has no W selector); no upper y cut.  The exact W > 3 GeV cut is
  applied in the analysis.  Sherpa's `INEL` y is built from the lepton
  momenta, so it is looser than, not equal to, the invariant y.
- EW: Gμ scheme with fixed α_em = 1/132.119; NLO virtuals from the internal
  `DY_QCD_Virtual` (DIS crossing); Amegic BVI + Comix RS.  Massless charm.
- Weighted events: the delivered cross-section is `sum(Weight)/sum(NTrials)`,
  never the mean weight.

## Known traps

1. **`RESPECT_MASSIVE_FLAG: true` is required.**  Without it the ME-to-shower
   interface silently discards 9% of muon and 23% of neutrino events and guts
   the charm channel, while the integrator's cross-section looks fine.  With
   it the shower is massless too, so g → cc̄ in the shower makes charm at
   12–16 times a massive shower's rate; the Sherpa charm numbers are quoted
   with that shower charm subtracted (`analyze.py`, `subtract_shower_charm`).
2. Sherpa 3.0.5's massive-charm options remove charm from the jet container,
   which is one reason the benchmark is massless (ZM).
3. **CC at NLO needs the full-CKM virtual**, which the stock install lacks:
   `patches/sherpa-dy-qcd-virtual-ckm.diff` (card `CKM: Order: 3`).
4. **The neutrino factor 2.**  Sherpa averages over two helicities of the
   incoming neutrino, so neutrino-beam cross-sections come out at half the
   physical value; `analysis/analyze_nu.py` applies
   `SHERPA_NU_SPIN_FACTOR = 2`.
5. **AHADIC and low W.**  AHADIC stalls on the very-low-W events a Q²-only cut
   allows at x → 1; the W > 3 GeV region avoids them.
6. **Neutron beams.**  Stock 3.0.5 segfaults on the first event of a 2112
   beam (`Hadron_Remnant::RemnantFlavour` builds an invalid diquark for u d d);
   the one-line fix is `patches/sherpa-hadron-remnant-neutron.diff`.  It is
   built into a separate copy at `$SHERPA_NPATCH` by
   `tools/build_sherpa_npatch.sh`, which copies the install, rebuilds
   `libRemnants`, rewrites every rpath (otherwise the copy silently loads the
   unpatched library) and re-signs.  **The build script is macOS-only**
   (`install_name_tool`, `codesign`); on Linux the equivalent is
   `patchelf --set-rpath`.  `tools/sherpa_production.sh` runs every neutron
   directory from the copy and checks, per job, the library actually mapped
   and the neutron valence content of every event.
7. A 2112 beam needs the neutron PDF set (`NNPDF40_nnlo_as_01180_n`, in
   `data/pdfs/lhapdf/`), because Sherpa's LHAPDF interface does no isospin
   swap, and `PDF_LIBRARY` must be named, because the default is set only for
   2212.
