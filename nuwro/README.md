# NuWro in the benchmark

NuWro (Wrocław) is a dedicated neutrino–nucleus event generator.  It enters
only the appendix comparison of GENIE with NuWro and GiBUU: ν_μ CC DIS at
E_ν = 200 GeV, in the benchmark region, on a free proton and on tungsten.

- Upstream: <https://github.com/NuWro/nuwro>, built version `NuWro_25_11_1`,
  installed at `$NUWRO_DIR` (`config.sh`: `BENCH_NUWRO_DIR`).  Only the
  driver, the converter and this file are ours.

## Building it

Conda's ROOT ships `libEGPythia8` but not `libEGPythia6`, and NuWro's DIS
fragmentation is `TPythia6`.  The plugin is built separately from
<https://github.com/luketpickering/ROOTEGPythia6> **with its own builtin
Pythia 6** (`BENCH_ROOTEGPYTHIA6`) — deliberately not GENIE's Pythia 6, so
the two generators do not share a library version.  Then `make` in the NuWro
tree.

## Running it

    make -C nuwro                                         # converter, once
    nuwro/run_nuwro.sh --current nu --target p 200 100000
    nuwro/run_nuwro.sh --current nu --target W 200 100000

writes `nuwro/nuwro_job_<cur>_<target>[_<E>]_N/events.hepmc`, read by
`analysis/analyze_nu.py` unchanged.

- **DIS only.**  Quasi-elastic, resonant, coherent, meson-exchange and hyperon
  channels are switched off: this is a DIS comparison, and non-DIS channels
  are per mille of the fiducial rate (non-DIS appendix).
- **Targets.**  `p`, `n`: a free nucleon, no Pauli blocking, no cascade — like
  for like with the other generators.  `W`: tungsten (74 p + 110 n), local
  Fermi gas, Pauli blocking and intranuclear cascade (`BENCH_NUWRO_FSI=0`
  keeps the nucleus and switches the cascade off).
- **The cascade at TeV energies.**  At 200 GeV the full cascade completes.
  At 1 TeV on tungsten it aborts on a negative cross-section
  (`kaskada7.cc`) about once in a few tens of thousands of events.

## The converter (`nuwro2hepmc.cc`)

Reads NuWro's ROOT tree and writes HepMC3 through HepMC3's own writer.

- **Beams**: the neutrino and the struck nucleon, status 4.  On tungsten the
  struck nucleon is bound and carries Fermi momentum, so the DIS invariants
  are per struck nucleon.
- **Final state**: the particles leaving the nucleus (after FSI).
- **Weight**: NuWro events are unweighted and each carries the total
  cross-section; the converter writes σ/N per event (and checks that every
  event carries the same σ), so summed weights give σ_fid.
- **D mesons are not carried over**, so NuWro has no charm observable here.

## Stability convention

The benchmark calls a particle stable if cτ > 10 mm.  A scan of NuWro's final
states finds no ρ, ω, η or Δ (Pythia 6 decays them).  Two differences remain,
neither affecting any observable used here: π⁰ is left undecayed (neutral,
same energy in `Ehad`), and neutral kaons appear as K⁰/K̄⁰ (311/−311) rather
than K_S/K_L, which are both stable under our convention anyway.
