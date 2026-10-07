# Event Generation for Neutrino and Muon DIS at the LHC

Code, generator configurations and results of a benchmark of Monte Carlo event
generators for deep-inelastic scattering at TeV energies: neutral-current muon
DIS (μ⁻ N → μ⁻ X) and charged-current neutrino DIS (ν_μ N → μ⁻ X), under one
common specification, compared with analytic structure-function calculations
and applied to FASER.

This repository accompanies the paper

> **Event Generation for Neutrino and Muon DIS at the LHC**
> Paulina Hernandez-Sainz, Felix Kling, Jelle Koorn and Juan Rojo
> arXiv:XXXX.XXXXX (to appear)

whose LaTeX source is in `paper/`; every number and figure in it is derived
from the result files tracked here.

**If you use this code, its generator configurations or its results, please
cite the paper above:**

```bibtex
@article{HernandezSainz:2026evgen,
    author  = "Hernandez-Sainz, Paulina and Kling, Felix and Koorn, Jelle and Rojo, Juan",
    title   = "{Event Generation for Neutrino and Muon DIS at the LHC}",
    eprint  = "XXXX.XXXXX",
    archivePrefix = "arXiv",
    primaryClass  = "hep-ph",
    year    = "2026"
}
```

The arXiv number will be added here on submission.  Please also cite the
event generators and tools you use (see the paper's references).

The code is released under the MIT licence (`LICENSE`).  External inputs in
`data/` (fluxes, measurements) remain under the terms of their sources, given
in each directory's README.

## The benchmark specification

| item | value |
|---|---|
| processes | μ⁻ NC DIS (γ*/Z exchange) and ν_μ / ν̄_μ CC DIS |
| beam energies | 400, 700, 1000, 2000, 4000 GeV (target rest frame); FASER ladders from 30 GeV to 6.8 TeV for the applications |
| target | tungsten per nucleon, as (74 p + 110 n)/184 from separate proton and genuine-neutron samples |
| fiducial region | Q² > 4 GeV² and W > 3 GeV, no cut on y (W > 5 GeV for some charm results) |
| PDF | NNPDF4.0 NNLO (`NNPDF40_nnlo_as_01180`), isospin-mirrored for the neutron |
| scales | μ_F = μ_R = Q |
| electroweak | Gμ scheme with fixed α = 1/132.119, sin²θ_W = 0.22291 |
| heavy quarks | massless charm and bottom (ZM) in the generators |
| QED, MPI | off (QED switched on only in the QED study) |
| stable particles | cτ > 10 mm |

## Generators and references

| code | version | role |
|---|---|---|
| Sherpa | 3.0.5 | MC@NLO, CS shower, AHADIC cluster hadronisation |
| Herwig | 7.3.0 (ThePEG 2.3.0) | POWHEG-matched NLO, angular-ordered shower, cluster hadronisation |
| POWHEG-RES (`DIS_v`, POWHEG-BOX-RES) | — | muon NC NLO, showered with Pythia 8.311 or Herwig |
| POWHEG-V2 (`nu-DIS`, POWHEG-BOX-V2) | — | neutrino CC NLO (massless; "V2mc" with massive charm) |
| GENIE | 3.06.02 | tunes G18_02a ("GRV98LO", FASER's default), GHE19_00a ("HEDIS"), and G18_02a with NNPDF4.0; hadronisation through Pythia 8 |
| Pythia | 8.311 | LO DIS, and the shower of the POWHEG samples |
| MG5_aMC | 3.7.2 | LO, at parton level and showered by Pythia 8 |
| NuWro, GiBUU | 25.11.1, 2017 | comparison of neutrino generators (appendix) |
| YADISM | 0.13.11 (EKO 0.15.5) | structure functions at LO/NLO/NNLO in ZM-VFNS and FONLL |

Two POWHEG codes are used and are always named: **POWHEG-RES** for muon NC
and **POWHEG-V2** for neutrino CC.  YADISM results always state their scheme,
**ZM-VFNS** (inclusive) or **FONLL** (charm).

## Layout

| path | content |
|---|---|
| `analysis/` | HepMC3 analysis (`analyze.py` muon, `analyze_nu.py` neutrino), selections, YADISM references, FASER applications |
| `analysis/paper_plots/` | one script per paper figure — see its README |
| `results/`, `results_nu/` | the tracked results (muon / neutrino), `histos_<generator>_<region>_<p\|n\|W>[_<E>].json` |
| `sherpa/`, `herwig7/`, `powheg/`, `genie/`, `pythia8/`, `mg5/`, `nuwro/`, `gibuu/` | cards and drivers per generator, each with a README |
| `tools/` | production drivers, ladders, consistency checks |
| `patches/` | local modifications of the external codes — see its README |
| `data/` | external inputs (FASER fluxes, the FASER electronic-detector measurement, NOMAD, DONUT), each with its provenance |
| `data/pdfs/lhapdf/` | derived PDF sets (neutron mirrors, tungsten averages) |
| `PRODUCTION.md`, `REGENERATE_*.md` | how the event samples were produced, and how to regenerate deleted ones from their seeds |
| `paper/` | the LaTeX source of the paper and its figures |
| `CONVENTIONS.md` | the working rules the code follows (cited in comments as "rule N") |

## Installing

The repository holds everything written for this work.  The generators and
libraries are installed separately under one external tree, `$PHYSICS24`
(default `~/physics24`, codes in `$PHYSICS24/software/`); the repository itself
can live anywhere.  Every external location is named once, in `config.sh`;
shell scripts source it, Makefiles include `config.mk`, Python does
`import paths`.  No tracked file contains an absolute path.

**1. Python environment and LaTeX.**

    conda env create -f environment.yml && conda activate evgen-bench

`environment.yml` pins numpy, scipy, matplotlib, ROOT, LHAPDF, HepMC3 and
YADISM/EKO/PineAPPL at the versions the results were produced with.  Keep the
environment active when running the drivers (a few small helpers call
`python3`).  The figures use LaTeX (`usetex`), so a TeX installation with
`type1cm` and `dvipng` must be on `PATH`.

**2. Machine settings go in `config.local.sh`** (untracked, sourced by
`config.sh`; never edit `config.sh`).  With the conda environment above:

    BENCH_DEPS=$CONDA_PREFIX          # ROOT, LHAPDF, HepMC3 and the analysis Python
    BENCH_LOG4CPP=$CONDA_PREFIX       # log4cpp for GENIE

On a Linux cluster an LCG view can supply ROOT, LHAPDF and HepMC3 instead:

    PHYSICS24=/data/faser/physics24
    BENCH_DEPS=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt
    BENCH_PYTHON=$HOME/miniforge3/envs/evgen-bench/bin/python3

All knobs are spelled `BENCH_*` (listed at the top of `config.sh`) because
conda already exports generic names such as `PYTHIA8DATA` and `CXX`, which
point at a different Pythia version and an unsuitable compiler.  Platform
differences (macOS/Linux) are decided once in `config.sh`
(`BENCH_CXX`, `BENCH_NOSLEEP`, `BENCH_LIBPATH_VAR`); scripts target bash 3.2.

**3. PDF sets.**  Install the base sets with LHAPDF
(`lhapdf install <set>`): `NNPDF40_nnlo_as_01180`, and for the PDF and
nuclear-PDF studies `CT18NNLO`, `MSHT20nnlo_as118`, `ABMP16_5_nnlo`,
`ATLASpdf21_T1`, `GRV98lo`, `CT18ANLO`, `NNPDF31sx_nlo_as_0118_LHCb_nf_6`,
`nNNPDF30_nlo_as_0118_A184_Z74` and `nCTEQ15HQ_FullNuc_184_74`.  The derived
sets in `data/pdfs/lhapdf/` are partly tracked (their `.info` files and the
central members) and are rebuilt in about a minute:

    tools/make_neutron_pdf.py <SET> --all-members      # isospin-mirrored neutron sets
    tools/make_isoscalar_pdf.py NNPDF40_nnlo_as_01180 --index 341000 --all-members

**4. Check the setup:** `tools/check_paths.sh` reports every configured path
that does not exist, every Python module the analysis interpreter cannot
import, and a missing `latex`.

**5. Rebuild the paper figures and run the checks** — no generator needed,
since every result is tracked:

    for f in analysis/paper_plots/pp*.py; do "$BENCH_PYTHON" "$f"; done
    tools/run_checks.sh --fast

`tools/run_checks.sh` re-derives every number quoted in the paper and in the
paper-plot scripts from the result files, and checks naming and path
conventions.  Checks that need event
samples skip rather than fail.

Three inputs are not distributed: the FASERν emulsion and electronic-detector
data digitised from FASER conference notes, and the NNLO SIDIS curves of
arXiv:2504.05376 provided by its authors.  Without them Figs. 7.1 and 7.2 are
drawn with the predictions only, the claims that compare with those data are
reported as skipped; the predictions themselves are all included.

**6. Regenerating event samples** needs the generators, built from source with
the diffs in `patches/` applied (each patch gives its own `patch` line).  The
productions and their drivers are described in `PRODUCTION.md`; how to
regenerate deleted samples is in `REGENERATE_*.md`, with the per-file lists
in `regenerate/*.csv`.  Drivers skip
work whose output exists, so they can be restarted, and refuse to start when
the disk holding `$PHYSICS24` is short of space.  Sherpa's genuine-neutron
samples need a patched second install (`tools/build_sherpa_npatch.sh`), whose
build script currently supports macOS only.

Generator cards are kept here only: the run directories in the install trees
hold symlinks back into the repository (`tools/link_cards.sh --check` reports
missing or diverged links, `tools/link_cards.sh` repairs them).

## Conventions in the analysis

- Every observable is computed from frame-invariant quantities in the target
  rest frame (E = p·P/m_N, angles from p_T and E), because generators write
  their events in different frames (lab, or the lepton–nucleon c.m.).
- Weighted samples pass a delivery closure test (delivered vs integrated
  cross-section, refused beyond 2%), and every result records the input files
  it was built from (`analyze.py: closure_check, input_manifest`).
- Jobs are combined with inverse-variance weights for repeated integrals and
  accepted-count weights for pooled event samples.
