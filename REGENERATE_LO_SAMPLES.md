# The LO event samples (a record)

The LO generators of the paper's generator table (Pythia 8 and MG5_aMC at
LO, Sherpa at LO and at fixed-order NLO, Herwig at LO) were run in the first
stage of the benchmark, which used a different setup from the paper's:
a free proton target, 1 TeV (and 400 GeV / 4 TeV), and the region
Q2 > 4 GeV2 with **0.2 < y < 0.9**, imposed at generation by Sherpa and
Herwig.  Their event files are deleted (Pythia 8 and Sherpa LO on 2026-09-07,
about 46 GB; Herwig LO and the showered MG5_aMC events with the rest of that
production), and **their results were deleted on 2026-09-19**: no LO result
is tracked, shown in the report, or used in a paper figure.  Git history
keeps the deleted results.

Because Sherpa and Herwig cut the y window at generation, their samples
could not be re-analysed in the benchmark region even if regenerated; a LO
comparison in that region needs new samples generated without the y cut.  This file records
how the original ones were made.

## Pythia 8 LO (`pythia8/`)

| job directories | events per job | seeds |
|---|---|---|
| `job_1`, `job_2` (mu NC, 1 TeV) | 250 000 | 4322, 4323 |
| `job_400GeV_1`, `_2` | 250 000 | 4722, 4723 |
| `job_4TeV_1`, `_2` | 250 000 | 8322, 8323 |
| `nu_job_1` .. `nu_job_4` (nu CC, 1 TeV) | 150 000 | 7322 .. 7325 |
| `nu_job_400GeV_1`, `_2` | 250 000 | 7722, 7723 |
| `nu_job_4TeV_1`, `_2` | 250 000 | 11322, 11323 |
| `shw{dire,vincia,simple,norecoil}_job_1` (mu NC, 1 TeV, shower models) | 250 000 | 4322 |

Seed = SEED0 + i, SEED0 = 4321 (mu) or 7321 (nu) at 1 TeV and SEED0 + E/GeV
elsewhere.  Commands, from `pythia8/`:

    ./run_pythia.sh    250000 2 [400|4000]       # mu NC
    ./run_pythia_nu.sh 150000 4                  # nu CC, 1 TeV (250000 2 at 400 / 4000)
    ./run_pythia.sh --shower <dire|vincia|simple|norecoil> 250000 1

`run_pythia_me.sh` is the matrix-element-level variant.  Each job directory
still holds its `pythia.log` and `events_xsec.json`.

## Sherpa LO (`$SHERPA_RUNS/<run>/job_N/evtfull`)

| run directory | jobs | events per job | seeds |
|---|---|---|---|
| `MuonDIS_LO` | 1 | 500 000 | 2001 |
| `MuonDIS_LO_400GeV`, `MuonDIS_LO_4TeV` | 2 | 250 000 | 1235, 1236 |
| `NuDIS_LO` | 8 | 25 000 | 1235 .. 1242 |
| `NuDIS_LO_400GeV`, `NuDIS_LO_4TeV` | 2 | 250 000 | 1235, 1236 |
| `NuDIS_LO_Dire` (shower model) | 8 | 25 000 | 1235 .. 1242 |

`sherpa/Runs/<run>/run_parallel.sh <events_per_job> <njobs>` (seed 1234 + i;
linked into the install tree by `tools/link_cards.sh`) and
`sherpa/run_sherpa_shower.sh nu lo Dire 25000 8`.  The integration grids
(`Results.zip`) and process libraries are kept in the run directories, so
generation needs no re-integration.  `RESPECT_MASSIVE_FLAG: true` must stay in
the card, or 9-23% of events are silently discarded.

**Sherpa fixed-order NLO** is parton level (`NLO_Mode: Fixed_Order`, shower
off), in `sherpa/Runs/MuonDIS_NLO_ME/`, with the LO counterparts
`MuonDIS_LO_ME` and `NuDIS_LO_ME`.

## Herwig 7 LO and MG5_aMC LO

Herwig LO: `herwig7/run_herwig.sh` and `run_herwig_nu.sh <events_per_job> <njobs> [E]`
with the cards `herwig7/DIS-mu.in` and `DIS-nu-POL.in` (and their 400 GeV /
4 TeV variants).  MG5_aMC LO, at parton level and showered by Pythia 8:
`mg5/run_mg5.sh --current <mu|nu> [--energy E] <nevents>`, with the DIS cut
and scale hooks of `mg5/dis_hooks.f`.  Their event counts are not recorded
here.

Regenerating any of these reproduces cross-sections to their statistical
error, not the events.
