# Regenerating the event files away from 1 TeV

The production samples at 300, 400, 700, 2000 and 4000 GeV (all generators,
both currents, p and n), and POWHEG-V2's 20, 50, 100 and 200 GeV points, are
deleted once analysed.  1 TeV is kept.  The per-file record of every deletion
is `regenerate/energies.csv`:

| column | meaning |
|---|---|
| `deleted` | date of the deletion |
| `group` | generator block of the first deletion, or the pion-pass prune of `tools/prune_energies.py` |
| `file` | the deleted event file (`$SHERPA_RUNS` as in `config.sh`) |
| `events` | events in the file, where recorded (empty for Herwig's first deletion and for Sherpa in the prunes) |
| `GB` | size |

A file appears once per deletion: the samples were regenerated for the
single-inclusive hadron passes and deleted again, so the CSV holds 1435 rows
for 495 distinct files.

## What they feed, and what is kept

Every result read from them is tracked: the integrated and charm results
`histos_<gen>[_charmfinal]_{q4w3,q4w5}_{p,n,W}_<TAG>.json` (the
cross-section-versus-energy figures, inclusive and charm) and the hadron spectra
`results{,_nu}/pions_*_q4w3_*.json` (the hadron-yield study).  Only the event
files went (`events.hepmc`, Sherpa `evtfull`).  Every job directory keeps its
sidecars (`events_xsec.json`, `V2_OK`, Herwig `.out` tables, Sherpa
`sherpa.log` and `Results.zip`), and the POWHEG-RES and POWHEG-V2 Les Houches files
and integration grids are kept.  So those regenerate by re-showering alone,
Sherpa by generation without re-integration, Herwig and GENIE from scratch.

The GENIE G18_02a neutrino samples (`nu_grv` at 300/400/700 GeV, `nu_nnpdf` at
400/700 GeV) were regenerated with the charm-propagator fix and are on disk
again; their rows in the CSV record the earlier deletions.

## How to regenerate

    tools/regenerate_energies.sh powheg|herwig|genie|sherpa|all|sidis-nu

which calls the production drivers of `PRODUCTION.md` with their own seeds
(the seed is a function of code, current, nucleon, energy and job), so the
regenerated samples are statistically equivalent to the deleted ones and,
for the POWHEG-RES and POWHEG-V2 re-showers, identical.  Each driver rebuilds only what is
missing.  `sidis-nu` brings back the neutrino samples of the whole hadron
ladder, POWHEG-V2's low points included.

| generator | driver | per point, per nucleon |
|---|---|---|
| POWHEG-RES (mu) | `powheg/production/production.sh` (re-shower, `run_res.sh`) | 200k showered, 25k per seed |
| POWHEG-V2, POWHEG-V2mc (nu) | `powheg/production/run_v2.sh`, `run_v2_mc.sh` (re-shower, matchInOut off) | 2 batches of 105k |
| Herwig 7 | `herwig7/production/production.sh` | 4 x 20k positive + 1 x 20k negative |
| Sherpa | `tools/sherpa_production.sh run` | 2 x 25k |
| GENIE | `genie/production/genie_job.sh run <cfg> <t> <E> 25000 8` | 8 x 25k (jobs 1-4 feed the integrated results, 5-8 were added for the hadron study) |

**POWHEG-V2 below 300 GeV** is not covered by the `powheg` lane: run
`powheg/production/run_v2.sh <E> <p|n>` for E = 20, 50, 100, 200.  The 50 GeV
proton point ran 7 batches (`V2_TARGET=630000`) to dilute a single
large-weight event in batch 1; pass the same target, or the spike returns at
2.5% of the weight (`DECLARED_DEFICITS` in `analysis/analyze.py`).

## Then

Do not re-run `tools/analyse_production.py` over the tracked integrated
results: they were analysed from the original samples, and a regenerated
sample reproduces them only to their statistical error.  A new per-event pass
reads the regenerated events (for the hadron study,
`tools/faser_pions_passes.sh all` with `ES="..."` and `NU_KEYS`/`MU_KEYS`), and
`tools/prune_energies.py` deletes them again once its results are written.

`tools/prune_energies.py` appends a Markdown table of what it deleted to the
end of this file; move such a table into `regenerate/energies.csv`.
