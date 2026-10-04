# Regenerating the shower-arm samples (QED and Herwig)

The 1 TeV shower arms re-shower a subset of the production Les Houches files
with one thing changed: QED radiation switched on in stages (`fsr`, `fsrisr`
for the muon only, `full`; overlays in `powheg/qed/`), or Herwig 7 in place of
Pythia 8.  The muon arms re-shower POWHEG-RES, the neutrino arms POWHEG-V2,
on p and n.  Each arm is compared with its own baseline jobs (the same job
numbers), never with the whole sample.

## What was deleted, and what is kept

The 92 arm `events.hepmc` files, 47.6 GB, listed per file in
`regenerate/arms.csv` (columns `sample`, `file`, `events`, `showered_from`,
`GB`; `showered_from` is the LHE, under `$POWHEG_RES` or `$POWHEG_V2`).  They
were deleted first in mid-September 2026 once the QED and shower-model
figures were made, regenerated (re-showered) when the POWHEG-V2 samples were
regenerated with `ubexcess_correct 1` and for the FASERnu track-threshold
re-analysis, and deleted again on 2026-09-20.  The CSV lists the first
deletion; the second covered the same 92 directories with the same event
counts to within the regeneration.

Every result they feed is tracked and was written before the deletion:
`histos_powheg[_nu]_{sub,qedfsr,qedfsrisr,qedfull,hw}_{q4w3,q4w3_faser_e}_{p,n,W}.json`
in `results/` and `results_nu/` (the QED and shower-model figures and the
FASERnu-selection arms).  Each job directory keeps `shower.cmnd` or
`shower.in`, `shower.log`, `events_xsec.json`, `V2_OK`/`HW_OK`,
`delivery.json` and the Herwig `.out` table, and the production LHE files
are kept, so every arm regenerates by re-showering alone.

## How to regenerate

    powheg/production/qed_arms.sh <mu|nu> <p|n> <fsr|fsrisr|full> [njobs] [max_parallel]
    powheg/production/herwig_arm.sh <mu|nu> <p|n> [njobs] [max_parallel]

`njobs` defaults to the arm subset: the 10 lowest-numbered complete muon
jobs (about 500k events; on the proton, jobs 1-5 and 7-11) and the first 2
neutrino batches (about 520k events) per nucleon; `analysis/sample_layout.py:arm_jobnums`
reads it back.  The QED arms reuse the production card with the overlay
appended; Pythia 8 sets no seed after `main_powheg`, so a re-shower
reproduces the deleted events exactly (`qed_arms.sh closure <mu|nu> <p|n>`
checks this on one job).  The Herwig arm uses `herwig7/LHE-{mu,nu}.in` with
seed 7000 + N and reproduces the results to their statistical error.

Then re-analyse the arm keys in `q4w3` and `q4w3_faser_e`
(`analysis/analyze.py v2_powheg_<arm>`, `analysis/analyze_nu.py v2_powheg_nu_<arm>`,
then `analysis/combine_target.py`), and delete the event files again once
their results are written.
