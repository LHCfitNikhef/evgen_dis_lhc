# Regenerating the GENIE non-DIS samples

The appendix on non-DIS processes compares, for GENIE's default tune, the full
event-generator list with the non-DIS processes alone, for both currents on
p and n at 100, 400 and 1000 GeV.  Its inputs are the `events.gst.root`
summary ntuples of two sample kinds:

| kind | directory | process list | events per point |
|---|---|---|---|
| full | `genie/v2ndfull_<mu\|nu>_<t>_job[_TAG]_N/` | EM (mu), CC (nu) | 1,000,000 (8 x 125k) |
| non-DIS only | `genie/v2ndonly_<mu\|nu>_<t>_job[_TAG]_N/` | EMNONDIS (mu), CCNONDIS (nu) | 200,000 (8 x 25k) |

with the FASER splines `genie/splines/faser/<pid>_<t>_{em,cc}_e8000.xml`.  The
per-file record (192 files, 14.2 GB) is `regenerate/nondis.csv`, whose
columns repeat each job's `events.meta`: `list`, `nev`, `beam_GeV`, `pid`,
`target`, `spline`, plus `file` and `GB`.

## State

All 192 files were deleted on 2026-09-30.  The **neutrino** samples (96 files,
all three energies) were then regenerated with the GENIE charm-propagator fix
(`patches/genie-aivazis-charm-propagator.diff`) and are on disk; the **muon**
samples (96 files) remain deleted.

Every result read from them is tracked:
`results{,_nu}/genie_nondis_diff_{p,n,W}[_<TAG>].json` from
`analysis/genie_nondis_diff.py`, which is all the appendix figures
(`analysis/paper_plots/ppA1_genie_nondis.py`, `ppA1b_genie_nondis_100GeV.py`)
read.  Each job directory keeps its sidecars (`events.meta`,
`events.ghep.status`, `rejections.json`, `genie.log.tail`).

## How to regenerate

    genie/run_genie_nondis.sh <mu|nu> <E> 1000000 200000 8 <p|n>

for E = 1000, 400, 100, then `analysis/genie_nondis_diff.py <mu|nu> <E>` for
each current and energy, with a Python that has ROOT (`$BENCH_DEPS/bin/python3`).  GENIE takes minutes per
point.  A regenerated sample reproduces the results to their statistical
error, not event by event.
