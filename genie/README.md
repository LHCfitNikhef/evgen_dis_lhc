# GENIE in the benchmark

GENIE is the default event generator of FASER, so it is shown in every
benchmark figure.  Three configurations, with these figure labels:

| label | tune / list | overlay (`GXMLPATH`) |
|---|---|---|
| GENIE (GRV98LO) | `G18_02a_00_000` — Bodek–Yang on GRV98 LO, FASER's default | `p8` |
| GENIE (NNPDF4.0) | `G18_02a_00_000` with GRV98 replaced by NNPDF4.0 NNLO | `nnpdf` |
| GENIE (HEDIS) | `GHE19_00a_00_000` — HEDIS, BGR18 NLO structure functions (via APFEL) | `p8` |

Muon NC uses the `EMDIS` list (photon exchange only, no Z); neutrino CC uses
`CCDISCHARM` (DIS-CC + DIS-CC-CHARM) for G18_02a and `CCHEDIS` for HEDIS.  The
build is the GENIE master branch (cited as 3.06.02 in the paper) at
`$GENIE_DIR`, with hadronisation and decays by Pythia 8 through the overlays.

## Running it

    source genie/genie_setup.sh          # environment, from config.sh
    make -C genie                        # gtohepmc3 converter
    genie/production/genie_job.sh spline <cfg> <p|n> <emax>
    genie/production/genie_job.sh run    <cfg> <p|n> <E> <nev_per_job> <jobs>
    genie/production/production.sh <cfg> ...   # the whole ladder, one lane per cfg

`<cfg>` is `mu_grv`, `mu_nnpdf`, `nu_grv`, `nu_nnpdf`, `nu_hedis` (and
`nubar_grv`, `nubar_hedis`).  Proton and neutron are separate free-nucleon
runs with their own splines (`genie/splines/v2/`) and disjoint seeds;
tungsten is combined per nucleon in the analysis.  Each job runs
`gevgen → gtohepmc3 → spline_to_json.py` into
`genie/v2g_<cfg>_<t>_job[_<E>]_N/`, verifies the event count and that the
overlay was read, and deletes the GHEP file and the full log.  The fiducial
region is applied in the analysis (`analyze.py`, `analyze_nu.py`).

The FASER ladders (Sect. 7) use `tools/genie_dimuon_ladder.sh` and
`tools/genie_faser_splines.sh`; the non-DIS appendix uses
`run_genie_nondis.sh`; the neutrino-generator appendix uses
`run_genie_nugen.sh` (free proton and tungsten).  The older single-config
drivers `run_genie{,_nnpdf,_hedis,_nucc_lo,_nnpdf_nu}.sh` are kept as
provenance.

## Config overlays (`genie/config/`)

- `p8/` — Pythia 8 hadroniser/decayer; charm hadrons decayed (GENIE's default
  leaves them stable) and π⁰ decayed, so the ctau > 10 mm convention holds.
- `nnpdf/` — `p8` plus `BYPDF.xml`/`LHAPDF6.xml`/`CommonParam.xml` swapping
  GRV98 LO for NNPDF40_nnlo_as_01180.
- `nu_cc_lo/`, `nondis/` — generated on every run (`make_nucc_lo_overlay.py`,
  `make_nondis_overlay.py`) and not tracked.
- `hienergy/`, `ultrahigh/` — raise the declared validity ceiling; see their
  READMEs.

## Local modifications (`patches/`)

1. `genie-em-q2-floor-and-pythia8-teardown.diff`: the EM Q² floor raised to
   4 GeV² (upstream 0.02), so muon generation and spline share the fiducial
   floor; also the Pythia 8 stability convention.  Every muon EM
   cross-section from this build, non-DIS included, is for Q² > 4 GeV².
2. `genie-aivazis-charm-propagator.diff`: upstream `AivazisCharmPXSecLO`
   multiplies by the W-propagator factor (1 + Q²/M_W²) instead of dividing,
   so charm was too large by (1 + Q²/M_W²)⁴.  Because `QPMDISPXSec` subtracts
   charm from inclusive DIS with `max(0, incl − charm)`, the bug also changed
   the total.  At 1 TeV on a proton the fix moves G18_02a charm by ×0.924 and
   total CC by ×0.998.  HEDIS is unaffected.  All G18_02a neutrino results use
   the patched build.

## Known traps

- **`GXMLPATH` fails open.**  A missing overlay directory falls through to
  `$GENIE/config` and silently runs the built-in tune (Pythia 6, charm
  stable).  `genie_setup.sh` and `genie_job.sh` refuse a missing directory,
  and each job checks the log for "Custom directory".
- **A spline must extend past the beam.**  Evaluating it at its last knot
  returns 0 and `gevgen` hangs in "Could not select interaction".
- **`gmkspl -e` is clamped to the tune's declared validity** (`GVLD-Emax`,
  1000 GeV for G18_02a) with only a WARN line, and the file keeps the name you
  asked for.  `genie_job.sh` checks the last knot of every spline.
  G18_02a is run in the benchmark ladder only up to 1 TeV; HEDIS (valid to
  10¹² GeV) runs at all five energies.  Above 1 TeV (FASER ladders) G18_02a is
  an extrapolation via `config/hienergy`.
- **`CCDIS` alone is charm-free.**  `QPMDISPXSec` subtracts the charm piece,
  so DIS-CC must run together with DIS-CC-CHARM (hence `CCDISCHARM`).
- GENIE enters "as is": its own EW parameters, Bodek–Yang low-Q² corrections
  and AGKY hadronisation.  No shower; for W above 3 GeV AGKY hands the system
  to Pythia string fragmentation.
- For nuclear targets `gtohepmc3` writes the struck nucleon as the beam, so x
  and W are per nucleon.
