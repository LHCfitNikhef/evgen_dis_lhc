# Production of the event samples

How the event samples behind the paper are produced: the benchmark settings,
the generation cuts, the targets, the statistics and seeds, the driver of
each generator, the directory and result-file layout, and the traps a
re-runner has to know.  Paths outside the repository are the variables of
`config.sh` (`$POWHEG_RES`, `$POWHEG_V2`, `$SHERPA_RUNS`, ...).  Deleted event
samples and how to bring them back are recorded in `REGENERATE_ENERGIES.md`
(energies away from 1 TeV, `tools/regenerate_energies.sh`), `REGENERATE_ARMS.md`
(the QED and Herwig shower arms), `REGENERATE_NONDIS.md` (GENIE non-DIS
processes) and `REGENERATE_LO_SAMPLES.md` (the LO studies), with per-file
lists in `regenerate/`.

## Settings

| setting | value |
|---|---|
| beams | fixed-energy mu- (NC) or nu_mu (CC) on a nucleon at rest |
| energies | 400, 700, 1000, 2000, 4000 GeV (`beams.BENCH_ENERGIES`); 1 TeV is the high-statistics point that carries every differential figure |
| target | tungsten by isospin, 74 p + 110 n; each generator runs a proton and a genuine neutron separately |
| normalisation | per nucleon, sigma_W = (74 sigma_p + 110 sigma_n) / 184 (`analysis/combine_target.py`) |
| PDF | `NNPDF40_nnlo_as_01180`, no nuclear modification |
| region (`q4w3`) | Q2 > 4 GeV2 and W > 3 GeV, no y cut, on both currents |
| charm region (`q4w5`) | Q2 > 4 GeV2 and W > 5 GeV (see "Traps", heavy-quark thresholds) |
| everything else | muF = muR = Q; Gmu EW scheme (alpha_em = 1/132.119, sin2thetaW = 0.22291); no QED radiation; massless c and b; stable if ctau > 10 mm |

W is always W2 = M_P^2 + 2 P.q - Q2 and E_lab = P.k / M_P, with the proton
mass, on both nucleons, and every cut is applied in the target rest frame.

## Generation cuts

The exact region is imposed by the analysis on the generated final state, and
sigma_fid = sigma_gen x (signed-weight fiducial fraction).  At generation each
code is cut looser than the region, never tighter, because the shower or the
event record moves the lepton after the hard process:

| generator | generation cut | why looser |
|---|---|---|
| POWHEG-RES | `Qmin 1.5` (Q2 > 2.25), ymin 0, ymax 1 | the LHE writer puts quarks on mass shells and reshuffles the lepton; 2.6% of events generated at Q2 > 4 land below 4 |
| POWHEG-V2, POWHEG-V2mc | `q2cut 2.25` | Pythia re-masses c and b and moves the lepton down in Q2 |
| Herwig 7 | `LeptonicDISCut` MinQ2 3.5, Miny 0, Maxy 1, MinW2 5; `{Neutral,Charged}CurrentCut` MinQ2 3.5, y in [0, 1] | the plugin's W2 drops the hadron mass, and Herwig moves the lepton in 6% (mu) / 2.5% (nu) of events |
| Sherpa | `[Q2, l, 13, 4, 1e12]` and `[INEL, l, 13, y_min, 1]` with y_min = (9 - M_P^2 + 4) / 2kP | Sherpa has no W selector: the y cut is the W > 3 edge at the Q2 floor, looser above it |
| GENIE | none beyond the spline range | |

Sherpa's INEL variable is not exactly the invariant y (it assumes a massless
target), but it is always larger, so the cut never removes an event inside
the region.  Sherpa treats the muon as massless, so its k.P exceeds E M_P by
1e-5 or less.

## Targets: a genuine neutron in every generator

A proton beam carrying a neutron PDF gives the right cross-section but a
proton remnant, one unit of charge off in the hadronic system, which the
hadron-level observables would inherit.  So the neutron run uses a 2112 beam
everywhere, and how the PDF is swapped depends on the code:

| generator | neutron run | PDF given |
|---|---|---|
| POWHEG-RES, POWHEG-V2, POWHEG-V2mc | `ih2 2`; LHE `<init>` beam 2112 | proton set: both codes isospin-swap it themselves, and Pythia 8 swaps `PDF:pSet` for a 2112 beam |
| Herwig 7 | `BeamB /Herwig/Particles/n0` | proton set: ThePEG's LHAPDF interface swaps it for id 2112 |
| Sherpa | 2112 beam, `PDF_LIBRARY: [None, LHAPDFSherpa]` | `NNPDF40_nnlo_as_01180_n` (Sherpa swaps nothing); needs the patched build below |
| GENIE | target 1000000010 | its own nucleon treatment, with the nucleon's own spline |

Giving the `_n` set to a code that swaps would generate a proton.  The
mirrored set `NNPDF40_nnlo_as_01180_n` (`tools/make_neutron_pdf.py`) and the
tungsten average `NNPDF40_nnlo_as_01180_W184free` (`tools/make_isoscalar_pdf.py`)
live in `data/pdfs/lhapdf/`, which the drivers prepend to `LHAPDF_DATA_PATH`.

Sherpa 3.0.5 segfaults on the first neutron event: `Hadron_Remnant::RemnantFlavour`
removes the first valence constituent rather than the struck one and builds a
non-existent diquark for udd.  `patches/sherpa-hadron-remnant-neutron.diff`
fixes it; `tools/build_sherpa_npatch.sh` builds a patched copy at
`$SHERPA_NPATCH` (it relinks `libRemnants` and rewrites the copy's rpaths, or
the copy would silently load the unpatched library), and only neutron jobs run
from it.  Proton events are byte-identical with and without the patch.

Each neutron sample is checked event by event: beams (l, 2112), dd/ud remnant
diquarks and no uu (`herwig7/production/inspect_sample.py`,
`powheg/production/inspect_sample.py`, `tools/sherpa_check_events.py`, which
also confirms the patched `libRemnants` was the one loaded).  The validation
runs (`powheg/production/validate.sh`, `herwig7/production/validate_neutron.sh`)
showed n beam + proton set and p beam + `_n` set to agree to every printed
digit for POWHEG-RES and POWHEG-V2 and within 0.2% for Herwig.

## Statistics per nucleon and current

| generator | 1 TeV | other energies |
|---|---|---|
| POWHEG-RES (mu) + Pythia 8 | 1M showered (50k per seed) | 200k (25k per seed) |
| POWHEG-V2 (nu) + Pythia 8 | 1M showered (4 batches of 260k) | 200k (2 batches of 105k) |
| POWHEG-V2mc (nu, charm only) + Pythia 8 | 1M | 200k |
| Herwig 7 POWHEG NLO | 300k positive + 200k negative (jobs of 50k) | 80k + 20k (jobs of 20k) |
| Sherpa MC@NLO | 200k (8 x 25k) | 50k (2 x 25k) |
| GENIE G18_02a (mu GRV98, mu NNPDF, nu GRV98, nu NNPDF) | 500k (8 x 62.5k) | 100k (4 x 25k) at 400 and 700 GeV; G18_02a is not run above its 1 TeV validity |
| GENIE HEDIS GHE19_00a (nu) | 400k (8 x 50k) | 100k (4 x 25k) at all four |

The POWHEG-RES and POWHEG-V2 targets count showered events, inside and outside the region (for
POWHEG-RES about 47% of them fall inside: the photon pole puts half the events
below Q2 = 4).  A 300 GeV point exists for every generator for the
single-inclusive hadron study (`beams.SIDIS_ENERGIES`), where the GENIE
samples away from 1 TeV run 8 x 25k; POWHEG-V2 also has 20, 50, 100 and
200 GeV (`beams.SIDIS_ENERGIES_LOW`).

## Seeds

Every (code, current, nucleon, energy, job) has its own stream, so the proton
and neutron samples combined into tungsten never share one:

| generator | seed |
|---|---|
| POWHEG-RES | `pwgseeds.dat` line i = base + i |
| POWHEG-V2 | batch j: iseed = base + j; POWHEG-V2mc: the same + 50000 |
| (both) | base = 1000 x (10 x [p 1, n 2] + energy index), index 400, 700, 1000, 2000, 4000, 300, 20, 50, 100, 200 -> 0..9 (`powheg/production/common.sh`) |
| Herwig 7 | 2000000 + 400000 [nu] + 200000 [n] + 100000 [neg] + 1000 x energy index + N, index over 400, 700, 1000, 2000, 4000, 300 (`herwig7/production/lib.sh`) |
| Sherpa | generation `RANDOM_SEED` 1234 + N (p), 5234 + N (n); integration at Sherpa's default seed |
| GENIE | SEED0 + 10 E/GeV + [p 0, n 500] + N; SEED0 = 18000 / 28000 / 38000 / 48000 / 58000 for mu_grv / mu_nnpdf / nu_grv / nu_nnpdf / nu_hedis (`genie/production/genie_job.sh`) |

Pythia 8 sets no `Random:setSeed` after `main_powheg`, so re-showering the same
LHE with the same card is byte-identical.  Regenerating any other sample
reproduces the results to their statistical error, not event by event.

## Drivers and layout

Every driver skips work whose output exists, checks the free disk before each
launch (it keeps a floor of free space on the volume holding `$PHYSICS24`), and can be killed and re-invoked.  `<t>` is
`p` or `n`, `[_TAG]` the energy tag of `beams.at_energy` (`400GeV`, `700GeV`,
`2TeV`, `4TeV`, `300GeV`; 1 TeV is untagged).

**POWHEG-RES and POWHEG-V2 + Pythia 8.**  `powheg/production/make_cards.py`
writes the cards in `powheg/cards/production/POWHEG-{RES,V2,V2mc}/`; the drivers
are `run_res.sh`, `run_v2.sh`, `run_v2_mc.sh`, with the two-lane launchers
`production.sh` and `production_mc.sh` and helpers `common.sh`, `set_keys.py`
(all in `powheg/production/`).  POWHEG-RES integrates on 8 seeds (stages 1-3),
then runs stage 4 seed by seed until the target is met; seeds that abort on an
ISR upper bound are listed in `<rundir>/DEAD_SEEDS` and not showered.  POWHEG-V2
generates batch 1 in `<rundir>` and batch j >= 2 in `<rundir>-b<j>` with batch
1's grids.  Each LHE is NaN-stripped (`powheg/strip_lhe_nan.py`, which writes
the `.nevents` sidecar) and showered with `powheg/powheg_mu1TeV.cmnd` (RES) or
`powheg/powheg_nu_v2.cmnd` (V2, `LesHouches:matchInOut = off`).  A job counts
only with LHE and HepMC beams as intended, >= 90% of the offered events
delivered, and its `V2_OK` marker.

| row | integrator | samples |
|---|---|---|
| POWHEG-RES mu | `$POWHEG_RES/v2-mu<TAG>-<t>/pwg-0001-st3-stat.dat`, "grand total total (pos.-\|neg.\|)" (copied to `INTEGRATED`) | `powheg/v2mu_<t>_job[_TAG]_<seed>/` |
| POWHEG-V2 nu | `$POWHEG_V2/v2-nu<TAG>-<t>/pwg-stat.dat`, "total (btilde+remnants)" | `powheg/v2nu_<t>_job[_TAG]_<j>/` |
| POWHEG-V2mc nu | `$POWHEG_V2/v2-numc<TAG>-<t>/pwg-stat.dat` | `powheg/v2numc_<t>_job[_TAG]_<j>/` |

(The run directories always carry the tag, `1TeV` included.)  POWHEG-V2mc is
the POWHEG-V2 card with `qmass 1.51d0`, `numflav 4` and no `iupperfsr`; it
generates charm production only, so its cross-section is compared with charm
rows and with the FONLL charm reference, never with an inclusive number.

**Herwig 7.3.0 POWHEG NLO.**  `herwig7/production/make_cards.py` writes the
cards in `herwig7/production/cards/` from `herwig7/DIS-{mu,nu}-POWHEG.in` and
`-PWGNEG.in` (only beams, target, PDF, cuts change); the driver is
`herwig7/production/production.sh` (functions in `lib.sh`), with the cut plugin
`herwig7/plugin/LeptonicDISCut`.  Samples
`herwig7/v2{mu,nu}pwg[neg]_<t>_job[_TAG]_N/` hold `events.hepmc`, the `.out`
cross-section table and `fingerprint.log` (md5 of card and plugin; the driver
refuses to mix builds).  The NLO result is positive minus negative half
(`herwig7/combine_nlo.py`); the cross-section bookkeeping, including the
events Herwig discards, is `herwig7/herwig_xsec.py` alone.  ThePEG is run
with `patches/thepeg-2.3.0-event-error-retry.diff`, which retries an event
whose remnant cannot be put on shell instead of discarding it.

**Sherpa 3.0.5 MC@NLO.**  `tools/make_sherpa_cards.py` writes
`sherpa/Runs/V2_{MuonDIS,NuDIS}_NLO_<t>[_TAG]/Sherpa.yaml` from `MuonDIS_NLO`
and `NuDIS_NLO_ckm3` (full CKM, with the virtual-correction patch); only the
beams, energies, PDF and selectors change.  The driver is
`tools/sherpa_production.sh run|check|status|ahadic`.  Run directories
`$SHERPA_RUNS/V2_{MuonDIS,NuDIS}_NLO_<t>[_TAG]/` hold `integ.log`,
`Results.zip` and `job_N/{evtfull,sherpa.log}`; the process libraries are
built once per current and linked.  The delivered cross-section is
sum(Weight)/sum(NTrials), not the mean weight.  The driver refuses a run whose
`integ.log` lacks "Massive PS flavours for Comix: (none)".

**GENIE 3.06.02.**  `genie/production/genie_job.sh spline|run <cfg> <t> <E> <nev> <njobs>`,
one lane per configuration through `genie/production/production.sh <cfg>...`:

| cfg | result key | tune, process list | `GXMLPATH` overlay | energies |
|---|---|---|---|---|
| mu_grv | `genie` (mu) | G18_02a, EMDIS | hienergy:p8 | 400, 700, 1000 |
| mu_nnpdf | `genie_nnpdf` (mu) | G18_02a, EMDIS | nnpdf | 400, 700, 1000 |
| nu_grv | `genie_lo` (nu, FASER's default tune) | G18_02a, CCDISCHARM | nu_cc_lo:p8 | 400, 700, 1000 |
| nu_nnpdf | `genie_nnpdf` (nu) | G18_02a, CCDISCHARM | nu_cc_lo:nnpdf | 400, 700, 1000 |
| nu_hedis | `genie` (nu) | GHE19_00a, CCHEDIS | p8 | all five |

Splines are made per nucleon (`genie/splines/v2/<cfg>_<t>_e{2000,8000}.xml`).
Samples `genie/v2g_<cfg>_<t>_job[_TAG]_<N>/` are complete with `V2_OK`; the
`.ghep.root` is removed and `genie.log` digested into `rejections.json` and
`genie.log.tail`.  The G18_02a neutrino charm channel is run with
`patches/genie-aivazis-charm-propagator.diff` (the W propagator multiplied
where it should divide); see `genie/README.md`.

## Analysis and result files

`tools/analyse_production.py [--energies ...] [--region q4w3|q4w5|q4w3_faser_e|...] [--tag charmfinal] [--only ...]`
runs the whole chain in dependency order: per nucleon (`analysis/analyze.py`
for the muon, `analysis/analyze_nu.py` for the neutrino, `herwig7/make_histos_nlo.py`
for Herwig), Herwig positive minus negative, then the tungsten combination.
It re-runs a step only when an input is newer than its result.

Results are `results/` (muon NC) and `results_nu/` (neutrino CC):

    histos_<gen>[_charmfinal]_<region>_{p,n,W}[_<TAG>].json

with `<region>` one of `q4w3`, `q4w5` (charm), `q4w3_faser_e` (FASERnu emulsion
selection, 1 TeV only), and `_charmfinal` the hadron-level charm tag.  `<gen>`:

| | muon NC (`results/`) | neutrino CC (`results_nu/`) |
|---|---|---|
| POWHEG-RES / POWHEG-V2 + Pythia 8 | `powheg` (RES) | `powheg_nu` (V2), `powheg_nu_mc` (V2mc) |
| Herwig 7 | `herwig_nlo_powheg_full` (halves `herwig_nlo_powheg`, `_neg`) | `herwig_nlo_full` (halves `herwig_nlo`, `_neg`) |
| Sherpa | `sherpa` | `sherpa_nlo` |
| GENIE | `genie`, `genie_nnpdf` | `genie_lo`, `genie_nnpdf`, `genie` (HEDIS) |

Every result stores its input manifest (files, sizes, mtimes) and the
delivery closure of the weighted samples (delivered event-level cross-section
against the integrator; beyond 2% no result is written).

## Analytic references

`tools/make_references.sh [grids|refs|combine|bands|all|list]`: PineAPPL grids,
YADISM in ZM-VFNS (inclusive), YADISM in FONLL (charm, and the FONLL inclusive
= ZM-VFNS inclusive - ZM-VFNS charm + FONLL charm), and the 7-point scale
bands, for the region on target W through `NNPDF40_nnlo_as_01180_W184free`
(structure functions are linear in the PDF, so this is the per-nucleon
tungsten result exactly; checked against the p/n combination to 1e-9).  Each
reference is made with and without target-mass corrections (`BENCH_TMC=3`,
the exact Georgi-Politzer form, set in `analysis/target.py`); the figures
show TMC on, an effect of at most -0.14% on sigma_fid.  muF and muR are
floored at the PDF's Q0^2 = 2.7225 GeV2.  There is no FONLL at NNLO for CC.

Names: `histos_yadism[_charm]_<order>[_fonll_damp]_q4w3_W[_tmc][_TAG].json`
(no suffix = ZM-VFNS, `_fonll_damp` = FONLL) and
`mhou_sigma_fonll{,_mu}_q4w3_W[_tmc][_pto2][_zm].json` for the bands.

## FASER ladders (paper Sect. 7)

The FASER predictions fold a fixed-energy sigma(E) ladder with the flux.
They use FASER's own cuts rather than the benchmark region, with no y window
at generation, for nu and nubar on p and n.  The emulsion, Herwig and Sherpa
ladders use a genuine neutron beam; the POWHEG-V2 rate and dimuon ladders run
a proton beam with the mirrored set, which leaves a proton remnant but the
right cross-section and muons:

| driver | what | energies [GeV] |
|---|---|---|
| `tools/powheg_v2_faser_ladder.sh` | POWHEG-V2 LHE only (Q2 and y are leptonic), event rates | 30-6800 (10 points) |
| `tools/faser_emulsion_ladder.sh`, `tools/faser_emulsion_ladder_ubexcess.sh` | POWHEG-V2 showered and Sherpa neutron points for the emulsion comparison; reduced to `faserdata_events.npz` | 30-6800 |
| `tools/sherpa_faser_ladder.sh` | Sherpa MC@NLO at FASER's cuts | 200-4000 by default; 100-6800 on disk |
| `tools/herwig_faser_ladder.sh` | Herwig 7 POWHEG NLO, cards `herwig7/production/cards/FL-*.in` | 30-6800 (11 points) |
| `tools/genie_dimuon_ladder.sh ladder\|emulsion\|flux` | GENIE default tune, full CC list: dimuon cut flow, emulsion ladder, and a flux-driven run as a check on the interpolation | 10-6800 |
| `tools/powheg_v2_dimuon_ladder.sh`, `..._pdf.sh` | POWHEG-V2 showered, reduced to `dimuon_events.npz` (every muon) | 30-6800 |
| `tools/genie_nubar_sidis.sh`, `tools/sherpa_nubar_sidis.sh` | antineutrino samples for the hadron yields | 300-4000 |

The showered ladders keep only the per-event tables; the HepMC is deleted
once extracted, and the LHE and grids are kept, so a re-shower is cheap.

## Traps a re-runner must know

- **The isospin swap.**  Give the `_n` set only to Sherpa.  POWHEG-RES, POWHEG-V2,
  Pythia 8 and ThePEG swap the proton set for a 2112 beam themselves.
- **Generation cuts at the region edge are not safe.**  The lepton moves after
  the hard process in every showered generator, so a coincident cut leaves a
  hole of about 1% in the first Q2 or W bin.
- **POWHEG-V2 needs `ubexcess_correct 1`** (set by `make_cards.py`).  Without it
  the upper-bound violations of the unweighting leave the LHE 6-8% low at
  x = 0.005-0.011, while the total still closes.  `powheg/production/regen_ubexcess.sh`
  regenerates existing samples in place with the same grids and seeds.
- **Pythia drops POWHEG-V2 charm events unless `LesHouches:matchInOut = off`.**
  About 2% of events, nearly all with an outgoing massless charm, fail with
  "setting mass failed", which biases charm-tagged observables by -8 to -10%.
  `analyze_nu.py` refuses a POWHEG-V2 job whose `shower.log` lacks the switch.
- **Heavy-quark thresholds in NC.**  In NC, an incoming c (b) leaves its
  antiquark in the remnant, so a hadronic final state needs W above about
  4.4 GeV (11 GeV).  The ZM calculation populates charm down to W = 3, and
  Pythia ("parton+hadronLevel failed") and Herwig ("Can't put the remnant
  on-shell") discard those events.  No shower setting cures it: it is energy
  conservation.  Charm-tagged POWHEG-RES is 3% (1 TeV) to 7% (400 GeV) low
  at W > 3, below 0.4% at W > 5, hence the `q4w5` charm region.  The charm tag
  excludes charm from b decays, so the untouched b threshold does not leak in.
- **Pythia re-masses c and b**, moving Q2 down in 6.9% of POWHEG-V2 events
  and raising the lowest x bin by up to 9% at hadron level against the LHE.
- **Sherpa hadronisation failures are vetoed.**  AHADIC failures become new
  events, with their trials carried over, so they lower sum(W)/sum(NTrials)
  and deplete low W.  Without a y cut the muon rate is 0.4%, the neutrino
  rate below 0.1% (`tools/sherpa_production.sh ahadic`).  `RESPECT_MASSIVE_FLAG: true`
  must stay in every card, or 9-23% of events are silently discarded.
- **Herwig discards events whose remnant cannot be put on shell**, mostly sea,
  charm and high-x events, which left its charm 5-26% low.  The ThePEG retry
  patch (`MaxEventErrorRetries`) rebuilds such an event on the same phase-space
  point; without it `herwig_xsec.py` restores the rate but not the shape.
  ThePEG must be rebuilt against the same macOS SDK as Herwig, or `Herwig read`
  fails (see the patch header).
- **GENIE overlays.**  `GXMLPATH` fails open: a wrong path silently falls back
  to the built-in tune.  The `p8` overlay has no `CommonParam` and inherits
  `GVLD-Emax = 1000`, which clamps the e2000 spline, so `config/hienergy` goes
  in front for G18_02a.  It must never go in front of HEDIS, whose
  `CommonParam` it would shadow.  G18_02a is not run above 1 TeV.
- **Glob the job directories exactly.**  The 1 TeV base `v2mu_p_job` must match
  `<base>_<integer>`; a `<base>_*` glob also catches `v2mu_p_job_400GeV_1`.
- **Neutron targets** carry m_n in the HepMC record; every observable uses M_P.
- **FASERnu tracks** are charged particles with E_lab > 1 GeV
  (`selection.TRACK_E_MIN`) and a signed lab p_z (`analyze.track_tan`).
