# POWHEG-BOX a-posteriori reweighting: scale variations and PDF sets

Event kinematics are not regenerated.  POWHEG-BOX stores enough per event
(`storeinfo_rwgt 1`, set in every production card) to re-evaluate the matrix
element with another scale or PDF; with `rwl_add 1` it reads the Les Houches
file and appends the extra weights.  The weights are then harvested into a
compact array (`harvest_weights.py`) and the reweighted LHE is deleted.

## Reweighting lists

| file | content | used by |
|---|---|---|
| `rwl_scale.xml` | the 7-point (ξ_R, ξ_F) scale variation | `powheg/production/reweight_scale.sh` (MHOU on hadron-level distributions), `tools/powheg_v2_faser_ladder_rwgt.sh` |
| `rwl_members.xml/.json` | every member of NNPDF4.0, CT18, MSHT20, ATLASpdf21, ABMP16, plus the scale points (`make_rwl_members.py`) | `run_reweight.sh`, `run_reweight_members.sh` |
| `rwl_nuclear.xml/.json` | nuclear-PDF sets for tungsten, with their free-nucleon baselines | `run_reweight_members.sh --rwl rwl_nuclear` |
| `rwl_dimuon_{p,n}.xml/.json` | the PDF sets above plus GRV98, on proton and on isospin-mirrored neutron sets, same weight ids (`make_rwl_dimuon.py`) | `tools/powheg_v2_dimuon_ladder_pdf.sh` |
| `rwl_npdf_hadron.xml/.json` | every member of nNNPDF3.0, EPPS21 and nCTEQ15HQ for tungsten, their free-nucleon baselines, NNPDF4.0 free tungsten nucleon and the closure (`make_rwl_npdf_hadron.py`) | `run_reweight_npdf.sh` (hadron-level nPDF bands, paper plot A2b, neutrino only) |
| `rwl_scale_pdf.xml` | scale points + central members only | `BENCH_RWL=rwl_scale_pdf.xml run_reweight.sh` |

## Scripts

- `run_reweight.sh [--current mu|nu] <E>` — one POWHEG-V2 sample (POWHEG-RES
  writes one LHE per seed and is reweighted per seed by
  `powheg/production/reweight_scale.sh`).
- `run_reweight_members.sh` — all members, in batches (`split_rwl.py`):
  `pwhg_main` aborts after roughly twenty PDF members in one pass.
- `run_reweight_selected.sh`, `subset_lhe.py` — reweight only the events that
  pass a selection.  This is exact, not an approximation: events outside the
  selection contribute nothing to the member dependence of the selected sum,
  and the normalisation is member-independent.
- `run_reweight_npdf.sh` — the nPDF members on the FASERν-selected events of
  the 1 TeV POWHEG-V2 proton sample, in chunks, six members per pass; joined
  to the showered events by `analysis/npdf_hadron.py`.  Only the PROTON
  sample is reweighted: the nPDF grids are the average tungsten nucleon, and
  the neutron sample (`ih2 2`) would swap u and d in them.
- `truncate_lhe.py` — reweight the first N events of a very large sample.

## Traps

- **The keyword is `lhapdf`, and only `lhapdf`.**  POWHEG-BOX's
  `rwl_setup_params_weights` tests exactly `lhapdf`, `facscfact` and
  `renscfact`.  Writing the production card's own `lhans1`/`lhans2` is
  accepted by the XML parser and **silently ignored**: every PDF weight then
  equals the nominal while the scale weights vary correctly.  Check that the
  PDF weights differ from the closure weight.
- **The closure weight** re-evaluates each event with the PDF it was
  generated with and must equal the nominal weight.
- **GRV98 cannot be addressed by number.**  POWHEG-BOX selects PDFs by numeric
  LHAPDF id, and `GRV98lo`'s declared `SetIndex` 80060 collides with another
  set's range in `pdfsets.index`, so LHAPDF tries to load a set that is not
  installed.  It is dropped from `rwl_members` (recorded in its JSON) and
  reached in the dimuon lists through the `pdfsets.index` of
  `data/pdfs/lhapdf/`, prepended to the LHAPDF search path.
- **POWHEG-BOX-RES cannot reweight a file mixing seeds** ("setrandom: failed"):
  each event restores its seed's random state, which only skips forward.
- **Parallel passes must not share the run directory's small files**: through
  symlinks pwhg_main rewrites FlavRegList, `*_equiv`, `pwhg_checklimits` and
  `pwgcounters` in the production directory, and concurrent writers crash.
  Copy them (`run_reweight_npdf.sh` does).
- **Different PDF sets in one pass abort sooner than members of one set**:
  about ten sets; six members per pass is safe.
- **nCTEQ15HQ_FullNuc_184_74 (SetIndex 122950) is not in conda's
  `pdfsets.index`**, and 122950 falls inside nCTEQ15WZSIH's range there; it is
  listed in `data/pdfs/lhapdf/pdfsets.index`.
