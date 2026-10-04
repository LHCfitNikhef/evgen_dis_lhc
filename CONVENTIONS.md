# Working conventions

The rules this code was written under.  Comments throughout the repository
cite them by number ("CONVENTIONS.md rule 2b").  Most exist because something
went wrong silently once: a run completed, the cross-section looked plausible,
and the result was wrong.

## 1. Everything lives in the repository; nothing names a path

- Everything written for the project (drivers, analysis code, generator cards,
  configuration overlays, environment setup) is in this repository.  External
  codes (generators, libraries) live under `$PHYSICS24/software/`; local
  changes to them are kept here as diffs in `patches/`.
- A generator's install tree gets a **symlink** back into the repository, never
  a copy (`tools/link_cards.sh`, `--check` to report).
- **No tracked file contains an absolute path.**  Every external location is
  named once in `config.sh`; shell scripts source it, Makefiles include
  `config.mk`, Python does `import paths`.  A machine's own settings go in the
  untracked `config.local.sh`.  `tools/check_no_abs_paths.sh` enforces it.
- Every knob is spelled `BENCH_*` and the working variables are assigned, not
  defaulted on generic names: conda exports `PYTHIA8_DIR`, `PYTHIA8DATA` and
  `CXX`, and a default on those names silently picks up a different Pythia
  version (8.312 instead of 8.311) or an unusable compiler.
- Scripts run on macOS and Linux and target bash 3.2 (no associative arrays,
  no zsh-isms; `"${ARR[@]}"`, never `$ARR`, which takes the first element
  only).  Platform differences are decided once in `config.sh`
  (`BENCH_CXX`, `BENCH_NOSLEEP`, `BENCH_LIBPATH_VAR`).
- **Intermediates do not outlive their results.**  An event file is deleted
  once every result read from it is written and checked
  (`tools/analyse_and_prune.py`), except the production samples a new
  selection would have to re-parse.  A step that re-derives a tracked result
  from event samples refuses to overwrite it when the samples are absent.
- `.gitignore` rule of thumb: if a generator wrote it, it is ignored; if we
  wrote it, it is tracked.

## 1b. Disk floor; GENIE is always shown

- Productions check the free space before starting (a floor of 150 GB on the
  volume holding `$PHYSICS24`): a generator that runs out of disk mid-write
  leaves a truncated file that still parses.
- GENIE, the generator FASER uses, appears in every comparison, in both its
  default tune (G18_02a, "GRV98LO") and HEDIS.  Below the statistics floor
  (`analyze.MIN_SELECTED_EVENTS`) the answer is more GENIE events, not a lower
  floor; where GENIE genuinely has no prediction (no charm in its NC DIS), that
  is stated in place.  The paper figures (`analysis/paper_plots/`) may omit it
  when a figure carries a different single message, and then say so.

## 1b'. One fiducial region; the paper figures are kept current

Only the benchmark region (Q² > 4 GeV², W > 3 GeV, no y cut, tungsten)
exists.  When the setup changes, the inputs of the paper figures and their
claims are recomputed; any other figure that would go stale is removed or says
in place that it predates the change.

## 1c. Never edit a shell script while it runs

Bash reads a script lazily from a byte offset, so editing a running driver
makes it resume in the middle of a different line.  Edit a copy, or stop the
driver, edit, and re-invoke it: every driver skips work whose output exists.

## 2. A plausible cross-section is not a validation

- **Delivery closure** (`analyze.py: closure_check`): every weighted sample
  compares its delivered event-level cross-section with the integrator's, and
  no result is written beyond 2%.  For Sherpa it is
  `sum(weight) / sum(trials)`, not the mean weight.
- **Input manifest** (`analyze.py: input_manifest`): every result records the
  files, sizes and times it came from, and warns when they span more than a
  day (the sign of a glob mixing old and new job directories).
- **One source per normalisation rule** (e.g. Herwig's cross-section
  bookkeeping only in `herwig7/herwig_xsec.py`).
- Jobs are combined with inverse-variance weights for repeated integrals and
  accepted-count weights for pooled event samples, never a plain mean.

## 2b. The two currents stay in sync

Whatever is done for muon NC DIS is done for neutrino CC DIS and vice versa
(cards, run scripts, analysis switches, energies, references, figures).  A
genuine asymmetry is stated where it occurs (e.g. no charm in GENIE's NC
DIS).  The SIDIS study of Sect. 7 shows the neutrino current only, by choice.

## 3. Naming and physics conventions

- **YADISM always names its scheme**: ZM-VFNS (inclusive references) or FONLL
  (charm).  FONLL always includes threshold damping.
  (`tools/check_yadism_scheme.py`)
- **POWHEG is always POWHEG-RES** (`DIS_v`, muon NC) **or POWHEG-V2**
  (`nu-DIS`, neutrino CC); "POWHEG" alone means the matching method.
  (`tools/check_powheg_naming.py`)
- GENIE legend labels: GENIE (GRV98LO), GENIE (NNPDF4.0), GENIE (HEDIS).
- Benchmark settings: NNPDF4.0 NNLO; μ_F = μ_R = Q; Gμ scheme with fixed
  α = 1/132.119 and sin²θ_W = 0.22291; no QED radiation; massless charm and
  bottom; stable if cτ > 10 mm; every cut in the target rest frame.

## 4. The paper's numbers are checked

Each paper figure is one script in `analysis/paper_plots/` whose quoted
numbers are re-derived from the result files (`tools/check_paper_plots.py`),
and every number in the paper text has an entry in
`tools/check_paper_claims.py`.  The figure PNGs used by the paper are tracked
in `paper/figures/`, so the paper builds from a clone.

## 5. The HTML report

`analysis/make_report.py` builds a self-contained page (`results/report.html`,
not tracked) with the figures embedded.  Its prose is sentence case (emphasis
in bold, not capitals; `tools/check_report_caps.py`) and every number it
quotes is checked by `tools/check_report_claims.py`.
