# Pythia 8 shower-model arms

Overlay cards of pure Pythia settings, appended to a base card by
`pythia8/run_pythia.sh --shower <arm>` (LO samples, `shw<arm>_job_N/`) or
`powheg/run_powheg_shower.sh --shower <arm>`.  The base card stays the single
source of every setting an arm does not change.

| arm        | setting                          | what it is |
|------------|----------------------------------|------------|
| `simple`   | none (model 1 is the default)    | the default p_T-ordered shower; must reproduce the baseline |
| `vincia`   | `PartonShowers:model = 2`        | antenna shower |
| `dire`     | `PartonShowers:model = 3`        | dipole shower |
| `norecoil` | `SpaceShower:dipoleRecoil = off` | diagnostic: separates the recoil choice from the shower algorithm |

These are the "also Vincia, Dire" alternatives listed for Pythia in the
paper's generator table.  The NLO shower comparison in the paper is made
differently: POWHEG-RES / POWHEG-V2 Les Houches files showered by Pythia 8 and by Herwig 7,
and Sherpa with CSS and with Dire.

## Two confounds

**1. The POWHEG-matched path is blocked.**  `main_powheg.cc` installs the
POWHEG-matching veto through `UserHooks` that intercept SimpleShower emission scales; Vincia
and Dire do not expose that interface.  On one seed of 59775 LHE events,
`simple` reproduced the production sample byte for byte, `vincia` delivered 0
events ("couldn't find Pythia FSR emission") and `dire` 225 ("couldn't find
Pythia ISR emission") — **and both runs exited 0**.  The delivery guard in
`run_powheg_shower.sh` exists to catch this.  Switching the veto off instead
would double-count the hardest emission, so the Pythia shower-model
comparison is made on the LO samples only.

**2. Vincia and Dire ship their own tunes.**  Changing the model does not
isolate the shower: at the shipped tunes it conflates shower and tune, at a
common tune it runs one shower with another's parameters.  The comparison is
therefore "shower model as shipped", not a shower-algorithm systematic.  The
base card's `SpaceShower:dipoleRecoil = on` is a SimpleShower setting the
other models never read, which is what `norecoil` tests.
