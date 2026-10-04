# QED arms for the POWHEG-RES and POWHEG-V2 shower

Pythia overlay cards, **appended** to the base Pythia steering of a
POWHEG-RES or POWHEG-V2 sample, that switch on QED radiation in stages.  Pythia's `readFile` applies
settings in order and the last one wins, so the base card stays the single
source of every non-QED setting.

| arm      | switches turned on                              | currents |
|----------|-------------------------------------------------|----------|
| `off`    | none — closure check                            | mu, nu   |
| `fsr`    | `TimeShower:QEDshowerByL`                       | mu, nu   |
| `fsrisr` | + `SpaceShower:QEDshowerByL`                    | mu only  |
| `full`   | + `QEDshowerByQ` in both showers                | mu, nu   |

Driver: `powheg/production/qed_arms.sh <mu|nu> <p|n> <arm>` (and `closure`).
It re-showers the **same 1 TeV Les Houches files** as a subset of the
production jobs — POWHEG-RES for the muon, POWHEG-V2 for the neutrino, proton
and neutron separately — with the same base card (`powheg_mu1TeV.cmnd`,
`powheg_nu_v2.cmnd`), into `powheg/qed<arm>_<baseline job>/`.  Each arm is
compared with its own baseline jobs, so the matrix element, matching, PDF and
scale cancel and only the QED radiation differs.  The older
`powheg/run_powheg_shower.sh --qed <arm>` reads the same overlays.

## Why `fsrisr` is muon-only

Initial-state radiation off the lepton line needs a charged incoming lepton.
The charged current starts from a neutrino, so the arm does not exist there;
the driver refuses it rather than silently producing a copy of `fsr`.

## The closure

`main_powheg.cc` never sets `Random:setSeed`, so Pythia runs from its fixed
default seed.  The `off` arm must therefore reproduce the baseline sample
byte for byte; `qed_arms.sh closure` checks this.  It also means every arm
shares the random sequence until QED first changes the event, so the arms are
correlated with the baseline and the error on a difference is far below the
error on either sample.

## Caveat

The benchmark fixes α_em = 1/132.119 (Gμ scheme) partly because there is no
QED radiation.  Turning radiation on without revisiting the EW scheme
double-counts part of the correction, so the arms measure **the shower-level
QED migration across the fiducial cuts**, not the full QED correction to DIS.
`TimeShower:QEDshowerByGamma`/`ByOther` stay at Pythia's defaults (on); they
only act once a photon has been radiated.
