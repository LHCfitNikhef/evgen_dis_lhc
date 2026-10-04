# Patches to the external codes

The generators are built from source outside this repository; these are the
local modifications the benchmark depends on.  Each diff states in its header
the exact `patch -pN -d <dir>` line to apply it (the `-p` level differs between
codes) and why it is needed.  Rebuild the code after patching.

| Patch | Code | What it fixes |
|---|---|---|
| `pythia8311-StandardModelParameters.diff` | Pythia 8.311 | widens the allowed `alphaEMmZ`/`sin2thetaW` ranges so the benchmark's Gμ-scheme values are accepted |
| `sherpa-dy-qcd-virtual-ckm.diff` | Sherpa 3.0.5 | serves the Cabibbo cross-family CC channels in the internal virtual ME, needed for CC MC@NLO with the full CKM matrix |
| `sherpa-hadron-remnant-neutron.diff` | Sherpa 3.0.5 | a neutron (2112) beam built a wrong remnant diquark and crashed; applied to a separate install (`tools/build_sherpa_npatch.sh`, macOS only) |
| `thepeg-2.3.0-event-error-retry.diff` | ThePEG 2.3.0 (Herwig 7.3.0) | retries an event that fails after its hard process was accepted on the same phase-space point (or LHE event) instead of discarding it, which had biased the charm and low-Q² rates |
| `powheg-v2-nu-dis-ckm.diff` | POWHEG-V2 `nu-DIS` | the CKM routine had no branch for b quarks: b-initiated Born points got an uninitialised weight and an outgoing b a Cabibbo weight |
| `powheg-res-dis-v-muflux-lepton-flux.diff` | POWHEG-RES `DIS_v` | lepton-flux hook, in a separate copy of the process (`DIS_v_muflux`), for the muon-flux comparisons |
| `genie-em-q2-floor-and-pythia8-teardown.diff` | GENIE 3.06.02 | raises the EM Q² generation floor to 4 GeV² (events generated inside the fiducial region) and makes Pythia 8 keep particles with cτ > 10 mm stable, the benchmark convention |
| `genie-aivazis-charm-propagator.diff` | GENIE 3.06.02 | the CC charm cross-section multiplied by the W-propagator factor instead of dividing (too large by (1+Q²/M_W²)⁴ at high Q²) |
| `mg5-3.7.2-dis-nlo-guard.diff` | MG5_aMC 3.7.2 | downgrades the refusal of lepton–hadron NLO to a warning, for the (unvalidated) NLO attempt only; the benchmark uses MG5_aMC at LO |
| `gibuu-dis-final-state-buffer.diff` | GiBUU 2017 | the DIS final-state buffer (20 particles) overflowed at TeV energies |
| `gibuu-masternbody-local-k.diff` | GiBUU 2017 | a loop variable made local in `ResetPosition` |
