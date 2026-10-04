# GiBUU in the benchmark

GiBUU (Giessen Boltzmann–Uehling–Uhlenbeck) is a transport model: it
propagates the hadrons out through the nucleus.  It enters only the appendix
comparison of GENIE with NuWro and GiBUU at E_ν = 200 GeV, **on the free
proton only** (see below).

- Upstream: <https://gibuu.hepforge.org>; installed at `$GIBUU_DIR` with its
  data package at `$BUUINPUT` (`config.sh`: `BENCH_GIBUU_DIR`,
  `BENCH_BUUINPUT`).

## Before using it

1. **This is the 2017 release.**  The hepforge download page answers
   automated requests with a proof-of-work HTML page, and the GitHub mirror
   (`github.com/gibuu`) stops at `GiBUU_2017` + `buuinput_2017`.  A newer
   release has to be downloaded by hand.
2. **Its neutrino module is documented for 1–50 GeV** (arXiv:1205.1061).
   Anything above is an extrapolation.
3. **Tungsten transport does not complete at 100–200 GeV**: the perturbative
   particle buffer overflows and the collision criterion reports
   probabilities above one (time step too coarse for these multiplicities);
   at 1 TeV it aborts before writing any event.  On tungsten it completes at
   50 GeV.  Hence no GiBUU tungsten row in the appendix.

## Building it

    PATH="$BENCH_REPO/tools/compat:$PATH" make -j 4      # in $GIBUU_DIR

- `tools/compat/gfind` supplies the `.` path GiBUU's Darwin Makefile omits
  in its GNU-only `gfind -maxdepth 1` call; on PATH for the build only.
- `patches/gibuu-masternbody-local-k.diff`: `ResetPosition` reused the host's
  DO variable `k`, collapsing a ten-attempt retry loop to one; gfortran 14
  rejects it.
- `patches/gibuu-dis-final-state-buffer.diff`: the DIS final-state buffer
  (20 particles) overflows at high energy after a few thousand events;
  raised to 200.

## Running it

    gibuu/run_gibuu.sh --current nu --target p 200 [ensembles] [runs] [job]

A jobcard with `eventtype = 5`, monoenergetic beam, every channel but DIS
switched off.  `p`, `n`: free nucleon, no transport (`numTimeSteps = 0`);
`W`: tungsten with transport.  Output is Les Houches (GiBUU 2017 has no
HepMC); `gibuu2hepmc.py` converts it into `gibuu/gibuu_job_<cur>_<t>[_<E>]_N/`.
Two features of that file a converter must handle:

- the outgoing lepton is not in the particle block — the neutrino and lepton
  four-momenta are on the `#` comment line closing each event, as
  `(E, px, py, pz)`;
- `XWGTUP` is in units of 10⁻³⁸ cm², so Σ XWGTUP × 0.01 is σ in pb.
