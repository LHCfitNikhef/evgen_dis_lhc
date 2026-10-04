#!/bin/bash
# GENIE environment setup (macOS arm64, conda ROOT/LHAPDF/PYTHIA8 + homebrew log4cpp)
#
# This file lives IN THE REPOSITORY (evgen-benchmark/genie/genie_setup.sh).
# ~/genie_setup.sh is a symlink to it, so `source ~/genie_setup.sh` keeps
# working; edit the repo copy, never the symlink target's old location.
#
# Every path outside the repo comes from ../config.sh -- the single source --
# so this file declares none of its own.  Relocate the tree with PHYSICS24, or
# an individual piece with its BENCH_* knob:
#     PHYSICS24=/data/physics24 source genie_setup.sh
#     BENCH_DEPS=/cvmfs/... BENCH_GENIE_DIR=/opt/genie source genie_setup.sh

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$_HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"

export GENIE=$GENIE_DIR
export ROOTSYS=$BENCH_DEPS

# NOTE: GENIE is built against $BENCH_DEPS's Pythia8, which is NOT the 8.311
# the rest of the benchmark uses -- so PYTHIA8DATA is deliberately overridden
# back to the dependency prefix here, unlike everywhere else.
export PYTHIA8_INC=$BENCH_DEPS/include
export PYTHIA8_LIB=$BENCH_DEPS/lib
export PYTHIA8DATA=$BENCH_DEPS/share/Pythia8/xmldoc

# LHAPDF6
export LHAPDF6_INC=$BENCH_DEPS/include
export LHAPDF6_LIB=$BENCH_DEPS/lib
export LHAPATH=$LHAPDF_DATA_PATH

# log4cpp
export LOG4CPP_INC=$BENCH_LOG4CPP/include
export LOG4CPP_LIB=$BENCH_LOG4CPP/lib

# NOTE: /usr/bin must come before $ROOTSYS/bin (conda's bin on this machine) so that
# clang++ is Apple clang. Conda's clang-18 ships an availability model that is
# incompatible with the macOS 26 SDK and fails on <charconv>/<chrono>.
# conda's libc++ headers mark <charconv> helpers as "introduced = 99.0" against
# the macOS 26 SDK, which cling/rootcling cannot parse. Disabling the libc++
# availability annotations fixes dictionary generation (and runtime cling).
export EXTRA_CLING_ARGS=-D_LIBCPP_DISABLE_AVAILABILITY

# /usr/bin first ONLY on macOS (Apple clang, above); on Linux it would put the
# system g++/python ahead of an LCG view's and mismatch the view's ROOT.
# The library path is extended on BOTH variables from their OWN previous
# value: copying DYLD_ into LD_ wiped an LCG setup's LD_LIBRARY_PATH on Linux.
_genie_libpath () {
    export DYLD_LIBRARY_PATH="$1${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
    export LD_LIBRARY_PATH="$1${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}
case "$(uname -s)" in
    Darwin) export PATH=/usr/bin:$ROOTSYS/bin:$GENIE/bin:$PATH ;;
    *)      export PATH=$ROOTSYS/bin:$GENIE/bin:$PATH ;;
esac
_genie_libpath "$GENIE/lib:$ROOTSYS/lib:$BENCH_LOG4CPP/lib"

# PYTHIA6 (built with $GENIE/src/scripts/build/ext/build_pythia6.sh)
export PYTHIA6=$PYTHIA6_LIB
_genie_libpath "$PYTHIA6"

# APFEL (built from github.com/scarrazza/apfel, needed for HEDIS NLO structure
# functions).  gevgen links libAPFEL.0.0.0.dylib by BARE NAME, so it is found
# only through DYLD_LIBRARY_PATH -- if APFEL_LIB is wrong, gevgen dies at
# load time rather than falling back.
export APFEL_INC=$APFEL_DIR/include
export APFEL_LIB=$APFEL_DIR/lib
_genie_libpath "$APFEL_LIB"

# ROOT 6.38 no longer ships the EGPythia6 module, so GENIE is built PYTHIA8-only.
# This overlay swaps the Pythia6 hadronizer/decayer algorithms for their Pythia8
# equivalents in the default tune, and enables charm decays (benchmark
# convention: stable if ctau > 10 mm).  It is BENCHMARK PHYSICS CONTENT, not
# just plumbing -- see genie/README.md.
export GXMLPATH=$BENCH/genie/config/p8

# GXMLPATH pointing at a missing directory does NOT error: GENIE silently falls
# back to the built-in tune (Pythia6 hadronizer, charm STABLE), which changes
# nch/Elead/D-meson observables without any warning.  Fail loudly instead.
if [ ! -d "$GXMLPATH" ]; then
    echo "genie_setup.sh: ERROR: GXMLPATH=$GXMLPATH does not exist." >&2
    echo "  GENIE would silently run the built-in tune instead of the" >&2
    echo "  benchmark overlay.  Refusing to set it." >&2
    unset GXMLPATH
    return 1 2>/dev/null || exit 1
fi

# ---------------------------------------------------------------------------
# THE SPLINE ENERGY RULE, in ONE place (all three run_genie*.sh read it).
#
# A GENIE spline evaluated AT its last knot reads past the knot array
# (Spline::FindClosestKnot -> GetKnot(iknot+1)) and returns 0.  That does NOT
# fail: gevgen accepts a zero cross-section and spins forever in "Could not
# select interaction".  So a spline is always built well past its beam energy
# -- the 1 TeV setup uses 2x, and that headroom is kept at every energy.
#
# The committed _e2000 splines cover 400 GeV and the 1 TeV ANCHOR; 4 TeV needs
# _e8000, built with `run_genie*.sh spline <E>`.  Returning e2000 for BOTH 400
# and 1000 is deliberate: the anchor must keep evaluating the exact spline its
# published sigma_gen came from, and a rebuild with different knot spacing
# would move a published number for no physics reason.
#
# Energies are enumerated rather than computed, so an unlisted one fails LOUDLY
# here instead of silently selecting a too-short spline and hanging gevgen.
genie_spline_emax() {
    case "${1:-1000}" in
        # 300 GeV is an OFF-SCAN comparison point (2026-09-07), the beam
        # energy of the FPF prediction in arXiv:2504.05376; it rides the same
        # e2000 spline as 400 and 1000, so it costs no spline build.  It is
        # deliberately NOT in beams.ENERGIES -- see tools/make_energy_cards.py.
        300|400|1000) echo 2000;;
        4000)     echo 8000;;
        *) echo "genie_spline_emax: no spline range defined for E=$1 GeV." >&2
           echo "  Add one here rather than guessing: a spline that is too" >&2
           echo "  short does not fail, it hangs gevgen forever." >&2
           return 1;;
    esac
}
