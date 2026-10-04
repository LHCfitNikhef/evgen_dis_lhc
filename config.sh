# ---------------------------------------------------------------------------
# evgen-benchmark -- the ONE place that knows where things are on this machine.
#
# Nothing in this repository may hardcode an absolute path.  Every location
# outside the repo is named here and is overridable from the environment, so
# the same checked-out tree runs on a laptop and on a cluster:
#
#     PHYSICS24=/project/faser ./pythia8/run_pythia.sh 50000 8
#     BENCH_DEPS=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9 make -C pythia8
#
# POSIX sh on purpose: sourced from sh, bash and zsh scripts, and read by
# config.mk (make) and analysis/paths.py (Python) through
# tools/print_config.sh -- so these defaults have exactly ONE definition.
#
# ---------------------------------------------------------------------------
# WHY EVERY KNOB IS SPELLED BENCH_*, AND WHY THE WORKING NAMES ARE ASSIGNED
# UNCONDITIONALLY.
#
# The first draft of this file used `: "${PYTHIA8_DIR:=...}"`.  It resolved to
# /opt/miniconda3, because conda already exports PYTHIA8_DIR and PYTHIA8DATA
# pointing at its Pythia 8.312 -- while this benchmark is built against 8.311.
# A "sensible default" on a generic name is therefore not sensible at all: it
# inherits whatever the ambient environment happens to hold and runs with the
# wrong particle data, silently.  That is the failure mode CONVENTIONS.md rule 2
# exists for.
#
# So: the ONLY things read from the environment are PHYSICS24 and the BENCH_*
# knobs, which nothing else defines.  Every working variable below is then set
# unconditionally from those.  To relocate something, set its BENCH_* knob.
# ---------------------------------------------------------------------------

# --- roots -----------------------------------------------------------------
: "${PHYSICS24:=$HOME/physics24}"                     # external tree + benchmark
: "${BENCH_SOFTWARE:=$PHYSICS24/software}"            # upstream codes (rule 1)
: "${BENCH_REPO:=$PHYSICS24/faser/evgen-benchmark}"   # this repository

PHYSICS24=$PHYSICS24
SOFTWARE=$BENCH_SOFTWARE
BENCH=$BENCH_REPO

# --- per-machine settings ----------------------------------------------------
# A clone on another machine sets its knobs in config.local.sh next to this
# file (NOT tracked: .gitignore), or in the environment -- never by editing
# this file.  Anything below may be set there; the roots above are read first,
# so a clone outside $PHYSICS24 sets BENCH_REPO in the environment (every
# driver does so itself, from its own location).
# shellcheck disable=SC1091
[ -f "$BENCH_REPO/config.local.sh" ] && . "$BENCH_REPO/config.local.sh"

# --- dependency prefixes ---------------------------------------------------
# BENCH_DEPS supplies HepMC3, LHAPDF, Pythia8 and ROOT (headers, libs, data).
# The default is the BASE conda installation, found from $CONDA_EXE, which
# conda sets once and does NOT change when an environment is activated --
# deliberately not $CONDA_PREFIX, which does.  Without conda it is /usr; on a
# cluster point it at the LCG view or module prefix.  tools/check_paths.sh
# reports a prefix that does not hold what is needed.
if [ -z "${BENCH_DEPS:-}" ]; then
  if [ -n "${CONDA_EXE:-}" ]; then
    BENCH_DEPS=$(dirname "$(dirname "$CONDA_EXE")")
  else
    BENCH_DEPS=/usr
  fi
fi
# log4cpp for GENIE, and an extra bin directory the run scripts keep on PATH:
# Homebrew's own prefix on macOS (HOMEBREW_PREFIX, else its standard
# /opt/homebrew), the system prefix elsewhere.
case "$(uname -s)" in
  Darwin) _bench_pkg=${HOMEBREW_PREFIX:-/opt/homebrew} ;;
  *)      _bench_pkg=/usr ;;
esac
: "${BENCH_LOG4CPP:=$_bench_pkg}"
: "${BENCH_EXTRA_BIN:=$_bench_pkg/bin}"
unset _bench_pkg
# THE PYTHON THE DRIVERS RUN.  It needs numpy, scipy, matplotlib, yadism, eko,
# pineappl and the lhapdf and ROOT bindings (environment.yml).  Not a bare
# `python3` and not $PYTHON: both resolve to whatever environment happens to be
# active, the trap this file exists to avoid.
if [ -z "${BENCH_PYTHON:-}" ]; then
  if [ -x "$BENCH_DEPS/bin/python3" ]; then BENCH_PYTHON=$BENCH_DEPS/bin/python3
  else BENCH_PYTHON=python3; fi
fi

# --- generator install trees ----------------------------------------------
# Each is independently overridable: a cluster rarely mirrors one tidy layout.
: "${BENCH_SHERPA_RUNS:=$SOFTWARE/sherpa/Runs}"
# The Sherpa install prefix.  Named here because run_sherpa_shower.sh
# needs the binary, and the older run_parallel.sh reaches it as
# ../../../install/bin/Sherpa from inside a job directory -- a relative
# path that only works at one nesting depth and breaks the moment a run
# directory moves.
: "${BENCH_SHERPA_INSTALL:=$SOFTWARE/sherpa/install}"
# A SECOND, PATCHED COPY of that install, used ONLY for neutron-beam (2112)
# runs (paper plots, 2026-09-13, user: "patched copy for neutron runs").
# Stock 3.0.5 segfaults on the first event of a 2112 beam in
# Hadron_Remnant::RemnantFlavour; patches/sherpa-hadron-remnant-neutron.diff
# fixes it.  Built by tools/build_sherpa_npatch.sh from the shared install
# (copied, libRemnants rebuilt, every rpath rewritten to the copy), so the
# shared install that the proton runs use is never modified.
: "${BENCH_SHERPA_NPATCH:=$SOFTWARE/sherpa-npatch/install}"
: "${BENCH_HERWIG_BIN:=$SOFTWARE/herwig7/local/bin}"
: "${BENCH_GENIE_DIR:=$SOFTWARE/genie-mc}"
: "${BENCH_PYTHIA6_LIB:=$SOFTWARE/pythia6/v6_428/lib}"
# NuWro, and the ROOT-Pythia6 plugin it needs (2026-09-08).  Conda's ROOT 6.38
# ships libEGPythia8 but NOT libEGPythia6, and NuWro's DIS fragmentation is
# TPythia6 -- so the plugin is built separately from
# github.com/luketpickering/ROOTEGPythia6 with its OWN builtin Pythia 6.
# >>> ITS OWN, NOT $BENCH_PYTHIA6_LIB. <<<  That library is GENIE's, built for
# GENIE, and sharing one Pythia 6 between two generators is precisely the kind
# of silent version coupling CONVENTIONS.md rule 2 is about.
: "${BENCH_NUWRO_DIR:=$SOFTWARE/nuwro}"
# GiBUU, and its input-data package (2026-09-08).  The RELEASE TARBALLS ON
# HEPFORGE CANNOT BE FETCHED AUTOMATICALLY -- the site is behind a
# proof-of-work challenge that answers every request, curl included, with an
# HTML page -- so what is built here is the GitHub mirror, which stops at the
# 2017 release.  If a newer one is wanted it has to be downloaded by hand.
: "${BENCH_GIBUU_DIR:=$SOFTWARE/gibuu}"
: "${BENCH_BUUINPUT:=$SOFTWARE/buuinput}"
: "${BENCH_ROOTEGPYTHIA6:=$SOFTWARE/ROOTEGPythia6/build/$(uname -s)}"
: "${BENCH_APFEL_DIR:=$SOFTWARE/apfel-install}"
: "${BENCH_MG5_DIR:=$SOFTWARE/mg5_amc/mg5_amc}"

# Pythia8 8.311, built for this benchmark.  NOTE it sits under faser/, not
# software/ -- an inconsistency with CONVENTIONS.md rule 1 preserved here rather
# than "tidied", because moving it is exactly the silent breakage that rule
# was written after.
: "${BENCH_PYTHIA8_DIR:=$PHYSICS24/faser/pythia8/pythia8311}"

# TWO DISTINCT POWHEG CODES -- never write a bare "POWHEG" (CONVENTIONS.md rule 3).
: "${BENCH_POWHEG_RES:=$SOFTWARE/powheg/powheg-dis-main/DIS_v}"
: "${BENCH_POWHEG_RES_MUFLUX:=$SOFTWARE/powheg/powheg-dis-main/DIS_v_muflux}"   # lepton-flux copy of DIS_v (patches/powheg-res-dis-v-muflux-lepton-flux.diff)
: "${BENCH_POWHEG_V2:=$SOFTWARE/powheg-cmass/nu-DIS-master}"

# FASER-FORMAT PRODUCTION (faser_format/, started 2026-09-23): benchmark events
# converted to FASER's GENIE flat ntuple (.gfaser.root) for calypso.
# BENCH_FASER_GENIE is FASER's GENIE fork (gitlab faser/offline/geniegenerator),
# read for its GDML geometry and flux-entry class -- external code, never edited.
# Its PRODUCTION branch is faser-R-3_04_00; files are read with `git show`, so
# the checkout's own branch does not matter.  BENCH_FASER_DATA holds the
# downloaded Kling flux ntuples and the produced samples (large, not in git).
: "${BENCH_FASER_GENIE:=$PHYSICS24/faser/geniegeneratorFASER}"
: "${BENCH_FASER_GENIE_BRANCH:=origin/faser-R-3_04_00}"
: "${BENCH_FASER_DATA:=$PHYSICS24/faser/faser-format-data}"
# FASER's GENIE fork BUILT HERE for the GENIE twin of each faser_format sample
# (faser_format/genie/build_genie_faser.sh): a local clone of $BENCH_FASER_GENIE
# at $BENCH_FASER_GENIE_BRANCH, so the user's checkout is never touched.
: "${BENCH_GENIE_FASER_DIR:=$SOFTWARE/genie-faser}"

# Shareable copy of the report, written outside the repo by make_report.py.
: "${BENCH_REPORT_COPY:=$PHYSICS24/faser/faser-evgen-benchmark-report.html}"

# THE PAPER'S OWN GIT REMOTE, for the Overleaf round trip (user, 2026-09-08:
# "I find the overleaf editor much nicer to work with than a local text
# editor").  Overleaf is linked to THAT repository, so the editor sees a repo
# containing only the document.  Named here because it is the one place things
# outside this repository are named; tools/sync_overleaf.sh is the only reader.
: "${BENCH_PAPER_REMOTE:=git@github.com:LHCfitNikhef/faser-evgen-bench-paper.git}"
# ...and an ORDINARY CLONE of it, which is how the round trip is made.
#
# >>> NOT A `git subtree` (2026-09-08, after it lost a section). <<<  A
# subtree push builds a SYNTHETIC history out of this repository's commits.
# A commit the editor makes on the remote is not in that history and cannot
# become part of it, because `git subtree pull` refuses ("unrelated
# histories") on a branch a `subtree push` created -- so every later push is
# a non-fast-forward, every push is refused, and the only way anything moves
# is a wholesale overwrite in one direction or the other.  A plain clone has
# the remote's real history in it, so a push is an ordinary fast-forward and
# a pull is an ordinary diff.  It lives outside the repository because it is
# a second checkout of a different repository, not part of this one.
: "${BENCH_PAPER_WORK:=$PHYSICS24/faser/faser-evgen-bench-paper}"

# --- working names: assigned, never defaulted (see header) -----------------
SHERPA_RUNS=$BENCH_SHERPA_RUNS
SHERPA_INSTALL=$BENCH_SHERPA_INSTALL
SHERPA_NPATCH=$BENCH_SHERPA_NPATCH
HERWIG_BIN=$BENCH_HERWIG_BIN
GENIE_DIR=$BENCH_GENIE_DIR
PYTHIA6_LIB=$BENCH_PYTHIA6_LIB
NUWRO_DIR=$BENCH_NUWRO_DIR
GIBUU_DIR=$BENCH_GIBUU_DIR
BUUINPUT=$BENCH_BUUINPUT
ROOTEGPythia6_ROOT=$BENCH_ROOTEGPYTHIA6
APFEL_DIR=$BENCH_APFEL_DIR
MG5_DIR=$BENCH_MG5_DIR
PYTHIA8_DIR=$BENCH_PYTHIA8_DIR
POWHEG_RES=$BENCH_POWHEG_RES
POWHEG_RES_MUFLUX=$BENCH_POWHEG_RES_MUFLUX
POWHEG_V2=$BENCH_POWHEG_V2
REPORT_COPY=$BENCH_REPORT_COPY
FASER_GENIE=$BENCH_FASER_GENIE
FASER_GENIE_BRANCH=$BENCH_FASER_GENIE_BRANCH
FASER_DATA=$BENCH_FASER_DATA
GENIE_FASER_DIR=$BENCH_GENIE_FASER_DIR
PYTHIA8DATA=$PYTHIA8_DIR/share/Pythia8/xmldoc
LHAPDF_DATA_PATH=$BENCH_DEPS/share/LHAPDF

# --- platform ---------------------------------------------------------------
# The benchmark runs on this macOS laptop AND on the Nikhef Linux cluster, so
# everything that differs between them is decided ONCE, here.

# Runtime library search path: Linux and macOS spell it differently.
case "$(uname -s)" in
  Darwin) BENCH_LIBPATH_VAR=DYLD_LIBRARY_PATH ;;
  *)      BENCH_LIBPATH_VAR=LD_LIBRARY_PATH   ;;
esac

# Keep the machine awake for the length of a long job.  macOS needs
# `caffeinate -i`; Linux has no equivalent and needs nothing.
#
# USE THE FUNCTION `bench_run`, NOT THE VARIABLE.  An unquoted $BENCH_NOSLEEP
# relies on the shell word-splitting it into "caffeinate" + "-i", which bash
# does and ZSH DOES NOT -- under zsh it becomes a single command name and every
# job dies with `command not found: caffeinate -i`.  config.sh is sourced from
# zsh (it is the interactive shell here), so the variable form is a trap.
# The function behaves identically in sh, bash and zsh.
#
# Probed rather than assumed from uname: caffeinate is absent in some macOS
# environments, and a missing command would fail every job under `set -e`.
# "${VAR+set}" distinguishes UNSET from SET-BUT-EMPTY, so that an explicit
# BENCH_NOSLEEP="" is honoured (that is how you turn the wrapper off on a
# machine that does have caffeinate) while an unset one is probed.  A plain
# ":=" default cannot express that, and an unconditional assignment would
# silently ignore the override.
if [ "${BENCH_NOSLEEP+set}" != set ]; then
    if command -v caffeinate >/dev/null 2>&1; then
        BENCH_NOSLEEP="caffeinate -i"
    else
        BENCH_NOSLEEP=""
    fi
fi

# Run a command, wrapped in the keep-awake prefix when there is one.
#   bench_run ./main_dis card.cmnd events 1000 42
bench_run () {
    if [ -n "$BENCH_NOSLEEP" ]; then
        # `eval set --` re-parses the prefix, which splits "caffeinate -i" into
        # two words in EVERY shell.  A bare `bench_run "$@"` would rely on
        # the shell doing that splitting, which zsh does not do at all -- see
        # the note above.  The caller's own arguments stay quoted and so
        # survive spaces intact.
        eval 'set -- '"$BENCH_NOSLEEP"' "$@"'
    fi
    "$@"
}

# ANALYSE-AND-DELETE.  With BENCH_PRUNE=1 a run script hands the finished
# sample to tools/analyse_and_prune.py, which analyses it and then deletes the
# event files -- but only after every result is written, non-empty and newer
# than the sample.  122 GB of HepMC backs ~20 MB of JSON here, so on a cluster
# the sample should not outlive its own analysis.  Off by default: on this
# laptop the samples are still live inputs for the FASER selection, which has
# to re-parse them.  BENCH_PRUNE_KEEP_ONE=1 leaves the first job dir.
: "${BENCH_PRUNE:=0}"
: "${BENCH_PRUNE_KEEP_ONE:=0}"

# Shell used for the inner `-c` wrappers in genie/run_genie*.sh.  bash is the
# one interactive-grade shell present on both platforms -- the cluster nodes
# are not guaranteed to have zsh, which these scripts used to hardcode.
: "${BENCH_SHELL:=bash}"

# Default C++ compiler (config.mk reads this).  Apple clang here, GNU on Linux.
case "$(uname -s)" in
  Darwin) : "${BENCH_CXX:=clang++}" ;;
  *)      : "${BENCH_CXX:=g++}"     ;;
esac

export PHYSICS24 SOFTWARE BENCH BENCH_DEPS BENCH_LOG4CPP BENCH_EXTRA_BIN BENCH_PYTHON
export SHERPA_RUNS SHERPA_INSTALL SHERPA_NPATCH BENCH_SHERPA_NPATCH HERWIG_BIN GENIE_DIR PYTHIA6_LIB APFEL_DIR MG5_DIR
export NUWRO_DIR ROOTEGPythia6_ROOT GIBUU_DIR BUUINPUT
export PYTHIA8_DIR POWHEG_RES POWHEG_RES_MUFLUX POWHEG_V2 REPORT_COPY PYTHIA8DATA
export LHAPDF_DATA_PATH BENCH_LIBPATH_VAR BENCH_REPORT_COPY \
    BENCH_PAPER_REMOTE BENCH_PAPER_WORK BENCH_NUWRO_DIR BENCH_ROOTEGPYTHIA6 \
    BENCH_GIBUU_DIR BENCH_BUUINPUT
export BENCH_NOSLEEP BENCH_SHELL BENCH_CXX BENCH_PRUNE BENCH_PRUNE_KEEP_ONE
export FASER_GENIE FASER_GENIE_BRANCH FASER_DATA GENIE_FASER_DIR

# The variables tools/print_config.sh exposes to make and to Python.
BENCH_CONFIG_VARS="PHYSICS24 SOFTWARE BENCH BENCH_DEPS BENCH_LOG4CPP BENCH_PYTHON \
BENCH_EXTRA_BIN SHERPA_RUNS SHERPA_INSTALL SHERPA_NPATCH HERWIG_BIN GENIE_DIR PYTHIA6_LIB APFEL_DIR \
MG5_DIR PYTHIA8_DIR POWHEG_RES POWHEG_RES_MUFLUX POWHEG_V2 REPORT_COPY PYTHIA8DATA \
NUWRO_DIR ROOTEGPythia6_ROOT GIBUU_DIR BUUINPUT \
LHAPDF_DATA_PATH BENCH_LIBPATH_VAR BENCH_NOSLEEP BENCH_SHELL \
BENCH_CXX BENCH_PRUNE \
BENCH_PRUNE_KEEP_ONE FASER_GENIE FASER_GENIE_BRANCH FASER_DATA GENIE_FASER_DIR"
