#!/bin/sh
# Verify that every path config.sh names actually exists on this machine.
#
# Run this FIRST on a new machine (a cluster, a fresh clone).  The benchmark's
# recurring failure mode is a wrong path that does not error: GXMLPATH falling
# back to GENIE's built-in tune, PYTHIA8DATA quietly resolving to conda's
# 8.312, a dropped card symlink.  Those cost days each.  This costs a second.
#
#   ./tools/check_paths.sh          # report, exit 1 if anything is missing
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$HERE/../config.sh"

# Variables that must name an existing directory or file.  BENCH_LIBPATH_VAR
# is a variable NAME, so it is not checked.
CHECK="PHYSICS24 SOFTWARE BENCH BENCH_DEPS BENCH_LOG4CPP BENCH_EXTRA_BIN \
SHERPA_RUNS SHERPA_INSTALL HERWIG_BIN GENIE_DIR PYTHIA6_LIB APFEL_DIR MG5_DIR PYTHIA8_DIR \
POWHEG_RES POWHEG_V2 PYTHIA8DATA LHAPDF_DATA_PATH"

missing=0
for v in $CHECK; do
    eval "val=\${$v-}"
    if [ -z "$val" ]; then
        printf 'UNSET    %-18s\n' "$v"; missing=$((missing + 1))
    elif [ -e "$val" ]; then
        printf 'ok       %-18s %s\n' "$v" "$val"
    else
        printf 'MISSING  %-18s %s\n' "$v" "$val"; missing=$((missing + 1))
    fi
done

# THE ANALYSIS PYTHON AND ITS PACKAGES (environment.yml).  A missing module
# fails loudly anyway, but hours into a production rather than here.
printf '\n'
if "$BENCH_PYTHON" -c "import sys" 2>/dev/null; then
    printf 'ok       %-18s %s\n' BENCH_PYTHON "$BENCH_PYTHON"
    for m in numpy scipy matplotlib yaml yadism eko pineappl lhapdf ROOT; do
        if "$BENCH_PYTHON" -c "import $m" 2>/dev/null; then
            printf 'ok         python: %s\n' "$m"
        else
            printf 'MISSING    python: %s  (environment.yml)\n' "$m"
            missing=$((missing + 1))
        fi
    done
else
    printf 'MISSING  %-18s %s\n' BENCH_PYTHON "$BENCH_PYTHON"; missing=$((missing + 1))
fi
command -v latex > /dev/null 2>&1 && printf 'ok       latex (figures use usetex)\n' \
    || { printf 'MISSING  latex (figures use usetex, analysis/plotstyle.py)\n'; missing=$((missing + 1)); }

printf '\n'
if [ "$missing" -gt 0 ]; then
    printf '%s check(s) failed.\n' "$missing"
    printf 'Set PHYSICS24, or the individual BENCH_* knob, and re-run.\n'
    printf 'Knobs are listed at the top of config.sh.\n'
    exit 1
fi
printf 'all paths resolve.\n'
