#!/bin/sh
# Print the benchmark path configuration as KEY=VALUE lines.
#
# This is how non-shell consumers read config.sh, so that config.sh stays the
# single definition of every default:  config.mk ($(shell ...)) and
# analysis/paths.py (subprocess) both go through here.
#
#   ./tools/print_config.sh          # all variables
#   ./tools/print_config.sh SHERPA_RUNS PYTHIA8_DIR   # just these
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}   # a clone works from anywhere
export BENCH_REPO
. "$HERE/../config.sh"
for v in ${*:-$BENCH_CONFIG_VARS}; do
    eval "printf '%s=%s\n' \"\$v\" \"\${$v-}\""
done
