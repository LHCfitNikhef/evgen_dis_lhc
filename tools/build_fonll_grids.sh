#!/usr/bin/env bash
# Pre-build the PineAPPL grids FONLL is assembled from, several at a time.
#
# Usage: tools/build_fonll_grids.sh [mu|nu|both] [max_parallel]
#
# WHY PRE-BUILD.  FONLL is  F = ZM + f_thr*(FFNS - FFN0), so a PDF study in
# FONLL needs THREE grids per observable instead of one.  The ZM run is a
# couple of minutes; the two MASSIVE runs are the expensive ones -- the
# coefficient functions are numerically integrated -- and they are what made
# the direct FONLL calculation a 6.75 h job before grids existed.
#
# But they are PDF-INDEPENDENT and beam-energy-specific, so they are built
# ONCE per (current, energy, scheme) and then convoluted with every member of
# every set at ~5 ms each.  Building them here, in parallel, turns a long
# serial wait inside pdf_dependence.py into one batch that uses the cores.
#
# THE GRIDS ARE A CACHE, NOT A RESULT (see analysis/pineappl_grids.py): they
# are keyed on current, energy, order, SCHEME and a hash of the node grid, and
# can always be rebuilt.  Re-running this is therefore always safe -- an
# already-built grid is skipped, not recomputed.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"

WHICH=${1:-both}
MAXPAR=${2:-5}
case "$WHICH" in
  mu|nu) CURRENTS="$WHICH" ;;
  both)  CURRENTS="nu mu" ;;
  *)     echo "usage: $0 [mu|nu|both] [max_parallel]" >&2; exit 1 ;;
esac

mkdir -p "$BENCH_REPO/logs"
n=0
for cur in $CURRENTS; do
  for e in 400 1000 4000; do
    # ZM first: it is cheap and every scheme needs it
    for spec in "ZM-VFNS 5" "FONLL-FFNS 3" "FONLL-FFN0 3"; do
      set -- $spec
      fns=$1; nfff=$2
      log="$BENCH_REPO/logs/grid_${cur}_${e}_$(echo "$fns" | tr -d '-').log"
      ( cd "$BENCH_REPO" && python3 analysis/pineappl_grids.py \
          build "$cur" "$e" 1 "$fns" "$nfff" > "$log" 2>&1 ) &
      n=$((n + 1))
      echo "  launched $cur $e GeV $fns  -> $(basename "$log")"
      # bash 3.2 has no `wait -n`, so throttle in whole batches
      if [ $((n % MAXPAR)) -eq 0 ]; then
        echo "  ... waiting for this batch of $MAXPAR"
        # `wait` UNDER `set -e` KILLS THIS SCRIPT if any child failed -- or was
        # killed by hand, which is how the first run of this died after one
        # batch with a bare "Terminated: 15" and no indication that ten builds
        # had never been launched.  A failed grid is worth reporting and worth
        # continuing past: the others are independent, and a missing grid is
        # caught later by the build that needs it.
        wait || echo "  (a build in that batch failed -- see its log)"
      fi
    done
  done
done
wait || echo "  (a build in the last batch failed -- see its log)"
echo "all $n grid build(s) finished"
