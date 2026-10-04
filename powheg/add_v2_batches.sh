#!/usr/bin/env bash
# Add statistics to a POWHEG-V2 production WITHOUT touching the existing one.
#
# WHY NOT JUST RAISE numevts.  run_powheg_v2.sh runs the code ONCE into one
# run directory and showers one pwgevents.lhe into <jobbase>_1 -- POWHEG-V2's
# cards have no parallelstage, so that is the shape of the code, not a choice.
# Raising numevts therefore REGENERATES the anchor LHE, and the 1 TeV neutrino
# LHE is what the PDF-member reweighting is joined to by lhe_index: rwgtnu_job_
# 1TeV_1 and rwgt-nu1TeV/weights.npz both point at it.  Replacing it silently
# invalidates that join -- the exact failure in the lhe-index episode, where a
# sample showered two minutes before a fix scrambled every reweighting number.
#
# So instead this ADDS independent batches: a new run directory per batch,
# carrying the SAME cached grids (use-old-grid / use-old-ubound are already 1
# in the cards) and a DIFFERENT iseed.  A fixed importance-sampling grid with
# independent seeds is exactly what POWHEG's own manyseeds mode does, so the
# batches are statistically independent draws from the same integral, and the
# integrator cross-section in pwg-stat.dat is unchanged and stays the right
# normalisation for all of them.
#
# The batches shower into <jobbase>_2, _3, ... which analyze_nu.py already
# pools: job_files() matches "<base>_<N>" with N an exact integer, so the new
# directories join the sample by existing under the name the analysis expects,
# and pooled_sigma_pb() combines them by accepted count.  Nothing that reads
# the old sample changes.
#
# Usage: ./add_v2_batches.sh [--mc] <ebeam> <first_job> <nbatches> <nev_each>
#
#   ./add_v2_batches.sh 1000 2 5 170000     # nu_job_2..6, 170k events each
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/powheg_variants.sh"

PV_CODE=V2
NEWARGS=()
for a in "$@"; do
    if [ "$a" = "--mc" ]; then PV_CODE=V2MC; else NEWARGS+=("$a"); fi
done
set -- "${NEWARGS[@]}"

powheg_parse_current nu "$@"
set -- "${PV_ARGS[@]}"
EBEAM=${1:?usage: add_v2_batches.sh [--mc] <ebeam> <first_job> <nbatches> <nev_each>}
FIRSTJOB=${2:?need a first job number}
NBATCH=${3:?need a batch count}
NEVEACH=${4:?need events per batch}

powheg_variant "$PV_CODE" "$PV_CURRENT" "$EBEAM"
SRCDIR=$PV_RUNDIR
JOBBASE=$PV_JOBBASE
PWHG=$POWHEG_V2/pwhg_main
[ -x "$PWHG" ] || { echo "no POWHEG-V2 binary at $PWHG" >&2; exit 1; }
[ -s "$SRCDIR/pwggrid.dat" ] || {
    echo "no cached grid in $SRCDIR -- run run_powheg_v2.sh $EBEAM first" >&2
    exit 1; }

# Refuse to clobber, per job AND per run directory, before anything runs.
for i in $(seq 0 $((NBATCH - 1))); do
    j=$((FIRSTJOB + i))
    [ -s "$HERE/${JOBBASE}_$j/events.hepmc" ] && [ "${FORCE:-0}" != "1" ] && {
        echo "refusing to overwrite $HERE/${JOBBASE}_$j/events.hepmc -- FORCE=1" >&2
        exit 1; }
done

echo "$PV_LABEL: adding $NBATCH batch(es) of $NEVEACH events at $EBEAM GeV"
echo "  grids cached from $SRCDIR; existing sample untouched"

pids=()
for i in $(seq 0 $((NBATCH - 1))); do
    j=$((FIRSTJOB + i))
    rd=$SRCDIR-b$j
    mkdir -p "$rd"
    # Everything the integration produced, so use-old-grid/use-old-ubound hit.
    # NOT the events or the logs: those are what this batch is here to make.
    for f in pwggrid.dat pwgubound.dat pwg-stat.dat pwgborngrid.top \
             FlavRegList bornequiv pwg-btlgrid.top pwg-rmngrid.top; do
        [ -f "$SRCDIR/$f" ] && cp -f "$SRCDIR/$f" "$rd/$f"
    done
    # A REAL FILE, not link_cards' symlink: the seed and the event count are
    # rewritten per batch and writing through the symlink would edit the
    # committed card (the same reason run_powheg_res.sh copies).
    python3 - "$SRCDIR/powheg.input" "$rd/powheg.input" "$NEVEACH" "$((12 + j))" <<'PY'
import re, sys
src, dst, nev, seed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
s = open(src).read()
s = re.sub(r"^numevts\s+\d+", f"numevts {nev}", s, flags=re.M)
s = re.sub(r"^iseed\s+\d+", f"iseed {seed}", s, flags=re.M)
open(dst, "w").write(s)
PY
    rm -f "$rd/pwgevents.lhe"
    ( cd "$rd" && bench_run "$PWHG" > run_v2.log 2>&1 ) &
    pids+=($!)
done
echo "  ${#pids[@]} generation process(es) running..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail + 1)); done
[ "$fail" -gt 0 ] && echo "  WARNING: $fail generation process(es) returned non-zero" >&2

# Shower each batch into its own job directory, in parallel.
for i in $(seq 0 $((NBATCH - 1))); do
    j=$((FIRSTJOB + i))
    rd=$SRCDIR-b$j
    [ -s "$rd/pwgevents.lhe" ] || { echo "  batch $j produced no LHE -- skipped" >&2; continue; }
    "$BENCH_REPO/powheg/strip_lhe_nan.py" "$rd/pwgevents.lhe"
    d=$HERE/${JOBBASE}_$j
    mkdir -p "$d"
    ( cd "$HERE" && bench_run ./main_powheg "$PV_CMND" \
        "$rd/pwgevents.lhe" "$d/events" > "$d/shower.log" 2>&1 ) &
done
wait

# DELIVERY GUARD, as in run_powheg_shower.sh: a shower job can fail on every
# event and still leave an empty events.hepmc behind an exit 0.
short=0
for i in $(seq 0 $((NBATCH - 1))); do
    j=$((FIRSTJOB + i))
    d=$HERE/${JOBBASE}_$j
    n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null); n=${n:-0}
    echo "  ${JOBBASE}_$j: $n events, $(cat "$d/events_xsec.json" 2>/dev/null)"
    [ "$n" -lt $((NEVEACH / 2)) ] && { echo "    *** under half of $NEVEACH ***" >&2; short=$((short+1)); }
done
[ "$short" -gt 0 ] && { echo "$short batch(es) short -- sample NOT usable" >&2; exit 1; }
echo "V2BATCHESDONE"
