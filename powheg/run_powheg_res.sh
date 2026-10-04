#!/usr/bin/env bash
# POWHEG-BOX-RES DIS at NLO: the full parallel run, all stages.
#
# Usage: ./run_powheg_res.sh [--current mu|nu] <beam_energy_gev> [nseeds|seed list]
#
#   ./run_powheg_res.sh 400                    # muon NC, 400 GeV, 10 seeds
#   ./run_powheg_res.sh 4000 10
#   ./run_powheg_res.sh --current nu 1000      # the neutrino CC cross-variant
#   ./run_powheg_res.sh 400 "11 12 13 14 15"   # ADD seeds 11-15 to what exists
#
# WHY A SEED LIST (2026-08-27).  POWHEG aborts a seed outright when the ISR
# radiation upper bound is exceeded -- gen_radiation.f calls exit(1) on any
# violation once t > 1.5*rad_ptsqmin, however tiny the violation is (0.1-0.2%
# in every case seen here).  It is deterministic: the same seed dies at the
# same event every time, so re-running a dead seed reproduces it exactly.  The
# events written before the abort are correctly distributed, so this costs
# STATISTICS and not accuracy, and the remedy is more seeds rather than
# repeated attempts at the same one.  Four of ten seeds died at 400 GeV and
# one of ten at 1 TeV; 4 TeV is clean.
#
# The t at which it happens sits just under the card's q2cut, so it is the
# generation cut the radiation is bumping into.  Raising q2cut would change
# the generation region and invalidate the comparison with every existing
# sample, so it is left alone.
#
# THE CURRENT IS A PARAMETER, THE CODE IS NOT.  This driver runs POWHEG-RES,
# whose native side is muon NC -- hence the default -- but the same five-stage
# sequence generates the neutrino CC cross-variant, whose only differences are
# the card and the run directory.  Those two live in powheg_variants.sh, which
# is also where the CC variant's INCLUSIVE-ONLY caveat is written down; it is
# echoed into the log below rather than left in a comment nobody reads.
#
# The 1 TeV run was driven BY HAND in August 2026 and never committed, which
# is exactly what CONVENTIONS.md rule 1 exists to prevent -- so the sequence is
# written down here.  From that run's Timings.txt the whole thing is about
# five minutes: st1 xg1 3 s, st1 xg2 4 s, st2 45 s, st3 11 s, st4 3.5 min.
#
# THE FIVE PASSES.  POWHEG-BOX-RES runs the same binary five times, each time
# with a different (parallelstage, xgriditeration) in powheg.input, and every
# seed is an independent process reading its index on stdin:
#
#   stage 1, xgriditeration 1   importance-sampling grid, first pass
#   stage 1, xgriditeration 2   the same grid, refined
#   stage 2                     upper bound for the btilde generation
#   stage 3                     upper bound for the radiation
#   stage 4                     event generation -> pwgevents-NNNN.lhe
#
# WHY THE CARD IS COPIED, NOT SYMLINKED.  tools/link_cards.sh normally makes
# the run directory's powheg.input a SYMLINK back into the repo, so there is
# one versioned copy.  This driver has to rewrite parallelstage between
# passes, and writing through that symlink would edit the committed card.  So
# it copies the card in as a real file and rewrites the copy.  After stage 4
# the copy is identical to the repo card again (stage 4 is the committed
# form), so link_cards.sh restores the symlink on its next run -- it reports
# that as "replace copy (identical to repo)".
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/powheg_variants.sh"     # powheg_variant / powheg_parse_current

powheg_parse_current mu "$@"     # default: POWHEG-RES's native current
set -- "${PV_ARGS[@]}"
EBEAM=${1:?usage: run_powheg_res.sh [--current mu|nu] <beam_energy_gev> [nseeds]}
NSEEDS=${2:-10}

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$PV_CURRENT" "$EBEAM")
powheg_variant RES "$PV_CURRENT" "$EBEAM"
CARD=$PV_CARD
RUNDIR=$PV_RUNDIR
PWHG=$POWHEG_RES/pwhg_main
[ -x "$PWHG" ] || { echo "no POWHEG binary at $PWHG" >&2; exit 1; }

# A seed LIST adds to what is already there; a seed COUNT means 1..N.
# NOTE "${SEEDLIST[@]}", never $SEEDLIST: bash 3.2 expands the latter to the
# first element only (CONVENTIONS.md rule 1 -- that once ran 1 GENIE job of 8).
case "$NSEEDS" in
    *' '*) SEEDLIST=($NSEEDS); ADDING=1 ;;
    *)     SEEDLIST=($(seq 1 "$NSEEDS")); ADDING=0 ;;
esac
MAXSEED=0
for s in "${SEEDLIST[@]}"; do [ "$s" -gt "$MAXSEED" ] && MAXSEED=$s; done

# Guard the seeds ACTUALLY BEING RUN, not seed 1.  When adding seeds 11+ to a
# finished production, pwgevents-0001.lhe exists by definition and a guard on
# it would refuse every time.
for s in "${SEEDLIST[@]}"; do
    f=$RUNDIR/pwgevents-$(printf "%04d" "$s").lhe
    if [ -s "$f" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "refusing to overwrite $f -- FORCE=1 to replace" >&2
        exit 1
    fi
done

mkdir -p "$RUNDIR"
powheg_install_card "$CARD" "$RUNDIR/powheg.input"   # a real file: see above

# POWHEG HAS ITS OWN CLOBBER GUARD, separate from this script's: stage 4 exits
# with "file pwgevents-0001.lhe exists! will not overwrite" and returns 0, so
# the stage LOOKS like it ran (one second, no error) and the previous LHE is
# silently reused.  That cost a full debugging cycle.  FORCE=1 means FORCE.
#
# ONLY THE SEEDS BEING RUN.  This used to be `rm -f pwgevents-*.lhe`, which
# with a seed list would delete a whole finished production to regenerate five
# of its files.
if [ "${FORCE:-0}" = "1" ]; then
    for s in "${SEEDLIST[@]}"; do
        rm -f "$RUNDIR/pwgevents-$(printf "%04d" "$s").lhe" \
              "$RUNDIR/pwgevents-$(printf "%04d" "$s").lhe.nevents"
    done
fi

# pwgseeds.dat must cover the highest seed.  Only ever GROWN: seq 1 N always
# reproduces lines 1..N, so extending it cannot change an existing seed's
# random stream, but shrinking it would strand finished seeds.
have=0
[ -f "$RUNDIR/pwgseeds.dat" ] && have=$(wc -l < "$RUNDIR/pwgseeds.dat" | tr -d ' ')
if [ "$MAXSEED" -gt "$have" ]; then
    seq 1 "$MAXSEED" > "$RUNDIR/pwgseeds.dat"
fi

# Rewrite the two staging keys in place.
set_stage () {   # $1 = parallelstage, $2 = xgriditeration
    python3 - "$RUNDIR/powheg.input" "$1" "$2" <<'PY'
import re, sys
path, stage, xg = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
s = re.sub(r"^parallelstage\s+\d+", f"parallelstage {stage}", s, flags=re.M)
s = re.sub(r"^xgriditeration\s+\d+", f"xgriditeration {xg}", s, flags=re.M)
open(path, "w").write(s)
PY
}

run_pass () {    # $1 = label, $2 = parallelstage, $3 = xgriditeration
    local label=$1
    set_stage "$2" "$3"
    echo "  $label: seeds ${SEEDLIST[*]} ..."
    echo "beg $label $(date)" >> "$RUNDIR/Timings.txt"
    for s in "${SEEDLIST[@]}"; do
        ( cd "$RUNDIR" && echo "$s" | bench_run "$PWHG" \
            > "run_${label}_$s.log" 2>&1 ) &
    done
    wait
    echo "end $label $(date)" >> "$RUNDIR/Timings.txt"
    # POWHEG exits 0 on some failures, so check the stage really produced
    # something rather than trusting the return code.
    if grep -qil "error\|not found\|abort" "$RUNDIR"/run_${label}_${SEEDLIST[0]}.log 2>/dev/null; then
        echo "  WARNING: $label log mentions an error -- check" \
             "$RUNDIR/run_${label}_${SEEDLIST[0]}.log" >&2
    fi
}

echo "$PV_LABEL at NLO: $EBEAM GeV [$TAG], seeds ${SEEDLIST[*]} -> $RUNDIR"
[ "$ADDING" = 1 ] && echo "  (adding to the seeds already there)"
run_pass st1xg1 1 1
run_pass st1xg2 1 2
run_pass st2    2 1
run_pass st3    3 1
run_pass st4    4 1

n=$(ls "$RUNDIR"/pwgevents-*.lhe 2>/dev/null | wc -l | tr -d ' ')
echo "done: $n LHE file(s) in $RUNDIR"
[ "$n" -gt 0 ] || { echo "NO LHE FILES PRODUCED -- check the stage logs" >&2; exit 1; }
if [ "$PV_CURRENT" = mu ]; then
    echo "next: powheg/run_powheg_shower.sh $EBEAM   (shower the LHE with Pythia8)"
else
    echo "next: powheg/run_powheg_shower.sh --current $PV_CURRENT $EBEAM   (shower the LHE with Pythia8)"
fi
