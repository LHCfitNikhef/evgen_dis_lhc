#!/usr/bin/env bash
# POWHEG-RES muon NC for the paper-plots production: integrate, generate
# seed by seed, shower, until the showered target is met.  One (energy,
# nucleon) per invocation.
#
# Usage: powheg/production/run_res.sh <energy_gev> <p|n>
#
#   powheg/production/run_res.sh 1000 n      # -> $POWHEG_RES/v2-mu1TeV-n,
#                                    #    powheg/v2mu_n_job_<seed>
#
# Environment knobs (all optional):
#   V2_PAR        processes at once (default 3; the whole production may use 4)
#   V2_INTSEEDS   seeds used for stages 1-3 (default 8)
#   V2_TARGET     showered events to reach (default: PRODUCTION.md)
#   V2_CARD / V2_RUNDIR / V2_JOBBASE / V2_JOBROOT / V2_SEEDBASE
#                 overrides used by powheg/production/validate.sh ONLY
#
# RE-INVOCABLE.  Everything already done is skipped: the integration (marker
# INTEGRATED), every complete LHE (closing tag present), every complete
# shower job (events_xsec.json written), every seed recorded in DEAD_SEEDS.
# An LHE left incomplete by a KILLED run is regenerated -- after deleting it,
# because POWHEG refuses to overwrite pwgevents-NNNN.lhe and RETURNS 0.
#
# THE STAGES.  Stages 1-3 (grids and upper bounds) run once, on seeds
# 1..V2_INTSEEDS, V2_PAR at a time; the integrated cross-section is written by
# stage 3 into every pwg-NNNN-st3-stat.dat ("grand total total (pos.-|neg.|)",
# the same number in each file -- checked below), which is what the analysis
# reads.  Stage 4 then runs seed by seed, INCLUDING seeds that took no part
# in stages 1-3: stage 4 only loads the combined grids.  (Validated in
# validate.sh by generating a seed beyond the integration set.)
#
# DEAD SEEDS ARE DROPPED.  POWHEG-RES aborts a seed outright when the ISR
# radiation upper bound is exceeded (gen_radiation.f, "upper bound lower than
# actual value", exit on any violation with t > 1.5 rad_ptsqmin), and it does
# so REPRODUCIBLY TO THE EVENT: re-running the seed dies at the same place
# (memory: powheg-truncated-samples).  So a dead seed is recorded in
# DEAD_SEEDS, its partial LHE is deleted without being showered, and the
# next seed is taken.  Dropping whole independent streams is unbiased.
#
# EVERY LHE PASSES powheg/strip_lhe_nan.py BEFORE SHOWERING: a NaN
# four-momentum is read by Pythia as end-of-file and silently truncates the
# job.  The script also writes pwgevents-NNNN.lhe.nevents, the count the
# closure gate divides by.  Seed N showers into <jobbase>_N, the convention
# analyze.lhe_events_offered_for relies on.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"

E=${1:?usage: run_res.sh <energy_gev> <p|n>}
T=${2:?usage: run_res.sh <energy_gev> <p|n>}
IH2=$(v2_target_ih2 "$T")
TAG=$(v2_tag mu "$E")
CARD=${V2_CARD:-$BENCH_REPO/powheg/cards/production/POWHEG-RES/mu$TAG-$T/powheg.input}
RUNDIR=${V2_RUNDIR:-$POWHEG_RES/v2-mu$TAG-$T}
JOBROOT=${V2_JOBROOT:-$V2_POWHEG_DIR}
JOBBASE=${V2_JOBBASE:-$(v2_jobbase mu "$T" "$E")}
TARGET=${V2_TARGET:-$(v2_target_events "$E")}
PAR=${V2_PAR:-3}
NINT=${V2_INTSEEDS:-8}
MAXSEED=400
SEEDBASE=${V2_SEEDBASE:-$(v2_seed_base "$T" "$E")}
PWHG=$POWHEG_RES/pwhg_main
CMND=$(v2_cmnd mu)

[ -x "$PWHG" ] || { echo "no POWHEG-RES binary at $PWHG" >&2; exit 1; }
[ -f "$CARD" ] || { echo "no card $CARD (powheg/production/make_cards.py)" >&2; exit 1; }
grep -Eq "^ih2[[:space:]]+$IH2([[:space:]]|$)" "$CARD" || {
    echo "card $CARD does not set ih2 $IH2 for target $T -- refusing" >&2; exit 1; }
NEV=$(awk '$1=="numevts" {print $2; exit}' "$CARD")

log () { echo "[$(date '+%F %T')] RES mu $TAG $T: $*"; }
log "card $CARD"
log "run dir $RUNDIR, samples $JOBROOT/${JOBBASE}_<seed>, target $TARGET showered events, $NEV per seed, $PAR at a time"

mkdir -p "$RUNDIR"

# pwgseeds.dat: line i holds the random seed of run index i.  Written once;
# a mismatch against this production's seed block is fatal, because it would
# mean the directory holds some other run.
want_first=$((SEEDBASE + 1))
if [ -f "$RUNDIR/pwgseeds.dat" ]; then
    got_first=$(head -1 "$RUNDIR/pwgseeds.dat" | tr -d ' ')
    [ "$got_first" = "$want_first" ] || {
        echo "$RUNDIR/pwgseeds.dat starts at $got_first, want $want_first" >&2
        exit 1; }
else
    seq "$want_first" "$((SEEDBASE + MAXSEED))" > "$RUNDIR/pwgseeds.dat"
fi

set_stage () {   # $1 parallelstage, $2 xgriditeration
    python3 "$HERE/set_keys.py" "$CARD" "$RUNDIR/powheg.input" \
        "parallelstage=$1" "xgriditeration=$2"
}

run_seed () {    # $1 label, $2 seed  (foreground)
    ( cd "$RUNDIR" && echo "$2" | bench_run "$PWHG" \
        > "run_$1_$2.log" 2>&1 ) || true
}

# ------------------------------------------------------------ stages 1-3
if [ -f "$RUNDIR/INTEGRATED" ]; then
    log "integration already done: $(cat "$RUNDIR/INTEGRATED")"
else
    v2_disk_check 1 "RES integration $TAG $T"
    for st in "st1xg1 1 1 xg1" "st1xg2 1 2 xg2" "st2 2 1 st2" "st3 3 1 st3"; do
        set -- $st
        label=$1; set_stage "$2" "$3"; statkey=$4
        log "stage $label on seeds 1..$NINT"
        echo "beg $label $(date)" >> "$RUNDIR/Timings.txt"
        for s in $(seq 1 "$NINT"); do
            v2_throttle "$PAR"
            run_seed "$label" "$s" &
        done
        wait
        echo "end $label $(date)" >> "$RUNDIR/Timings.txt"
        # POWHEG exits 0 on some failures: check the product, not the status
        for s in $(seq 1 "$NINT"); do
            f=$RUNDIR/pwg-$(printf %04d "$s")-$statkey-stat.dat
            [ -s "$f" ] || {
                echo "stage $label seed $s produced no $f -- see" \
                     "$RUNDIR/run_${label}_$s.log" >&2; exit 1; }
        done
    done
    # every seed's st3 file must carry the SAME combined total -- compared to
    # 10 significant digits, because the seeds sum the error in different
    # orders (400 GeV, first production: 182.41792809074411 vs ...408)
    ntot=$(grep -h "grand total total" "$RUNDIR"/pwg-*-st3-stat.dat \
        | awk '{printf "%.10g %.10g\n", $(NF-2), $NF}' | sort -u | wc -l | tr -d ' ')
    [ "$ntot" = 1 ] || {
        echo "the st3 stat files disagree on the grand total ($ntot values)" >&2
        exit 1; }
    grep -h "grand total total" "$RUNDIR/pwg-0001-st3-stat.dat" \
        | sed 's/^ *//' > "$RUNDIR/INTEGRATED"
    log "integrated: $(cat "$RUNDIR/INTEGRATED")"
fi
[ -s "$RUNDIR/pwg-0001-st3-stat.dat" ] || {
    echo "no $RUNDIR/pwg-0001-st3-stat.dat -- the analysis needs it" >&2; exit 1; }

# ------------------------------------------------------------ stage 4 + shower
touch "$RUNDIR/DEAD_SEEDS"
is_dead () { grep -q "^seed $1 " "$RUNDIR/DEAD_SEEDS"; }

delivered () { v2_delivered "$JOBROOT" "$JOBBASE"; }   # see common.sh

set_stage 4 1
next=1
while :; do
    have=$(delivered)
    [ "$have" -ge "$TARGET" ] && break
    # how many seeds the remaining events need, at most PAR, assuming 3% loss
    need=$(( (TARGET - have + NEV * 97 / 100 - 1) / (NEV * 97 / 100) ))
    [ "$need" -gt "$PAR" ] && need=$PAR
    batch=()
    while [ "${#batch[@]}" -lt "$need" ] && [ "$next" -le "$MAXSEED" ]; do
        s=$next; next=$((next + 1))
        is_dead "$s" && continue
        v2_shower_done "$JOBROOT/${JOBBASE}_$s" && continue
        batch+=("$s")
    done
    [ "${#batch[@]}" -gt 0 ] || {
        echo "ran out of seeds (MAXSEED $MAXSEED) at $have of $TARGET events" >&2
        exit 1; }
    # LHE ~4.6 kB/event + HepMC ~5.2 kB/event, rounded up
    gb=$(( ${#batch[@]} * NEV * 10 / 1000000 + 1 ))
    v2_disk_check "$gb" "RES st4+shower $TAG $T seeds ${batch[*]}"
    log "have $have/$TARGET; stage 4 on seeds ${batch[*]}"
    for s in "${batch[@]}"; do
        lhe=$RUNDIR/pwgevents-$(printf %04d "$s").lhe
        v2_lhe_complete "$lhe" && continue
        rm -f "$lhe" "$lhe.nevents"
        run_seed st4 "$s" &
    done
    wait
    for s in "${batch[@]}"; do
        lhe=$RUNDIR/pwgevents-$(printf %04d "$s").lhe
        if v2_lhe_complete "$lhe"; then
            "$BENCH_REPO/powheg/strip_lhe_nan.py" "$lhe"
        else
            nw=$(grep -c "<event>" "$lhe" 2>/dev/null); nw=${nw:-0}
            why=$(grep -m1 "upper bound lower than actual value" "$RUNDIR/run_st4_$s.log" 2>/dev/null | tr -s ' ' | cut -c1-120)
            echo "seed $s died after $nw events: ${why:-no ISR-bound message, see run_st4_$s.log}" >> "$RUNDIR/DEAD_SEEDS"
            log "seed $s DEAD after $nw events (dropped): ${why:-see run_st4_$s.log}"
            rm -f "$lhe"
        fi
    done
    for s in "${batch[@]}"; do
        lhe=$RUNDIR/pwgevents-$(printf %04d "$s").lhe
        [ -s "$lhe" ] || continue
        d=$JOBROOT/${JOBBASE}_$s
        rm -rf "$d"
        v2_throttle "$PAR"
        v2_shower "$lhe" "$d" "$CMND" "$T" 13 &
    done
    wait
    for s in "${batch[@]}"; do
        d=$JOBROOT/${JOBBASE}_$s
        [ -s "$RUNDIR/pwgevents-$(printf %04d "$s").lhe" ] || continue
        v2_shower_done "$d" || { echo "shower of seed $s failed ($d)" >&2; exit 1; }
    done
done

ndead=$(grep -c "^seed " "$RUNDIR/DEAD_SEEDS" || true)
{
    echo "POWHEG-RES mu NC $TAG target $T (ih2 $IH2)"
    echo "integrator: $(cat "$RUNDIR/INTEGRATED")  [pwg-0001-st3-stat.dat]"
    echo "showered events: $(delivered) (target $TARGET) in $JOBROOT/${JOBBASE}_<seed>"
    echo "dead seeds dropped: ${ndead:-0}"
    echo "finished $(date)"
} > "$RUNDIR/V2_SUMMARY.txt"
log "DONE: $(tr '\n' ';' < "$RUNDIR/V2_SUMMARY.txt")"
