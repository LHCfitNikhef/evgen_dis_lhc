#!/usr/bin/env bash
# POWHEG-V2 neutrino CC for the paper-plots production: integrate and
# generate batch 1, add batches with the cached grids, shower, until the
# showered target is met.  One (energy, nucleon) per invocation.
#
# Usage: powheg/production/run_v2.sh <energy_gev> <p|n>
#
#   powheg/production/run_v2.sh 1000 n   # -> $POWHEG_V2/v2-nu1TeV-n[-b<j>],
#                                #    powheg/v2nu_n_job_<j>
#
# Environment knobs (all optional):
#   V2_PAR        processes at once for batches >= 2 and showers (default 1)
#   V2_TARGET     showered events to reach (default: PRODUCTION.md)
#   V2_CARD / V2_RUNDIR / V2_JOBBASE / V2_JOBROOT / V2_SEEDBASE
#                 overrides used by powheg/production/validate.sh ONLY
#
# THE SHAPE OF THE CODE.  POWHEG-V2's cards have no parallelstage: one process
# integrates and generates.  Batch 1 in <rundir> does both and is the ONLY
# place grids are made -- a fresh directory, so use-old-grid can never pick up
# another region's grids.  Batch j >= 2 in <rundir>-b<j> gets a copy of batch
# 1's grids and a different iseed, which is what add_v2_batches.sh does for the earlier production
# and what manyseeds does in POWHEG-RES.  A copied grid is CHECKED to have
# been loaded, not re-integrated: after the run the batch's pwggrid.dat and
# pwgubound.dat must still be byte-identical to batch 1's.  (Comparing
# pwg-stat.dat would prove nothing: a loading run does not rewrite it, so the
# copy always matches -- measured on the earlier production's prod-nu1TeV-ckm-b2, whose
# pwg-stat.dat keeps the copy's mtime.  It is therefore NOT copied.)
# <rundir>-b<j> is the sibling pattern analyze.lhe_events_offered globs.
#
# The integrated cross-section is batch 1's pwg-stat.dat,
# "total (btilde+remnants) cross section in pb".
#
# RE-INVOCABLE.  A batch whose LHE is complete is not regenerated; a shower
# job with V2_OK is not redone.  An incomplete batch-1 directory (a killed
# integration) is WIPED before rerunning, since a half-written grid under
# use-old-grid 1 is exactly the silent reuse this layout exists to prevent.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"

E=${1:?usage: run_v2.sh <energy_gev> <p|n>}
T=${2:?usage: run_v2.sh <energy_gev> <p|n>}
IH2=$(v2_target_ih2 "$T")
TAG=$(v2_tag nu "$E")
CARD=${V2_CARD:-$BENCH_REPO/powheg/cards/production/POWHEG-V2/nu$TAG-$T/powheg.input}
RUNDIR=${V2_RUNDIR:-$POWHEG_V2/v2-nu$TAG-$T}
JOBROOT=${V2_JOBROOT:-$V2_POWHEG_DIR}
JOBBASE=${V2_JOBBASE:-$(v2_jobbase nu "$T" "$E")}
TARGET=${V2_TARGET:-$(v2_target_events "$E")}
PAR=${V2_PAR:-1}
MAXBATCH=40
SEEDBASE=${V2_SEEDBASE:-$(v2_seed_base "$T" "$E")}
PWHG=$POWHEG_V2/pwhg_main
CMND=$(v2_cmnd nu)

[ -x "$PWHG" ] || { echo "no POWHEG-V2 binary at $PWHG" >&2; exit 1; }
[ -f "$CARD" ] || { echo "no card $CARD (powheg/production/make_cards.py)" >&2; exit 1; }
grep -Eq "^ih2[[:space:]]+$IH2([[:space:]]|$)" "$CARD" || {
    echo "card $CARD does not set ih2 $IH2 for target $T -- refusing" >&2; exit 1; }
NEV=$(awk '$1=="numevts" {print $2; exit}' "$CARD")

log () { echo "[$(date '+%F %T')] V2 nu $TAG $T: $*"; }
log "card $CARD"
log "run dir $RUNDIR[-b<j>], samples $JOBROOT/${JOBBASE}_<j>, target $TARGET showered events, $NEV per batch"

bdir () { if [ "$1" = 1 ]; then echo "$RUNDIR"; else echo "$RUNDIR-b$1"; fi; }

stat_total () {
    grep "total (btilde+remnants) cross section in pb" "$1/pwg-stat.dat" 2>/dev/null \
        | sed 's/^ *//'
}

run_batch () {   # $1 = batch index (foreground)
    local j=$1 rd
    rd=$(bdir "$j")
    ( cd "$rd" && bench_run "$PWHG" > run_v2.log 2>&1 ) || true
    v2_lhe_complete "$rd/pwgevents.lhe" || {
        echo "batch $j: no complete pwgevents.lhe -- see $rd/run_v2.log" >&2
        return 1; }
}

# ------------------------------------------------------------ batch 1
rd1=$(bdir 1)
if v2_lhe_complete "$rd1/pwgevents.lhe" && [ -n "$(stat_total "$rd1")" ]; then
    log "batch 1 already generated: $(stat_total "$rd1")"
else
    v2_disk_check $(( NEV * 12 / 1000000 + 1 )) "V2 batch 1 $TAG $T"
    rm -rf "$rd1"
    mkdir -p "$rd1"
    python3 "$HERE/set_keys.py" "$CARD" "$rd1/powheg.input" "iseed=$((SEEDBASE + 1))"
    log "batch 1: integrating and generating $NEV events (iseed $((SEEDBASE + 1)))"
    run_batch 1
    [ -n "$(stat_total "$rd1")" ] || { echo "no pwg-stat.dat total in $rd1" >&2; exit 1; }
    log "batch 1 integrated: $(stat_total "$rd1")"
fi
REF_TOTAL=$(stat_total "$rd1")

delivered () { v2_delivered "$JOBROOT" "$JOBBASE"; }   # see common.sh

prepare_batch () {   # $1 = j >= 2: fresh dir with batch 1's grids
    local j=$1 rd f
    rd=$(bdir "$j")
    rm -rf "$rd"
    mkdir -p "$rd"
    for f in pwggrid.dat pwgubound.dat pwgborngrid.top pwgxgrid.dat \
             FlavRegList bornequiv virtequiv realequivregions-btl \
             realequivregions-rad pwg-btlgrid.top pwg-rmngrid.top; do
        [ -f "$rd1/$f" ] && cp -f "$rd1/$f" "$rd/$f"
    done
    python3 "$HERE/set_keys.py" "$rd1/powheg.input" "$rd/powheg.input" \
        "iseed=$((SEEDBASE + j))"
}

# ------------------------------------------------------------ batches + showers
j=1
while :; do
    have=$(delivered)
    [ "$have" -ge "$TARGET" ] && break
    need=$(( (TARGET - have + NEV * 97 / 100 - 1) / (NEV * 97 / 100) ))
    [ "$need" -gt "$PAR" ] && need=$PAR
    batch=()
    while [ "${#batch[@]}" -lt "$need" ] && [ "$j" -le "$MAXBATCH" ]; do
        v2_shower_done "$JOBROOT/${JOBBASE}_$j" || batch+=("$j")
        j=$((j + 1))
    done
    [ "${#batch[@]}" -gt 0 ] || { echo "ran out of batches at $have/$TARGET" >&2; exit 1; }
    gb=$(( ${#batch[@]} * NEV * 12 / 1000000 + 1 ))    # LHE 5.7 + HepMC 6.4 kB/ev
    v2_disk_check "$gb" "V2 batches ${batch[*]} $TAG $T"
    log "have $have/$TARGET; batches ${batch[*]}"
    for b in "${batch[@]}"; do
        rd=$(bdir "$b")
        v2_lhe_complete "$rd/pwgevents.lhe" && continue
        [ "$b" = 1 ] && { echo "batch 1 LHE vanished" >&2; exit 1; }
        prepare_batch "$b"
        v2_throttle "$PAR"
        run_batch "$b" &
    done
    wait
    for b in "${batch[@]}"; do
        rd=$(bdir "$b")
        v2_lhe_complete "$rd/pwgevents.lhe" || { echo "batch $b failed" >&2; exit 1; }
        if [ "$b" != 1 ] && ! { cmp -s "$rd/pwggrid.dat" "$rd1/pwggrid.dat" \
                && cmp -s "$rd/pwgubound.dat" "$rd1/pwgubound.dat"; }; then
            echo "batch $b: pwggrid.dat/pwgubound.dat differ from batch 1's --" \
                 "it RE-INTEGRATED instead of loading the grids; refusing to" \
                 "pool it" >&2
            exit 1
        fi
        "$BENCH_REPO/powheg/strip_lhe_nan.py" "$rd/pwgevents.lhe"
    done
    for b in "${batch[@]}"; do
        d=$JOBROOT/${JOBBASE}_$b
        rm -rf "$d"
        v2_throttle "$PAR"
        v2_shower "$(bdir "$b")/pwgevents.lhe" "$d" "$CMND" "$T" 14 &
    done
    wait
    for b in "${batch[@]}"; do
        v2_shower_done "$JOBROOT/${JOBBASE}_$b" || {
            echo "shower of batch $b failed" >&2; exit 1; }
    done
done

{
    echo "POWHEG-V2 nu CC $TAG target $T (ih2 $IH2)"
    echo "integrator: $REF_TOTAL  [$(basename "$rd1")/pwg-stat.dat]"
    echo "showered events: $(delivered) (target $TARGET) in $JOBROOT/${JOBBASE}_<batch>"
    echo "finished $(date)"
} > "$rd1/V2_SUMMARY.txt"
log "DONE: $(tr '\n' ';' < "$rd1/V2_SUMMARY.txt")"
