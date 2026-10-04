#!/usr/bin/env bash
# SEVEN-POINT SCALE WEIGHTS FOR THE  1 TeV POWHEG SAMPLES, for the hadron-level
# MHOU bands of paper plots 8, 9 and 11 (analysis/mhou_hadron.py).
#
# Usage:  powheg/production/reweight_scale.sh <mu|nu> <p|n> [job ...] [max_parallel=4]
#           mu: POWHEG-RES, one pass per SEED -- the job numbers ARE the seeds
#               (default: the first 10 complete jobs, the QED/Herwig subset)
#           nu: POWHEG-V2, one pass per BATCH (job 1 = the base directory,
#               job j = <rundir>-b<j>; default: job 1)
#
# The hard process is not rerun: `storeinfo_rwgt 1` (every card sets it)
# keeps what pwhg_main needs to re-evaluate each event at another scale, and
# `rwl_add 1` appends the weights of powheg/reweight/rwl_scale.xml.  Each pass
# works in its OWN directory, $root/v2rwgt-<cur>1TeV-<t>/<job>/, holding
# symlinks to the integration grids and a COPY-free read of the LHE the
# production job actually showered (its shower.log names it), and the
# weights are harvested at once into
#     results{,_nu}/mhou_weights/scale_<cur>_<t>_job<N>.npz
# after which the reweighted LHE is deleted (harvest_weights.py's own rule:
# the numbers are kept, the duplicated event records are not).
#
# >>> POWHEG-RES REWEIGHTS PER SEED (new, 2026-09-15). <<<  the earlier production's
# run_reweight.sh refused RES because `manyseeds` writes one LHE per seed.
# pwhg_main in manyseeds mode reads the seed index from stdin and then, with
# rwl_add, ASKS for the event file name on stdin too (it looks for a bare
# pwgevents.lhe otherwise); it writes pwgevents-rwgt-<NNNN>.lhe, ~4 min/50k:
# one pass per seed is exactly one pass per showered job, and the join to the
# shower is then per job as well (lhe_index counts within its own file).
#
# Safe to re-invoke: a job whose .npz exists is skipped.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
RWL=$BENCH_REPO/powheg/reweight/rwl_scale.xml

one () {   # $1 cur, $2 t, $3 job
    local cur=$1 t=$2 n=$3 job lhe rd dest out npz resdir pwhg seedidx
    job=$V2_POWHEG_DIR/v2${cur}_${t}_job_$n
    v2_shower_done "$job" || { echo "  $job is not a complete job" >&2; return 1; }
    lhe=$(grep -o 'Beams:LHEF *| *[^ |]*' "$job/shower.log" | awk '{print $NF}')
    [ -s "$lhe" ] || { echo "  no LHE for $job ($lhe)" >&2; return 1; }
    rd=$(dirname "$lhe")
    resdir=$BENCH_REPO/results_nu; [ "$cur" = mu ] && resdir=$BENCH_REPO/results
    mkdir -p "$resdir/mhou_weights"
    npz=$resdir/mhou_weights/scale_${cur}_${t}_job$n.npz
    [ -s "$npz" ] && { echo "  skip $cur $t job $n (weights exist)"; return 0; }
    if [ "$cur" = mu ]; then
        root=$POWHEG_RES; pwhg=$POWHEG_RES/pwhg_main
    else
        root=$POWHEG_V2; pwhg=$POWHEG_V2/pwhg_main
    fi
    dest=$root/v2rwgt-${cur}1TeV-$t/job$n
    rm -rf "$dest"; mkdir -p "$dest"
    # every integration file by symlink (grids, upper bounds, seeds, flavours);
    # the input card and the LHE are the only things the pass may rewrite
    for f in "$rd"/*; do
        case "$(basename "$f")" in
            powheg.input|pwgevents*|*.log|DEAD_SEEDS|V2_SUMMARY.txt|INTEGRATED) ;;
            *) ln -s "$f" "$dest/" ;;
        esac
    done
    cp "$rd/powheg.input" "$dest/"
    grep -q '^storeinfo_rwgt 1' "$dest/powheg.input" || {
        echo "  $rd/powheg.input has no storeinfo_rwgt 1 -- cannot reweight" >&2; return 1; }
    cp "$RWL" "$dest/"
    cat >> "$dest/powheg.input" <<KEYS

! ---- a-posteriori scale reweighting (powheg/production/reweight_scale.sh) ----
rwl_file 'rwl_scale.xml'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
    if [ "$cur" = mu ]; then
        seedidx=$(basename "$lhe" | sed -E 's/pwgevents-0*([0-9]+)\.lhe/\1/')
        ln -s "$lhe" "$dest/$(basename "$lhe")"
        ( cd "$dest" && printf "%s\n%s\n" "$seedidx" "$(basename "$lhe")" | bench_run "$pwhg" > reweight.log 2>&1 ) || {
            echo "  pwhg_main failed for $cur $t job $n (see $dest/reweight.log)" >&2; return 1; }
        out=$(ls "$dest"/pwgevents-rwgt-*.lhe 2>/dev/null | head -1)
    else
        ln -s "$lhe" "$dest/pwgevents.lhe"
        ( cd "$dest" && bench_run "$pwhg" > reweight.log 2>&1 ) || {
            echo "  pwhg_main failed for $cur $t job $n (see $dest/reweight.log)" >&2; return 1; }
        out=$dest/pwgevents-rwgt.lhe
    fi
    [ -s "$out" ] || { echo "  no reweighted LHE in $dest" >&2; ls "$dest" | head >&2; return 1; }
    python3 "$BENCH_REPO/powheg/reweight/harvest_weights.py" "$out" "$npz" > "$dest/harvest.log" 2>&1 || {
        echo "  harvest failed for $out" >&2; return 1; }
    rm -f "$out"
    echo "  ok $cur $t job $n -> $(basename "$npz")"
}

if [ "$1" = "--one" ]; then shift; one "$@"; exit $?; fi
cur=${1:?usage: reweight_scale.sh <mu|nu> <p|n> [job ...] [max_parallel]}
t=${2:?target}
shift 2
maxpar=${MAXPAR:-4}
jobs=()
for a in "$@"; do jobs+=("$a"); done
if [ ${#jobs[@]} -eq 0 ]; then
    if [ "$cur" = mu ]; then
        while read -r j; do jobs+=("$j"); done < <(
            for d in "$V2_POWHEG_DIR/v2mu_${t}_job"_*; do
                r=${d##*_}; case "$r" in (''|*[!0-9]*) continue ;; esac   # '(' form: bash 3.2 cannot parse a bare ')' pattern inside <( )
                v2_shower_done "$d" && echo "$r"
            done | sort -n | head -10)
    else
        jobs=(1)
    fi
fi
echo "scale reweighting $cur $t: jobs ${jobs[*]}"
for j in "${jobs[@]}"; do echo "$cur $t $j"; done | xargs -n 3 -P "$maxpar" "$HERE/reweight_scale.sh" --one
echo "== reweight_scale $cur $t finished =="
