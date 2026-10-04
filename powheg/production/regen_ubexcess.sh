#!/usr/bin/env bash
# REGENERATE THE  POWHEG-V2 AND POWHEG-V2mc SAMPLES WITH ubexcess_correct 1
# (user, 2026-09-19: "Lets do 1 now").
#
#   powheg/production/regen_ubexcess.sh [--list] [max_parallel=10]
#
# Why: without the switch POWHEG under-produces every point where the true
# weight exceeds its unweighting bound; for POWHEG-V2 at FASER energies those
# points cluster at low x and low Q2, and the tungsten x spectrum sat 6-9% low
# at x = 0.004-0.011 (pp02) while the total stayed right.  Evidence and the
# same-grid test (0.92-0.94 -> 0.98-1.02): PRODUCTION.md.
#
# What it does, per EXISTING shower job powheg/v2nu_<t>_job[_TAG]_<j> and
# powheg/v2numc_<t>_job[_TAG]_<j> (the job's shower.log names its LHE, so the
# run directory is read, never re-derived):
#   1. in the run directory: keep the integration grids and the upper bounds,
#      add `ubexcess_correct 1` to its powheg.input, keep its iseed -- the SAME
#      random stream as the original batch, so old and new are maximally
#      correlated and the difference is the correction alone;
#   2. regenerate the LHE (pwhg_main reads the grids: use-old-grid/-ubound 1)
#      and REFUSE if pwggrid.dat or pwgubound.dat changed (a re-integration);
#   3. strip NaN events (powheg/strip_lhe_nan.py), as run_v2.sh does;
#   4. re-shower into the SAME job directory with the production card
#      (v2_shower: beams checked, matchInOut off, V2_OK).
# A run directory that has done 1-3 carries UBEXCESS_OK; a job whose run dir
# has it and whose shower is complete and newer is skipped -- safe to kill and
# re-invoke.  CONVENTIONS.md 1c: do not edit while running.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
PWHG=$POWHEG_V2/pwhg_main
CMND=$(v2_cmnd nu)
MAXPAR=${2:-${1:-10}}
[ "${1:-}" = "--list" ] && MAXPAR=0

jobs_list () {
    for d in "$V2_POWHEG_DIR"/v2nu_[pn]_job* "$V2_POWHEG_DIR"/v2numc_[pn]_job*; do
        [ -d "$d" ] || continue
        case "$(basename "$d")" in *_[0-9]|*_[0-9][0-9]) ;; *) continue ;; esac
        echo "$d"
    done
}
lhe_of () { grep -o 'Beams:LHEF *| *[^ |]*' "$1/shower.log" | awk '{print $NF}'; }
md5f () { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }

one () {   # $1 = job dir
    local d=$1 lhe rd t g0 u0 base
    base=$(basename "$d")
    lhe=$(lhe_of "$d")
    [ -n "$lhe" ] || { echo "  $base: no LHE named in shower.log" >&2; return 1; }
    rd=$(dirname "$lhe")
    case "$base" in *_p_*) t=p ;; *_n_*) t=n ;; *) echo "  $base: target?" >&2; return 1 ;; esac
    if [ -f "$rd/UBEXCESS_OK" ] && v2_shower_done "$d" && [ "$d/V2_OK" -nt "$rd/UBEXCESS_OK" ]; then
        echo "  skip $base (regenerated and re-showered)"; return 0
    fi
    if [ ! -f "$rd/UBEXCESS_OK" ]; then
        [ -s "$rd/pwggrid.dat" ] && [ -s "$rd/pwgubound.dat" ] || {
            echo "  $base: $rd lacks its grids -- refusing" >&2; return 1; }
        grep -q '^use-old-grid *1' "$rd/powheg.input" && grep -q '^use-old-ubound *1' "$rd/powheg.input" || {
            echo "  $base: $rd/powheg.input does not reuse grid and bound -- refusing" >&2; return 1; }
        g0=$(md5f "$rd/pwggrid.dat"); u0=$(md5f "$rd/pwgubound.dat")
        grep -q '^ubexcess_correct' "$rd/powheg.input" || \
            printf 'ubexcess_correct 1  ! regen_ubexcess.sh, 2026-09-19\n' >> "$rd/powheg.input"
        rm -f "$rd/pwgevents.lhe" "$rd/pwgevents.lhe.nevents" "$rd/pwgboundviolations.dat"
        echo "  gen   $base  ($(basename "$rd"), iseed $(awk '$1=="iseed" {print $2}' "$rd/powheg.input"))  $(date '+%T')"
        ( cd "$rd" && bench_run "$PWHG" > run_ubexcess.log 2>&1 ) || true
        v2_lhe_complete "$rd/pwgevents.lhe" || { echo "  $base: no complete LHE (run_ubexcess.log)" >&2; return 1; }
        [ "$(md5f "$rd/pwggrid.dat")" = "$g0" ] && [ "$(md5f "$rd/pwgubound.dat")" = "$u0" ] || {
            echo "  $base: grids CHANGED -- pwhg_main re-integrated; refusing" >&2; return 1; }
        "$BENCH_REPO/powheg/strip_lhe_nan.py" "$rd/pwgevents.lhe" > /dev/null
        date '+%F %T' > "$rd/UBEXCESS_OK"
    fi
    rm -rf "$d"
    echo "  show  $base  $(date '+%T')"
    v2_shower "$lhe" "$d" "$CMND" "$t" 14 || { echo "  $base: shower failed" >&2; return 1; }
    v2_shower_done "$d" || { echo "  $base: shower incomplete" >&2; return 1; }
    echo "  done  $base  $(date '+%T')"
}

if [ "$MAXPAR" = 0 ]; then jobs_list; exit 0; fi
if [ "${1:-}" = "--one" ]; then one "$2"; exit $?; fi
N=$(jobs_list | wc -l | tr -d ' ')
echo "regen_ubexcess: $N jobs, $MAXPAR at a time, start $(date '+%F %T')"
v2_disk_check 20 "regen_ubexcess (in place)"
jobs_list | xargs -n 1 -P "$MAXPAR" "$0" --one
echo "regen_ubexcess: end $(date '+%F %T')"
