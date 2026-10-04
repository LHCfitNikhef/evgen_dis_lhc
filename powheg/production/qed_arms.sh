#!/usr/bin/env bash
# THE QED STUDY IN  (paper plots 10 and 11): re-shower  1 TeV POWHEG
# samples with QED radiation switched on in stages.
#
# Usage:  powheg/production/qed_arms.sh <mu|nu> <p|n> <arm> [njobs] [max_parallel=4]
#         powheg/production/qed_arms.sh closure <mu|nu> <p|n>
#
#   arm    fsr | fsrisr (mu only) | full      -- the overlays in powheg/qed/
#   njobs  how many baseline jobs to re-shower, lowest job numbers first
#          (default: mu 10 of ~21 seeds = ~500k events, nu 2 of 4 batches =
#          ~520k events, per nucleon)
#
# WHAT IS HELD FIXED.  Each arm re-showers the SAME LHE file its baseline job
# showered (read from that job's shower.log, Beams:LHEF), with the SAME base
# card (v2_cmnd: powheg_mu1TeV.cmnd for POWHEG-RES, powheg_nu_v2.cmnd --
# matchInOut off -- for POWHEG-V2) and the overlay APPENDED, so only the QED
# radiation changes.  Output: powheg/qed<arm>_<baseline job name>, e.g.
# qedfsr_v2mu_p_job_3.  The arm is compared with ITS OWN baseline jobs (the
# same job numbers), never with the whole baseline sample.
#
# THE NEUTRINO ARMS RE-SHOWER POWHEG-V2 (user, 2026-09-14), not POWHEG-RES as
# in the earlier production: the production has no POWHEG-RES neutrino sample, and a ratio
# taken against its own QED-off baseline cancels the generator.
#
# `closure` re-showers ONE baseline job with the empty `off` overlay into
# powheg/qedoff_<job> and compares its events.hepmc with the baseline's.
# main_powheg sets no Random:setSeed, so the two must be byte-identical; if
# they are not, the arms are not event-by-event comparable to the baseline
# and the study's small errors on a difference do not hold.  The closure
# directory is deleted after a successful comparison (a diagnostic).
#
# Safe to kill and re-invoke (CONVENTIONS.md 1c): a job with V2_OK is skipped.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
cd "$V2_POWHEG_DIR"

lep_of () { case "$1" in mu) echo 13 ;; nu) echo 14 ;; esac; }

baseline_jobs () {    # $1 cur, $2 target: 1 TeV complete jobs, numeric order
    local d rest
    for d in "$V2_POWHEG_DIR/v2$1_$2_job"_*; do
        rest=${d#"$V2_POWHEG_DIR/v2$1_$2_job"_}
        case "$rest" in ''|*[!0-9]*) continue ;; esac
        v2_shower_done "$d" && echo "$rest"
    done | sort -n
}

lhe_of () { grep -o 'Beams:LHEF *| *[^ |]*' "$1/shower.log" | awk '{print $NF}'; }

one () {    # $1 cur, $2 target, $3 arm, $4 baseline job number
    local cur=$1 t=$2 arm=$3 n=$4 base lhe out card
    base=$V2_POWHEG_DIR/v2${cur}_${t}_job_$n
    out=$V2_POWHEG_DIR/qed${arm}_v2${cur}_${t}_job_$n
    v2_shower_done "$out" && { echo "  skip $(basename "$out") (complete)"; return 0; }
    lhe=$(lhe_of "$base")
    [ -s "$lhe" ] || { echo "  no LHE for $(basename "$base"): '$lhe'" >&2; return 1; }
    rm -rf "$out"; mkdir -p "$out"
    card=$out/shower.cmnd
    cat "$(v2_cmnd "$cur")" "$BENCH_REPO/powheg/qed/$arm.cmnd" > "$card"
    v2_shower "$lhe" "$out" "$card" "$t" "$(lep_of "$cur")"
    echo "arm $arm on baseline $(basename "$base"), LHE $lhe" >> "$out/V2_OK"
}

if [ "$1" = "--one" ]; then shift; one "$@"; exit $?; fi

if [ "$1" = closure ]; then
    cur=$2; t=$3
    n=$(baseline_jobs "$cur" "$t" | head -1)
    one "$cur" "$t" off "$n"
    a=$V2_POWHEG_DIR/v2${cur}_${t}_job_$n/events.hepmc
    b=$V2_POWHEG_DIR/qedoff_v2${cur}_${t}_job_$n/events.hepmc
    if cmp -s "$a" "$b"; then
        echo "CLOSURE OK: qedoff re-shower of v2${cur}_${t}_job_$n is byte-identical to the baseline"
        rm -rf "$V2_POWHEG_DIR/qedoff_v2${cur}_${t}_job_$n"
    else
        echo "CLOSURE FAILED: $b differs from $a -- kept for inspection" >&2
        exit 1
    fi
    exit 0
fi

cur=${1:?usage: qed_arms.sh <mu|nu> <p|n> <arm> [njobs] [max_parallel]}
t=${2:?target p|n}
arm=${3:?arm fsr|fsrisr|full}
case "$arm" in fsr|fsrisr|full) ;; *) echo "arm must be fsr, fsrisr or full" >&2; exit 1 ;; esac
if [ "$cur" = nu ] && [ "$arm" = fsrisr ]; then
    echo "arm fsrisr does not exist in the charged current (neutral beam)" >&2; exit 1
fi
case "$cur" in mu) def=10 ;; nu) def=2 ;; *) echo "current mu|nu" >&2; exit 1 ;; esac
njobs=${4:-$def}
maxpar=${5:-4}
list=$(baseline_jobs "$cur" "$t" | head -"$njobs")
# disk: each arm job writes about what its baseline holds
gb=0
for n in $list; do
    k=$(du -sk "$V2_POWHEG_DIR/v2${cur}_${t}_job_$n/events.hepmc" | awk '{print $1}')
    gb=$(( gb + k / 1048576 + 1 ))
done
v2_disk_check "$gb" "QED $arm $cur $t"
echo "QED arm $arm, $cur on $t: baseline jobs $(echo $list | tr '\n' ' ')"
for n in $list; do echo "$cur $t $arm $n"; done \
    | xargs -n 4 -P "$maxpar" "$HERE/qed_arms.sh" --one
echo "== QED arm $arm $cur $t finished =="
