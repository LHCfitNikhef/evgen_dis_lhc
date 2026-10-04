#!/usr/bin/env bash
# THE HERWIG ARM OF THE SHOWER STUDY IN  (paper plot 9): the  1 TeV POWHEG
# Les Houches files, showered and hadronised by Herwig 7 instead of Pythia 8.
#
# Usage:  powheg/production/herwig_arm.sh <mu|nu> <p|n> [njobs] [max_parallel=4] [nev]
#
#   njobs  baseline jobs to re-shower, lowest job numbers first -- the SAME
#          subset the QED arms use (default mu 10, nu 2), so the Pythia and
#          Herwig arms compare the same events
#   nev    events per job; default the whole LHE file.  A small number is a
#          test run: it goes to hwtest_* and is never picked up as a sample.
#
# Output: powheg/hw_v2<cur>_<t>_job_<N>, one per baseline job <cur>_<t>_job_<N>,
# reading the LHE that job's shower.log names.  The card is herwig7/LHE-<cur>.in,
# the earlier study's card (the matching, the PDF, 0.5 GeV intrinsic pT, no QED,
# ctau > 10 mm), with the target handled as the  Herwig production does
# (herwig7/production/make_cards.py):
#
#   NEUTRON: the LHE declares beam B as 2112, so ThePEG builds the beam from
#   /Herwig/Particles/n0, which is given the PROTON set.  ThePEG's LHAPDF
#   interface isospin-swaps the set itself for id 2112 (LHAPDF6.cc::xfx);
#   giving it the neutron set would swap twice and shower a proton.
#   herwig7/production/inspect_sample.py checks each job's remnant diquarks and GenPdfInfo.
#
# The muon LHE (POWHEG-RES, fixed target) is boosted to the c.m. frame by
# tools/lhe_boost_cm.py; the neutrino one (POWHEG-V2) is already collinear.
# POWHEG-RES writes beam B at 0.938 GeV for the neutron as well: its
# kinematics use one nucleon mass, and the boost follows the file.
#
# Delivery guard as run_powheg_herwig.sh: Herwig throws when it exhausts the
# reader, after writing everything, so the exit code is ignored and >= 90% of
# the offered events must arrive.  Safe to kill and re-invoke: HW_OK skips.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
cd "$V2_POWHEG_DIR"

baseline_jobs () {
    local d rest
    for d in "$V2_POWHEG_DIR/v2$1_$2_job"_*; do
        rest=${d#"$V2_POWHEG_DIR/v2$1_$2_job"_}
        case "$rest" in ''|*[!0-9]*) continue ;; esac
        v2_shower_done "$d" && echo "$rest"
    done | sort -n
}
lhe_of () { grep -o 'Beams:LHEF *| *[^ |]*' "$1/shower.log" | awk '{print $NF}'; }

one () {    # $1 cur, $2 target, $3 job number, $4 nev (or "all")
    local cur=$1 t=$2 n=$3 nev=$4 base src d stem run nlhe nout card
    base=$V2_POWHEG_DIR/v2${cur}_${t}_job_$n
    stem=hw_v2${cur}_${t}_job_$n
    [ "$nev" = all ] || stem=hwtest_v2${cur}_${t}_job_$n
    d=$V2_POWHEG_DIR/$stem
    [ -f "$d/HW_OK" ] && { echo "  skip $stem (complete)"; return 0; }
    src=$(lhe_of "$base")
    [ -s "$src" ] || { echo "  no LHE for $(basename "$base")" >&2; return 1; }
    rm -rf "$d"; mkdir -p "$d"; cd "$d"
    python3 "$BENCH_REPO/tools/lhe_boost_cm.py" "$src" "$d/events_cm.lhe" > boost.log 2>&1
    nlhe=$(grep -c '<event>' "$d/events_cm.lhe")
    [ "$nev" = all ] && nev=$nlhe
    run=LHE-$stem
    card=$d/shower.in
    sed -e "s|@LHEFILE@|$d/events_cm.lhe|" -e "s|@RUNNAME@|$run|" \
        "$BENCH_REPO/herwig7/LHE-$cur.in" > "$card"
    # ERROR CAP.  LHE-<cur>.in allows 10000 errors, written for 50k-event
    # jobs.  Every "Can't put the remnant on-shell" discard counts as one, and
    # a 260k-event POWHEG-V2 batch on a proton discards ~5%: the first run
    # stopped at 183k of 260k events on the cap (2026-09-14).  The cap is a
    # tripwire for a broken configuration, which the >= 90% delivery guard
    # below already is, so it is raised to the job size.
    grep -q '^set EventGenerator:MaxErrors 10000$' "$card" || {
        echo "  LHE-$cur.in no longer sets MaxErrors 10000 -- refusing" >&2; return 1; }
    # -i.bak, not -i '': the one in-place form BSD and GNU sed both accept
    sed -i.bak "s|^set EventGenerator:MaxErrors 10000\$|set EventGenerator:MaxErrors $nev|" "$card" && rm -f "$card.bak"
    if [ "$t" = n ]; then
        grep -q '^set /Herwig/Particles/p+:PDF NNPDF40$' "$card" || {
            echo "  LHE-$cur.in no longer sets p+:PDF -- refusing" >&2; return 1; }
        # the neutron beam takes the PROTON set; ThePEG swaps u<->d for 2112
        sed -i.bak 's|^set /Herwig/Particles/p+:PDF NNPDF40$|&\
#  NEUTRON (powheg/production/herwig_arm.sh): n0 gets the PROTON set, ThePEG\
# isospin-swaps it for id 2112 itself (LHAPDF6.cc::xfx)\
set /Herwig/Particles/n0:PDF NNPDF40|' "$card" && rm -f "$card.bak"
    fi
    "$HERWIG_BIN/Herwig" read "$card" > read.log 2>&1 || {
        echo "  Herwig read failed for $stem (read.log)" >&2; return 1; }
    bench_run "$HERWIG_BIN/Herwig" run "$run.run" -N "$nev" -s "$((7000 + n))" \
        > shower.log 2>&1 || true
    nout=$(grep -c '^E ' events.hepmc 2>/dev/null); nout=${nout:-0}
    printf '{"n_lhe": %s, "n_requested": %s, "n_delivered": %s, "baseline": "%s", "lhe": "%s"}\n' \
        "$nlhe" "$nev" "$nout" "$(basename "$base")" "$src" > delivery.json
    rm -f "$d/events_cm.lhe" "$d/$run.run" "$d"/*-S*.log
    if [ "$nout" -lt $(( nev * 9 / 10 )) ]; then
        echo "  SHORT: $stem delivered $nout of $nev (under 90%) -- see shower.log" >&2
        return 1
    fi
    echo "$nout of $nev events, target $t, $(date)" > HW_OK
    echo "  ok $stem: $nout/$nev"
}

if [ "$1" = "--one" ]; then shift; one "$@"; exit $?; fi

cur=${1:?usage: herwig_arm.sh <mu|nu> <p|n> [njobs] [max_parallel] [nev]}
t=${2:?target p|n}
case "$cur" in mu) def=10 ;; nu) def=2 ;; *) echo "current mu|nu" >&2; exit 1 ;; esac
case "$t" in p|n) ;; *) echo "target p|n" >&2; exit 1 ;; esac
njobs=${3:-$def}
maxpar=${4:-4}
nev=${5:-all}
list=$(baseline_jobs "$cur" "$t" | head -"$njobs")
gb=0
for n in $list; do
    k=$(du -sk "$V2_POWHEG_DIR/v2${cur}_${t}_job_$n/events.hepmc" | awk '{print $1}')
    gb=$(( gb + k / 1048576 + 1 ))
done
[ "$nev" = all ] && v2_disk_check "$gb" "Herwig arm $cur $t"
echo "Herwig arm, $cur on $t: baseline jobs $(echo $list | tr '\n' ' ') (nev $nev)"
for n in $list; do echo "$cur $t $n $nev"; done \
    | xargs -n 4 -P "$maxpar" "$HERE/herwig_arm.sh" --one
echo "== Herwig arm $cur $t finished =="
