#!/usr/bin/env bash
# THE HERWIG FASER LADDER (user, 2026-10-04): Herwig 7 POWHEG-NLO for the
# comparisons with FASER's data (paper plots 12, 12b) and the antineutrino
# SIDIS yields (paper plot 15).
#
# Usage: nohup /bin/bash tools/herwig_faser_ladder.sh > logs/herwig_faser_ladder.log 2>&1 &
#        tools/herwig_faser_ladder.sh --list
#
# Cards: herwig7/production/cards/FL-<nu|nubar>-PWG[NEG]-<p|n>-<tag>.in
# (make_cards.py, region "faser": Q2 > 4 at generation, no y cut, MinW2 0.3),
# at the eleven energies of make_cards.FASER_ENERGIES, nu and nubar, p and a
# genuine neutron.  Per point: 2 x 50k positive-weight events and 1 x 10k
# negative-weight ones (the negative half is ~0.1% of the CC rate).
#
# AFTER EACH JOB, AT ONCE (disk, CONVENTIONS.md rule 1): the per-event FASER
# table is extracted (analysis/faser_emulsion_shapes.py extract ->
# faserdata_events.npz) and events.hepmc is DELETED -- except for the
# positive-weight antineutrino jobs at the SIDIS energies (300 GeV-4 TeV),
# whose hadron spectra analysis/faser_pions.py still has to read; those are
# deleted by that pass (tools/nubar_sidis.py).  The .out/.log files, which
# carry the cross-section, stay.
#
# Re-invoke freely: a complete job (lib.sh v2_job_state, or a job whose
# hepmc was already turned into faserdata_events.npz) is skipped.
# DO NOT EDIT WHILE IT RUNS (CONVENTIONS.md rule 1c).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
. "$HERE/../herwig7/production/lib.sh"
MAXPAR=${MAXPAR:-10}
KB_PER_EVENT=13
PY=${BENCH_PYTHON:-python3}
ES="1000 400 4000 700 2000 300 200 100 6800 60 30"   # launch order
SIDIS_ES=" 300 400 700 1000 2000 4000 "

fl_tag() { python3 "$V2_REPO/analysis/beams.py" tag nu "$1"; }
fl_stem() {   # <cur> <pos|neg> <t> <E>
    local h=PWG; [ "$2" = neg ] && h=PWGNEG
    echo "FL-$1-$h-$3-$(fl_tag "$4")"
}
fl_dir() {    # <cur> <pos|neg> <t> <E> <N>
    local h=pwg; [ "$2" = neg ] && h=pwgneg
    echo "$V2_HDIR/fl$1${h}_$3_job_$(fl_tag "$4")_$5"
}
fl_seed() {   # distinct from every production seed (those are 2xxxxxx)
    local ci=0 ti=0 hi=0 ei=0 e
    [ "$1" = nubar ] && ci=1; [ "$3" = n ] && ti=1; [ "$2" = neg ] && hi=1
    for e in 30 60 100 200 300 400 700 1000 2000 4000 6800; do
        [ "$e" = "$4" ] && break; ei=$((ei + 1)); done
    echo $((3000000 + ci*400000 + ti*200000 + hi*100000 + ei*1000 + $5))
}

tasks() {   # cur half t E N nev
    local E cur t
    for E in $ES; do for cur in nu nubar; do for t in p n; do
        echo "$cur pos $t $E 1 50000"
        echo "$cur pos $t $E 2 50000"
        echo "$cur neg $t $E 1 10000"
    done; done; done
}

if [ "${1:-}" = "--list" ]; then
    tasks | awk '{n+=$6; t++} END {print t, "jobs,", n, "events"}'; exit 0
fi

postprocess() {   # <cur> <half> <t> <E> <dir>
    local cur=$1 half=$2 t=$3 E=$4 d=$5 beam=14
    [ "$cur" = nubar ] && beam=-14
    "$PY" "$V2_REPO/analysis/faser_emulsion_shapes.py" extract \
        --beam="$beam" --target="$t" --energy="$E" --generator=herwig_nlo \
        "$d/events.hepmc" "$d/faserdata_events.npz" > "$d/extract.log" 2>&1 || {
        echo "  EXTRACT FAILED $(basename "$d")" >&2; return 1; }
    if [ "$cur" = nubar ] && [ "$half" = pos ] && [[ "$SIDIS_ES" == *" $E "* ]]; then
        echo "  kept  $(basename "$d")/events.hepmc for the SIDIS pass"
    else
        rm -f "$d/events.hepmc"
    fi
}

LOCK=$V2_BUILD/faser_ladder.lock
mkdir -p "$V2_BUILD"
mkdir "$LOCK" 2>/dev/null || { echo "another ladder driver holds $LOCK" >&2; exit 1; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT
echo "Herwig FASER ladder start $(date)  (MAXPAR $MAXPAR, floor $V2_DISK_FLOOR_GB GB)"
if python3 "$V2_REPO/herwig7/production/make_cards.py" --check | grep -q "would write"; then
    echo "the cards on disk differ from make_cards.py -- run it first" >&2; exit 1
fi
FAILED=$V2_BUILD/faser_ladder.failed.log; : > "$FAILED"
TASKS=$V2_BUILD/faser_ladder.tasks.log; tasks > "$TASKS"
while read -r cur half t E N nev <&3; do
    stem=$(fl_stem "$cur" "$half" "$t" "$E")
    d=$(fl_dir "$cur" "$half" "$t" "$E" "$N")
    seed=$(fl_seed "$cur" "$half" "$t" "$E" "$N")
    if [ -s "$d/faserdata_events.npz" ] && [ -s "$d/$stem-S$seed.out" ]; then
        echo "  skip  $(basename "$d") (extracted)"; continue; fi
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXPAR" ]; do sleep 10; done
    running=$(jobs -pr | wc -l | tr -d ' ')
    need=$(( (running + 1) * nev * KB_PER_EVENT / 1048576 + 1 ))
    avail=$(v2_avail_gb)
    if [ $((avail - need)) -lt "$V2_DISK_FLOOR_GB" ]; then
        echo "DISK: ${avail} GB free, need ${need} GB -- stopping launches" | tee -a "$FAILED"
        break
    fi
    v2_compile "$stem" || { echo "  compile failed: $stem" | tee -a "$FAILED"; continue; }
    ( v2_run_in "$stem" "$d" "$seed" "$nev" \
        && postprocess "$cur" "$half" "$t" "$E" "$d" \
        || echo "  FAILED $cur $half $t $E job $N" >> "$FAILED" ) < /dev/null &
done 3< "$TASKS"
wait
echo "Herwig FASER ladder end $(date)"
[ -s "$FAILED" ] && { echo "PROBLEMS:"; cat "$FAILED"; exit 2; }
echo "all tasks complete"
