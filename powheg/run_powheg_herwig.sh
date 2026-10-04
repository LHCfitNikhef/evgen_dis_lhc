#!/usr/bin/env bash
# Shower a POWHEG Les Houches sample with HERWIG 7 instead of Pythia 8.
#
# Usage: ./run_powheg_herwig.sh [--current mu|nu] [beam_energy_gev] [job list]
#
#   ./run_powheg_herwig.sh                  # 1 TeV muon NC  -> hw_job_N
#   ./run_powheg_herwig.sh --current nu     # 1 TeV nu CC     -> hw_nu_job_N
#   ./run_powheg_herwig.sh 400 "1 2"        # two seeds at 400 GeV
#
# WHAT THIS IS FOR (user, 2026-09-08): "what we want is to compare one POWHEG
# LHE with Pythia8 and with Herwig 7, varying shower + hadronisation at fixed
# matrix element."  The events are the SAME events the published `powheg`
# (muon) and `powheg_nu` (neutrino) samples were showered from -- same file,
# same weights -- so the matrix element, the matching, the parton distribution
# and the scale are held fixed by construction and the difference is
#
#     Pythia 8.311   p_T-ordered dipole shower + Lund string
#     Herwig 7.3.0   angular-ordered shower    + cluster model
#
# IT IS ALSO THE SHOWER STUDY AT NLO.  Switching PartonShowers:model inside
# Pythia does not work on a POWHEG sample -- the veto hooks are a SimpleShower
# feature and Vincia and Dire deliver 0 and 225 events of 59775
# (pythia8/shower/README.md, logs/shower_probe_powheg.log).  Herwig applies the
# POWHEG veto through its own MaxPtIsMuF/RestrictPhasespace machinery, so this
# path carries an NLO-matched shower comparison that the model switch cannot.
#
# >>> THE MUON LHE MUST BE BOOSTED FIRST, AND THE NEUTRINO ONE MUST NOT. <<<
# POWHEG-RES runs the muon side with `fixed_target 1`: the init line declares
# beam B as a proton of energy 0.938 GeV, i.e. AT REST WITH ZERO MOMENTUM.
# ThePEG builds its parton bins from momentum fractions along the beam, so a
# beam with no momentum is a division by zero -- the run dies at `saverun`
# with "Tried to write a NaN or Inf double to a persistent stream", AFTER
# writing a .run file that then fails to read back.  POWHEG-V2 writes the
# neutrino side already in the c.m. frame and needs nothing.
# tools/lhe_boost_cm.py decides which case it has, boosts or symlinks, and
# says so; the boost is exact and the analysis is covariant, so this changes
# no number.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/powheg_variants.sh"

powheg_parse_current mu "$@"
set -- "${PV_ARGS[@]}"
EBEAM=${1:-1000}
JOBS=${2:-}

# WHICH CODE PER CURRENT, and it is each current's OWN published NLO row:
# POWHEG-RES on the muon neutral current, POWHEG-V2 on the neutrino charged
# current.  Using one code for both would make the Herwig arm comparable to a
# row that is not on any figure.
case "$PV_CURRENT" in
    mu) CODE=RES ;;
    nu) CODE=V2  ;;
esac
powheg_variant "$CODE" "$PV_CURRENT" "$EBEAM"
CARD=$BENCH_REPO/herwig7/LHE-$PV_CURRENT.in
[ -f "$CARD" ] || { echo "no Herwig card $CARD" >&2; exit 1; }

# hw_job_N (muon, one job per POWHEG-RES seed file) and hw_nu_job_N (neutrino,
# one job for POWHEG-V2's single pwgevents.lhe).  The prefix keeps them clear
# of `job_*` and of each other: `hw_job_*` does not match `hw_nu_job_*`, which
# is the glob trap CONVENTIONS.md rule 1 records.
case "$PV_CURRENT" in
    mu) JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name hw_job "$EBEAM") ;;
    nu) JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name hw_nu_job "$EBEAM") ;;
esac

# THE JOB INDEX MUST MATCH THE PYTHIA ARM'S, or the comparison is between
# two different subsets of the same integral and every difference carries an
# extra statistical fluctuation nobody accounted for.
#
#   POWHEG-RES writes one pwgevents-000N.lhe per seed in ONE directory, and
#   run_powheg_shower.sh showers seed N into job_N.  Index = seed.
#
#   POWHEG-V2 writes ONE pwgevents.lhe per run directory, and
#   add_v2_batches.sh adds statistically independent batches as
#   <rundir>-b<j>, showering batch j into nu_job_j (job 1 being the base
#   directory itself).  Index = batch, and the ORDER below reproduces it.
LHES=$(ls "$PV_RUNDIR"/pwgevents-*.lhe 2>/dev/null || true)
if [ -z "$LHES" ]; then
    LHES=$(ls "$PV_RUNDIR"/pwgevents.lhe 2>/dev/null || true)
    j=2
    while [ -s "$PV_RUNDIR-b$j/pwgevents.lhe" ]; do
        LHES="$LHES $PV_RUNDIR-b$j/pwgevents.lhe"
        j=$((j + 1))
    done
fi
[ -n "$LHES" ] || { echo "no LHE files in $PV_RUNDIR" >&2; exit 1; }
set -- $LHES
NLHE=$#
if [ -n "$JOBS" ]; then JOBLIST=($JOBS); else JOBLIST=($(seq 1 "$NLHE")); fi

for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc -- FORCE=1 to replace" >&2
    exit 1
  fi
done

# The same NaN guard the Pythia path runs: POWHEG-RES writes a NaN
# four-momentum into a handful of CC events, and a reader that answers an
# unparseable token with "end of file" truncates the rest of the job silently.
"$BENCH_REPO/powheg/strip_lhe_nan.py" $LHES

echo "$PV_LABEL + Herwig 7: ${#JOBLIST[@]} of $NLHE LHE file(s) -> ${JOBBASE}_N"
for i in "${JOBLIST[@]}"; do
  (
    eval "src=\${$i}"
    d=$HERE/${JOBBASE}_$i
    mkdir -p "$d"
    cd "$d"
    python3 "$BENCH_REPO/tools/lhe_boost_cm.py" "$src" "$d/events_cm.lhe" \
        > boost.log 2>&1
    nev=$(grep -c '<event>' "$d/events_cm.lhe")
    # The run name carries the job so two jobs cannot share a repository file.
    sed -e "s|@LHEFILE@|$d/events_cm.lhe|" -e "s|@RUNNAME@|LHE-${JOBBASE}_$i|" \
        "$CARD" > "$d/shower.in"
    "$HERWIG_BIN/Herwig" read "$d/shower.in" > read.log 2>&1
    # >>> `|| true` IS DELIBERATE, AND THE DELIVERY GUARD BELOW IS WHY. <<<
    # Herwig ENDS this run by throwing: asking for the file's own event count
    # exhausts the reader, because every event Herwig discards in the shower
    # consumes one more Les Houches event.  It throws AFTER writing everything
    # it made, so the sample is complete and the exit code is not.  Asking for
    # fewer events instead would throw away good statistics to avoid a
    # message.  What must not be swallowed is a REAL failure, so the guard
    # below counts what arrived rather than trusting the exit code -- and it
    # is set at 90%, not the other drivers' 50%, because the only expected
    # loss here is Herwig's few-per-cent shower discard.
    bench_run "$HERWIG_BIN/Herwig" run "LHE-${JOBBASE}_$i.run" \
        -N "$nev" -s "$i" > shower.log 2>&1 || true
    n_out=$(grep -c '^E ' events.hepmc 2>/dev/null); n_out=${n_out:-0}
    printf '{"n_lhe": %s, "n_delivered": %s}\n' "$nev" "$n_out" > delivery.json
    # NOTHING BUT THE EVENTS AND THE BOOKKEEPING OUTLIVES THE SHOWER
    # (CONVENTIONS.md rule 1).  The boosted copy is 44-130 MB per job and
    # tools/lhe_boost_cm.py rebuilds it in seconds; the ThePEG repository
    # dump is 8 MB and is rebuilt by `Herwig read`; the per-event log is
    # 3.8 MB of progress lines.  The .out file stays: it carries Herwig's own
    # cross-section summary.
    rm -f "$d/events_cm.lhe" "$d/LHE-${JOBBASE}_$i.run" "$d"/*-S*.log
  ) &
done
wait
echo "ALL HERWIG SHOWER JOBS DONE"

# DELIVERY GUARD, the same one run_powheg_shower.sh carries and for the same
# reason: a shower that fails on every event still leaves this script exiting
# 0, with an empty events.hepmc and a count printed to a log nobody re-reads.
short=0
for i in "${JOBLIST[@]}"; do
    eval "src=\${$i}"
    d=$HERE/${JOBBASE}_$i
    n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null); n=${n:-0}
    lhe=$(grep -c '<event>' "$src" 2>/dev/null); lhe=${lhe:-0}
    echo "  ${JOBBASE}_$i: $n events of $lhe offered" \
         "($(python3 -c "print('%.2f%%' % (100.0*(1-$n/max($lhe,1))))") discarded)"
    if [ "$lhe" -gt 0 ] && [ "$n" -lt $((lhe * 9 / 10)) ]; then
        echo "    *** delivered $n of $lhe LHE events (under 90%) ***" >&2
        short=$((short + 1))
    fi
done
if [ "$short" -gt 0 ]; then
    echo "" >&2
    echo "$short job(s) delivered under 90% of the events offered.  Herwig" >&2
    echo "discards a few per cent in the shower; ten is a broken" >&2
    echo "broken configuration, not statistics -- read the shower.log." >&2
    exit 1
fi
