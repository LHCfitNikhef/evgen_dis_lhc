#!/usr/bin/env bash
# Shower the POWHEG-RES LHE files in parallel -> <jobdir>_N/events.hepmc
#
# Usage: ./run_powheg_shower.sh [--current mu|nu] [--qed ARM] \
#                               [beam_energy_gev] [job list]
#
#   ./run_powheg_shower.sh                    # 1 TeV muon NC -> job_N
#   ./run_powheg_shower.sh 400                # 400 GeV -> job_400GeV_N
#   ./run_powheg_shower.sh --current nu 1000  # the CC cross-variant -> resnu_job_N
#   FORCE=1 ./run_powheg_shower.sh 400 "5 6 7 9"   # redo those four only
#   ./run_powheg_shower.sh --qed fsr 1000     # the QED study -> qedfsr_job_N
#   ./run_powheg_shower.sh --shower vincia 1000   # -> shwvincia_job_N
#
# --qed RE-SHOWERS THE SAME LHE FILES WITH QED RADIATION ON.  The arm names an
# overlay card in powheg/qed/ which is APPENDED to the base Pythia steering,
# so the matrix element, the matching, the PDF and the scale are all held
# fixed and only the radiation changes.  Output goes to qed<ARM>_<jobbase>_N
# so it never collides with the published sample.  See powheg/qed/README.md;
# `--qed off` must reproduce the published sample exactly, which is what makes
# the study self-validating.
#
# THE JOB LIST EXISTS BECAUSE RE-SHOWERING IS A REPAIR OPERATION.  Three seeds
# across the three muon energies were truncated by the NaN trap below, and four
# more at 400 GeV had to be regenerated outright (their LHEs were cut off
# mid-write, with no closing tag).  Re-running the whole energy to fix four
# seeds means re-showering 600000 events to replace 240000, so the list is
# spelled the same way the GENIE drivers spell theirs.
#
# THIS IS THE POWHEG-RES SHOWER STEP, for either current.  POWHEG-V2 produces
# a single pwgevents.lhe rather than one per seed, so run_powheg_v2.sh showers
# its own output inline and does not come here -- see the note in that script
# about why the two drivers have different shapes.
#
# The Pythia card is used at EVERY energy despite its name: it sets no beams at
# all.  Beams:frameType and the LHEF input are set by the driver
# (main_powheg.cc) from the LHE file itself, so the beam energy travels with
# the events rather than the card.  Which card it is DOES depend on the
# current, and powheg_variants.sh is what picks it.
#
# PYTHIA8DATA must point at the 8.311 xmldoc, not conda's 8.312 -- config.sh
# assigns it unconditionally for exactly that reason.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo
. "$HERE/powheg_variants.sh"       # powheg_variant / powheg_parse_current

# --qed is stripped BEFORE powheg_parse_current, which knows only about the
# current and the positional arguments.  bash 3.2: ARR+=("x") is fine (3.1+),
# associative arrays and ${var:t} are not.
QED_ARM=""
SHW_ARM=""
QED_ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --qed)    QED_ARM=$2; shift 2 ;;
        --shower) SHW_ARM=$2; shift 2 ;;
        *)        QED_ARGS+=("$1"); shift ;;
    esac
done
# ONE variant at a time.  Both mechanisms compose an overlay onto the base
# card and both rename the output; allowing them together would produce a
# directory whose name says one thing and whose card says two, and the
# resulting sample would answer neither question.
if [ -n "$QED_ARM" ] && [ -n "$SHW_ARM" ]; then
    echo "--qed and --shower cannot be combined: each is a controlled" \
         "variation against the SAME published baseline, and a sample" \
         "carrying both varies two things at once." >&2
    exit 1
fi
[ ${#QED_ARGS[@]} -gt 0 ] && set -- "${QED_ARGS[@]}" || set --

powheg_parse_current mu "$@"       # default: POWHEG-RES's native current
set -- "${PV_ARGS[@]}"
EBEAM=${1:-1000}
JOBS=${2:-}

powheg_variant RES "$PV_CURRENT" "$EBEAM"
JOBBASE=$PV_JOBBASE

# The QED arm renames the output and supplies an overlay card.  Both checks
# below fail LOUDLY: a missing overlay would otherwise shower the baseline
# under a QED name, and a neutrino `fsrisr` would silently duplicate `fsr`.
QED_OVERLAY=""
if [ -n "$QED_ARM" ]; then
    QED_OVERLAY=$HERE/qed/$QED_ARM.cmnd
    if [ ! -f "$QED_OVERLAY" ]; then
        echo "no such QED arm '$QED_ARM' -- powheg/qed/ has:" >&2
        ls "$HERE"/qed/*.cmnd | sed 's|.*/||; s|\.cmnd$|  |' >&2
        exit 1
    fi
    if [ "$PV_CURRENT" = nu ] && [ "$QED_ARM" = fsrisr ]; then
        echo "arm 'fsrisr' does not exist in the charged current: ISR off" \
             "the lepton line needs a CHARGED incoming lepton, and the" \
             "beam here is a neutrino.  Use --qed fsr or --qed full." >&2
        exit 1
    fi
    JOBBASE=qed${QED_ARM}_$PV_JOBBASE
fi
if [ -n "$SHW_ARM" ]; then
    # ONE copy of the shower overlays, under pythia8/, since they are pure
    # Pythia settings and both drivers use them.  NOTE vincia and dire do NOT
    # work on this path -- the POWHEG veto is a SimpleShower feature -- which
    # the delivery guard at the end of this script will catch; see
    # pythia8/shower/README.md.
    QED_OVERLAY=$BENCH_REPO/pythia8/shower/$SHW_ARM.cmnd
    if [ ! -f "$QED_OVERLAY" ]; then
        echo "no such shower arm '$SHW_ARM' -- pythia8/shower/ has:" >&2
        ls "$BENCH_REPO"/pythia8/shower/*.cmnd | sed 's|.*/||; s|\.cmnd$|  |' >&2
        exit 1
    fi
    JOBBASE=shw${SHW_ARM}_$PV_JOBBASE
fi
LHEDIR=${LHEDIR:-$PV_RUNDIR}

nlhe=$(ls "$LHEDIR"/pwgevents-*.lhe 2>/dev/null | wc -l | tr -d ' ')
[ "$nlhe" -gt 0 ] || {
    echo "no LHE files in $LHEDIR -- run powheg/run_powheg_res.sh" \
         "--current $PV_CURRENT $EBEAM first" >&2
    exit 1; }

# NOTE "${JOBLIST[@]}" below, not $JOBLIST: bash 3.2 expands the latter to the
# FIRST element only, which would silently shower one job instead of N.  That
# exact mistake once ran 1 GENIE job instead of 8 (CONVENTIONS.md rule 1).
if [ -n "$JOBS" ]; then JOBLIST=($JOBS); else JOBLIST=($(seq 1 "$nlhe")); fi

for i in "${JOBLIST[@]}"; do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc -- FORCE=1 to replace" >&2
    exit 1
  fi
done

# GUARD THE READER BEFORE IT READS.  POWHEG-RES writes a NaN four-momentum
# into a handful of CC events, and Pythia's LHEF reader answers an unparseable
# token with "reached end of Les Houches Events File" -- so ONE bad event
# silently truncates the rest of that job and nothing anywhere returns
# non-zero.  On the first nu CC production that cost 14.5% of the sample.
# One pass over a file that is about to be read anyway; a clean file is not
# rewritten.  See powheg/strip_lhe_nan.py.
"$BENCH_REPO/powheg/strip_lhe_nan.py" "$LHEDIR"/pwgevents-*.lhe

echo "$PV_LABEL: showering ${#JOBLIST[@]} of $nlhe LHE files from $LHEDIR" \
     "-> ${JOBBASE}_N"
for i in "${JOBLIST[@]}"; do
    tag=$(printf "%04d" "$i")
    d=$HERE/${JOBBASE}_$i
    mkdir -p "$d"
    # Compose base + overlay per job.  Appending means the overlay's "on"
    # overrides the base card's "off" (Pythia applies settings in order, last
    # wins), so the base card remains the one source of every other setting.
    card=$PV_CMND
    if [ -n "$QED_OVERLAY" ]; then
        card=$d/shower.cmnd
        cat "$PV_CMND" "$QED_OVERLAY" > "$card"
    fi
    ( bench_run ./main_powheg "$card" "$LHEDIR/pwgevents-$tag.lhe" \
        "$d/events" > "$d/shower.log" 2>&1 ) &
done
wait
echo "ALL SHOWER JOBS DONE"
# DELIVERY GUARD.  A shower job can fail on EVERY event and still leave this
# script exiting 0: main_powheg writes an empty events.hepmc, the loop
# completes, and the summary below prints "0 events" to a log nobody re-reads.
# That is how the Vincia arm was first run -- zero events, exit 0 -- and Dire
# delivered 225 of ~60000 for the same reason, the POWHEG veto hooks being a
# SimpleShower feature that the other two shower models do not expose.
# So the count is CHECKED, not merely printed.  The threshold is deliberately
# loose: this is a tripwire for a broken configuration, not a statistics test.
short=0
for i in "${JOBLIST[@]}"; do
    d=$HERE/${JOBBASE}_$i
    # NOT `grep -c ... || echo 0`: grep prints "0" AND returns 1 when there
    # are no matches, so the fallback fires too and the variable becomes the
    # two-line string "0\n0", which `[` then rejects as not an integer.
    n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null); n=${n:-0}
    echo "  ${JOBBASE}_$i: $n events, $(cat "$d/events_xsec.json" 2>/dev/null)"
    lhe=$(grep -c '<event>' "$LHEDIR/pwgevents-$(printf "%04d" "$i").lhe" 2>/dev/null); lhe=${lhe:-0}
    if [ "$lhe" -gt 0 ] && [ "$n" -lt $((lhe / 2)) ]; then
        echo "    *** delivered $n of $lhe LHE events (under half) ***" >&2
        short=$((short + 1))
    fi
done
if [ "$short" -gt 0 ]; then
    echo "" >&2
    echo "$short job(s) delivered under half the events offered.  Something is" >&2
    echo "wrong with the configuration, not with statistics -- check the" >&2
    echo "shower.log for aborts.  The sample is NOT usable." >&2
    exit 1
fi
