#!/usr/bin/env bash
# Parallel Pythia8 production for the muon-DIS benchmark: N seeded jobs.
# Usage: ./run_pythia.sh <events_per_job> <njobs> [beam_energy_gev]
#
# The energy defaults to the 1 TeV anchor, whose job dirs and card keep their
# original untagged names so nothing published moves.  Other energies use
# dis_mu<tag>.cmnd and job_<tag>_N -- both names come from analysis/beams.py,
# which is the one place that knows the rule.
#
#   ./run_pythia.sh 250000 8          # 1 TeV, the published sample
#   ./run_pythia.sh 250000 8 400      # 400 GeV -> job_400GeV_1..8
#   ./run_pythia.sh --shower vincia 250000 1   # -> shwvincia_job_1
#
# --shower RE-RUNS THE SAME CARD WITH A DIFFERENT PARTON SHOWER.  The arm names
# an overlay in pythia8/shower/ which is APPENDED to the base card, so the
# matrix element, the PDF, the scales and the cuts are all held fixed and only
# the shower changes.  Output goes to shw<ARM>_<jobbase>_N.
#
# THE COMPARISON IS DONE AT LO FOR A REASON.  On the POWHEG samples the veto
# hooks are a SimpleShower feature that Vincia and Dire do not expose, and both
# arms deliver essentially nothing (logs/shower_probe_powheg.log).  At LO there
# is no matching veto at all, so the shower is the only thing that differs.
set -e
SHW_ARM=""
PY_ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --shower) SHW_ARM=$2; shift 2 ;;
        *)        PY_ARGS+=("$1"); shift ;;
    esac
done
[ ${#PY_ARGS[@]} -gt 0 ] && set -- "${PY_ARGS[@]}" || set --
NEV=${1:-50000}
NJOBS=${2:-8}
EBEAM=${3:-1000}
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag mu "$EBEAM")
JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name job "$EBEAM")
CARD=$HERE/dis_mu$TAG.cmnd
[ -f "$CARD" ] || {
    echo "no card $CARD -- run tools/make_energy_cards.py first" >&2; exit 1; }

# The shower arm renames the output and supplies an overlay.  A missing overlay
# would otherwise run the BASELINE shower under a variant name, which is the
# quiet kind of wrong this repository keeps finding, so it fails loudly.
SHW_OVERLAY=""
if [ -n "$SHW_ARM" ]; then
    SHW_OVERLAY=$HERE/shower/$SHW_ARM.cmnd
    if [ ! -f "$SHW_OVERLAY" ]; then
        echo "no such shower arm '$SHW_ARM' -- pythia8/shower/ has:" >&2
        ls "$HERE"/shower/*.cmnd | sed 's|.*/||; s|\.cmnd$|  |' >&2
        exit 1
    fi
    JOBBASE=shw${SHW_ARM}_$JOBBASE
fi

# Seeds must differ between jobs; the offset keeps different energies from
# reusing the same stream.  The anchor keeps 4321 so it reproduces exactly.
if [ "$EBEAM" = "1000" ]; then SEED0=4321; else SEED0=$((4321 + EBEAM)); fi

# A job dir that already holds events is a PRODUCTION SAMPLE -- job_1 alone is
# 1.3 GB.  Refuse rather than silently overwrite days of running; FORCE=1 to
# override deliberately.
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $d/events.hepmc ($(du -h "$d/events.hepmc" | cut -f1))." >&2
    echo "  Re-run with FORCE=1 to replace it." >&2
    exit 1
  fi
done

echo "Pythia8: $NJOBS x $NEV events, mu- $EBEAM GeV [$TAG] -> ${JOBBASE}_N"
for i in $(seq 1 "$NJOBS"); do
  d=$HERE/${JOBBASE}_$i
  mkdir -p "$d"
  # base + overlay, composed per job.  Appending means the overlay wins:
  # Pythia applies settings in order.  The SEED IS UNCHANGED by the arm, so
  # the `simple` arm must reproduce the published job exactly.
  card=$CARD
  if [ -n "$SHW_OVERLAY" ]; then
      card=$d/run.cmnd
      cat "$CARD" "$SHW_OVERLAY" > "$card"
  fi
  (cd "$d" && bench_run "$HERE/main_dis" "$card" events \
      "$NEV" $((SEED0 + i)) > pythia.log 2>&1) &
done
wait
echo "all $NJOBS Pythia jobs done"
grep -h "sigma_gen_mb" "$HERE"/${JOBBASE}_*/events_xsec.json

# Analyse and delete, if asked.  The sample is an intermediate: nothing reads
# an event file once analyze.py has run.  The pruner refuses to delete unless
# every result was written, is non-empty and is newer than the sample, and
# refuses the 1 TeV anchor outright while the FASER selection still has to
# re-parse it.  See config.sh for the knobs.
if [ "${BENCH_PRUNE:-0}" = "1" ]; then
    echo "BENCH_PRUNE=1: analysing, then deleting the event files"
    keep=""
    [ "${BENCH_PRUNE_KEEP_ONE:-0}" = "1" ] && keep="--keep-one"
    "$BENCH_REPO/tools/analyse_and_prune.py" --gen pythia --energy "$EBEAM" $keep
fi
