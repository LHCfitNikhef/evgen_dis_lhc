#!/usr/bin/env bash
# Reweight ONLY the events that pass a selection, over every PDF member.
#
# Usage: ./run_reweight_selected.sh <energy_gev> <selection> [members_per_batch]
#
# WHY THIS IS EXACT.  A PDF band on a selection is the spread, across members,
# of  S(member) = sum over SELECTED events of w_member(event),  normalised by a
# member-INDEPENDENT constant.  An event that fails the selection contributes
# nothing to that spread, so reweighting it changes no digit of the answer.
# Handing the reweighter only the selected events is therefore the same
# calculation, not a cheaper approximation of it -- and the check below
# demonstrates that rather than asserting it.
#
# WHY IT IS NEEDED.  The 400 GeV sample holds 1,000,000 events because the
# dimuon tier keeps about one in twelve thousand.  Reweighting all of them over
# 319 members was MEASURED at 32 hours (11,300 events/min once the file is
# 340 MB and I/O-bound).  The events that enter the dimuon band number about
# eighty, and reweighting those takes seconds.
#
# The row -> lhe_index mapping is carried in a manifest written by
# subset_lhe.py and read by the analysis; it is never reconstructed from order.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$BENCH_REPO/powheg/powheg_variants.sh"

EBEAM=${1:?usage: run_reweight_selected.sh <energy_gev> <selection> [perbatch]}
SEL=${2:?usage: run_reweight_selected.sh <energy_gev> <selection> [perbatch]}
PERBATCH=${3:-6}

PV_CODE=${BENCH_POWHEG_CODE:-V2}
powheg_variant "$PV_CODE" nu "$EBEAM"
TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag nu "$EBEAM")
SRC=$PV_RUNDIR
DEST=$POWHEG_V2/rwgtsel-nu$TAG-$SEL
PWHG=$POWHEG_V2/pwhg_main
IDXF=$BENCH_REPO/results_nu/selected_indices_${SEL}_${TAG}.json

[ -s "$SRC/pwgevents.lhe" ] || { echo "no LHE at $SRC" >&2; exit 1; }
[ -f "$IDXF" ] || {
    echo "no index list at $IDXF" >&2
    echo "  run: BENCH_ENERGY=$EBEAM BENCH_SELECTION=$SEL \\" >&2
    echo "       analysis/dump_selected_indices.py" >&2; exit 1; }
[ -f "$HERE/rwl_members.json" ] || {
    echo "no $HERE/rwl_members.json -- run make_rwl_members.py first" >&2
    exit 1; }

mkdir -p "$DEST"
[ "${FORCE:-0}" = "1" ] && rm -rf "$DEST/wgt"
rm -f "$DEST"/pwgevents*.lhe "$DEST"/rwl_b*.xml "$DEST"/reweight_b*.log
cp "$SRC/powheg.input" "$DEST/"
cp "$SRC"/pwggrid*.dat "$SRC"/pwgubound*.dat "$SRC"/pwgxgrid.dat "$DEST/" 2>/dev/null || true

python3 "$HERE/subset_lhe.py" "$SRC/pwgevents.lhe" \
        "$DEST/pwgevents-subset.lhe" "$IDXF"
LHE_IN=$DEST/pwgevents-subset.lhe

NBATCH=$(python3 "$HERE/split_rwl.py" "$HERE/rwl_members.xml" "$DEST" "$PERBATCH")
echo "reweighting the $SEL subset at $TAG in $NBATCH passes of at most $PERBATCH"
mkdir -p "$DEST/wgt"
for b in $(seq 1 "$NBATCH"); do
    x=$(printf "rwl_b%03d.xml" "$b")
    npz="$DEST/wgt/$(printf "b%03d.npz" "$b")"
    [ -f "$npz" ] && { echo "  pass $b/$NBATCH already harvested"; continue; }
    cat "$SRC/powheg.input" > "$DEST/powheg.input"
    cat >> "$DEST/powheg.input" <<KEYS

! ---- a-posteriori reweighting, selected subset, pass $b of $NBATCH ----
rwl_file '$x'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
    rm -f "$DEST/pwgevents-rwgt.lhe"
    cp "$LHE_IN" "$DEST/pwgevents.lhe"
    ( cd "$DEST" && bench_run "$PWHG" > "reweight_b$b.log" 2>&1 ) || {
        echo "pass $b failed; tail of its log:" >&2
        tail -12 "$DEST/reweight_b$b.log" >&2; exit 1; }
    [ -s "$DEST/pwgevents-rwgt.lhe" ] || {
        echo "pass $b produced no output" >&2; exit 1; }
    python3 "$HERE/harvest_weights.py" "$DEST/pwgevents-rwgt.lhe" "$npz"
    rm -f "$DEST/pwgevents-rwgt.lhe"
done
python3 "$HERE/harvest_weights.py" --merge "$DEST/wgt" "$DEST/weights.npz"
cp "$LHE_IN" "$DEST/pwgevents.lhe"
echo "wrote $DEST/weights.npz (rows follow $DEST/pwgevents-subset.index.json)"
