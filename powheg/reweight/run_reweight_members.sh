#!/usr/bin/env bash
# Reweight a POWHEG-V2 sample over EVERY member of every PDF set, in batches.
#
# Usage: ./run_reweight_members.sh [--current mu|nu] [--mc] [--rwl FILE] \
#                                  <energy_gev> [members_per_batch]
#
#   run_reweight_members.sh 1000                       # neutrino CC, the six
#                                                      # proton sets (as before)
#   run_reweight_members.sh --current mu --rwl rwl_nuclear 1000
#                                                      # muon NC, the nuclear
#                                                      # sets -> rwgtnuclear-mu1TeV
#
# --rwl names a reweight file in this directory without its extension, and a
# non-default one writes into ITS OWN destination directory: two studies over
# the same sample must never share a weights.npz.
#
# WHY BATCHES.  Asking pwhg_main for all 309 members at once ABORTS: it loads
# them one by one and dies with a bare SIGABRT after about the twentieth, with
# no error message of its own -- something internal overflows.  A batch of ten
# runs cleanly, so the members are added ten at a time and `rwl_add 1` makes
# each pass APPEND to the weights already in the file.  The limit was found by
# running into it, not from documentation, so the batch size is deliberately
# well under where it broke.
#
# WHY IT IS AFFORDABLE AT ALL.  The README beside this script says PDF error
# sets "would cost about ten times this whole run", which was an estimate and
# never a measurement.  Measured: 87 microseconds per (event, weight), so the
# 309 members over 50000 events is roughly half an hour of compute per beam
# energy.  What makes it slower than that in practice is I/O -- every pass
# rewrites the whole LHE, which grows as weights accumulate.
#
# THE RESULT is a Les Houches file carrying the seven scale points, the four
# central members under their original ids, and every member of all six sets,
# so a PDF UNCERTAINTY can be formed from the Monte Carlo and compared with
# the analytic one from YADISM.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$BENCH_REPO/powheg/powheg_variants.sh"

# --mc / BENCH_POWHEG_CODE=V2MC selects the massive-charm entry, which has
# its own run directory and must never be reweighted into the massless one's.
PV_CODE=${BENCH_POWHEG_CODE:-V2}
# THE CURRENT IS A PARAMETER, as on run_reweight.sh and the three generation
# drivers (CONVENTIONS.md rule 2b: what one current gets, the other is owed).  The
# V2 muon cross-variant is the same CODE, so the reweighting machinery reaches
# it unchanged; POWHEG-RES still cannot be reweighted here, for the structural
# reason in run_reweight.sh's header.
CURRENT=nu
# --rwl selects the reweight file, WITHOUT its extension.  A file other than
# the default writes into its OWN destination directory, so a study can never
# overwrite the weights.npz another study's published results were made from.
RWL=rwl_members
ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --mc)      PV_CODE=V2MC; shift ;;
        --current) CURRENT=${2:?--current takes mu or nu}; shift 2 ;;
        --rwl)     RWL=${2:?--rwl takes a reweight file base name}
                   RWL=${RWL%.xml}; shift 2 ;;
        *)         ARGS+=("$1"); shift ;;
    esac
done
set -- "${ARGS[@]}"
EBEAM=${1:-1000}
# SIX, not ten.  The ceiling is on PDFs loaded in ONE pwhg_main process, not
# on weights in the file: a pass that loaded 10 survived and one that loaded 12
# aborted.  Six leaves margin, at the cost of more passes.
PERBATCH=${2:-6}
case "$CURRENT" in mu|nu) ;; *) echo "--current takes mu or nu" >&2; exit 1 ;; esac
[ "$PV_CODE" = V2MC ] && [ "$CURRENT" != nu ] && {
    echo "POWHEG-V2mc is neutrino CC only" >&2; exit 1; }

# The derived PDF sets this repository builds -- the isospin-mirrored neutrons
# and the free-nucleon tungsten averages -- live here and are reached by
# PREPENDING the data path, exactly as the ladders do it.
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH

powheg_variant "$PV_CODE" "$CURRENT" "$EBEAM"
TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$CURRENT" "$EBEAM")
SRC=$PV_RUNDIR
# THE DESTINATION CARRIES THE CODE AND THE REWEIGHT FILE, so a V2mc harvest can
# never land in the massless entry's directory and be picked up as if it were
# the same sample, and a second study's weights can never replace the first's.
DEST=$POWHEG_V2/rwgt-$CURRENT$TAG
[ "$PV_CODE" = V2MC ] && DEST=$POWHEG_V2/rwgt-numc$TAG
[ "$RWL" != "rwl_members" ] && DEST=$POWHEG_V2/rwgt${RWL#rwl_}-$CURRENT$TAG
PWHG=$POWHEG_V2/pwhg_main

[ -s "$SRC/pwgevents.lhe" ] || { echo "no LHE at $SRC" >&2; exit 1; }

# CAP THE EVENTS HANDED TO THE REWEIGHTER (default 50000, matching the 1 TeV
# and 4 TeV samples).  The 400 GeV sample carries 1,000,000 events because the
# dimuon tier needs them, and at the measured ~5900 events/minute one pass over
# it takes ~170 min -- 54 passes is six and a half days.  The PDF band is a
# spread across members on the SAME events, so the per-event MC error largely
# cancels and 50000 is what the other two energies give it anyway.
# BENCH_RWGT_MAXEV=0 reweights everything.
MAXEV=${BENCH_RWGT_MAXEV:-50000}
[ -f "$HERE/$RWL.json" ] || {
    echo "no $HERE/$RWL.json -- run make_${RWL#rwl_}... first" >&2
    exit 1; }

# KEEP THE HARVESTED WEIGHTS.  This used to `rm -rf "$DEST"`, which deleted
# wgt/ along with everything else and so made the per-pass resume check below
# dead code: an interrupted run -- and this one takes hours -- started again
# from nothing.  Only the scratch is cleared now, so a rerun picks up where it
# stopped.  Pass FORCE=1 to discard the harvest and redo every pass.
mkdir -p "$DEST"
[ "${FORCE:-0}" = "1" ] && rm -rf "$DEST/wgt"
rm -f "$DEST"/pwgevents*.lhe "$DEST"/rwl_b*.xml "$DEST"/reweight_b*.log
cp "$SRC/powheg.input" "$DEST/"
# the capped copy is made ONCE, here, and every pass then copies from it
if [ "$MAXEV" != "0" ]; then
    python3 "$HERE/truncate_lhe.py" "$SRC/pwgevents.lhe" \
            "$DEST/pwgevents-capped.lhe" "$MAXEV"
    LHE_IN=$DEST/pwgevents-capped.lhe
else
    LHE_IN=$SRC/pwgevents.lhe
fi
cp "$LHE_IN" "$DEST/pwgevents.lhe"
cp "$SRC"/pwggrid*.dat "$SRC"/pwgubound*.dat "$SRC"/pwgxgrid.dat "$DEST/" 2>/dev/null || true

# INDEPENDENT PASSES, NOT ACCUMULATION.  pwhg_main aborts once the file it is
# reading carries about twenty weights: batch 1 wrote ten and batch 2 then
# died partway through its own ten.  So each pass reads the ORIGINAL
# pwgevents.lhe, writes its own output, and the weights are harvested into a
# compact array afterwards -- nothing accumulates in any LHE.
NBATCH=$(python3 "$HERE/split_rwl.py" "$HERE/$RWL.xml" "$DEST" "$PERBATCH")
echo "reweighting [$PV_LABEL] $TAG with $RWL in $NBATCH independent passes of at most $PERBATCH"

mkdir -p "$DEST/wgt"
for b in $(seq 1 "$NBATCH"); do
    x=$(printf "rwl_b%03d.xml" "$b")
    npz="$DEST/wgt/$(printf "b%03d.npz" "$b")"
    if [ -f "$npz" ]; then
        echo "  pass $b/$NBATCH already harvested"
        continue
    fi
    cat "$SRC/powheg.input" > "$DEST/powheg.input"
    cat >> "$DEST/powheg.input" <<KEYS

! ---- a-posteriori reweighting, independent pass $b of $NBATCH ----
rwl_file '$x'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
    rm -f "$DEST/pwgevents-rwgt.lhe"
    # every pass starts from the pristine file, so none of them ever sees
    # more than $PERBATCH weights
    cp "$LHE_IN" "$DEST/pwgevents.lhe"
    ( cd "$DEST" && bench_run "$PWHG" > "reweight_b$b.log" 2>&1 ) || {
        echo "pass $b failed; tail of its log:" >&2
        tail -12 "$DEST/reweight_b$b.log" >&2; exit 1; }
    [ -s "$DEST/pwgevents-rwgt.lhe" ] || {
        echo "pass $b produced no output" >&2; exit 1; }
    # harvest into a compact array and throw the 60 MB LHE away
    python3 "$HERE/harvest_weights.py" "$DEST/pwgevents-rwgt.lhe" "$npz"
    rm -f "$DEST/pwgevents-rwgt.lhe"
    echo "  pass $b/$NBATCH harvested -> $(basename "$npz")"
done

python3 "$HERE/harvest_weights.py" --merge "$DEST/wgt" "$DEST/weights.npz"
cp "$LHE_IN" "$DEST/pwgevents.lhe"
echo "wrote $DEST/weights.npz"
