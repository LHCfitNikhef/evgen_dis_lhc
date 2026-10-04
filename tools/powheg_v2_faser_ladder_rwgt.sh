#!/usr/bin/env bash
# Add the 7-point scale weights to every point of the POWHEG-V2 FASER ladder
# (user, 2026-09-04: "add the MHOUs on top of the POWHEG prediction").
#
# The ladder points were generated with `storeinfo_rwgt 1`, so the matrix
# element can be re-evaluated at each scale a posteriori -- the same
# reweighting pass as powheg/reweight/run_reweight.sh, applied to the 40
# short runs of tools/powheg_v2_faser_ladder.sh instead of a production.
# Scale weights ONLY (powheg/reweight/rwl_scale.xml): the neutron points run
# on the isospin-swapped set, for which the proton PDF ids of the fuller file
# would be meaningless.
#
# Works in a `rwgt/` subdirectory of each point, never in place, and leaves
# `rwgt/pwgevents-rwgt.lhe`, which analysis/faser_powheg_rates.py prefers
# over the plain LHE when it exists.  About 20 s per point.
#
# Usage: tools/powheg_v2_faser_ladder_rwgt.sh            (MAXJOBS=4 parallel)
#        FORCE=1 tools/powheg_v2_faser_ladder_rwgt.sh    (redo existing)
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PWHG=$POWHEG_V2/pwhg_main
BASEDIR=$POWHEG_V2/ladder-faser
RWL=$BENCH_REPO/powheg/reweight/rwl_scale.xml
MAXJOBS=${MAXJOBS:-4}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 2; done; }

one () {
  R=$1
  D=$R/rwgt
  if [ -s "$D/pwgevents-rwgt.lhe" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "exists: $D"; return; fi
  rm -rf "$D"; mkdir -p "$D"
  cp "$R/powheg.input" "$R/pwgevents.lhe" "$D/"
  cp "$R"/pwggrid*.dat "$R"/pwgubound*.dat "$R"/pwgxgrid.dat "$D/" 2>/dev/null || true
  cp "$RWL" "$D/"
  # rwl_add 1 is what SELECTS the reweighting path (see run_reweight.sh)
  cat >> "$D/powheg.input" <<KEYS

! ---- a-posteriori scale reweighting (tools/powheg_v2_faser_ladder_rwgt.sh) ----
rwl_file 'rwl_scale.xml'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
  if ( cd "$D" && bench_run "$PWHG" > reweight.log 2>&1 ) && [ -s "$D/pwgevents-rwgt.lhe" ]; then
    echo "done: $D ($(grep -c '<wgt' "$D/pwgevents-rwgt.lhe") weights)"
  else
    echo "FAILED: $D -- see $D/reweight.log" >&2
  fi
}

for R in "$BASEDIR"/nu_*_E* "$BASEDIR"/nubar_*_E*; do
  [ -s "$R/pwgevents.lhe" ] || continue
  limit
  one "$R" &
done
wait
echo "ladder reweighting finished"
