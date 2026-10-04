#!/usr/bin/env bash
# Reclaim the intermediates whose PRODUCT is kept -- rule 1, applied to what
# the old production directories accumulated before the drivers learned to
# clean up after themselves.
#
# This is not a prune of samples.  Nothing here is an event sample and nothing
# here is read by any analysis; every item is a by-product sitting beside the
# file that superseded it.  `tools/analyse_and_prune.py` is still the only
# thing that deletes a sample, and it still refuses the anchor.
#
# THREE CLASSES, and why each is safe:
#
#   1. genie.log -- GENIE's per-event "Interaction Summary" dump, 40 GB over
#      180 files, and the single largest thing on this disk.  It is 3.7x the
#      size of the events it describes: nucclo_job_400GeV_9 carries a 1.29 GB
#      log beside a 344 MB HepMC.  NOTHING READS IT.  The only *.log any
#      analysis opens are Sherpa's integ.log and sherpa.log.
#
#      THE DIGEST IS NOT NEW HERE.  tools/genie_dimuon_ladder.sh met this
#      first and tools/donut_genie.sh repeated it: grep the errors, keep the
#      last 200 lines, write genie_log_digest.txt, drop the log.  The old
#      production directories simply predate that code.  donut_genie.sh's own
#      comment records what not reusing it cost -- 21 GB, and the volume to
#      within 3 GB of the 150 GB floor.  This applies the SAME recipe to every
#      directory that never got it, so the diagnostic content survives.
#
#   2. events.ghep.root -- GENIE's native output, converted ONCE to
#      events.hepmc by gtohepmc3, which is what every analysis parses.  Six GB
#      over 140 files, read by nothing: the only GENIE ROOT file any analysis
#      opens is events.gst.root, in genie_nondis_diff.py, and those live in
#      the ndfull_*/ndonly_* directories alone and are NOT touched here.
#      A ghep.root is deleted only when a NON-EMPTY events.hepmc sits beside
#      it -- the conversion's product present is the whole licence to drop it.
#
#   3. Herwig *.run -- the compiled run card, 1.3 GB over 171 files,
#      regenerated from the tracked *.in by `Herwig read` in seconds.
#
# WHAT IT DELIBERATELY DOES NOT TOUCH: events.hepmc, evtfull, pwgevents.lhe
# (the ladder LHE are re-read by powheg_v2_dimuon_ladder_pdf.sh and
# powheg_v2_faser_ladder_rwgt.sh -- live inputs, not leftovers), events.gst.root,
# events_xsec.json, shower.log, integ.log, sherpa.log, cards or results.
#
#   tools/reclaim_intermediates.sh            # report only, changes nothing
#   tools/reclaim_intermediates.sh --apply    # do it
#
# Run tools/check_manifests.py afterwards: results record the sizes and mtimes
# of their inputs, and nothing here is an input, so it must stay green.
set -u

resolve () {
  p=$1
  while [ -L "$p" ]; do p=$(readlink "$p"); done
  printf '%s/%s\n' "$(cd "$(dirname "$p")" 2>/dev/null && pwd -P)" "$(basename "$p")"
}

SELF=$(resolve "$0")
REPO=$(cd "$(dirname "$SELF")/.." && pwd)
APPLY=0
[ "${1-}" = "--apply" ] && APPLY=1

BENCH_REPO=${BENCH_REPO:-$REPO}
. "$REPO/config.sh"   # GENIE_DIR and friends; BENCH is this repo

TOTAL_KB=0
NFILES=0

note () {   # note <file>  -- tally it, and say so when only reporting
  kb=$(du -sk "$1" 2>/dev/null | awk '{print $1}')
  TOTAL_KB=$((TOTAL_KB + kb))
  NFILES=$((NFILES + 1))
  [ "$APPLY" = "1" ] || printf '  %8s MB  %s\n' "$((kb / 1024))" "${1#$BENCH/}"
}

# --- 1. genie.log -> genie_log_digest.txt ----------------------------------
# The recipe is donut_genie.sh:116-122 verbatim, so the digests this writes
# are the same artefact the newer drivers produce.  The rm is chained onto the
# digest with && : a failed digest keeps the log.
echo "== genie.log (digest kept, log dropped) =="
find "$BENCH" -name genie.log -type f | sort | while read -r f; do
  d=$(dirname "$f")
  note "$f"
  [ "$APPLY" = "1" ] || continue
  ( grep -n -E "FATAL|ERROR|CRIT|WARN|Warning|rror" "$f" \
      | grep -v "cling::AutoLoading" | head -300
    echo "---- last 200 lines of genie.log ($(wc -l < "$f") lines, \
$(wc -c < "$f") bytes) ----"
    tail -n 200 "$f" ) > "$d/genie_log_digest.txt" \
    && rm -f "$f" && echo "  digested $(echo "${d#$BENCH/}")"
done

# --- 2. events.ghep.root, ONLY beside a non-empty events.hepmc -------------
echo "== events.ghep.root (converted; events.hepmc present) =="
find "$BENCH" -name '*.ghep.root' -type f | sort | while read -r f; do
  d=$(dirname "$f")
  if [ ! -s "$d/events.hepmc" ]; then
    echo "  KEPT, no events.hepmc beside it: ${d#$BENCH/}"
    continue
  fi
  note "$f"
  [ "$APPLY" = "1" ] && rm -f "$f" && echo "  removed ${f#$BENCH/}"
done

# --- 3. Herwig *.run ------------------------------------------------------
echo "== Herwig *.run (rebuilt from the tracked *.in) =="
find "$BENCH" -name '*.run' -type f | sort | while read -r f; do
  note "$f"
  [ "$APPLY" = "1" ] && rm -f "$f" && echo "  removed ${f#$BENCH/}"
done

# The while loops above run in subshells, so TOTAL_KB does not survive them
# (bash 3.2: no lastpipe).  Recompute the standing total instead of carrying
# it -- after --apply this is what is LEFT, which is the honest number.
echo
echo "-- still on disk after this run --"
for pat in genie.log '*.ghep.root' '*.run'; do
  find "$BENCH" -name "$pat" -type f -exec du -sk {} + 2>/dev/null \
    | awk -v p="$pat" '{s+=$1; n++} END {printf "  %7.1f GB  %4d  %s\n", s/1048576, n, p}'
done
df -h "$PHYSICS24" | tail -1 | awk '{printf "\n  volume: %s free of %s\n", $4, $2}'
