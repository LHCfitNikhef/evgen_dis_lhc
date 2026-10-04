#!/usr/bin/env bash
# RE-SHOWER THE  POWHEG-V2 (nu CC) SAMPLES WITH LesHouches:matchInOut = off.
#
# Usage:  powheg/production/reshower_nu_matchinout.sh [max_parallel=6] [job ...]
#         (no job names: all powheg/v2nu_{p,n}_job*; a name is e.g. v2nu_p_job_1)
#
# Why: see powheg/powheg_nu_v2.cmnd.  The LHE files are the production's own
# (each job's shower.log names its Beams:LHEF) -- nothing is regenerated, only
# showered again, with the one card change.
#
# SAFE TO KILL AND RE-INVOKE (CONVENTIONS.md 1c).  Each job showers into
# <job>.reshower/ and is verified there by v2_shower (LHE and HepMC beams, at
# least 90% of the offered events, V2_OK); only then is the old events.hepmc
# replaced, and the job is stamped MATCHINOUT_OFF.  A stamped job is skipped.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"
CMND=$(v2_cmnd nu)
cd "$V2_POWHEG_DIR"

jobs_list () {
    if [ $# -gt 0 ]; then for j in "$@"; do echo "$V2_POWHEG_DIR/$j"; done
    else
        for d in "$V2_POWHEG_DIR"/v2nu_p_job* "$V2_POWHEG_DIR"/v2nu_n_job*; do
            case "$d" in *.reshower) continue ;; esac
            [ -d "$d" ] && echo "$d"
        done
    fi
}

one () {   # $1 = job dir
    local d=$1 name lhe tgt new
    name=$(basename "$d")
    [ -f "$d/MATCHINOUT_OFF" ] && { echo "  skip $name (already re-showered)"; return 0; }
    lhe=$(grep -o 'Beams:LHEF *| *[^ |]*' "$d/shower.log" | awk '{print $NF}')
    [ -s "$lhe" ] || { echo "  no LHE for $name ($lhe)" >&2; return 1; }
    case "$name" in v2nu_p_*) tgt=p ;; v2nu_n_*) tgt=n ;; *) return 1 ;; esac
    new="$d.reshower"
    rm -rf "$new"
    v2_shower "$lhe" "$new" "$CMND" "$tgt" 14 || return 1
    mv -f "$new/events.hepmc" "$d/events.hepmc"
    mv -f "$new/events_xsec.json" "$d/events_xsec.json"
    mv -f "$new/shower.log" "$d/shower.log"
    mv -f "$new/V2_OK" "$d/V2_OK"
    rmdir "$new" 2>/dev/null || rm -rf "$new"
    echo "matchInOut off, $(date), card $(basename "$CMND")" > "$d/MATCHINOUT_OFF"
    echo "  replaced $name"
}
# One job per process: the driver re-invokes itself as "--one <dir>", so every
# worker sources common.sh (and config.sh) itself.  bash 3.2 has no wait -n;
# xargs keeps max_parallel slots full.
if [ "$1" = "--one" ]; then one "$2"; exit $?; fi
MAXPAR=${1:-6}; shift || true
jobs_list "$@" | xargs -n 1 -P "$MAXPAR" "$HERE/reshower_nu_matchinout.sh" --one
echo "== re-shower finished =="
