#!/usr/bin/env bash
# Nuclear-PDF weights for the SELECTED events of the 1 TeV POWHEG-V2 proton
# sample, for the hadron-level nPDF bands of paper plot A2b
# (analysis/npdf_hadron.py).  NEUTRINO CURRENT ONLY, by user decision
# (2026-10-07): the figure is an illustration, and one current suffices.
#
# Usage: powheg/reweight/run_reweight_npdf.sh [weights_per_pass=6] [max_parallel=6]
#
# Input: the selected-event lists results_nu/npdf_weights/indices_nu_p_job<N>[c<k>].json
# written by `analysis/npdf_hadron.py dump <N>` (and split into chunks
# <N>c<k> by hand: one POWHEG-V2 file holds a single seed, read forward, so
# chunks of it reweight independently and in parallel).  Each list's events
# are cut out of the job's own Les Houches file (named in its shower.log;
# subset_lhe.py, exact: see its docstring) and reweighted with
# rwl_npdf_hadron.xml in passes of a few members.  Each pass's reweighted LHE
# is harvested at once and deleted.
#
# THREE TRAPS, all met on 2026-10-07:
#  * SIX PDF members per pass, not twelve: pwhg_main aborts (SIGABRT, or
#    SIGSEGV) while loading roughly its tenth PDF set (README: "about
#    twenty" for the members of one set; fewer when the sets differ).
#  * SMALL FILES ARE COPIED, NOT SYMLINKED.  pwhg_main rewrites FlavRegList,
#    the *_equiv tables, pwhg_checklimits and pwgcounters at start-up; through
#    a symlink those writes land in the PRODUCTION directory, and parallel
#    passes writing the one shared FlavRegList crash.
#  * NEVER CONCATENATE SUBSETS OF DIFFERENT SEEDS into one file (POWHEG-RES
#    manyseeds): each event's #rwgt line restores the random state, and
#    POWHEG-BOX-RES/random.f can only skip forward within one seed --
#    "setrandom: failed".  Irrelevant for the single-seed POWHEG-V2 file, kept
#    here for whoever extends this to the muon current.
#
# Output per list: results_nu/npdf_weights/weights_nu_p_job<N>[c<k>].npz and
# its .index.json (row k -> lhe_index[k], subset_lhe.py's manifest), never
# reconstructed from order.
#
# Safe to re-invoke: harvested passes and finished lists are skipped.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$BENCH_REPO/powheg/production/common.sh"
# the nCTEQ15HQ id lives in the repo's own pdfsets.index (make_rwl_npdf_hadron.py)
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
PY=${BENCH_PYTHON:-python3}
RWL=$HERE/rwl_npdf_hadron.xml
WD=$BENCH_REPO/results_nu/npdf_weights

one () {   # $1 list name <N>[c<k>], $2 perpass
    local n=$1 per=$2 job lhe rd dest nb b x npz out f
    [ -s "$WD/weights_nu_p_job$n.npz" ] && { echo "  skip job $n (done)"; return 0; }
    job=$V2_POWHEG_DIR/v2nu_p_job_${n%%c*}
    lhe=$(grep -o 'Beams:LHEF *| *[^ |]*' "$job/shower.log" | awk '{print $NF}')
    [ -s "$lhe" ] || { echo "  no LHE for $job" >&2; return 1; }
    rd=$(dirname "$lhe")
    grep -q '^storeinfo_rwgt 1' "$rd/powheg.input" || { echo "$rd: no storeinfo_rwgt 1" >&2; return 1; }
    dest=$POWHEG_V2/npdfrwgt-nu1TeV-p/job$n
    mkdir -p "$dest/wgt"
    for f in "$rd"/*; do
        case "$(basename "$f")" in
            powheg.input|pwgevents*|*.log|DEAD_SEEDS|V2_SUMMARY.txt|INTEGRATED) ;;
            *) [ -e "$dest/$(basename "$f")" ] && continue
               if [ "$(wc -c < "$f")" -lt 1000000 ]; then cp "$f" "$dest/"
               else ln -s "$f" "$dest/"; fi ;;
        esac
    done
    [ -s "$dest/subset.index.json" ] ||
        "$PY" "$HERE/subset_lhe.py" "$lhe" "$dest/subset.lhe" "$WD/indices_nu_p_job$n.json" > /dev/null
    nb=$("$PY" "$HERE/split_rwl.py" "$RWL" "$dest" "$per")
    for b in $(seq 1 "$nb"); do
        x=$(printf "rwl_b%03d.xml" "$b")
        npz=$dest/wgt/$(printf "b%03d.npz" "$b")
        [ -s "$npz" ] && continue
        cp "$rd/powheg.input" "$dest/powheg.input"
        cat >> "$dest/powheg.input" <<KEYS

! ---- a-posteriori nPDF reweighting (powheg/reweight/run_reweight_npdf.sh), pass $b of $nb ----
rwl_file '$x'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS
        rm -f "$dest"/pwgevents-rwgt*.lhe
        ln -sf subset.lhe "$dest/pwgevents.lhe"
        ( cd "$dest" && bench_run "$POWHEG_V2/pwhg_main" > "reweight_b$b.log" 2>&1 ) || {
            echo "  job $n pass $b failed:" >&2; tail -5 "$dest/reweight_b$b.log" >&2; return 1; }
        out=$dest/pwgevents-rwgt.lhe
        [ -s "$out" ] || { echo "  job $n pass $b: no reweighted LHE" >&2; return 1; }
        "$PY" "$HERE/harvest_weights.py" "$out" "$npz" > /dev/null
        rm -f "$out" "$dest/reweight_b$b.log"
    done
    "$PY" "$HERE/harvest_weights.py" --merge "$dest/wgt" "$WD/weights_nu_p_job$n.npz"
    cp "$dest/subset.index.json" "$WD/weights_nu_p_job$n.index.json"
    rm -f "$dest/subset.lhe"
    echo "  ok job $n $(date +%H:%M)"
}

if [ "$1" = "--one" ]; then shift; one "$@"; exit $?; fi
PERPASS=${1:-6}
MAXPAR=${2:-6}
[ -s "$RWL" ] || { echo "no $RWL -- run make_rwl_npdf_hadron.py" >&2; exit 1; }
lists=()
for idxf in "$WD"/indices_nu_p_job*.json; do
    n=$(basename "$idxf" .json); lists+=("${n##*job}")
done
echo "nPDF reweighting nu: lists ${lists[*]}, $PERPASS weights per pass, $MAXPAR in parallel"
for j in "${lists[@]}"; do echo "$j $PERPASS"; done |
    xargs -n 2 -P "$MAXPAR" bash "$0" --one
echo "== run_reweight_npdf finished =="
