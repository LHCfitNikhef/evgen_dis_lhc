#!/bin/bash
# ONE ENTRY POINT FOR EVERY GATE (TODO item 7, 2026-09-19).
#
#   tools/run_checks.sh            run every checker, print a summary, exit 1 on any FAIL
#   tools/run_checks.sh --build    rebuild results/report.html first (make_report.py --no-copy)
#   tools/run_checks.sh --pull     pull the paper from Overleaf first (CONVENTIONS.md rule 4)
#   tools/run_checks.sh --fast     skip the checks that read event files (parser drift)
#
# Before 2026-09-19 there was no single way to run the gates: the README named
# three of the twelve, one of them failed with no arguments (its default
# sample had been deleted) and another reported "cannot run" as a pass.  Each
# checker here is run with no arguments, and its exit code decides:
#   0 = PASS, 77 = SKIP (the checker says it has nothing to test), else FAIL.
#
# ORDER MATTERS for two groups, and the script says so rather than guessing:
#   * the PAGE checks read results/report.html -- if the page is older than a
#     result JSON or than make_report.py, they are checking a stale page; the
#     script warns (use --build);
#   * the PAPER checks read paper/*.tex -- if Overleaf is ahead of the local
#     copy they check a superseded text; the script prints the sync status
#     (use --pull).
#
# bash 3.2 (macOS /bin/bash): no associative arrays.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/.." && pwd)
cd "$REPO" || exit 2
# the analysis interpreter (config.sh: $BENCH_PYTHON, not a bare python3)
. "$REPO/config.sh"

BUILD=0; PULL=0; FAST=0
for a in "$@"; do
    case "$a" in
        --build) BUILD=1 ;;
        --pull) PULL=1 ;;
        --fast) FAST=1 ;;
        -h|--help) sed -n 2,26p "$0"; exit 0 ;;
        *) echo "unknown option $a" >&2; exit 2 ;;
    esac
done

if [ $PULL = 1 ]; then
    echo "== pulling the paper from Overleaf"
    tools/sync_overleaf.sh pull || { echo "pull failed" >&2; exit 2; }
fi
if [ $BUILD = 1 ]; then
    echo "== rebuilding results/report.html"
    _log=$(mktemp "${TMPDIR:-/tmp}/run_checks_build.XXXXXX")
    "$BENCH_PYTHON" analysis/make_report.py --no-copy > "$_log" 2>&1 \
        || { tail -20 "$_log"; rm -f "$_log"; exit 2; }
    rm -f "$_log"
fi

# stale-page warning: any result JSON or the page builder newer than the page
PAGE=results/report.html
if [ ! -s "$PAGE" ]; then
    echo "WARNING: $PAGE does not exist -- the page checks will fail (use --build)"
elif [ -n "$(find results results_nu analysis/make_report.py analysis/paper_plots \
            -newer "$PAGE" \( -name '*.json' -o -name '*.py' \) 2>/dev/null | head -1)" ]; then
    echo "WARNING: $PAGE is older than some results or code -- page checks see a stale page (use --build)"
fi
echo "== paper sync status"
tools/sync_overleaf.sh status 2>/dev/null | tail -3

NAMES=""; STATES=""; NFAIL=0
run() {   # run <name> <command...>
    local name=$1; shift
    local out rc
    out=$("$@" 2>&1); rc=$?
    local state=PASS
    if [ $rc = 77 ]; then state=SKIP
    elif [ $rc != 0 ]; then state=FAIL; NFAIL=$((NFAIL + 1)); fi
    printf '  %-26s %-4s  %s\n' "$name" "$state" "$(printf '%s\n' "$out" | grep -v '^\s*$' | tail -1 | cut -c1-110)"
    if [ $state = FAIL ]; then printf '%s\n' "$out" | tail -15 | sed 's/^/      | /'; fi
}

echo "== results and code"
run check_no_abs_paths       bash tools/check_no_abs_paths.sh
run check_manifests          "$BENCH_PYTHON" tools/check_manifests.py
run check_figure_labels      "$BENCH_PYTHON" tools/check_figure_labels.py
run check_paper_plots        "$BENCH_PYTHON" tools/check_paper_plots.py
run check_paper_figures_used "$BENCH_PYTHON" tools/check_paper_figures_used.py
run check_report_claims      "$BENCH_PYTHON" tools/check_report_claims.py
run check_paper_claims       "$BENCH_PYTHON" tools/check_paper_claims.py
echo "== the rendered page (results/report.html)"
run check_report_complete    "$BENCH_PYTHON" tools/check_report_complete.py
run check_report_caps        "$BENCH_PYTHON" tools/check_report_caps.py
run check_report_js          "$BENCH_PYTHON" tools/check_report_js.py
run check_yadism_scheme      "$BENCH_PYTHON" tools/check_yadism_scheme.py
run check_powheg_naming      "$BENCH_PYTHON" tools/check_powheg_naming.py
echo "== parsers (read event files)"
if [ $FAST = 1 ]; then
    printf '  %-26s %-4s  %s\n' check_parser_drift SKIP "--fast"
else
    run check_parser_drift   "$BENCH_PYTHON" tools/check_parser_drift.py --max=3000
fi
run check_subset_reweight    "$BENCH_PYTHON" tools/check_subset_reweight.py

echo
if [ $NFAIL = 0 ]; then echo "all gates pass"; exit 0; fi
echo "$NFAIL gate(s) FAIL"; exit 1
