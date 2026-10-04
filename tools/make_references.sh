#!/usr/bin/env bash
# THE ANALYTIC REFERENCES AND SCALE BANDS FOR PAPER PLOTS (PRODUCTION.md).
#
# Usage:  tools/make_references.sh [phase] [max_parallel]
#           phase: grids | refs | combine | bands | all | list  (default all;
#                  list prints the grid and ref tasks without running them)
#           max_parallel: default 3 -- generator productions share the machine
#         Knobs: REGION (default q4w3; q4w5 is the charm region, W > 5 GeV,
#                user 2026-09-14), TMCS (default "3 0": TMC on, then off).
#                The grids do not depend on the region, so a new region needs
#                only refs, combine and bands.
#
# Region q4w3 (Q2 > 4, W > 3 GeV, no y cut), target W (per-nucleon tungsten
# through NNPDF40_nnlo_as_01180_W184free: structure functions are linear in the
# PDF, so this IS (74 sigma_p + 110 sigma_n)/184), the five BENCH_ENERGIES, and
# every reference TWICE: target-mass corrections off (BENCH_TMC=0) and on
# (BENCH_TMC=3, the exact Georgi-Politzer form; see analysis/target.py).
#
# WHAT EACH PHASE DOES, in dependency order:
#   grids    the PDF-independent PineAPPL grids, one per (current, energy,
#            order, scheme, TMC).  The massive FONLL components are charm
#            only.  Already-built grids are skipped, and the TMC-off grids at
#            400/1000/4000 GeV exist from the earlier production -- a grid does not know the PDF
#            or the region, so they serve tungsten and q4w3 unchanged.
#   refs     the direct YADISM runs in ZM-VFNS: inclusive NLO and NNLO and
#            charm NLO (and NNLO on the muon side), both currents.
#   combine  charm in FONLL (convolved from the grids), then the FONLL
#            inclusive = inclusive ZM - charm ZM + charm FONLL.
#   bands    the 7-point bands, one run per (current, order, scheme, TMC)
#            over all five energies.
# grids and refs are independent and share one pool of max_parallel slots.
#
# NO FONLL AT NNLO ON THE NEUTRINO SIDE, and it is not asked for: the massive
# charged-current coefficient functions are known but not implemented in YADISM, and
# mhou_sigma_fonll.py refuses the request outright.
#
# Every task SKIPS work whose output already exists (CONVENTIONS.md 1c), so the
# driver can be killed and re-invoked at the cost of the unfinished tasks.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
cd "$BENCH_REPO"

PHASE=${1:-all}
MAXPAR=${2:-3}
PY=${BENCH_PYTHON:-python3}
ES=${ES:-"400 700 1000 2000 4000"}   # knob: ES=300 adds the SIDIS-ladder point (beams.SIDIS_ENERGIES)
TMCS=${TMCS:-"3 0"}
REGION=${REGION:-q4w3}
case "$REGION" in q4w3|q4w5) ;; *) echo "REGION must be q4w3 or q4w5" >&2; exit 1 ;; esac
LOG="$BENCH_REPO/logs/v2refs"
[ "$REGION" = q4w3 ] || LOG="$BENCH_REPO/logs/v2refs_$REGION"
mkdir -p "$LOG"

# KEEP THE MACHINE AWAKE FOR THE WHOLE DRIVER.  The first run of this script
# (2026-09-13/14) had no keep-awake and the Mac idle-slept from 01:01, pausing
# every task; `caffeinate -i -w` added by hand did NOT stop the "Sleep Service
# Back to Sleep" cycles, `-s` (no system sleep on AC) did.  -w ties the
# assertion to this shell, so it ends with the driver.  config.sh's
# BENCH_NOSLEEP is empty on Linux, which skips this.
case "$BENCH_NOSLEEP" in
  caffeinate*) $BENCH_NOSLEEP -s -w $$ > /dev/null 2>&1 & ;;
esac

export BENCH_SELECTION=$REGION BENCH_TARGET=W
unset BENCH_Q2MIN
# one core per task: yadism's numerics must not fan out behind the throttle
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
       MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

tsuf() { [ "$1" = "0" ] && echo "" || echo "_tmc"; }
# result file for (resdir, stem, energy, tmc), named as the calculators name it
rfile() {
  echo "$1/$($PY analysis/beams.py name "histos_$2_${REGION}_W$(tsuf "$4")" "$3").json"
}

run_pool() {   # stdin: one "logname@@command" per line
  # bash 3.2 has no `wait -n`; xargs -P keeps max_parallel slots full.  NUL
  # delimited, so xargs does no quote processing on the commands.
  tr '\n' '\0' | xargs -0 -n 1 -P "$MAXPAR" sh -c '
    [ -n "$1" ] || exit 0
    name=${1%%@@*}; cmd=${1#*@@}
    echo "  start $name"
    if sh -c "$cmd" > "$0/$name.log" 2>&1; then echo "  done  $name"
    else echo "  FAILED $name -- see $0/$name.log"; fi' "$LOG"
}

grid_tasks() {
  for T in $TMCS; do for E in $ES; do
    # the massive components first: they are the long ones
    for spec in "mu 1 FONLL-FFNS 3 charm" "mu 1 FONLL-FFN0 3 charm" \
                "mu 2 FONLL-FFNS 3 charm" "mu 2 FONLL-FFN0 3 charm" \
                "nu 1 FONLL-FFNS 3 charm" "nu 1 FONLL-FFN0 3 charm" \
                "mu 1 ZM-VFNS 5 all" "mu 2 ZM-VFNS 5 all" \
                "nu 1 ZM-VFNS 5 all" "nu 2 ZM-VFNS 5 all"; do
      set -- $spec
      printf 'grid_%s_%s_pto%s_%s_tmc%s@@%s\n' "$1" "$E" "$2" "$3" "$T" \
        "BENCH_TMC=$T $PY analysis/pineappl_grids.py build $1 $E $2 $3 $4 $5"
    done
  done; done
}

ref_tasks() {
  for T in $TMCS; do for E in $ES; do
    env="BENCH_TMC=$T BENCH_ENERGY=$E"
    for o in nlo nnlo; do
      f=$(rfile results "yadism_$o" "$E" "$T")
      [ -s "$f" ] || printf 'mu_%s_%s_tmc%s@@%s\n' "$o" "$E" "$T" \
        "$env $PY analysis/yadism_calc.py $o"
      f=$(rfile results "yadism_charm_$o" "$E" "$T")
      [ -s "$f" ] || printf 'mu_charm_%s_zm_%s_tmc%s@@%s\n' "$o" "$E" "$T" \
        "$env $PY analysis/yadism_charm_calc.py charm $o zm"
      f=$(rfile results_nu "yadism_$o" "$E" "$T")
      [ -s "$f" ] || printf 'nu_%s_%s_tmc%s@@%s\n' "$o" "$E" "$T" \
        "$env $PY analysis/yadism_cc_calc.py $o"
    done
    f=$(rfile results_nu "yadism_charm_nlo" "$E" "$T")
    [ -s "$f" ] || printf 'nu_charm_nlo_zm_%s_tmc%s@@%s\n' "$E" "$T" \
      "$env $PY analysis/yadism_cc_charm_calc.py nlo zm"
  done; done
}

combine() {
  for T in $TMCS; do for E in $ES; do
    export BENCH_TMC=$T BENCH_ENERGY=$E
    for o in nlo nnlo; do
      f=$(rfile results "yadism_charm_${o}_fonll_damp" "$E" "$T")
      [ -s "$f" ] || $PY analysis/yadism_charm_calc.py charm $o fonll_damp \
        > "$LOG/mu_charm_${o}_fonll_${E}_tmc$T.log" 2>&1
    done
    f=$(rfile results_nu "yadism_charm_nlo_fonll_damp" "$E" "$T")
    [ -s "$f" ] || $PY analysis/yadism_cc_charm_calc.py nlo fonll_damp \
      > "$LOG/nu_charm_nlo_fonll_${E}_tmc$T.log" 2>&1
    # both currents in one call; the neutrino NNLO has no FONLL charm input
    # and is skipped there by construction
    $PY analysis/yadism_fonll_inclusive.py nlo  > "$LOG/fonll_incl_nlo_${E}_tmc$T.log" 2>&1
    $PY analysis/yadism_fonll_inclusive.py nnlo > "$LOG/fonll_incl_nnlo_${E}_tmc$T.log" 2>&1 || true
    echo "  combined $E GeV, TMC $T"
    unset BENCH_TMC BENCH_ENERGY
  done; done
}

bands() {
  EL=$(echo $ES | tr ' ' ',')
  for T in $TMCS; do
    for args in "--current mu --zm" "--current mu" "--current mu --pto 2" \
                "--current nu --zm" "--current nu" "--current nu --pto 2 --zm"; do
      tag=$(echo "$args" | tr -d '-' | tr ' ' '_')
      printf 'band%s_tmc%s@@%s\n' "$tag" "$T" \
        "BENCH_TMC=$T $PY analysis/mhou_sigma_fonll.py $args --energies $EL"
    done
  done | run_pool
}

case "$PHASE" in
  grids)   grid_tasks | run_pool ;;
  refs)    ref_tasks | run_pool ;;
  combine) combine ;;
  bands)   bands ;;
  list)    grid_tasks; ref_tasks ;;
  all)     { grid_tasks; ref_tasks; } | run_pool; combine; bands ;;
  *) echo "usage: $0 [grids|refs|combine|bands|all] [max_parallel]" >&2; exit 1 ;;
esac
echo "== references: phase $PHASE finished =="
