#!/usr/bin/env bash
# Sherpa shower-model arms: the SAME card and the SAME matrix element run
# through a different parton shower.
#
# Usage: ./run_sherpa_shower.sh <mu|nu> <lo|nlo> <CSS|Dire> [nev_per_job] [njobs]
#
#   ./run_sherpa_shower.sh mu nlo Dire 25000 8      # -> MuonDIS_NLO_Dire/job_N
#   ./run_sherpa_shower.sh nu lo  Dire 25000 8      # -> NuDIS_LO_Dire/job_N
#
# WHY SHERPA AND NOT PYTHIA, for the NLO half of the shower study.  The LO
# study (pythia8/shower/) switches PartonShowers:model, which also switches
# hadronisation tune, so it conflates shower with tune; and it cannot be done
# at NLO at all, because the POWHEG veto is a Pythia UserHooks object that
# intercepts SIMPLE-shower emission scales and Vincia and Dire do not expose
# it (logs/shower_probe_powheg.log).  Sherpa has neither problem: it owns its
# MC@NLO matching, and its two showers share one cluster hadronisation model.
# So the difference between CSS and Dire here IS the shower.
#
# ONLY TWO SHOWERS QUALIFY.  This build also has DIM, which refuses to start
# under MC@NLO -- "Shower needs to be set for MC@NLO" -- rather than running
# with the matching quietly dropped.  It is therefore rejected by name below
# rather than being allowed to fail per-run.
#
# NO SECOND COPY OF THE CARD.  Sherpa takes YAML overrides on the command
# line, so the shower is selected with `"SHOWER_GENERATOR: <mod>"` and the
# committed card stays the single source of every other setting.  The variant
# run directory symlinks that same card, exactly as the base directory does.
#
# EACH ARM INTEGRATES ITSELF, and this was MEASURED, not assumed.  The
# MuonDIS_NLO card integrated from scratch with each shower gives:
#
#            BVI [pb]           RS [pb]          TOTAL [nb]
#   CSS    34850.00 +- 43.74    234.87 +- 6.40   35.0849 +- 0.0442
#   Dire   27624.40 +- 35.28   7488.39 +- 6.05   35.1128 +- 0.0358
#
# The TOTAL agrees to 0.08% (0.5 sigma), as it must -- it is the NLO
# cross-section and cannot depend on the shower.  That is the consistency
# check this whole comparison rests on, and it passes.
#
# THE SPLIT DOES NOT AGREE AT ALL: RS is 0.7% of the total under CSS and 21%
# under Dire, a factor of 32.  That is the MC@NLO subtraction being built from
# the shower's splitting kernels.  So an arm must NOT inherit a Results.zip
# integrated with a different shower -- and note HOW it would fail: the total
# is pinned by the integrator, so the cross-section would still come out
# right, while the S-event/H-event population belonged to one shower and the
# events were showered with another.  Nothing would look wrong.
#
# The compiled Amegic libraries under Process/ ARE shower-independent (they
# are matrix-element code), so those are shared by symlink rather than
# recompiled.
set -e
CUR=${1:-mu}
ORDER=${2:-nlo}
SHOWER=${3:-Dire}
NEV=${4:-25000}
NJOBS=${5:-8}

HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"

case "$CUR" in mu|nu) ;; *) echo "current must be mu or nu" >&2; exit 1;; esac
case "$ORDER" in lo|nlo) ;; *) echo "order must be lo or nlo" >&2; exit 1;; esac
case "$SHOWER" in
  CSS|Dire) ;;
  DIM) echo "DIM has no MC@NLO interface in Sherpa 3.0.5 -- it exits with" \
            "'Shower needs to be set for MC@NLO'.  Use CSS or Dire." >&2
       exit 1 ;;
  *) echo "shower must be CSS or Dire" >&2; exit 1 ;;
esac

# (current, order) -> the committed base run directory.  A TABLE, not a rule:
# the neutrino NLO directory is NuDIS_NLO_ckm3, because NuDIS_NLO is the
# superseded unit-CKM run kept as documentation.  Deriving the name would
# silently pick the wrong one.
case "$CUR:$ORDER" in
  mu:lo)   BASE=MuonDIS_LO      ;;
  mu:nlo)  BASE=MuonDIS_NLO     ;;
  nu:lo)   BASE=NuDIS_LO        ;;
  nu:nlo)  BASE=NuDIS_NLO_ckm3  ;;
esac
VARIANT=${BASE}_${SHOWER}
SRC=$SHERPA_RUNS/$BASE
DST=$SHERPA_RUNS/$VARIANT
REPOCARD=$BENCH_REPO/sherpa/Runs/$BASE/Sherpa.yaml

[ -d "$SRC" ]        || { echo "no base run dir $SRC" >&2; exit 1; }
[ -f "$REPOCARD" ]   || { echo "no committed card $REPOCARD" >&2; exit 1; }
[ -d "$SRC/Process" ] || { echo "no compiled Process/ in $SRC -- integrate the" \
                                "baseline first" >&2; exit 1; }

mkdir -p "$DST"
ln -sfn "$REPOCARD" "$DST/Sherpa.yaml"      # the SAME committed card
[ -e "$DST/Process" ] || ln -s "$SRC/Process" "$DST/Process"

SHERPA_BIN=$SHERPA_INSTALL/bin/Sherpa
[ -x "$SHERPA_BIN" ] || { echo "no Sherpa at $SHERPA_BIN" >&2; exit 1; }

# ---- integrate this arm, once ------------------------------------------
if [ ! -f "$DST/Results.zip" ]; then
    echo "$VARIANT: integrating (no Results.zip yet)"
    ( cd "$DST" && bench_run "$SHERPA_BIN" -e 0 -R 4321 \
        "SHOWER_GENERATOR: $SHOWER" > integ.log 2>&1 ) || {
        echo "integration failed; see $DST/integ.log" >&2; exit 1; }
    echo "$VARIANT: integration done"
else
    echo "$VARIANT: reusing its own Results.zip"
fi

# ---- generate ------------------------------------------------------------
echo "$VARIANT: $NJOBS x $NEV events with SHOWER_GENERATOR: $SHOWER"
for i in $(seq 1 "$NJOBS"); do
  d=$DST/job_$i
  if [ -s "$d/evtfull" ] && [ "${FORCE:-0}" != "1" ]; then
      echo "refusing to overwrite $d/evtfull -- FORCE=1 to replace" >&2
      exit 1
  fi
done
for i in $(seq 1 "$NJOBS"); do
  d=$DST/job_$i
  mkdir -p "$d"
  ln -sfn "$REPOCARD" "$d/Sherpa.yaml"
  [ -e "$d/Process" ] || ln -s "$DST/Process" "$d/Process"
  cp -f "$DST/Results.zip" "$d/"
  ( cd "$d" && bench_run "$SHERPA_BIN" -e "$NEV" \
      "RANDOM_SEED: $((1234 + i))" \
      "SHOWER_GENERATOR: $SHOWER" \
      'EVENT_OUTPUT: ["HepMC3_GenEvent[evtfull]"]' \
      > sherpa.log 2>&1 ) &
done
wait
echo "$VARIANT: all $NJOBS jobs done"

# DELIVERY GUARD, the same lesson as run_powheg_shower.sh: a job can fail on
# every event and still leave this script exiting 0.
short=0
for i in $(seq 1 "$NJOBS"); do
  d=$DST/job_$i
  n=$(grep -c '^E ' "$d/evtfull" 2>/dev/null); n=${n:-0}
  echo "  job_$i: $n events"
  [ "$n" -lt $((NEV / 2)) ] && short=$((short + 1))
done
if [ "$short" -gt 0 ]; then
    echo "$short job(s) delivered under half of $NEV events -- check" \
         "$DST/job_*/sherpa.log" >&2
    exit 1
fi
