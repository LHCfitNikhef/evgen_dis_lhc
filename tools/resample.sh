#!/usr/bin/env bash
# "THE SAMPLE CHANGED": re-derive EVERYTHING that depends on one sample.
#
# Usage:
#   tools/resample.sh powheg_nu   [energies...]     # default: all three
#   tools/resample.sh powheg_nu_mc 1000 4000
#   tools/resample.sh --list                        # what it knows about
#   tools/resample.sh --dry-run powheg_nu           # print, do not run
#
# WHY THIS EXISTS (user, 2026-08-29): "make sure that 'the sample changed' is
# one command rather than many rediscoveries".
#
# Tripling the POWHEG-V2 statistics set off a cascade that was found ONE LINK
# AT A TIME, each by a different safety gate rather than by planning:
#
#   the indexed twin went stale        -> caught by the lhe_index join
#   the member weights went stale      -> caught by the closure gate
#   faser_s and faser_dimuon were never re-analysed (only `inclusive` was)
#                                      -> caught by the closure gate, -72%
#   the dimuon SUBSET indices went stale, leaving 5 events instead of ~120
#                                      -> caught by the statistics floor
#
# Every one was caught and none reached a figure, which is the system working.
# But four rediscoveries of one fact is three too many, and the next person --
# or the next sample -- may not be so well covered.  The dependency is written
# down here ONCE.
#
# "ally" -- Q2 > 4 with no y window -- IS IN EVERY SELECTION LOOP, even though
# it is not the benchmark's region.  It is what the arXiv:2402.13318
# comparison is made on, and a comparison whose generators were re-analysed
# and whose foreign region was not is exactly the stale-by-one-region failure
# this script exists to stop.  Sherpa and Herwig have no events outside the y
# window and simply produce nothing for it.
#
# THE ORDER MATTERS and is not alphabetical: a twin must be re-showered before
# the weights are harvested against it, the weights before the bands, the
# selections before the subsets that are chosen from them.  Each stage below
# depends only on the ones above it.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$REPO/config.sh"
cd "$REPO"

DRY=0
[ "${1:-}" = "--dry-run" ] && { DRY=1; shift; }
run() { echo "  \$ $*"; [ "$DRY" = 1 ] || "$@"; }

SAMPLE=${1:-}
case "$SAMPLE" in
  --list|"")
    cat <<'LIST'
samples this knows how to re-derive:

  powheg_nu      POWHEG-V2 neutrino CC (massless charm)
                 twin -> member weights -> pdf members -> all four selections
                 -> selection dependence -> dimuon subset -> dimuon bands
  powheg_nu_mc   POWHEG-V2mc neutrino CC (massive charm)
                 twin -> member weights -> pdf members -> selections
  powheg         POWHEG-RES muon NC
                 shower -> all four selections (no reweighting chain)
  genie_nu       GENIE neutrino CC, both tunes: all four selections
  mg5_ps         MG5_aMC LO + Pythia 8, both currents (1 TeV only)
                 parton level from the same LHE -> every selection

after ANY of them, finish with:
  tools/check_manifests.py             (every result reads its WHOLE sample)
  tools/make_all_plots.sh && analysis/make_report.py
  tools/sync_paper_figures.sh          (if a paper plot moved)
  tools/check_paper_claims.py          (the paper's numbers)
  tools/export_faser_emulsion_csv.py --check
  tools/export_sidis_share.py --check
                                       (the CSVs shared outside the repo)
LIST
    exit 0 ;;
esac

shift
ES=${*:-"1000 4000 400"}
echo "== re-deriving everything downstream of $SAMPLE, energies: $ES =="

case "$SAMPLE" in
powheg_nu|powheg_nu_mc)
  MC=""; [ "$SAMPLE" = powheg_nu_mc ] && MC="--mc"
  for E in $ES; do
    echo "-- $E GeV: indexed twin (the join needs lhe_index on every event)"
    run ./powheg/reshower_v2_indexed.sh "$E"
  done
  for E in $ES; do
    echo "-- $E GeV: member weights (tied to the LHE, so they follow it)"
    run env FORCE=1 ./powheg/reweight/run_reweight_members.sh "$E" 6 $MC
  done
  for E in $ES; do
    echo "-- $E GeV: inclusive PDF bands"
    run env BENCH_ENERGY="$E" "$BENCH_PYTHON" analysis/powheg_pdf_members.py
  done
  KEYS="$SAMPLE ${SAMPLE}_charmfinal"
  for E in $ES; do
    for S in inclusive faser_s faser_e faser_dimuon ally; do
      for K in $KEYS; do
        echo "-- $E GeV $S: $K"
        run env BENCH_ENERGY="$E" BENCH_SELECTION="$S" \
            "$BENCH_PYTHON" analysis/analyze_nu.py "$K" || true
      done
    done
  done
  if [ "$SAMPLE" = powheg_nu ]; then
    # THE NUCLEAR-PDF STUDY reads its own weights, harvested from the same LHE
    # into their own directory, so it goes stale with the sample exactly as the
    # member weights do.  Anchor energy only -- that is where the study is
    # defined.  ITS MUON HALF IS NOT REACHED FROM HERE: it runs on the POWHEG-V2
    # muon cross-variant, which no resample target owns; rerun it by hand with
    #   powheg/reweight/run_reweight_members.sh --current mu --rwl rwl_nuclear 1000
    if echo "$ES" | grep -qw 1000; then
      echo "-- 1000 GeV: nuclear PDF weights and bands"
      run env FORCE=1 ./powheg/reweight/run_reweight_members.sh \
          --current nu --rwl rwl_nuclear 1000 6
      run "$BENCH_PYTHON" analysis/npdf_impact.py --current nu
    fi

    for E in $ES; do
      for S in inclusive faser_s faser_dimuon; do
        echo "-- $E GeV $S: selection dependence"
        run env BENCH_ENERGY="$E" BENCH_SELECTION="$S" \
            "$BENCH_PYTHON" analysis/powheg_selection_pdf.py || true
      done
    done
    # THE SUBSET IS CHOSEN FROM THE SELECTION, so it must be re-chosen after
    # the twin changes -- its old indices point at events that are no longer
    # the selected ones, which showed up as 5 events instead of ~120.
    for E in $ES; do
      echo "-- $E GeV: dimuon subset indices and reweighting"
      run env BENCH_ENERGY="$E" BENCH_SELECTION=faser_dimuon \
          "$BENCH_PYTHON" analysis/dump_selected_indices.py
      run env FORCE=1 ./powheg/reweight/run_reweight_selected.sh \
          "$E" faser_dimuon 6
      run env BENCH_ENERGY="$E" BENCH_SELECTION=faser_dimuon \
          "$BENCH_PYTHON" analysis/powheg_selection_members.py || true
    done
  fi
  ;;
powheg)
  for E in $ES; do
    echo "-- $E GeV: shower, then every selection"
    run ./powheg/run_powheg_shower.sh "$E"
    for S in inclusive faser_s faser_e faser_dimuon ally; do
      for K in powheg powheg_charmfinal; do
        run env BENCH_ENERGY="$E" BENCH_SELECTION="$S" \
            "$BENCH_PYTHON" analysis/analyze.py "$K" || true
      done
    done
  done
  ;;
mg5_ps)
  # MG5_aMC LO + Pythia 8, both currents.  The sample itself is remade by
  # mg5/run_mg5.sh; everything below is what reads it.  ANCHOR ENERGY ONLY --
  # the MG5 arm exists at 1 TeV, and passing another energy here would analyse
  # directories that do not exist and write empty results.
  echo "-- parton level from the same Les Houches files"
  run "$BENCH_PYTHON" analysis/mg5_lhe_histos.py mg5/PROC_MU_LO mg5_me
  run "$BENCH_PYTHON" analysis/mg5_lhe_histos.py mg5/PROC_NU_LO mg5_me --nu
  for S in inclusive faser_s faser_e faser_dimuon ally; do
    for K in mg5_ps mg5_ps_charmfinal; do
      echo "-- 1000 GeV $S: $K (mu)"
      run env BENCH_SELECTION="$S" "$BENCH_PYTHON" analysis/analyze.py "$K" || true
      echo "-- 1000 GeV $S: $K (nu)"
      run env BENCH_SELECTION="$S" "$BENCH_PYTHON" analysis/analyze_nu.py "$K" || true
    done
  done
  ;;
genie_nu)
  for E in $ES; do
    for S in inclusive faser_s faser_e faser_dimuon ally; do
      for K in genie genie_charmfinal genie_lo genie_lo_charmfinal; do
        run env BENCH_ENERGY="$E" BENCH_SELECTION="$S" \
            "$BENCH_PYTHON" analysis/analyze_nu.py "$K" || true
      done
    done
  done
  ;;
*) echo "unknown sample '$SAMPLE' -- try --list" >&2; exit 1 ;;
esac

echo
echo "== downstream re-derivation done.  Now: =="
echo "   tools/check_manifests.py"
echo "   tools/make_all_plots.sh && analysis/make_report.py"
echo "   tools/check_paper_claims.py"
echo "   tools/export_faser_emulsion_csv.py --check"
echo "   tools/export_sidis_share.py --check"
