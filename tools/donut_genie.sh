#!/usr/bin/env bash
# GENIE for the DONUT comparison (user, 2026-09-10): the charged-particle
# multiplicity at the primary vertex of neutrino interactions in DONUT's
# emulsion-and-steel target, arXiv:0711.0728 FIG. 9.
#
# THE ENERGY WEIGHT IS DONUT'S OWN FIG. 2, which is the spectrum of neutrinos
# INTERACTING in the target -- already folded with the cross-section and the
# acceptance.  So gevgen is handed that spectrum as a flux histogram and the
# events come out distributed as the experiment's were, with no ladder and no
# second convolution.  tools/digitise_donut.py extracts it from the paper.
#
# >>> nu_mu AND nubar_mu ONLY, AND THAT IS THE USER'S POINT (2026-09-10): <<<
# "neutrino cross-sections are flavour independent, so if we have interaction
# cross-sections for muon neutrinos, we can also use them for tau neutrinos".
# At the primary vertex the multiplicity is flavour-blind for the same reason:
# every charged-current event contributes exactly one lepton track, whether it
# is an electron, a muon or a tau (the tau decays microns away, and DONUT
# counts the parent).  So the three flavour curves of FIG. 2 are SUMMED into
# one spectrum.  What is not carried over is the tau mass suppression of the
# y distribution, which acts on the 5% of the sample that is nu_tau; said in
# place in the analysis.
#
# THE BEAM IS HALF ANTINEUTRINO.  The paper measures nubar_mu/nu_mu = 1.05 +-
# 0.13 (Section IX B) and takes the fluxes as equal; the two are generated
# separately and combined 50/50, so the assumption is visible rather than
# buried in a single sample.
#
# THE TARGET IS IRON.  DONUT's Emulsion Cloud Chambers interleave 1 mm
# stainless-steel sheets with ~100 um emulsion layers on a plastic base
# (Section IV), so the steel carries about 86% of the mass and Fe-56 is the
# nucleus almost every interaction happens on.  The emulsion's own silver and
# bromine are heavier and would give slightly more final-state interaction;
# that is a stated approximation, not a modelled one.
#
# >>> AND THE EVENT-GENERATOR LIST IS Default, NOT CC. <<<  FIG. 9 is "all the
# located events", which is every channel the beam produces: charged and
# neutral current, quasi-elastic, resonance, DIS and coherent.  Running the CC
# list here -- the list every other GENIE driver in this repository uses --
# would drop a quarter of the sample and shift the multiplicity down by one
# track over that quarter, which looks exactly like a hadronisation effect.
#
# Usage: tools/donut_genie.sh [splines|events|analyse|all]
#        NEV=200000 NJOBS=4 by default
# Output: genie/donut_job_<nu|nubar>_<i>/events.hepmc
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
G=$BENCH_REPO/genie
STAGE=${1:-all}
TUNE=G18_02a_00_000
EMAX=${EMAX:-400}
TGT=1000260560                     # Fe-56
LIST=Default
NEV=${NEV:-200000}
NJOBS=${NJOBS:-4}
SPLDIR=$G/splines/donut
FLUXDIR=$G/donut_flux_hist
MSG=$G/config/p8/Messenger_quiet.xml

mkdir -p "$SPLDIR" "$FLUXDIR"

do_splines () {
  for pid in 14 -14; do
    f=$SPLDIR/${pid}_Fe56_all_e$EMAX.xml
    if [ -s "$f" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $f"; continue; fi
    # the environment is set INSIDE the wrapped shell: macOS strips
    # DYLD_LIBRARY_PATH across the caffeinate wrapper, as every other GENIE
    # driver here does it
    ( cd "$SPLDIR" && bench_run "$BENCH_SHELL" -c "
        source $G/genie_setup.sh > /dev/null
        export GXMLPATH=$G/config/hienergy:$G/config/p8
        gmkspl -p $pid -t $TGT -e $EMAX --tune $TUNE \
               --event-generator-list $LIST -o $f" \
        > "$SPLDIR/gmkspl_${pid}.log" 2>&1 \
      && echo "done: $f" || echo "FAILED: see $SPLDIR/gmkspl_${pid}.log" ) &
  done
  wait
}

make_flux_hist () {
  local out=$FLUXDIR/donut_flux.root
  [ -s "$out" ] && [ "${FORCE:-0}" != "1" ] && return 0
  # ALL THREE FLAVOURS SUMMED, for the reason in the header: at the primary
  # vertex the multiplicity does not know which lepton it made.
  python3 - "$BENCH_REPO/data/donut/flux.json" "$FLUXDIR/mk_donut.C" "$out" <<'PY'
import json, os, sys
src, macro, out = sys.argv[1:4]
d = json.load(open(src))
edges = d["e_edges_gev"]
tot = [sum(d["spectra"][k][i] for k in d["spectra"]) for i in range(len(edges) - 1)]
body = ("void mk_donut(){ double ed[]={%s}; double c[]={%s}; "
        "TH1D h(\"flux\",\"DONUT interacting-neutrino spectrum\",%d,ed); "
        "for(int i=0;i<%d;i++) h.SetBinContent(i+1,c[i]); "
        "TFile f(\"%s\",\"recreate\"); h.Write(); f.Close(); }") % (
    ",".join("%.6f" % x for x in edges), ",".join("%.6e" % x for x in tot),
    len(tot), len(tot), os.path.basename(out))   # ROOT runs in $FLUXDIR
open(macro, "w").write(body)
print(f"flux histogram {out}: {len(tot)} bins, "
      f"{edges[0]:.0f}-{edges[-1]:.0f} GeV, {sum(tot):.1f} events in FIG. 2")
PY
  ( cd "$FLUXDIR" && bench_run "$BENCH_SHELL" -c "
      source $G/genie_setup.sh > /dev/null
      root -l -b -q $FLUXDIR/mk_donut.C" > "$FLUXDIR/mk_donut.log" 2>&1 )
  [ -s "$out" ] || { echo "flux histogram not written: see $FLUXDIR/mk_donut.log" >&2
                     exit 1; }
}

# THE GENIE LOG IS 2.7 GB PER JOB, and it is the Pythia8 decayer talking.
# tools/genie_dimuon_ladder.sh met this first and solved it the same way; not
# reusing that solution here cost 21 GB and took the volume to within 3 GB of
# the 150 GB floor CONVENTIONS.md rule 1b sets.  A digest is kept, the log is not.
reduce_log () {
  local d=$1
  [ -s "$d/genie.log" ] || return 0
  ( grep -n -E "FATAL|ERROR|CRIT|WARN|Warning|rror" "$d/genie.log" \
      | grep -v "cling::AutoLoading" | head -300
    echo "---- last 200 lines of genie.log ($(wc -l < "$d/genie.log") lines, \
$(wc -c < "$d/genie.log") bytes) ----"
    tail -n 200 "$d/genie.log" ) > "$d/genie_log_digest.txt" \
    && rm -f "$d/genie.log"
}

do_events () {
  make_flux_hist
  for cur in nu:14 nubar:-14; do
    c=${cur%%:*}; pid=${cur#*:}
    spline=$SPLDIR/${pid}_Fe56_all_e$EMAX.xml
    # A MISSING SPLINE HANGS gevgen rather than failing it, so it is checked
    # here (memory: genie-layout-and-setup).
    [ -s "$spline" ] || { echo "no spline $spline -- run 'splines' first" >&2
                          exit 1; }
    for i in $(seq 1 "$NJOBS"); do
      d=$G/donut_job_${c}_$i
      if [ -s "$d/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "exists: $d"; continue
      fi
      mkdir -p "$d"
      seed=$(( 4100 + 10 * i + $([ "$pid" = 14 ] && echo 0 || echo 1) ))
      ( cd "$d" && bench_run "$BENCH_SHELL" -c "
          source $G/genie_setup.sh > /dev/null
          export GXMLPATH=$G/config/hienergy:$G/config/p8
          gevgen -n $NEV -p $pid -t $TGT -e 3,300 \
            -f $FLUXDIR/donut_flux.root,flux,NOWIDTH \
            --tune $TUNE --event-generator-list $LIST \
            --cross-sections $spline --seed $seed \
            --message-thresholds $MSG -o events.ghep.root > genie.log 2>&1 \
          && $G/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
          && rm -f events.ghep.root" \
        && reduce_log "$d" \
        && rm -f "$d/input-flux.root" \
        && echo "done: $(basename "$d")" || echo "FAILED: $d" ) &
    done
  done
  wait
}

do_analyse () {
  # THE WHOLE CHAIN IN ONE PLACE (CONVENTIONS.md rule 1): the multiplicity, the
  # merge with the data, the figure and the cross-section comparison.  It was
  # briefly a scratch wrapper, which is exactly the thing that rule exists to
  # stop.
  python3 "$BENCH_REPO/analysis/donut_nch.py" genie "$G"/donut_job_*/events.hepmc
  python3 "$BENCH_REPO/analysis/donut_nch.py" --merge
  python3 "$BENCH_REPO/analysis/plot_donut_nch.py"
  python3 "$BENCH_REPO/analysis/donut_xsec.py"
}

case "$STAGE" in
  splines) do_splines ;;
  events)  do_events ;;
  analyse) do_analyse ;;
  all)     do_splines; do_events; do_analyse ;;
  *) echo "usage: donut_genie.sh [splines|events|analyse|all]" >&2; exit 1 ;;
esac
echo "donut_genie.sh $STAGE done"
