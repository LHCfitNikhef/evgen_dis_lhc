#!/usr/bin/env bash
# GENIE event samples for the FASER DIMUON CUT FLOW (user, 2026-09-07):
# the FASER collaboration's own GENIE estimate of the opposite-sign dimuon
# signal at 300 fb^-1, reproduced with our GENIE -- see
# analysis/faser_dimuon_cutflow.py for the cut flow and the physics.
#
# TWO KINDS OF RUN, both with the DEFAULT tune G18_02a_00_000 (= FASER's), the
# full CC event-generator list (every charged-current process, as the slide's
# GENIE would have run) and the benchmark's Pythia8 decay overlay, on a free
# proton and a free neutron, nu_mu and nubar_mu:
#
#   ladder   fixed beam energies, 12 points from 10 GeV to 6.8 TeV, so that
#            the cut flow is convolved with the flux exactly as the POWHEG-V2
#            one is (tools/powheg_v2_dimuon_ladder.sh) -- the common method;
#   flux     gevgen handed the FASERnu flux histogram itself (-f), so its
#            events are distributed as Phi(E) sigma(E) and every step's
#            flux-averaged efficiency is a plain count: the check on the
#            ladder's interpolation, with no ladder anywhere.
#
# Each job is reduced to a compact per-event table (dimuon_events.npz) by
# faser_dimuon_cutflow.py the moment it finishes, and the GHEP and HepMC
# files are deleted: 9.6M ladder events would be 37 GB of event files
# backing a few MB of tables that hold every muon of every event (CONVENTIONS.md
# rule 1).  Regeneration is ~1000 events/s per core.
#
#   emulsion the FASERnu EMULSION ladder of paper plot 12 (fdladder_job_*):
#            10 energies 30 GeV-6.8 TeV, 60k events per point, the same tune,
#            list and splines, reduced by faser_emulsion_shapes.py extract to
#            faserdata_events.npz.  Its first production (2026-09-08) ran
#            outside any tracked driver and recorded no seeds; this stage is
#            how it was regenerated on 2026-10-01 after the charm fix
#            (patches/genie-aivazis-charm-propagator.diff), seeds below.
#
# Usage: tools/genie_dimuon_ladder.sh [ladder|flux|emulsion|all]
#        NEV_LADDER=200000 NEV_FLUX=250000 NJOBS_FLUX=4 MAXJOBS=8 by default
#        NEV_EMULSION=60000 ES_EMULSION="30 60 ... 6800"
# Output: genie/dmuladder_job_<nu|nubar>_<p|n>_E<GeV>_1/dimuon_events.npz
#         genie/dmuflux_job_<nu|nubar>_<p|n>_<i>/dimuon_events.npz
#         genie/fdladder_job_<nu|nubar>_<p|n>_E<GeV>_1/faserdata_events.npz
#         (both match .gitignore's *job_*/ rule; the flux histograms go to
#          genie/faser_flux_hist/ and are gitignored too)
#
# THE SPLINES are those of tools/genie_faser_splines.sh: the CC list to
# 8 TeV, validity cap raised through config/hienergy.  A missing spline HANGS
# gevgen rather than failing it (memory: genie-layout-and-setup), so it is
# checked up front.
#
# THE FLUX HISTOGRAM.  gevgen's text-file flux driver rebuilds the spectrum
# by accept-reject into 300 LINEAR bins from a spline through the points,
# which for a spectrum falling four orders of magnitude over 10 GeV-6 TeV
# leaves the tail as Poisson noise; a ROOT TH1D with the flux file's own
# log-spaced bins and its counts as contents is sampled exactly, and
# ",NOWIDTH" tells gevgen the contents are already counts per bin.  The
# histogram is built here with a ROOT macro (this environment has no pyROOT).
#
# THE LOG.  The Pythia8 decayer prints three lines per decay and the stock
# whisper thresholds do not know its stream; config/p8/Messenger_quiet.xml,
# passed after whisper, silences it.  500 MB of log per point otherwise.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh" > /dev/null
G=$BENCH_REPO/genie
TUNE=G18_02a_00_000
LIST=CC
STAGE=${1:-all}
ES=${ES:-"10 20 30 60 100 200 400 700 1000 2000 4000 6800"}
NEV_LADDER=${NEV_LADDER:-200000}
NEV_FLUX=${NEV_FLUX:-250000}
NJOBS_FLUX=${NJOBS_FLUX:-4}
ES_EMULSION=${ES_EMULSION:-"30 60 100 200 400 700 1000 2000 4000 6800"}
NEV_EMULSION=${NEV_EMULSION:-60000}
# which per-event table a job is reduced to: dimuon (faser_dimuon_cutflow.py)
# or emulsion (faser_emulsion_shapes.py)
TABLE=dimuon
MAXJOBS=${MAXJOBS:-8}
MSG="Messenger_whisper.xml:Messenger_quiet.xml"
FLUXDIR=$G/faser_flux_hist
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }

free_gb () { df -Pk "$PHYSICS24" | awk 'NR==2 {print int($4/1048576)}'; }   # -k: macOS and Linux alike
if [ "$(free_gb)" -lt 158 ]; then
  echo "only $(free_gb) GB free on the $PHYSICS24 volume; the floor is 150 GB" >&2
  echo "and $MAXJOBS parallel jobs hold up to ~7 GB of transient event files" >&2
  exit 1
fi

for pid in 14 -14; do for t in p n; do
  f=$G/splines/faser/${pid}_${t}_cc_e8000.xml
  [ -s "$f" ] || { echo "missing spline $f -- run tools/genie_faser_splines.sh" >&2; exit 1; }
done; done
[ -x "$G/gtohepmc3" ] || { echo "no $G/gtohepmc3 (make -C genie)" >&2; exit 1; }

# The Pythia8 hadroniser prints a particle listing per DIS event straight
# to stdout, outside GENIE's Messenger: ~1 kB per event even with every
# GENIE stream quiet.  Once a job has been extracted the log is reduced to
# its complaints and its tail; a job that FAILS keeps the whole log.
reduce_log () {
  local d=$1
  ( grep -n -E "FATAL|ERROR|CRIT|WARN|Warning|rror" "$d/genie.log" | grep -v "cling::AutoLoading" | head -300
    echo "---- last 200 lines of genie.log ($(wc -l < "$d/genie.log") lines, $(wc -c < "$d/genie.log") bytes) ----"
    tail -n 200 "$d/genie.log" ) > "$d/genie_log_digest.txt" && rm -f "$d/genie.log"
}

# one gevgen job: <jobdir> <pid> <target p|n> <energy-arg> <flux-arg> <seed> <nev> <energy-tag>
run_job () {
  local d=$1
  local pid=$2
  local t=$3
  local earg=$4
  local farg=$5
  local seed=$6
  local nev=$7
  local etag=$8
  local code spline
  case "$t" in p) code=1000010010;; n) code=1000000010;; esac
  spline=$G/splines/faser/${pid}_${t}_cc_e8000.xml
  local script=faser_dimuon_cutflow.py out=dimuon_events.npz
  if [ "$TABLE" = emulsion ]; then script=faser_emulsion_shapes.py; out=faserdata_events.npz; fi
  if [ -s "$d/$out" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $d"; return 0; fi
  mkdir -p "$d"
  limit
  # the GENIE environment is re-sourced INSIDE the wrapper shell: macOS strips
  # DYLD_LIBRARY_PATH across caffeinate (as every run_genie*.sh does)
  ( cd "$d" && bench_run "$BENCH_SHELL" -c "
      source $G/genie_setup.sh > /dev/null
      export GXMLPATH=$G/config/hienergy:$G/config/p8
      gevgen -n $nev -p $pid -t $code -e $earg $farg \
        --tune $TUNE --event-generator-list $LIST \
        --cross-sections $spline --seed $seed \
        --message-thresholds $MSG -o events.ghep.root > genie.log 2>&1 \
      && $G/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
      && rm -f events.ghep.root" \
    && python3 "$BENCH_REPO/analysis/$script" extract \
         --beam="$pid" --target="$t" --energy="$etag" --generator=genie \
         "$d/events.hepmc" "$d/$out" > "$d/extract.log" 2>&1 \
    && [ -s "$d/$out" ] \
    && rm -f "$d/events.hepmc" "$d/input-flux.root" \
    && reduce_log "$d" \
    && echo "done: $(basename "$d")  $(tail -1 "$d/extract.log")" \
    || echo "FAILED: $d (files kept)" ) &
}

do_ladder () {
  for cur in nu:14 nubar:-14; do
    for t in p n; do
      for E in $ES; do
        c=${cur%%:*}; pid=${cur#*:}
        # seed: unique per (beam, nucleon, energy), reproducible
        seed=$(( 7000 + (E % 977) * 4 + $([ "$pid" = 14 ] && echo 0 || echo 2) + $([ "$t" = p ] && echo 0 || echo 1) ))
        run_job "$G/dmuladder_job_${c}_${t}_E${E}_1" "$pid" "$t" "$E" "" "$seed" "$NEV_LADDER" "$E"
      done
    done
  done
  wait
}

make_flux_hist () {   # <pid> -> $FLUXDIR/flux_<pid>.root with TH1D "flux"
  # TWO `local` statements, not one: bash expands every word of a `local`
  # line BEFORE assigning, so `local pid=$1 out=...${pid}...` reads the
  # CALLER'S pid -- which wrote the nu_mu flux into flux_-14.root and then
  # ran the antineutrino jobs on it (2026-09-07, caught by the nu_mu jobs
  # failing to find their own file).
  local pid=$1
  local out=$FLUXDIR/flux_${pid}.root
  [ -s "$out" ] && [ "${FORCE:-0}" != "1" ] && return 0
  mkdir -p "$FLUXDIR"
  # conda python (numpy) writes the macro; ROOT from the GENIE environment runs it
  python3 - "$BENCH_REPO/data/faser_flux/FASER_${pid}.txt" "$FLUXDIR/mk_${pid}.C" "$out" <<'PY'
import sys
import numpy as np
src, macro, out = sys.argv[1:4]
d = np.loadtxt(src)
e, n = d[:, 0], d[:, 1]
r = e[1] / e[0]                       # the file's log spacing
edges = np.concatenate(([e[0] / np.sqrt(r)], np.sqrt(e[:-1] * e[1:]), [e[-1] * np.sqrt(r)]))
body = ("void mk_%s(){ double ed[]={%s}; double c[]={%s}; "
        "TH1D h(\"flux\",\"FASERnu flux, counts per bin\",%d,ed); "
        "for(int i=0;i<%d;i++) h.SetBinContent(i+1,c[i]); "
        "TFile f(\"%s\",\"recreate\"); h.Write(); f.Close(); }") % (
    sys.argv[1].split("FASER_")[1].split(".")[0].replace("-", "m"),
    ",".join("%.6f" % x for x in edges), ",".join("%.6e" % x for x in n),
    len(e), len(e), out)
open(macro, "w").write(body)
print(f"flux histogram {out}: {len(e)} bins, {edges[0]:.2f}-{edges[-1]:.0f} GeV, {n.sum():.4e} neutrinos at 150 fb^-1")
PY
  # ROOT wants the macro's function named after the file, minus the sign
  local mac=$FLUXDIR/mk_${pid}.C
  if [ "$pid" = -14 ]; then mv "$mac" "$FLUXDIR/mk_m14.C"; mac=$FLUXDIR/mk_m14.C; fi
  ( cd "$FLUXDIR" && bench_run "$BENCH_SHELL" -c "
      source $G/genie_setup.sh > /dev/null
      root -l -b -q $mac" > "$FLUXDIR/mk_${pid}.log" 2>&1 )
  [ -s "$out" ] || { echo "flux histogram not written: see $FLUXDIR/mk_${pid}.log" >&2; exit 1; }
}

do_flux () {
  make_flux_hist 14
  make_flux_hist -14
  for cur in nu:14 nubar:-14; do
    for t in p n; do
      c=${cur%%:*}; pid=${cur#*:}
      for i in $(seq 1 $NJOBS_FLUX); do
        seed=$(( 9000 + 10 * i + $([ "$pid" = 14 ] && echo 0 || echo 2) + $([ "$t" = p ] && echo 0 || echo 1) ))
        run_job "$G/dmuflux_job_${c}_${t}_$i" "$pid" "$t" "10,6403" \
                "-f $FLUXDIR/flux_${pid}.root,flux,NOWIDTH" "$seed" "$NEV_FLUX" flux
      done
    done
  done
  wait
}

do_emulsion () {
  TABLE=emulsion
  for cur in nu:14 nubar:-14; do
    for t in p n; do
      for E in $ES_EMULSION; do
        c=${cur%%:*}; pid=${cur#*:}
        # seed: unique per (beam, nucleon, energy), disjoint from the ladder's
        # 7000-10900 and the flux mode's 9000-9050
        seed=$(( 11000 + (E % 977) * 4 + $([ "$pid" = 14 ] && echo 0 || echo 2) + $([ "$t" = p ] && echo 0 || echo 1) ))
        run_job "$G/fdladder_job_${c}_${t}_E${E}_1" "$pid" "$t" "$E" "" "$seed" "$NEV_EMULSION" "$E"
      done
    done
  done
  wait
}

case "$STAGE" in
  ladder)   do_ladder ;;
  flux)     do_flux ;;
  emulsion) do_emulsion ;;
  all)      do_ladder; do_flux; do_emulsion ;;
  *) echo "stage must be ladder, flux, emulsion or all" >&2; exit 2 ;;
esac
echo "stage $STAGE done"
