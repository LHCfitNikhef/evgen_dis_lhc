#!/usr/bin/env bash
# Sherpa MC@NLO ANTINEUTRINO samples for the SIDIS yields (paper plot 15),
# 2026-10-04.  The FASER-ladder run directories already hold the integration
# (Results.zip) at 400 GeV-4 TeV:
#   p: $SHERPA_RUNS/NuDIS_NLO_FASER_nubar_p_E<E>      (2212 beam)
#   n: $SHERPA_RUNS/V2_NuDIS_NLO_FASER_nubar_n_E<E>   (genuine 2112, patched copy)
# 300 GeV is built first with the ladder's own tools (tools/sherpa_faser_ladder.sh
# ONE=nubar:p:300, then tools/faser_emulsion_ladder.sh sherpa-n 300).
# Here each point generates NEV new events with its own seed into `evtsidis`,
# kept for analysis/faser_pions.py (key sherpa_nlo_nubar) and deleted after
# the extraction by tools/nubar_sidis.py.
# Usage: nohup /bin/bash tools/sherpa_nubar_sidis.sh > logs/sherpa_nubar_sidis.log 2>&1 &
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"
export LHAPDF_DATA_PATH=$BENCH_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH
NEV=${NEV:-50000}
MAXJOBS=${MAXJOBS:-4}
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 10; done; }

gen () {   # t E
  local t=$1 E=$2
  local d=$SHERPA_RUNS/NuDIS_NLO_FASER_nubar_p_E$E sh=$SHERPA_INSTALL/bin/Sherpa
  [ "$t" = n ] && { d=$SHERPA_RUNS/V2_NuDIS_NLO_FASER_nubar_n_E$E; sh=$SHERPA_NPATCH/bin/Sherpa; }
  if [ -s "$d/SIDIS_GEN_OK" ]; then echo "skip nubar_${t}_E$E"; return 0; fi
  [ -s "$d/Results.zip" ] || { echo "no integration in $d" >&2; return 1; }
  local toff=0; [ "$t" = n ] && toff=1
  local seed=$(( 7100 + E / 10 + toff ))
  (
    cd "$d" || exit 1
    if [ "$t" = n ]; then
      export SHERPA_LIBRARY_PATH="$SHERPA_NPATCH/lib/SHERPA-MC/" SHERPA_SHARE_PATH="$SHERPA_NPATCH/share/SHERPA-MC/"
    fi
    bench_run "$sh" -e "$NEV" "RANDOM_SEED: $seed" \
        'EVENT_OUTPUT: ["HepMC3_GenEvent[evtsidis]"]' > gen_sidis.log 2>&1 \
      && [ "$(grep -c '^E ' evtsidis)" = "$NEV" ] \
      && grep -q "Time: " gen_sidis.log \
      && date '+%F %T' > SIDIS_GEN_OK \
      && echo "gen done nubar_${t}_E$E" \
      || echo "gen FAILED nubar_${t}_E$E (see $d/gen_sidis.log)" >&2
  )
}

echo "Sherpa nubar SIDIS start $(date)"
# 300 GeV integrations first (proton, then the neutron that reuses its grid)
if [ ! -s "$SHERPA_RUNS/NuDIS_NLO_FASER_nubar_p_E300/Results.zip" ]; then
  ONE=nubar:p:300 /bin/bash "$HERE/sherpa_faser_ladder.sh" > "$BENCH_REPO/logs/sherpa_nubar_p300.log" 2>&1 &
fi
for E in 1000 400 700 2000 4000; do for t in p n; do limit; gen "$t" "$E" & done; done
wait
if [ ! -s "$SHERPA_RUNS/V2_NuDIS_NLO_FASER_nubar_n_E300/Results.zip" ]; then
  /bin/bash "$HERE/faser_emulsion_ladder.sh" sherpa-n 300 > "$BENCH_REPO/logs/sherpa_nubar_n300.log" 2>&1
fi
gen p 300 & gen n 300 & wait
echo "Sherpa nubar SIDIS end $(date)"
