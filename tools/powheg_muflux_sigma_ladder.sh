#!/usr/bin/env bash
# sigma(E) for the arXiv:2506.13889 setup at FIXED beam energies, so that the
# flux-mode POWHEG-RES total can be checked against a direct convolution
#     N = integral dx f(x) sigma(x * 7000 GeV)
# of the same card (their alpha_em, Qmin 1.65, no y window) with the same
# flux.  Integrator only (stages 1-3, one seed): no events are needed.
#
# Usage: tools/powheg_muflux_sigma_ladder.sh [E1 E2 ...]   (GeV; default ladder)
# Output: $POWHEG_RES_MUFLUX/ladder-muflux-2506.13889/E<GeV>/pwg-0001-st3-stat.dat
#         and a summary line per energy on stdout; read by
#         analysis/faser_muflux_check.py --ladder.
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"
CARD=$BENCH_REPO/powheg/cards/POWHEG-RES/muflux-2506.13889/powheg.input
PWHG=$POWHEG_RES_MUFLUX/pwhg_main
BASEDIR=$POWHEG_RES_MUFLUX/ladder-muflux-2506.13889
ES=${*:-"30 60 100 200 400 700 1000 2000 4000 6800"}
for E in $ES; do
    R=$BASEDIR/E$E
    if [ -s "$R/pwg-0001-st3-stat.dat" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "E=$E: exists"; continue
    fi
    mkdir -p "$R"
    # the flux card with a fixed beam of energy E: no LEPpdf, fewer calls
    python3 - "$CARD" "$R/powheg.input" "$E" <<'PY'
import re, sys
src, dst, e = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(src).read()
s = re.sub(r"^ebeam1 .*", f"ebeam1 {e}d0        ! fixed beam (sigma ladder)", s, flags=re.M)
s = re.sub(r"^fixed_lepton_beam .*", "fixed_lepton_beam 1", s, flags=re.M)
s = re.sub(r"^LEPpdf .*\n", "", s, flags=re.M)
s = re.sub(r"^ncall1 .*", "ncall1 30000", s, flags=re.M)
s = re.sub(r"^ncall2 .*", "ncall2 30000", s, flags=re.M)
open(dst, "w").write(s)
PY
    echo 1 > "$R/pwgseeds.dat"
    stage () {
        python3 - "$R/powheg.input" "$1" "$2" <<'PY'
import re, sys
p, st, xg = sys.argv[1:4]
s = open(p).read()
s = re.sub(r"^parallelstage\s+\d+", f"parallelstage {st}", s, flags=re.M)
s = re.sub(r"^xgriditeration\s+\d+", f"xgriditeration {xg}", s, flags=re.M)
open(p, "w").write(s)
PY
        ( cd "$R" && echo 1 | bench_run "$PWHG" > "run_st$1_xg$2.log" 2>&1 )
    }
    stage 1 1; stage 1 2; stage 2 1; stage 3 1
    tot=$(grep "grand total total" "$R/pwg-0001-st3-stat.dat" | awk '{print $5, $7}')
    echo "E=$E GeV: sigma (pos-|neg|) = $tot pb"
done
