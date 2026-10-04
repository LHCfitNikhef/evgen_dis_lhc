#!/usr/bin/env bash
# The single-inclusive pion/kaon extraction over the samples (analysis/
# faser_pions.py): both currents (CONVENTIONS.md rule 2b), GENIE always (rule 1b),
# each nucleon separately at every energy of the SIDIS ladder
# (beams.SIDIS_ENERGIES), then the tungsten combination per energy.
# Each extraction is one single-threaded parse of one sample and writes one
# JSON, so they run detached and in parallel; a pass whose JSON already
# exists is skipped unless FORCE=1.
#
# The combination needs the region's tracked cross-section on each nucleon,
# histos_<key>_q4w3_{p,n}[_TAG].json -- at 300 GeV that is
#     tools/analyse_production.py --energies 300
# on the 300 GeV samples, run once before this.
#
# Usage: tools/faser_pions_passes.sh [extract|combine|all]     (MAXJOBS=6)
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh" > /dev/null
MAXJOBS=${MAXJOBS:-6}
mkdir -p "$BENCH_REPO/logs"
limit () { while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 5; done; }
NU_KEYS=${NU_KEYS:-"powheg_nu genie_lo genie herwig_nlo_full sherpa_nlo"}
MU_KEYS=${MU_KEYS:-"powheg genie herwig_nlo_powheg_full sherpa"}
ES=${ES:-$(python3 -c "import sys; sys.path.insert(0, '$BENCH_REPO/analysis'); import beams; print(' '.join('%g' % e for e in beams.SIDIS_ENERGIES))")}
WHAT=${1:-all}
FP=$BENCH_REPO/analysis/faser_pions.py

out_of () {   # $1 cur $2 key $3 E $4 target -> the JSON the pass writes
    python3 - "$BENCH_REPO" "$@" <<'EOF'
import sys
sys.path.insert(0, sys.argv[1] + "/analysis")
import faser_pions as fp
print(fp.spectra_path(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]))
EOF
}

tasks () {
    for k in $NU_KEYS; do for E in $ES; do echo "nu $k $E"; done; done
    for k in $MU_KEYS; do for E in $ES; do echo "mu $k $E"; done; done
}

if [ "$WHAT" = extract ] || [ "$WHAT" = all ]; then
    while read -r cur key E; do
        for t in p n; do
            # GENIE's default tune is declared valid to 1 TeV (genie/production/genie_job.sh)
            case "$cur $key" in "nu genie_lo"|"mu genie") [ "${E%.*}" -le 1000 ] || continue ;; esac
            o=$(out_of "$cur" "$key" "$E" "$t")
            if [ -s "$o" ] && [ "${FORCE:-0}" != 1 ]; then continue; fi
            limit
            log=$BENCH_REPO/logs/pions_${cur}_${key}_${E}_${t}.log
            ( python3 "$FP" extract --current "$cur" --key "$key" --energy "$E" --target "$t" \
                  > "$log" 2>&1 \
              && echo "done: $cur $key $E $t  $(tail -1 "$log")" \
              || echo "FAILED: $cur $key $E $t  $(tail -1 "$log")" ) &
        done
    done < <(tasks)
    wait
    echo "extraction passes finished"
fi

if [ "$WHAT" = combine ] || [ "$WHAT" = all ]; then
    while read -r cur key E; do
        case "$cur $key" in "nu genie_lo"|"mu genie") [ "${E%.*}" -le 1000 ] || continue ;; esac
        p=$(out_of "$cur" "$key" "$E" p); n=$(out_of "$cur" "$key" "$E" n)
        [ -s "$p" ] && [ -s "$n" ] || { echo "SKIP combine $cur $key $E: a nucleon pass is missing"; continue; }
        python3 "$FP" combine --current "$cur" --key "$key" --energy "$E" \
            > "$BENCH_REPO/logs/pions_${cur}_${key}_${E}_W.log" 2>&1 \
          && echo "combined: $cur $key $E" \
          || echo "FAILED combine: $cur $key $E  $(tail -1 "$BENCH_REPO/logs/pions_${cur}_${key}_${E}_W.log")"
    done < <(tasks)
    echo "combinations finished"
fi
