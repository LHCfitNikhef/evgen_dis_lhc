#!/usr/bin/env bash
# GENIE for the neutrino-generator appendix: GENIE beside NuWro and GiBUU, on
# a free proton AND on tungsten with GENIE's full nuclear model.
#
# Usage:
#   genie/run_genie_nugen.sh spline <cfg> <p|W> <emax>
#   genie/run_genie_nugen.sh run    <cfg> <p|W> <energy> <nev_per_job> <jobs>
#
#   cfg  nu_grv    numu CCDISCHARM  G18_02a_00_000   (FASER's default tune)
#        nu_hedis  numu CCHEDIS     GHE19_00a_00_000 (HEDIS, BGR18)
#   the same tunes, lists and overlays as genie/production/genie_job.sh; only
#   the target and the energy differ.
#
# WHY A SEPARATE DRIVER (user, 2026-10-01): the production driver runs one
# FREE nucleon at a time and builds tungsten by isospin.  Here tungsten is the
# NUCLEUS, 1000741840, so GENIE applies its own nuclear model (Fermi motion,
# binding, and for G18_02a the hA2018 intranuclear cascade) -- which is what
# the comparison with NuWro and GiBUU is about.  Energy 200 GeV (user): the
# lowest HEDIS spline knot is 100 GeV, and none of the three generators'
# nuclear models is trustworthy far above it.
#
# >>> PER NUCLEON. <<<  A nuclear spline is the cross-section PER NUCLEUS;
# NuWro and GiBUU quote per nucleon.  events_xsec.json is divided by A here,
# and says so (per_nucleon_of), so a figure cannot mix the two silently.
#
# >>> THE BEAM RECORD IS THE STRUCK NUCLEON. <<<  genie/gtohepmc3 writes the
# hit nucleon, not the nucleus, as the status-4 target on a nuclear run, so x
# and W are per nucleon -- as for NuWro's and GiBUU's converters.
#
# Jobs whose NUGEN_OK exists are skipped, so the driver can be re-invoked.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
G=$HERE
BENCH_REPO=${BENCH_REPO:-$(cd "$G/.." && pwd)}
. "$BENCH_REPO/config.sh" > /dev/null
SPL=$G/splines/nugen
mkdir -p "$SPL"

cfg () {   # sets PID LIST TUNE OVERLAY SEED0
    case "$1" in
        nu_grv)   PID=14; LIST=CCDISCHARM; TUNE=G18_02a_00_000;   OVERLAY=$G/config/nu_cc_lo:$G/config/p8; SEED0=68000 ;;
        nu_hedis) PID=14; LIST=CCHEDIS;    TUNE=GHE19_00a_00_000; OVERLAY=$G/config/p8;                SEED0=78000 ;;
        *) echo "unknown cfg '$1' (nu_grv | nu_hedis)" >&2; exit 1 ;;
    esac
    [ "$1" = nu_grv ] && python3 "$G/make_nucc_lo_overlay.py" > /dev/null
    local d
    for d in $(echo "$OVERLAY" | tr ':' ' '); do
        # GXMLPATH FAILS OPEN: a missing overlay silently runs the built-in tune
        [ -d "$d" ] || { echo "missing overlay $d -- refusing" >&2; exit 1; }
    done
}
tcode () { case "$1" in p) echo 1000010010 ;; W) echo 1000741840 ;; *) echo "target p|W" >&2; exit 1 ;; esac; }
anuc () { case "$1" in p) echo 1 ;; W) echo 184 ;; esac; }
toff () { case "$1" in p) echo 0 ;; W) echo 700 ;; esac; }
spline_file () { echo "$SPL/$1_$2_e$3.xml"; }
last_knot () { grep -o '<E> *[0-9.]*' "$1" | tail -1 | awk '{print $2}'; }

do_spline () {
    local c=$1 t=$2 emax=$3 f k
    cfg "$c"; f=$(spline_file "$c" "$t" "$emax")
    if [ -s "$f" ]; then echo "exists: $f (last knot $(last_knot "$f"))"; return 0; fi
    echo "building $f ($LIST, pid $PID on $t, $TUNE, to $emax GeV)"
    ( cd "$SPL" && bench_run "$BENCH_SHELL" -c "
        source $G/genie_setup.sh > /dev/null
        export GXMLPATH=$OVERLAY
        gmkspl -p $PID -t $(tcode "$t") -e $emax --tune $TUNE \
               --event-generator-list $LIST -o $f.part" \
      > "$SPL/gmkspl_${c}_${t}_e$emax.log" 2>&1 )
    [ -s "$f.part" ] || { echo "gmkspl produced nothing: $SPL/gmkspl_${c}_${t}_e$emax.log" >&2; exit 1; }
    k=$(last_knot "$f.part")
    # evaluating a spline AT its last knot returns 0 and gevgen hangs, and a
    # clamped spline keeps the requested name with a shorter range
    awk -v k="$k" -v e="$emax" 'BEGIN{exit !(k+0 >= e-1e-6)}' || {
        echo "spline $f.part ends at $k GeV, asked for $emax" >&2; exit 1; }
    mv "$f.part" "$f"
    echo "spline done: $f (last knot $k)"
}

do_run () {
    local c=$1 t=$2 e=$3 nev=$4 jobs=$5 f base i d emax=$6
    cfg "$c"
    f=$(spline_file "$c" "$t" "$emax")
    [ -s "$f" ] || { echo "missing spline $f -- run: $0 spline $c $t $emax" >&2; exit 1; }
    awk -v e="$e" -v k="$(last_knot "$f")" 'BEGIN{exit !(e+0 < k+0)}' || {
        echo "beam $e GeV is not below the last knot of $f" >&2; exit 1; }
    base=$(python3 "$BENCH_REPO/analysis/beams.py" name "nugen_${c}_${t}_job" "$e")
    echo "GENIE $c on $t at $e GeV: $jobs x $nev -> genie/${base}_N"
    for i in $(seq 1 "$jobs"); do
        d=$G/${base}_$i
        [ -f "$d/NUGEN_OK" ] && { echo "  skip ${base}_$i (complete)"; continue; }
        mkdir -p "$d"
        ( cd "$d" && bench_run "$BENCH_SHELL" -c "
            source $G/genie_setup.sh > /dev/null
            export GXMLPATH=$OVERLAY
            gevgen -n $nev -p $PID -t $(tcode "$t") -e $e \
              --tune $TUNE --event-generator-list $LIST \
              --cross-sections $f --seed $((SEED0 + ${e%.*} + $(toff "$t") + i)) \
              -o events.ghep.root > genie.log 2>&1 \
            && $G/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
            && python3 $G/spline_to_json.py $f $e $nev > events_xsec_nucleus.json" \
          && verify "$d" "$nev" "$c" "$t" ) &
    done
    wait
}

verify () {   # $1 dir, $2 nev, $3 cfg, $4 target
    local d=$1 nev=$2 n a
    n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null || true)
    if [ "${n:-0}" -ne "$nev" ]; then
        echo "  SHORT: $d has ${n:-0} of $nev events -- kept for inspection" >&2
        return 1
    fi
    grep -q "Custom directory" "$d/genie.log" || {
        echo "  $d/genie.log does not report a custom config directory" >&2; return 1; }
    a=$(anuc "$4")
    python3 - "$d" "$a" <<'EOF'
import json, sys
d, a = sys.argv[1], int(sys.argv[2])
x = json.load(open(f"{d}/events_xsec_nucleus.json"))
x["sigma_nucleus_mb"] = x["sigma_gen_mb"]
x["sigma_gen_mb"] = x["sigma_gen_mb"] / a
x["per_nucleon_of"] = a
json.dump(x, open(f"{d}/events_xsec.json", "w"))
EOF
    python3 "$G/nondis_log_digest.py" "$d" > /dev/null
    rm -f "$d/events.ghep.root"
    echo "$n events, cfg $3, target $4, $(date)" > "$d/NUGEN_OK"
    echo "  ok $(basename "$d"): $n events"
}

case "$1" in
    spline) do_spline "$2" "$3" "$4" ;;
    run)    do_run "$2" "$3" "$4" "$5" "$6" "${7:-400}" ;;
    *) sed -n '2,12p' "$0"; exit 1 ;;
esac
