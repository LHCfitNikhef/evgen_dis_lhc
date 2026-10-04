#!/usr/bin/env bash
# GENIE FOR PAPER PLOTS V2: the five earlier GENIE rows, one nucleon at a time.
#
# Usage:
#   genie/production/genie_job.sh spline <cfg> <p|n> <emax>
#   genie/production/genie_job.sh run    <cfg> <p|n> <energy> <nev_per_job> <jobs>
#   genie/production/genie_job.sh list
#
# <cfg> is one of the earlier configurations, unchanged apart from the target:
#   mu_grv    mu-  EMDIS       G18_02a_00_000  overlay p8           (earlier run_genie.sh)
#   mu_nnpdf  mu-  EMDIS       G18_02a_00_000  overlay nnpdf        (earlier run_genie_nnpdf.sh)
#   nu_grv    numu CCDISCHARM  G18_02a_00_000  nu_cc_lo:p8          (earlier run_genie_nucc_lo.sh)
#   nu_nnpdf  numu CCDISCHARM  G18_02a_00_000  nu_cc_lo:nnpdf       (earlier run_genie_nnpdf_nu.sh)
#   nu_hedis  numu CCHEDIS     GHE19_00a_00_000 p8                  (earlier run_genie_hedis.sh CC)
# and <p|n> is the target nucleon, 1000010010 or 1000000010.  The final region
# (Q2 > 4, W > 3, no y cut) is applied in the analysis: GENIE generates the
# full list phase space, except the muon EM cross-sections, which this build
# floors at Q2 > 4 (patches/genie-em-q2-floor-and-pythia8-teardown.diff).
#
# BOTH NUCLEONS GET FRESH SPLINES, in genie/splines/v2/, rather than the proton
# reusing the earlier production's: p and n then come from the same build and knot grid, and the
# earlier proton spline is a cross-check of the new one (a GENIE rebuild that
# changed the cross-section would show up there, not in a figure).
#
# THE VALIDITY RULE (user, 2026-09-01) is enforced here: G18_02a is declared
# valid to 1000 GeV, so its configurations refuse a beam above that.  HEDIS
# (declared to 1e12) runs at every energy.
#
# THE SPLINE MUST RUN PAST THE BEAM: evaluating a GENIE spline AT its last knot
# returns 0 and gevgen then hangs.  400/700/1000 use e2000, 2000/4000 use
# e8000, and every spline's last knot is checked after the build.
#
# DISK: genie.log runs ~3.7x the HepMC it describes and nothing reads it
# (memory: reclaim-intermediates), and the ghep.root is only the converter's
# input.  After a verified conversion both are removed; the log is reduced to
# genie/nondis_log_digest.py's digest (rejections.json, genie.log.tail) first.
#
# Jobs are skipped when complete (V2_OK), so the driver can be killed and
# re-invoked (CONVENTIONS.md 1c).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
G=$(cd "$HERE/.." && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$G/.." && pwd)}
. "$BENCH_REPO/config.sh" > /dev/null
SPL=$G/splines/v2
mkdir -p "$SPL"

cfg () {   # sets PID LIST TUNE OVERLAY SEED0 EVALID
    # OVERLAY is the FULL GXMLPATH.  config/hienergy (stock CommonParam with
    # GVLD-Emax = 1e4, nothing else) goes in front of p8 for mu_grv ONLY: p8
    # has no CommonParam, so it inherits GENIE's 1000 GeV ceiling and gmkspl
    # CLAMPS an e2000 spline at 1000 (seen 2026-09-14) -- a 1 TeV beam would
    # sit on the last knot, read 0 and hang gevgen.  nnpdf and nu_cc_lo already
    # carry 5000.  It must NEVER go in front of HEDIS: its CommonParam shadows
    # p8/GHE19_00c's and gmkspl dies with "No Common parameters available for
    # Param list HEDIS-SF" (also seen 2026-09-14).
    case "$1" in
        mu_grv)   PID=13; LIST=EMDIS;      TUNE=G18_02a_00_000;   OVERLAY=$G/config/hienergy:$G/config/p8; SEED0=18000; EVALID=1000 ;;
        mu_nnpdf) PID=13; LIST=EMDIS;      TUNE=G18_02a_00_000;   OVERLAY=$G/config/nnpdf;                SEED0=28000; EVALID=1000 ;;
        nu_grv)   PID=14; LIST=CCDISCHARM; TUNE=G18_02a_00_000;   OVERLAY=$G/config/nu_cc_lo:$G/config/p8;    SEED0=38000; EVALID=1000 ;;
        nu_nnpdf) PID=14; LIST=CCDISCHARM; TUNE=G18_02a_00_000;   OVERLAY=$G/config/nu_cc_lo:$G/config/nnpdf; SEED0=48000; EVALID=1000 ;;
        nu_hedis) PID=14; LIST=CCHEDIS;    TUNE=GHE19_00a_00_000; OVERLAY=$G/config/p8;                   SEED0=58000; EVALID=1000000 ;;
        # THE ANTINEUTRINO ROWS (2026-10-04, SIDIS yields with nubar): the
        # same tunes and overlays with the beam turned round, own seed blocks
        nubar_grv)   PID=-14; LIST=CCDISCHARM; TUNE=G18_02a_00_000;   OVERLAY=$G/config/nu_cc_lo:$G/config/p8; SEED0=68000; EVALID=1000 ;;
        nubar_hedis) PID=-14; LIST=CCHEDIS;    TUNE=GHE19_00a_00_000; OVERLAY=$G/config/p8;                SEED0=78000; EVALID=1000000 ;;
        *) echo "unknown cfg '$1'" >&2; exit 1 ;;
    esac
    case "$1" in nu_grv|nu_nnpdf|nubar_grv) python3 "$G/make_nucc_lo_overlay.py" > /dev/null ;; esac
    local d
    for d in $(echo "$OVERLAY" | tr ':' ' '); do
        # GXMLPATH FAILS OPEN: a missing overlay silently runs the built-in tune
        [ -d "$d" ] || { echo "missing overlay $d -- refusing" >&2; exit 1; }
    done
}
# DISJOINT RANDOM STREAMS for the two nucleons: without the offset p and n
# of the same (cfg, energy, job) would share a seed and their statistical
# fluctuations would be correlated in the tungsten combination.
toff () { case "$1" in p) echo 0 ;; n) echo 500 ;; esac; }
tcode () { case "$1" in p) echo 1000010010 ;; n) echo 1000000010 ;; *) echo "target p|n" >&2; exit 1 ;; esac; }
emax_for () { case "$1" in 300|400|700|1000) echo 2000 ;; 2000|4000) echo 8000 ;; *) echo "no spline range for $1 GeV" >&2; exit 1 ;; esac; }
spline_file () { echo "$SPL/$1_$2_e$3.xml"; }
last_knot () { grep -o '<E> *[0-9.]*' "$1" | tail -1 | awk '{print $2}'; }

# GXMLPATH for a build or a run: the overlay as cfg() states it, at every
# energy (the HEDIS e8000 spline needs no raised ceiling: validity 1e12).
gxml () { echo "$1"; }

do_spline () {
    local c=$1 t=$2 emax=$3 f
    cfg "$c"; f=$(spline_file "$c" "$t" "$emax")
    if [ -s "$f" ]; then echo "exists: $f (last knot $(last_knot "$f"))"; return 0; fi
    echo "building $f ($LIST, pid $PID on $t, $TUNE, to $emax GeV)"
    ( cd "$SPL" && bench_run "$BENCH_SHELL" -c "
        source $G/genie_setup.sh > /dev/null
        export GXMLPATH=$(gxml "$OVERLAY" "$emax")
        gmkspl -p $PID -t $(tcode "$t") -e $emax --tune $TUNE \
               --event-generator-list $LIST -o $f.part" \
      > "$SPL/gmkspl_${c}_${t}_e$emax.log" 2>&1 )
    [ -s "$f.part" ] || { echo "gmkspl produced nothing: $SPL/gmkspl_${c}_${t}_e$emax.log" >&2; exit 1; }
    local k; k=$(last_knot "$f.part")
    # a clamped spline carries the requested name and a shorter range
    awk -v k="$k" -v e="$emax" 'BEGIN{exit !(k+0 >= e-1e-6)}' || {
        echo "spline $f.part ends at $k GeV, asked for $emax (validity clamp?)" >&2; exit 1; }
    mv "$f.part" "$f"
    echo "spline done: $f (last knot $k)"
}

do_run () {
    local c=$1 t=$2 e=$3 nev=$4 jobs=$5 emax f base i d
    cfg "$c"
    awk -v e="$e" -v v="$EVALID" 'BEGIN{exit !(e+0 <= v+0)}' || {
        echo "$c ($TUNE) is declared valid to $EVALID GeV; refusing $e GeV" >&2; exit 1; }
    emax=$(emax_for "$e"); f=$(spline_file "$c" "$t" "$emax")
    [ -s "$f" ] || { echo "missing spline $f -- run: $0 spline $c $t $emax" >&2; exit 1; }
    base=$(python3 "$BENCH_REPO/analysis/beams.py" name "v2g_${c}_${t}_job" "$e")
    echo "GENIE $c on $t at $e GeV: $jobs x $nev -> genie/${base}_N (spline e$emax)"
    for i in $(seq 1 "$jobs"); do
        d=$G/${base}_$i
        [ -f "$d/V2_OK" ] && { echo "  skip ${base}_$i (complete)"; continue; }
        mkdir -p "$d"
        ( cd "$d" && bench_run "$BENCH_SHELL" -c "
            source $G/genie_setup.sh > /dev/null
            export GXMLPATH=$(gxml "$OVERLAY" "$emax")
            gevgen -n $nev -p $PID -t $(tcode "$t") -e $e \
              --tune $TUNE --event-generator-list $LIST \
              --cross-sections $f --seed $((SEED0 + ${e%.*} * 10 + $(toff "$t") + i)) \
              -o events.ghep.root > genie.log 2>&1 \
            && $G/gtohepmc3 events.ghep.root events.hepmc >> genie.log 2>&1 \
            && python3 $G/spline_to_json.py $f $e $nev > events_xsec.json" \
          && v2_verify "$d" "$nev" "$c" "$t" ) &
    done
    wait
}

v2_verify () {   # $1 dir, $2 nev, $3 cfg, $4 target
    local d=$1 nev=$2 n
    n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null || true)
    if [ "${n:-0}" -ne "$nev" ]; then
        echo "  SHORT: $d has ${n:-0} of $nev events -- kept for inspection" >&2
        return 1
    fi
    # the overlay really was used: GENIE prints the custom config directory
    grep -q "Custom directory" "$d/genie.log" || {
        echo "  $d/genie.log does not report a custom config directory" >&2; return 1; }
    # rejections.json + genie.log.tail, and the full log deleted
    python3 "$G/nondis_log_digest.py" "$d" > /dev/null
    rm -f "$d/events.ghep.root"
    echo "$n events, cfg $3, target $4, $(date)" > "$d/V2_OK"
    echo "  ok $(basename "$d"): $n events"
}

case "$1" in
    spline) do_spline "$2" "$3" "$4" ;;
    run)    do_run "$2" "$3" "$4" "$5" "$6" ;;
    list)   echo "mu_grv mu_nnpdf nu_grv nu_nnpdf nu_hedis" ;;
    *) sed -n 2,8p "$0"; exit 1 ;;
esac
