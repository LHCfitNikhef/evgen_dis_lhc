# Shared functions for the Herwig 7 production (PRODUCTION.md).
# SOURCED, never executed:. "$HERE/lib.sh"
#
# bash 3.2 (macOS /bin/bash) and Linux: no associative arrays, no wait -n.
#
# Layout
#   herwig7/production/cards/<STEM>.in            written by make_cards.py (tracked)
#   herwig7/production/build/<STEM>/              one `Herwig read` per stem: the .run,
#                                         its Herwig-cache, the plugin copy,
#                                         read.log and fingerprint.log (ignored)
#   herwig7/v2<cur>pwg[neg]_<t>_job[_TAG]_N/   one job (ignored, *job_*/)
#
# THE .run CARRIES THE CUT OBJECT SERIALISED FIELD BY FIELD, so it is stale
# when EITHER the card OR LeptonicDISCut.so changes -- run_herwig_pwg.sh:58
# tests only the card.  Here the fingerprint is md5(card + plugin), stored
# with the .run AND in every job directory; a completed job whose fingerprint
# differs from the current build is refused, never silently mixed or
# overwritten.

V2_REPO=${BENCH_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
. "$V2_REPO/config.sh"
V2_HDIR=$V2_REPO/herwig7
V2_DIR=$V2_HDIR/production   # was .../v2 until the 09-30 rename missed it (fixed 2026-10-04)
V2_CARDS=$V2_DIR/cards
V2_BUILD=$V2_DIR/build
V2_PLUGIN=$V2_HDIR/plugin/LeptonicDISCut.so
V2_ENERGIES="400 700 1000 2000 4000 300"  # beams.SIDIS_ENERGIES: 300 GeV is the
                                          # pion-study point, appended so the seed
                                          # blocks of the five energies do not move
# THE NEUTRON SET LIVES IN THIS REPOSITORY (tools/make_neutron_pdf.py), not on
# config.sh's LHAPDF_DATA_PATH.  Prepended, never replacing.  Production cards
# use only the proton set; the validation arms need the neutron one.
V2_LHAPDF_PATH="$V2_REPO/data/pdfs/lhapdf:$LHAPDF_DATA_PATH"
# Rule 1b is 150 GB; the brief for this production asks for 155.
V2_DISK_FLOOR_GB=155

v2_md5() {
    if command -v md5sum >/dev/null 2>&1; then md5sum | cut -d' ' -f1
    else md5 -q; fi
}

hw_env() {
    # LD_LIBRARY_PATH passes through: macOS ignores it, and on Linux an LCG
    # compiler's libstdc++ is found only through it ("GLIBCXX not found")
    bench_run env -i HOME="$HOME" LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" \
        PATH="$HERWIG_BIN:$BENCH_EXTRA_BIN:/usr/bin:/bin" \
        LHAPDF_DATA_PATH="$V2_LHAPDF_PATH" "$@" < /dev/null
}

v2_tag() {   # v2_tag <E>  -> "" at the anchor, "_400GeV" etc. (beams.at_energy)
    if [ "$1" = 1000 ]; then echo ""
    else echo "_$(python3 "$V2_REPO/analysis/beams.py" tag mu "$1")"; fi
}

v2_stem() {  # v2_stem <cur> <pos|neg> <t> <E>   (make_cards.stem)
    local h=PWG; [ "$2" = neg ] && h=PWGNEG
    local t; t=$(v2_tag "$4")
    echo "V2-$1-$h-$3${t:+-${t#_}}"
}

v2_jobdir() {  # v2_jobdir <cur> <pos|neg> <t> <E> <N>
    local h=pwg; [ "$2" = neg ] && h=pwgneg
    echo "$V2_HDIR/v2$1${h}_$3_job$(v2_tag "$4")_$5"
}

# Distinct per current, target, half, energy and job:
#   2 000 000 + 400 000 [nu] + 200 000 [n] + 100 000 [neg] + 1000*e_index + N
v2_seed() {  # v2_seed <cur> <pos|neg> <t> <E> <N>
    local ci=0 ti=0 hi=0 ei=0 e
    [ "$1" = nu ] && ci=1
    [ "$3" = n ] && ti=1
    [ "$2" = neg ] && hi=1
    for e in $V2_ENERGIES; do [ "$e" = "$4" ] && break; ei=$((ei + 1)); done
    [ "$ei" -le 5 ] || { echo "v2_seed: energy $4 not in $V2_ENERGIES" >&2; return 1; }
    [ "$5" -ge 1 ] && [ "$5" -le 999 ] || { echo "v2_seed: job $5 out of range" >&2; return 1; }
    echo $((2000000 + ci*400000 + ti*200000 + hi*100000 + ei*1000 + $5))
}

v2_fingerprint() {  # v2_fingerprint <stem>
    cat "$V2_CARDS/$1.in" "$V2_PLUGIN" | v2_md5
}

# Compile one stem into build/<stem>/ if missing or stale.  Safe to call from
# several background jobs for DIFFERENT stems (each has its own directory).
v2_compile() {  # v2_compile <stem>
    local stem=$1 b=$V2_BUILD/$1 fp
    [ -f "$V2_CARDS/$stem.in" ] || {
        echo "no card $V2_CARDS/$stem.in -- run herwig7/production/make_cards.py" >&2
        return 1; }
    fp=$(v2_fingerprint "$stem")
    if [ -f "$b/$stem.run" ] && [ "$(cat "$b/fingerprint.log" 2>/dev/null)" = "$fp" ]; then
        return 0
    fi
    echo "  compiling $stem"
    rm -rf "$b"; mkdir -p "$b"
    cp "$V2_CARDS/$stem.in" "$V2_PLUGIN" "$b/"
    (cd "$b" && hw_env Herwig read "$stem.in" > read.log 2>&1) || {
        echo "Herwig read FAILED for $stem; see $b/read.log" >&2; tail -5 "$b/read.log" >&2
        return 1; }
    [ -f "$b/$stem.run" ] || { echo "no $stem.run after read -- check saverun" >&2; return 1; }
    echo "$fp" > "$b/fingerprint.log"
}

# A job is COMPLETE when the end-of-run table reports exactly NEV generated
# events, events.hepmc holds NEV events and its closing tag, and the stored
# fingerprint matches.  Return: 0 complete, 1 absent/incomplete, 2 complete
# but built from a DIFFERENT card or plugin (refuse).
v2_job_state() {  # v2_job_state <dir> <stem> <seed> <nev>
    local d=$1 out="$1/$2-S$3.out" n
    [ -s "$d/events.hepmc" ] && [ -s "$out" ] || return 1
    grep -Eq "Total \(from generated events\):[[:space:]]+$4[[:space:]]" "$out" || return 1
    tail -c 200 "$d/events.hepmc" | grep -q "END_EVENT_LISTING" || return 1
    n=$(grep -c '^E ' "$d/events.hepmc")
    [ "$n" = "$4" ] || return 1
    [ "$(cat "$d/fingerprint.log" 2>/dev/null)" = "$(v2_fingerprint "$2")" ] || return 2
    return 0
}

# df -P: POSIX output, one line per filesystem (Linux wraps long device names)
v2_avail_gb() { df -Pk "${PHYSICS24:-$HOME}" | awk 'NR==2 {printf "%d", $4/1048576}'; }

# Run one job to completion (foreground).  Compiles on demand.
v2_run_job() {  # v2_run_job <cur> <pos|neg> <t> <E> <N> <nev>
    local stem d seed
    stem=$(v2_stem "$1" "$2" "$3" "$4")
    d=$(v2_jobdir "$1" "$2" "$3" "$4" "$5")
    seed=$(v2_seed "$1" "$2" "$3" "$4" "$5") || return 1
    v2_run_in "$stem" "$d" "$seed" "$6"
}

v2_run_in() {  # v2_run_in <stem> <dir> <seed> <nev>
    local stem=$1 d=$2 seed=$3 nev=$4 st
    v2_job_state "$d" "$stem" "$seed" "$nev"; st=$?
    if [ $st = 0 ]; then echo "  skip  $(basename "$d") (complete)"; return 0; fi
    if [ $st = 2 ]; then
        echo "  REFUSE $(basename "$d"): complete but built from a different card/plugin" >&2
        return 1
    fi
    v2_compile "$stem" || return 1
    rm -rf "$d"; mkdir -p "$d"
    cp "$V2_BUILD/$stem/$stem.run" "$V2_BUILD/$stem/LeptonicDISCut.so" \
       "$V2_BUILD/$stem/fingerprint.log" "$d/"
    [ -d "$V2_BUILD/$stem/Herwig-cache" ] && cp -R "$V2_BUILD/$stem/Herwig-cache" "$d/"
    echo "  start $(basename "$d")  seed $seed  $nev events  $(date '+%H:%M:%S')"
    (cd "$d" && hw_env Herwig run "$stem.run" -N "$nev" -s "$seed" -d 0 \
        > herwig_run.log 2>&1)
    v2_job_state "$d" "$stem" "$seed" "$nev"; st=$?
    if [ $st = 0 ]; then
        echo "  done  $(basename "$d")  $(date '+%H:%M:%S')"
    else
        echo "  FAILED $(basename "$d") (state $st); see $d/herwig_run.log" >&2
        return 1
    fi
}
