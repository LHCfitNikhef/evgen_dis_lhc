# Shared helpers for the POWHEG part of the paper-plots production.
#
# Sourced (bash 3.2) by run_res.sh, run_v2.sh, validate.sh and production.sh,
# AFTER config.sh.  The spec is PRODUCTION.md; this file is the one place
# that turns (code, energy, target) into its names, so the drivers and the
# analysis branches cannot drift apart.
#
# NAMES (PRODUCTION.md, "Layout"):
#   POWHEG-RES mu NC   run dir  $POWHEG_RES/v2-mu<TAG>-<t>
#                      samples  powheg/v2mu_<t>_job[_TAG]_<seed>
#   POWHEG-V2  nu CC   run dir  $POWHEG_V2/v2-nu<TAG>-<t>        (batch 1)
#                               $POWHEG_V2/v2-nu<TAG>-<t>-b<j>   (batch j>1)
#                      samples  powheg/v2nu_<t>_job[_TAG]_<j>
# <TAG> is ALWAYS present in the run directories (1TeV included, as in the earlier production);
# the sample directories use beams.at_energy, i.e. untagged at 1 TeV.
#
# The earlier names (parallel-mu*-wide, prod-nu*, job_*, nu_job_*, v2mu_job_*) are
# never produced here.  NB `v2mu_job_1` already exists in powheg/ -- it is the
# earlier POWHEG-V2 muon NC cross-variant -- and differs from every name below by
# the explicit target field.

# EXTRA LHAPDF DIRECTORY, prepended AFTER config.sh.  config.sh assigns
# LHAPDF_DATA_PATH unconditionally, so a caller's `LHAPDF_DATA_PATH=... driver`
# is silently overwritten -- validate.sh's first proton-beam + neutron-PDF run
# died on exactly that: id 339900 resolved through conda's index to
# NNPDF40_nlo_nf_4_pdfas, not installed, and pwhg_main aborted.  (Had that
# id resolved to an INSTALLED set, it would have run on the wrong PDF without
# a word -- memory: lhapdf-index-collision.)  The production never sets this:
# its neutron uses the proton set, which config.sh's path already holds.
if [ -n "${V2_LHAPDF_PREPEND:-}" ]; then
    LHAPDF_DATA_PATH=$V2_LHAPDF_PREPEND:$LHAPDF_DATA_PATH
    export LHAPDF_DATA_PATH
fi

V2_HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
V2_POWHEG_DIR=$(cd "$V2_HERE/.." && pwd)

# The floor this production refuses to cross (CONVENTIONS.md rule 1b says 150 GB;
# the task allows no less than 155 GB left after a batch).
V2_DISK_FLOOR_GB=${V2_DISK_FLOOR_GB:-155}

v2_tag () {        # $1 = current, $2 = energy  ->  1TeV / 400GeV
    python3 "$BENCH_REPO/analysis/beams.py" tag "$1" "$2"
}

v2_jobbase () {    # $1 = mu|nu, $2 = p|n, $3 = energy  ->  v2mu_p_job[_TAG]
    python3 "$BENCH_REPO/analysis/beams.py" name "v2$1_$2_job" "$3"
}

v2_target_ih2 () { # POWHEG's hadron-2 code: 1 proton, 2 neutron
    case "$1" in
        p) echo 1 ;;
        n) echo 2 ;;
        *) echo "target must be p or n, not '$1'" >&2; return 1 ;;
    esac
}

v2_target_pdg () {
    case "$1" in p) echo 2212 ;; n) echo 2112 ;; esac
}

# Seed bases.  Every (code-independent) energy x target gets its own block of
# 1000 seeds, so the proton and neutron runs -- and the five energies -- draw
# DISJOINT random streams.  Sharing a stream would correlate the p and n
# integrator errors, which then ADD in the tungsten combination instead of
# averaging, and the propagated error would be understated.
v2_seed_base () {  # $1 = p|n, $2 = energy
    local ti ei
    case "$1" in p) ti=1 ;; n) ti=2 ;; esac
    case "$2" in
        400) ei=0 ;; 700) ei=1 ;; 1000) ei=2 ;; 2000) ei=3 ;; 4000) ei=4 ;;
        300) ei=5 ;;   # the SIDIS-ladder point (beams.SIDIS_ENERGIES), its own block
        20) ei=6 ;; 50) ei=7 ;; 100) ei=8 ;; 200) ei=9 ;;   # POWHEG-V2 only (beams.SIDIS_ENERGIES_LOW)
        *) echo "energy $2 is not on the energy ladder (beams.SIDIS_ENERGIES)" >&2
           return 1 ;;
    esac
    echo $(( 1000 * (10 * ti + ei) ))
}

# The statistics targets (PRODUCTION.md): showered events per target.
v2_target_events () {
    if [ "$1" = 1000 ]; then echo 1000000; else echo 200000; fi
}

# Free space on the volume holding $PHYSICS24 (the samples), in GB.
v2_free_gb () {
    df -Pk "${PHYSICS24:-$HOME}" | awk 'NR==2 {printf "%d", $4 / 1048576}'
}

# Refuse to start work that would leave less than the floor.
v2_disk_check () {   # $1 = GB the next batch may write, $2 = label
    local free need
    free=$(v2_free_gb)
    need=$(( V2_DISK_FLOOR_GB + $1 ))
    if [ "$free" -lt "$need" ]; then
        echo "DISK: $free GB free, batch '$2' may write $1 GB and the floor" \
             "is $V2_DISK_FLOOR_GB GB -- STOPPING (CONVENTIONS.md rule 1b)" >&2
        return 1
    fi
    echo "  disk: $free GB free (batch '$2' <= $1 GB, floor $V2_DISK_FLOOR_GB GB)"
}

# Wait until fewer than $1 background jobs of THIS shell are running.
v2_throttle () {
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$1" ]; do sleep 2; done
}

# Is an LHE file complete (closing tag present)?
v2_lhe_complete () {
    [ -s "$1" ] && tail -c 200 "$1" | grep -q "</LesHouchesEvents>"
}

# The two beam ids on the <init> line of an LHE file: "13 2112".
v2_lhe_beams () {
    awk '/<init>/ {getline; print $1, $2; exit}' "$1"
}

# A shower job is complete when main_powheg has written its JSON sidecar
# (only after closing the HepMC writer) AND v2_shower has verified it and
# left V2_OK.  The marker matters: a job that wrote its JSON but failed the
# beam or event-count check must not be counted as delivered.
v2_shower_done () {
    [ -s "$1/events.hepmc" ] && [ -s "$1/events_xsec.json" ] && [ -f "$1/V2_OK" ]
}

# Showered events over the COMPLETE job directories <root>/<base>_<N>, N an
# exact integer.
#
# >>> "<base>_*" IS NOT ENOUGH. <<<  At 1 TeV the base is untagged
# (beams.at_energy), so v2nu_p_job_* also matches v2nu_p_job_400GeV_1, whose
# trailing field is an integer too.  The first production counted the 400 GeV
# sample into the 1 TeV total ("have 205450/1000000" before a single 1 TeV
# job existed) and would have stopped ~200k events short, silently.  The
# WHOLE remainder after "<base>_" must be digits -- the rule analyze.py's
# job_files() applies.
v2_delivered () {   # $1 = job root, $2 = job base
    local tot=0 d n rest
    for d in "$1/$2"_*; do
        rest=${d#"$1/$2"_}
        case "$rest" in ''|*[!0-9]*) continue ;; esac
        v2_shower_done "$d" || continue
        n=$(sed -E 's/.*"n_accepted": ([0-9]+).*/\1/' "$d/events_xsec.json")
        tot=$((tot + n))
    done
    echo "$tot"
}

v2_count_events () {   # events in a HepMC file
    local n
    n=$(grep -c '^E ' "$1" 2>/dev/null); echo "${n:-0}"
}

# Shower one LHE into one job directory and verify it.  Runs in the
# foreground; callers background it.
#   v2_shower <lhe> <jobdir> <cmnd> <target p|n> <lepton pdg>
v2_shower () {
    local lhe=$1 d=$2 card=$3 tgt=$4 lep=$5 want got n nlhe
    want="$lep $(v2_target_pdg "$tgt")"
    got=$(v2_lhe_beams "$lhe")
    if [ "$got" != "$want" ]; then
        echo "  REFUSING to shower $lhe: LHE <init> beams '$got', want '$want'" >&2
        return 1
    fi
    mkdir -p "$d"
    ( cd "$V2_POWHEG_DIR" && bench_run ./main_powheg "$card" "$lhe" \
        "$d/events" > "$d/shower.log" 2>&1 ) || {
        echo "  main_powheg failed for $d (see shower.log)" >&2; return 1; }
    [ -s "$d/events_xsec.json" ] || { echo "  no events_xsec.json in $d" >&2; return 1; }
    n=$(v2_count_events "$d/events.hepmc")
    nlhe=$(cat "$lhe.nevents" 2>/dev/null) || {
        echo "  no $lhe.nevents -- strip_lhe_nan.py was not run" >&2; return 1; }
    got=$(python3 "$V2_HERE/inspect_sample.py" beam "$d")
    if [ "$got" != "$want" ]; then
        echo "  BAD BEAMS in $d/events.hepmc: '$got', want '$want'" >&2
        return 1
    fi
    if [ "$n" -lt $(( nlhe * 9 / 10 )) ]; then
        echo "  SHORT: $d delivered $n of $nlhe offered LHE events" >&2
        return 1
    fi
    echo "$n events of $nlhe offered, beams $got, $(date)" > "$d/V2_OK"
    echo "  shower ok: $(basename "$d") $n/$nlhe events, beams $got"
}

# The Pythia steering is the earlier card, unchanged: it sets no beams (they come
# from the LHE) and PDF:pSet MUST stay the PROTON set for a neutron beam too,
# because Pythia 8.311 isospin-conjugates a 2112 beam's PDF itself
# (PDF::xf, beamType -1).  A neutron set there would swap TWICE and shower a
# neutron-labelled beam with proton PDFs.
#
# >>> THE NEUTRINO CARD IS powheg_nu_v2.cmnd (2026-09-14): the earlier card plus
# LesHouches:matchInOut = off, which stops Pythia dropping ~2% of POWHEG-V2
# events, nearly all charm (see the card).  The muon (POWHEG-RES) card is the earlier production's:
# its losses are different ("parton+hadronLevel failed") and the switch does
# not cure them.
v2_cmnd () {       # $1 = mu|nu
    local c="$V2_POWHEG_DIR/powheg_${1}1TeV.cmnd"
    [ "$1" = nu ] && c="$V2_POWHEG_DIR/powheg_nu_v2.cmnd"
    [ "$1" = nu ] && { grep -q '^LesHouches:matchInOut = off$' "$c" || {
        echo "$c no longer sets LesHouches:matchInOut = off -- refusing" >&2
        return 1; }; }
    grep -q '^PDF:pSet = LHAPDF6:NNPDF40_nnlo_as_01180$' "$c" || {
        echo "$c no longer sets the proton PDF set -- refusing" >&2; return 1; }
    echo "$c"
}
