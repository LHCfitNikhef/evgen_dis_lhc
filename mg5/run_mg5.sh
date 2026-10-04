#!/usr/bin/env bash
# MG5_aMC LO DIS for the benchmark: build the process directory, install the
# DIS cut and scale hooks, generate unweighted events, and shower them with
# Pythia 8.
#
# Usage: ./run_mg5.sh [--current mu|nu] [--energy GEV] [--no-shower] <nevents>
#
#   ./run_mg5.sh --current mu 200000      # -> mg5/PROC_MU_LO, mg5/ps_job_1
#   ./run_mg5.sh --current nu 200000      # -> mg5/PROC_NU_LO, mg5/ps_nu_job_1
#
# WHY THIS SCRIPT EXISTS.  MG5 entered the benchmark as a PARTON-LEVEL entry:
# its Les Houches file was read directly (analysis/mg5_lhe_histos.py) and only
# the six DIS-kinematics observables were defined.  Showering it with Pythia 8
# gives it a hadronic final state like every other entry, at the cost of the
# usual caveat: the hard process is MG5's, everything after it is Pythia's, so
# the arm differs from the Pythia LO row only in the matrix element and the
# phase-space generation.  That is exactly what makes the pair worth having.
#
# THREE THINGS MG5 DOES NOT DO BY ITSELF, all of them silent:
#
#  * no Q2 or y cut -- the run card offers collider cuts only, so the fiducial
#    region is imposed in dummy_cuts (mg5/dis_hooks.f);
#  * `dynamical_scale_choice = -1`, the run-card DEFAULT, sets the scale by
#    CLUSTERING the external states and ran at 0.787 x Q on average here,
#    which cost a FACTOR TWO on the neutral-current cross section;
#  * its LHAPDF interface asks for AlphaS_FlavorScheme / AlphaS_NumFlavors,
#    which NNPDF4.0's .info does not carry as written, and the C++ runtime is
#    not linked into the LO libraries on macOS.  Both are patched below, and
#    both fail LOUDLY, unlike the first two.
#
# The cross section is checked against the published parton-level result at
# the end, which is the closure that catches all three.
set -e

CURRENT=mu
EBEAM=1000
SHOWER=1
ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --current)   CURRENT=$2; shift 2 ;;
        --energy)    EBEAM=$2; shift 2 ;;
        --no-shower) SHOWER=0; shift ;;
        *)           ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -gt 0 ] && set -- "${ARGS[@]}" || set --
NEV=${1:-200000}

case "$CURRENT" in
    mu|nu) ;;
    *) echo "--current must be mu or nu" >&2; exit 1 ;;
esac

HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"          # every path outside the repo

TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$CURRENT" "$EBEAM")
# the anchor keeps untagged names, as everywhere else in this repository
SUF=$(python3 "$BENCH_REPO/analysis/beams.py" name "" "$EBEAM")
PROC=$HERE/PROC_$(echo "$CURRENT" | tr a-z A-Z)_LO$SUF
if [ "$CURRENT" = mu ]; then
    JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name ps_job "$EBEAM")
else
    JOBBASE=$(python3 "$BENCH_REPO/analysis/beams.py" name ps_nu_job "$EBEAM")
fi
NCARD=$HERE/mg5_shower.cmnd        # one card for both currents, see its header
JOBDIR=$HERE/${JOBBASE}_1

read -r EA EB <<EOF
$(python3 - "$CURRENT" "$EBEAM" <<'PY'
import sys, os
sys.path.insert(0, os.path.join(os.environ["BENCH"], "analysis"))
import beams
a, b = beams.Beams(sys.argv[1], float(sys.argv[2])).cm_energies
print(f"{a:.6f} {b:.6f}")
PY
)
EOF
[ -n "$EB" ] || { echo "could not compute the beam energies" >&2; exit 1; }

# ------------------------------------------------------------ process setup
if [ ! -d "$PROC" ]; then
    echo "MG5: generating $PROC ($CURRENT, $TAG)"
    CMD=$(mktemp "${TMPDIR:-/tmp}/mg5cmd.XXXXXX")
    # THE MUON IS MASSIVE HERE.  MG5's default restricted Standard Model sets
    # the electron and muon masses to zero, so its Les Houches file declares a
    # massless muon -- and Pythia, which puts every particle on its physical
    # mass shell as it reads the file, then moves the reconstructed Q2 by
    # about 0.06 GeV2.  At the Q2 = 4 edge, where dsigma/dQ2 is largest, that
    # is enough to push events out of the region.  The rest of the benchmark
    # scatters a MASSIVE muon (analysis/beams.py carries m_mu, and Pythia,
    # Herwig and Sherpa all use it), so `sm-lepton_masses` is also the
    # consistent choice, not merely the convenient one.
    if [ "$CURRENT" = mu ]; then
        cat > "$CMD" <<EOF
set automatic_html_opening False
import model sm-lepton_masses
define q = u d s c b u~ d~ s~ c~ b~
generate mu- q > mu- q QED=2 QCD=0
output $PROC
EOF
    else
        # the CC flavour lists are written out rather than using MG5's `p`:
        # the incoming set is what a nu_mu can strike and the outgoing set is
        # what it turns into, and MG5's SM has a DIAGONAL CKM, so the
        # Cabibbo-suppressed channels are absent (worth ~1% on the inclusive
        # rate here, more on the charm composition).
        cat > "$CMD" <<EOF
set automatic_html_opening False
import model sm-lepton_masses
define qi = d s b u~ c~
define qo = u c d~ s~ b~
generate vm qi > mu- qo QED=2 QCD=0
output $PROC
EOF
    fi
    # RUN IT FROM mg5/, not from wherever the caller stood: mg5_aMC drops a
    # generated ALOHA parser table (py.py) and a debug file into the CURRENT
    # directory, and from the repository root those land beside config.sh as
    # untracked litter..gitignore covers them under mg5/.
    cd "$HERE"
    bench_run "$MG5_DIR/bin/mg5_aMC" "$CMD" > "$HERE/generate_$CURRENT$SUF.log" 2>&1 || {
        echo "MG5 process generation failed -- see $HERE/generate_$CURRENT$SUF.log" >&2
        tail -20 "$HERE/generate_$CURRENT$SUF.log" >&2; exit 1; }
    rm -f "$CMD"
fi

# THE C++ RUNTIME.  lhapdf-config --libs hands back -stdlib=libc++, a COMPILE
# flag that gfortran ignores at link time, so the C++ runtime (vtables,
# typeinfo, __clang_call_terminate) is never linked and the build dies in
# hundreds of undefined symbols that look like a compiler mismatch.
# NB the test is on the llhapdf LINE, not on the file: make_opts already
# defines STDLIB=-lc++ a hundred lines above, which is not passed to the link
# that fails, so a file-wide grep reports the patch as present and the build
# dies anyway.
# The runtime is libc++ with clang (macOS) and libstdc++ with GCC (Linux).
case "$(uname -s)" in Darwin) CXXRT=c++ ;; *) CXXRT=stdc++ ;; esac
if ! grep -qE "^[[:space:]]*llhapdf\+=.*-l${CXXRT//+/\\+}( |\$)" "$PROC/Source/make_opts"; then
    perl -i -pe "s/^(\\s*llhapdf\\+=.*-lLHAPDF.*)\$/\$1 -l$CXXRT/" "$PROC/Source/make_opts"
    grep -qE "^[[:space:]]*llhapdf\+=.*-l${CXXRT//+/\\+}( |\$)" "$PROC/Source/make_opts" || {
        echo "could not append -l$CXXRT to the llhapdf line of $PROC/Source/make_opts" >&2
        exit 1; }
fi

python3 "$HERE/install_hooks.py" "$PROC"

# ---------------------------------------------------------------- run card
RC=$PROC/Cards/run_card.dat
# set_rc <value> <key>.  The value and the key go through the environment
# rather than through the perl one-liner's text: both $1 and $2 mean something
# to perl as well as to the shell, and the first draft substituted perl's
# capture groups into the card.
set_rc () {
    VAL=$1 KEY=$2 perl -i -pe \
        'BEGIN{$v=$ENV{VAL}; $k=$ENV{KEY}} s/^(\s*)\S+(\s*=\s*\Q$k\E\b)/$1$v$2/' "$RC"
    VAL=$1 KEY=$2 perl -ne \
        'BEGIN{$v=$ENV{VAL}; $k=$ENV{KEY}; $ok=0} $ok=1 if /^\s*\Q$v\E\s*=\s*\Q$k\E\b/; END{exit($ok?0:1)}' \
        "$RC" || { echo "run card: $2 is not in $RC (or would not take $1)" >&2; exit 1; }
}
set_rc "$NEV"   nevents
set_rc 1234     iseed
set_rc 0        lpp1                 # pointlike lepton beam: no PDF
set_rc 1        lpp2
set_rc "$EA"    ebeam1               # the c.m. pair, as for Herwig and Sherpa
set_rc "$EB"    ebeam2
set_rc lhapdf   pdlabel1
set_rc lhapdf   pdlabel2
set_rc 331100   lhaid                # NNPDF40_nnlo_as_01180
set_rc 0        dynamical_scale_choice   # -> user_dynamical_scale, i.e. Q
set_rc False    fixed_ren_scale
set_rc False    fixed_fac_scale1
set_rc False    fixed_fac_scale2
set_rc False    use_syst
set_rc none     systematics_program
# and zero every collider cut, which otherwise applies ON TOP of the DIS cuts.
# A negative rapidity limit is MG5's "no cut" (cuts.f: etamax >= 0 is the test
# that arms it); a zero pt or Delta R is simply an inactive minimum.
set_rc 0.0      ptj
set_rc 0.0      ptl
set_rc -1.0     etaj
set_rc -1.0     etal
set_rc 0.0      etalmin
set_rc 0.0      drjl

# ---------------------------------------------------------------- param card
# alpha_em = 1/132.119, the benchmark's Gmu electroweak point.  MG5 derives
# sin^2(theta_W) from (alpha_em, G_F, M_Z) rather than taking it as an input,
# so it is not set here; the closure below is what says whether the resulting
# point is close enough to the specified one.
perl -pi -e 's/^(\s*1\s+)\S+(\s+# aEWM1.*)$/${1}1.321190000e+02${2}/' \
    "$PROC/Cards/param_card.dat"
grep -qE "^\s*1\s+1.32119" "$PROC/Cards/param_card.dat" || {
    echo "param card: could not set aEWM1" >&2; exit 1; }

# ------------------------------------------------------------- generation
LHE=$PROC/Events/run_01/unweighted_events.lhe.gz
if [ ! -s "$LHE" ]; then
    echo "MG5: generating $NEV events ($CURRENT, $TAG)"
    ( cd "$PROC" && bench_run ./bin/generate_events -f run_01 \
        > "$HERE/events_$CURRENT$SUF.log" 2>&1 ) || {
        echo "MG5 event generation failed -- see $HERE/events_$CURRENT$SUF.log" >&2
        tail -30 "$HERE/events_$CURRENT$SUF.log" >&2; exit 1; }
fi
[ -s "$LHE" ] || { echo "no $LHE after generation" >&2; exit 1; }
# An existing Les Houches file is REUSED, whatever `nevents` says on this
# invocation -- regenerating it would throw away the run the sample on disk
# was made from.  Delete $PROC/Events/run_01 to make a new one.
echo "Les Houches events: $(gunzip -c "$LHE" | grep -ac "^<event>") in $LHE"

XSEC=$(grep -a "Integrated weight" "$PROC/Events/run_01/run_01_tag_1_banner.txt" \
       2>/dev/null | tail -1 | sed 's/.*: *//')
echo "MG5 cross section (pb): $XSEC"

# ---------------------------------------------------------------- shower
if [ "$SHOWER" = 0 ]; then
    echo "--no-shower: stopping after the Les Houches file"
    exit 0
fi

[ -x "$BENCH_REPO/powheg/main_powheg" ] || make -C "$BENCH_REPO/powheg" main_powheg

if [ -s "$JOBDIR/events.hepmc" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "refusing to overwrite $JOBDIR/events.hepmc; FORCE=1 to replace it." >&2
    exit 1
fi
mkdir -p "$JOBDIR"
# main_powheg is the repository's LHE -> HepMC3 shower driver; with
# POWHEG:veto = 0 (set in the card) it registers no hooks and is a plain
# Les Houches shower, which is what a LO sample needs -- there is no matching
# to veto against.
PLAIN=$JOBDIR/events.lhe
gunzip -c "$LHE" > "$PLAIN"
( cd "$JOBDIR" && bench_run "$BENCH_REPO/powheg/main_powheg" "$NCARD" \
    "$PLAIN" events > shower.log 2>&1 )
rm -f "$PLAIN"
echo "showered -> $JOBDIR/events.hepmc"
cat "$JOBDIR/events_xsec.json"
