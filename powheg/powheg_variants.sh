# The POWHEG variant table -- the ONE place that maps (code, current) to paths.
#
# Sourced (not executed) by run_powheg_res.sh, run_powheg_v2.sh and
# run_powheg_shower.sh, AFTER config.sh, which is where $POWHEG_RES and
# $POWHEG_V2 come from.
#
# WHY THERE IS A TABLE AND NOT A RULE.  There are four (code, current)
# combinations and their four card-directory suffixes are all different --
# -wide, -prod, -NC, -CC.  Nothing derives them: they are the names the
# existing samples, the symlinks and tools/link_cards.sh already use, and
# renaming one would orphan a published sample.  So they are listed here once
# and an unknown combination RAISES, rather than falling through to a
# plausible-looking name -- the same reason SHERPA_NU_NLO in analyze_nu.py is
# a table keyed by beam tag.
#
# WHICH COMBINATION IS WHICH, and why each exists:
#
#   RES + mu   the published muon NC NLO sample.  POWHEG-RES's native side.
#   V2  + nu   the published neutrino CC NLO sample.  POWHEG-V2's native side,
#              and the code whose point is massive charm.
#   V2  + mu   muon NC through the V2 code: a THIRD independent NLO matching
#              on the muon side.  Its card REQUIRES `q2cut 2.25d0`, which the
#              CC cards do not -- NC has a photon pole and POWHEG's default
#              q2cut of 5d-2 integrates it and returns 15 microbarn against
#              the ~38 nb it should give, with no error message.
#   RES + nu   neutrino CC through the RES code.  >>> INCLUSIVE SELECTION
#              ONLY <<<: DIS_v hardcodes CKM_diag = .true., the CKM matrix
#              never reaches the matrix elements and d -> c is absent.  CKM
#              unitarity protects the inclusive rate (0.5% against POWHEG-V2)
#              but not the charm composition -- the identical defect cost
#              Sherpa 3.4 points on the charm fraction.  The label carries the
#              warning so it appears in every run log.
#
# CROSS-VARIANT CARDS EXIST AT 1 TeV ONLY.  tools/make_energy_cards.py knows
# the two native variants and not these, deliberately: they are matching
# cross-checks at the benchmark point, not scan members.  Asking for another
# energy therefore fails on the missing card, which is the honest outcome.
#
# powheg_variant sets, for one (code, current, energy):
#
#   PV_CARD     the committed card in this repo
#   PV_RUNDIR   the generation directory inside the generator's install tree
#   PV_JOBBASE  the showered-sample directory base, per beams.py naming (the
#               anchor keeps its untagged name, so nothing published moves)
#   PV_CMND     the Pythia8 steering card, keyed by CURRENT and not by code:
#               main_powheg.cc takes the beams from the LHE, so the card only
#               has to know which lepton it is showering, not which code or
#               which energy produced it
#   PV_LABEL    a human label for the log lines
#
powheg_variant () {   # $1 = RES|V2   $2 = mu|nu   $3 = beam energy in GeV
    local code=$1 current=$2 ebeam=$3
    local here repo tag jobbase

    here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    repo=${BENCH_REPO:-$(cd "$here/.." && pwd)}

    # The combination is checked BEFORE beams.py is called, so a typo reports
    # itself here rather than as a Python traceback out of beams.py.
    case "$code:$current" in
        RES:mu) jobbase=job       ;;
        RES:nu) jobbase=resnu_job ;;
        V2:nu)  jobbase=nu_job    ;;
        V2:mu)  jobbase=v2mu_job  ;;
        # POWHEG-V2mc: POWHEG-V2 with MASSIVE charm, CC only (user,
        # 2026-08-28).  A SEPARATE code name throughout -- separate cards,
        # separate run directories, separate job base, separate result keys --
        # so it can never be confused with the massless POWHEG-V2 results,
        # which stay exactly as they are.
        V2MC:nu) jobbase=numc_job ;;
        *) echo "powheg_variant: no such variant '$code' + '$current'" \
                "(want RES|V2|V2MC and mu|nu; V2MC is nu only)" >&2
           return 1 ;;
    esac

    tag=$(python3 "$repo/analysis/beams.py" tag "$current" "$ebeam") || return 1
    PV_JOBBASE=$(python3 "$repo/analysis/beams.py" name "$jobbase" "$ebeam") \
        || return 1

    case "$code:$current" in
        RES:mu) PV_CARD=$here/cards/POWHEG-RES/mu$tag-wide/powheg.input
                PV_RUNDIR=$POWHEG_RES/parallel-mu$tag-wide
                PV_LABEL="POWHEG-RES mu- NC" ;;
        RES:nu) PV_CARD=$here/cards/POWHEG-RES/nu$tag-CC/powheg.input
                PV_RUNDIR=$POWHEG_RES/parallel-nu$tag-CC
                PV_LABEL="POWHEG-RES nu_mu CC [INCLUSIVE ONLY: diagonal CKM, no d->c]" ;;
        V2:nu)  PV_CARD=$here/cards/POWHEG-V2/nu$tag-prod/powheg.input
                # 1 TeV LIVES IN prod-nu1TeV-ckm, NOT prod-nu1TeV.  The
                # CKM bug (patches/powheg-v2-nu-dis-ckm.diff: a b-initiated
                # Born returned an uninitialised flavour weight, and an
                # outgoing b was given |V_cs|^2) was fixed on 2026-08-21 and
                # the sample REGENERATED into a new directory; the superseded
                # 2026-08-13 tree was left in place beside it.
                # This table kept pointing at the OLD one, so anything that
                # resolved its LHE through here -- the member reweighting did
                # -- silently used a sample the published results do not come
                # from.  Verified against powheg/nu_job_1/shower.log, which
                # names prod-nu1TeV-ckm as its Beams:LHEF.
                # 400 GeV and 4 TeV were produced after the fix and need no
                # such special case.
                PV_RUNDIR=$POWHEG_V2/prod-nu$tag
                [ "$tag" = 1TeV ] && [ -d "$POWHEG_V2/prod-nu1TeV-ckm" ] \
                    && PV_RUNDIR=$POWHEG_V2/prod-nu1TeV-ckm
                PV_LABEL="POWHEG-V2 nu_mu CC" ;;
        V2:mu)  PV_CARD=$here/cards/POWHEG-V2/mu$tag-NC/powheg.input
                PV_RUNDIR=$POWHEG_V2/prod-mu$tag
                PV_LABEL="POWHEG-V2 mu- NC" ;;
        V2MC:nu) PV_CARD=$here/cards/POWHEG-V2mc/nu$tag-prod/powheg.input
                PV_RUNDIR=$POWHEG_V2/prod-numc$tag
                PV_LABEL="POWHEG-V2mc nu_mu CC [MASSIVE CHARM]" ;;
    esac

    case "$current" in
        mu) PV_CMND=powheg_mu1TeV.cmnd ;;
        nu) PV_CMND=powheg_nu1TeV.cmnd ;;
    esac

    [ -f "$PV_CARD" ] || {
        echo "no card $PV_CARD" >&2
        echo "  the two native variants take new energies from" \
             "tools/make_energy_cards.py;" >&2
        echo "  the two cross-variants (V2+mu, RES+nu) exist at 1 TeV only." >&2
        return 1; }
    return 0
}

# Read a leading `--current mu|nu` off a driver's argument list.  Every driver
# defaults to its own native current, so every invocation that predates the
# cross-variants keeps working unchanged:
#
#     ./run_powheg_res.sh 400                 # muon NC, as before
#     ./run_powheg_res.sh --current nu 1000   # the CC cross-variant
#
# Sets PV_CURRENT and leaves the remaining arguments in PV_ARGS, which the
# caller re-expands with "${PV_ARGS[@]}" -- NOT $PV_ARGS, which in bash 3.2
# takes only the first element.
powheg_parse_current () {   # $1 = default current, then the driver's "$@"
    PV_CURRENT=$1; shift
    while [ $# -gt 0 ]; do
        case $1 in
            --current) PV_CURRENT=$2; shift 2 ;;
            --current=*) PV_CURRENT=${1#--current=}; shift ;;
            --) shift; break ;;
            *) break ;;
        esac
    done
    case $PV_CURRENT in
        mu|nu) ;;
        *) echo "--current must be mu or nu, not '$PV_CURRENT'" >&2; return 1 ;;
    esac
    PV_ARGS=("$@")
}

# Copy the committed card into a run directory as a REAL FILE.
#
# `cp -f` alone is not enough, and the failure is not the obvious one.
# tools/link_cards.sh makes a run directory's powheg.input a SYMLINK back into
# this repo, so the destination frequently points AT the source -- and cp then
# refuses with "dst and src are identical (not copied)" and returns 1, which
# under `set -e` aborts the driver before POWHEG is ever called.  The repo card
# is NOT damaged; the cost is a dead run with a confusing message.  Removing
# the destination first makes it a plain copy in every case.
#
# It has to be a real file because run_powheg_res.sh rewrites parallelstage
# between its five passes, and writing through the symlink would edit the
# committed card.  After stage 4 the copy is identical to the repo card again,
# so link_cards.sh restores the symlink on its next run and reports it as
# "replace copy (identical to repo)".
powheg_install_card () {   # $1 = source card, $2 = destination path
    rm -f "$2"
    cp "$1" "$2"
}
