#!/usr/bin/env bash
# Append scale-variation and PDF weights to an EXISTING POWHEG LHE.
#
# Usage: powheg/reweight/run_reweight.sh [--current mu|nu] [--code V2|RES] \
#                                        <beam_energy_gev>
#
#   run_reweight.sh 1000                      # neutrino CC, POWHEG-V2 (as before)
#   run_reweight.sh --current mu 1000         # muon NC through POWHEG-V2
#
# The hard process is NOT rerun.  POWHEG's `storeinfo_rwgt 1` (set in the
# production cards) keeps, per event, everything needed to re-evaluate the
# matrix element at a different scale or with a different PDF; this reads
# pwgevents.lhe and writes pwgevents-rwgt.lhe with the extra weights.
#
# It works on a COPY of the run directory, never in place: the reweighting
# pass rewrites powheg.input and pwgevents.lhe, and the originals are what
# every published number was produced from.
#
# THE CURRENT IS A PARAMETER, as it is on the three generation drivers
# (powheg/powheg_variants.sh, 2026-08-26).  CONVENTIONS.md rule 2b: an uncertainty
# band offered on one current is owed on the other, and the muon side is
# reachable here because the V2 muon cross-variant is the same CODE -- the
# reweighting machinery travels with the code, not with the current.
#
# >>> POWHEG-RES IS NOT SUPPORTED YET, AND THE REASON IS STRUCTURAL. <<<
# RES runs with `manyseeds` and writes TEN files, pwgevents-0001.lhe through
# -0010.lhe, where V2 writes one.  Reweighting them means one pwhg_main pass
# per seed with the matching seed index, then a combination that respects the
# per-seed statistics -- not the single pass below.  The binary does carry the
# rwl_* symbols, so it is a driver question rather than a capability one.
# Until then the muon NLO band comes from the V2 cross-variant at 1 TeV, which
# is stated wherever it is shown rather than left looking like the primary
# muon entry.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$BENCH_REPO/powheg/powheg_variants.sh"

CODE=V2
ARGS=()
CURRENT=nu
while [ $# -gt 0 ]; do
    case "$1" in
        --current) CURRENT=${2:?--current takes mu or nu}; shift 2 ;;
        --code)    CODE=${2:?--code takes V2 or RES}; shift 2 ;;
        *)         ARGS+=("$1"); shift ;;
    esac
done
set -- "${ARGS[@]}"
EBEAM=${1:?usage: run_reweight.sh [--current mu|nu] [--code V2|RES] <energy_gev>}
case "$CURRENT" in mu|nu) ;; *) echo "--current takes mu or nu" >&2; exit 1 ;; esac

if [ "$CODE" != "V2" ]; then
    echo "run_reweight.sh: --code $CODE is not supported -- see the header." >&2
    echo "  POWHEG-RES writes ten manyseeds LHEs and needs one pass per seed." >&2
    exit 1
fi

powheg_variant "$CODE" "$CURRENT" "$EBEAM"
TAG=$(python3 "$BENCH_REPO/analysis/beams.py" tag "$CURRENT" "$EBEAM")

# The neutrino anchor's events live in the CKM-fixed directory, which is NOT
# what powheg_variant names -- that table points at the generation directory
# for a NEW run.  Prefer the -ckm tree when it exists, exactly as before.
# analysis/powheg_nlo_uncertainty.py closes the result against the published
# showered sigma_fid, so picking the wrong source here is caught rather than
# quietly reweighting a sample nobody published.
SRC=$PV_RUNDIR-ckm
[ -d "$SRC" ] || SRC=$PV_RUNDIR
[ -s "$SRC/pwgevents.lhe" ] || { echo "no LHE at $SRC" >&2; exit 1; }
DEST=$POWHEG_V2/rwgt-$CURRENT$TAG
PWHG=$POWHEG_V2/pwhg_main

rm -rf "$DEST"; mkdir -p "$DEST"
cp "$SRC/powheg.input" "$DEST/"
cp "$SRC/pwgevents.lhe" "$DEST/"
# BENCH_RWL selects the reweight file.  The default now carries every member
# of every PDF set, so the Monte Carlo can produce PDF UNCERTAINTIES and not
# just central values -- the README's "about ten times this whole run" was an
# estimate, never a measurement, and the measurement is 87 us per (event,
# weight), i.e. about half an hour per energy.  Set BENCH_RWL=rwl_scale_pdf.xml
# for the old central-only file.
RWL=${BENCH_RWL:-rwl_members.xml}
[ -f "$HERE/$RWL" ] || { echo "no reweight file $HERE/$RWL" >&2; exit 1; }
cp "$HERE/$RWL" "$DEST/"
# grids are copied so POWHEG does not try to re-integrate
cp "$SRC"/pwggrid*.dat "$SRC"/pwgubound*.dat "$SRC"/pwgxgrid.dat "$DEST/" 2>/dev/null || true

# the reweighting keys.  `rwl_add 1` is what actually SELECTS the reweighting
# path: pwhg_main.f tests `flg_newweight .or. flg_rwl_add` and, when neither is
# set, falls through to ordinary event generation -- which then refuses to
# overwrite the very pwgevents.lhe it was meant to read, and exits.  Setting
# rwl_file alone is not enough.
# rwl_group_events batches events through the variations so the PDFs are not
# re-initialised per event, which otherwise dominates the runtime.
cat >> "$DEST/powheg.input" <<KEYS

! ---- a-posteriori reweighting (powheg/reweight/) ----
rwl_file '$RWL'
rwl_add 1
rwl_group_events 2000
rwl_format_rwgt 1
KEYS

echo "reweighting [$PV_LABEL] $SRC -> $DEST"
( cd "$DEST" && bench_run "$PWHG" > reweight.log 2>&1 ) || {
    echo "pwhg_main failed; tail of the log:" >&2
    tail -25 "$DEST/reweight.log" >&2; exit 1; }

OUT=$DEST/pwgevents-rwgt.lhe
[ -s "$OUT" ] || OUT=$DEST/pwgevents.lhe
echo "wrote $OUT"
echo "weights per event:"
grep -m1 -c "<wgt" "$OUT" 2>/dev/null || true
grep -m1 -A20 "<initrwgt" "$OUT" | head -25
echo
echo "next: BENCH_SELECTION=inclusive analysis/powheg_nlo_uncertainty.py \\"
echo "        --current $CURRENT     (BENCH_ENERGY=$EBEAM)"
