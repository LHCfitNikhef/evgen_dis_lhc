#!/usr/bin/env bash
# Validation runs for the POWHEG part of the paper-plots production, at
# 1 TeV, BEFORE any production.  Settles, with numbers:
#
#  (1) the NEUTRON: a genuine neutron beam (ih2 = 2, proton PDF set, the code
#      swaps u<->d) must reproduce a proton beam carrying the isospin-mirrored
#      set NNPDF40_nnlo_as_01180_n (LHAPDF 339900, data/pdfs/lhapdf).  Both
#      runs use the SAME random seeds, so they differ only in how the neutron
#      PDF is obtained.  The shower of the neutron beam must show a 2112 beam
#      and a charge-neutral target (final-state charge = lepton charge).
#  (2) the GENERATION Q2 CUT: the same card generated at Q2 > 2.25 (the earlier production) and at
#      Q2 > 4 must give the same sigma(Q2 > 4, W > 3) = sigma_int x fraction,
#      at LHE level and after the shower, if nothing moves the lepton across
#      Q2 = 4.  A seed beyond the integration set is also generated, which
#      checks that stage 4 of POWHEG-RES works from the combined grids alone.
#
# Usage: powheg/production/validate.sh            (re-invocable; the drivers skip work)
#        powheg/production/validate.sh summary    (print the comparison only)
#
# Output: $POWHEG_RES/v2check-<name>/ and $POWHEG_V2/v2check-<name>/, each with
# the run in run/ and the showered jobs beside it.  These are DIAGNOSTICS
# (CONVENTIONS.md rule 1): once the numbers are recorded in PRODUCTION.md,
# `powheg/production/validate.sh clean` deletes the event files and keeps the logs and
# stat files.
#
# Cores: the RES lane runs 3 processes, the V2 lane 1 -- 4 in total.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/../.." && pwd)}
. "$BENCH_REPO/config.sh"
. "$HERE/common.sh"

CARDS=$BENCH_REPO/powheg/cards/v2
NPDF_DIR=$BENCH_REPO/data/pdfs/lhapdf

# name | nucleon for the driver | card | key overrides | numevts | showered target | intseeds (RES) | seed base
RES_RUNS="
res-p-q225|p|mu1TeV-p||50000|550000|6|90000
res-p-q4|p|mu1TeV-p|Qmin=2d0|50000|300000|6|90000
res-n-beam|n|mu1TeV-n||10000|30000|4|91000
res-p-npdf|p|mu1TeV-p|lhans1=339900 lhans2=339900|10000|30000|4|91000
"
V2_RUNS="
v2-p-q225|p|nu1TeV-p||200000|190000|90000
v2-p-q4|p|nu1TeV-p|q2cut=4d0|200000|190000|90000
v2-n-beam|n|nu1TeV-n||50000|45000|91000
v2-p-npdf|p|nu1TeV-p|lhans1=339900 lhans2=339900|50000|45000|91000
"

res_lane () {
    echo "$RES_RUNS" | while IFS='|' read -r name t card keys nev target nint sb; do
        [ -n "$name" ] || continue
        chk=$POWHEG_RES/v2check-$name
        mkdir -p "$chk"
        # shellcheck disable=SC2086
        python3 "$HERE/set_keys.py" "$CARDS/POWHEG-RES/$card/powheg.input" \
            "$chk/card.input" "numevts=$nev" $keys
        env_pdf=""
        case "$keys" in *339900*) env_pdf=$NPDF_DIR ;; esac
        # NOT LHAPDF_DATA_PATH=...: the driver sources config.sh, which
        # reassigns it; common.sh prepends V2_LHAPDF_PREPEND afterwards
        V2_LHAPDF_PREPEND=$env_pdf \
        V2_CARD=$chk/card.input V2_RUNDIR=$chk/run V2_JOBROOT=$chk V2_JOBBASE=job \
        V2_TARGET=$target V2_INTSEEDS=$nint V2_SEEDBASE=$sb V2_PAR=3 \
            "$HERE/run_res.sh" 1000 "$t" > "$chk/driver.log" 2>&1 \
            || { echo "RES validation $name FAILED -- $chk/driver.log" >&2; return 1; }
        echo "RES validation $name done"
    done
}

v2_lane () {
    echo "$V2_RUNS" | while IFS='|' read -r name t card keys nev target sb; do
        [ -n "$name" ] || continue
        chk=$POWHEG_V2/v2check-$name
        mkdir -p "$chk"
        # shellcheck disable=SC2086
        python3 "$HERE/set_keys.py" "$CARDS/POWHEG-V2/$card/powheg.input" \
            "$chk/card.input" "numevts=$nev" $keys
        env_pdf=""
        case "$keys" in *339900*) env_pdf=$NPDF_DIR ;; esac
        # NOT LHAPDF_DATA_PATH=...: the driver sources config.sh, which
        # reassigns it; common.sh prepends V2_LHAPDF_PREPEND afterwards
        V2_LHAPDF_PREPEND=$env_pdf \
        V2_CARD=$chk/card.input V2_RUNDIR=$chk/run V2_JOBROOT=$chk V2_JOBBASE=job \
        V2_TARGET=$target V2_SEEDBASE=$sb V2_PAR=1 \
            "$HERE/run_v2.sh" 1000 "$t" > "$chk/driver.log" 2>&1 \
            || { echo "V2 validation $name FAILED -- $chk/driver.log" >&2; return 1; }
        echo "V2 validation $name done"
    done
}

summary () {
    python3 - "$POWHEG_RES" "$POWHEG_V2" "$HERE" <<'PY'
import glob, json, os, re, subprocess, sys
res, v2, here = sys.argv[1:4]
sys.path.insert(0, here)
import inspect_sample as I

def res_sigma(d):
    t = open(f"{d}/run/INTEGRATED").read()
    m = re.search(r"\(pos\.-\|neg\.\|\):\s+([0-9.eE+-]+)\s+\+-\s+([0-9.eE+-]+)", t)
    return float(m.group(1)), float(m.group(2))

def v2_sigma(d):
    t = open(f"{d}/run/pwg-stat.dat").read()
    m = re.search(r"total \(btilde\+remnants\) cross section in pb\s+([0-9.eE+-]+)\s+\+-\s+([0-9.eE+-]+)", t)
    return float(m.group(1)), float(m.group(2))

def frac_lhe(lhes):
    num = den = 0.0; nev = 0
    for lhe in lhes:
        for rec in I.read_lhe(lhe):
            if rec[0] != "event": continue
            _, w, k, kp, fl = rec
            Q2, y, _ = I.q2_y(k, kp)
            W2 = I.M_P**2 + y * 2 * 1000.0 * I.M_P - Q2
            den += w; nev += 1
            if Q2 > 4 and W2 > 9: num += w
    return num / den, nev

def frac_hepmc(files):
    num = den = 0.0; nev = 0
    for f in files:
        for ev in I.hepmc_events(f):
            r = I.event_kinematics(ev)
            if r is None: continue
            Q2, W2 = r[0], r[1]
            w = ev["weight"]; den += w; nev += 1
            if Q2 > 4 and W2 > 9: num += w
    return num / den, nev

rows = {}
for code, base, sig, lhepat in (("RES", res, res_sigma, "run/pwgevents-*.lhe"),
                                ("V2", v2, v2_sigma, "run*/pwgevents.lhe")):
    for d in sorted(glob.glob(f"{base}/v2check-*")):
        name = os.path.basename(d)[len("v2check-"):]
        try:
            s, e = sig(d)
        except Exception as ex:
            print(f"{name}: no integrator yet ({ex})"); continue
        jobs = sorted(glob.glob(f"{d}/job_*/events.hepmc"))
        if code == "RES":
            lhes = [f"{d}/run/pwgevents-%04d.lhe" % int(j.split('_')[-1].split('/')[0]) for j in jobs]
        else:
            lhes = sorted(glob.glob(f"{d}/run/pwgevents.lhe") + glob.glob(f"{d}/run-b*/pwgevents.lhe"))
        fl, nl = frac_lhe(lhes)
        fh, nh = frac_hepmc(jobs)
        beams = I.hepmc_summary(jobs[0], 3000)
        rows[name] = (s, e, fl, nl, fh, nh)
        print(f"{name:12s} sigma_int = {s:.6g} +- {e:.3g} ({100*e/s:.2f}%)  "
              f"LHE f(Q2>4,W>3) = {fl:.5f} [{nl} ev] -> {s*fl:.6g}   "
              f"HepMC f = {fh:.5f} [{nh} ev] -> {s*fh:.6g}")
        print(f"{'':12s} first job: beams {beams['beams']} charge*3 {beams['charge_sum_times3']} remnants {beams['remnants']}")

def cmp(a, b, what, idx_s=0, frac=None):
    if a not in rows or b not in rows: return
    ra, rb = rows[a], rows[b]
    if frac is None:
        va, vb = ra[0], rb[0]; ea, eb = ra[1], rb[1]
    else:
        fa, fb = ra[frac], rb[frac]; na, nb = ra[frac+1], rb[frac+1]
        va, vb = ra[0]*fa, rb[0]*fb
        ea = va*((ra[1]/ra[0])**2 + (1-fa)/(fa*na))**0.5
        eb = vb*((rb[1]/rb[0])**2 + (1-fb)/(fb*nb))**0.5
    r = va/vb; er = r*((ea/va)**2+(eb/vb)**2)**0.5
    print(f"  {what}: {a}/{b} = {r:.5f} +- {er:.5f}  ({(r-1)/er:+.1f} sigma)")

print()
cmp("res-n-beam", "res-p-npdf", "RES neutron beam vs proton beam + neutron PDF (integrator)")
cmp("v2-n-beam", "v2-p-npdf", "V2  neutron beam vs proton beam + neutron PDF (integrator)")
cmp("res-n-beam", "res-p-q225", "RES n/p integrator (NC mu: expect < 1)")
cmp("v2-n-beam", "v2-p-q225", "V2  n/p integrator (CC nu: expect > 1)")
cmp("res-p-q4", "res-p-q225", "RES gen Q2>4 vs Q2>2.25, sigma(Q2>4,W>3) LHE  ", frac=2)
cmp("res-p-q4", "res-p-q225", "RES gen Q2>4 vs Q2>2.25, sigma(Q2>4,W>3) HepMC", frac=4)
cmp("v2-p-q4", "v2-p-q225", "V2  gen Q2>4 vs Q2>2.25, sigma(Q2>4,W>3) LHE  ", frac=2)
cmp("v2-p-q4", "v2-p-q225", "V2  gen Q2>4 vs Q2>2.25, sigma(Q2>4,W>3) HepMC", frac=4)
PY
}

clean () {
    for d in "$POWHEG_RES"/v2check-* "$POWHEG_V2"/v2check-*; do
        [ -d "$d" ] || continue
        find "$d" \( -name "*.lhe" -o -name "*.hepmc" \) -delete
        echo "cleaned $d"
    done
}

case "${1:-run}" in
    summary) summary; exit 0 ;;
    clean)   clean; exit 0 ;;
    run) ;;
    *) echo "usage: validate.sh [run|summary|clean]" >&2; exit 1 ;;
esac

v2_disk_check 12 "validation"
res_lane & pr=$!
v2_lane & pv=$!
fail=0
wait "$pr" || fail=1
wait "$pv" || fail=1
[ "$fail" = 0 ] || { echo "VALIDATION RUNS FAILED" >&2; exit 1; }
echo "VALIDATION RUNS DONE"
summary
