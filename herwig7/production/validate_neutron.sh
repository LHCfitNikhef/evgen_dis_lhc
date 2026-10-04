#!/usr/bin/env bash
# Validate the NEUTRON target of the Herwig production before it runs.
#
# Usage: herwig7/production/validate_neutron.sh [nev_xsec] [nev_inspect]
#        (defaults 200000 and 5000; at most 4 Herwig processes at a time)
#
# THE QUESTION.  Does a genuine neutron beam (/Herwig/Particles/n0) with the
# PROTON set give the neutron -- or does ThePEG want the neutron set, and
# would the pair then swap twice?  make_cards.py records what the source
# says (ThePEG::LHAPDF swaps u<->d itself for id 2112); this measures it.
#
# Arms, 1 TeV, BOTH currents (CONVENTIONS.md rule 2b), positive half:
#   p    proton beam  + proton set        the production proton card
#   n    neutron beam + proton set        the production neutron card
#   pn   proton beam  + neutron set       the fallback route
#   nn   neutron beam + neutron set       the double swap
# Expected if ThePEG swaps:  sigma(n) = sigma(pn) within statistics, and
# sigma(nn) = sigma(p); n/p well above 1 for nu CC, below 1 for mu NC.
# Plus the NEGATIVE half of p and n, so that n/p is also quoted for the full
# NLO (pos - neg) and the neg fraction that sets the pos/neg job split is
# measured in the new region.
#
# Two passes per arm:
#   xsec     nev_xsec events with events.hepmc pointed at /dev/null: only the
#            sampler integral in the .log is wanted, and 200k events of HepMC
#            per arm would be ~2 GB of disk for nothing (rule 1)
#   inspect  nev_inspect events WITH HepMC, read by inspect_sample.py: beams as
#            written, remnant diquarks, the PDF used event by event
#
# Directories herwig7/v2val_job_<arm>_<pass>/ are diagnostics: delete them
# once the numbers are recorded in PRODUCTION.md.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
. "$HERE/lib.sh"
NX=${1:-200000}
NI=${2:-5000}
MAXPAR=4

ARMS="V2-mu-PWG-p V2-mu-PWG-n V2VAL-mu-PWG-pn V2VAL-mu-PWG-nn \
V2-nu-PWG-p V2-nu-PWG-n V2VAL-nu-PWG-pn V2VAL-nu-PWG-nn \
V2-mu-PWGNEG-p V2-mu-PWGNEG-n V2-nu-PWGNEG-p V2-nu-PWGNEG-n"

avail=$(v2_avail_gb)
echo "disk: ${avail} GB free (floor $V2_DISK_FLOOR_GB)"
[ "$avail" -gt $((V2_DISK_FLOOR_GB + 5)) ] || { echo "too little disk; stopping" >&2; exit 1; }

# compile serially first: jobs of one stem must never race on its build dir
for s in $ARMS; do v2_compile "$s"; done

run_arm() {  # run_arm <stem> <pass> <nev> <seed>
    local stem=$1 d=$V2_HDIR/v2val_job_$1_$2
    rm -rf "$d"; mkdir -p "$d"
    cp "$V2_BUILD/$stem/$stem.run" "$V2_BUILD/$stem/LeptonicDISCut.so" "$d/"
    [ -d "$V2_BUILD/$stem/Herwig-cache" ] && cp -R "$V2_BUILD/$stem/Herwig-cache" "$d/"
    [ "$2" = xsec ] && ln -s /dev/null "$d/events.hepmc"
    (cd "$d" && hw_env Herwig run "$stem.run" -N "$3" -s "$4" -d 0 > herwig_run.log 2>&1) \
        && echo "  done $stem $2 $(date '+%H:%M:%S')" \
        || echo "  FAILED $stem $2 -- $d/herwig_run.log" >&2
}

i=0
for pass in inspect xsec; do
    nev=$NI; [ "$pass" = xsec ] && nev=$NX
    for s in $ARMS; do
        i=$((i + 1))
        while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$MAXPAR" ]; do sleep 5; done
        run_arm "$s" "$pass" "$nev" $((3000000 + i)) &
    done
done
wait
echo "all validation arms finished $(date '+%H:%M:%S')"

python3 - "$V2_HDIR" $ARMS <<'PY'
import os, sys
sys.path.insert(0, sys.argv[1])
from herwig_xsec import read_job
base, arms = sys.argv[1], sys.argv[2:]
res = {}
for a in arms:
    d = os.path.join(base, f"v2val_job_{a}_xsec")
    outs = [f for f in os.listdir(d) if f.endswith(".out")]
    j = read_job(os.path.join(d, outs[0]))
    res[a] = j
    print(f"  {a:18s} sigma = {j.sigma_pb:12.5f} +- {j.err_pb:9.5f} pb   "
          f"survival {j.survival:.4f}")
def ratio(a, b):
    x, y = res[a], res[b]
    r = x.sigma_pb / y.sigma_pb
    e = r * ((x.err_pb/x.sigma_pb)**2 + (y.err_pb/y.sigma_pb)**2) ** 0.5
    return r, e
for cur in ("mu", "nu"):
    print(f"\n  {cur}:")
    r, e = ratio(f"V2-{cur}-PWG-n", f"V2VAL-{cur}-PWG-pn")
    print(f"    n-beam+p-set / p-beam+n-set = {r:.5f} +- {e:.5f}  ({(r-1)/e:+.2f} sigma)")
    r, e = ratio(f"V2VAL-{cur}-PWG-nn", f"V2-{cur}-PWG-p")
    print(f"    n-beam+n-set / p-beam+p-set = {r:.5f} +- {e:.5f}  ({(r-1)/e:+.2f} sigma)  [double swap = proton?]")
    r, e = ratio(f"V2-{cur}-PWG-n", f"V2-{cur}-PWG-p")
    print(f"    n/p, positive half          = {r:.5f} +- {e:.5f}")
    P = res[f"V2-{cur}-PWG-p"]; N = res[f"V2-{cur}-PWG-n"]
    Pn = res[f"V2-{cur}-PWGNEG-p"]; Nn = res[f"V2-{cur}-PWGNEG-n"]
    sp = P.sigma_pb - Pn.sigma_pb; sn = N.sigma_pb - Nn.sigma_pb
    ep = (P.err_pb**2 + Pn.err_pb**2) ** 0.5; en = (N.err_pb**2 + Nn.err_pb**2) ** 0.5
    print(f"    full NLO (gen region): p {sp:.5f} +- {ep:.5f}, n {sn:.5f} +- {en:.5f} pb;"
          f"  n/p = {sn/sp:.5f} +- {sn/sp*((ep/sp)**2+(en/sn)**2)**0.5:.5f}")
    print(f"    neg/pos: p {Pn.sigma_pb/P.sigma_pb:.4f}, n {Nn.sigma_pb/N.sigma_pb:.4f}")
PY

echo
echo "== inspection (HepMC passes) =="
for s in $ARMS; do
    pset=NNPDF40_nnlo_as_01180
    case "$s" in *-pn|*-nn) pset=NNPDF40_nnlo_as_01180_n ;; esac
    python3 "$HERE/inspect_sample.py" "$V2_HDIR/v2val_job_${s}_inspect" --pdfset "$pset"
done
