#!/usr/bin/env bash
# Fetch the POWHEG forward heavy-hadron NEUTRINO event files of
# arXiv:2309.12793 (Buonocore, Kling, Rottoli, Sominka).
#
# >>> WHY THESE FILES AND NOT THE ONES IN data/faser_flux/ (user, 2026-09-04).
# <<<  "The nu_e charm flux comes from POWHEG, this is the one that we should
# be using."  The vendored histograms in data/faser_flux/ are the combination
# row of arXiv:2105.08270, which AVERAGES four hadronic-interaction models for
# the charm component -- one of them DPMJET, whose nu_e prediction is 3457
# against SIBYLL's 901.  That average sits 55% above the charm flux
# arXiv:2402.13318 Table I actually uses, which is POWHEG + Pythia 8.3, i.e.
# exactly the calculation these files come from: 2309.12793 is Ref. [46] of
# that table.
#
# WHAT THEY ARE.  Neutrinos from charm and bottom hadron decay crossing a
# 1 m x 1 m square centred on the beam collision axis at z = 480 m, at
# sqrt(s) = 13.6 TeV, one row per neutrino:
#
#   vpid  neutrino PDG id          x0, y0   position at z = 480 m [m]
#   hpid  parent hadron PDG id     thx, thy px/pz, py/pz
#   en    energy [GeV]             w        weight [pb]
#   iEvent  key into the weights file, which carries the 7-point scale
#           variation as ratios to the central setup
#
# The transverse position is IN THE FILE, so the FASERnu aperture is applied
# exactly rather than inferred -- which is the whole reason the vendored
# histograms needed an aperture archaeology in the first place
# (faser_rates.FLUX_APERTURE_CM2).
#
# WHERE THEY GO, AND WHY NOT IN GIT.  ~300 MB of generator output, so
# data/forward_charm/ is gitignored under the repository's usual rule ("if a
# generator wrote it, it is ignored").  What IS tracked is this script, the
# README beside the data, and the small per-flavour histograms
# analysis/forward_charm_flux.py derives from them -- a clone can rebuild the
# flux from one command and does not need the 300 MB to reproduce a number.
#
# Usage: tools/fetch_forward_charm.sh [charm|bottom|both]     (default: both)
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT=$HERE/data/forward_charm
BASE=https://raw.githubusercontent.com/KlingFelix/forward_heavy_hadrons_NLONLLx/main/NeutrinoFluxMC
TAG=13.6TeV_POWHEG-smallxLHCb_P8-Monash
WHAT=${1:-both}
mkdir -p "$OUT"

get () {  # $1 = filename
  if [ -s "$OUT/$1" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "exists: $1"; return
  fi
  echo "fetching $1"
  curl -sfL -o "$OUT/$1.part" "$BASE/$1"
  mv "$OUT/$1.part" "$OUT/$1"
}

case "$WHAT" in charm|both)
  for i in 0 1 2 3 4; do get "events_charm_${TAG}_${i}.csv.gzip"; done
  get "weights_charm_${TAG}.csv.gzip" ;;
esac
case "$WHAT" in bottom|both)
  for i in 0 1 2; do get "events_bottom_${TAG}_${i}.csv.gzip"; done
  get "weights_bottom_${TAG}.csv.gzip" ;;
esac

# The authors' own GENIE tungsten cross-sections, shipped beside the events.
# We do NOT use them for a rate -- the benchmark's whole contribution is the
# cross-section, and ours comes from our own splines -- but they are a second,
# independent copy of the number faser_genie_rates.py closes against, so they
# are worth having on disk to compare with.
for pid in 12 -12 14 -14 16 -16; do
  f=xs_GENIE_W_${pid}.txt
  if [ -s "$OUT/$f" ] && [ "${FORCE:-0}" != "1" ]; then echo "exists: $f"; continue; fi
  echo "fetching $f"
  curl -sfL -o "$OUT/$f" "$BASE/GENIE/$f"
done

echo "done -- $(du -sh "$OUT" | cut -f1) in $OUT"
