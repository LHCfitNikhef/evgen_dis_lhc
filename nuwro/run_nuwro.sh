#!/usr/bin/env bash
# Generate one NuWro sample and convert it to HepMC3.
#
# Usage: nuwro/run_nuwro.sh [--current nu|nubar] [--target p|n|W] \
#                           <beam energy in GeV> [n events] [job number]
#
# WHY NuWro IS HERE (user, 2026-09-08): a "Comparison with Neutrino
# Generators" study, setting GENIE beside the dedicated neutrino generators.
# NuWro is the one of those that reaches these energies -- its README claims
# "threshold to TeV" and the FPF whitepaper (arXiv:2203.05090, Fig. 7.20) runs
# it at 1 TeV on tungsten, which is this benchmark's own anchor point.
#
# >>> DIS ONLY, AND THAT IS THE POINT OF THE COMPARISON. <<<  Every other
# channel NuWro has -- quasi-elastic, resonant, coherent, meson-exchange,
# hyperon -- is switched off, because this benchmark is a DIS benchmark and
# the non-DIS channels were measured (appendix: non-DIS processes) to be per
# mille of the fiducial rate.  Leaving them on would compare NuWro's total
# against everyone else's DIS.
#
# >>> THE TARGET IS A CHOICE WITH PHYSICS IN IT. <<<
#   p, n   a FREE nucleon (nucleus_target = 0, no Pauli blocking, no cascade).
#          Like for like with every other generator in this benchmark, which
#          all scatter off a free proton -- so a difference here is DIS and
#          hadronisation modelling and nothing else.
#   W      tungsten, 74 protons and 110 neutrons, with the local Fermi gas and
#          the intranuclear cascade ON -- set BENCH_NUWRO_FSI=0 to keep the
#          nucleus but switch the cascade off, which is how its TeV failure
#          was isolated (see nuwro/README.md).  The difference between this and the
#          free-nucleon run is the nuclear physics this benchmark otherwise
#          omits, and it is what NuWro is in the comparison FOR.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$BENCH_REPO/config.sh"

CURRENT=nu
TARGET=p
while [ $# -gt 0 ]; do
  case "$1" in
    --current) CURRENT=$2; shift 2 ;;
    --target)  TARGET=$2;  shift 2 ;;
    *) break ;;
  esac
done
EBEAM=${1:?usage: run_nuwro.sh [--current nu|nubar] [--target p|n|W] <E_GeV> [nev] [job]}
NEV=${2:-200000}
JOB=${3:-1}

case "$CURRENT" in
  nu)    PID=14  ;;
  nubar) PID=-14 ;;
  *) echo "--current takes nu or nubar" >&2; exit 1 ;;
esac
case "$TARGET" in
  p) NP=1;  NN=0;   TGT=0; FSI=0; PAULI=0 ;;
  n) NP=0;  NN=1;   TGT=0; FSI=0; PAULI=0 ;;
  W) NP=74; NN=110; TGT=2; FSI=${BENCH_NUWRO_FSI:-1}; PAULI=1 ;;
  *) echo "--target takes p, n or W" >&2; exit 1 ;;
esac

[ -x "$HERE/nuwro2hepmc" ] || { echo "build the converter first: make -C nuwro" >&2; exit 1; }
[ -x "$NUWRO_DIR/bin/nuwro" ] || { echo "no $NUWRO_DIR/bin/nuwro" >&2; exit 1; }

# >>> THE ANCHOR KEEPS ITS UNTAGGED NAME. <<<  beams.py's `name` is the one
# place that convention lives, and job_files() in the analysis matches
# "<base>_<integer>" exactly -- so a directory called nuwro_job_nu_p_1TeV_1
# at the anchor energy would simply never be found, silently, which is the
# trap analyze.py's job_files docstring was written after.
BASE_NAME=$(python3 "$BENCH_REPO/analysis/beams.py" name \
                "nuwro_job_${CURRENT}_${TARGET}" "$EBEAM")
d=$HERE/${BASE_NAME}_${JOB}
mkdir -p "$d"
echo "NuWro: $CURRENT on $TARGET at $EBEAM GeV, $NEV events -> $(basename "$d")"

# NuWro writes its scratch files into the CURRENT directory, so each job runs
# in its own.  The energy is given in MeV (src/jednostki.h sets MeV = 1).
(
  cd "$d" && bench_run "$NUWRO_DIR/bin/nuwro" \
      -i "$NUWRO_DIR/data/params.txt" -o events.root \
      -p "number_of_test_events = $((NEV * 10))" \
      -p "number_of_events = $NEV" \
      -p "random_seed = $JOB" \
      -p "beam_type = 0" \
      -p "beam_energy = $(python3 -c "print(int(round(float('$EBEAM') * 1000)))")" \
      -p "beam_particle = $PID" \
      -p "nucleus_p = $NP" -p "nucleus_n = $NN" \
      -p "nucleus_target = $TGT" \
      -p "FSI_on = $FSI" -p "pauli_blocking = $PAULI" \
      -p "dyn_qel_cc = 0" -p "dyn_qel_nc = 0" \
      -p "dyn_res_cc = 0" -p "dyn_res_nc = 0" \
      -p "dyn_dis_cc = 1" -p "dyn_dis_nc = 0" \
      -p "dyn_coh_cc = 0" -p "dyn_coh_nc = 0" \
      -p "dyn_mec_cc = 0" -p "dyn_mec_nc = 0" \
      -p "dyn_hyp_cc = 0" -p "dyn_lep = 0" \
      > gen.log 2>&1
)
# THE PROGRESS BAR IS 40 kB OF CARRIAGE RETURNS.  Kept only in digest form,
# the same treatment the GENIE logs get.
tr '\r' '\n' < "$d/gen.log" | grep -v "% of DIScc\|% of real events" > "$d/gen_digest.log"
mv -f "$d/gen_digest.log" "$d/gen.log"

"$HERE/nuwro2hepmc" "$d/events.root" "$d/events.hepmc" 2>/dev/null | tee "$d/convert.log"
# The ROOT file is an intermediate and nothing reads it once the HepMC exists
# (CONVENTIONS.md rule 1: intermediates do not outlive their results).
rm -f "$d/events.root"
n=$(grep -c '^E ' "$d/events.hepmc" 2>/dev/null); n=${n:-0}
echo "  $n events in $(basename "$d")/events.hepmc"
[ "$n" -eq "$NEV" ] || { echo "  !! wrote $n events, asked for $NEV" >&2; exit 1; }
