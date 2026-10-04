#!/usr/bin/env bash
# Reproduce the YADISM NLO reference with GENIE HEDIS's OWN parton densities.
#
# WHY.  GENIE's GHE19_00a tune builds its structure functions with APFEL from
# NNPDF31sx_nlo_as_0118_LHCb_nf_6 (its CommonParam.xml, param_set HEDIS-SF),
# while this benchmark's reference uses NNPDF40_nnlo_as_01180.  Running the
# SAME calculation with the SAME PDF is what turns "HEDIS is 5% low" into a
# number with a cause; analysis/genie_hedis_pdf.py collects the result.
#
# It builds the tungsten free-nucleon average of each set -- the reference is
# computed on (74 p + 110 n)/184, not on a proton -- and then convolves the
# cached PineAPPL grids, so nothing here re-runs YADISM.
#
# Usage: tools/genie_hedis_pdf.sh          (about ten minutes, all cached)
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
. "$REPO/config.sh"

HEDIS_PDF=NNPDF31sx_nlo_as_0118_LHCb_nf_6
BENCH_PDF=NNPDF40_nnlo_as_01180
# LHAPDF ids for the two averages; 341000 is the tracked NNPDF4.0 one, so the
# HEDIS set takes the next free block.
# The .dat is gitignored (rebuilt in seconds); the .info and pdfsets.index are
# tracked, so a clone has the id and not the grid -- test for the MEMBER.
if [ ! -f "$REPO/data/pdfs/lhapdf/${HEDIS_PDF}_W184free/${HEDIS_PDF}_W184free_0000.dat" ]; then
    python3 "$REPO/tools/make_isoscalar_pdf.py" "$HEDIS_PDF" --index 341900
fi

ENERGIES=400,700,1000,2000,4000
for P in "${HEDIS_PDF}_W184free" "${BENCH_PDF}_W184free"; do
    # ZM-VFNS across the ladder: the PDF effect, in the scheme the inclusive
    # benchmark is defined in.
    BENCH_SELECTION=q4w3 BENCH_TMC=3 python3 \
        "$REPO/analysis/mhou_sigma_fonll.py" --zm --energies "$ENERGIES" \
        --pdf "$P" 2>&1 | grep -E '\(1, 1\)' | sed 's/^/  /'
    # FONLL at the anchor: the scheme GENIE HEDIS itself runs (APFEL FONLL-B),
    # so the residual after the PDF is like-for-like.
    BENCH_SELECTION=q4w3 BENCH_TMC=3 python3 \
        "$REPO/analysis/mhou_sigma_fonll.py" --energies 1000 \
        --pdf "$P" 2>&1 | grep -E '\(1, 1\)' | sed 's/^/  /'
done

# The proton-level decomposition: which PDF difference the shift actually is.
for P in NNPDF40_nlo_as_01180 NNPDF31_nnlo_as_0118 "$HEDIS_PDF" \
         NNPDF31sx_nlonllx_as_0118_LHCb_nf_6 "$BENCH_PDF"; do
    BENCH_SELECTION=q4w3 BENCH_TMC=3 python3 \
        "$REPO/analysis/mhou_sigma_fonll.py" --zm --energies 1000 \
        --pdf "$P" 2>&1 | grep -E '\(1, 1\)' | sed 's/^/  /'
done

python3 "$REPO/analysis/genie_hedis_pdf.py"
