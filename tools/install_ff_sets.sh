#!/usr/bin/env bash
# Install the LHAPDF fragmentation-function sets of the report's
# "Fragmentation functions" tab (analysis/ff_comparison.py, 2026-09-21) into
# $LHAPDF_DATA_PATH as config.sh sets it.  The list is ff_comparison.SETS,
# the one source; a set already present is skipped.  ~0.5 GB of tarballs,
# ~2 GB unpacked.
#
# Usage: tools/install_ff_sets.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=$(cd "$HERE/.." && pwd)
. "$BENCH_REPO/config.sh" > /dev/null
DEST=${LHAPDF_DATA_PATH%%:*}
URL=https://lhapdfsets.web.cern.ch/current
SETS=$(cd "$BENCH_REPO/analysis" && python3 -c \
    'import ff_comparison as f; print(" ".join(n for v in f.SETS.values() for n, *_ in v))')
for s in $SETS; do
    if [ -f "$DEST/$s/$s.info" ]; then
        echo "have    $s"
        continue
    fi
    echo "install $s -> $DEST"
    curl -sfL "$URL/$s.tar.gz" | tar -xz -C "$DEST"
done
