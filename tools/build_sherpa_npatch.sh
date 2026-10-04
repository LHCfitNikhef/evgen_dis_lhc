#!/usr/bin/env bash
# Build the PATCHED COPY of the Sherpa install used for neutron-beam runs.
#
# Usage: tools/build_sherpa_npatch.sh [--force]
#
# WHY (user decision 2026-09-13, "patched copy for neutron runs").  Stock
# Sherpa 3.0.5 segfaults on the first event of any 2112 (neutron) beam:
# Hadron_Remnant::RemnantFlavour builds the diquark kf 21101 / 201, which is
# not a particle.  patches/sherpa-hadron-remnant-neutron.diff is the one-line
# fix.  It is applied to a COPY so the shared install -- which every proton
# run of this benchmark uses, some of them running while this is built -- is
# never modified.
#
# WHAT IT DOES
#   1. copies $SHERPA_INSTALL to $SHERPA_NPATCH (config.sh, BENCH_SHERPA_NPATCH);
#   2. applies the patch to a scratch copy of REMNANTS/Main/Hadron_Remnant.C
#      (the source tree is not touched either), compiles it with the exact
#      compile flags CMake used for the install (read from the build tree's
#      flags.make), and relinks libRemnants from the build tree's other object
#      files with CMake's own link line (link.txt), into the COPY;
#   3. REWRITES EVERY LC_RPATH in the copy's binary and dylibs from the shared
#      lib dir to the copy's.  This is the step that matters: every Sherpa
#      dylib carries an absolute rpath to the install it was built in, so a
#      copied binary would otherwise load the UNPATCHED libRemnants from the
#      shared install -- silently, and with -flat_namespace possibly a mix of
#      both trees;
#   4. re-signs (ad hoc) everything it changed, which arm64 macOS requires.
#
# RUNNING THE COPY.  $SHERPA_NPATCH/bin/Sherpa, with SHERPA_LIBRARY_PATH and
# SHERPA_SHARE_PATH pointed at the copy (the compiled-in defaults name the
# shared install; Sherpa's Library_Loader dlopens by those paths).  Both are
# ordinary environment variables and survive caffeinate, unlike DYLD_*, which
# macOS strips when a SIP-protected binary sits in the exec chain.
# tools/sherpa_production.sh does this and VERIFIES it inside every neutron job
# (lsof on the running process: libRemnants from the copy, nothing from the
# shared lib dir).
#
# Portability: macOS only as written (install_name_tool, codesign).  On Linux
# the equivalent is patchelf --set-rpath; not needed yet.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH_REPO=${BENCH_REPO:-$(cd "$HERE/.." && pwd)}
export BENCH_REPO
. "$BENCH_REPO/config.sh"

SRC_INSTALL=$SHERPA_INSTALL
DST=$SHERPA_NPATCH
SRC_TREE=$(cd "$SRC_INSTALL/.." && pwd)          # the Sherpa source tree
BUILD=$SRC_TREE/build
PATCH=$BENCH_REPO/patches/sherpa-hadron-remnant-neutron.diff
LIBSUB=lib/SHERPA-MC
OLD_RPATH=$SRC_INSTALL/$LIBSUB
NEW_RPATH=$DST/$LIBSUB

[ "$(uname -s)" = Darwin ] || { echo "macOS only (install_name_tool/codesign)" >&2; exit 1; }
[ -x "$SRC_INSTALL/bin/Sherpa" ] || { echo "no Sherpa at $SRC_INSTALL" >&2; exit 1; }
[ -f "$PATCH" ] || { echo "no $PATCH" >&2; exit 1; }
OBJDIR=$BUILD/REMNANTS/Main/CMakeFiles/Remnants.dir
[ -f "$OBJDIR/flags.make" ] && [ -f "$OBJDIR/link.txt" ] || {
    echo "no CMake build tree at $BUILD -- cannot rebuild libRemnants with the install's flags" >&2; exit 1; }

if [ -e "$DST" ]; then
    if [ "${1:-}" = --force ]; then rm -rf "$DST"
    else echo "$DST exists; --force rebuilds it" >&2; exit 1; fi
fi

echo "== copying $SRC_INSTALL -> $DST"
mkdir -p "$(dirname "$DST")"
cp -Rp "$SRC_INSTALL" "$DST"

echo "== patching a scratch copy of Hadron_Remnant.C"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/REMNANTS/Main"
cp "$SRC_TREE/REMNANTS/Main/Hadron_Remnant.C" "$WORK/REMNANTS/Main/"
( cd "$WORK" && patch -p0 < "$PATCH" )
grep -q 'if (!taken && flav==flit) { taken = true; continue; }' \
    "$WORK/REMNANTS/Main/Hadron_Remnant.C" || { echo "patch did not apply" >&2; exit 1; }

echo "== compiling with the install's CMake flags"
flag () { sed -n "s/^$1 = //p" "$OBJDIR/flags.make" | head -1; }
CXX=$(awk '{print $1; exit}' "$OBJDIR/link.txt")
# word splitting of the flag strings is intended
# shellcheck disable=SC2046
$CXX $(flag CXX_DEFINES) $(flag CXX_INCLUDES) $(flag CXX_FLAGS) \
    -c "$WORK/REMNANTS/Main/Hadron_Remnant.C" -o "$WORK/Hadron_Remnant.C.o"

echo "== relinking libRemnants into the copy"
# CMake's link line, run from its own directory, with two substitutions: the
# output goes into the copy, and the patched object replaces the stock one.
LINK=$(cat "$OBJDIR/link.txt")
LINK=${LINK//CMakeFiles\/Remnants.dir\/Hadron_Remnant.C.o/$WORK\/Hadron_Remnant.C.o}
LINK=$(printf '%s' "$LINK" | sed -E "s#-o [^ ]*libRemnants\.0\.0\.0\.dylib#-o $NEW_RPATH/libRemnants.0.0.0.dylib#")
case "$LINK" in *"$WORK/Hadron_Remnant.C.o"*) ;; *) echo "link line: object substitution failed" >&2; exit 1 ;; esac
case "$LINK" in *"-o $NEW_RPATH/libRemnants.0.0.0.dylib"*) ;; *) echo "link line: output substitution failed" >&2; exit 1 ;; esac
rm -f "$NEW_RPATH/libRemnants.0.0.0.dylib"
( cd "$BUILD/REMNANTS/Main" && eval "$LINK" )

echo "== rewriting rpaths $OLD_RPATH -> $NEW_RPATH"
n=0
for f in "$DST/bin/Sherpa" "$NEW_RPATH"/*.dylib; do
    [ -L "$f" ] && continue
    if otool -l "$f" | grep -A2 LC_RPATH | grep -q "path $OLD_RPATH "; then
        install_name_tool -rpath "$OLD_RPATH" "$NEW_RPATH" "$f" 2>/dev/null
        n=$((n + 1))
    fi
    codesign -f -s - "$f" 2>/dev/null
done
echo "   rewrote $n rpaths"

echo "== checks"
bad=0
for f in "$DST/bin/Sherpa" "$NEW_RPATH"/*.dylib; do
    [ -L "$f" ] && continue
    if otool -l "$f" | grep -A2 LC_RPATH | grep -q "path $OLD_RPATH "; then
        echo "   still points at the shared install: $f"; bad=1; fi
done
if cmp -s "$SRC_INSTALL/$LIBSUB/libRemnants.0.0.0.dylib" "$NEW_RPATH/libRemnants.0.0.0.dylib"; then
    echo "   libRemnants in the copy is identical to the shared one -- not patched"; bad=1; fi
[ "$bad" = 0 ] || exit 1
{
  echo "Sherpa install copy with patches/sherpa-hadron-remnant-neutron.diff"
  echo "built $(date '+%Y-%m-%d %H:%M:%S') by tools/build_sherpa_npatch.sh from $SRC_INSTALL"
  echo "patch md5 $(md5 -q "$PATCH")"
  echo "libRemnants md5 $(md5 -q "$NEW_RPATH/libRemnants.0.0.0.dylib")"
} > "$DST/NPATCH_STAMP"
cat "$DST/NPATCH_STAMP"
echo "== done"
