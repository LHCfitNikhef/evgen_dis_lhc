#!/usr/bin/env bash
# No tracked file may name a location on one particular machine (CONVENTIONS.md
# rule 1): a home directory, a user's temporary directory, a macOS private
# path.  Everything outside the repository goes through config.sh, and a
# clone must run on any machine -- the 2026-09-30 audit found such paths in
# card comments, patch headers and a generated ROOT macro.
#
# What is allowed: config.sh's documented DEFAULTS (Homebrew's standard
# prefix) and the literal placeholder "/Users/..." in prose.  Binary files are
# not scanned.
#
#   tools/check_no_abs_paths.sh      # exit 1 if any tracked text file offends
set -u
cd "$(dirname "$0")/.."
hits=$(git grep -n -I -E '(/Users/[A-Za-z0-9_]|/home/[a-z][A-Za-z0-9_-]*/|/private/(tmp|var)/|/tmp/[A-Za-z])' \
        -- ':!tools/check_no_abs_paths.sh' | grep -v '/Users/\.\.\.' || true)
if [ -n "$hits" ]; then
    echo "$hits"
    echo
    echo "tracked files name a machine-specific path -- route it through config.sh"
    exit 1
fi
echo "no machine-specific absolute paths in tracked files"
