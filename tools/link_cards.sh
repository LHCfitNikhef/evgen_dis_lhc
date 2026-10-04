#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Re-create the symlinks from the generator install trees into this repo.
#
# The benchmark's generator cards (Sherpa.yaml, powheg.input) are VERSIONED
# HERE and only here.  The install trees under physics24/software (and
# physics24/faser/pythia8) hold symlinks pointing back, so there is exactly
# ONE copy of every card: what you ran is what is committed, always.
#
# Run this after reinstalling/rebuilding a generator, or after moving a run
# directory, either of which silently removes the symlink.  It is idempotent
# and REFUSES to overwrite a real file (that would destroy an edit made in
# the install tree) -- it reports those instead so you can merge by hand.
#
#   ./tools/link_cards.sh          # create/repair links
#   ./tools/link_cards.sh --check  # report only, change nothing
#
# Written for BASH 3.2, which is what /bin/bash is on macOS: no associative
# arrays, no ${var:A}/${var:h}/${var:t}, no glob qualifiers.  It used to be a
# zsh script; the cluster nodes are not guaranteed to have zsh.
# ---------------------------------------------------------------------------
shopt -s nullglob

# Resolve a path to its physical location, following symlinks in the final
# component too -- this is zsh's ${var:A}, which bash has no equivalent for.
resolve () {
  local p=$1 t n=0
  while [ -L "$p" ] && [ "$n" -lt 32 ]; do
    t=$(readlink "$p")
    case $t in
      /*) p=$t ;;
       *) p=$(dirname "$p")/$t ;;
    esac
    n=$((n + 1))
  done
  printf '%s/%s\n' "$(cd "$(dirname "$p")" 2>/dev/null && pwd -P)" "$(basename "$p")"
}

SELF=$(resolve "$0")
REPO=$(cd "$(dirname "$SELF")/.." && pwd)
CHECK=0
[ "${1-}" = "--check" ] && CHECK=1

BENCH_REPO=${BENCH_REPO:-$REPO}
. "$REPO/config.sh"   # SHERPA_RUNS, POWHEG_RES, POWHEG_V2

# ---------------------------------------------------------------------------
# TWO DISTINCT POWHEG CODES -- always name them explicitly, never just "POWHEG".
#
#   POWHEG-RES  powheg/powheg-dis-main/DIS_v   (POWHEG-BOX-RES)
#               inclusive AND charm production, in both CC and NC DIS;
#               mass effects NEGLECTED.  Supplies the muon NC NLO entry.
#
#   POWHEG-V2   powheg-cmass/nu-DIS-master     (POWHEG-BOX-V2, arXiv:2406.05115)
#               neutrino DIS; its purpose is MASSIVE charm production
#               (qmass > 0, CC s->c only).  Cards live under POWHEG-V2/.
#
# Repo card paths mirror this split: powheg/cards/POWHEG-{RES,V2}/<name>/.
# ---------------------------------------------------------------------------

# "<repo card dir>|<install run dir>", one per line.  A plain list rather than
# an associative array, which bash 3.2 does not have.
POWHEG_CARDS="\
POWHEG-RES/mu1TeV-wide|$POWHEG_RES/parallel-mu1TeV-wide
POWHEG-V2/nu1TeV-cmass|$POWHEG_V2/testrun-nu1TeV-cmass
POWHEG-V2/nu1TeV-prod|$POWHEG_V2/prod-nu1TeV
POWHEG-V2/mu1TeV-NC|$POWHEG_V2/prod-mu1TeV
POWHEG-RES/nu1TeV-CC|$POWHEG_RES/parallel-nu1TeV-CC"

nlink=0; nok=0; nconflict=0; nmissing=0

link_one () {  # $1 = repo file (source of truth), $2 = install-tree path
  local src=$1 dst=$2
  if [ ! -f "$src" ]; then
    printf '%s\n' "  MISSING IN REPO  $src"; nmissing=$((nmissing + 1)); return
  fi
  if [ -L "$dst" ]; then
    if [ "$(resolve "$dst")" = "$(resolve "$src")" ]; then nok=$((nok + 1)); return; fi
    printf '%s\n' "  RELINK           $dst"
    [ "$CHECK" -eq 1 ] || { rm -f "$dst"; ln -s "$src" "$dst"; }
    nlink=$((nlink + 1)); return
  fi
  if [ -e "$dst" ]; then
    # a REAL file here means someone edited the install tree directly
    if diff -q "$src" "$dst" >/dev/null 2>&1; then
      printf '%s\n' "  replace copy     $dst (identical to repo)"
      [ "$CHECK" -eq 1 ] || { rm -f "$dst"; ln -s "$src" "$dst"; }
      nlink=$((nlink + 1))
    else
      printf '%s\n' "  !! CONFLICT      $dst differs from repo copy - merge by hand"
      nconflict=$((nconflict + 1))
    fi
    return
  fi
  printf '%s\n' "  create           $dst"
  [ "$CHECK" -eq 1 ] || { mkdir -p "$(dirname "$dst")"; ln -s "$src" "$dst"; }
  nlink=$((nlink + 1))
}

printf '%s\n' "Sherpa cards:"
for d in "$REPO"/sherpa/Runs/*; do
  [ -d "$d" ] || continue                      # was the zsh (/) glob qualifier
  name=$(basename "$d")
  for f in "$d"/*; do
    [ -f "$f" ] || continue
    link_one "$f" "$SHERPA_RUNS/$name/$(basename "$f")"
  done
done

printf '%s\n' "POWHEG cards:"
# Fed by here-doc, NOT a pipe: a pipe would run the loop in a subshell and the
# counters would be lost.
while IFS='|' read -r name dir; do
  [ -n "$name" ] || continue
  link_one "$REPO/powheg/cards/$name/powheg.input" "$dir/powheg.input"
done <<EOF
$POWHEG_CARDS
EOF

printf '\n'
printf '%s\n' "linked/repaired: $nlink   already correct: $nok   conflicts: $nconflict   missing in repo: $nmissing"
[ "$CHECK" -eq 1 ] && printf '%s\n' "(--check: nothing was changed)"
[ "$nconflict" -gt 0 ] && exit 1
exit 0
