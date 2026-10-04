#!/usr/bin/env bash
# The sigma_fid(E) ladder for the arXiv:2506.13889 DIS region (muon NC only).
#
# ONE LINE, because the ladder itself is now tools/yadism_ladder_region.sh --
# there are two foreign regions and there was one copy of this loop per region,
# which is how the two would drift.  The region's own caveats live in that
# script's header and in analysis/selection.py.
exec "$(dirname "$0")/yadism_ladder_region.sh" dis2506 mu "$@"
