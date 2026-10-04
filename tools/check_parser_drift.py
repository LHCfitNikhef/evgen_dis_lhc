#!/usr/bin/env python3
"""Run BOTH HepMC parsers over the same file and require they agree.

WHY THIS EXISTS.  `parse_hepmc3` is written out twice -- 168 lines in
analyze.py and 142 in analyze_nu.py -- and the two were written independently.
It is the function that turns a HepMC record into the weight, the beams, the
final-state particles, the D mesons and the heavy-flavour tag, so every number
in this benchmark passes through one copy or the other, and the muon and
neutrino sides never pass through the same one.

CONVENTIONS.md rule 2b asks that whatever is done for one current is done for the
other.  Rule 2 says a result is not trustworthy because it looks right.  Both
apply here with force: if the copies drift, each side stays perfectly
self-consistent, both cross-sections stay plausible, and the comparison the
whole benchmark exists to make quietly stops being like-for-like.  That is not
hypothetical -- the DIS kinematics WERE duplicated the same way, byte for
byte, until 2026-08-27 (commit e2fc9c9).

WHY A TEST RATHER THAN A MERGE.  Merging the two is the better end state, but
it is surgery on the most safety-critical function in the project and would
have to be validated by re-running every sample.  This gets most of the safety
for none of the risk, and it keeps working as a regression test afterwards.

HOW IT COMPARES THEM.  A NEUTRINO sample can be read by both: analyze.py's
parser takes beam_pid and supports 14 (it is used for nu NC), analyze_nu.py's
hardcodes it.  Every field both parsers return is then compared BIT for BIT,
event by event -- they are reading the same text, so anything less than exact
agreement is a difference in the code and not in the input.

The reverse direction is not available: analyze_nu.py's parser cannot read a
muon sample, because the beam PDG is not a parameter there.  That asymmetry is
itself part of what a merge would fix.

Usage:
  tools/check_parser_drift.py [hepmc file] [--max N]
Exit status 1 on any disagreement.
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")

import analyze                                            # noqa: E402
import analyze_nu                                         # noqa: E402

# The kept  1 TeV POWHEG-V2 proton sample (the earlier default, powheg/nu_job_1,
# was deleted with earlier on 2026-09-19 and left this check failing with no arguments).
DEFAULT_SAMPLE = f"{BASE}/powheg/v2nu_p_job_1/events.hepmc"
DEFAULT_MAX = 20000


def describe(ev):
    """The comparable content of one yielded event, in a stable form."""
    w, k, P, parts, dmes, hard = ev
    return {
        "weight": w,
        "beam_lepton": k,
        "beam_proton": P,
        "n_particles": len(parts),
        "particles": sorted(parts),
        "dmesons": sorted(dmes),
        "hard_in": sorted(hard.get("in", [])),
        "hard_out": sorted(hard.get("out", [])),
        "fs_heavy": sorted(hard.get("fs_heavy", frozenset())),
        # the reweighting join key: a drift here silently mispairs weights
        "lhe_index": hard.get("lhe_index"),
    }


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    nmax = DEFAULT_MAX
    for a in sys.argv[1:]:
        if a.startswith("--max"):
            nmax = int(a.split("=", 1)[1]) if "=" in a else DEFAULT_MAX
    path = args[0] if args else DEFAULT_SAMPLE
    if not os.path.exists(path):
        sys.exit(f"no sample at {path}\n"
                 f"  pass a neutrino HepMC file -- both parsers must be able "
                 f"to read it, which\n  means a nu_mu beam (analyze_nu's "
                 f"parser hardcodes the beam PDG).")

    print(f"comparing both parse_hepmc3 implementations on\n  {path}")
    print(f"  (up to {nmax} events)")

    stats_mu, stats_nu = {}, {}
    gen_mu = analyze.parse_hepmc3(path, beam_pid=14, run_stats=stats_mu,
                                  hard_graph=True)
    gen_nu = analyze_nu.parse_hepmc3(path, run_stats=stats_nu,
                                     hard_graph=True)

    bad, n = [], 0
    while n < nmax:
        a = next(gen_mu, None)
        b = next(gen_nu, None)
        if a is None and b is None:
            break
        if a is None or b is None:
            bad.append(f"event {n}: one parser stopped and the other did not "
                       f"(analyze.py {'ended' if a is None else 'continued'})")
            break
        n += 1
        da, db = describe(a), describe(b)
        for key in da:
            if da[key] != db[key]:
                if len(bad) < 20:
                    bad.append(f"event {n}, field {key!r}:\n"
                               f"      analyze.py    {da[key]!r}\n"
                               f"      analyze_nu.py {db[key]!r}")

    for label, sa, sb in (("run_stats", stats_mu, stats_nu),):
        for key in sorted(set(sa) | set(sb)):
            if sa.get(key) != sb.get(key):
                bad.append(f"{label}[{key!r}]: analyze.py {sa.get(key)!r} "
                           f"vs analyze_nu.py {sb.get(key)!r}")

    if bad:
        print(f"\n*** PARSER DRIFT: {len(bad)} disagreement(s) over "
              f"{n} events ***\n", file=sys.stderr)
        for b in bad:
            print("  " + b, file=sys.stderr)
        print("\n  The two copies of parse_hepmc3 have diverged. Each side of "
              "the benchmark\n  stays self-consistent when this happens, so "
              "no cross-section will look wrong.",
              file=sys.stderr)
        return 1
    print(f"\n{n} events: the two parsers agree bit for bit on every field "
          f"they both return")
    print(f"  run_stats agree too: {sorted(stats_mu)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
