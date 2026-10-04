#!/usr/bin/env python3
"""YADISM -> PineAPPL interpolation grids, and convolution against any PDF.

WHY.  The PDF-dependence study spends essentially all of its time convoluting:
the YADISM operator is built once per energy (~90-150 s, PDF-independent) and
then applied once per PDF member at ~9 s each.  With ~280 members over three
energies that is about seven hours, and every new PDF set costs another two.

A PineAPPL grid stores the same PDF-independent information in a form built for
exactly this: convolution becomes a weighted sum over interpolation nodes.
Measured here, on an identical calculation, a convolution drops from 225 ms to
4.9 ms -- 46x -- and the two agree to 1.2e-14, which is machine precision
rather than "close enough".

WHAT IT DOES NOT CHANGE.  The grid is a faster route to the SAME numbers, not a
different calculation: same theory card, same x/Q2 nodes, same observables.
Everything downstream -- the spline interpolation, the fiducial integration,
the charm fraction -- is untouched, because this replaces only the step that
turns an operator plus a PDF into structure-function values.

SCALE VARIATION COMES FREE.  `Grid.convolve` takes `xi=[(xiR, xiF)]`, so the
7-point envelope can be evaluated from the same grid rather than re-running
YADISM per scale point.  That is the MHOU study as well as the PDF one.

THE GRIDS ARE A CACHE, NOT A RESULT.  They are derived from the theory card and
the node grid and can always be rebuilt, so they live under
`analysis/pineappl_cache/` and are gitignored like the other caches here.  The
cache key includes the current, the beam energy, the perturbative order and a
hash of the node grid -- the last because a grid silently reused across beam
energies is exactly the bug the FONLL cache had, where a 1 TeV grid was served
for 400 GeV and produced plausible wrong numbers.

Usage:
  pineappl_grids.py build nu 1000 1        # build the NLO CC grids at 1 TeV
  pineappl_grids.py check nu 1000 1        # rebuild and verify against yadism
  pineappl_grids.py check nu 1000 1 FONLL-FFNS 3   # ... a FONLL component
  pineappl_grids.py build mu 700 2 FONLL-FFNS 3 charm   # charm observables only
  BENCH_TMC=3 pineappl_grids.py build ...    # the target-mass-corrected grids
"""
import hashlib
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "pineappl_cache")

# the structure functions the studies need, per flavour group
NAMES_TOTAL = ["F2_total", "FL_total", "F3_total"]
NAMES_CHARM = ["F2_charm", "FL_charm", "F3_charm"]


def _calc_module(current):
    """The YADISM driver for one current -- the only thing that differs."""
    return "yadism_cc_calc" if current == "nu" else "yadism_calc"


def grid_key(current, energy, pto, x_nodes, q2_nodes,
             fns="ZM-VFNS", nfff=5, tmc=None):
    """Cache key, including a hash of the NODES and the FLAVOUR SCHEME.

    The nodes are hashed rather than assumed: a grid built on one beam
    energy's nodes and served for another gives plausible, wrong numbers, and
    that has happened here before with the FONLL cache.

    THE SCHEME IS IN THE KEY FOR THE SAME REASON, and it is the more dangerous
    omission of the two.  FONLL is assembled from three runs that differ ONLY
    in (FNS, NfFF) -- ZM-VFNS/5, FONLL-FFNS/3, FONLL-FFN0/3 -- so with the
    scheme left out of the key all three would collide on one filename, the
    first one built would be served for the other two, and
    F = ZM + f_thr*(FFNS - FFN0) would evaluate to exactly the ZM answer.
    The study would then report "FONLL" while computing ZM-VFNS, with no
    error and a perfectly plausible plot.  Adding this invalidated the
    pre-existing ZM grids, which is the intended cost: they rebuild in ~2 min
    per energy and a cache is not a result.

    THE TARGET-MASS CORRECTION IS IN THE KEY (2026-09-13), for the reason the
    scheme is: TMC is a theory-card switch on the same observables on the same
    nodes, so a TMC-on grid and a TMC-off one would otherwise share a filename
    and the first one built would be served for both.  The PDF is NOT in the
    key and must not be -- a grid is PDF-independent, which is its point.
    Mode 0 adds nothing, so every pre-existing grid keeps its name.
    """
    import target
    h = hashlib.sha1(
        np.concatenate([np.asarray(x_nodes, dtype=float),
                        np.asarray(q2_nodes, dtype=float)]).tobytes()
    ).hexdigest()[:12]
    tag = f"{fns.replace('-', '').lower()}{nfff}"
    m = target.tmc() if tmc is None else int(tmc)
    return (f"{current}_{energy:g}_pto{pto}_{tag}_{h}"
            + ("" if m == 0 else f"_tmc{m}"))


def grid_file(key, name):
    return os.path.join(CACHE, f"{key}_{name}.pineappl.lz4")


def build(current, energy, pto, names, x_nodes=None, q2_nodes=None,
          force=False, fns="ZM-VFNS", nfff=5):
    """Build (or reuse) one PineAPPL grid per observable.  Returns the key.

    `fns`/`nfff` override the theory card's flavour scheme, which is what the
    FONLL assembly needs: it wants the SAME observables on the SAME nodes in
    three different schemes.  They go into the cache key -- see grid_key().
    """
    import importlib
    os.environ["BENCH_ENERGY"] = f"{energy:g}"
    for m in ("yadism_calc", "yadism_cc_calc", "beams", "analyze",
              "runmeta", "selection"):
        sys.modules.pop(m, None)
    calc = importlib.import_module(_calc_module(current))
    from yadism_calc import sf_grid_nodes
    if x_nodes is None:
        x_nodes, q2_nodes = sf_grid_nodes()

    key = grid_key(current, energy, pto, x_nodes, q2_nodes, fns, nfff)
    missing = [n for n in names if not os.path.exists(grid_file(key, n))]
    if not missing and not force:
        print(f"[grids] {key}: all {len(names)} grid(s) cached")
        return key, x_nodes, q2_nodes

    import yadism
    from yadbox.export import dump_pineappl_to_file
    theory, obs = calc.make_cards(pto, x_nodes, q2_nodes)
    theory["FNS"], theory["NfFF"] = fns, nfff
    # the card must carry the TMC mode the key names.  Both read target.tmc()
    # today; a card builder that stopped doing so would build a TMC-off grid
    # under a TMC-on name, so it is checked rather than trusted.
    import target
    if theory.get("TMC", 0) != target.tmc():
        raise SystemExit(f"[grids] card TMC {theory.get('TMC')} != "
                         f"target.tmc() {target.tmc()} for {key}")
    kins = obs["observables"][list(obs["observables"])[0]]
    obs["observables"] = {n: kins for n in names}
    print(f"[grids] {key}: running yadism ({fns}, NfFF={nfff}) for "
          f"{len(names)} observable(s) ...", flush=True)
    out = yadism.run_yadism(theory, obs)
    os.makedirs(CACHE, exist_ok=True)
    for n in names:
        f = grid_file(key, n)
        # WRITE THEN RENAME.  Two builds can legitimately be in flight for the
        # same (current, energy, scheme) -- the prebuild driver runs a batch in
        # parallel, and a `check` run rebuilds with force -- and writing the
        # final path directly lets them interleave into a torn .lz4 that fails
        # to load, or worse, loads short.  os.replace is atomic on one
        # filesystem, so the loser of a race simply overwrites with identical
        # content and no reader ever sees a partial file.
        tmp = f"{f}.tmp{os.getpid()}"
        dump_pineappl_to_file(out, tmp, n)
        os.replace(tmp, f)
        print(f"[grids]   wrote {os.path.basename(f)} "
              f"({os.path.getsize(f)/1e3:.0f} kB)")
    return key, x_nodes, q2_nodes


def convolve(key, name, pdf, xir=1.0, xif=1.0):
    """Structure-function values for one observable, as yadism would give them.

    Returns a plain array in the grid's bin order, which is the same order
    `Output[name]` uses -- so a caller can substitute this for
    `[p["result"] for p in res[name]]` without changing anything else.
    """
    from pineappl.grid import Grid
    f = grid_file(key, name)
    if not os.path.exists(f):
        raise FileNotFoundError(f"no grid at {f} -- build it first")
    g = Grid.read(f)
    # THREE entries, not two: this PineAPPL wants (xiR, xiF, xiA), the last
    # being the fragmentation scale, which DIS has no use for and which stays
    # at 1.  Passing a 2-tuple raises "expected tuple of length 3" -- an
    # honest error rather than a silent misreading, but only because the
    # binding checks; it is worth pinning down here rather than at each call.
    xi = None if (xir == 1.0 and xif == 1.0) else [(xir, xif, 1.0)]
    vals = g.convolve(pdg_convs=g.convolutions,
                      xfxs=[lambda pid, x, q2: pdf.xfxQ2(pid, x, q2)],
                      alphas=lambda q2: pdf.alphasQ2(q2),
                      xi=xi)
    return np.asarray(vals)


def res_like(key, names, pdf, xir=1.0, xif=1.0):
    """Grid results shaped exactly like `yadism.Output.apply_pdf...` returns.

    `{name: [{"result": v}, ...]}` -- the shape `splines_from()` consumes -- so
    a caller swaps the SOURCE of the numbers and changes nothing else.  That is
    deliberate: the fiducial integration, the splines and the charm fraction
    stay on the code path they were validated on.
    """
    return {n: [{"result": float(v)} for v in convolve(key, n, pdf, xir, xif)]
            for n in names}


def check(current, energy, pto, fns="ZM-VFNS", nfff=5):
    """Rebuild and verify the grids against yadism's own convolution.

    THE POINT OF THE WHOLE EXERCISE is that the grid is the same calculation,
    so that has to be measured rather than asserted.  A grid that is merely
    close would be worse than no grid: it would shift every number on the page
    by a little, consistently, and nothing would look wrong.
    """
    import importlib
    import lhapdf
    lhapdf.setVerbosity(0)
    names = NAMES_TOTAL + NAMES_CHARM
    key, x_nodes, q2_nodes = build(current, energy, pto, names, force=True,
                                   fns=fns, nfff=nfff)

    calc = importlib.import_module(_calc_module(current))
    import yadism
    theory, obs = calc.make_cards(pto, x_nodes, q2_nodes)
    theory["FNS"], theory["NfFF"] = fns, nfff
    kins = obs["observables"][list(obs["observables"])[0]]
    obs["observables"] = {n: kins for n in names}
    out = yadism.run_yadism(theory, obs)

    pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
    res = out.apply_pdf_alphas_alphaqed_xir_xif(
        pdf, lambda muR: pdf.alphasQ(muR), lambda muR: calc.ALPHA, 1.0, 1.0)
    worst = 0.0
    for n in names:
        yv = np.array([p["result"] for p in res[n]])
        gv = convolve(key, n, pdf)
        m = np.abs(yv) > 1e-12
        rel = np.abs(gv[m] / yv[m] - 1.0)
        worst = max(worst, float(rel.max()) if rel.size else 0.0)
        print(f"[check] {n:10s} max |grid/yadism - 1| = "
              f"{rel.max() if rel.size else 0:.3e}")
    print(f"[check] WORST over all observables: {worst:.3e}")
    if worst > 1e-9:
        sys.exit("[check] FAILED: the grid is not reproducing yadism.\n"
                 "  A grid that is merely close is worse than none -- it would "
                 "shift every\n  number consistently and look fine.")
    print("[check] the grid reproduces yadism to machine precision")


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__.strip().splitlines()[-2])
    mode, current, energy = sys.argv[1], sys.argv[2], float(sys.argv[3])
    pto = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    if current not in ("mu", "nu"):
        sys.exit("current must be mu or nu")
    # optional 6th/7th argument: the flavour scheme, so the FONLL component
    # grids can be verified against yadism exactly like the ZM ones are
    fns = sys.argv[5] if len(sys.argv) > 5 else "ZM-VFNS"
    nfff = int(sys.argv[6]) if len(sys.argv) > 6 else (5 if fns == "ZM-VFNS"
                                                       else 3)
    # optional 8th argument "charm": build the charm observables only.  The
    # massive FONLL components are only ever asked for charm (the inclusive
    # correction IS the charm one -- pdf_dependence.fonll_combine), and a
    # massive run over the total structure functions costs more for nothing.
    groups = sys.argv[7] if len(sys.argv) > 7 else "all"
    if groups not in ("all", "charm"):
        sys.exit("the 8th argument, if given, must be all or charm")
    if mode == "build":
        build(current, energy, pto,
              NAMES_CHARM if groups == "charm" else NAMES_TOTAL + NAMES_CHARM,
              fns=fns, nfff=nfff)
    elif mode == "check":
        check(current, energy, pto, fns, nfff)
    else:
        sys.exit("mode must be build or check")


if __name__ == "__main__":
    main()
