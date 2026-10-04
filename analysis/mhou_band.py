#!/usr/bin/env python3
"""The 7-point MHOU band of the matched calculation, ready to draw.

`analysis/mhou_hadron.py` produces it -- the POWHEG scale weights joined to
the SHOWERED events by `lhe_index`, so the band survives a cut on the
charged-hadron multiplicity, which is what the FASER$\\nu$ selection is.  This
module is the reader: it applies the same window and the same rebinning a
figure applies to its curves, and only then takes the envelope.

>>> MERGE FIRST, ENVELOPE SECOND. <<<  Two bins merge by ADDING their
contents at each of the seven scales; the band of the merged bin is the
envelope of those sums.  Enveloping first and averaging the ratios afterwards
is a different number wherever the band varies across the merged bins, and it
is the one a stored relative band would force -- which is why
`mhou_hadron.py` stores the seven histograms themselves.

>>> THE MUON BAND IS MEASURED ON THE POWHEG-V2 CROSS-VARIANT. <<<  POWHEG-RES
runs with `manyseeds` over a hundred-odd Les Houches files and is not
reweightable in one pass, so the muon column's band comes from the same code
as the neutrino column's, run on the muon neutral current.  Every figure that
draws it says so; `LABEL` below is the phrase they use.
"""
import json
import os

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# WHICH SAMPLE THE BAND IS MEASURED ON, per column.  Both are POWHEG-V2: the
# neutrino entry is the published one, the muon entry is the cross-variant.
SOURCE = {"results": "powheg_v2", "results_nu": "powheg_nu"}


# The current each column is, for the legend.  Both columns carry a band and
# the two are DIFFERENT NUMBERS -- the muon neutral current's is several times
# the neutrino charged current's -- so a legend with two entries reading
# "7-pt scale" would leave the reader unable to tell which panel each belongs
# to.  Every entry names its current.
CURRENT = {"results": r"\mu", "results_nu": r"\nu"}


def label(resdir, denom_key):
    """What the legend calls the band, given the curve it is drawn on.

    THE PROVENANCE IS NAMED WHENEVER IT DIFFERS FROM THE CURVE.  A band
    measured on POWHEG-V2 and drawn silently on a POWHEG-RES curve puts one
    code's number on another code's line -- the ambiguity CONVENTIONS.md rule 3
    forbids in prose, and no better on a figure.  Where the curve IS the
    sample the band was measured on, there is nothing to disambiguate and the
    current alone is enough.
    """
    cur = CURRENT.get(resdir, "")
    if denom_key == SOURCE.get(resdir):
        return rf"7-pt scale (${cur}$)"
    return rf"7-pt scale (${cur}$, POWHEG-V2)"


def path(resdir, selection, esuffix=""):
    return f"{BASE}/{resdir}/mhou_hadron_{selection}{esuffix}.json"


def load(resdir, selection, esuffix=""):
    p = path(resdir, selection, esuffix)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        d = json.load(f)
    sel = (d.get("selection") or {}).get("name")
    if sel != selection:
        raise SystemExit(f"{p}: selection is {sel!r}, not {selection!r}")
    return d


# HOW WELL THE BAND HAS TO BE RESOLVED BEFORE IT IS DRAWN AS ITS OWN BIN.
# The envelope must exceed its OWN Monte Carlo error by this factor; a bin
# that does not is merged with its neighbours until it does.  What such a bin
# would otherwise draw is sampling noise wearing the shape of a theory
# uncertainty, and it is not a small effect: POWHEG's reweighting factor can
# be enormous for an event near the edge of the Born phase space -- ONE event
# in the neutrino sample carries f = 156 with a negative nominal weight, and
# on its own it turns a one per cent band into a twenty-three per cent one in
# the bin that holds it.  The error stored by mhou_hadron.py sees that
# immediately, because a single event dominating a weighted mean makes the
# error as large as the mean's own displacement.
RESOLVE = 3.0
# and a bin is never merged past this many of the figure's own bins.  The cap
# matters: an outlier is diluted by merging, not removed, so without it a bin
# holding one f = 156 event would keep merging until it had swallowed the
# whole axis.  Four bins takes that spike from twenty-three per cent to two,
# which is the size of the honest band around it.
MAX_MERGE = 4


def _err(sw, swf, sw2, sw2f, sw2f2):
    """The Monte Carlo error on each varied-to-nominal ratio, per bin.

    r = S(w f) / S(w), and the events are the SAME at every scale, so this is
    the error on a weighted mean of f:
        Var(r) = [ S(w^2 f^2) - 2 r S(w^2 f) + r^2 S(w^2) ] / S(w)^2
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(sw != 0, swf / np.where(sw != 0, sw, 1.0), 1.0)
        var = (sw2f2 - 2.0 * r * sw2f + r * r * sw2) / np.where(sw != 0,
                                                                sw * sw, 1.0)
    return r, np.sqrt(np.maximum(var, 0.0))


def rel(resdir, obs, selection, esuffix="", trim_lo=None, rebin=None,
        resolve=RESOLVE):
    """(edges, rel_lo, rel_hi) for one observable, or None if absent.

    `trim_lo` and `rebin` are the figure's own, applied in the figure's own
    order -- drop the bins below the window, then merge -- so the band's
    edges land ON the curve's edges.  Adjacent bins are then merged further,
    left to right, until each one's envelope is larger than its own Monte
    Carlo error by the factor `resolve`; the returned edges are therefore a
    SUBSET of the figure's, never a different grid.
    """
    d = load(resdir, selection, esuffix)
    if d is None:
        return None
    h = (d.get("hists") or {}).get(obs)
    if not h or "curves" not in h or "sw2" not in h:
        return None
    e = np.asarray(h["edges"], dtype=float)
    c = np.asarray(h["curves"], dtype=float)          # (7, nbins)
    sw2 = np.asarray(h["sw2"], dtype=float)
    sw2f = np.asarray(h["sw2f"], dtype=float)
    sw2f2 = np.asarray(h["sw2f2"], dtype=float)
    if trim_lo is not None:
        k = int(np.searchsorted(e, trim_lo - 1e-9))
        e, c = e[k:], c[:, k:]
        sw2, sw2f, sw2f2 = sw2[k:], sw2f[:, k:], sw2f2[:, k:]

    def merge(n, arr):
        m = (arr.shape[-1] // n) * n
        return arr[..., :m].reshape(arr.shape[:-1] + (-1, n)).sum(axis=-1)

    if rebin and rebin > 1:
        m = (c.shape[1] // rebin) * rebin
        c, sw2, sw2f, sw2f2 = (merge(rebin, c), merge(rebin, sw2),
                               merge(rebin, sw2f), merge(rebin, sw2f2))
        e = e[:m + 1:rebin]

    # THE ADAPTIVE PASS.  Accumulate bins from the left; close the group as
    # soon as its envelope is resolved, and close whatever is left at the end
    # by folding it into the previous group rather than drawing it alone.
    groups, start = [], 0
    n = c.shape[1]
    while start < n:
        stop = start
        while stop < n:
            stop += 1
            sl = slice(start, stop)
            C = c[:, sl].sum(axis=1)
            r, err = _err(C[0], C, sw2[sl].sum(), sw2f[:, sl].sum(axis=1),
                          sw2f2[:, sl].sum(axis=1))
            env = float(np.max(np.abs(r - 1.0)))
            if C[0] != 0 and env > resolve * float(np.max(err)):
                break
            if stop - start >= MAX_MERGE:
                break
        groups.append((start, stop))
        start = stop
    if len(groups) > 1:
        # the trailing group may be unresolved; fold it back
        a, b = groups[-1]
        C = c[:, a:b].sum(axis=1)
        r, err = _err(C[0], C, sw2[a:b].sum(), sw2f[:, a:b].sum(axis=1),
                      sw2f2[:, a:b].sum(axis=1))
        if C[0] == 0 or float(np.max(np.abs(r - 1.0))) <= \
                resolve * float(np.max(err)):
            groups[-2] = (groups[-2][0], b)
            groups.pop()

    edges = np.array([e[a] for a, _b in groups] + [e[groups[-1][1]]])
    lo = np.zeros(len(groups))
    hi = np.zeros(len(groups))
    for i, (a, b) in enumerate(groups):
        C = c[:, a:b].sum(axis=1)
        if C[0] == 0:
            continue
        r = C / C[0] - 1.0
        lo[i], hi[i] = float(r.min()), float(r.max())
    return edges, lo, hi
