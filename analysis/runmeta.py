#!/usr/bin/env python3
"""The (selection, energy) stamp every result carries.

ONE definition, called by every script that writes a results/*.json --
analyze.py, analyze_nu.py, me_nlo_samples.py, the yadism_* calculators,
herwig7/make_histos*.py and mg5_lhe_histos.py.  Without it a JSON does not say
which fiducial region or which beam produced it, and the energy scan makes that
a real hazard: two files whose names differ by a suffix would otherwise be
indistinguishable once loaded, and a plot could combine 400 GeV with 4 TeV, or
the inclusive region with the FASER one, and look entirely reasonable.

    from runmeta import stamp
    out = {...}
    stamp(out)                    # muon, the configured energy
    stamp(out, current="nu")
"""
import json
import os
import tempfile

import beams
import selection


def energy():
    """The configured beam energy in GeV -- $BENCH_ENERGY, else the anchor."""
    return float(os.environ.get("BENCH_ENERGY") or beams.ANCHOR_ENERGY)


def stamp(out, current="mu", sel=None, e_lab=None):
    """Add the selection and beam fields to a result dict, and return it."""
    sel = sel or selection.get()
    e_lab = energy() if e_lab is None else float(e_lab)
    b = beams.Beams(current, e_lab)
    out["selection"] = sel.as_dict()
    out["energy_gev"] = e_lab
    out["beam_tag"] = b.tag
    out["two_kP"] = b.two_kP
    return out


def _plain(o):
    """numpy scalars -> Python scalars, recursively.

    json.dump cannot serialise numpy types and fails PART WAY THROUGH, leaving
    a truncated file behind.  That is exactly how results/histos_yadism*.json
    came to be 300-byte fragments ending at `"lo_closure_ok": ` -- the value
    was an np.bool_ from `abs(closure - 1) <= tol`.  Coerce before writing.
    """
    if isinstance(o, dict):
        return {k: _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    if hasattr(o, "item") and getattr(o, "shape", None) == ():
        return o.item()                      # numpy scalar
    if hasattr(o, "tolist") and hasattr(o, "shape"):
        return o.tolist()                    # numpy array
    return o


def write_result(path, out):
    """Write a result JSON atomically, coercing numpy types first.

    Atomic because a crash mid-dump previously left a truncated JSON that was
    then COMMITTED and only discovered when a later run tried to read it back.
    Writing to a temp file and renaming means a reader sees either the old
    file or the complete new one, never a fragment.
    """
    out = _plain(out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(out, f)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return path
