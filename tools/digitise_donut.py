#!/usr/bin/env python3
"""Read DONUT's Fig. 2 and Fig. 9 off arXiv:0711.0728 -- EXACTLY, not by eye.

WHAT IS TAKEN, AND WHY THOSE TWO FIGURES (user, 2026-09-10).  DONUT (E872) is
the experiment that first observed the tau neutrino, in a beam-dump beam on an
emulsion-and-steel target at Fermilab.  Two of its figures are what a
generator benchmark can be compared against:

  * FIG. 2, the energy spectra of the neutrinos INTERACTING in the target, per
    flavour.  Not a flux: it already carries the cross-section and the
    detector acceptance, so it is exactly the energy weight our prediction
    needs and nothing has to be folded with a cross-section twice.
  * FIG. 9, the charged-particle multiplicity n_ch at the primary vertex of
    all 578 located events, data and DONUT's own LEPTO-based Monte Carlo.

>>> THE FIGURES ARE VECTOR GRAPHICS, SO THIS IS AN EXTRACTION AND NOT A
DIGITISATION. <<<  The arXiv PDF carries the ROOT output as paths, so the
histogram outlines and the marker centres come out as the numbers that were
plotted, to the precision the PDF stores them (six decimals).  There is no
pixel step, no reading of a log axis by eye, and no operator judgement.  What
IS judgement -- which path is which curve -- is settled by the dash pattern
and then CHECKED against numbers published in the paper's text and tables.

THE THREE CHECKS, all of which must pass:

  1. FIG. 9's data points sum to 578, the located-event count in Section VII.
  2. FIG. 9's Monte Carlo sums to the same, since it is normalised to the data.
  3. FIG. 2's flavour composition agrees with TABLE III, which gives the
     Monte Carlo fractions of the 578 located events: nu_e CC 0.181, nu_mu CC
     0.199 + 0.159, nu_tau CC 0.018.  Table III's fourth category, "NCeff", is
     a CLASSIFICATION and not a flavour, so the comparison is between flavours
     of the CC part.

     >>> AND THE TWO ARE NOT EXPECTED TO BE EQUAL. <<<  Fig. 2 counts
     INTERACTIONS; Table III counts LOCATED events, and the location
     efficiency is 0.31 for the shower-like "LP" events against 0.77 for the
     rest (Section VII D).  Electron showers make nu_e CC the flavour most
     often lost, so the located sample must carry a SMALLER nu_e fraction than
     the interacting one.  Measured here: 38.4% against 32.5%, in that
     direction.  The check therefore allows eight points and would still catch
     the failure it is for -- two curves swapped, which moves nu_mu and nu_e by
     twenty.

WHY THE AXIS CALIBRATION NEEDS NO TICK LABEL.  Every tick is a path segment
with its own coordinates, and the major ticks are longer than the minor ones,
so the axis transformation is read off the geometry.  The one thing that has
to be identified is which major tick carries which value, and that is fixed by
the labels' glyph positions -- checked, in Fig. 2, by requiring the three
decades to be equally spaced to better than half a per cent.

Usage:
  tools/digitise_donut.py [path/to/0711.0728.pdf]
The PDF is downloaded to a scratch file if not given.  Writes
data/donut/flux.json and data/donut/nch.json, which ARE tracked -- they are a
few kilobytes and a clone must not need the paper to build the report.
"""
import json
import math
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
OUT = os.path.join(BASE, "data", "donut")
PDF_URL = "https://arxiv.org/pdf/0711.0728"

# Page numbers in the arXiv PDF (1-based), not the printed page numbers.
PAGE_FLUX, PAGE_NCH = 24, 31

# Published numbers this extraction is checked against.
N_LOCATED = 578                      # Section VII C: 578 of 866 located
TABLE_III = {"nu_e": 0.181, "nu_mu": 0.199 + 0.159, "nu_tau": 0.018}

NUM = re.compile(r"-?\d+\.?\d*")


def die(msg):
    sys.exit(f"digitise_donut: {msg}")


def page_svg(pdf, page):
    """One page of the PDF as SVG text."""
    tmp = os.path.join(tempfile.gettempdir(), f"donut_p{page}.svg")
    r = subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page),
                        pdf, tmp], capture_output=True, text=True)
    if r.returncode or not os.path.exists(tmp):
        die(f"pdftocairo failed on page {page}: {r.stderr.strip()}\n"
            f"  (poppler's pdftocairo is required; brew install poppler)")
    return open(tmp).read()


def paths(svg):
    """[(attributes, [(cmd, x, y), ...]), ...] for every <path> on the page."""
    out = []
    # THE WHOLE ELEMENT, not the part before `d=`.  poppler writes `transform`
    # AFTER the path data on some elements and before it on others, and a
    # transform that is not seen is a marker silently placed in the wrong
    # coordinate frame -- which lands it at a perfectly plausible height.
    for el in re.findall(r"<path\b[^>]*/?>", svg):
        m = re.search(r'd="([^"]*)"', el)
        if not m:
            continue
        attr, d = el, m.group(1)
        toks = (d.replace("M", " M ").replace("L", " L ")
                 .replace("C", " C ").replace("Z", " Z ").split())
        pts, i, cmd = [], 0, None
        while i < len(toks):
            t = toks[i]
            if t in "MLCZ":
                cmd, i = t, i + 1
                continue
            try:
                pts.append((cmd, float(toks[i]), float(toks[i + 1])))
            except (ValueError, IndexError):
                break
            i += 2
        out.append((attr, pts))
    return out


def to_frame(attr, x, y):
    """Page coordinates -> the frame the axes are drawn in.

    pdftocairo emits the plot under `matrix(s, 0, 0, -s, tx, ty)`, so the
    frame's y grows UPWARD while the page's grows downward.  Elements that
    already carry that transform are in frame coordinates; the markers of
    Fig. 9 do not carry it and have to be mapped.  Getting this backwards
    silently mirrors the figure, which is why it is done from the attribute
    rather than assumed.
    """
    m = re.search(r"matrix\(([-0-9.]+), 0, 0, ([-0-9.]+), ([-0-9.]+), "
                  r"([-0-9.]+)\)", attr)
    if m:
        return x, y                          # already in frame coordinates
    s, tx, ty = 0.83023, 54.0, 709.616       # the one transform this file uses
    return (x - tx) / s, (ty - y) / s


def ticks(pts, want_vertical):
    """(position, length) of every tick in an axis path."""
    out = []
    for i in range(len(pts) - 1):
        (c1, x1, y1), (c2, x2, y2) = pts[i], pts[i + 1]
        if c1 != "M" or c2 != "L":
            continue
        if want_vertical and abs(x1 - x2) < 1e-6:
            out.append((x1, abs(y2 - y1)))
        elif not want_vertical and abs(y1 - y2) < 1e-6:
            out.append((y1, abs(x2 - x1)))
    return out


def frame_of(ps):
    """The plot frame as (xmin, xmax, ymin, ymax), from the widest short path.

    ROOT draws the frame as a four- or five-point polyline, so it is the only
    short path that spans the whole plot in both directions.
    """
    best = None
    for _attr, pts in ps:
        if not pts or len(pts) > 6:
            continue
        xs = [x for _, x, _ in pts]
        ys = [y for _, _, y in pts]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if w > 300 and h > 200 and (best is None or w * h > best[0]):
            best = (w * h, min(xs), max(xs), min(ys), max(ys))
    if best is None:
        die("could not find the plot frame")
    return best[1:]


def major_ticks(ps, frame, vertical):
    """Every major tick on one axis, gathered from the paths that ARE axes.

    AN AXIS PATH IS RECOGNISED BY ITS SHAPE, not by its length or its order in
    the file: it is a run of segments, all perpendicular to the axis, all
    starting exactly on the frame edge.  A data curve fails that test even
    when it touches the axis, and poppler is free to split one axis over
    several paths or to draw a histogram with more segments than the axis has.

    MAJOR AND MINOR ARE THEN SEPARATED BY LENGTH, which inside an axis path is
    a clean two-valued distribution -- ROOT draws minors at half the major.
    The split is taken at 1.4x the shortest rather than at a written
    threshold, so it needs no maintenance per figure.
    """
    x0, x1, y0, y1 = frame
    ticks_found = []
    for _attr, pts in ps:
        seg = []
        for i in range(len(pts) - 1):
            (c1, xa, ya), (c2, xb, yb) = pts[i], pts[i + 1]
            if c1 != "M" or c2 != "L":
                continue
            if vertical and abs(xa - xb) < 1e-6:
                on = abs(min(ya, yb) - y0) < 1e-3 or abs(max(ya, yb) - y1) < 1e-3
                seg.append((round(xa, 3), abs(ya - yb), on))
            elif (not vertical) and abs(ya - yb) < 1e-6:
                on = abs(min(xa, xb) - x0) < 1e-3 or abs(max(xa, xb) - x1) < 1e-3
                seg.append((round(ya, 3), abs(xa - xb), on))
        if len(seg) >= 5 and all(s[2] for s in seg):
            ticks_found.extend(seg)
    if not ticks_found:
        die("no axis path found: none has five or more segments starting on "
            "the frame edge")
    shortest = min(ln for _p, ln, _o in ticks_found)
    return sorted(set(p for p, ln, _o in ticks_found if ln > 1.4 * shortest))


def read_flux(svg):
    """The three interacting-neutrino spectra of FIG. 2, per 6 GeV bin."""
    ps = paths(svg)
    fr = frame_of(ps)
    xmaj = major_ticks(ps, fr, True)
    if len(xmaj) != 7:
        die(f"FIG. 2: found {len(xmaj)} major x ticks, expected 7 "
            f"(0 to 300 GeV in steps of 50)")
    x0, x1 = xmaj[0], xmaj[-1]

    # the y axis is logarithmic: two major ticks plus the frame bottom, which
    # carries the "0.1" label.  Their spacing must be one decade each.
    ymaj = major_ticks(ps, fr, False)
    if len(ymaj) < 2:
        die("FIG. 2: could not find the y major ticks")
    decades = sorted([fr[2]] + ymaj[:2])
    d1, d2 = decades[1] - decades[0], decades[2] - decades[1]
    if abs(d1 / d2 - 1.0) > 0.005:
        die(f"FIG. 2: the two decades of the log axis measure {d1:.3f} and "
            f"{d2:.3f} device units -- they must be equal, so the axis has "
            f"been misidentified")
    dec = 0.5 * (d1 + d2)
    v_of_y = lambda y: 0.1 * 10.0 ** ((y - decades[0]) / dec)  # noqa: E731

    # the three curves, told apart by their dash pattern, exactly as the
    # figure's own legend does: dotted nu_mu, dashed nu_e, solid nu_tau
    def style(attr):
        m = re.search(r'stroke-dasharray="([^"]*)"', attr)
        return m.group(1) if m else "solid"
    curves = {}
    for attr, pts in ps:
        if len(pts) < 80 or "stroke-width=\"1.5\"" not in attr:
            continue
        xs = [x for _, x, _ in pts]
        if max(xs) - min(xs) < 400:
            continue                              # a legend line, not a curve
        curves.setdefault(style(attr), pts)
    name = {"1 2": "nu_mu", "3 3": "nu_e", "solid": "nu_tau"}
    if set(curves) != set(name):
        die(f"FIG. 2: found curve styles {sorted(curves)}, expected "
            f"{sorted(name)} (dotted nu_mu, dashed nu_e, solid nu_tau)")

    edges = [6.0 * i for i in range(51)]
    out = {}
    for st, pts in curves.items():
        # a ROOT histogram outline is a staircase: read the height of the
        # horizontal segment covering each bin centre
        segs = []
        for i in range(len(pts) - 1):
            (c1, xa, ya), (c2, xb, yb) = pts[i], pts[i + 1]
            if c1 in "ML" and c2 == "L" and abs(ya - yb) < 1e-6:
                segs.append((min(xa, xb), max(xa, xb), ya))
        # THE FRAME IS APPENDED TO EACH CURVE'S PATH by ROOT, so its top edge
        # is a horizontal segment spanning every bin.  Taking the highest
        # segment at a bin centre therefore returns the top of the plot for
        # all fifty bins -- a perfectly smooth, perfectly wrong spectrum.  A
        # staircase step is at most a few bins wide, so the narrowest segment
        # covering the centre is the one that belongs to the bin.
        bw = (x1 - x0) * 6.0 / 300.0
        vals = []
        for i in range(50):
            xc = x0 + (x1 - x0) * (edges[i] + 3.0) / 300.0
            hit = [(b - a, y) for a, b, y in segs
                   if a - 1e-6 <= xc <= b + 1e-6 and (b - a) < 4.0 * bw]
            vals.append(v_of_y(min(hit)[1]) if hit else 0.0)
        out[name[st]] = vals
    return {"e_edges_gev": edges, "spectra": out}


def read_nch(svg):
    """FIG. 9: the data and Monte Carlo n_ch distributions."""
    ps = paths(svg)
    fr = frame_of(ps)
    xmaj = major_ticks(ps, fr, True)
    ymaj = major_ticks(ps, fr, False)
    if len(xmaj) != 11 or len(ymaj) != 7:
        die(f"FIG. 9: found {len(xmaj)} x and {len(ymaj)} y major ticks, "
            f"expected 11 (0..20 by 2) and 7 (0..120 by 20)")
    nch_of_x = lambda x: 20.0 * (x - xmaj[0]) / (xmaj[-1] - xmaj[0])  # noqa
    n_of_y = lambda y: 120.0 * (y - ymaj[0]) / (ymaj[-1] - ymaj[0])   # noqa

    # A MARKER IS A CIRCLE, so it is drawn with Bezier segments and its
    # bounding box is square.  Error bars are straight lines with straight
    # caps and fail both tests -- which matters, because a cap read as a
    # marker lands at the end of an error bar and so at a plausible height.
    filled, open_ = [], []
    for attr, pts in ps:
        if sum(1 for c, _x, _y in pts if c == "C") < 4:
            continue
        loc = [to_frame(attr, x, y) for _, x, y in pts]
        xs = [q[0] for q in loc]
        ys = [q[1] for q in loc]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if not (3.0 < w < 12.0 and abs(w - h) < 0.3):
            continue
        c = (0.5 * (max(xs) + min(xs)), 0.5 * (max(ys) + min(ys)))
        # inside the frame: the axis titles are set in a font whose round
        # glyphs pass every other test here
        if not (fr[0] <= c[0] <= fr[1] and fr[2] <= c[1] <= fr[3]):
            continue
        (filled if "fill-rule" in attr else open_).append(c)
    # THE LAST MARKER OF EACH KIND IS THE LEGEND KEY, not a data point.  They
    # are told apart by count -- thirteen points and one key -- and the count
    # is checked, because a legend marker read as data would add ~100 events
    # at n_ch = 13 and still look like a plot.
    for lst, what in ((filled, "data"), (open_, "Monte Carlo")):
        if len(lst) != 14:
            die(f"FIG. 9: found {len(lst)} {what} markers, expected 14 "
                f"(13 points plus the legend key)")
    out = {}
    for lst, what in ((filled, "data"), (open_, "mc")):
        pts = sorted(lst[:13])
        nch = [int(round(nch_of_x(x) - 0.5)) for x, _y in pts]
        if nch != list(range(1, 14)):
            die(f"FIG. 9: the {what} markers sit at n_ch = {nch}, expected "
                f"1 through 13")
        out[what] = [n_of_y(y) for _x, y in pts]
    return {"nch": list(range(1, 14)), **out}


def main():
    pdf = sys.argv[1] if len(sys.argv) > 1 else None
    if pdf is None:
        pdf = os.path.join(tempfile.gettempdir(), "donut_0711.0728.pdf")
        if not os.path.exists(pdf):
            print(f"fetching {PDF_URL}")
            r = subprocess.run(["curl", "-sfL", "--max-time", "180", "-o",
                                pdf, PDF_URL])
            if r.returncode:
                die(f"could not fetch {PDF_URL}; pass the PDF as an argument")
    print(f"reading {pdf}")

    flux = read_flux(page_svg(pdf, PAGE_FLUX))
    nch = read_nch(page_svg(pdf, PAGE_NCH))

    # ---- the three checks ----
    tot = {k: sum(v) for k, v in flux["spectra"].items()}
    cc = sum(tot.values())
    frac = {k: v / cc for k, v in tot.items()}
    print("  FIG. 2 interacting-neutrino spectra, summed over energy:")
    for k in ("nu_e", "nu_mu", "nu_tau"):
        want = TABLE_III[k] / sum(TABLE_III.values())
        print(f"    {k:8s} {tot[k]:8.2f} events, {100*frac[k]:5.2f}% of the "
              f"CC total   (Table III: {100*want:5.2f}%)")
        if abs(frac[k] - want) > 0.08:
            die(f"FIG. 2's {k} fraction is {100*frac[k]:.1f}%, against Table "
                f"III's {100*want:.1f}%.  A location-efficiency difference "
                f"cannot be this large: the curves have been mixed up or the "
                f"log axis is wrong.")
    s_data, s_mc = sum(nch["data"]), sum(nch["mc"])
    print(f"  FIG. 9 sums: data {s_data:.1f}, Monte Carlo {s_mc:.1f}  "
          f"(the paper locates {N_LOCATED} events)")
    for s, what in ((s_data, "data"), (s_mc, "Monte Carlo")):
        if abs(s / N_LOCATED - 1.0) > 0.08:
            die(f"FIG. 9's {what} sums to {s:.1f} against {N_LOCATED} located "
                f"events -- the y axis or the marker list is wrong")

    os.makedirs(OUT, exist_ok=True)
    flux.update({
        "source": "arXiv:0711.0728 (DONUT, PRD 78 052002) FIG. 2",
        "what": ("energy spectra of neutrinos INTERACTING in the DONUT "
                 "emulsion target, per flavour, in bins of 6 GeV.  Already "
                 "carries the cross-section and the acceptance, so it is an "
                 "event weight and must not be folded with a cross-section "
                 "again."),
        "units": "events per 6 GeV bin",
        "table_III_cc_fractions": TABLE_III,
        "extracted_by": "tools/digitise_donut.py",
    })
    nch.update({
        "source": "arXiv:0711.0728 (DONUT, PRD 78 052002) FIG. 9",
        "what": ("charged-particle multiplicity at the primary vertex of the "
                 "578 located events; `mc` is DONUT's own LEPTO-based "
                 "simulation, normalised to the data"),
        "n_located": N_LOCATED,
        "data_err": [math.sqrt(max(v, 0.0)) for v in nch["data"]],
        "extracted_by": "tools/digitise_donut.py",
    })
    for name, obj in (("flux", flux), ("nch", nch)):
        p = os.path.join(OUT, f"{name}.json")
        with open(p, "w") as f:
            json.dump(obj, f, indent=1)
        print(f"  wrote {os.path.relpath(p, BASE)}")


if __name__ == "__main__":
    main()
