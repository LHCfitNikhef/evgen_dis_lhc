#!/usr/bin/env python3
"""Sum the channel splines of a gmkspl XML file at a given energy and emit
the events_xsec.json record analyze.py expects for unweighted samples.

gmkspl stores xsec knots in natural units (GeV^-2); 1 GeV^-2 = 0.3894 mb.
Stdlib only: this also runs under the GENIE conda python (no numpy there).

Usage: spline_to_json.py <spline.xml> <E_GeV> <n_accepted>
"""
import json
import math
import sys
import xml.etree.ElementTree as ET

GEV2_TO_MB = 0.3893793722

def interp_loglog(e0, knots):
    """Linear interpolation in (log E, log xsec) between bracketing knots."""
    lx = math.log(e0)
    for (e1, x1), (e2, x2) in zip(knots, knots[1:]):
        if e1 <= e0 <= e2:
            l1, l2 = math.log(e1), math.log(e2)
            t = 0.0 if l2 == l1 else (lx - l1) / (l2 - l1)
            return math.exp(math.log(x1) + t*(math.log(x2) - math.log(x1)))
    raise ValueError(f"E = {e0} outside spline knot range "
                     f"[{knots[0][0]}, {knots[-1][0]}]")

def main():
    fname, e0, nacc = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    root = ET.parse(fname).getroot()
    total = 0.0
    for spl in root.iter("spline"):
        knots = []
        for knot in spl.iter("knot"):
            e = float(knot.find("E").text)
            x = float(knot.find("xsec").text)
            if e > 0 and x > 0:
                knots.append((e, x))
        # channels that are zero everywhere (or dead at E0, e.g. below a
        # heavy-quark threshold) contribute nothing
        if not knots or not knots[0][0] <= e0 <= knots[-1][0]:
            continue
        total += interp_loglog(e0, knots)
    print(json.dumps({"sigma_gen_mb": total * GEV2_TO_MB,
                      "sigma_err_mb": 0.0,
                      "n_accepted": nacc}))

if __name__ == "__main__":
    main()
