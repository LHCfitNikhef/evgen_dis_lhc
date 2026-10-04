#!/usr/bin/env python3
"""Write the mu- + mu+ SUM of the FASERnu muon flux as one LHAPDF set.

WHY.  POWHEG-RES's lepton-flux mode reads the flux for the beam flavour
`ih1` alone, and its LHAPDF interface sets f(-13) = f(13) whatever the grid
holds, so a single run cannot see both charges.  FASERnu muon DIS is photon
exchange at Q2 << m_Z^2, which does not know the lepton charge, so the
mu+ + mu- rate of arXiv:2506.13889 (its Table 2.1 sums the two) is one run
of a mu- beam on the SUMMED flux.  The set is written beside the vendored
originals under data/faser_muon_flux/lhapdf/, with index 111222334 (the
original mu-/mu+ grid is 111222333).

Usage: tools/make_muon_flux_sum.py   (idempotent; writes two files)
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data", "faser_muon_flux")
DST = os.path.join(SRC, "lhapdf")
NAME = "muon_flux_FASERv_Run3_var2"
OUT = NAME + "_sum"
INDEX = 111222334


def main():
    with open(os.path.join(SRC, NAME + "_0000.dat")) as f:
        parts = f.read().split("---")
    head, block = parts[0], parts[1].strip().split("\n")
    xs, qs, fl = block[0], block[1], [int(v) for v in block[2].split()]
    im, ip = fl.index(13), fl.index(-13)
    rows = []
    for ln in block[3:]:
        if not ln.strip():
            continue
        v = [float(t) for t in ln.split()]
        s = v[im] + v[ip]
        rows.append(f"{s:.8e} {s:.8e}")      # f(-13) = f(13) = the sum
    odir = os.path.join(DST, OUT)
    os.makedirs(odir, exist_ok=True)
    with open(os.path.join(odir, OUT + "_0000.dat"), "w") as f:
        f.write(head + "---\n" + xs + "\n" + qs + "\n-13 13\n"
                + "\n".join(rows) + "\n---\n")
    with open(os.path.join(SRC, NAME + ".info")) as f:
        info = f.read()
    info = re.sub(r"SetDesc: .*", 'SetDesc: "Muon flux 25x30 cm, mu- plus '
                  'mu+ summed (tools/make_muon_flux_sum.py)"', info)
    info = re.sub(r"SetIndex: .*", f"SetIndex: {INDEX}", info)
    with open(os.path.join(odir, OUT + ".info"), "w") as f:
        f.write(info)
    idx = os.path.join(DST, "pdfsets.index")
    lines = open(idx).read().splitlines() if os.path.exists(idx) else []
    lines = [l for l in lines if not l.startswith(str(INDEX))]
    lines.append(f"{INDEX} {OUT} 1")
    with open(idx, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {odir} (index {INDEX}), {len(rows)} x nodes")


if __name__ == "__main__":
    main()
