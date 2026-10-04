#!/usr/bin/env python3
"""Write the BOUND-PROTON PDF of a nucleus from its average-nucleon LHAPDF set.

WHY.  The FASER-format productions (faser_format/) generate on the struck
nucleon -- a proton or a neutron of the tungsten nucleus, drawn by their
columns -- and POWHEG-V2 builds the neutron from the set it is given by
swapping u <-> d itself (ih2 2).  So what a generator needs from a nuclear PDF
is the BOUND PROTON f^{p/A}.  But the nNNPDF3.0 (and EPPS21) LHAPDF grids are
the AVERAGE NUCLEON of the nucleus, the isospin mix already done -- measured,
see tools/make_isoscalar_pdf.py and memory nuclear-pdf-impact:

    A f_u^A = Z u_p + N d_p        (n = p with u <-> d, ubar <-> dbar)
    A f_d^A = Z d_p + N u_p

nNNPDF3.0 is built exactly this way from its bound-proton parametrisation
(isospin symmetry of the bound nucleons), so the mix inverts exactly:

    u_p = A (Z f_u^A - N f_d^A) / (Z^2 - N^2)      and the same for ubar,
    d_p = A (Z f_d^A - N f_u^A) / (Z^2 - N^2)      dbar; s, c, b, g unchanged.

Fed to a generator on a proton (ih2 1) and a neutron (ih2 2), in the ratio Z:N,
this reproduces the average nucleon exactly (DIS is linear in the PDF) AND
gives each event a remnant of the right charge -- which one run on the
average-nucleon set with a proton beam would not.

The inversion amplifies node noise by |A Z / (Z^2 - N^2)| ~ 2 for tungsten,
harmless on the central member; the check below re-mixes the written set
through LHAPDF and requires the source back, and the valence sum rules of the
bound proton (u_v = 2, d_v = 1) are printed.

Usage:
  tools/make_bound_proton_pdf.py SET --index N [--Z 74 --N 110] [--name OUT]
  tools/make_bound_proton_pdf.py nNNPDF30_nlo_as_0118_A184_Z74 --index 342000
"""
import argparse
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from make_isoscalar_pdf import lhapdf_path  # noqa: E402

SWAP = {1: 2, 2: 1, -1: -2, -2: -1}


def invert_block(block, Z, N):
    lines = block.strip("\n").split("\n")
    if not "".join(lines).strip():
        return block.strip("\n")
    xs, qs, flav_line = lines[0], lines[1], lines[2]
    fl = [int(v) for v in flav_line.split()]
    mirror = [fl.index(SWAP.get(pid, pid)) for pid in fl]
    A, D = Z + N, Z * Z - N * N
    rows = []
    for ln in lines[3:]:
        if not ln.strip():
            continue
        v = [float(t) for t in ln.split()]
        out = []
        for j, pid in enumerate(fl):
            if pid in SWAP:
                out.append(A * (Z * v[j] - N * v[mirror[j]]) / D)
            else:
                out.append(v[j])
        rows.append(" ".join(f"{o: .9E}" for o in out))
    return "\n".join([xs, qs, flav_line] + rows)


def check(out_name, src_name, Z, N):
    import lhapdf
    lhapdf.setVerbosity(0)
    new, old = lhapdf.mkPDF(out_name, 0), lhapdf.mkPDF(src_name, 0)
    A = Z + N
    worst, where = 0.0, None
    for x in (1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.6):
        for q in (2.0, 5.0, 10.0, 100.0):
            for pid in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21):
                got = (Z * new.xfxQ(pid, x, q)
                       + N * new.xfxQ(SWAP.get(pid, pid), x, q)) / A
                want = old.xfxQ(pid, x, q)
                d = abs(got - want) / max(abs(want), 1e-3)
                if d > worst:
                    worst, where = d, (pid, x, q, got, want)
    if worst > 2e-4:
        raise SystemExit(f"make_bound_proton_pdf: re-mixing the written set "
                         f"does not give the source back: {where} "
                         f"(relative {worst:.2e})")
    print(f"  re-mixed through LHAPDF: worst relative deviation {worst:.2e} "
          f"over 264 (flavour, x, Q) points")
    # valence sum rules of the bound proton, at Q = 10 GeV
    import numpy as np
    lx = np.linspace(np.log(1e-6), 0.0, 4001)[:-1]
    xs = np.exp(lx)
    for q, name, want in ((2, "u_v", 2.0), (1, "d_v", 1.0)):
        f = np.array([new.xfxQ(q, x, 10.0) - new.xfxQ(-q, x, 10.0) for x in xs])
        s = np.trapezoid(f, lx) if hasattr(np, "trapezoid") else np.trapz(f, lx)
        print(f"  bound proton {name} = {s:.3f} (a proton: {want:.0f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set")
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--Z", type=int, default=74)
    ap.add_argument("--N", type=int, default=110)
    ap.add_argument("--name")
    a = ap.parse_args()
    out_name = a.name or f"{a.set}_boundp"
    src = lhapdf_path(a.set)
    odir = os.path.join(BASE, "data", "pdfs", "lhapdf", out_name)
    os.makedirs(odir, exist_ok=True)
    mfile = sorted(glob.glob(os.path.join(src, a.set + "_[0-9][0-9][0-9][0-9].dat")))[0]
    with open(mfile) as f:
        parts = f.read().split("---")
    blocks = [invert_block(b, a.Z, a.N) for b in parts[1:]]
    with open(os.path.join(odir, out_name + "_0000.dat"), "w") as f:
        f.write(parts[0] + "---\n" + "\n---\n".join(blocks) + "\n---\n")
    info = open(os.path.join(src, a.set + ".info")).read()
    desc = (f"{a.set} member 0 inverted to the BOUND PROTON of the Z = {a.Z}, "
            f"A = {a.Z + a.N} nucleus (isospin symmetry of the bound nucleons; "
            f"tools/make_bound_proton_pdf.py)")
    info = re.sub(r'SetDesc: "(?:[^"\\]|\\.)*"', f'SetDesc: "{desc}"', info,
                  count=1, flags=re.S)
    info = re.sub(r"SetIndex: .*", f"SetIndex: {a.index}", info)
    info = re.sub(r"NumMembers: .*", "NumMembers: 1", info)
    info = re.sub(r"ErrorType: .*", "ErrorType: replicas", info)
    with open(os.path.join(odir, out_name + ".info"), "w") as f:
        f.write(info)
    idx = os.path.join(BASE, "data", "pdfs", "lhapdf", "pdfsets.index")
    lines = [ln for ln in open(idx).read().splitlines()
             if not ln.startswith(f"{a.index} ") and ln.split()[1:2] != [out_name]]
    lines.append(f"{a.index} {out_name} 1")
    with open(idx, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {odir} (LHAPDF id {a.index})")
    check(out_name, a.set, a.Z, a.N)


if __name__ == "__main__":
    main()
