#!/usr/bin/env python3
"""Write the FREE-NUCLEON tungsten average of a proton PDF as an LHAPDF set.

WHY.  The nuclear-PDF study compares the benchmark's own prediction for a
tungsten target -- NNPDF4.0 with nuclear effects ASSUMED TO VANISH -- against
nNNPDF3.0 and EPPS21, which are fitted for tungsten.  The comparison is only
like-for-like if the baseline is the same TARGET: tungsten is 74 protons and
110 neutrons, so the object to compare a nuclear PDF against is the average
free nucleon

    f_i^{W,free} = (Z f_i^p + N f_i^n) / A ,   n = p with u <-> d, ubar <-> dbar

and NOT the free proton.  Comparing against the proton would fold the isospin
of the target into what the figure calls a nuclear effect, and for neutrino
charged currents the isospin term is the LARGER of the two by a wide margin
(nu d -> mu- u runs on the valence d, and tungsten has 1.5 neutrons per
proton).  Both lines are shown in the study, and this tool is what makes the
first one exist.

>>> AND BOTH nPDF SETS ARE ALREADY THE AVERAGE NUCLEON, NOT THE BOUND PROTON.
The LHAPDF grids of nNNPDF30_nlo_as_0118_A184_Z74 and of
EPPS21nlo_CT18Anlo_W184 both carry the PER-NUCLEON PDF of the nucleus, with
the isospin mix already done -- measured here, not assumed: at x = 0.1,
Q = 10 GeV both give u/d = 0.88-0.91, against 1.59 for their own free proton
and 0.92 for the free-nucleon average this tool builds.  A bound-PROTON grid
would have needed the mix applying to it as well; these do not, and applying
it twice would be a silent 10% error on the neutrino rate.

WHY A PDF SET AND NOT A COMBINATION OF TWO CROSS-SECTIONS.  Deep-inelastic
scattering is linear in the PDF at every order -- sigma = sum_i C_i (x) f_i --
so the cross-section of the average nucleon IS the average of the two
cross-sections, exactly.  Doing it in the PDF costs one POWHEG reweighting
pass per member instead of two, and the generator then needs no knowledge of
the target at all.  The same argument is why tools/make_neutron_pdf.py exists.

The mixing is exact on the grid nodes and LHAPDF's interpolation is linear in
the node values, so it is exact between them too.  alpha_s and every other
piece of metadata come from the source set unchanged.

Usage:
  tools/make_isoscalar_pdf.py SET --index N [--all-members] [--name OUT]
                              [--Z 74 --N 110]
  tools/make_isoscalar_pdf.py NNPDF40_nnlo_as_01180 --index 341000 --all-members
  tools/make_isoscalar_pdf.py NNPDF40_nnlo_as_01180 --index 341400 --Z 26 --N 30
"""
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "analysis"))

import target as target_mod                                   # noqa: E402

# Tungsten by default -- the FASERnu target and what this study is about --
# but the mix is a parameter, because DONUT's target is steel and the same
# construction is what a free-nucleon iron prediction needs.
Z, N = target_mod.Z_W, target_mod.N_W
A = Z + N
SWAP = {1: 2, 2: 1, -1: -2, -2: -1}


def opt(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    return default


def lhapdf_path(setname):
    out = subprocess.run(["lhapdf-config", "--datadir"], capture_output=True,
                         text=True).stdout.strip()
    cands = [os.path.join(p, setname) for p in
             os.environ.get("LHAPDF_DATA_PATH", "").split(":") + [out] if p]
    for c in cands:
        if os.path.isdir(c):
            return c
    raise SystemExit(f"cannot find the set {setname} in LHAPDF_DATA_PATH")


def mix_block(block):
    """One '---'-delimited grid block, columns mixed into the average nucleon."""
    lines = block.strip("\n").split("\n")
    if not "".join(lines).strip():
        return block.strip("\n")
    xs, qs, flav_line = lines[0], lines[1], lines[2]
    fl = [int(v) for v in flav_line.split()]
    # column j of the OUTPUT is (Z * col_j + N * col_of_mirror(j)) / A
    mirror = [fl.index(SWAP.get(pid, pid)) for pid in fl]
    rows = []
    for ln in lines[3:]:
        if not ln.strip():
            continue
        v = [float(t) for t in ln.split()]
        rows.append(" ".join(
            f"{(Z * v[j] + N * v[mirror[j]]) / A: .7E}" for j in range(len(v))))
    return "\n".join([xs, qs, flav_line] + rows)


def check(out_name, src_name, n_members):
    """Reload both sets through LHAPDF and verify the mix, node by node.

    CONVENTIONS.md rule 2: a grid file that parses and gives plausible numbers is
    exactly the failure mode this benchmark keeps meeting.  A permuted column,
    a dropped row or a format that LHAPDF re-interpolates differently would
    all pass "it loaded".
    """
    import lhapdf
    lhapdf.setVerbosity(0)
    new, old = lhapdf.mkPDF(out_name, 0), lhapdf.mkPDF(src_name, 0)
    worst, where = 0.0, None
    for x in (1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.6):
        for q in (2.0, 5.0, 10.0, 100.0):
            for pid in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21):
                got = new.xfxQ(pid, x, q)
                want = (Z * old.xfxQ(pid, x, q)
                        + N * old.xfxQ(SWAP.get(pid, pid), x, q)) / A
                scale = max(abs(want), 1e-3)
                d = abs(got - want) / scale
                if d > worst:
                    worst, where = d, (pid, x, q, got, want)
    if worst > 2e-4:
        pid, x, q, got, want = where
        raise SystemExit(
            f"make_isoscalar_pdf: the written set does not reproduce the "
            f"mix: flavour {pid} at x={x:g}, Q={q:g} gives {got:.6g}, "
            f"expected {want:.6g} (relative {worst:.2e})")
    print(f"  checked against LHAPDF: worst relative deviation {worst:.2e} "
          f"over 264 (flavour, x, Q) points")
    if n_members > 1:
        # the member count must be what the .info claims, or POWHEG will
        # address a member that is not there and LHAPDF will hand back the
        # nearest set instead -- see powheg/reweight/make_rwl_members.py
        got_n = lhapdf.getPDFSet(out_name).size
        if got_n != n_members:
            raise SystemExit(f"make_isoscalar_pdf: {out_name} reports "
                             f"{got_n} members, wrote {n_members}")


def main():
    global Z, N, A
    if "--Z" in sys.argv or "--N" in sys.argv:
        Z = int(opt("--Z", Z))
        N = int(opt("--N", N))
        A = Z + N
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        raise SystemExit(__doc__)
    src_set = args[0]
    all_members = "--all-members" in sys.argv
    index = opt("--index")
    if index is None:
        raise SystemExit("--index is required: LHAPDF addresses a member as "
                         "id + member, so the id ranges must not overlap")
    index = int(index)
    sym = {184: "W", 56: "Fe"}.get(A, f"A{A}")
    out_name = opt("--name", f"{src_set}_{sym}{A}free")

    src = lhapdf_path(src_set)
    odir = os.path.join(BASE, "data", "pdfs", "lhapdf", out_name)
    os.makedirs(odir, exist_ok=True)

    members = sorted(glob.glob(os.path.join(src, src_set + "_[0-9][0-9][0-9][0-9].dat")))
    if not all_members:
        members = members[:1]
    for mfile in members:
        with open(mfile) as f:
            parts = f.read().split("---")
        head, blocks = parts[0], parts[1:]
        new_blocks = [mix_block(b) for b in blocks]
        suffix = os.path.basename(mfile)[len(src_set):]        # _NNNN.dat
        with open(os.path.join(odir, out_name + suffix), "w") as f:
            f.write(head + "---\n" + "\n---\n".join(new_blocks) + "\n---\n")
    n_members = len(members)

    with open(os.path.join(src, src_set + ".info")) as f:
        info = f.read()
    what = "every member" if all_members else "member 0"
    desc = (f"{src_set} {what} combined as ({Z} p + {N} n) / {A}: the average "
            f"FREE nucleon of the Z = {Z}, A = {A} nucleus, i.e. the target "
            f"with nuclear effects assumed to vanish "
            f"(tools/make_isoscalar_pdf.py)")
    info = re.sub(r'SetDesc: "(?:[^"\\]|\\.)*"', f'SetDesc: "{desc}"', info,
                  count=1, flags=re.S)
    info = re.sub(r"SetDesc: '(?:[^'\\]|\\.)*'", f'SetDesc: "{desc}"', info,
                  count=1, flags=re.S)
    info = re.sub(r"SetIndex: .*", f"SetIndex: {index}", info)
    info = re.sub(r"NumMembers: .*", f"NumMembers: {n_members}", info)
    if not all_members:
        # one member is not an error band; say replicas so LHAPDF's own
        # uncertainty() cannot be fooled into reading a Hessian pair
        info = re.sub(r"ErrorType: .*", "ErrorType: replicas", info)
    with open(os.path.join(odir, out_name + ".info"), "w") as f:
        f.write(info)

    idx = os.path.join(BASE, "data", "pdfs", "lhapdf", "pdfsets.index")
    lines = open(idx).read().splitlines() if os.path.exists(idx) else []
    lines = [ln for ln in lines
             if not ln.startswith(str(index) + " ") and
             not ln.split()[1:2] == [out_name]]
    lines.append(f"{index} {out_name} {n_members}")
    with open(idx, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"wrote {odir} ({n_members} member(s), LHAPDF id {index})")
    check(out_name, src_set, n_members)


if __name__ == "__main__":
    main()
