#!/usr/bin/env python3
"""Write an isospin-swapped copy of the benchmark PDF as a 'neutron' LHAPDF set.

WHY.  Tungsten is 74 protons and 110 neutrons, and the neutron cross-section
is the larger one for neutrino charged currents (nu n -> mu- p scatters off
the d quark).  Every generator here takes its target through LHAPDF, so the
generator-independent way to scatter on a neutron is a PDF set in which
u <-> d and ubar <-> dbar are exchanged -- exact isospin symmetry, the same
construction 2402.13318's flux authors and every DIS fit use.  A set rather
than a generator switch, so POWHEG-V2, POWHEG-RES and Pythia all get the
same neutron.

Member 0 only by default: the rates need the central value.  Written under
data/pdfs/lhapdf/<set>_n/ with a pdfsets.index entry (id 339900), to be
PREPENDED to LHAPDF_DATA_PATH like the muon-flux grid.

--all-members mirrors EVERY member (2026-09-07, for the PDF bands on the
FASER dimuon rate: POWHEG-V2's neutron ladder points can only be reweighted
to a set that exists as a neutron), keeping the source's ErrorType so
LHAPDF's own uncertainty() applies to the mirror as to the original.  The
id must then be given, and the ids of the mirrored sets must not overlap
member ranges: id + member is how POWHEG addresses a member.

Usage: tools/make_neutron_pdf.py [SET] [--all-members] [--index N]
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
_args = [a for a in sys.argv[1:] if not a.startswith("--")]
SET = _args[0] if _args else "NNPDF40_nnlo_as_01180"
ALL_MEMBERS = "--all-members" in sys.argv
if "--index" in sys.argv:
    _INDEX_ARG = int(sys.argv[sys.argv.index("--index") + 1])
else:
    _INDEX_ARG = None
INDEX = 339900 if _INDEX_ARG is None else _INDEX_ARG
SWAP = {1: 2, 2: 1, -1: -2, -2: -1}


def lhapdf_path(setname):
    out = subprocess.run(["lhapdf-config", "--datadir"], capture_output=True,
                         text=True).stdout.strip()
    cands = [os.path.join(p, setname) for p in
             os.environ.get("LHAPDF_DATA_PATH", "").split(":") + [out] if p]
    for c in cands:
        if os.path.isdir(c):
            return c
    raise SystemExit(f"cannot find the set {setname} in LHAPDF_DATA_PATH")


def main():
    src = lhapdf_path(SET)
    out_name = SET + "_n"
    odir = os.path.join(BASE, "data", "pdfs", "lhapdf", out_name)
    os.makedirs(odir, exist_ok=True)
    # --- the member grids, columns permuted (member 0, or every member)
    import glob as _glob
    members = sorted(_glob.glob(os.path.join(src, SET + "_[0-9][0-9][0-9][0-9].dat")))
    if not ALL_MEMBERS:
        members = members[:1]
    for mfile in members:
        with open(mfile) as f:
            txt = f.read()
        parts = txt.split("---")
        head, blocks = parts[0], parts[1:]
        new_blocks = []
        for b in blocks:
            b = b.strip("\n")
            if not b.strip():
                new_blocks.append(b)
                continue
            lines = b.split("\n")
            xs, qs, fl = lines[0], lines[1], [int(v) for v in lines[2].split()]
            perm = [fl.index(SWAP.get(pid, pid)) for pid in fl]
            rows = []
            for ln in lines[3:]:
                if not ln.strip():
                    continue
                v = ln.split()
                rows.append(" ".join(v[i] for i in perm))
            new_blocks.append("\n".join([xs, qs, lines[2]] + rows))
        suffix = os.path.basename(mfile)[len(SET):]          # _NNNN.dat
        with open(os.path.join(odir, out_name + suffix), "w") as f:
            f.write(head + "---\n" + "\n---\n".join(new_blocks) + "\n---\n")
    n_members = len(members)
    # --- the info file: one member, its own index, a description that says so
    with open(os.path.join(src, SET + ".info")) as f:
        info = f.read()
    what = "every member" if ALL_MEMBERS else "member 0"
    # the description may be a MULTI-LINE quoted string (ATLASpdf21's is);
    # replacing its first line alone leaves an unterminated YAML value
    info = re.sub(r'SetDesc: "(?:[^"\\]|\\.)*"', f'SetDesc: "{SET} {what} with u<->d and '
                  f'ubar<->dbar exchanged: a NEUTRON target by isospin '
                  f'(tools/make_neutron_pdf.py)"', info, count=1, flags=re.S)
    info = re.sub(r"SetIndex: .*", f"SetIndex: {INDEX}", info)
    info = re.sub(r"NumMembers: .*", f"NumMembers: {n_members}", info)
    if not ALL_MEMBERS:
        info = re.sub(r"ErrorType: .*", "ErrorType: replicas", info)
    with open(os.path.join(odir, out_name + ".info"), "w") as f:
        f.write(info)
    idx = os.path.join(BASE, "data", "pdfs", "lhapdf", "pdfsets.index")
    lines = open(idx).read().splitlines() if os.path.exists(idx) else []
    lines = [l for l in lines if not l.startswith(str(INDEX) + " ")]
    lines.append(f"{INDEX} {out_name} {n_members}")
    with open(idx, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {odir} ({n_members} member(s), LHAPDF id {INDEX})")


if __name__ == "__main__":
    main()
