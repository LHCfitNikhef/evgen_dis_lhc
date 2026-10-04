#!/usr/bin/env python3
"""Reweighting XML for the FASER dimuon ladder: six PDF sets with every
member, on a PROTON and on a NEUTRON, plus the seven scale points.

WHY A SECOND FILE beside rwl_members.xml.  The dimuon ladder of
tools/powheg_v2_dimuon_ladder.sh has neutron points, generated on the
isospin-mirrored NNPDF4.0 (id 339900), and a neutron point can only be
reweighted to a set that exists AS A NEUTRON: tools/make_neutron_pdf.py
--all-members writes mirrored copies of all six sets under
data/pdfs/lhapdf/ with their own ids (pdfsets.index there).  So two XMLs
are written, `rwl_dimuon_p.xml` and `rwl_dimuon_n.xml`, with THE SAME
WEIGHT IDS mapping to the proton set and to its mirror respectively, so the
analysis can read a p point and an n point by one id.

GRV98 IS IN (user, 2026-09-07: "GRV98, the latter without PDF error band").
It was unresolvable by id -- LHAPDF's global index hands 80060 to
METAv10LHC -- and data/pdfs/lhapdf/pdfsets.index now carries the line
"80060 GRV98lo 1", which the prepended data path makes win.  Every id is
still resolved back to its set name before it is written, as in
make_rwl_members.py, for the same reason: a wrong resolution would produce
plausible weights.

The id scheme keeps rwl_members.py's: 1001-1007 scales, 3000+N members in
the order of SETS.  The central member of each set is its member 0.

Usage: make_rwl_dimuon.py            (writes rwl_dimuon_{p,n}.{xml,json})
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = ["NNPDF40_nnlo_as_01180", "CT18NNLO", "MSHT20nnlo_as118",
        "ATLASpdf21_T1", "ABMP16_5_nnlo", "GRV98lo"]
SCALES = [("1001", "1d0", "1d0"), ("1002", "2d0", "1d0"),
          ("1003", "0.5d0", "1d0"), ("1004", "1d0", "2d0"),
          ("1005", "1d0", "0.5d0"), ("1006", "2d0", "2d0"),
          ("1007", "0.5d0", "0.5d0")]


def build(target):
    import lhapdf
    lhapdf.setVerbosity(0)
    suffix = "" if target == "p" else "_n"
    lines = ["<initrwgt>",
             "<weightgroup name='scale_variation' combine='envelope'>"]
    for wid, r, f in SCALES:
        lines.append(f"<weight id='{wid}'> renscfact={r} facscfact={f} </weight>")
    lines.append("</weightgroup>")
    lines.append("<weightgroup name='PDF_members' combine='none'>")
    mapping, n, dropped = {}, 0, []
    for name in SETS:
        pname = name + suffix
        pset = lhapdf.getPDFSet(pname)
        base, size = pset.lhapdfID, pset.size
        ids = []
        for imem in range(size):
            pid = base + imem
            try:
                got = lhapdf.mkPDF(pid)
                gname, gmem = got.set().name, got.memberID
            except Exception as exc:                          # noqa: BLE001
                dropped.append((pname, imem, pid, f"unresolvable: {exc}"))
                n += 1          # keep the ids aligned between p and n
                continue
            if gname != pname or gmem != imem:
                dropped.append((pname, imem, pid, f"resolves to {gname}/{gmem}"))
                n += 1
                continue
            wid = f"{3000 + n}"
            lines.append(f"<weight id='{wid}'> lhapdf={pid} </weight>")
            ids.append(wid)
            n += 1
        mapping[name] = {"weight_ids": ids, "lhapdf_base": base,
                         "lhapdf_set": pname, "members": len(ids),
                         "members_in_set": size,
                         "error_type": lhapdf.getPDFSet(name).errorType}
        print(f"  {target}: {pname:28s} ids {ids[0]}-{ids[-1]}  ({len(ids)} of {size})")
    lines += ["</weightgroup>", "</initrwgt>"]
    if dropped:
        print(f"  {target}: {len(dropped)} member(s) DROPPED:")
        for a, b, c, d in dropped:
            print(f"    {a} member {b} (lhapdf={c}): {d}")
    with open(f"{HERE}/rwl_dimuon_{target}.xml", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(f"{HERE}/rwl_dimuon_{target}.json", "w") as f:
        json.dump({"target": target, "sets": mapping,
                   "scales": [s[0] for s in SCALES],
                   "dropped": [{"set": a, "member": b, "lhapdf_id": c, "why": d}
                               for a, b, c, d in dropped]}, f, indent=1)
    return n


def main():
    n_p = build("p")
    n_n = build("n")
    if n_p != n_n:
        sys.exit(f"proton and neutron member counts differ ({n_p} vs {n_n})")
    print(f"  {n_p} member weights + 7 scales per file")


if __name__ == "__main__":
    main()
