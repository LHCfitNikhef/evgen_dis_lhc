#!/usr/bin/env python3
"""Generate a POWHEG rwl file carrying EVERY member of each PDF set.

WHY THIS SUPERSEDES THE CENTRAL-ONLY FILE.  `rwl_scale_pdf.xml` carries the
seven scale points and the CENTRAL member of four PDF sets, and the README
next to it records why the error members were left out: "each extra member
costs a full matrix-element re-evaluation per event, so the ~100 NNPDF
replicas would cost about ten times this whole run".

That estimate was pessimistic, and it was never measured.  Timed on 2000
events with ten variations, the cost is 87 microseconds per (event, weight),
so all 309 members of the six sets over 50000 events comes to roughly half an
hour per beam energy -- an hour and a half for the whole scan.  Affordable, so
the PDF uncertainty can now be computed from the Monte Carlo as well as from
YADISM, and the two compared.

WHAT THE IDs MEAN.  POWHEG selects a PDF by its LHAPDF id, and members of a
set are consecutive from the set's own base id.  Those bases are read from
LHAPDF here rather than written down, because a hardcoded id that drifts
against the installed set would silently reweight to the wrong PDF -- the
weights would still vary, and still look reasonable.

THE WEIGHT ID SCHEME is deliberate and must stay stable, because the analysis
maps ids back to (set, member):

    1001-1007   the seven scale points, unchanged from the central-only file
    2001-2004   the four central members, unchanged, so existing results and
                the published nlo_unc_* JSONs keep their meaning
    3000 + N    member N of the concatenated member list, in the order the
                sets appear below; the mapping is written alongside as JSON

Usage: make_rwl_members.py [out.xml] [out.json]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# The sets, in the order their members are numbered.  Same six as
# analysis/pdf_dependence.py, so the two studies compare like with like.
SETS = ["NNPDF40_nnlo_as_01180", "CT18NNLO", "MSHT20nnlo_as118",
        "ATLASpdf21_T1", "ABMP16_5_nnlo", "GRV98lo"]

SCALES = [("1001", "1d0", "1d0"), ("1002", "2d0", "1d0"),
          ("1003", "0.5d0", "1d0"), ("1004", "1d0", "2d0"),
          ("1005", "1d0", "0.5d0"), ("1006", "2d0", "2d0"),
          ("1007", "0.5d0", "0.5d0")]


def main():
    out_xml = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/rwl_members.xml"
    out_map = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/rwl_members.json"
    import lhapdf
    lhapdf.setVerbosity(0)

    lines = ["<initrwgt>",
             "<weightgroup name='scale_variation' combine='envelope'>"]
    for wid, r, f in SCALES:
        lines.append(f"<weight id='{wid}'> renscfact={r} facscfact={f} "
                     f"</weight>")
    lines.append("</weightgroup>")

    # the four central members keep their original ids so nothing already
    # produced changes meaning
    central = {"NNPDF40_nnlo_as_01180": "2001", "CT18NNLO": "2002",
               "MSHT20nnlo_as118": "2003", "ATLASpdf21_T1": "2004"}
    lines.append("<weightgroup name='PDF_sets' combine='none'>")
    for name, wid in central.items():
        base = lhapdf.getPDFSet(name).lhapdfID
        lines.append(f"<weight id='{wid}'> lhapdf={base} </weight>")
    lines.append("</weightgroup>")

    mapping, n, dropped = {}, 0, []
    lines.append("<weightgroup name='PDF_members' combine='none'>")
    for name in SETS:
        pset = lhapdf.getPDFSet(name)
        base, size = pset.lhapdfID, pset.size
        ids = []
        for imem in range(size):
            # EVERY ID IS RESOLVED BACK BEFORE IT IS WRITTEN.  POWHEG selects a
            # PDF by NUMBER through LHAPDF's global index, and a set's own
            # `SetIndex` is not guaranteed to be registered in that index: the
            # GRV98lo shipped here declares 80060, which pdfsets.index does not
            # list, so the lookup walks back to the nearest lower entry
            # (80000 = METAv10LHC) and pwhg_main aborts on a set that is not
            # installed.  That killed pass 54 of 54 after five hours of
            # reweighting, taking the merge with it.
            #
            # A resolution that lands on the WRONG set would be worse than the
            # abort, since the weights would vary plausibly, so the check is
            # that the id comes back as the set we asked for -- not merely that
            # it comes back.
            pid = base + imem
            try:
                got = lhapdf.mkPDF(pid).set().name
            except Exception as exc:                          # noqa: BLE001
                dropped.append((name, imem, pid, f"unresolvable: {exc}"))
                continue
            if got != name:
                dropped.append((name, imem, pid,
                                f"resolves to {got}, not {name}"))
                continue
            wid = f"{3000 + n}"
            lines.append(f"<weight id='{wid}'> lhapdf={pid} "
                         f"</weight>")
            ids.append(wid)
            n += 1
        if not ids:
            print(f"  {name:24s} SKIPPED ENTIRELY -- no usable member id")
            continue
        mapping[name] = {"weight_ids": ids, "lhapdf_base": base,
                         "members": len(ids), "members_in_set": size,
                         "error_type": pset.errorType}
        print(f"  {name:24s} ids {ids[0]}-{ids[-1]}  ({len(ids)} of {size} "
              f"members)")
    lines.append("</weightgroup>")
    # NEVER SILENTLY.  A dropped member is a set the Monte Carlo arm cannot
    # cover, which the figures have to say rather than imply by absence.
    if dropped:
        print(f"\n  {len(dropped)} member(s) DROPPED as unusable:")
        for name, imem, pid, why in dropped:
            print(f"    {name} member {imem} (lhapdf={pid}): {why}")
    lines.append("</initrwgt>")

    with open(out_xml, "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(out_map, "w") as f:
        json.dump({"sets": mapping, "scales": [s[0] for s in SCALES],
                   "central": central,
                   "dropped": [{"set": a, "member": b, "lhapdf_id": c,
                                "why": d} for a, b, c, d in dropped],
                   "note": "weight id -> (set, member); see the module "
                           "docstring for the id scheme"}, f, indent=1)
    print(f"  wrote {out_xml} with {7 + len(central) + n} weights")
    print(f"  wrote {out_map}")


if __name__ == "__main__":
    main()
