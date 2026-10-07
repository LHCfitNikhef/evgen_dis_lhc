#!/usr/bin/env python3
"""Write rwl_npdf_hadron.{xml,json}: the nuclear-PDF reweighting list for the
hadron-level nPDF bands (analysis/npdf_hadron.py, paper plot A2b).

The sets are those of Fig. C.1 (analysis/npdf_impact.py), every member of the
three nuclear sets and the central member of each free-nucleon baseline:

  NNPDF40_nnlo_as_01180            closure: the generation PDF, must return the nominal
  NNPDF40_nnlo_as_01180_W184free   the benchmark's free tungsten nucleon
  nNNPDF30_nlo_as_0118_A184_Z74    200 MC replicas      + own baseline nNNPDF30_..._p_W184free
  EPPS21nlo_CT18Anlo_W184          106 Hessian, 90% CL  + own baseline CT18ANLO_W184free
  nCTEQ15HQ_FullNuc_184_74         38 Hessian, 90% CL   + own baseline nCTEQ15HQ_1_1_W184free

Ids are LHAPDF numbers, read from each set's .info SetIndex.  POWHEG-BOX
addresses PDFs by number only, and the number is checked here to resolve
BACK to the intended set and member -- memory lhapdf-index-collision: a
SetIndex need not.  nCTEQ15HQ_FullNuc_184_74 is absent from the conda
pdfsets.index (its 122950 would fall inside nCTEQ15WZSIH's range), so it is
listed in data/pdfs/lhapdf/pdfsets.index, which the reweighting driver
prepends to LHAPDF_DATA_PATH.

Usage: LHAPDF_DATA_PATH=data/pdfs/lhapdf:$LHAPDF_DATA_PATH make_rwl_npdf_hadron.py
"""
import json
import os

import lhapdf

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = [  # (name, members used)
    ("NNPDF40_nnlo_as_01180", 1),
    ("NNPDF40_nnlo_as_01180_W184free", 1),
    ("nNNPDF30_nlo_as_0118_A184_Z74", None),
    ("nNNPDF30_nlo_as_0118_p_W184free", 1),
    ("EPPS21nlo_CT18Anlo_W184", None),
    ("CT18ANLO_W184free", 1),
    ("nCTEQ15HQ_FullNuc_184_74", None),
    ("nCTEQ15HQ_1_1_W184free", 1),
]
FIRST_ID = 5000


def main():
    lhapdf.setVerbosity(0)
    sets, lines, wid = {}, [], FIRST_ID
    for name, nuse in SETS:
        ps = lhapdf.getPDFSet(name)
        base = int(ps.lhapdfID)
        nmem = ps.size
        nuse = nmem if nuse is None else nuse
        ids = []
        for m in range(nuse):
            p = lhapdf.mkPDF(base + m)
            if p.set().name != name or p.memberID != m:
                raise SystemExit(f"LHAPDF id {base + m} resolves to "
                                 f"{p.set().name}/{p.memberID}, not {name}/{m}")
            lines.append(f"<weight id='{wid}'> lhapdf={base + m} </weight>")
            ids.append(str(wid))
            wid += 1
        sets[name] = {"weight_ids": ids, "lhapdf_base": base, "members": nuse,
                      "members_in_set": nmem, "error_type": ps.errorType,
                      "conf_level": ps.errorConfLevel}
    xml = ("<initrwgt>\n<weightgroup name='npdf_hadron' combine='none'>\n"
           + "\n".join(lines) + "\n</weightgroup>\n</initrwgt>\n")
    with open(os.path.join(HERE, "rwl_npdf_hadron.xml"), "w") as f:
        f.write(xml)
    with open(os.path.join(HERE, "rwl_npdf_hadron.json"), "w") as f:
        json.dump({"sets": sets, "closure_set": SETS[0][0],
                   "closure_id": sets[SETS[0][0]]["weight_ids"][0],
                   "note": "weight id -> (set, member); the closure id "
                           "re-evaluates each event with the PDF it was "
                           "generated with and must return the nominal"},
                  f, indent=1)
    print(f"{wid - FIRST_ID} weights in rwl_npdf_hadron.xml")


if __name__ == "__main__":
    main()
