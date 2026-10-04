#!/usr/bin/env python3
"""Display names for the samples, in ONE place.

WHY THIS EXISTS.  A sample's name was written down in five different places:
the `label` field each analysis stamps into its result JSON (which the plot
legends printed raw), a rename table inside make_report.py that only the
HTML tables consulted, and hand-written row labels in make_report.SCAN_ROWS,
plot_sigma_vs_E.py and plot_charm_ratio.py.  They drifted, as duplicated
naming rules in this repository always have: the same GENIE sample appeared
as "GENIE LO (GRV98LO, charm)" in a figure legend and
"GENIE G18_02a (default tune = FASER)" in the table directly beneath it, and
the muon scan table called it simply "GENIE".

GENIE is where this actually hurts, because it enters the benchmark FOUR
times under two modules and two tunes, and two of those four are easy to
confuse with each other in exactly the way that matters:

  * G18_02a EMDIS      muon NC, Bodek-Yang on GRV98LO, LO.  FASER's default.
  * G18_02a EMDIS with NNPDF4.0 swapped into BYPDF -- a diagnostic variant,
    NOT a FASER configuration.
  * G18_02a            neutrino CC, the same classic tune.  FASER's default.
  * HEDIS GHE19_00a    neutrino CC, BGR18 structure functions at NLO through
    APFEL.  NOT the FASER tune, and at a different perturbative order, so it
    belongs against the NLO reference rather than the LO one.

The canonical names below say the tune, the module where it disambiguates,
and whether the row is the configuration FASER actually runs -- the last
being the one a reader is most likely to want and least able to infer.

    from labels import display
    display(result_json["label"])        # -> the canonical name

`display()` passes anything it does not recognise straight through, so
adding a generator needs no edit here; only names that are genuinely
ambiguous are worth centralising.
"""

import re as _re

# The four GENIE entries.  Kept as named constants because the row tables in
# plot_sigma_vs_E.py and plot_charm_ratio.py refer to them directly -- they
# build their rows from generator KEYS and never see a result's own label.
# SHORTENED (user, 2026-08-28).  These ran to 39 characters and crowded every
# legend they appeared in.  The distinguishing fact between the two G18_02a
# rows is the PARTON DISTRIBUTION -- Bodek-Yang on GRV98 LO as shipped, or
# NNPDF4.0 substituted into BYPDF as a diagnostic -- so naming the PDF says
# what the row is in a third of the width.
#
# WHAT THE SHORT FORM DROPS, and where it went.  "(FASER default tune)" marked
# the row FASER actually runs, which CONVENTIONS.md rule 1b calls the thing a FASER
# reader most wants and is least able to infer.  It is now carried by the
# report prose and the figure footnotes instead of by every legend entry --
# see FASER_DEFAULT_NOTE below, which those call sites use.  The module and
# tune are unchanged; only the parenthesis is shorter.
#
# The muon and neutrino names are now the SAME string, which is safe because
# no table or legend mixes the two currents -- and it is honest, since it is
# in fact the same tune run on two different beams.
#
# >>> STANDING RULE (user, 2026-10-04): EXACTLY THESE THREE NAMES, IN EVERY
# FIGURE. <<<  "GENIE (GRV98LO)", "GENIE (NNPDF4.0)", "GENIE (HEDIS)" -- the
# minimum that tells the three predictions apart.  No tune names, no
# "(FASER tune)", no "(default)": the tunes are stated once in the paper
# (Sect. 2) and in FASER_DEFAULT_NOTE.  Every hand-written legend tuple in
# analysis/ uses these constants, so the rule has one place to change.
GENIE_MU = "GENIE (GRV98LO)"
GENIE_MU_NNPDF = "GENIE (NNPDF4.0)"
GENIE_NU = "GENIE (GRV98LO)"
GENIE_NU_HEDIS = "GENIE (HEDIS)"
# the CC counterpart of GENIE_MU_NNPDF: same tune, NNPDF4.0
# substituted into Bodek-Yang in place of GRV98 LO
GENIE_NU_NNPDF = "GENIE (NNPDF4.0)"

# The fact the long labels used to carry, for the figures and prose that need
# it.  Kept here so it has ONE wording, like the names themselves.
FASER_DEFAULT_NOTE = ("GENIE G18_02a on GRV98 is the tune FASER runs; "
                      "the NNPDF4.0 row is a diagnostic variant, not a "
                      "FASER configuration.")



def _charm(name):
    """The charm-tagged variant of a canonical name.

    The marker goes INSIDE an existing parenthesis -- "(FASER default tune,
    charm)" -- rather than being appended as a second group, so these read
    like their neighbours in the same table: "Pythia 8.311 LO (ME, charm)".

    A name with NO trailing group gets its own.  That used to be an assertion,
    on the assumption that every canonical name ends in one; dropping the
    "(not the FASER tune)" suffix from the HEDIS label falsified it and the
    module stopped importing. Failing loudly was the right behaviour -- a
    silent "GENIE HEDIS GHE19_00, charm" would have been worse -- but the
    assumption itself was never necessary.
    """
    return (name[:-1] + ", charm)") if name.endswith(")") else (name + " (charm)")

# Keyed on the sample's OWN label, NOT on the generator key.  The key "genie"
# means the muon EMDIS sample under results/ but the HEDIS one under
# results_nu/, so a key-based table silently stamped "default tune = FASER"
# onto the HEDIS rows -- announcing precisely the confusion this file exists
# to prevent.  A label is unique across both directories.
# The QED study's shower arms.  Their STORED labels carry the full
# provenance -- which code, which arm, and on the CC side the diagonal-CKM
# caveat -- because that is what a result JSON is for.  A figure legend
# carrying three of them at once was wider than the axes and covered the data,
# so the legend gets these short forms instead.  The caveat is not lost: it is
# in the view subtitle, in the report text and in the stored label itself.
_QED_SHORT = {
    "POWHEG-RES NLO + Pythia8, QED FSR (lepton line)":
        "POWHEG-RES, QED FSR",
    "POWHEG-RES NLO + Pythia8, QED FSR+ISR (lepton line)":
        "POWHEG-RES, QED FSR+ISR",
    "POWHEG-RES NLO + Pythia8, QED FSR+ISR + quark line":
        "POWHEG-RES, QED FSR+ISR+quark",
    "POWHEG-RES NLO + Pythia8, QED FSR (lepton line) "
    "(diagonal CKM, inclusive only)": "POWHEG-RES, QED FSR",
    "POWHEG-RES NLO + Pythia8, QED FSR+ISR + quark line "
    "(diagonal CKM, inclusive only)": "POWHEG-RES, QED FSR+ISR+quark",
    "POWHEG-RES NLO + Pythia8 (diagonal CKM, inclusive only)":
        "POWHEG-RES, QED off",
}

# YADISM ALWAYS NAMES ITS SCHEME (user rule, 2026-08-27).  The stored labels
# predate the rule and say only the order, which does not identify a
# calculation: ZM-VFNS and FONLL differ by 7-13% on charm.  Mapped here rather
# than restamped into 542 result files, because the JSON label is provenance --
# it records what the run called itself -- while THIS table is what a reader
# sees.  A label that already names FONLL is passed through unchanged.
_YADISM_SCHEME = {
    "YADISM LO":              "YADISM LO (ZM-VFNS)",
    "YADISM NLO":             "YADISM NLO (ZM-VFNS)",
    "YADISM NNLO":            "YADISM NNLO (ZM-VFNS)",
    "YADISM CC LO":           "YADISM CC LO (ZM-VFNS)",
    "YADISM CC NLO":          "YADISM CC NLO (ZM-VFNS)",
    "YADISM CC NNLO":         "YADISM CC NNLO (ZM-VFNS)",
    "YADISM LO (charm)":      "YADISM LO (ZM-VFNS, charm)",
    "YADISM NLO (charm)":     "YADISM NLO (ZM-VFNS, charm)",
    "YADISM NNLO (charm)":    "YADISM NNLO (ZM-VFNS, charm)",
    "YADISM CC LO (charm)":   "YADISM CC LO (ZM-VFNS, charm)",
    "YADISM CC NLO (charm)":  "YADISM CC NLO (ZM-VFNS, charm)",
    "YADISM LO (bottom)":     "YADISM LO (ZM-VFNS, bottom)",
    "YADISM NLO (charm, FONLL)":     "YADISM NLO (FONLL, charm)",
    "YADISM NNLO (charm, FONLL)":    "YADISM NNLO (FONLL, charm)",
    "YADISM CC NLO (charm, FONLL)":  "YADISM CC NLO (FONLL, charm)",
    "YADISM, FONLL":                 "YADISM (FONLL)",
}

BY_LABEL = {
    # POWHEG-V2mc: the massive-charm CC entry.  Named distinctly everywhere so
    # it can never be read as the massless POWHEG-V2 -- they are the same code
    # at different settings, which is exactly the confusion a shared label
    # would create.
    "POWHEG-V2mc NLO + Pythia8 (massive charm, charm)":
        "POWHEG-V2mc (massive charm)",
    "POWHEG-V2mc NLO + Pythia8 (massive charm)":
        "POWHEG-V2mc (massive charm)",
    "GENIE G18_02a (EMDIS)": GENIE_MU,
    "GENIE G18_02a (EMDIS) (charm)": _charm(GENIE_MU),
    "GENIE G18_02a (EMDIS, NNPDF4.0)": GENIE_MU_NNPDF,
    "GENIE LO (GRV98LO)": GENIE_NU,
    "GENIE LO (GRV98LO, charm)": _charm(GENIE_NU),
    "GENIE LO (NNPDF4.0)": GENIE_NU_NNPDF,
    "GENIE LO (NNPDF4.0, charm)": _charm(GENIE_NU_NNPDF),
    "GENIE HEDIS (BGR18 NLO)": GENIE_NU_HEDIS,
    "GENIE HEDIS (BGR18 NLO, charm)": _charm(GENIE_NU_HEDIS),
    "GENIE HEDIS BGR18 NLO (nu CC)": GENIE_NU_HEDIS,
    "GENIE HEDIS BGR18 NLO (nu NC)": GENIE_NU_HEDIS,
    # stored labels that spell out the tune (analyze_nu's nugen samples), and
    # the display names used before 2026-10-04
    "GENIE G18_02a (GRV98)": GENIE_NU,
    "GENIE G18_02a (NNPDF4.0)": GENIE_NU_NNPDF,
    "GENIE HEDIS GHE19_00a": GENIE_NU_HEDIS,
    **_QED_SHORT,
    **_YADISM_SCHEME,
}


def display(label):
    """The canonical display name for a sample's stored label."""
    return BY_LABEL.get(label, label)


# A figure whose TITLE already says the order should not repeat it in every
# legend entry (user, 2026-08-28).
_ORDER_RE = {
    # `(?<![A-Za-z@])` keeps MC@NLO intact -- dropping the order there would
    # leave "Sherpa 3.0.5 MC@", since for that generator NLO is part of the
    # matching's NAME rather than a statement of the order.
    # `\b` on the left already prevents NNLO from matching NLO, because the
    # inner "NLO" of "NNLO" is preceded by a word character.
    "LO": _re.compile(r"(?<![A-Za-z@])\bLO\b"),
    "NLO": _re.compile(r"(?<![A-Za-z@])\bNLO\b"),
    "NNLO": _re.compile(r"(?<![A-Za-z@])\bNNLO\b"),
}


def strip_order(name, order):
    """Drop the perturbative order from a legend entry that repeats the title.

    ONLY the figure's OWN order is removed.  A row at a DIFFERENT order than
    the figure's -- the NNLO reference drawn on an NLO figure, or GENIE HEDIS
    which is NLO on a figure of LO generators -- keeps its label, because
    there the order is what distinguishes the row rather than what repeats the
    title.  Removing it everywhere would have quietly made those rows look
    like the rest.
    """
    rx = _ORDER_RE.get(order)
    if not rx:
        return name
    out = rx.sub("", name)
    # tidy what removal leaves behind: "YADISM CC NLO, FONLL" -> "YADISM CC,
    # FONLL", "Pythia 8.311 LO (Lund)" -> "Pythia 8.311 (Lund)", and a name
    # that was ONLY the order is left alone rather than blanked.
    out = _re.sub(r"\s+([,)])", r"\1", out)
    out = _re.sub(r"\(\s*,\s*", "(", out)
    out = _re.sub(r"\s{2,}", " ", out).strip(" ,")
    return out or name


if __name__ == "__main__":
    for k, v in sorted(BY_LABEL.items()):
        print(f"{k!r:40s} -> {v!r}")
