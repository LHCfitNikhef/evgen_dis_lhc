#!/usr/bin/env python3
"""sigma(charm) / sigma(inclusive) versus beam energy, for both currents.

SIX PANELS: the perturbative ORDER on the rows (LO, NLO), the charm
observable on the columns (total, >=1 charged D, >=1 neutral D).  The orders
are drawn separately rather than overlaid because they are not competing
predictions of the same quantity -- the ZM-VFNS charm series does not converge
at these kinematics (sigma_charm runs 3.879 -> 2.947 -> 1.915 nb from LO to
NNLO in the muon NC case), so an LO and an NLO generator sitting on different
lines is the expected physics, not a disagreement to be read off by eye.
Overlaying them invited exactly that misreading.  Each row carries its OWN
same-order YADISM reference, which is the convention every other comparison in
this benchmark uses.

WHY THIS RATIO IS WORTH PLOTTING SEPARATELY.  It is the one charm observable
that is free of the normalisation conventions this benchmark keeps tripping
over.  Numerator and denominator come from the SAME sample, so the overall
cross-section normalisation cancels exactly -- the Sherpa/Herwig neutrino
helicity factor 2, the Herwig attempted-vs-generated choice, the
integrator-vs-event-level choice, all of it.  What is left is the physics:
what fraction of deep-inelastic scattering produces charm.

THE THREE COLUMNS are NOT a partition: an event can contain both D species,
and charm can hadronise into a D_s or a charmed baryon instead, so the two D
columns do not add up to the first.  What they separate is the perturbative
channel (column 1, which tracks the matrix element) from the hadronisation
model (columns 2 and 3, which track the fragmentation).

DEFINITIONS -- and they do not all agree.
  * Generators: an event counts if it contains at least one PROMPT charm quark
    in the final state after shower and hadronisation, i.e. at least one charm
    hadron that is not a decay product of a bottom hadron (analyze.py).
  * YADISM: the charm-flavour structure-function contribution in the neutral
    current, and every term proportional to V_cd or V_cs in the charged
    current.  A fixed-order parton-level object, so it has no D mesons and
    appears only in the first column, drawn hollow.

THREE RULES THIS FILE FOLLOWS, each of them a bug that has already happened
here at least once:

  * NO ENERGY LITERAL IN A RESULT NAME.  Every per-energy name is built by
    `beams.at_energy()`, the same function analyze.py and the run scripts use.
    The previous version of this file hardcoded one dict of names per sample,
    and the neutrino table was still pinned at 1 TeV alone long after the scan
    had filled in 400 GeV and 4 TeV -- so the neutrino figure silently claimed
    "one energy so far" while the samples sat on disk.  The muon table had
    gone stale the same way for Sherpa MC@NLO, Herwig NLO and YADISM NNLO.
  * VERIFY THE ENERGY STAMP.  Every JSON carries `energy_gev`; a result whose
    stamp disagrees with the energy it is being plotted at is a hard error, not
    a warning.  That stamp is what caught the 1 TeV anchor being overwritten.
  * HONOUR `stat_insufficient`.  A result the analysis has already flagged as
    unmeasurable is skipped and listed, never drawn.  Storing a flag is not
    honouring it.

Usage: plot_charm_ratio.py [mu|nu]      (default: both)
Writes results{,_nu}/cmp_charm_ratio_E_{LO,NLO}.png and charm_ratio_vs_E.json
"""
import collections
import json
import math
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{BASE}/analysis")
# the house style, applied BEFORE pyplot is imported (it sets rcParams and
# selects the Agg backend)
import plotstyle  # noqa: E402
plotstyle.apply()
from plotstyle import tex  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402

import beams  # noqa: E402
import labels  # noqa: E402  -- the one display-name table
# bound directly: the plotting function has a local `labels` list for the
# legend, which shadows the module inside it
from labels import strip_order as _strip_order  # noqa: E402

# A sample is declared ONCE, with its BASE (anchor-energy) result names; the
# per-energy names come from beams.at_energy().  `energies=None` means "every
# configured energy", and a point whose files are absent is reported as
# missing rather than silently dropped -- so a newly generated energy appears
# in the figure with no edit here at all.
Sample = collections.namedtuple(
    "Sample", "key label row order color marker incl charm energies")


def S(key, label, row, order, color, marker, incl, charm, energies=None):
    return Sample(key, label, row, order, color, marker, incl, charm, energies)


# Rows: "LO" and "NLO".  GENIE is placed by the structure functions its tune
# actually uses -- HEDIS is BGR18 at NLO, the classic tunes are Bodek-Yang at
# LO -- and labelled "tune" so it is never mistaken for a matched calculation.
SAMPLES_NU = [
    S("sherpa_lo", "Sherpa 3.0.5 LO", "LO", "LO", "#3fa66a", "o",
      "sherpa_lo", "sherpa_lo_charmfinal"),
    S("pythia", "Pythia 8.311 LO", "LO", "LO", "#eb6834", "o",
      "pythia", "pythia_charmfinal"),
    S("herwig", "Herwig 7.3.0 LO", "LO", "LO", "#8e5bd0", "o",
      "herwig", "herwig_charmfinal"),
    S("genie_lo", labels.GENIE_NU, "LO", "tune", "#c2317b", "^",
      "genie_lo", "genie_lo_charmfinal"),
    S("yadism", "YADISM CC LO, ZM-VFNS (V_cd/V_cs)", "LO", "analytic", "#33383f", "D",
      "yadism", "yadism_charm"),

    S("sherpa_nlo", "Sherpa 3.0.5 MC@NLO", "NLO", "NLO", "#2a78d6", "s",
      "sherpa_nlo", "sherpa_nlo_charmfinal"),
    S("powheg_nu", "POWHEG-V2 NLO + Pythia8", "NLO", "NLO", "#0f9b9b", "s",
      "powheg_nu", "powheg_nu_charmfinal"),
    S("herwig_nlo", "Herwig 7.3.0 NLO (POWHEG)", "NLO", "NLO", "#6a3fb0", "s",
      "herwig_nlo_full", "herwig_nlo_full_charmfinal"),
    S("genie", labels.GENIE_NU_HEDIS, "NLO", "tune", "#e88bb8", "^",
      "genie", "genie_charmfinal"),
    S("yadism_nlo", "YADISM CC NLO, FONLL", "NLO", "analytic",
      "#7a8087", "D", "yadism_nlo_fonll_damp",
      "yadism_charm_nlo_fonll_damp"),
]

SAMPLES_MU = [
    S("sherpa_lo", "Sherpa 3.0.5 LO", "LO", "LO", "#3fa66a", "o",
      "sherpa_lo", "sherpa_lo_charmfinal"),
    S("pythia", "Pythia 8.311 LO", "LO", "LO", "#eb6834", "o",
      "pythia", "pythia_charmfinal"),
    S("herwig", "Herwig 7.3.0 LO", "LO", "LO", "#8e5bd0", "o",
      "herwig", "herwig_charmfinal"),
    S("yadism", "YADISM LO (ZM-VFNS)", "LO", "analytic", "#33383f", "D",
      "yadism", "yadism_charm"),

    S("sherpa", "Sherpa 3.0.5 MC@NLO", "NLO", "NLO", "#2a78d6", "s",
      "sherpa", "sherpa_charmfinal"),
    S("powheg", "POWHEG-RES NLO + Pythia8", "NLO", "NLO", "#0f9b9b", "s",
      "powheg", "powheg_charmfinal"),
    S("herwig_nlo_powheg", "Herwig 7.3.0 NLO (POWHEG)", "NLO", "NLO",
      "#6a3fb0", "s", "herwig_nlo_powheg_full",
      "herwig_nlo_powheg_full_charmfinal"),
    S("yadism_nlo", "YADISM NLO, FONLL", "NLO", "analytic",
      "#7a8087", "D", "yadism_nlo_fonll_damp",
      "yadism_charm_nlo_fonll_damp"),
    S("yadism_nnlo", "YADISM NNLO, FONLL", "NLO", "analytic",
      "#aab0b7", "D", "yadism_nnlo_fonll_damp",
      "yadism_charm_nnlo_fonll_damp"),
]

# THIS FOOTNOTE IS DELIBERATELY EMPTY, and the history is why it is not simply
# deleted.  Sherpa's CC MC@NLO once ran a UNIT CKM matrix -- the stock install
# has no virtual matrix element for the Cabibbo cross-family channels -- which
# left its charm composition 3.4 points low and the row stuck at 1 TeV.  Both
# were fixed by the two-line patch to DY_QCD_Virtual.C (see patches/), and all
# three energies now come from the full-CKM `NuDIS_NLO_ckm3*` runs: the charm
# fraction moved 18.69% -> 22.10% at 1 TeV, matching its own LO (22.81%).
#
# The caveat outlived the fix on this figure and had to be reported by the
# user, so: a caveat on a plot is a claim about the CURRENT sample, and when
# the sample changes it is as wrong as a stale number would be.
FOOTNOTE_NU = ""

# GENIE's muon EM DIS has no charm at all -- a measured zero, not a gap -- so
# it is deliberately absent from the muon figure rather than drawn as a flat
# line at the origin.  The charm TABLE reports it, with the reason.
FOOTNOTE_MU = ("GENIE is absent here by construction: its EM/NC DIS has no "
               "charm matrix element and AGKY makes no charm in "
               "fragmentation, so its charm fraction is a measured ZERO "
               "(0 of 200,000 events) rather than a point on this axis.")

ROWS = ["LO", "NLO"]

# THE DENOMINATOR OF THE RATIO PANEL, named explicitly per row.
# It cannot be inferred from order == "analytic": the muon NLO row carries TWO
# analytic curves, YADISM NLO and YADISM NNLO, and picking "the analytic one"
# would silently normalise to whichever happened to come first in the list.
# NNLO is context for the NLO row, not its reference -- so it appears IN the
# ratio panel, as the NNLO/NLO K-factor on the charm fraction, which is worth
# seeing in its own right.
REFERENCE = {"LO": "yadism", "NLO": "yadism_nlo"}

# THE CHARM REFERENCE IS FONLL (DAMPED) WHEREVER IT EXISTS (user, 2026-08-27).
# ZM-VFNS treats the charm as massless, which is wrong exactly where this
# benchmark sits: the ZM charm series does not converge at FASER kinematics
# (K = 0.760 then 0.650 on the muon side), and that is a threshold artefact of
# the scheme rather than a physical effect.  FONLL restores the mass and damps
# the mass correction near threshold, so it is the reference a charm result
# should be judged against.
#
# >>> THE FRACTION TAKES BOTH HALVES FROM ONE SCHEME. <<<
# The charm fraction is sigma_charm / sigma_inclusive, so moving the numerator
# to FONLL without the denominator plots FONLL charm over ZM inclusive -- two
# schemes in one ratio, which is what the analytic rows did until the user
# caught it.  The FONLL inclusive is light-ZM + charm-FONLL and is built by
# analysis/yadism_fonll_inclusive.py.
#
# LO IS THE ONE EXCEPTION, AND IT IS NOT AN OVERSIGHT: FONLL DOES NOT EXIST AT
# LO.  The massive contribution is gamma* g -> c cbar, which first appears at
# O(alpha_s), so there is nothing to combine with the massless description at
# O(alpha_s^0).  The LO rows therefore keep their ZM-VFNS reference, which is
# stated wherever they are shown rather than left looking like an omission.

PANELS = [
    ("total", "total charm", None),
    ("d_ch", "charm seen as $D^{\\pm}$", "frac_ge1_d_ch"),
    ("d_0", "charm seen as $D^{0}/\\bar{D}^{0}$", "frac_ge1_d_0"),
]


def _wrap(text, width):
    """Hard-wrap a footnote to `width` characters.

    fig.text(wrap=True) wraps to the figure box but does not report the
    resulting height, so the space beneath it cannot be reserved and the
    footnote ends up printed over the panels below.  Wrapping here means the
    line count is known before the layout is computed.
    """
    import textwrap
    return "\n".join(textwrap.wrap(text, width))


def load(resdir, name, e_gev):
    """One result JSON, with its energy stamp CHECKED against the energy asked.

    A name built for one energy that returns a result stamped with another is
    the anchor-overwrite failure mode, and it is silent in every other
    respect -- so it raises here rather than warning.
    """
    path = f"{resdir}/histos_{name}.json"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        d = json.load(f)
    stamp = d.get("energy_gev")
    if stamp is not None and abs(float(stamp) - e_gev) > 1e-6:
        raise SystemExit(
            f"ENERGY STAMP MISMATCH: {path} is stamped energy_gev={stamp} "
            f"but was resolved as the {e_gev:g} GeV result. One of the name "
            f"and the output path was made energy-aware without the other.")
    return d


def ratio_points(resdir, incl_name, charm_name, e_gev):
    """{panel key: (ratio, error)} for one generator at one energy.

    Returns (points, skip_reason); points is None when the sample is absent or
    has been flagged unmeasurable.

    The error is BINOMIAL on the tagged fraction, using the inclusive sample's
    fiducial event count: the charm events are a subset of the inclusive ones,
    drawn from the same run, so quadrature over two independent cross-section
    errors would be wrong (and much too large).  Analytic entries carry no
    event count and get no error bar.
    """
    incl = load(resdir, incl_name, e_gev)
    charm = load(resdir, charm_name, e_gev)
    if incl is None or charm is None:
        return None, "no sample"
    # A result the analysis has already judged unmeasurable must not be drawn.
    for tag, d in (("inclusive", incl), ("charm", charm)):
        if d.get("stat_insufficient"):
            why = d.get("stat_insufficient_why", "insufficient statistics")
            return None, f"{tag} flagged: {why}"
    if incl["sigma_fid_pb"] <= 0.0:
        return None, "inclusive sigma_fid <= 0"
    r_tot = charm["sigma_fid_pb"] / incl["sigma_fid_pb"]
    n = incl.get("n_fiducial")
    out = {}
    for key, _, frackey in PANELS:
        if frackey is None:
            r = r_tot
        elif charm.get("d_supported") and charm.get(frackey) is not None:
            # the species fraction is measured WITHIN the charm sample, so
            # multiplying by the charm fraction puts it back on sigma_incl
            r = r_tot * charm[frackey]
        else:
            continue                    # parton-level sample: no D mesons
        err = math.sqrt(max(r * (1.0 - r), 0.0) / n) if n else 0.0
        out[key] = (r, err)
    return out, None


def make(mode):
    """One figure per (current, perturbative order).

    LAYOUT: the TOTAL charm fraction gets a full-width panel of its own on
    top, and the two D-species panels sit side by side beneath it.  The three
    are not equals -- the total is the physics result and the species split is
    the breakdown of it -- and a flat 1x3 row gave them equal billing while
    making each one narrow.  Giving the total the full width also lets the
    generator spread, which is a few percent, actually be seen.

    LO and NLO go to SEPARATE FIGURES rather than separate rows of one figure,
    and each figure says which it is in its title.  They are not competing
    predictions of the same quantity: in the neutral current the ZM-VFNS charm
    series does not converge at these kinematics, so an LO and an NLO curve
    sitting apart is expected physics rather than a disagreement.  Each figure
    carries its OWN same-order YADISM reference.
    """
    nu = mode == "nu"
    resdir = f"{BASE}/results_nu" if nu else f"{BASE}/results"
    samples = SAMPLES_NU if nu else SAMPLES_MU
    beam = ("CC $\\nu_\\mu$ DIS" if nu else "NC $\\mu^-$ DIS")
    xlab = (r"$E_{\nu}$ [GeV]" if nu else r"$E_{\mu}$ [GeV]")

    records, missing = [], []
    written = []
    for row in ROWS:
        in_row = [s for s in samples if s.row == row]
        # EXPLICIT geometry, not tight_layout().  tight_layout ignores
        # fig.text, so the footnote landed on top of the lower panels; and
        # with a gridspec it silently declines to honour `rect` at all, so
        # reserving the space that way did not work either.  Fixed fractions
        # are predictable and this figure's content is a fixed shape.
        foot = (FOOTNOTE_NU if (nu and row == "NLO")
                else FOOTNOTE_MU if (not nu and row == "LO") else None)
        bottom = 0.235 if foot else 0.185
        fig = plt.figure(figsize=(11.0, 11.0))
        # THREE rows now: the total charm fraction, its ratio to the same-order
        # YADISM reference directly beneath it, then the two D-species panels.
        # The ratio panel is deliberately attached to the TOTAL and not to the
        # species panels -- the species fractions are a breakdown of the total
        # and YADISM has no D-meson prediction to normalise them against.
        # NESTED gridspecs, not one 3x2.  A single grid has ONE hspace, and
        # the two gaps here want opposite things: the ratio panel belongs to
        # the total panel and should sit flush beneath it, while the D-species
        # panels are a separate block and need room for their own titles.  With
        # a uniform hspace the ratio panel floated between the two and its
        # x-label landed on top of the D titles.
        outer = fig.add_gridspec(2, 1, height_ratios=[1.92, 1.0],
                                 left=0.085, right=0.975, top=0.845,
                                 bottom=bottom, hspace=0.34)
        top_gs = outer[0].subgridspec(2, 1, height_ratios=[1.30, 0.62],
                                      hspace=0.07)
        bot_gs = outer[1].subgridspec(1, 2, wspace=0.24)
        axes = {"total": fig.add_subplot(top_gs[0]),
                "d_ch": fig.add_subplot(bot_gs[0]),
                "d_0": fig.add_subplot(bot_gs[1])}
        rax = fig.add_subplot(top_gs[1], sharex=axes["total"])

        # The reference series, gathered BEFORE the sample loop so every
        # generator can be divided by it as it is drawn.
        ref_key = REFERENCE.get(row)
        ref_smp = next((x for x in in_row if x.key == ref_key), None)
        ref_tot = {}
        if ref_smp is not None:
            for e_ref in (ref_smp.energies or beams.ENERGIES):
                got_r, _why = ratio_points(
                    resdir, beams.at_energy(ref_smp.incl, e_ref),
                    beams.at_energy(ref_smp.charm, e_ref), e_ref)
                if got_r and got_r.get("total"):
                    ref_tot[e_ref] = got_r["total"][0]
        ratio_vals = []
        handles, labels = [], []
        drew = False
        for isamp, smp in enumerate(in_row):
            # a small horizontal fan so overlapping markers and their error
            # bars stay separable; the RECORDED energy is untouched.
            jitter = 1.014 ** (isamp - 0.5 * (len(in_row) - 1))
            series = {k: ([], [], []) for k, _, _ in PANELS}
            ratio_series = ([], [], [])
            for e_gev in (smp.energies or beams.ENERGIES):
                incl_name = beams.at_energy(smp.incl, e_gev)
                charm_name = beams.at_energy(smp.charm, e_gev)
                got, why = ratio_points(resdir, incl_name, charm_name, e_gev)
                if got is None:
                    missing.append(f"{smp.label} @ {e_gev:g} GeV ({why})")
                    continue
                rec = {"generator": smp.key, "label": smp.label,
                       "row": smp.row, "order": smp.order, "E_GeV": e_gev,
                       "inclusive": incl_name, "charm": charm_name}
                for pkey, (r, err) in got.items():
                    series[pkey][0].append(e_gev * jitter)
                    series[pkey][1].append(r)
                    series[pkey][2].append(err)
                    rec[pkey] = r
                    rec[pkey + "_err"] = err
                    # ratio to the same-order reference.  The reference itself
                    # is skipped -- it would be a row of exact 1.0 -- and the
                    # unity line stands in for it.
                    if (pkey == "total" and smp.key != ref_key
                            and ref_tot.get(e_gev)):
                        rr = r / ref_tot[e_gev]
                        ratio_series[0].append(e_gev * jitter)
                        ratio_series[1].append(rr)
                        ratio_series[2].append(err / ref_tot[e_gev])
                        ratio_vals.append(rr)
                        rec["total_over_reference"] = rr
                records.append(rec)
            hollow = smp.order == "analytic"
            # THE TITLE ALREADY SAYS THE ORDER (user, 2026-08-28), so the
            # legend does not repeat it: the "[LO]"/"[NLO]" tag went on every
            # generator row of a figure headed "-- LO generators", and the
            # order also sat inside most of the names themselves.  A row at a
            # DIFFERENT order from the figure's keeps it -- strip_order()
            # removes only the figure's own -- because there it is what tells
            # the row apart rather than what repeats the heading.
            label = _strip_order(smp.label, row)
            if smp.order not in ("analytic", "tune") and smp.order != row:
                label = f"{label} [{smp.order}]"
            for pkey, _, _ in PANELS:
                es, rs, errs = series[pkey]
                if not es:
                    continue
                drew = True
                h = axes[pkey].errorbar(
                    es, rs, yerr=errs, color=smp.color, marker=smp.marker,
                    markersize=7.5,
                    linestyle="-" if len(es) > 1 else "none",
                    markerfacecolor="none" if hollow else smp.color,
                    markeredgewidth=1.6, capsize=3, elinewidth=1.2,
                    label=label)
                if label not in labels:
                    handles.append(h)
                    labels.append(label)
            if ratio_series[0]:
                rax.errorbar(ratio_series[0], ratio_series[1],
                             yerr=ratio_series[2], color=smp.color,
                             marker=smp.marker, markersize=7.0,
                             linestyle="-" if len(ratio_series[0]) > 1
                             else "none",
                             markerfacecolor="none" if hollow else smp.color,
                             markeredgewidth=1.6, capsize=3, elinewidth=1.2)
        if not drew:
            plt.close(fig)
            continue

        # THE TWO D PANELS SHARE A Y RANGE.  They are the same quantity split
        # by species, so a reader compares them against each other; separate
        # autoscales make D+- and D0 look equally large when one is routinely
        # twice the other.  The total panel keeps its own scale, being a
        # different quantity.
        d_lims = [ax.get_ylim() for k, ax in axes.items() if k != "total"]
        if d_lims:
            lo = min(l[0] for l in d_lims)
            hi = max(l[1] for l in d_lims)
            for k, ax in axes.items():
                if k != "total":
                    ax.set_ylim(lo, hi)

        for pkey, title, _ in PANELS:
            ax = axes[pkey]
            ax.set_xscale("log")
            ax.set_xlim(250.0, 6.5e3)
            # Explicit ticks AT the scanned energies.  Matplotlib's default
            # log minor labels ("3x10^2, 4x10^2, ...") overprint each other
            # into an unreadable smear, and none of them is an energy the
            # benchmark was actually run at.
            ax.set_xticks(list(beams.ENERGIES))
            ax.set_xticklabels([f"{e:g}" for e in beams.ENERGIES])
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            plotstyle.ticks(ax)
            # tex() keeps the $...$ maths verbatim and wraps only the prose
            ax.set_xlabel(tex(xlab), fontsize=plotstyle.FS_XLABEL)
            ax.set_title(tex(title), fontsize=plotstyle.FS_PANEL_TITLE)
            ax.set_ylabel(r"$\sigma_{\rm charm}/\sigma_{\rm inclusive}$",
                          fontsize=plotstyle.FS_YLABEL)
        axes["total"].set_title(tex("total charm fraction"),
                                fontsize=plotstyle.FS_PANEL_TITLE)
        # The total panel hands its x-axis to the ratio panel beneath it.
        axes["total"].set_xlabel("")
        plt.setp(axes["total"].get_xticklabels(), visible=False)

        # ---- the ratio panel -------------------------------------------
        # SHORT, because this panel is only ~0.6 of the main one and a long
        # label overflows upward into the main panel's own y-label: the full
        # "YADISM CC NLO, FONLL" rendered across both, so the left
        # edge read "FONLL sigma_charm/sigma_inclusive".
        ref_label = ("FONLL"
                     if ref_smp is not None and "FONLL" in ref_smp.label
                     else (ref_smp.label.split(" (")[0]
                           if ref_smp is not None else f"YADISM {row}, ZM-VFNS"))
        rax.axhline(1.0, color=(ref_smp.color if ref_smp is not None
                                else "#33383f"),
                    linestyle="--", linewidth=1.3, zorder=1)
        rax.set_xscale("log")
        rax.set_xlim(250.0, 6.5e3)
        rax.set_xticks(list(beams.ENERGIES))
        rax.set_xticklabels([f"{e:g}" for e in beams.ENERGIES])
        rax.xaxis.set_minor_formatter(ticker.NullFormatter())
        plotstyle.ticks(rax)
        # no x-label here: the D panels beneath repeat this axis and label it,
        # and a label in this position collides with their titles.
        rax.set_ylabel(tex("ratio to") + "\n" + tex(ref_label),
                       fontsize=plotstyle.FS_YLABEL)
        # AUTO-RANGED, like the ratio panels in plot_compare.py.  A fixed
        # window either clips a generator that is genuinely far off -- GENIE
        # sits ~0.6 of the CC charm reference -- or squashes the others into a
        # line when they all agree.  The floor on the padding keeps a panel
        # where everything lands on 1.000 from being blank.
        if ratio_vals:
            lo_r, hi_r = min(ratio_vals), max(ratio_vals)
            pad = max(0.06 * (hi_r - lo_r), 0.02)
            rax.set_ylim(lo_r - pad, hi_r + pad)
        else:
            rax.text(0.5, 0.5, tex(f"no {ref_label} reference at this order"),
                     transform=rax.transAxes, ha="center", va="center",
                     fontsize=9, color="#8a9098", style="italic")
            rax.set_yticks([])

        # The title names the ORDER and lists the generators in it, so a
        # reader never has to infer from the legend which block they are
        # looking at -- the two figures are otherwise near-identical in form.
        gen_names = ", ".join(x.label.split(" (")[0]
                              for x in in_row if x.order != "analytic")
        fig.suptitle(tex(f"{beam}: charm fraction vs beam energy"
                         f"  --  {row} generators"),
                     fontsize=plotstyle.FS_SUPTITLE, y=0.975)
        # The generator list goes on its own WRAPPED line rather than inside
        # suptitle: as a second title line it ran past both edges of the
        # canvas and was clipped, which is worse than not naming them at all.
        fig.text(0.5, 0.935,
                 _wrap(gen_names + "   (+ same-order YADISM reference, "
                       "drawn hollow)", 112),
                 ha="center", va="top", fontsize=10.0, color="#2c3038")
        ncol = 3 if len(labels) <= 6 else 4
        fig.legend(handles, [tex(x) for x in labels], loc="lower center",
                   ncol=ncol, fontsize=plotstyle.FS_LEGEND, frameon=True,
                   bbox_to_anchor=(0.5, 0.012))
        note = ("$Q^2>4$ GeV$^2$, $0.2<y<0.9$;  charm = $\\geq 1$ prompt "
                "charm hadron in the final state;  points fanned out "
                "horizontally for legibility")
        # bottom-up: legend, then the footnote if any, then the note, then
        # the axes (whose x-labels hang ~0.05 below the gridspec bottom).
        fig.text(0.5, bottom - 0.068, note, ha="center", va="top",
                 fontsize=8.4, style="italic", color="#5d6470")
        if foot:
            fig.text(0.5, bottom - 0.098, _wrap(foot, 120), ha="center",
                     va="top", fontsize=8.0, style="italic", color="#5d6470")

        os.makedirs(resdir, exist_ok=True)
        png = f"{resdir}/cmp_charm_ratio_E_{row}.png"
        fig.savefig(png, dpi=150)
        plt.close(fig)
        written.append(png)

    with open(f"{resdir}/charm_ratio_vs_E.json", "w") as f:
        json.dump({"note": "charm fraction of the fiducial cross-section; "
                           "total, and the part with >=1 charged / neutral D. "
                           "Charm = >=1 PROMPT charm hadron in the final "
                           "state (analytic rows: fixed-order charm terms)",
                   "current": "CC nu_mu" if nu else "NC mu-",
                   "energies_gev": list(beams.ENERGIES),
                   "points": records}, f, indent=1)
    for rec in sorted(records, key=lambda r: (r["row"] != "LO", r["label"],
                                              r["E_GeV"])):
        print(f"{rec['row']:3s} {rec['label']:32s} "
              f"E = {rec['E_GeV']:7.1f} GeV   total {rec['total']:.4f}"
              + (f"   D+- {rec['d_ch']:.4f}   D0 {rec['d_0']:.4f}"
                 if "d_ch" in rec else "   (parton level)"))
    for m in missing:
        print(f"MISSING (not plotted): {m}")
    for png in written:
        print(f"wrote {png}")


def main():
    modes = sys.argv[1:] or ["mu", "nu"]
    for m in modes:
        if m not in ("mu", "nu"):
            sys.exit(f"unknown mode {m} (use mu|nu)")
        print(f"--- {m} ---")
        make(m)


if __name__ == "__main__":
    main()
