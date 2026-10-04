#!/usr/bin/env python3
"""Fiducial selections for the benchmark.

ONE PLACE defining what "in the fiducial region" means, so that the muon
analysis, the neutrino analysis and the ME-level samples cannot drift apart --
they used to carry `Q2_MIN, Y_MIN, Y_MAX` as module constants in three files.

Two things need this, and they are the same problem:

  * the ENERGY SCAN, where the cuts stay put and the beam moves;
  * the FASER SELECTION, a second fiducial region beside the inclusive one.

A selection is deliberately a plain object rather than a set of globals: the
analysis is run once per (selection, energy) pair and the result carries the
names of both, so a plot can never silently mix two regions.

Pick one from the environment:

    BENCH_SELECTION=inclusive analysis/analyze.py pythia
"""
import math
import os


# ------------------------------------------------------------- the Q2 floor
#
# BENCH_Q2MIN RAISES THE FIDUCIAL Q2 FLOOR ON EVERY SELECTION AT ONCE.  It was
# a knob on the histogram binning alone (analyze.bins_for) and on the analytic
# integrated band (mhou_sigma_fonll), which meant "the Q2 > 11 study" existed
# for the reference calculation and NOT for the generator samples -- the
# asymmetry CONVENTIONS.md rule 2 is about, since both sides printed perfectly
# plausible numbers.  Setting it here makes one floor, and every consumer that
# reads a Selection gets it: the event loops, the YADISM references, the
# differential scale bands.
#
# EVERY RESULT COMPUTED WITH IT SET CARRIES THE VALUE IN ITS FILENAME
# (q2_suffix below), so a Q2 > 11 number can never be mistaken for a Q2 > 4
# one -- the same rule the energy tag and the selection suffix follow, and for
# the same reason: both traps have been sprung here before.
#
# THE DEFAULT IS 4.0 AND THAT IS STILL THE BENCHMARK'S REGION on the neutrino
# charged current.  The muon neutral current is now shown at 11 in the paper
# figures (user, 2026-08-31), because the neutral-current scale variation and
# the NNLO correction are both dominated by the bottom of the Q2 range; the
# two currents therefore carry DIFFERENT floors on purpose, and each figure
# states its own.  Going back is setting one variable.
Q2_FLOOR = float(os.environ.get("BENCH_Q2MIN", 4.0))
Q2_FLOOR_DEFAULT = 4.0


def q2_suffix(q2_min=None):
    """Filename suffix for a non-default Q2 floor: "" or "_q2min11"."""
    v = Q2_FLOOR if q2_min is None else float(q2_min)
    return "" if v == Q2_FLOOR_DEFAULT else f"_q2min{v:g}"


class Selection:
    """A fiducial region.

    Always a cut on the leptonic invariants; the FASER tiers add cuts on the
    scattered-lepton angle and, for Tier E, on the charged hadronic final
    state.  The two are applied at different points in the event loop -- the
    leptonic ones before the hadrons are walked, the hadronic ones after -- so
    they are separate predicates rather than one.
    """

    def __init__(self, name, q2_min, y_min, y_max, label=None,
                 e_lep_min=None, theta_max=None, theta_min=None,
                 n05_min=None, n01_min=None, dphi_x_min=None,
                 dimuon=False, w2_min=None, x_min=None, note=""):
        self.name = name
        self.q2_min = float(q2_min)
        self.y_min = float(y_min)
        self.y_max = float(y_max)
        self.label = label or name
        # FASER tiers.  None means "not applied", which is what keeps the
        # inclusive selection bit-for-bit what it always was.
        self.e_lep_min = e_lep_min        # GeV, scattered lepton, lab
        self.theta_max = theta_max        # rad
        self.theta_min = theta_min        # rad
        self.n05_min = n05_min            # charged hadrons, tan(theta) < 0.5
        self.n01_min = n01_min            # charged hadrons, tan(theta) < 0.1
        self.dphi_x_min = dphi_x_min      # rad, lepton vs summed hadron system
        # dimuon signal: require a SECOND muon of OPPOSITE SIGN passing the
        # same energy and angle cuts as the primary.  False keeps every other
        # selection bit-for-bit what it was.
        self.dimuon = bool(dimuon)
        # HADRONIC INVARIANT MASS, W^2 = m_p^2 + 2 P.q - Q^2, in GeV^2.
        # None means "not applied", which is what keeps every selection that
        # existed before this field bit-for-bit what it was -- the benchmark's
        # own region is defined on (Q2, y) and stays so.  It exists for the
        # ONE comparison that needs a foreign region: arXiv:2506.13889 cuts on
        # (Q, W) rather than (Q2, y), and a ratio against their number is only
        # a measurement of anything if both sides use the same region.
        self.w2_min = None if w2_min is None else float(w2_min)
        # BJORKEN x FLOOR, same discipline as w2_min: None means "not
        # applied", so every selection that predates this field is
        # bit-for-bit what it was.  It exists for ONE comparison, the SIDIS
        # regions of arXiv:2504.05376, which cut on x > 0.1 and W > 3 GeV
        # rather than on y -- and a multiplicity quoted against theirs is a
        # comparison of anything only if both sides use the same region.
        self.x_min = None if x_min is None else float(x_min)
        self.note = note

    def passes(self, q2, y, w2=None, x=None):
        """Leptonic invariants only, as the original test.

        Kept identical in form to `Q2 < Q2_MIN or y < Y_MIN or y > Y_MAX`, so
        the boundary treatment does not shift.

        w2 is OPTIONAL and only consulted by a selection that asks for it, so
        a caller that does not compute W^2 is unaffected.  A selection that
        DOES set w2_min and is handed w2=None raises rather than silently
        passing every event: a cut that quietly does not apply is exactly the
        failure mode CONVENTIONS.md rule 2 is about.
        """
        if self.w2_min is not None:
            if w2 is None:
                raise ValueError(
                    f"selection {self.name!r} cuts on W2 > {self.w2_min:g} but "
                    "was handed w2=None -- the caller must compute W^2")
            if w2 < self.w2_min:
                return False
        if self.x_min is not None:
            if x is None:
                raise ValueError(
                    f"selection {self.name!r} cuts on x > {self.x_min:g} but "
                    "was handed x=None -- the caller must compute x")
            if x < self.x_min:
                return False
        return not (q2 < self.q2_min or y < self.y_min or y > self.y_max)

    def passes_lepton(self, e_lep, theta):
        """The scattered-lepton energy and angle requirements of a tier."""
        if self.e_lep_min is not None and e_lep < self.e_lep_min:
            return False
        if self.theta_max is not None and theta > self.theta_max:
            return False
        if self.theta_min is not None and theta < self.theta_min:
            return False
        return True

    def passes_hadrons(self, n05, n01, dphi_x):
        """Tier E's requirements on the charged hadronic final state.

        dphi_x is measured against the VECTOR SUM of the charged hadrons, which
        is what the emulsion analysis' back-to-back topology cut means.  It is
        NOT the existing `dphi` observable, which is the minimum over
        individual hadrons; the two differ and conflating them would silently
        change the selection.
        """
        if self.n05_min is not None and n05 < self.n05_min:
            return False
        if self.n01_min is not None and n01 < self.n01_min:
            return False
        if self.dphi_x_min is not None and (dphi_x is None
                                            or dphi_x < self.dphi_x_min):
            return False
        return True

    def passes_dimuon(self, n_os_mu):
        """The opposite-sign second-muon requirement.

        n_os_mu counts muons of charge OPPOSITE to the primary that already
        pass passes_lepton() -- i.e. the same E and theta cuts the primary
        had to satisfy, which is what "a second muon satisfying the same
        selection criteria" means.

        In muon NC the primary is a mu-, so the partner is a mu+; in neutrino
        CC the primary is also a mu-, so again a mu+.  Either way the pair is
        opposite-sign, and at these energies an opposite-sign muon comes
        overwhelmingly from the semileptonic decay of a charm hadron.  That
        makes this tier a CHARM-TAGGED selection: it is the classic dimuon
        signal, and its rate is set by the charm fraction times a semileptonic
        branching ratio times the acceptance for a decay muon above 100 GeV.
        Expect it to be small and, on the smaller samples, statistics limited.
        """
        return not self.dimuon or n_os_mu >= 1

    @property
    def has_hadronic_cuts(self):
        return any(v is not None for v in
                   (self.n05_min, self.n01_min, self.dphi_x_min))

    def as_dict(self):
        """Stamped into every result, so a JSON says which region produced it."""
        d = {"name": self.name, "label": self.label, "q2_min": self.q2_min,
             "y_min": self.y_min, "y_max": self.y_max}
        for k in ("e_lep_min", "theta_max", "theta_min", "n05_min",
                  "n01_min", "dphi_x_min", "w2_min", "x_min"):
            v = getattr(self, k)
            if v is not None:
                d[k] = v
        if self.dimuon:
            d["dimuon"] = True
        if self.note:
            d["note"] = self.note
        return d

    def __repr__(self):
        extra = "".join(f", {k}={getattr(self, k):g}" for k in
                        ("e_lep_min", "theta_max", "theta_min", "w2_min",
                         "x_min")
                        if getattr(self, k) is not None)
        return (f"Selection({self.name!r}, Q2 > {self.q2_min:g}, "
                f"{self.y_min:g} < y < {self.y_max:g}{extra})")


# The published region: Q2 > 4 GeV^2, 0.2 < y < 0.9, leptonic invariants only.
INCLUSIVE = Selection("inclusive", Q2_FLOOR, 0.2, 0.9, label="Inclusive selection")

# --------------------------------------------------------------- FASER tiers
# Truth-level proxies for FASER's two published event selections, per
# paper/bench-sec-settings.tex Sec. "FASER-inspired selections".  TWO tiers because the experiment has two detector concepts
# with genuinely different acceptances, so "the FASER cuts" is ambiguous.
#
# BOTH TIERS ARE NESTED INSIDE THE INCLUSIVE REGION (user decision,
# 2026-08-25): they keep Q2 > 4 AND 0.2 < y < 0.9, and add their own cuts on
# top.  A FASER cross-section is therefore always SMALLER than the inclusive
# one and the ratio is a genuine efficiency.
#
# The alternative -- keeping only the perturbativity cut, so each tier's
# energy cut replaces the y window -- was implemented first and rejected: it
# made Tier S three times LARGER than the inclusive region (111 nb against
# 37.5), because dropping y > 0.2 opens the low-y region that the cut exists
# to avoid, and it is precisely where hadronisation models are least reliable.
FASER_S = Selection(
    "faser_s", Q2_FLOOR, 0.2, 0.9, label="FASER Tier S (spectrometer-like)",
    e_lep_min=100.0, theta_max=0.025,
    note="Inclusive region plus E' > 100 GeV and theta < 25 mrad, after "
         "the spectrometer muon-appearance analyses. At 1 TeV the energy "
         "cut coincides with the y < 0.9 already imposed, so the ANGULAR "
         "cut is what acts: it removes the high-Q2, high-y corner.")

# A CHARGED TRACK in the emulsion needs E > TRACK_E_MIN in the lab, as in
# FASER's reconstruction.  Every Tier E multiplicity count applies it (the
# n05/n01 requirements below); until 2026-09-19 Tier E counted charged hadrons
# of any energy, which raised its rate by 1.7% (Sherpa) to ~4% (GENIE), i.e.
# generator-dependently -- the emulsion-shape extraction (pp12) always had it.
TRACK_E_MIN = 1.0   # GeV

FASER_E = Selection(
    "faser_e", Q2_FLOOR, 0.2, 0.9, label="FASER Tier E (emulsion-like)",
    e_lep_min=200.0, theta_min=0.005,
    n05_min=5, n01_min=4, dphi_x_min=math.pi / 2.0,
    note="E' > 200 GeV, theta > 5 mrad, at least 5 charged tracks (E > 1 "
         "GeV) with tan(theta) < 0.5 of which at least 4 with tan(theta) < "
         "0.1, and "
         "Delta phi > pi/2 between the lepton and the summed hadron system. "
         "Its efficiency is generator dependent through the multiplicity "
         "requirement, which is the point rather than a nuisance.")

# The DIMUON SIGNAL tier (user request, 2026-08-25): Tier S exactly, plus a
# second muon of OPPOSITE SIGN passing the same E and theta cuts as the
# primary.  Nested inside Tier S, which is nested inside the inclusive region,
# so its efficiency is bounded by one against either.
#
# Physically this is a CHARM TAG rather than another acceptance cut.  The
# primary muon is a mu- in both currents; an opposite-sign mu+ at these
# energies comes overwhelmingly from a charm hadron decaying semileptonically.
# So the tier measures charm through the same channel the neutrino experiments
# historically used, and it is the one selection here whose rate is set by a
# BRANCHING RATIO as well as by the hadronisation model.
#
# THE ctau > 10 mm STABILITY CONVENTION IS WHAT MAKES THIS A CLEAN TAG, and it
# is worth being explicit about why.  The obvious background to a charm-decay
# muon is an in-flight pi+ or K+ decay, and at 100 GeV those would be copious.
# They cannot contribute here: c*tau is 7.8 m for the pion and 3.7 m for the
# kaon, both far above the 10 mm convention, so the benchmark holds them STABLE
# and they never decay to muons at all.  Charm hadrons sit at c*tau ~ 0.1-0.3 mm,
# well below it, so they do decay.  The convention therefore separates the
# signal from its main background by construction -- which also means this tier
# would mean something different under a different stability convention.
#
# Consequence to keep in mind when reading it: the rate is roughly the charm
# fraction (6-15% over this scan) times a ~10% semileptonic branching ratio
# times the acceptance for a decay muon above 100 GeV, so it is orders of
# magnitude below the other tiers and the smaller samples will be statistics
# limited.  An empty result here is a real answer, not a bug -- analyze.py
# writes one rather than dividing by zero.
FASER_DIMUON = Selection(
    "faser_dimuon", Q2_FLOOR, 0.2, 0.9, label="FASER dimuon signal",
    e_lep_min=100.0, theta_max=0.025, dimuon=True,
    note="Tier S (E' > 100 GeV, theta < 25 mrad) plus a second, "
         "opposite-sign muon passing those same two cuts. The opposite-sign "
         "partner tags charm through its semileptonic decay, so this is the "
         "classic dimuon charm signal rather than a further acceptance cut.")

# ------------------------------------------- the arXiv:2506.13889 DIS region
#
# A FOREIGN REGION, USED FOR ONE COMPARISON AND NOTHING ELSE (user, 2026-09-03).
#
# The "Predictions for FASER" tab compares our muon rate against Table 2.1 of
# arXiv:2506.13889.  Their DIS region is Q > 1.65 GeV and W > 2 GeV; ours is
# Q2 > 4 GeV2 and 0.2 < y < 0.9.  Quoting a ratio between the two measures the
# difference in fiducial region far more than it measures the physics: the y
# window alone discards 66% of the muon rate above Q2 > 4 (measured off the
# full-y POWHEG-RES LHE), which is why that comparison has been reading 16-25%.
#
# So THIS SELECTION EXISTS SO THAT ONE COMPARISON CAN BE LIKE-FOR-LIKE, and it
# is deliberately NOT wired into the benchmark: the published fiducial region
# is unchanged, every other result keeps 0.2 < y < 0.9, and the question of
# whether the benchmark itself should move to a W cut is a separate decision
# the user has parked until the dust settles.  Nothing but the FASER-tab muon
# reproduction may read it.
#
# WHY A W CUT IS THE PHYSICALLY BETTER SHAPE, for whenever that decision comes:
# y > 0.2 exists to keep AHADIC away from the low-W, x -> 1 corner (see the
# note on the FASER tiers above, and sherpa/README.md).  But at Q2 > 4 only
# 0.21% of the muon rate lies below W = 2 GeV and 2.6% below W = 3 GeV, so
# y > 0.2 is discarding two thirds of the rate to remove a couple of per cent.
# A W cut removes the same corner and keeps the rest.
DIS_2506 = Selection(
    "dis2506", 1.65**2, 0.0, 1.0, label="arXiv:2506.13889 DIS region",
    w2_min=4.0,
    note="Q > 1.65 GeV and W > 2 GeV, the DIS cuts of arXiv:2506.13889 "
         "Table 2.1, with NO y window. Used only for the FASER-tab muon "
         "reproduction, so that our rate and theirs are measured on the same "
         "region; the benchmark's own fiducial region is unchanged.")

# ------------------------------------------- the arXiv:2402.13318 comparison
#
# A SECOND FOREIGN REGION, AND FOR THE SAME REASON AS DIS_2506 (user,
# 2026-09-04): "in comparing with the numbers of 2402.13318, we should not be
# imposing the 0.2 < y < 0.9 cut at all, but rather going fully inclusive.  We
# can still keep the Q2 > 4 GeV2 cut to make sure MCs are well behaved."
#
# Table I of arXiv:2402.13318 is every CC interaction at every Q2 and every y.
# Dividing OUR fiducial rate by it measures the y window far more than it
# measures the physics -- the window discards 21% of the neutrino rate above
# Q2 > 4 and 66% of the muon rate -- which is precisely the trap the muon-side
# comparison against arXiv:2506.13889 had to be rebuilt to escape.  So this
# region drops the y window and keeps the perturbativity floor.
#
# APPLIED TO BOTH CURRENTS (CONVENTIONS.md rule 2b).  The neutrino side is what the
# 2402.13318 comparison needs; the muon side is owed it, and is the honest
# partner to the dis2506 region already there.
#
# WHAT IT DOES NOT DO: it is NOT the benchmark's region and nothing else may
# read it.  The published fiducial region is still Q2 > 4 GeV2 with
# 0.2 < y < 0.9 on both currents, and the standing question of whether the
# benchmark should move to a W cut instead is the user's to settle later.
#
# THE X REACH IS THE THING TO WATCH.  With no lower y cut, x runs to the 0.999
# cap of the analytic integrators and W2 falls to the proton mass, where the
# calculation carries no target-mass corrections and the hadronisation models
# are least trustworthy.  It is a small piece of the rate -- 0.21% of the muon
# rate above Q2 > 4 lies below W = 2 GeV -- but it is the part of this region
# that a number quoted from it should be read against.
ALL_Y = Selection(
    "ally", Q2_FLOOR, 0.0, 1.0, label="Q2 > 4 GeV2, all y",
    note="The benchmark's Q2 floor with NO y window, for the comparison "
         "against arXiv:2402.13318 Table I, which is fully inclusive in y. "
         "Used only there; the benchmark's own region is unchanged.")

# ------------------------------------------------- the SIDIS comparison regions
# arXiv:2504.05376 (Bonino, Gehrmann, Loechner, Schoenwald, Stagnitto),
# "Identified Hadron Production in Deeply Inelastic Neutrino-Nucleon
# Scattering": charged pion production in (anti-)neutrino SIDIS at NNLO, the
# calculation a FASER single-inclusive pion measurement would be compared
# with.  Their region -- for the ABCMO comparison AND for their Fig. 4, the
# FPF-type prediction at E_nu ~ 300 GeV -- is x > 0.1 and W > 3 GeV, and
# nothing else.
#
# >>> THESE REGIONS CARRY NO y WINDOW (user, 2026-09-07): "we should not use
#     the 0.2 < y < 0.9 cut here, we already have a W cut and this is
#     sufficient". <<<  That is the right call and it is the same move the
#     paper is eventually to make everywhere (TODO.md): the y window is a
#     Sherpa AHADIC workaround rather than physics, and W > 3 GeV protects
#     the thing it was protecting -- with no lower y cut x runs to the
#     integrator cap and W falls to the proton mass, which is exactly what a
#     W floor forbids.  What stays is the benchmark's Q2 > 4 GeV2, the
#     perturbativity floor, which is also what keeps these regions inside
#     the generation cuts of the samples that can serve them.
#
# >>> ONLY POWHEG AND GENIE CAN SERVE THEM. <<<  Herwig generates with
#     Miny/Maxy = 0.2/0.9 on its LeptonicDISCut and its
#     {Neutral,Charged}CurrentCut, and Sherpa with an INEL selector of the
#     same window, on BOTH currents -- so those samples contain no event
#     outside it and a no-y-window number from them would be silently
#     truncated, not merely uncertain.  POWHEG-V2/RES cut on q2cut alone and
#     GENIE is fully inclusive.  analysis/faser_pions.py measures the
#     fraction of each sample below y = 0.2 and analysis/faser_sidis.py
#     refuses to quote a generator that has none, rather than drawing it.
#
# >>> AND THE FINAL RATES CARRY NO x CUT EITHER (user, 2026-09-07): "we don't
#     cut in x, it is fine to account for this when comparing with their
#     numbers but for the final event rate predictions we should not cut in
#     x". <<<  So the yields come from sidis_e, which has none; sidis_ex
#     exists only to measure what their x floor would cost, and sidis_paper
#     only to compare with their curves on their own region.
SIDIS_E_W3 = Selection(
    "sidis_e", Q2_FLOOR, 0.0, 1.0, label="FASER Tier E + W > 3 GeV, no y cut",
    e_lep_min=200.0, theta_min=0.005,
    n05_min=5, n01_min=4, dphi_x_min=math.pi / 2.0, w2_min=9.0,
    note="The region the pion yields in z are quoted for (user, 2026-09-07): "
         "the emulsion tier with the benchmark's Q2 floor and W > 3 GeV in "
         "place of the y window, and NO cut on x. The companion z > 0.1 is a "
         "cut on the hadron and is applied by binning.")

SIDIS_E_W3_X = Selection(
    "sidis_ex", Q2_FLOOR, 0.0, 1.0,
    label="FASER Tier E + W > 3 GeV + x > 0.1",
    e_lep_min=200.0, theta_min=0.005,
    n05_min=5, n01_min=4, dphi_x_min=math.pi / 2.0, w2_min=9.0, x_min=0.1,
    note="Diagnostic only: what the x > 0.1 of arXiv:2504.05376 would cost a "
         "FASER measurement. No published yield here uses it.")

SIDIS_W = Selection(
    "sidis_w", Q2_FLOOR, 0.0, 1.0, label="Q2 > 4 GeV2 + W > 3 GeV",
    w2_min=9.0,
    note="Our own region with no detector tier and no x cut: the fiducial "
         "region the paper is eventually to adopt (TODO.md). The pion "
         "multiplicity is quoted here as well as on the paper's region, so "
         "the effect of their x floor on a multiplicity is a measured "
         "number.")

SIDIS_PAPER = Selection(
    "sidis_paper", Q2_FLOOR, 0.0, 1.0,
    label="Q2 > 4 GeV2 + W > 3 GeV + x > 0.1",
    w2_min=9.0, x_min=0.1,
    note="The region of arXiv:2504.05376 (x > 0.1, W > 3 GeV) with no "
         "detector tier and no y window, for the like-for-like comparison of "
         "the pion multiplicity dM/dz against their Fig. 4. The one thing it "
         "still carries that theirs does not is the benchmark's Q2 > 4 GeV2 "
         "floor, which is stated wherever the comparison is made.")

# --------------------------------------------- the FINAL region (paper plots)
# >>> THE BENCHMARK REGION OF THE FINAL PAPER (user, 2026-09-11/13): <<<
#     tungsten, Q2 > 4 GeV2 and W > 3 GeV, NO y cut at all.
#
# The same cuts as sidis_w, under its OWN NAME, and that is the point: the
# histos_*_sidis_w.json files already on disk are PROTON results from the
# W-cut study, and a tungsten result written under that name would
# replace them with a plausible number and no error.  Every result is
# also tagged with its target (target.tag: _p, _n, _W), so the name alone
# says which nucleus it is.
#
# Unlike sidis_w this region is served by FRESH samples generated inside it
# (PRODUCTION.md) -- Herwig and Sherpa included, whose earlier samples carry the
# y window as a generation cut and could never serve it.
Q4W3 = Selection(
    "q4w3", Q2_FLOOR, 0.0, 1.0, label="Q2 > 4 GeV2, W > 3 GeV, no y cut",
    w2_min=9.0,
    note="The final benchmark region (user, 2026-09-13): tungsten by isospin "
         "(74p + 110n, NNPDF4.0, no nuclear effects), generated inside the "
         "region. Paper plots.")

# ---------------------------------------- the CHARM region (paper plots)
# >>> CHARM-PRODUCTION OBSERVABLES USE W > 5 GeV (user, 2026-09-14). <<<
# In the neutral current a massless charm-initiated event leaves its cbar in
# the target remnant, and Pythia and Herwig cannot build a hadronic final
# state below W ~ 4.4 GeV: they discard those events, depleting the POWHEG-RES
# charm-tagged rate by 3.3% (1 TeV) and 7.3% (400 GeV).  At W > 5 the loss is
# -0.03..-0.36% (powheg/production/lost_events.py --w2min 25, PRODUCTION.md).
# Applied to BOTH currents (rule 2b) and to the inclusive rates that
# normalise a charm fraction, so numerator and denominator share one region.
# A SUBSET of q4w3, so every sample (generated for q4w3) serves it.
Q4W5 = Selection(
    "q4w5", Q2_FLOOR, 0.0, 1.0, label="Q2 > 4 GeV2, W > 5 GeV, no y cut",
    w2_min=25.0,
    note="The charm-production region of the final paper (user, "
         "2026-09-14): the benchmark region with W raised from 3 to 5 GeV, "
         "which removes the heavy-quark remnant threshold the showered "
         "generators cannot hadronise. Charm figures and their inclusive "
         "normalisations only.")

# --------------------------------- the FASER tiers NESTED IN THE FINAL region
# User ruling 2026-09-14 (paper Sec. 3): "Tier S, Tier E and the dimuon
# selection are nested in the region" -- the dimuon one in the charm
# region W > 5 (below) -- the earlier tiers above keep
# 0.2 < y < 0.9 and stay as they are for the earlier results; these carry W > 3 and
# no y window instead, with the tier cuts unchanged.  Subsets of q4w3.
Q4W3_FASER_S = Selection(
    "q4w3_faser_s", Q2_FLOOR, 0.0, 1.0,
    label="FASER Tier S (spectrometer-like), final region",
    e_lep_min=100.0, theta_max=0.025, w2_min=9.0,
    note="Q2 > 4 GeV2 and W > 3 GeV, plus E' > 100 GeV and theta < 25 mrad.")

Q4W3_FASER_E = Selection(
    "q4w3_faser_e", Q2_FLOOR, 0.0, 1.0,
    label="FASER Tier E (emulsion-like), final region",
    e_lep_min=200.0, theta_min=0.005,
    n05_min=5, n01_min=4, dphi_x_min=math.pi / 2.0, w2_min=9.0,
    note="Q2 > 4 GeV2 and W > 3 GeV, plus E' > 200 GeV, theta > 5 mrad, "
         "at least 5 charged tracks (E > 1 GeV) with tan(theta) < 0.5 of which "
         "at least 4 with tan(theta) < 0.1, and Delta phi > pi/2 between the lepton "
         "and the summed hadron system.")

# THE DIMUON TIER IS A CHARM TAG, so it takes the CHARM region's W > 5 GeV
# (user, 2026-09-14: "for dimuon tier we should also use W > 5 GeV for
# consistency") -- hence q4w5 in its name, not q4w3.  Still nested in q4w3.
Q4W5_FASER_DIMUON = Selection(
    "q4w5_faser_dimuon", Q2_FLOOR, 0.0, 1.0,
    label="FASER dimuon signal, charm region (W > 5 GeV)",
    e_lep_min=100.0, theta_max=0.025, dimuon=True, w2_min=25.0,
    note="Tier S with W raised to 5 GeV, as for every charm observable, plus "
         "a second, opposite-sign muon passing the same energy and angle "
         "cuts.")

SELECTIONS = {"inclusive": INCLUSIVE, "faser_s": FASER_S, "faser_e": FASER_E,
              "faser_dimuon": FASER_DIMUON, "dis2506": DIS_2506,
              "ally": ALL_Y, "sidis_e": SIDIS_E_W3, "sidis_ex": SIDIS_E_W3_X,
              "sidis_w": SIDIS_W, "sidis_paper": SIDIS_PAPER, "q4w3": Q4W3,
              "q4w5": Q4W5, "q4w3_faser_s": Q4W3_FASER_S,
              "q4w3_faser_e": Q4W3_FASER_E,
              "q4w5_faser_dimuon": Q4W5_FASER_DIMUON}

DEFAULT = "inclusive"


def get(name=None):
    """The named selection, defaulting to $BENCH_SELECTION then to inclusive."""
    name = name or os.environ.get("BENCH_SELECTION") or DEFAULT
    try:
        return SELECTIONS[name]
    except KeyError:
        raise SystemExit(
            f"unknown selection {name!r}; known: {', '.join(sorted(SELECTIONS))}"
        ) from None


if __name__ == "__main__":
    for s in SELECTIONS.values():
        print(s)
