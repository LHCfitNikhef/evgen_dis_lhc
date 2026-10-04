#!/usr/bin/env python3
"""Herwig cross-section bookkeeping, shared by the make_histos*.py scripts.

WHICH NUMBER IS THE CROSS-SECTION OF A HERWIG SAMPLE
====================================================

Herwig's end-of-run table quotes two totals:

    Total (from attempted events): including vetoed events        38.3(1)e+00
    Total (from generated events):     25000     25271            37.8(1)e+00

They differ by the survival fraction (25000/25271 here, 1.1%).  Events are
lost during showering/hadronisation -- "Can't put the remnant on-shell in
HwRemDecayer::setRemMasses()" and "Failed to generate the shower after 100
attempts in ShowerHandler::cascade()" -- and each failure is silently
discarded and regenerated, so the file still contains the full N events.

**Use the ATTEMPTED (= integrated) number.**  The vetoed events are not a
physical acceptance: they passed the generation cuts and were then dropped
because Herwig could not finish them.  In this benchmark the generation cuts
ARE the fiducial cuts (DIS-mu.in sets Q2 > 4, 0.2 < y < 0.9, and the leptons
do not radiate, so hard-level y is the y the analysis reconstructs), which
means every vetoed event was a FIDUCIAL event.  Normalising with the
"generated events" number therefore reports a fiducial cross-section that is
low by exactly the veto rate.

Cross-check that settles it: the ME-level sample (DIS-mu-ME.in, shower and
hadronisation off) has NO veto at all -- attempts == generated -- and closes
at 0.9975 x YADISM LO.  The hadron-level sample normalised the old way gave
0.9861; normalised this way it gives 0.9977.  The shower does not change the
fiducial rate (POWHEG shows 0.03%), so the two must agree, and only the
"attempted" convention makes them agree.

CAVEAT — THIS FIXES THE RATE, NOT THE SHAPE
===========================================

The vetoed events are NOT a random subset.  Measured on the neutrino sample
by comparing the delivered hadron-level initial-state composition against the
veto-free ME-level one (NEUTRINO_NOTES.md), the discard rate is strongly
flavour dependent:

    d 0.1%   |   ubar 3.5%   |   s 7.4%   |   cbar 17.4%

which is physical -- a valence d costs the remnant nothing, while a sea quark
forces its flavour partner into the remnant, hardest of all for charm -- and
closes on the measured overall rate exactly.  The losses also concentrate at
high x, where little energy is left for the remnant.

So using the integrated cross-section restores the TOTAL correctly, because
the missing events belong in it.  It does NOT repair the delivered sample:
Herwig's charm-producing events are still depleted ~4% net, its D-meson
yields are correspondingly low, and its bin-by-bin ratio still falls to
0.85-0.90 in the highest x bins while its total is right.  There is no card
setting that removes this; it is a limitation of Herwig's remnant treatment
at fixed-target energies and has to be stated wherever Herwig shapes are used.

>>> FIXED FOR  (2026-09-19): REBUILD, DO NOT REPLACE. <<<  The paragraph
above is the earlier state.  In's region (no y cut, W down to 3 GeV) the same
discards cost Herwig 5-26% of its charm (pp06/pp07).  The mechanism: after a
failed remnant, ThePEG draws a NEW phase-space point, so the loss is not a
random subset.  patches/thepeg-2.3.0-event-error-retry.diff makes ThePEG
rebuild the event on the SAME point (sampled runs) or the SAME Les Houches
event (LHE showering) up to MaxEventErrorRetries times; the cards set 100.
They also set RemnantDecayer:DISRemnantOption NoLepton, so the scattered
lepton is never moved to absorb the remnant's recoil (that migrated high-x
charm into x ~ 0.1-0.2 by up to 58%); such an event fails and is retried.
Survival is then 99.4-99.9% and the delivered flavour composition matches the
matrix-element level.  The attempted-events normalisation above still holds
for what is left.  read_job()/read_lhe_job() REFUSE a job whose .out does
not carry the retry line, so a pre-fix sample cannot be analysed by mistake.

The .out file quotes the integrated total to 3 significant figures only, so
the precise value and its error are read from the .log's sampler summary
("Total integrated xsec" / "error in xsec"); the .out is used for the event
counts.  This is the rule analysis/analyze_nu.py already applied on the
neutrino side; these helpers exist so that both sides state it once.

COMBINING JOBS
==============

Each job is an independent estimate of the SAME integral, so they are
combined by inverse variance, not by a plain mean.  A plain mean is only
correct for jobs of equal precision, and the drivers permit partial reruns
and different event counts.  The fiducial fraction is pooled over jobs.
"""
import glob
import math
import os
import re
from typing import List, NamedTuple

# ".out": "Total (from generated events):  25000  25271  3.78(1)e+01"
RE_OUT_GEN = re.compile(
    r"Total \(from generated events\):\s+(\d+)\s+(\d+)\s+"
    r"(-?[0-9.]+)\(([0-9]+)\)e([+-]\d+)")
# ".out": "Total (from attempted events): including vetoed events  0.108(1)e+03"
# Same quantity as the sampler integral, but only 3 significant figures.
# NOTE the -? and the anchoring: a lazy [^\n]*? in front of an unsigned
# number will happily skip over a minus sign, silently turning a negative
# cross-section into a positive one.  Matchbox NLO jobs CAN come out negative
# when the weight tail is bad, so this matters.
RE_OUT_ATT = re.compile(
    r"Total \(from attempted events\)[^\n]*?\s"
    r"(-?[0-9.]+)\(([0-9]+)\)e([+-]\d+)\s*$", re.M)
# ".log": the ThePEG sampler's integral over the generation cuts, full precision
# ".out": printed by the patched ThePEG when MaxEventErrorRetries > 0
RE_RETRY = re.compile(
    r"evgen-benchmark event-error retry \(MaxEventErrorRetries (\d+)\): "
    r"events retried (\d+), retries (\d+), recovered (\d+), abandoned (\d+)")


def _require_retry(outfile: str, txt: str) -> int:
    """Abandoned-event count of a job; REFUSE a job run without retries.

    job directories are herwig7/v2<cur>pwg[neg]_<t>_job... (sampled) and
    powheg/hw_v2<cur>_<t>_job_N (the POWHEG LHE showered by Herwig).  See the
    module docstring: without the retry a  Herwig sample is charm-depleted.
    """
    m = RE_RETRY.search(txt)
    d = os.path.basename(os.path.dirname(os.path.abspath(outfile)))
    if m is None and (d.startswith("v2") or d.startswith("hw_v2")):
        raise RuntimeError(
            f"{outfile}: a production Herwig job without the event-error retry line.  "
            f"It was produced before patches/thepeg-2.3.0-event-error-retry."
            f"diff and discards remnant failures by drawing NEW points, which "
            f"depletes charm by 5-26%.  Regenerate it (herwig7/production/production.sh"
            f" or powheg/production/herwig_arm.sh); do not analyse it.")
    return int(m.group(5)) if m else -1


RE_LOG_INT = re.compile(
    r"Total integrated xsec:\s+([0-9.eE+-]+)\s*\n\s*error in xsec:\s+"
    r"([0-9.eE+-]+)")


def _bracket_value_nb(mant: str, err_digits: str, expo: int):
    """Decode Herwig's "3.78(1)e+01" notation: error is in the last digits."""
    val = float(mant) * 10.0**expo
    ndec = len(mant.split(".")[1]) if "." in mant else 0
    err = float(err_digits) * 10.0**(expo - ndec)
    return val, err


class JobXsec(NamedTuple):
    """One Herwig job's cross-section bookkeeping.  All sigmas in pb."""
    path: str
    sigma_pb: float        # integrated over the generation cuts (USE THIS)
    err_pb: float          # its integration error
    n_written: int         # events actually in events.hepmc
    n_attempted: int       # incl. events vetoed during shower/hadronisation
    sigma_delivered_pb: float   # Herwig's "from generated events" (NOT used)
    n_abandoned: int = -1  # events dropped after MaxEventErrorRetries; -1 = no retry

    @property
    def survival(self) -> float:
        return self.n_written / self.n_attempted if self.n_attempted else 1.0


def read_job(outfile: str) -> JobXsec:
    """Parse one job's cross-section from its .out (counts) and .log (value)."""
    with open(outfile) as f:
        txt = f.read()
    m = RE_OUT_GEN.search(txt)
    if not m:
        raise RuntimeError(f"no end-of-run cross-section table in {outfile}")
    n_written, n_attempted = int(m.group(1)), int(m.group(2))
    n_abandoned = _require_retry(outfile, txt)

    mant, expo = m.group(3), int(m.group(5))
    sigma_delivered_nb = float(mant) * 10.0**expo

    # Preferred source: the sampler summary in the .log, full precision.
    sigma_nb = err_nb = None
    logfile = outfile[:-4] + ".log"
    if os.path.exists(logfile):
        with open(logfile) as f:
            mi = RE_LOG_INT.search(f.read())
        if mi:
            sigma_nb, err_nb = float(mi.group(1)), float(mi.group(2))

    if sigma_nb is None:
        # Fallback: the .out's "attempted events" row.  Same quantity, but
        # only 3 significant figures.  Matchbox (NLO) runs write no sampler
        # summary to the .log, so this is the only source there.
        ma = RE_OUT_ATT.search(txt)
        if not ma:
            raise RuntimeError(
                f"no integrated cross-section in {logfile} or {outfile}.  "
                f"Do NOT fall back to the 'generated events' number -- it is "
                f"scaled down by the shower veto; see this module's docstring.")
        sigma_nb, err_nb = _bracket_value_nb(ma.group(1), ma.group(2),
                                             int(ma.group(3)))

    return JobXsec(path=outfile, sigma_pb=sigma_nb * 1e3, err_pb=err_nb * 1e3,
                   n_written=n_written, n_attempted=n_attempted,
                   sigma_delivered_pb=sigma_delivered_nb * 1e3,
                   n_abandoned=n_abandoned)


# ------------------------------------------- the Les Houches event handler
# A run that SHOWERS an external Les Houches file prints a different table:
# there is no sampler and no "attempted events" row, because ThePEG did not
# integrate anything.  One line carries everything -- events generated,
# events attempted (= Les Houches events consumed, the difference being
# Herwig's shower discards) and the cross-section it delivers:
#
#     Total:      57846     60000    0.2114(5)e+03
#
# THIS NUMBER IS THE DELIVERED CROSS-SECTION, NOT THE INTEGRATED ONE, and
# that is the point of reading it.  The sample's normalisation is pinned to
# POWHEG's integrator, exactly as the Pythia-showered arm's is; what the
# analysis needs from Herwig is the INDEPENDENT measurement to check that
# pinning against.  Reading Herwig's own delivered sigma is the only closure
# available here, because the HepMC weights ThePEG writes are normalised to
# the largest weight in the file and carry no absolute scale at all -- a mean
# weight over them is a pure number, and comparing it to a cross-section in
# picobarns produces "-100%", which is what this reader exists to replace.
RE_LHE_TOTAL = re.compile(
    r"^Total:\s+(\d+)\s+(\d+)\s+([0-9.]+)\((\d+)\)e([+-]\d+)",
    re.MULTILINE)


def read_lhe_job(outfile: str) -> JobXsec:
    """Parse one POWHEG-LHE-showered Herwig job's .out file.

    sigma_pb here is Herwig's DELIVERED cross-section (see above), so the
    field name is shared with read_job() but the meaning is not -- which is
    why this is a separate function rather than a flag on that one.
    """
    with open(outfile) as f:
        txt = f.read()
    m = RE_LHE_TOTAL.search(txt)
    if not m:
        raise RuntimeError(
            f"no Les Houches event handler total in {outfile}.  A run that "
            f"showers an external file prints 'Total: <generated> "
            f"<attempts> <sigma>'; if that line is missing the run did not "
            f"finish and the sample must not be normalised.")
    n_written, n_attempted = int(m.group(1)), int(m.group(2))
    n_abandoned = _require_retry(outfile, txt)
    sigma_nb, err_nb = _bracket_value_nb(m.group(3), m.group(4),
                                        int(m.group(5)))
    return JobXsec(path=outfile, sigma_pb=sigma_nb * 1e3, err_pb=err_nb * 1e3,
                   n_written=n_written, n_attempted=n_attempted,
                   sigma_delivered_pb=sigma_nb * 1e3, n_abandoned=n_abandoned)


def combine_lhe_jobs(outfiles) -> tuple:
    """(sigma_delivered_pb, n_written, n_attempted) over several such jobs.

    WEIGHTED BY ATTEMPTS, not a plain mean over jobs: each job's sigma is
    sum(weights)/attempts, so the pooled estimate is the attempt-weighted
    average.  CONVENTIONS.md rule 2 names the plain mean over jobs of unequal size
    as a thing not to do.
    """
    tot_w = tot_a = 0
    acc = 0.0
    for fn in sorted(outfiles):
        j = read_lhe_job(fn)
        acc += j.sigma_delivered_pb * j.n_attempted
        tot_w += j.n_written
        tot_a += j.n_attempted
    if not tot_a:
        return None, 0, 0
    return acc / tot_a, tot_w, tot_a


def read_jobs(pattern: str, verbose: bool = True) -> List[JobXsec]:
    """Read every job matching a glob of .out files, newest-consistent order."""
    outs = sorted(glob.glob(pattern))
    if not outs:
        raise SystemExit(f"no Herwig .out files matching {pattern}")
    jobs = [read_job(o) for o in outs]
    if verbose:
        # neutrino CC sits at ~3 pb and muon NC at ~38 nb: pick the unit that
        # keeps the error visible rather than printing "+- 0.0000"
        nb = all(j.sigma_pb >= 1e3 for j in jobs)
        scale, unit = (1e3, "nb") if nb else (1.0, "pb")
        for j in jobs:
            print(f"  {os.path.relpath(j.path)}: sigma = {j.sigma_pb/scale:.4f}"
                  f" +- {j.err_pb/scale:.4f} {unit}, {j.n_written} events, "
                  f"survival {j.survival:.4f}")
    return jobs


def combine(jobs: List[JobXsec]):
    """Inverse-variance combination of the per-job integrals.

    Returns (sigma_pb, err_pb, n_written_total, survival).  Falls back to an
    n_written-weighted mean if any job failed to quote an error.
    """
    n_written = sum(j.n_written for j in jobs)
    n_attempted = sum(j.n_attempted for j in jobs)
    survival = n_written / n_attempted if n_attempted else 1.0

    if all(j.err_pb > 0 for j in jobs):
        wts = [1.0 / j.err_pb**2 for j in jobs]
        sigma = sum(w * j.sigma_pb for w, j in zip(wts, jobs)) / sum(wts)
        err = math.sqrt(1.0 / sum(wts))
    else:  # no usable errors: weight by the statistics each job carries
        tot = sum(j.n_written for j in jobs)
        sigma = sum(j.n_written * j.sigma_pb for j in jobs) / tot
        err = 0.0
    return sigma, err, n_written, survival


def fiducial_sigma(sigma_pb: float, err_pb: float, n_fid: int, n_written: int):
    """sigma_fid and its error for an UNWEIGHTED sample.

    sigma_fid = sigma_gen_region * (n_fid / n_written).  Two uncertainties
    enter and they are independent: the integration error on sigma_gen_region,
    and the BINOMIAL error on the fiducial fraction.  The fraction is binomial,
    not Poisson: when it is close to 1 (generation cuts == fiducial cuts) the
    Poisson form sqrt(n_fid) overstates the error by more than an order of
    magnitude, because n_fid cannot fluctuate above n_written.
    """
    f = n_fid / n_written
    sigma_fid = sigma_pb * f
    rel_int = (err_pb / sigma_pb) if sigma_pb else 0.0
    var_f = f * (1.0 - f) / n_written
    rel_f = math.sqrt(var_f) / f if f > 0 else 0.0
    return sigma_fid, sigma_fid * math.hypot(rel_int, rel_f)
