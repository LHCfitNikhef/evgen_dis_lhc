#!/usr/bin/env python3
"""Where the paper-plots samples live, and what their results are called.

ONE PLACE, and the drivers follow it rather than the other way round: the
layout is written down in PRODUCTION.md, and every analysis branch that
reads a sample asks this module for the path.  A second copy of a sample
path is how this benchmark has repeatedly analysed the wrong directory
(POWHEG_NU_LHE_DIR, 2026-08-28), so there is not one.

The samples are generated in the final region (selection q4w3) on ONE
nucleon at a time, p or n.  An analysis of one therefore needs BOTH knobs set
explicitly, and refuses otherwise:

    BENCH_SELECTION=q4w3 BENCH_TARGET=n BENCH_ENERGY=700 \\
        analysis/analyze.py v2_powheg

A sample analysed under another selection would be a proton-style result
from a tungsten-style sample, and one analysed without a target would be
written under a name that does not say which nucleon it is.
"""
import os

import beams
import paths
import selection
import target

REGION = "q4w3"
# The selections a sample may be analysed under: the generation region
# itself and the regions NESTED in it -- the charm region (W > 5, user
# 2026-09-14) and the FASER tiers (user ruling, paper Sec. 3).  Each is
# checked to be a subset of q4w3 below, so a looser one cannot be added here
# by mistake and silently read outside what the samples contain.
BENCH_SELECTIONS = ("q4w3", "q4w5", "q4w3_faser_s", "q4w3_faser_e",
                 "q4w5_faser_dimuon")


def _nested_in_region(sel):
    base = selection.SELECTIONS[REGION]
    return (sel.q2_min >= base.q2_min and sel.w2_min is not None
            and sel.w2_min >= base.w2_min and sel.y_min >= base.y_min
            and sel.y_max <= base.y_max)


def require():
    """(selection, nucleon) for an analysis, refusing anything implicit."""
    sel = selection.get()
    if sel.name not in BENCH_SELECTIONS or not _nested_in_region(sel):
        raise SystemExit(f"the samples are generated in the {REGION} region; "
                         f"set BENCH_SELECTION to one of {BENCH_SELECTIONS} "
                         f"(got {sel.name!r})")
    t = os.environ.get("BENCH_TARGET")
    if t not in ("p", "n"):
        raise SystemExit("a sample is ONE nucleon: set BENCH_TARGET=p or n "
                         "explicitly (the tungsten result is the combination, "
                         "analysis/combine_target.py)")
    if sel.q2_min != selection.Q2_FLOOR_DEFAULT:
        raise SystemExit("the benchmark region is Q2 > 4 on both currents (user, 2026-09-13); "
                         "unset BENCH_Q2MIN")
    return sel, t


def tag(energy):
    return beams.Beams("mu", energy).tag


# ------------------------------------------------------------ sample dirs --
# Job-directory BASES, to be passed through beams.at_energy() and then to
# analyze.job_files(parent, base, fname), which appends _<N>.
def powheg_job_base(current, t, mc=False):
    """mc=True: POWHEG-V2mc, the MASSIVE-charm neutrino entry (charm only),
    powheg/production/run_v2_mc.sh -> powheg/v2numc_<t>_job[_TAG]_<j>."""
    if mc:
        assert current == "nu", "POWHEG-V2mc is charged current only"
        return f"v2numc_{t}_job"
    return f"v2{current}_{t}_job"


def herwig_job_base(current, t, neg=False):
    return f"v2{current}pwg{'neg' if neg else ''}_{t}_job"


def herwig_runname(current, t, energy, neg=False):
    """The Herwig run name, as herwig7/production/make_cards.py writes it:
    V2-<cur>-PWG[NEG]-<t>[-<tag>], 1 TeV untagged.  The .out files a job
    leaves are <runname>-S<seed>.out, and globbing on the full name is what
    keeps the untagged 1 TeV glob off the tagged energies' files."""
    stem = f"V2-{current}-PWG{'NEG' if neg else ''}-{t}"
    return stem if energy == beams.ANCHOR_ENERGY else f"{stem}-{tag(energy)}"


# GENIE: genie/production/genie_job.sh writes genie/v2g_<cfg>_<t>_job[_TAG]_<N>/.
# The generator key -> the driver's configuration, per current.  The keys
# are the earlier production's result keys, so a figure reads the same names: on the NEUTRINO
# side "genie" is HEDIS/BGR18 and "genie_lo" the classic G18_02a, as in the earlier production.
GENIE_CFG = {
    "mu": {"genie": "mu_grv", "genie_nnpdf": "mu_nnpdf"},
    "nu": {"genie_lo": "nu_grv", "genie_nnpdf": "nu_nnpdf", "genie": "nu_hedis"},
}


def genie_job_base(current, key, t):
    return f"v2g_{GENIE_CFG[current][key]}_{t}_job"


def genie_files(current, key, t, energy, fname):
    """<fname> of every COMPLETE  GENIE job (V2_OK): the driver leaves a
    short job in place for inspection, and a glob must not pick it up."""
    import glob
    import re
    base = beams.at_energy(genie_job_base(current, key, t), energy)
    pat = re.compile(rf"^{re.escape(base)}_\d+$")
    out = []
    for d in sorted(glob.glob(f"{paths.REPO}/genie/{base}_*")):
        if (os.path.isdir(d) and pat.match(os.path.basename(d))
                and os.path.exists(f"{d}/V2_OK")
                and os.path.exists(f"{d}/{fname}")):
            out.append(f"{d}/{fname}")
    return out


def sherpa_rundir(current, t, energy):
    stem = "V2_MuonDIS_NLO" if current == "mu" else "V2_NuDIS_NLO"
    return f"{paths.SHERPA_RUNS}/{beams.at_energy(f'{stem}_{t}', energy)}"


def powheg_rundir(current, t, energy, mc=False):
    """The POWHEG integration directory; its tag is ALWAYS present, as in the earlier production.
    mc=True: POWHEG-V2mc's $POWHEG_V2/v2-numc<TAG>-<t>."""
    root = paths.POWHEG_RES if current == "mu" else paths.POWHEG_V2
    if mc:
        assert current == "nu", "POWHEG-V2mc is charged current only"
        return f"{root}/v2-numc{tag(energy)}-{t}"
    return f"{root}/v2-{current}{tag(energy)}-{t}"


# ---------------------------------------------------------------- results --
def result_key(gen, t, region=REGION):
    """<gen>_<region>_<t>, before the energy tag: e.g. powheg_q4w3_n,
    powheg_charmfinal_q4w5_W.  The region is EXPLICIT: the analyses pass the
    selection they ran under, and a reader asks for the one it wants."""
    if region not in BENCH_SELECTIONS:
        raise SystemExit(f"{region!r} is not a benchmark selection {BENCH_SELECTIONS}")
    return f"{gen}_{region}{target.tag(t)}"


# ------------------------------------------------------------ shower arms --
# Paper plots 9-11 : the SAME 1 TeV POWHEG Les Houches files re-showered
# (powheg/production/qed_arms.sh, powheg/production/herwig_arm.sh) on a SUBSET of the
# production jobs.  Every arm and its Pythia baseline ("sub") must cover the
# SAME job numbers, or a ratio between them carries a fluctuation nobody
# accounted for; arm_files() refuses otherwise.
ARMS = {  # arm -> (directory prefix before "v2<cur>_<t>_job_N", OK marker)
    "sub": ("", "V2_OK"),
    "qedfsr": ("qedfsr_", "V2_OK"),
    "qedfsrisr": ("qedfsrisr_", "V2_OK"),
    "qedfull": ("qedfull_", "V2_OK"),
    "hw": ("hw_", "HW_OK"),
}


def arm_jobnums(current, t):
    """The job numbers of the arm subset: those the QED fsr arm showered."""
    import glob
    import re
    out = []
    for d in glob.glob(f"{paths.REPO}/powheg/qedfsr_v2{current}_{t}_job_*"):
        m = re.fullmatch(rf"qedfsr_v2{current}_{t}_job_(\d+)", os.path.basename(d))
        if m and os.path.exists(f"{d}/V2_OK"):
            out.append(int(m.group(1)))
    if not out:
        raise SystemExit(f"no complete qedfsr_v2{current}_{t}_job_* arm jobs")
    return sorted(out)


def arm_files(current, arm, t):
    """events.hepmc of arm `arm` over the subset; refuses a partial arm."""
    if arm not in ARMS:
        raise SystemExit(f"unknown arm {arm!r}; known {sorted(ARMS)}")
    pre, ok = ARMS[arm]
    if arm == "qedfsrisr" and current == "nu":
        raise SystemExit("no fsrisr arm in the charged current (neutral beam)")
    files = []
    for n in arm_jobnums(current, t):
        d = f"{paths.REPO}/powheg/{pre}v2{current}_{t}_job_{n}"
        if not (os.path.exists(f"{d}/{ok}") and os.path.exists(f"{d}/events.hepmc")):
            raise SystemExit(f"{d}: arm {arm} is not complete on the subset "
                             f"{arm_jobnums(current, t)}")
        files.append(f"{d}/events.hepmc")
    return files


def powheg_offered(current, t, energy, files):
    """LHE events offered to the shower by the jobs behind `files`, from the
    .nevents sidecars: POWHEG-RES seed N is pwgevents-000N.lhe; POWHEG-V2
    batch j is <rundir>[-b<j>]/pwgevents.lhe.  None if a sidecar is missing."""
    import re
    rd = powheg_rundir(current, t, energy)
    tot = 0
    for f in files:
        n = int(re.search(r"_(\d+)$", os.path.basename(os.path.dirname(f))).group(1))
        side = (f"{rd}/pwgevents-{n:04d}.lhe.nevents" if current == "mu" else
                f"{rd if n == 1 else f'{rd}-b{n}'}/pwgevents.lhe.nevents")
        if not os.path.exists(side):
            return None
        tot += int(open(side).read().split()[0])
    return tot

