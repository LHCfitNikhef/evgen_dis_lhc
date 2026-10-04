# `hienergy` — raising GENIE's validity ceiling

One file, one changed line: `GVLD-Emax` 1000 → 10000 GeV in `CommonParam.xml`,
otherwise a copy of the GENIE install's.

## Why it exists

`gmkspl -e <E>` does not build a spline beyond the tune's declared validity.
GENIE clamps the range and says so only in a WARN line:

    WARN GEVGDriver : Refusing to exceed validity range: Emax = 1000

The file is still written under the requested name, ending at 1000 GeV.  A
beam at or above the last knot then reads 0 from the spline and `gevgen`
hangs in "Could not select interaction".  Always check the last `<E>` knot of
a new spline (`genie/production/genie_job.sh` does).

## Where it is used

- `mu_grv` in `genie/production/genie_job.sh` (`hienergy:p8`): `p8` carries no
  `CommonParam.xml`, so without this the 2000 GeV spline needed for a 1 TeV
  beam would be clamped at 1000.  The benchmark still refuses G18_02a beams
  above 1 TeV.
- The FASER ladders up to 6.8 TeV (`tools/genie_faser_splines.sh`,
  `tools/genie_dimuon_ladder.sh`).  There G18_02a — a Bodek–Yang tune its
  authors declare valid to 1000 GeV — is an **extrapolation**, done because
  it is the configuration FASER uses.

Never put it in front of the HEDIS overlay: its `CommonParam.xml` shadows
`p8/GHE19_00c/` and `gmkspl` fails with "No Common parameters available for
Param list HEDIS-SF".  HEDIS needs no override (validity to 10¹² GeV).
`ultrahigh/` is the same override at 10⁹ GeV for the neutrino-telescope
ladders.
