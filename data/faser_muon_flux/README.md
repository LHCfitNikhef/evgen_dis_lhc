# The LHC muon flux at FASERν

Vendored from **LHCfitNikhef/DIS_with_LHC_muons**, `muon_fluxes/`, the data
release of arXiv:2506.13889.  Kept in the repository (data, not code) so a
clone can reproduce the muon event rates.

## Files

Each flux is an **LHAPDF grid**, the form in which it is handed to POWHEG-RES
(`powheg/run_powheg_muflux.sh`):

    muon_flux_FASERv_Run3_var1{.info,_0000.dat}   r < 9 cm (33% of the flux)
    muon_flux_FASERv_Run3_var2{.info,_0000.dat}   full FASERν area, 25 × 30 cm (used here)

Flavours 13 and −13 are μ⁻ and μ⁺.  The grid variable is
x_μ = 2E_μ/√s_pp with √s_pp = 13.6 TeV (x_μ = 1 is 6.8 TeV); 22 nodes from
8.1 GeV.  The Q dependence is flat, and as in every LHAPDF grid the stored
numbers are x·f(x).

`lhapdf/` lays the same grids out as LHAPDF looks them up (relative
symlinks + `pdfsets.index`), for use with `LHAPDF_DATA_PATH` prepended.  It
also holds `muon_flux_FASERv_Run3_var2_sum` (μ⁻ + μ⁺ in one flavour, written
by `tools/make_muon_flux_sum.py`), because POWHEG-RES reads a single beam charge
and photon exchange at Q² ≪ M_Z² does not distinguish them.

## Normalisation

Eq. (2.1) of arXiv:2506.13889 defines f_μ(x_μ) = n_T L_T dN_μ/dx_μ in pb⁻¹,
so **the grid already carries the target**: n_T is the tungsten nucleon
density and L_T = 50 cm (the fiducial length the paper uses, shorter than
FASERν), and N = ∫ dx_μ f_μ σ with σ in pb.  The published Table 2.1 of that
paper corresponds instead to the full FASERν (1.1 t, 76 cm of tungsten), a
factor 1.520; `analysis/faser_muflux_check.py` computes both.  The built-in
luminosity is **250 fb⁻¹**, not the 150 fb⁻¹ of `data/faser_flux/`.
Integrated: 1.90 pb⁻¹ for 25 × 30 cm (μ⁻ 1.06, μ⁺ 0.84), 0.63 pb⁻¹ for
r < 9 cm.

## Provenance and citation

FLUKA simulations of the muon flux reaching FASER, interpolated into LHAPDF
grids by R. Francener, V. P. Goncalves, F. Kling, P. Krack and J. Rojo,
"Deep-Inelastic Scattering at TeV Energies with LHC Muons",
[arXiv:2506.13889](https://arxiv.org/abs/2506.13889).  Their caveat: the
flux changes every data-taking year, and the geometry is that of the start of
Run 3, with the detector axis at (x, y) = (1, −3.3) cm from the nominal line
of sight.
