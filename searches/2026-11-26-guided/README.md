# Guided search, round 1 (2026-11-26)

Experiment for the paper: searches chosen from the catalog of fragile numerical decisions, each logged with target,
pattern, oracle, effort and outcome, including negatives (see `../LOG.md`).

## S1. XCFun vs libxc in PySCF (pattern: transcription slips / inconsistent configuration; oracle: independent implementation)
All 140 XCFun names in PySCF (+ aliases); 104 compared (same xc string, or explicit component pairs), 36 without a
libxc counterpart. Realistic densities: H2O (RKS), OH and stretched N2 (UKS); deriv 1-3; integrated E_xc difference
and max relative difference of v, f, k (rho > 1e-6); hybrid/RSH coefficients as used by PySCF.
- New: TW (OSTRA-027), SPBEC (OSTRA-028).
- Known, re-detected: REVTPSSX/REVTPSS (012), REVTPSSC (013), R4SCANX (014), TF/TFK spin (016).
- Artifacts / not bugs: M05X, M052XX... vs HYB_ libxc components (the XCFun component carries no HF fraction);
  VWN3C vs LDA_C_VWN_3 (XCFun VWN3 is the RPA parametrization; the `VWN3` string agrees in both backends);
  rsh/alpha differences for global hybrids (alpha unused when omega = 0); PBE-family 6.7e-6 (constant precision);
  zvPBEint/zvPBEsol (documented polynomial fit of |zeta|^omega).
- Candidates: P86C (0.08 mEh, rel 2e-4 in exc, also in BP86/B3P86); PW91C (v rel 2e-3, E identical);
  BRX and CSC listed in PySCF but fail to evaluate (need jp / laplacian); libxc d3 not available for r2SCAN family.
## S2. Equivalent-cell invariance (pattern: convergence fixed for the typical case; oracle: invariance)
8 crystals x 6-7 equivalent cells.
- simple-dftd3 periodic D3(BJ): invariant to 1e-12 kJ/mol; gradients and virial vs FD to 1e-11 on triclinic cells
  -> the old candidate "triclinic FD deviations in periodic D3" is refuted.
- PySCF PBC Ewald (energy_nuc): invariant.
- tblite GFN2-xTB: invariant among cells with the same number of atoms (<= 0.01 kJ/mol/f.u.); supercells differ by
  hundreds of kJ/mol because tblite samples only the Gamma point (by design: k-point convergence, not a bug).
- dftd4 (extends OSTRA-005): in the sheared MgO cell the analytic gradient (3.8e-5 Eh/bohr) is not the derivative
  of the energy (FD ~ 0); in sheared CaF2 the forces are consistent but nonzero, breaking the crystal symmetry.
