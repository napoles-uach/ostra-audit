# Search log

Every search, including those that find nothing. This is the data for the paper's comparison between catalog-guided
and free exploration. Effort is agent wall time (approximate), excluding human review.

| Date | Agent | Target | Pattern (catalog) | Oracle | Effort | Outcome |
|---|---|---|---|---|---|---|
| 2026-11-26 | Explorer 1 | XCFun vs libxc in PySCF, all 140 XCFun names (104 comparable) | transcription-slip, inconsistent-configuration | independent-implementation | ~1.5 h | 2 new (OSTRA-027, 028); 4 known re-detected; 3 candidates; 5 artifacts explained |
| 2026-11-26 | Explorer 1 | simple-dftd3 periodic D3 (8 crystals, equivalent cells, gradient and virial vs FD) | convergence-fixed | invariance, energy-gradient consistency | ~0.5 h | 0 (old candidate refuted) |
| 2026-11-26 | Explorer 1 | PySCF PBC Ewald nuclear energy (equivalent cells) | convergence-fixed | invariance | ~0.2 h | 0 |
| 2026-11-26 | Explorer 1 | tblite GFN2-xTB (8 crystals, equivalent cells) | convergence-fixed | invariance | ~0.5 h | 0 (supercell differences = Gamma-point sampling, by design) |
| 2026-11-26 | Explorer 1 | dftd4 gradient/virial vs FD on sheared cells | convergence-fixed | energy-gradient consistency | ~0.2 h | extension of OSTRA-005 (gradient inconsistent with energy in sheared MgO) |
| 2026-11-26 | Explorer 2 | tblite, dftd4, simple-dftd3, mctc-lib optional outputs | implicit-api-contract | api-contract | running | — |
