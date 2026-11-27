# OSTRA-005: EEQ Ewald sums with a fixed +-2 image range: charges depend on the equivalent unit cell

| Field | Value |
|---|---|
| library | multicharge / dftd4 |
| location | multicharge Ewald kernels (rep=[2,2,2]) |
| versions | dftd4 4.3.1; multicharge 4bb6c96 |
| detected | 2026-11-20 |
| detected_by | Explorer 1 |
| oracle | Invariance under equivalent cells (a' = a + k b) |
| libraries | multicharge, dftd4 |
| pattern | convergence-fixed |
| oracle_class | invariance |
| class | Convergence, known |
| impact_level | I for non-reduced cells of compact ionic solids (MgO: 34 % in E_D4, dq = 0.6); III-IV for standard cells (<= 0.05 kJ/mol/f.u.) |
| status | Resolved upstream in multicharge main (aa2daad, PR #81, 2026-11-27; not yet released, dftd4 still bundles the old version) |
| report | Commented on known issue: [multicharge#72](https://github.com/grimme-lab/multicharge/issues/72) |

Coordination numbers are invariant (1e-15); the discrepancy comes from the charge model. LLL lattice reduction
before the Ewald sum is an effective workaround. Full comment text: `report/issue.md`.

## Impact on observables (2026-11-25)
`studies/impact-d4-ewald-2026-11-25/`: 8 crystals, equivalent cells. Standard cells agree to < 0.001 kJ/mol/f.u.
except the MgO primitive cell (0.046 kJ/mol/f.u., 0.15 %). Sheared equivalent cells (a' = a + 2b): MgO -10.2
kJ/mol/f.u. (-34 %), q(Mg) +0.61; CaF2 -2.1 %; LiF -1.9 %. Silent.

## Gradient consistency on sheared cells (guided search, 2026-11-26)
dftd4 4.3.1, MgO a' = a + 2b: analytic gradient 3.8e-5 Eh/bohr while central FD of the energy is ~0 (forces must
vanish by symmetry): the gradient is not the derivative of the energy. Sheared CaF2: forces 5.7e-5 Eh/bohr,
consistent with FD but nonzero, breaking the crystal symmetry. Standard cells: consistent to 1e-12. Virial consistent
in all cases. `searches/2026-11-26-guided/cell_invariance/d4_virial.py`.

## Upstream resolution check (2026-11-27)
multicharge PR #81 (T. Zakharov, merged 2026-11-27, aa2daad) replaces the fixed [2,2,2] image boxes by
tolerance-based real/reciprocal cutoffs. The PR reports supercell charge errors going from 0.12 to 1.1e-13.
Independent check with our invariance oracle (multicharge CLI, rocksalt MgO a = 4.212 A, primitive cell vs
equivalent cells a' = a + 2b and a' = a + 5b - 3c):

| build | primitive q(Mg) | a+2b q(Mg) | a+5b-3c q(Mg) | E_EEQ primitive / a+2b / a+5b-3c (Eh) |
|---|---|---|---|---|
| 4bb6c96 (before #81) | 0.7049 | 1.3140 | -2.9564 | -0.12313 / -0.22953 / +0.51641 |
| aa2daad (after #81) | 0.7011 | 0.7011 | 0.7011 | -0.12247 (all three identical to 14 digits) |

The primitive cell itself was off by 0.004 e before the fix. Residual exposure: dftd4 4.3.x and tblite bundle
the earlier multicharge until they update the subproject. OSTRA-026 patch re-applied on aa2daad: applies cleanly
(offset only), meson test 50 OK, 0 Fail.

