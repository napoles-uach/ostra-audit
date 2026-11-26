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
| status | Known issue; data contributed |
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

