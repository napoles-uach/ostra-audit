# OSTRA-005: EEQ Ewald sums with a fixed +-2 image range: charges depend on the equivalent unit cell

| Field | Value |
|---|---|
| library | multicharge / dftd4 |
| location | multicharge Ewald kernels (rep=[2,2,2]) |
| versions | dftd4 4.3.1; multicharge 4bb6c96 |
| detected | 2026-11-20 |
| detected_by | Explorer 1 |
| oracle | Invariance under equivalent cells (a' = a + k b) |
| class | Convergence, known |
| impact_level | I (up to 18 % in charges, 1.25 % in E_D4 for skewed cells) |
| status | Known issue; data contributed |
| report | grimme-lab/multicharge#72 (comment with skewed-cell data and LLL workaround) |

Coordination numbers are invariant (1e-15); the discrepancy comes from the charge model. LLL lattice reduction
before the Ewald sum is an effective workaround. Full comment text: `report/issue.md`.
