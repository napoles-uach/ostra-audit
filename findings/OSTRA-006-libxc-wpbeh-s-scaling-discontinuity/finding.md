# OSTRA-006: HSE03/06/12, LC-wPBE and wPBEh energy is discontinuous at s = 1 (s_scaling_2)

| Field | Value |
|---|---|
| library | libxc |
| location | maple/scaling.mpl (s_scaling_2), used by maple/gga_exc/gga_x_wpbeh.mpl |
| versions | 7.1.2 and devel 3c7794c50 |
| detected | 2026-11-21 |
| detected_by | Explorer 1 (spin-polarization oracle, by chance) |
| oracle | Spin-polarized vs unpolarized identity; energy scan; FD vs analytic forces |
| class | Known design (Henderson et al. 2009, eq. 8 vs eq. 9) |
| impact_level | II |
| status | Classified; not to be reported (deliberate fidelity to HSE06) |
| report | Not to be reported (known design, Henderson et al. 2009) |

## Cause
eq. (8) of Henderson, Izmaylov, Scalmani, Scuseria, J. Chem. Phys. 131, 044108 (2009) is piecewise: identity for
s <= 1 and s - ln(1 + e^(s - s_max)) above; the second branch is 1 - 5.1e-4 at s = 1. The same paper notes the
discontinuity and proposes eq. (9), s - (1 - e^-s) ln(1 + e^(s - s_max)), which libxc does not offer.

## Measured impact (H2O, def2-SVP, grid level 3, H displaced in 40 steps of 0.0005 A)
- Energy noise around a quartic fit: HSE06 4.6e-8 Eh vs PBE0 1.8e-12 Eh.
- Subtracting only the s = 1 jump post hoc: 9.8e-12 Eh (full attribution).
- Analytic vs finite-difference force: HSE06 up to 1.4e-4 Eh/bohr vs PBE0 1.2e-7.
- Switching to eq. (9) would shift E_xc by +9.7e-5 Eh (H2O) and +3.3e-4 Eh (benzene), non-self-consistent.

## Notes from explorer 2 (verified in part)
- Second cut at s = 15 (clamp to s_max) and inconsistent inequalities w < 14 / w > 14 in wpbeh.mpl (measure-zero point).
- TDDFT analytic gradients with HSE06, LC-wPBE and pure GGA_X_WPBEH (omega = 0) do not converge to FD4 with grid
  refinement (CO, def2-SVP: ~2-4e-5 Eh/bohr vs 3e-9 for CAM-B3LYP); consistent with this s = 1 jump.
- Its 0.10-0.20 kcal/mol "vs CP2K" figure compares two different rescaling definitions and is NOT a defect measure.
