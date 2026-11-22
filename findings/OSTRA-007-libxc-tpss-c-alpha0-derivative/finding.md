# OSTRA-007: Derivative of TPSS/TM-family correlation at tau = tau_W (z = 1) depends on rounding

| Field | Value |
|---|---|
| library | libxc |
| location | TPSS-family correlation (z clamp) |
| versions | 7.1.2 and devel |
| detected | 2026-11-21 |
| detected_by | Explorer 1 (spin-polarization oracle) |
| oracle | Spin-polarized vs unpolarized identity; RKS vs UKS |
| class | No impact |
| impact_level | IV |
| status | Documented; not reported |
| report | Not to be reported (no impact) |

At alpha = 0 exactly, v_tau is 0 on half of the He grid points (clamped branch) and finite on the other half.
RKS and UKS energies of He and H2 nevertheless agree to <1e-15 Eh for TPSS, SCAN and M06-L: in one-orbital
systems every orbital variation keeps tau = tau_W, so delta z = 0 to first order and the ambiguous derivative
always multiplies zero. Could matter only for second derivatives (f_xc).
