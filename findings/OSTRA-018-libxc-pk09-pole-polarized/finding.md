# OSTRA-018: lda_c_pk09: NaN and huge values in a band around r_s = 22.739 for spin-polarized densities

| Field | Value |
|---|---|
| library | libxc |
| location | maple/lda_exc/lda_c_pk09.mpl (alpha_z with exponent s_s(r_s); DETOL guard lines 22-25) |
| versions | 7.1.2 and devel 3c7794c50 |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | Edge scan; zeta = 0 vs zeta != 0 comparison |
| class | Numerical bug or model singularity (pending the paper) |
| impact_level | II-III (open-shell densities near rho ~ 2e-5) |
| status | Verified (anomaly reproduced); cause: model pole vs implementation pending |
| report | Pending |

Explorer 2 (batch 1) reported a NaN band; it was not reproduced at zeta = 0. In batch 3 it gave the exact
condition: spin-polarized with 0 < |zeta| < 1. Verified in 7.1.2 (zeta = 0.4): NaN for |r_s - 22.7394121| < ~5e-5 and
|eps_c| more than twice the zeta = 0 value over +-5e-4 in r_s (e.g. eps_c = -1287 Ha at -1e-4, normal -0.0104).
At zeta = 0 nothing happens because 1^s_s = 1. The DETOL = 1e-12 guard only acts within ~4e-16 of the pole.
Pending: whether the pole of s_s(r_s) is in the published PK09 form (Proynov-Kong 2009) or an implementation detail;
estimate how often molecular grids hit the band.
