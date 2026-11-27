# OSTRA-018: lda_c_pk09: NaN and huge values in a band around r_s = 22.739 for spin-polarized densities

| Field | Value |
|---|---|
| library | libxc |
| location | maple/lda_exc/lda_c_pk09.mpl (alpha_z with exponent s_s(r_s); DETOL guard lines 22-25) |
| versions | 7.1.2 and devel 3c7794c50 |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | Edge scan; zeta = 0 vs zeta != 0 comparison |
| libraries | libxc |
| pattern | none |
| oracle_class | branch-boundary, exact-identity |
| class | Model singularity in the published fit (design); residual NaN band in the implementation |
| impact_level | IV on molecular grids (band rarely hit, no measurable energy effect); III for direct UEG use near r_s = 22.74 |
| status | Attributed (paper checked) |
| report | Not to be reported (pole of the published fit; libxc devel already screens it, 2026-11-19); optional comment on the residual NaN band |

Explorer 2 (batch 1) reported a NaN band; it was not reproduced at zeta = 0. In batch 3 it gave the exact
condition: spin-polarized with 0 < |zeta| < 1. Verified in 7.1.2 (zeta = 0.4): NaN for |r_s - 22.7394121| < ~5e-5 and
|eps_c| more than twice the zeta = 0 value over +-5e-4 in r_s (e.g. eps_c = -1287 Ha at -1e-4, normal -0.0104).
At zeta = 0 nothing happens because 1^s_s = 1. The DETOL = 1e-12 guard only acts within ~4e-16 of the pole.
Pending: whether the pole of s_s(r_s) is in the published PK09 form (Proynov-Kong 2009) or an implementation detail;
estimate how often molecular grids hit the band.

## Attribution against the paper (2026-11-27)
Proynov and Kong, Phys. Rev. A 79, 014103 (2009), Eqs. (22)-(25): alpha_xi = 2/[(1+xi)^s + (1-xi)^s] with
s(r_s, xi) = f_r(r_s) f_s(xi) and f_r = P1/P2 a fitted rational function. With the published coefficients P2 has a
real root at r_s = 22.73941213 (np.roots of b_0..b_4), so s -> +-infinity there. The pole is in the published model,
inside its stated range of validity (0.1 <= r_s <= 30; the paper tabulates r_s = 20 and 30 only). At xi = 0, 1^s = 1
and nothing happens; for 0 < |xi| < 1 alpha_xi -> 0 as |s| -> infinity, so the energy has a finite limit but a sharp
unphysical feature (eps_c -0.0316 at r_s = 22.739 vs -0.0099 nearby, xi = 0.4; potential up to 1.5e3) and floating-
point NaN within ~5e-5 of the pole. libxc reproduces the paper faithfully (the extra factor 1.28 comes from the authors'
code, erratum); libxc devel added a DETOL screen of the exact zero of P2 on 2026-11-19, so the maintainers are aware
of the pole; the NaN band of width ~1e-4 remains.

## Frequency on molecular grids (2026-11-27)
UKS densities of OH, O2, CH3 and NO (def2-TZVP, grid levels 3, 5, 7): at most 4 points within 5e-4 of the pole (CH3),
no non-finite values, PK09 correlation energy unchanged to 1e-6 Eh between grid levels. Practical impact negligible.

