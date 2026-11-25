# OSTRA-014: R4SCANX uses d_p4 = d_a4 = 0.232 instead of 0.802 / 0.178

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/SCAN_like_eps.hpp:156-157 |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Differential vs libxc with exact parameter attribution |
| libraries | XCFun |
| pattern | stale-threshold |
| oracle_class | independent-implementation |
| class | New bug (parameters) or outdated preprint values |
| impact_level | I (silent; atomization energies MAD 0.63, max 1.99 kcal/mol, impact study 2026-11-24) |
| status | Confirmed (level 3): published values 0.802/0.178; minimal fix matches libxc to 3e-14 |
| report | Pending |

XCFun vs libxc defaults: 1.0e-2 max relative difference; libxc with dp4 = da4 = 0.232: 1.5e-14 (exact attribution).
XCFun added r4SCAN in 2020-09 (#140), before the r4SCAN paper; 0.232 may be a preliminary value.

## Confirmation phase (2026-11-21) — level 3: defect scientifically confirmed
- **Publication check.** r4SCAN preprint (arXiv:2110.00647v1, Furness et al., = JCP 156, 034109): d_p2 = 0.361,
  d_p4 = 0.802, d_alpha4 = 0.178 (eqs. C40, C51, C52; Fig. 5 shows the optimisation). 0.232 does not appear in it.
- **Minimal fix:** DX_DAMP4_P = 0.802, DX_DAMP4_A = 0.178.
- **Result** vs libxc 7.1.2: eps 2.4e-2 -> 3.2e-14, derivatives 1.3e2 -> 1.4e-10. Controls unchanged.
- **Status:** confirmed defect (parameters). Origin of 0.232 unknown (r4SCAN entered XCFun in 2020-09, before the preprint).

## Impact on observables (2026-11-24)
`studies/impact-xcfun-2026-11-24/`: self-consistent r4SCAN, released vs fixed XCFun, 26 molecules, def2-TZVP.
Atomization energies: MAD 0.63, max 1.99 kcal/mol (C2H6); reaction energies: MAD 0.58, max 2.03 kcal/mol
(2 H2O2 -> 2 H2O + O2). Fixed build matches libxc. Silent. Impact level I confirmed.
