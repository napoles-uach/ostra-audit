# OSTRA-012: REVTPSSX registered without XC_KINETIC and with two formula errors

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/revtpssx.cpp, revtpssx_eps.hpp |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Code comparison with TPSSX in the same library; differential vs libxc |
| libraries | XCFun |
| pattern | inconsistent-configuration, transcription-slip |
| oracle_class | independent-implementation |
| class | New bug (registration + formula) |
| impact_level | III (fails loudly: NaN as registered) |
| status | Confirmed (level 3): minimal fix matches the publication and libxc to 1e-15 |
| report | Pending |

(a) flags lack XC_KINETIC, so PySCF types it as GGA and it returns NaN on all test points; (b) denominator
(1 + sqrt(e) p) instead of (1 + sqrt(e) p)^2; (c) inside the square root `p` is used where TPSSX uses `p0 = p/sigma`,
so the term carries an extra factor sigma. (b) and (c) only matter once (a) is fixed.

## Confirmation phase (2026-11-21) — level 3: defect scientifically confirmed
- **Publication check.** revTPSS preprint (arXiv:0906.5029, = PRL 103, 026403): exchange keeps the TPSS form with
  c z^2 p/(1+z^2)^2 replaced by c z^3 p/(1+z^2)^2 and new constants (mu = 0.14, c = 2.35204, e = 2.1677); the TPSS
  form has the denominator (1 + sqrt(e) p)^2. The 2011 erratum (PRL 106, 179902) only corrects cohesive-energy tables.
- **Minimal fix** (`xcfun_minimal_fixes.patch`): add XC_KINETIC; use p0 inside the square root; square the denominator.
- **Result** (`validation_results.txt`, 400 random points vs libxc 7.1.2): type GGA/NaN -> MGGA; eps 1.0e-8, derivatives
  3.6e-6. The remainder is a constant-precision convention: libxc uses c = 2.35203946, e = 2.16769874. With those
  (`xcfun_fixes_plus_libxc_constant_conventions.patch`): eps 9.6e-16, derivatives 1.0e-14. Controls TPSSX, TPSSC,
  R2SCANX, SCANX unchanged.
- **Status:** confirmed defect (registration + two formula errors). Convention difference in constants is not a defect.

## Impact on observables (2026-11-24)
`studies/impact-xcfun-2026-11-24/`: with the released XCFun, revTPSS SCF in PySCF does not finish for 7 of 26
molecules (singular matrix or > 15 min); loud failure, no silent wrong results observed.
