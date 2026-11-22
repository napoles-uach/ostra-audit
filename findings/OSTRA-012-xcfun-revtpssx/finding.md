# OSTRA-012: REVTPSSX registered without XC_KINETIC and with two formula errors

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/revtpssx.cpp, revtpssx_eps.hpp |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Code comparison with TPSSX in the same library; differential vs libxc |
| class | New bug (registration + formula) |
| impact_level | III (fails loudly: NaN as registered) |
| status | Verified by code reading; NaN reproduced |
| report | none yet (XCFun last commit 2022-11) |

(a) flags lack XC_KINETIC, so PySCF types it as GGA and it returns NaN on all test points; (b) denominator
(1 + sqrt(e) p) instead of (1 + sqrt(e) p)^2; (c) inside the square root `p` is used where TPSSX uses `p0 = p/sigma`,
so the term carries an extra factor sigma. (b) and (c) only matter once (a) is fixed.
