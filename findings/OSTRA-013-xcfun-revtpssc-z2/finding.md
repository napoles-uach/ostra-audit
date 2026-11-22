# OSTRA-013: REVTPSSC uses z^2 instead of z^3 in eps_c = eps_PKZB (1 + d eps_PKZB z^3)

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/revtpssc_eps.hpp:107-109 |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Code comparison with TPSSC; differential vs libxc |
| class | New bug (formula), silent |
| impact_level | I-II (explorer: SCF shifts 1.5-2.1 kcal/mol; not re-measured) |
| status | Verified |
| report | none yet |

`revtpssc_eps.hpp` uses `pow2(...)` where `tpssc_eps.hpp` (same formula, same d = 2.8) uses `pow3(...)`.
Max relative difference vs libxc 2 % on 200 random points (control TPSSC vs libxc: 3e-11). Silent: passes the
uniform gas (z = 0) and one-electron densities (z = 1).
