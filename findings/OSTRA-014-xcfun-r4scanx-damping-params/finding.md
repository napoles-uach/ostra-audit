# OSTRA-014: R4SCANX uses d_p4 = d_a4 = 0.232 instead of 0.802 / 0.178

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/SCAN_like_eps.hpp:156-157 |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Differential vs libxc with exact parameter attribution |
| class | New bug (parameters) or outdated preprint values |
| impact_level | I-II (explorer: 0.5-0.9 kcal/mol; not re-measured) |
| status | Verified vs libxc; published values to be confirmed in Furness et al., JCP 156, 034109 (2022) |
| report | none yet |

XCFun vs libxc defaults: 1.0e-2 max relative difference; libxc with dp4 = da4 = 0.232: 1.5e-14 (exact attribution).
XCFun added r4SCAN in 2020-09 (#140), before the r4SCAN paper; 0.232 may be a preliminary value.
