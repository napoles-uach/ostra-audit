# OSTRA-015: SCANX/SCANC return +-inf at alpha = 1 (uniform gas) and SCANX v_sigma is halved at s = 0

| Field | Value |
|---|---|
| library | XCFun |
| location | src/functionals/SCAN_like_eps.hpp (interpolation 1/(1-alpha)) |
| versions | 2.1.1 (PySCF 2.15.1 bundled; also macOS PySCF 2.14.0); same code on master a89b783 |
| detected | 2026-11-21 |
| detected_by | Explorer 2 (XCFun vs libxc differential audit); verified by explorer 1 |
| oracle | Edge cases; FD of v_sigma |
| class | Unguarded edge case |
| impact_level | III-IV (exact points rarely on molecular grids) |
| status | Verified |
| report | none yet |

At alpha = 1 exactly: SCANX exc = +inf, SCANC exc = -inf, v_sigma NaN (libxc finite). At s = 0, alpha = 1.2:
SCANX v_sigma = -0.0508 vs libxc -0.0993; finite differences give -0.0993, so XCFun is wrong there.
