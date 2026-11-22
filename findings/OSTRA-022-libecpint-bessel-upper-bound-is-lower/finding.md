# OSTRA-022: libecpint BesselFunction::upper_bound returns a lower bound on the rising side of K_L(z)

| Field | Value |
|---|---|
| library | libecpint |
| location | src/lib/bessel.cpp:145-152; include/libecpint/bessel.hpp:107 |
| versions | v1.1 (d7966f1, 2026-11-18) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | Comparison with exact values (explorer); code reading (verifier) |
| class | Documented contract violated, new |
| impact_level | IV inside libecpint (no caller); III for external users of the public header |
| status | Verified (code); impact downgraded |
| report | none yet |

upper_bound takes the tabulated point to the left of z, which bounds a decreasing function from above but
is a lower bound where K_L(z) rises (L >= 1 below its maximum). Explorer reports underestimates up to 1958x
(L = 24, z = 34.975). No function inside libecpint calls upper_bound.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/04-libecpint-bessel-upper-bound/finding.md` (copied below as `explorer_finding.md`).
