# OSTRA-024: libecpint BesselFunction::calculate at z <= 0 writes only K_0; scalar overload returns 1 for all L

| Field | Value |
|---|---|
| library | libecpint |
| location | src/lib/bessel.cpp:162-163, 223-224 |
| versions | v1.1 (d7966f1) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | API contract (comment states K_l = 0 for l > 0) |
| libraries | libecpint |
| pattern | unverified-comment, implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | IV inside libecpint (callers guard x > 0 or rely on zero-initialized buffers); III for external users |
| status | Verified (code); impact downgraded |
| report | Pending |

The vector overload leaves values[1..maxL] untouched at z <= 0, contradicting its own comment; the scalar
overload returns 1.0 instead of 0.0 for L >= 1. radial.cpp guards x > 0 and sets the orders by hand (with a comment
explaining why); shared_grid.cpp is safe because its buffer starts at zero and x is fixed per call.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/06-libecpint-bessel-contract/finding.md` (copied below as `explorer_finding.md`).
