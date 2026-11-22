# OSTRA-023: simple-dftd3 gCP: B_k(x) series vs closed form chosen by abs(za - zb) instead of abs(x)

| Field | Value |
|---|---|
| library | simple-dftd3 |
| location | src/dftd3/gcp.f90:806, 1106-1132 ('12 terms seem accurate enough'), 1303, 1378-1380 |
| versions | v1.7.0 (84cb7cd) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | mpmath reference; FD gradient vs Hessian (explorer) |
| class | Numerical precision, new |
| impact_level | II (gCP derivatives for unlike-element pairs; energy impact not re-measured) |
| status | Partly verified |
| report | none yet |

Verifier (mpmath): the 12-term series has relative errors 3.5e-4 to 1.8e-3 at x = 5; the closed-form
upward recursion has 4e-3 to 0.92 at x = 0.05-0.1 for k = 7-8. The switch uses |za - zb| < 0.1, while x = (zb - za)R/2
depends on R, so each formula can be used in the regime where it fails. Explorer: Cl-Ge (r2SCAN-3c) FD gradient vs
Hessian diverges to 4.7e-3 at h = 1e-5, 1e-11 with the patch.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/05-simple-dftd3-gcp-bint-baux/finding.md` (copied below as `explorer_finding.md`).
