# OSTRA-025: simple-dftd3 gCP: gradient silently zero unless the virial is also requested

| Field | Value |
|---|---|
| library | simple-dftd3 |
| location | src/dftd3/gcp.f90:104 (`present(gradient) .and. present(sigma)`); api.f90:1398-1418 |
| versions | v1.7.0 (84cb7cd) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | API contract; consistency with sibling routine (disp.f90:102 uses .or.) |
| class | API contract, new |
| impact_level | III (C API / Fortran callers requesting gradients without sigma) |
| status | Verified (code) |
| report | none yet |

With sigma absent (C API with NULL sigma: the unallocated Fortran array is passed as absent), grad is false
and only energies are computed; the gradient stays at zero without a warning. The dispersion routine uses
`present(gradient) .or. present(sigma)` and handles each output separately.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/07-simple-dftd3-gcp-grad-contract/finding.md` (copied below as `explorer_finding.md`).
