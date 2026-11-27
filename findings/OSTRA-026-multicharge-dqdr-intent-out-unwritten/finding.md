# OSTRA-026: multicharge get_charges: intent(out) dqdr left undefined when dqdL is not requested

| Field | Value |
|---|---|
| library | multicharge |
| location | src/multicharge/charge.f90:67; model/type.F90:219-221 |
| versions | 4bb6c96 (2026-11-17); still present in aa2daad (2026-11-27), patch re-tested there |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | API contract (intent(out) must be defined) |
| libraries | multicharge |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | III (callers requesting only Cartesian charge derivatives) |
| status | Verified (run), attributed, fix tested; issue drafted |
| report | Pending (issue drafted in OSTRA-026/report/, together with 021) |

grad = (present(dqdr) .and. present(dqdL)) .or. present(cache); with only dqdr present no derivative is
computed and the intent(out) array is never written.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/08-multicharge-dqdr-contract/finding.md` (copied below as `explorer_finding.md`).

## Dynamic confirmation and tested fix (explorer 1, 2026-11-25)
multicharge main 4bb6c96 built with meson; Fortran test (verification/test_multicharge.f90): dqdr requested alone is
not written (sentinel 999 kept). Patch in report/: equal to dqdr of the full call; meson test 52 OK.
