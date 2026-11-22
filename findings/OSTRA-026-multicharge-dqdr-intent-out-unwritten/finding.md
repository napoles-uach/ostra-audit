# OSTRA-026: multicharge get_charges: intent(out) dqdr left undefined when dqdL is not requested

| Field | Value |
|---|---|
| library | multicharge |
| location | src/multicharge/charge.f90:67; model/type.F90:219-221 |
| versions | 4bb6c96 (2026-11-17) |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | API contract (intent(out) must be defined) |
| class | API contract, new |
| impact_level | III (callers requesting only Cartesian charge derivatives) |
| status | Verified (code) |
| report | Pending |

grad = (present(dqdr) .and. present(dqdL)) .or. present(cache); with only dqdr present no derivative is
computed and the intent(out) array is never written.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/08-multicharge-dqdr-contract/finding.md` (copied below as `explorer_finding.md`).
