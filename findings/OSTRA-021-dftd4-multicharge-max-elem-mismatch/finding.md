# OSTRA-021: dftd4 accepts Z = 112-118 but multicharge EEQ parameters stop at Z = 103: unphysical charges without error

| Field | Value |
|---|---|
| library | dftd4 x multicharge |
| location | dftd4 src/dftd4/model/d4.f90:106-111; multicharge param/eeq2019.f90 (max_elem = 103) |
| versions | dftd4 4a5bbf0 and released 4.3.1; multicharge 4bb6c96 |
| detected | 2026-11-22 |
| detected_by | Explorer 2 (Gemini, Antigravity), batch 3; verified by explorer 1 (Claude) |
| oracle | Chemical sanity of charges; API contract |
| libraries | dftd4, multicharge |
| pattern | fixed-size-table |
| oracle_class | api-contract |
| class | Fixed-size table across libraries, new |
| impact_level | III (superheavy elements only) |
| status | Verified (released dftd4 4.3.1) |
| report | Pending |

d4.f90 rejects only 104-111; EEQ getters return -1.0 beyond 103 and new_eeq_model does not validate them.
Verifier with dftd4 4.3.1 (Python): Hg-O gives q(Hg) = +0.48; Cn-O and Og-O give q(Cn) = -0.87, q(Og) = -0.87 and
oxygen +0.87, no error raised.

Explorer's original record: `agents/explorer2/exploracion3-snapshot-2026-11-22/03-dftd4-multicharge-maxelem/finding.md` (copied below as `explorer_finding.md`).
