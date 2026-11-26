# OSTRA-030: tblite: xtb_singlepoint and several get_engrad routines skip the gradient unless sigma is also requested

| Field | Value |
|---|---|
| library | tblite |
| location | src/tblite/xtb/singlepoint.f90:151; repulsion/effective.f90:125; classical/halogen.f90:138; disp/d3, disp/d4; ceh.f90:1308 |
| versions | tblite c6f66b9 (main) |
| detected | 2026-11-26 |
| detected_by | Explorer 2 (Gemini 4 Neon, Antigravity), guided search round 1 (pattern implicit-api-contract); verified by explorer 1 (Claude) |
| oracle | API contract: every combination of optional outputs vs the full call, sentinel-initialized outputs |
| libraries | tblite |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | III-IV (Fortran API only; the C API, Python and CLI allocate both) |
| status | Verified (code, by explorer 1) and run (explorer 2); fix tested by explorer 2 |
| report | Pending |

`grad = present(gradient) .and. present(sigma)` (same pattern as OSTRA-025 in simple-dftd3). Explorer 2: gradient or
sigma requested alone is left unwritten (sentinel kept); with CDS solvation the gradient accumulates on uninitialized
memory. Verifier confirmed the three .and. conditions by reading main. Patch: `patches/tblite-singlepoint.patch`
(7/7 unit suites pass, bit-identical when both are requested).

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.
