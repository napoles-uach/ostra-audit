# OSTRA-031: mctc-lib: get_coordination_number leaves dcndr unset unless dcndL is also requested; with cut > 0 it rescales the unset memory

| Field | Value |
|---|---|
| library | mctc-lib |
| location | src/mctc/ncoord/type.f90:216-230 (present(dcndr) .and. present(dcndL)); cut_coordination_number |
| versions | mctc-lib 71a21af (main) |
| detected | 2026-11-26 |
| detected_by | Explorer 2 (Gemini 4 Neon, Antigravity), guided search round 1 (pattern implicit-api-contract); verified by explorer 1 (Claude) |
| oracle | API contract: every combination of optional outputs vs the full call, sentinel-initialized outputs |
| libraries | mctc-lib |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | III (Fortran API of mctc-lib and its users in dftd4/simple-dftd3/tblite/multicharge; silent) |
| status | Verified (run by explorer 1), attributed; fix tested by explorer 2 |
| report | Pending |

Verifier (explorer 1, own Fortran program against mctc-lib 71a21af): H2O, erf counting function. dcndr requested
alone: left at the sentinel 999 (error 1e3); with cut = 1.5 the unset values are rescaled to 444-654 by
cut_coordination_number. Explorer 2: same for dcndL and dcndrlist; patch `patches/mctc-lib-ncoord.patch`
(38/38 tests, FD 7.5e-12).

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.
