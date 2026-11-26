# OSTRA-029: dftd4: get_dispersion with sigma but no gradient crashes (C API: dftd4_get_dispersion with gradient = NULL)

| Field | Value |
|---|---|
| library | dftd4 |
| location | src/dftd4/disp.f90:89,126-127,207; api.f90:725-744 (comment: 'Would need to refactor ... to make sigma truly optional') |
| versions | dftd4 4a5bbf0 (main) and released 4.3.1 |
| detected | 2026-11-26 |
| detected_by | Explorer 2 (Gemini 4 Neon, Antigravity), guided search round 1 (pattern implicit-api-contract); verified by explorer 1 (Claude) |
| oracle | API contract: every combination of optional outputs vs the full call, sentinel-initialized outputs |
| libraries | dftd4 |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | III (C API and Fortran callers; Python and CLI always pass both; loud: segfault) |
| status | Verified (run on released 4.3.1 by explorer 1), attributed; fix tested by explorer 2 |
| report | Pending |

Explorer 2: SIGBUS in Fortran with only gradient or only sigma; ASan BUS WRITE in the C API with gradient = NULL,
sigma != NULL. Verifier (explorer 1, independent script with the released Python package calling the C function
directly): gradient only -> correct; sigma only -> segmentation fault (exit 139); both -> correct. The C wrapper
already allocates sigma when only the gradient is requested, with a comment acknowledging that sigma is not truly
optional; the opposite case is not handled. Patch: explorer 2, `patches/dftd4-disp-capi.patch` (meson test 54 OK).

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.
