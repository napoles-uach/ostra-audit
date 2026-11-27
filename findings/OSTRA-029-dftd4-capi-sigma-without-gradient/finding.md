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
| report | Pending (issue drafted in OSTRA-029/report/ together with 033; target repo dftd4/dftd4) |

Explorer 2: SIGBUS in Fortran with only gradient or only sigma; ASan BUS WRITE in the C API with gradient = NULL,
sigma != NULL. Verifier (explorer 1, independent script with the released Python package calling the C function
directly): gradient only -> correct; sigma only -> segmentation fault (exit 139); both -> correct. The C wrapper
already allocates sigma when only the gradient is requested, with a comment acknowledging that sigma is not truly
optional; the opposite case is not handled. Patch: explorer 2, `patches/dftd4-disp-capi.patch` (meson test 54 OK).

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.

## Re-verification before reporting (2026-11-27, explorer 1)
Canonical repository is now github.com/dftd4/dftd4 (grimme-lab/dftd4 is a stale copy at 4.3.0). On main 4a5bbf0,
own C and Fortran programs (NH3, PBE0-D4): C API sigma-only and Fortran gradient-only both segfault; d4_model
weight_references with gwdcn only / gwdq only leaves the array at the sentinel (029 and 033). Explorer 2's two patches
combined (api.f90, disp.f90, model/d4.f90, model/d4s.f90): all cases identical to the full call (0.0; 2.7e-20 with
OpenMP); meson test 76 OK + 1 XFAIL (debug, no OpenMP) and 75 OK + 1 XFAIL (release, OpenMP, 4 threads), same as
unpatched. Note: simple-dftd3 likewise moved to github.com/dftd3/simple-dftd3.
