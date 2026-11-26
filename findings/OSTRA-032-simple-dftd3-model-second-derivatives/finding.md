# OSTRA-032: simple-dftd3: d3_model weight_references / get_atomic_c6 ignore requested second CN derivatives unless first derivatives are also requested

| Field | Value |
|---|---|
| library | simple-dftd3 |
| location | src/dftd3/model.f90:245-254, 312-360 |
| versions | simple-dftd3 84cb7cd (main, v1.7.0) |
| detected | 2026-11-26 |
| detected_by | Explorer 2 (Gemini 4 Neon, Antigravity), guided search round 1 (pattern implicit-api-contract); verified by explorer 1 (Claude) |
| oracle | API contract: every combination of optional outputs vs the full call, sentinel-initialized outputs |
| libraries | simple-dftd3 |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | IV-III (public Fortran model API only; get_dispersion and C/Python/CLI unaffected) |
| status | Verified (code, explorer 1) and run (explorer 2); fix tested by explorer 2 |
| report | Pending |

gwd2cn is computed inside `if (present(gwdcn))`; second derivatives of C6 likewise. Explorer 2: gwd2cn alone returns
0 (reference 1.129); get_atomic_c6 leaves d2c6dcn2 / d2c6dcnij unset or crashes (SIGBUS) when second derivatives are
requested without first ones. Patch: `patches/simple-dftd3-model.patch`.

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.
