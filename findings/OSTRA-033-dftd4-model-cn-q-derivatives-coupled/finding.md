# OSTRA-033: dftd4: d4_model / d4s_model couple CN and charge derivatives with .and.; requesting only one leaves the other intent(out) unset

| Field | Value |
|---|---|
| library | dftd4 |
| location | src/dftd4/model/d4.f90:270 (present(gwdcn) .and. present(gwdq)) and 5 more methods in d4.f90 / d4s.f90 |
| versions | dftd4 4a5bbf0 (main) and released 4.3.1 |
| detected | 2026-11-26 |
| detected_by | Explorer 2 (Gemini 4 Neon, Antigravity), guided search round 1 (pattern implicit-api-contract); verified by explorer 1 (Claude) |
| oracle | API contract: every combination of optional outputs vs the full call, sentinel-initialized outputs |
| libraries | dftd4 |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | IV-III (public Fortran model API; not reachable from C/Python/CLI) |
| status | Verified (code, explorer 1) and run (explorer 2); fix tested by explorer 2 |
| report | Pending |

Explorer 2: in weight_references, get_atomic_c6 and get_polarizabilities of both models, asking only for derivatives
with respect to CN or only with respect to q leaves the requested intent(out) arrays unwritten (12 partial cases).
Patch: `patches/dftd4-model.patch` (meson test 54 OK, bit-identical).

Explorer's record: `agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/reporte_contratos.md`; inventory of the 56 routines checked (positives and negatives):
`agents/explorer2/exploracion5-contratos-snapshot-2026-11-26/rutinas_revisadas.csv`.
