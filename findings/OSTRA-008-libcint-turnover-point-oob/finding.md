# OSTRA-008: Out-of-bounds read of TURNOVER_POINT[40] for Boys order m >= 40 corrupts SR integrals at high l

| Field | Value |
|---|---|
| library | libcint |
| location | src/fmt.c:42 (TURNOVER_POINT), callers in rys_roots.c / rys_wheeler.c |
| versions | v6.2.1 (ab88a54) and PySCF 2.15.1 bundled; table unchanged since 8fdc4d8 (2021-07) |
| detected | 2026-11-21 |
| detected_by | Explorer 1 (threshold mining) |
| oracle | Boys values vs 40-digit quadrature; erfc + erf = full Coulomb |
| class | Memory safety, new |
| impact_level | III (SR integrals with l >~ 10; errors up to 3e8, build-dependent) |
| status | Verified and attributed; issue drafted |
| report | Pending (issue drafted in report/) |

MXRYSROOTS = 32 lets the Rys code request m up to 2*nroots (Schmidt fallback) or 2*nroots-1 (SR Jacobi moments),
but the table has 40 entries. Reading past it is undefined behaviour; the failing orders differ between binaries
(PySCF bundled: from m = 41; v6.2.1 release: m = 43, 45, 47). Extending the table with the formula in its own
comment (reproduces the 38 non-zero entries to 8e-12) removes every error. Patch: `report/turnover_point_extend.patch`.
