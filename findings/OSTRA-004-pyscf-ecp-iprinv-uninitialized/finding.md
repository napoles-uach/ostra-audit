# OSTRA-004: ECP-centre derivative kernels return uninitialized memory when the rinv atom has no ECP

| Field | Value |
|---|---|
| library | PySCF (libcgto) |
| location | pyscf/lib/gto/nr_ecp_deriv.c |
| versions | PySCF 2.15.1 |
| detected | 2026-11-20 |
| detected_by | Explorer 1 (harness bug exposed it) |
| oracle | API contract (outputs must be written) |
| libraries | PySCF |
| pattern | implicit-api-contract |
| oracle_class | api-contract |
| class | API contract, new |
| impact_level | III (zero inside PySCF; garbage up to 1e256 for direct API users) |
| status | Reported |
| report | Reported: [pyscf#3669](https://github.com/pyscf/pyscf/issues/3669) |

`ECPscalar_iprinv` (and `_ipiprinv`, `_iprinvip`) return 0 without writing the output when `AS_RINV_ORIG_ATOM`
points to an atom without ECP; `getints2c` allocates with `numpy.ndarray` (uninitialized).
Full report: `report/issue.md`.
