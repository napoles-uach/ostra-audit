# OSTRA-017: Segmentation fault for VV10/rVV10 functionals with the XCFun backend

| Field | Value |
|---|---|
| library | PySCF |
| location | pyscf/dft/xcfun.py:368 (placeholder ids 5000+), reaching XCFUN_xc_type via xc_type (xcfun.py:338) and XCFUN_eval_xc |
| versions | PySCF 2.14.0 (macOS) and 2.15.1 (Linux); code unchanged on master 6b288fd (2026-11-23) |
| detected | 2026-11-21 |
| detected_by | Explorer 2; verified by explorer 1 |
| oracle | Crash / exit code |
| libraries | PySCF |
| pattern | inconsistent-configuration |
| oracle_class | api-contract |
| class | New bug (crash) |
| impact_level | III (loud: process killed, exit 139) |
| status | Verified |
| report | Pending |

`mf._numint.libxc = dft.xcfun` with `wb97m_v` or `b97m_v` on H2/STO-3G terminates with SIGSEGV (exit 139);
the same runs with libxc succeed.

## Re-check and cause (explorer 1, 2026-11-23)
See `evidence/2026-11-23-master/results.md`. The segfault happens in `xc_type`, before any evaluation:
`XCFUN_xc_type(5000+)` -> `xcfun_enumerate_parameters` returns NULL -> `xcfun_set(fun, NULL, ...)`.
The supported route (`mf.nlc = 'vv10'`) works with XCFun and matches libxc. Important for the fix: named VV10
functionals map only to the placeholder in the XCFun backend, so simply skipping ids >= 5000 would replace the
crash with a silent ~0.6 Eh error (no XC at all). Pending: blind verification by explorer 2; GitHub issue search.
