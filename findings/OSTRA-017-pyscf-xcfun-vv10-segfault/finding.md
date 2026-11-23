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
| status | Reported (verified, attributed, fix tested) |
| report | Reported: [pyscf#3679](https://github.com/pyscf/pyscf/issues/3679) (2026-11-23) |

`mf._numint.libxc = dft.xcfun` with `wb97m_v` or `b97m_v` on H2/STO-3G terminates with SIGSEGV (exit 139);
the same runs with libxc succeed.

## Re-check and cause (explorer 1, 2026-11-23)
See `evidence/2026-11-23-master/results.md`. The segfault happens in `xc_type`, before any evaluation:
`XCFUN_xc_type(5000+)` -> `xcfun_enumerate_parameters` returns NULL -> `xcfun_set(fun, NULL, ...)`.
The supported route (`mf.nlc = 'vv10'`) works with XCFun and matches libxc. Important for the fix: named VV10
functionals map only to the placeholder in the XCFun backend, so simply skipping ids >= 5000 would replace the
crash with a silent ~0.6 Eh error (no XC at all). Pending: blind verification by explorer 2; GitHub issue search.

## Prior art and history (2026-11-23)
GitHub search by the human reviewer: #919 "Bug in VV10 with XCFun backend" (2021, wb97m_v + XCFun, ZeroDivisionError),
closed by #928 (parser refactor, VV10_XC dict). The placeholder ids 5000+ were added later in #1681 (eb7d9043,
2023-04-11, "Refactor NLC parser and integral code"). Loading the pre-#1681 xcfun.py, `parse_xc('wb97m_v')` raised
KeyError 'Unknown XC functional WB97M_V': the change turned a clear error into a segfault. The #919 script segfaults on
2.15.1. Candidate fix (report/xcfun_vv10_guard.patch) applied to 2.15.1: named VV10 strings raise NotImplementedError;
mf.nlc route, b3lyp, pbe unchanged; test_xcfun.py 11 passed.
