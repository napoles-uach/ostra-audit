# OSTRA-017: Segmentation fault for VV10/rVV10 functionals with the XCFun backend

| Field | Value |
|---|---|
| library | PySCF |
| location | pyscf/dft/xcfun.py:368 (placeholder ids 5000+) passed unfiltered to XCFUN_eval_xc |
| versions | PySCF 2.14.0 (macOS) and 2.15.1 (Linux) |
| detected | 2026-11-21 |
| detected_by | Explorer 2; verified by explorer 1 |
| oracle | Crash / exit code |
| class | New bug (crash) |
| impact_level | III (loud: process killed, exit 139) |
| status | Verified |
| report | none yet (check pyscf issues first) |

`mf._numint.libxc = dft.xcfun` with `wb97m_v` or `b97m_v` on H2/STO-3G terminates with SIGSEGV (exit 139);
the same runs with libxc succeed.
