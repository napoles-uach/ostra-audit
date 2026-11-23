# OSTRA-008: Out-of-bounds read of TURNOVER_POINT[40] for Boys order m >= 40 corrupts SR integrals at high l

| Field | Value |
|---|---|
| library | libcint |
| location | src/fmt.c:42 (TURNOVER_POINT), callers in rys_roots.c / rys_wheeler.c |
| versions | v6.2.1 (ab88a54) and PySCF 2.15.1 bundled; table unchanged since 8fdc4d8 (2021-07) |
| detected | 2026-11-21 |
| detected_by | Explorer 1 (Claude, threshold mining); independently verified by explorer 2 (Gemini, Antigravity), 2026-11-23 |
| oracle | Boys values vs 40-digit quadrature; erfc + erf = full Coulomb |
| class | Memory safety, new |
| impact_level | III (SR integrals with l >~ 10; errors up to 3e8, build-dependent) |
| status | Reported (verified independently, attributed, minimal fix verified) |
| report | Reported: [libcint#143](https://github.com/sunqm/libcint/issues/143) (2026-11-23) |

MXRYSROOTS = 32 lets the Rys code request m up to 2*nroots (Schmidt fallback) or 2*nroots-1 (SR Jacobi moments),
but the table has 40 entries. Reading past it is undefined behaviour; the failing orders differ between binaries
(PySCF bundled: from m = 41; v6.2.1 release: m = 43, 45, 47). Extending the table with the formula in its own
comment (reproduces the 38 non-zero entries to 8e-12) removes every error. Patch: `report/turnover_point_extend.patch`.

## Independent verification (explorer 2, Gemini, 2026-11-23)
Blind protocol: the verifier received only the location and the claim; the patch was given after phase 1.
Material: `verification/gemini-2026-11-23/` (verification.md, scripts, logs, patched fmt.c). Platform: macOS arm64.

- **Confirmed:** 40-entry table; OOB for m >= 40; SR path reads indices 41/43/45/47 for nroots = 21..24 (nroots >= 25
  exits in `CINTsr_rys_roots`); Schmidt fallback with m = 2*nroots (observed with PySCF's libcint 6.1.3, omega > 0).
  ASan global-buffer-overflow at fmt.c:214/256/363/404 in four binaries (clang -O0/-O2, gcc-16 + libquadmath, PySCF
  libcint 6.1.3); error size depends on what follows the table in memory (clang -O0 happened to be correct in one
  case). Silent example: 3.77 % error, SR int2e_sph (11,11|10,10), gcc-16, stderr empty.
- **Patch:** formula reproduces entries m = 2..39 to 8.2e-12 abs; ASan clean; SR integrals 5e-12 to 7e-11 vs full - LR;
  bit-identical for m <= 39 and int2e_sph l <= 5.
- **Corrections to our draft:** 71 is a loose upper bound, not the exact maximum (reachable: 47 in SR, 62-63 via
  MXRYSROOTS and the fmt[64] buffer); direct API calls with m >= 72 still overflow. Issue text updated accordingly,
  plus a macOS arm64 confirmation paragraph and the signature (Claude and Gemini).
- **Side observations, not part of this finding (moved to candidates):** (1) double-precision Boys values lose digits
  just above t0, already for existing entries (re-checked by explorer 1 with the patched Linux build:
  gamma_inc_like 1e-8 at m = 39, 8e-8 at m = 47, 1e-4 at m = 63, 1.5e-3 at m = 71; long double 3e-12 / 1e-10 / 7e-8
  at m = 39 / 47 / 63). Relevant where long double is 64-bit and quadmath is absent (macOS arm64), since
  CINTqrys_* then map to CINTlrys_*. (2) A singularity in lrys_wheeler at R = 1.0 on macOS arm64 without quadmath,
  present also in the m = 39 control.
