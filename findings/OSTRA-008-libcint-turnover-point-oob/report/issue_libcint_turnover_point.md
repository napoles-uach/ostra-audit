TITLE: Out-of-bounds read of TURNOVER_POINT in fmt.c for Boys order m >= 40 corrupts SR integrals with high angular momentum

## Summary

`TURNOVER_POINT[]` in `src/fmt.c` has 40 entries (m = 0..39), but the Boys-like functions that index it are called with larger m. The Rys-root code asks for orders up to m = 2*nroots - 1 (Jacobi moments for SR integrals) or m = 2*nroots (Schmidt fallback), and `MXRYSROOTS = 32` allows nroots up to about 31. For m >= 40 the code reads past the end of the array. That is undefined behaviour: the branch between the power series and upward recursion then depends on whatever happens to sit in memory after the table. When the value read is small, upward recursion is used at small t, where it is unstable, and F_m comes out wrong by many orders of magnitude.

Checked on master (v6.2.1, ab88a54). The table has had 40 entries since 8fdc4d8 ("Improve fmt function accuracy", July 2021).

## Where

- Table: https://github.com/sunqm/libcint/blob/ab88a54/src/fmt.c#L42 (40 entries)
- Indexed in `gamma_inc_like`, `lgamma_inc_like`, `qgamma_inc_like`, `fmt_erfc_like`, `fmt_lerfc_like`, `fmt_qerfc_like`
- Reached with m = 2*nroots - 1 from `lnaive_jacobi_moments` / `qnaive_jacobi_moments` (rys_wheeler.c, called by `CINTlrys_jacobi` / `CINTqrys_jacobi` with lower != 0, i.e. SR integrals). `CINTsr_rys_roots` uses these for nroots up to 24, so nroots = 21..24 read entries 41..47.
- Reached with m = 2*nroots from `CINTlrys_schmidt` / `CINTqrys_schmidt` (rys_roots.c), which `segment_solve` uses as a fallback.

## Evidence

1) Boys values. `libcint_turnover_repro.py` (below) calls the exported `fmt_lerfc_like(f, t, lower=0.2, m)` and compares F_m with a 40-digit quadrature. m <= 41 is correct in both builds below and is omitted.

| m | t | v6.2.1 as released | v6.2.1, table extended |
|---|---|---|---|
| 43 | 0.5 | 1.0e+47 | 5.9e-17 |
| 43 | 5.0 | 3.2e+06 | 5.7e-17 |
| 45 | 5.0 | 2.6e+08 | 4.0e-17 |
| 47 | 5.0 | 2.3e+10 | 5.2e-17 |

The failing orders differ between binaries, as expected for an out-of-bounds read. In the libcint bundled with PySCF 2.15.1, m = 41 already fails (1.4e+43 at t = 0.5).

2) Integrals. `libcint_sr_identity_repro.py` (below) checks erfc + erf = full Coulomb for an (AA|BB) quartet of l = 11 shells on two centres (nroots = 23, m = 45), omega = 0.4, exponent 1.0. The table shows max |SR + LR - full| / max |full|.

| d (Angstrom) | 0.3 | 0.5 | 0.6 | 0.7 | 0.9 | other d in 0.3-1.3 |
|---|---|---|---|---|---|---|
| v6.2.1 as released | 4e-12 | 9e-12 | 5e-12 | 8e-12 | 3.1e-05 | <= 3e-11 |
| PySCF 2.15.1 bundled libcint | 3.2e+08 | 5.0e-03 | 8.7e-07 | 9.3e-09 | 3.1e-05 | <= 3e-11 |
| v6.2.1, table extended | 4e-12 | 9e-12 | 5e-12 | 8e-12 | 5e-12 | <= 2e-11 |

With l <= 9 (m <= 37) all three agree at 1e-12. No warning is printed in the failing cases.

<details><summary>libcint_turnover_repro.py</summary>

```python
"""Out-of-bounds read of TURNOVER_POINT[m] in libcint src/fmt.c for m >= 40.

Usage:  python3 libcint_turnover_repro.py /path/to/libcint.so
Needs:  numpy, mpmath.  Calls the exported long-double Boys-like function
fmt_lerfc_like(f, t, lower, m), which the SR Rys-root code calls with
m = 2*nroots - 1 (rys_wheeler.c, lnaive_jacobi_moments), and compares
F_m(t) = int_lower^1 u^(2m) exp(-t u^2) du with a 40-digit quadrature.
"""
import sys, ctypes, mpmath as mp
mp.mp.dps = 40
L = ctypes.CDLL(sys.argv[1])
LD = ctypes.c_longdouble
L.fmt_lerfc_like.argtypes = [ctypes.POINTER(LD), LD, LD, ctypes.c_int]
lower = 0.2
print(" m    t     rel. error of F_m")
for m in (37, 39, 41, 43, 45, 47):
    for t in (0.5, 5.0, 30.0):
        f = (LD * (m + 1))()
        L.fmt_lerfc_like(f, LD(t), LD(lower), m)
        ref = mp.quad(lambda u: u**(2*m) * mp.exp(-t*u*u), [lower, 1])
        err = abs((mp.mpf(float(f[m])) - ref) / ref)
        print("%2d  %5.1f   %.1e" % (m, t, float(err)))
```
</details>

<details><summary>libcint_sr_identity_repro.py</summary>

```python
"""End-to-end check with PySCF: erfc (SR) + erf (LR) = full Coulomb for an
(AA|BB) shell quartet of l = 11 functions (nroots = 23, Boys order m = 45).
Run as-is, or with LD_PRELOAD=/path/to/libcint.so.6 to test a given build."""
import numpy as np
from pyscf import gto
omega, l = 0.4, 11
for d in np.arange(0.3, 1.31, 0.1):
    mol = gto.M(atom=f"He 0 0 0; He 0 0 {d}", basis={"He": [[l, (1.0, 1.0)]]}, verbose=0)
    sl = (0, 1, 0, 1, 1, 2, 1, 2)
    full = mol.intor("int2e_sph", shls_slice=sl)
    with mol.with_range_coulomb(omega):
        lr = mol.intor("int2e_sph", shls_slice=sl)
    with mol.with_range_coulomb(-omega):
        sr = mol.intor("int2e_sph", shls_slice=sl)
    print("d = %.1f A   max|SR + LR - full| / max|full| = %.1e"
          % (d, np.abs(sr + lr - full).max() / np.abs(full).max()))
```
</details>

## Practical impact

It only matters for very high angular momentum: SR integrals need roughly sum of l (plus derivative order) >= 40, i.e. shells with l around 10 or more. Standard basis sets never get there. The users who do are those running basis-set-limit or partial-wave-type calculations, which libcint supports with ANG_MAX = 15. When it does trigger, the integrals are silently wrong.

## Suggested fix

Extend the table to cover every m that can be requested. The formula in the comment above the table, t0 = 0.5*((2m-1)!!/(2m-1)^0.5)^(1/(m-1)), reproduces the existing 38 non-zero entries to 8e-12. The patch below extends the table to m = 71 with values generated from that formula. With it, both reproducers give 1e-17 (Boys values) and 1e-11 (integrals) at every m and distance tested. A guard against m exceeding the table size would also prevent the problem from reappearing if MXRYSROOTS grows.

<details><summary>Patch (src/fmt.c)</summary>

```diff
@@ -79,7 +79,39 @@ static double TURNOVER_POINT[] = {
         13.81992759110,
         14.19235654675,
         14.56467047710,
-        14.93687515212
+        14.93687515212,
+        15.30897591434,
+        15.68097772026,
+        16.05288517693,
+        16.42470257417,
+        16.79643391324,
+        17.16808293229,
+        17.53965312906,
+        17.91114778117,
+        18.2825699643,
+        18.65392256854,
+        19.02520831311,
+        19.39642975963,
+        19.76758932414,
+        20.13868928797,
+        20.50973180762,
+        20.88071892376,
+        21.25165256942,
+        21.62253457746,
+        21.99336668739,
+        22.36415055167,
+        22.73488774143,
+        23.10557975173,
+        23.47622800644,
+        23.8468338627,
+        24.21739861506,
+        24.58792349925,
+        24.95840969574,
+        25.328858333,
+        25.69927049048,
+        26.06964720149,
+        26.43998945575,
+        26.81029820184,
 };
```
</details>

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by an AI model (Claude, Anthropic); all results were checked and the report reviewed by a human before submission.
