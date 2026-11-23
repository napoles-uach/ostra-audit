TITLE: [BUG REPORT] Segmentation fault with VV10 functionals (wb97m_v, b97m_v, scan_vv10, ...) and the XCFun backend: placeholder ids 5000+ reach libxcfun

**Describe the bug**
With the XCFun backend (`mf._numint.libxc = dft.xcfun`), any `xc` string that names a VV10 functional (`wb97m_v`, `b97m_v`, `wb97x_v`, `lc_vv10`, `scan_vv10`, ..., or `X+vv10`) kills the process with SIGSEGV. The script from #919 now segfaults as well. The same inputs work with the libxc backend, and the documented route `mf.xc = 'b3lyp'; mf.nlc = 'vv10'` works with XCFun and gives the same energy as libxc.

Cause: `xcfun.py:368` assigns temporary ids 5000+ to the `VV10_XC` names so that `parse_xc` can parse NLC functionals (added in #1681). These ids are not filtered before being passed to the C interface. The crash happens in `xc_type` (`xcfun.py:338`), before any evaluation:

```
rks.get_veff -> numint.nr_rks -> numint._xc_type -> xcfun.xc_type -> _itrf.XCFUN_xc_type(5001)
```

`XCFUN_xc_type` (`lib/dft/xcfun_itrf.c:170`) calls `xcfun_enumerate_parameters(5001)`, which returns NULL for an out-of-range id, and then `xcfun_set(fun, NULL, 1.)`. The same ids would also reach `XCFUN_eval_xc` from `eval_xc1`.

**To Reproduce**
```python
from pyscf import gto, dft
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
mf = dft.RKS(mol)
mf.xc = "wb97m_v"                 # also b97m_v, wb97x_v, lc_vv10, scan_vv10, b3lyp+vv10
mf._numint.libxc = dft.xcfun
mf.kernel()                       # Segmentation fault (exit 139)
```

**Expected behavior**
A Python exception with a clear message (as before #1681, when `parse_xc` raised `KeyError: 'Unknown XC functional WB97M_V'`), or a correct result.

**Error message**
```
Fatal Python error: Segmentation fault
  File ".../pyscf/dft/xcfun.py", line 338 in <genexpr>
  File ".../pyscf/dft/xcfun.py", line 338 in xc_type
  File ".../pyscf/dft/numint.py", line 2815 in _xc_type
  File ".../pyscf/dft/numint.py", line 1117 in nr_rks
  File ".../pyscf/dft/rks.py", line 81 in get_veff
```

**Version and hardware (please complete the following information):**
- PySCF version: 2.15.1 (`dft/xcfun.py`, `dft/numint.py` and `lib/dft/xcfun_itrf.c` unchanged on master 6b288fd); first seen with 2.14.0
- PySCF install type: PyPI/pip
- Operating system: Linux x86_64 (2.15.1) and macOS arm64 (2.14.0)

**Additional context**
Skipping the ids >= 5000 is not enough. In the XCFun backend a named VV10 functional maps *only* to its placeholder; its semilocal part is not represented. With the placeholders simply filtered out, the calculation runs without any exchange-correlation and silently returns a wrong energy (H2/cc-pVDZ):

| xc | libxc backend | XCFun, ids >= 5000 skipped |
|---|---|---|
| `scan_vv10` | -1.16323561 | -0.52497908 |
| `lc_vv10` | -1.16626798 | -0.52497908 |
| `b97m_v` | -1.16513535 | -0.52497908 |
| `wb97m_v` | -1.15338673 | -0.52497908 |

So either the semilocal parts are mapped in the XCFun backend, or the backend should refuse these strings. A minimal patch for the second option:

<details><summary>Patch (pyscf/dft/xcfun.py)</summary>

```diff
@@ -324,6 +324,13 @@ VV10_XC = {
 }
 VV10_XC.update([(key.replace('_', ''), val) for key, val in VV10_XC.items()])
 
+def _check_no_nlc(xc_code, fn_facs):
+    '''NLC placeholder ids (>= 5000) are not XCFun functionals and must not reach libxcfun_itrf.'''
+    if any(xid >= 5000 for xid, fac in fn_facs):
+        raise NotImplementedError(
+            f'{xc_code} contains a VV10 functional. The semilocal part of VV10 functionals is not '
+            'available in the XCFun backend. Use the libxc backend, or a semilocal xc with mf.nlc.')
+
 @lru_cache(100)
 def xc_type(xc_code):
     if xc_code is None:
@@ -332,6 +339,7 @@ def xc_type(xc_code):
         hyb, fn_facs = parse_xc(xc_code)
     else:
         fn_facs = [(xc_code, 1)]  # mimic fn_facs
+    _check_no_nlc(xc_code, fn_facs)
 
     if not fn_facs:
         return 'HF'
@@ -1037,6 +1045,7 @@ def eval_xc1(xc_code, rho, spin=0, deriv=1, omega=None):
         rho = numpy.asarray(rho[...,[0,1,2,3,5],:], order='C')
 
     hyb, fn_facs = parse_xc(xc_code)
+    _check_no_nlc(xc_code, fn_facs)
     if omega is not None:
         hyb = hyb[:2] + (float(omega),)
```
</details>

Tested on 2.15.1 with this patch (H2/cc-pVDZ, XCFun backend):

| Input | Before | After |
|---|---|---|
| `wb97m_v`, `scan_vv10`, `b3lyp+vv10`; script from #919 | segfault | `NotImplementedError` |
| `b3lyp` + `mf.nlc = 'vv10'` | -1.16416231 | -1.16416231 (= libxc) |
| `scan` + `mf.nlc = 'vv10'` | -1.15674014 | -1.15674014 (= libxc) |
| `b3lyp`, `pbe` | unchanged | unchanged |
| `pyscf/dft/test/test_xcfun.py` | — | 11 passed |

Related: #919 (same symptom area, fixed in #928); the placeholder ids were introduced later in #1681.

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
