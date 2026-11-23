TITLE: ECP nuclear gradients crash (stack overflow) with i functions: omega_nuc[CART_CUM] in type1_rad_ang is too small

## Summary

`ECPscalar_ipnuc` and `ECPscalar_iprinv`, the integrals used by ECP nuclear gradients, overflow a stack array when the basis contains i functions (l = 6) and the ECP has a local channel (l = -1). The process aborts with `*** stack smashing detected ***` (or segfaults). As a result, `mf.nuc_grad_method()` cannot run with built-in i-function basis sets (e.g. cc-pV5Z-PP, cc-pwCV5Z-PP) combined with ECPs that have a local channel (def2, LANL2DZ, SBKJC, CRENBL...).

Seen in PySCF 2.15.1 (pip). The same code is on master (484fd7c).

## Cause

`CART_CUM` is `455+1` ("upto l = 12"), and `_offset_cart` has 15 entries. `type1_rad_ang` (local channel) declares

```c
double omega_nuc[CART_CUM];
for (i = 0; i <= lmax; i++) {
        pnuc = omega_nuc + _offset_cart[i];
        ang_nuc_in_cart(pnuc, i, unitr);
}
```

It is called from `ECPtype1_cart` with `lmax = li + lj`, which does not depend on the ECP channels. For gradients, `_deriv1_cart` raises `li` by one, so a diagonal block of l = 6 shells gives `lmax = 7 + 6 = 13`, past the end of `omega_nuc`. This is why l = 6 fails for any ECP. The same limit applies to `type2_facs_ang` (non-local channels, `lmax = li + lc`), which can be reached with l = 5 shells and lc = 7.

## Reproducer

```python
from pyscf import gto, scf
mol = gto.M(atom="Ag 0 0 0; H 0 0 1.62",
            basis={"Ag": "cc-pv5z-pp", "H": "sto-3g"},   # cc-pv5z-pp contains i functions
            ecp={"Ag": "def2-svp"},                       # ECP with a local channel
            verbose=0)
mol.intor("ECPscalar")          # works
mol.intor("ECPscalar_ipnuc")    # *** stack smashing detected *** / segfault
# equivalently: scf.RHF(mol).run().nuc_grad_method().kernel()
```

Control: the same ECP with def2-QZVPP (l <= 4) works.

## Suggested fix

Size the table to libcint's `ANG_MAX = 15`:

```diff
--- a/pyscf/lib/gto/nr_ecp.h
+++ b/pyscf/lib/gto/nr_ecp.h
-#define CART_CUM        (455+1) // upto l = 12
+#define CART_CUM        (815+1) // upto l = 15 (ANG_MAX in libcint)
--- a/pyscf/lib/gto/nr_ecp.c
+++ b/pyscf/lib/gto/nr_ecp.c
 static int _offset_cart[] = {0, 1, 4, 10, 20, 35, 56, 84, 120,
-                             165, 220, 286, 364, 455, 560};
+                             165, 220, 286, 364, 455, 560, 680, 816};
```

Tested by rebuilding `libcgto` from the v2.15.1 sources with and without this change:

| Test | Unpatched | Patched |
|---|---|---|
| AgH (cc-pV5Z-PP / def2-SVP), `ECPscalar_ipnuc`, `ECPscalar_iprinv` | stack smashing | runs |
| Synthetic l = 6 shell, local-only ECP and lc = 4; l = 5 shell, lc = 7 (`ECPscalar_ipnuc`) | stack smashing | runs |
| `ECPscalar`, `_ipnuc`, `_iprinv` for (l, lc) = (3,2), (4,3), (5,5), (5,6) | — | bit-identical to unpatched (12/12) |
| AgH RHF gradient on H (z), analytic vs 4-point FD (h = 1e-3 bohr) | aborts | agree to 2.3e-9 Eh/bohr; sum of forces 1.6e-14 |

## Related, minor

`ECPsph_ine_opt` leaves the [0, 16) tables with `z > 16` (nr_ecp.c:4717); for `z == 16.0` exactly, `entry = floor(16.0/K_TAB_INTERVAL) = 400` is one row past the 400-row table. Using `z >= 16` there (and at nr_ecp.c:4663) avoids it. It only affects that exact value.

This is in the same ECP derivative code as #3669 (uninitialized output of the `iprinv` kernels for atoms without ECP).

---
Found and verified as part of OSTRA (Open-Source elecTRonic-structure Audit), a physics-based audit of electronic-structure libraries at UACH. Analysis assisted by AI models (Claude, Anthropic; Gemini, Google); all results were checked and the report reviewed by a human before submission.
