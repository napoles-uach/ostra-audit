**Title:** `ECPscalar_iprinv` (and `_ipiprinv`, `_iprinvip`) return uninitialized memory when the rinv atom carries no ECP

### Summary

The ECP-centre derivative kernels in `pyscf/lib/gto/nr_ecp_deriv.c` return `0` without writing to the output buffer when `AS_RINV_ORIG_ATOM` points to an atom that has no ECP. `getints2c` allocates the result with `numpy.ndarray` (uninitialized), so the caller receives whatever was in that memory instead of zeros.

PySCF itself is not affected: `grad/rhf.py` and `hessian/rhf.py` only call these integrals for atoms listed in `ecp_atoms`. Any direct caller that loops over all atoms gets silently wrong, non-reproducible derivatives. This includes custom gradient code and consumers of the vendored copy of this file (e.g. the Rust crate `libcint-src`, `src/cecp/nr_ecp_deriv.c`).

Affected functions (6): `ECPscalar_iprinv_{cart,sph}`, `ECPscalar_ipiprinv_{cart,sph}`, `ECPscalar_iprinvip_{cart,sph}`.

### Reproducer (PySCF 2.15.1; `nr_ecp_deriv.c` identical on current `master`)

```python
import numpy as np
from pyscf import gto

mol = gto.M(atom="I 0 0 0; H 0 0 3.04", unit="bohr", basis="def2-svp", ecp="def2-svp", verbose=0)
# atom 1 (H) has no ECP -> d U_H / d R_H must be zero
with mol.with_rinv_at_nucleus(1):
    v = mol.intor("ECPscalar_iprinv", comp=3)
print(np.abs(v).max())            # 7.147e+256 here; varies between runs

for fill in (0.0, 7.0):           # the kernel leaves the buffer untouched
    buf = np.full((3, mol.nao, mol.nao), fill)
    with mol.with_rinv_at_nucleus(1):
        print(np.unique(mol.intor("ECPscalar_iprinv", comp=3, out=buf)))   # [0.] then [7.]
```

### Cause

For example, lines 344–349 of `ECPscalar_iprinv_cart`:

```c
int atm_id = (int)env[AS_RINV_ORIG_ATOM];
...
int shl_id = _one_shell_ecpbas(&necpbas, atm_id, atm, natm, bas, nbas, env);
if (shl_id < 0) {
        return 0;          // output block is never written
}
```

Every other path that has no value calls `ECPscalar_distribute0(out, dims, comp, di, dj)` (e.g. `_cart_factory`, line 324), so the drivers rely on the kernel always filling its block.

### Proposed fix

Zero the block before the early return, as `_cart_factory` does. For the `_cart` variants:

```c
if (shl_id < 0) {
        const int li = bas[ANG_OF+shls[0]*BAS_SLOTS];
        const int lj = bas[ANG_OF+shls[1]*BAS_SLOTS];
        const int di = (li+1)*(li+2)/2 * bas[NCTR_OF+shls[0]*BAS_SLOTS];
        const int dj = (lj+1)*(lj+2)/2 * bas[NCTR_OF+shls[1]*BAS_SLOTS];
        ECPscalar_distribute0(out, dims, 3, di, dj);    /* 9 for ipiprinv / iprinvip */
        return 0;
}
```

The `_sph` variants need the same change with `(2*l+1)` in place of `(l+1)(l+2)/2`. A regression test can check that `ECPscalar_iprinv` with `with_rinv_at_nucleus` on a non-ECP atom returns zeros when the output buffer is pre-filled with a non-zero value.

---
*Found during a systematic, AI-assisted audit (Claude) of numerical libraries used in DFT codes. The finding was verified by the reproducer above and by reading the source.*
