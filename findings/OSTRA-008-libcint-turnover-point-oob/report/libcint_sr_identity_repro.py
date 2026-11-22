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
