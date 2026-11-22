"""ECPscalar_iprinv returns uninitialized memory when the rinv atom carries no ECP."""
import numpy as np
import pyscf
from pyscf import gto

mol = gto.M(atom="I 0 0 0; H 0 0 3.04", unit="bohr", basis="def2-svp", ecp="def2-svp", verbose=0)
print("pyscf", pyscf.__version__, "| ECP on atoms:", sorted(set(mol._ecpbas[:, gto.ATOM_OF])))

# Atom 1 (H) has no ECP, so d U_H / d R_H is identically zero.
for trial in range(3):
    np.full(200000, 1e10 * (trial + 1)).sum()          # touch some memory between calls
    with mol.with_rinv_at_nucleus(1):
        v = mol.intor("ECPscalar_iprinv", comp=3)
    print(f"call {trial}: max|v| = {np.abs(v).max():.3e}  (expected 0)")

# The kernel does not write to the output buffer at all:
for fill in (0.0, 7.0):
    buf = np.full((3, mol.nao, mol.nao), fill)
    with mol.with_rinv_at_nucleus(1):
        v = mol.intor("ECPscalar_iprinv", comp=3, out=buf)
    print(f"output buffer pre-filled with {fill}: values returned = {np.unique(v)}")
