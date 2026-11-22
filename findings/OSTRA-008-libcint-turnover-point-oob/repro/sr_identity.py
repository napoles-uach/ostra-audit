"""Oracle: erfc (SR) + erf (LR) = full Coulomb, shell quartet (A A|B B) with high angular momentum l.
nroots = 2l + 1 ; the SR path builds Boys values F_m with m = 2*nroots - 1 from TURNOVER_POINT[m] (40 entries)."""
import numpy as np, sys
from pyscf import gto
omega=0.4
for l in (8,9,10,11):
    for d in (0.5,1.0,2.0):
        mol=gto.M(atom=f"He 0 0 0; He 0 0 {d}",basis={"He":[[l,(1.0,1.0)]]},cart=False,verbose=0)
        sl=(0,1,0,1,1,2,1,2)
        full=mol.intor("int2e_sph",shls_slice=sl)
        with mol.with_range_coulomb(omega): lr=mol.intor("int2e_sph",shls_slice=sl)
        with mol.with_range_coulomb(-omega): sr=mol.intor("int2e_sph",shls_slice=sl)
        err=np.abs(sr+lr-full).max()/np.abs(full).max()
        print("l=%2d nroots=%2d m_max=%2d d=%.1f A   max|SR+LR-full|/max|full| = %.1e   max|SR|=%.2e max|full|=%.2e"%(l,2*l+1,2*(2*l+1)-1,d,err,np.abs(sr).max(),np.abs(full).max()),flush=True)
