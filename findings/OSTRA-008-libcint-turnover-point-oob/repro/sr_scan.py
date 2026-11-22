import numpy as np, sys
from pyscf import gto
omega=0.4; l=int(sys.argv[1])
for d in np.arange(0.3,1.31,0.1):
    mol=gto.M(atom=f"He 0 0 0; He 0 0 {d}",basis={"He":[[l,(1.0,1.0)]]},verbose=0)
    sl=(0,1,0,1,1,2,1,2)
    full=mol.intor("int2e_sph",shls_slice=sl)
    with mol.with_range_coulomb(omega): lr=mol.intor("int2e_sph",shls_slice=sl)
    with mol.with_range_coulomb(-omega): sr=mol.intor("int2e_sph",shls_slice=sl)
    print("l=%d d=%.1f  err=%.1e"%(l,d,np.abs(sr+lr-full).max()/np.abs(full).max()),flush=True)
