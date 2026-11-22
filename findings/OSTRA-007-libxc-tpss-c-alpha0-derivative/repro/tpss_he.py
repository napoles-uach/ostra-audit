"""One-orbital systems have tau = tau_W exactly (z = 1). Compare RKS vs UKS (same physics, different libxc branches)."""
import numpy as np
from pyscf import gto, dft
for atom,sp in (("He 0 0 0",0),("H 0 0 0; H 0 0 0.74",0)):
    for xc in ("TPSS","PBE","SCAN","M06-L"):
        mol=gto.M(atom=atom,basis="def2-qzvp",verbose=0,spin=sp)
        out=[]
        for M in (dft.RKS,dft.UKS):
            mf=M(mol); mf.xc=xc; mf.grids.level=4; mf.conv_tol=1e-12; out.append(mf.kernel())
        print("%-22s %-6s E_RKS=%.10f  E_UKS-E_RKS=%+.2e Eh"%(atom[:20],xc,out[0],out[1]-out[0]))
