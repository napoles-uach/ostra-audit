"""RHF nuclear gradient of AgH (Ag cc-pV5Z-PP with i functions, def2-SVP ECP; H sto-3g): analytic vs central FD on H z."""
import numpy as np
from pyscf import gto, scf
def mol_at(z): return gto.M(atom=f"Ag 0 0 0; H 0 0 {z}", basis={"Ag":"cc-pv5z-pp","H":"sto-3g"}, ecp={"Ag":"def2-svp"}, unit="Bohr", verbose=0)
z0=3.06; m=mol_at(z0); mf=scf.RHF(m); mf.conv_tol=1e-12; e0=mf.kernel()
g=mf.nuc_grad_method().kernel()
h=1e-3; es=[]
for d in (-2,-1,1,2):
    mf2=scf.RHF(mol_at(z0+d*h)); mf2.conv_tol=1e-12; es.append(mf2.kernel(dm0=mf.make_rdm1()))
fd=(es[0]-8*es[1]+8*es[2]-es[3])/(12*h)
print("E=%.10f  analytic dE/dz(H)=%.10e  FD4=%.10e  |diff|=%.2e  sum of forces=%.1e"%(e0,g[1,2],fd,abs(g[1,2]-fd),abs(g.sum(0)).max()))
