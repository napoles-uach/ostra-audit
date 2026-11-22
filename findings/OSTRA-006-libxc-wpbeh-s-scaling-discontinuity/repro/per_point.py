"""Size of the s=1 energy jump per grid point in a real molecule (H2O, HSE06)."""
import numpy as np, sys
sys.path.insert(0,'/home/claude/xcprobe'); import lxc
from pyscf import gto, dft
mol=gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587",basis="def2-tzvp",verbose=0)
mf=dft.RKS(mol); mf.xc="HSE06"; mf.grids.level=3; e=mf.kernel()
ao=dft.numint.eval_ao(mol,mf.grids.coords,deriv=1)
rho=dft.numint.eval_rho(mol,ao,mf.make_rdm1(),xctype="GGA")
r,g=rho[0],np.linalg.norm(rho[1:4],axis=0); w=mf.grids.weights
kf=(3*np.pi**2*r)**(1/3); s=g/(2*kf*r+1e-300)
m=r>1e-8
f=lxc.Func(428,1)
def eps(rr,ss):
    sig=(2*(3*np.pi**2*rr)**(1/3)*rr*ss)**2
    return f.exc_vxc({"rho":rr[:,None],"sigma":sig[:,None]})["zk"]
d=1e-9
jump=w[m]*r[m]*(eps(r[m],np.full(m.sum(),1+d))-eps(r[m],np.full(m.sum(),1-d)))  # energy jump if point i sat at s=1
print("E(HSE06)=%.8f  points=%d"%(e,m.sum()))
print("per-point jump |w rho d_eps| at s=1:  max %.2e  median %.2e  sum over all points %.2e Eh"%(np.abs(jump).max(),np.median(np.abs(jump)),np.abs(jump).sum()))
near=np.abs(s[m]-1)<0.01
print("points with |s-1|<0.01: %d  (fraction %.3f)"%(near.sum(),near.mean()))
print("weighted jump of points with |s-1|<0.01: %.2e Eh"%np.abs(jump[near]).sum())
