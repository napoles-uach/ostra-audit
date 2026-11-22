import numpy as np, lxc
from pyscf import gto, dft
mol=gto.M(atom="He 0 0 0",basis="def2-qzvp",verbose=0)
mf=dft.RKS(mol); mf.xc="TPSS"; mf.grids.level=4; mf.kernel()
ao=dft.numint.eval_ao(mol,mf.grids.coords,deriv=1)
rho=dft.numint.eval_rho(mol,ao,mf.make_rdm1(),xctype="MGGA")
r,sig,tau=rho[0],(rho[1:4]**2).sum(0),rho[5]
m=r>1e-8; r,sig,tau,w=r[m],sig[m],tau[m],mf.grids.weights[m]
z=sig/(8*r*tau)
print("points",m.sum()," z>1:",(z>1).sum()," z==1:",(z==1).sum()," |z-1|<1e-12:",(abs(z-1)<1e-12).sum(), " min z %.6f"%z.min())
fu,fp=lxc.Func(231,1),lxc.Func(231,2)
xu={"rho":r[:,None],"sigma":sig[:,None],"tau":tau[:,None],"lapl":np.zeros((len(r),1))}
xp={"rho":np.stack([r/2]*2,1),"sigma":np.stack([sig/4]*3,1),"tau":np.stack([tau/2]*2,1),"lapl":np.zeros((len(r),2))}
ou,op=fu.exc_vxc(xu),fp.exc_vxc(xp)
dvt=np.abs(ou["vtau"][:,0]-op["vtau"].mean(1)); 
print("points where vtau differs (rel>1e-6):",(dvt>1e-6*np.abs(ou["vtau"][:,0])+1e-300).sum())
print("vtau==0 in unp:",(ou["vtau"][:,0]==0).sum()," in pol:",(op["vtau"]==0).all(1).sum())
