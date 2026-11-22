"""HSE06 scan + post-hoc removal of ONLY the s=1 jump: E_cont = E - sum_{s_i>=1} w_i rho_i J(rho_i),
J(rho) = eps(rho,1+) - eps(rho,1-).  If the residual collapses, the s=1 switch is the cause."""
import numpy as np, sys, json
sys.path.insert(0,'/home/claude/xcprobe'); import lxc
from pyscf import gto, dft
f=lxc.Func(428,1)
def eps(rr,ss):
    sig=(2*(3*np.pi**2*rr)**(1/3)*rr*ss)**2
    return f.exc_vxc({"rho":rr[:,None],"sigma":sig[:,None]})["zk"]
n,h=41,0.0005
base=np.array([[0,0,0],[0,0.757,0.587],[0,-0.757,0.587]])
E=[];C=[];N=[];dm=None
for k in range(n):
    c=base.copy(); c[1,1]+=k*h
    mol=gto.M(atom=[["O",c[0]],["H",c[1]],["H",c[2]]],basis="def2-svp",verbose=0,unit="Ang")
    mf=dft.RKS(mol); mf.xc="HSE06"; mf.grids.level=3; mf.conv_tol=1e-12
    e=mf.kernel(dm0=dm); dm=mf.make_rdm1()
    ao=dft.numint.eval_ao(mol,mf.grids.coords,deriv=1); rho=dft.numint.eval_rho(mol,ao,dm,xctype="GGA")
    r=rho[0]; g=np.linalg.norm(rho[1:4],axis=0); m=r>1e-10
    s=np.zeros_like(r); s[m]=g[m]/(2*(3*np.pi**2*r[m])**(1/3)*r[m])
    up=m&(s>=1)
    J=eps(r[up],np.full(up.sum(),1+1e-9))-eps(r[up],np.full(up.sum(),1-1e-9))
    E.append(e); C.append(e-np.sum(mf.grids.weights[up]*r[up]*J)); N.append(int(up.sum()))
E,C=np.array(E),np.array(C); x=np.arange(n)*h
for lab,y in (("HSE06 as is",E),("HSE06 minus s=1 jump",C)):
    res=y-np.polyval(np.polyfit(x,y,4),x); print("%-22s rms residual %.2e Eh  max %.2e"%(lab,res.std(),np.abs(res).max()))
print("grid points with s>=1 along the scan: min %d max %d"%(min(N),max(N)))
json.dump({"E":E.tolist(),"C":C.tolist(),"N":N,"h":h},open("scan_corr.json","w"))
