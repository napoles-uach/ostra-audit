"""Energy along a small H displacement: HSE06 (has s=1 jump) vs PBE0 (control). Steps from points crossing s=1."""
import numpy as np, sys, json
from pyscf import gto, dft
xc=sys.argv[1]; n=int(sys.argv[2]); h=float(sys.argv[3])
base=np.array([[0,0,0],[0,0.757,0.587],[0,-0.757,0.587]])
E=[]; dm=None
for k in range(n):
    c=base.copy(); c[1,1]+=k*h
    mol=gto.M(atom=[["O",c[0]],["H",c[1]],["H",c[2]]],basis="def2-svp",verbose=0,unit="Ang")
    mf=dft.RKS(mol); mf.xc=xc; mf.grids.level=3; mf.conv_tol=1e-12
    E.append(mf.kernel(dm0=dm)); dm=mf.make_rdm1()
E=np.array(E); x=np.arange(n)*h
p=np.polyfit(x,E,4); res=E-np.polyval(p,x)
json.dump({"xc":xc,"h":h,"E":E.tolist(),"res":res.tolist()},open(f"scan_{xc}.json","w"))
print(xc,"rms residual %.2e Eh, max |res| %.2e, max |2nd diff| %.2e"%(res.std(),np.abs(res).max(),np.abs(np.diff(E,2)).max()))
