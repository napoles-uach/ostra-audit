import numpy as np, json
from pyscf import gto, dft
B=1/0.52917721092
for xc in ("HSE06","PBE0"):
    d=json.load(open(f"scan_{xc}.json")); E=np.array(d["E"]); h=d["h"]
    fd=(E[2:]-E[:-2])/(2*h/ 0.52917721092)   # Eh/bohr along H1 y
    devs=[]
    for k in (5,15,25,35):
        c=np.array([[0,0,0],[0,0.757+k*h,0.587],[0,-0.757,0.587]])
        mol=gto.M(atom=[["O",c[0]],["H",c[1]],["H",c[2]]],basis="def2-svp",verbose=0,unit="Ang")
        mf=dft.RKS(mol); mf.xc=xc; mf.grids.level=3; mf.conv_tol=1e-12; mf.kernel()
        g=mf.nuc_grad_method(); g.grid_response=True; gy=g.kernel()[1,1]
        devs.append(abs(gy-fd[k-1]))
    print(xc,"|analytic - FD(h=0.0005 A)| along scan: ",["%.1e"%x for x in devs],"Eh/bohr")
