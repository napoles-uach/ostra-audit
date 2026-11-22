import numpy as np
from pyscf import gto
from scipy.spatial.transform import Rotation
L="SPDFGHIKLM"
def mk(geom, shells):
    atoms=[[f"X{i}",c] for i,c in enumerate(geom)]
    basis={f"X{i}":gto.basis.parse("".join(f"X {L[l]}\n {e:.17g} 1.0\n" for l,e in sh)) for i,sh in enumerate(shells)}
    m=gto.Mole(atom=atoms,basis=basis,unit="bohr",verbose=0); m.build(dump_input=False,parse_arg=False); return m
for lmax in (4,5,6,7):
    worst=(0,None)
    for R in np.linspace(0.25,12,48):
        for a,b in ((0.9,1.7),(0.3,3.0)):
            g=np.array([[0,0,0],[0.3*R,-0.2*R,0.93*R]])
            sh=[[(lmax,a)],[(lmax,b)]]
            m=mk(g,sh); n=m.nao
            try: e=m.intor("int2e",aosym="s1")
            except Exception as ex: print("fail",lmax,ex); break
            v=np.linalg.eigvalsh(e.reshape(n*n,n*n))
            Rm=Rotation.random(random_state=int(R*10)).as_matrix()
            v2=np.linalg.eigvalsh(mk(g@Rm.T,sh).intor("int2e",aosym="s1").reshape(n*n,n*n))
            d=max(-v[0]/v[-1], np.max(np.abs(v2-v))/v[-1])
            if d>worst[0]: worst=(d,(R,a,b))
    print(f"l={lmax}: worst (neg-eig or rotation dev)/max eig = {worst[0]:.2e} at R,a,b={worst[1]}",flush=True)
