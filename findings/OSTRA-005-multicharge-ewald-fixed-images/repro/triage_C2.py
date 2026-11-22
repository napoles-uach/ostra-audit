import numpy as np
from dftd4.interface import DispersionModel as M, DampingParam as P
P3=np.array([True]*3)
a=10.66
prim=0.5*a*np.array([[0,1,1],[1,0,1],[1,1,0]]); n=np.array([11,17]); p=np.array([[0,0,0],[0.5,0.5,0.5]])@prim
def q_e(n,p,lat):
    m=M(n,p,lattice=lat,periodic=P3); return m.get_properties()["partial charges"], m.get_dispersion(P(method="pbe0"),grad=False)["energy"]/len(n)
q0,e0=q_e(n,p,prim)
# conventional cubic cell (8 atoms)
conv=a*np.eye(3); fr=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
nc=np.r_[[11]*4,[17]*4]; pc=np.r_[fr,fr+[.5,0,0]]@conv
qc,ec=q_e(nc,pc,conv)
print(f"primitive: q(Na)={q0[0]:.8f}  E/atom={e0:.10e}")
print(f"conventional cubic: q(Na)={qc[nc==11].mean():.8f} (spread {np.ptp(qc[nc==11]):.1e})  E/atom={ec:.10e}")
for k in range(1,6):
    T=np.array([[1,k,0],[0,1,0],[0,0,1]])  # increasingly skewed equivalent cells
    lat=T@prim; ang=np.degrees(np.arccos(lat[0]@lat[1]/np.linalg.norm(lat[0])/np.linalg.norm(lat[1])))
    q,e=q_e(n,p,lat)
    print(f"skew k={k}: |a'|={np.linalg.norm(lat[0]):6.1f} bohr, angle(a',b)={ang:5.1f} deg  q(Na)={q[0]:.8f}  dq={q[0]-qc[nc==11].mean():+.2e}  dE/atom rel={abs(e-ec)/abs(ec):.2e}")
