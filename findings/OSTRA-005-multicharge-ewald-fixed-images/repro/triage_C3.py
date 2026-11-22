import numpy as np
from dftd4.interface import DispersionModel as M, DampingParam as P
P3=np.array([True]*3)
def lll(B, delta=0.75):
    """LLL-reduce lattice rows (bohr); returns reduced basis spanning the same lattice"""
    B=B.copy().astype(float); n=len(B); k=1
    def gs(B):
        Bs=np.zeros_like(B); mu=np.zeros((n,n))
        for i in range(n):
            Bs[i]=B[i]
            for j in range(i):
                mu[i,j]=B[i]@Bs[j]/(Bs[j]@Bs[j]); Bs[i]-=mu[i,j]*Bs[j]
        return Bs,mu
    Bs,mu=gs(B)
    while k<n:
        for j in range(k-1,-1,-1):
            q=round(mu[k,j])
            if q: B[k]-=q*B[j]; Bs,mu=gs(B)
        if Bs[k]@Bs[k] >= (delta-mu[k,k-1]**2)*(Bs[k-1]@Bs[k-1]): k+=1
        else: B[[k,k-1]]=B[[k-1,k]]; Bs,mu=gs(B); k=max(k-1,1)
    return B
a=10.66; prim=0.5*a*np.array([[0,1,1],[1,0,1],[1,1,0]]); n=np.array([11,17]); p=np.array([[0,0,0],[0.5,0.5,0.5]])@prim
ref=0.34806642
for k in (1,2,3,5):
    lat=np.array([[1,k,0],[0,1,0],[0,0,1]])@prim
    q_raw=M(n,p,lattice=lat,periodic=P3).get_properties()["partial charges"][0]
    red=lll(lat)
    q_red=M(n,p,lattice=red,periodic=P3).get_properties()["partial charges"][0]
    print(f"k={k}: raw q(Na)={q_raw:.8f} (dq {q_raw-ref:+.1e})   after LLL reduction q(Na)={q_red:.8f} (dq {q_red-ref:+.1e})")
