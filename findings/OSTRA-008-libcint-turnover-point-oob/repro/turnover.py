import ctypes, glob, os, pyscf, numpy as np, mpmath as mp
from math import lgamma, exp, log
mp.mp.dps=40
# 1) table formula from the comment: t0 = .5*((2m-1)!!/(2m-1)**.5)**(1/(m-1))
def t0(m):
    df=mp.mpf(1)
    for k in range(1,2*m,2): df*=k
    return float(0.5*(df/mp.sqrt(2*m-1))**(mp.mpf(1)/(m-1)))
print("formula m=2,10,39:",t0(2),t0(10),t0(39),"  (table: 0.866025403784 4.05033123037 14.93687515212)")
print("formula would give m=40..64:",[round(t0(m),3) for m in (40,44,48,56,64)])
# 2) long-double Boys function with m beyond the table (lower != 0 path, as used by SR Jacobi moments)
import sys; L=ctypes.CDLL(sys.argv[1] if len(sys.argv)>1 else glob.glob(os.path.dirname(pyscf.__file__)+'/lib/deps/lib/libcint.so*')[0])
LD=ctypes.c_longdouble
L.fmt_lerfc_like.argtypes=[ctypes.POINTER(LD),LD,LD,ctypes.c_int]
L.lgamma_inc_like.argtypes=[ctypes.POINTER(LD),LD,ctypes.c_int]
def ref(t,low,m):
    t=mp.mpf(t); low=mp.mpf(low)
    return mp.quad(lambda u: u**(2*m)*mp.exp(-t*u*u),[low,1])
for m in (37,39,41,43,45,47):
    for t in (0.5,5.0,30.0):
        f=(LD*(m+1))(); L.fmt_lerfc_like(f,LD(t),LD(0.2),m)
        r=ref(t,0.2,m); print("erfc-like m=%d t=%5.1f  rel err %.1e"%(m,t,abs(float((mp.mpf(float(f[m]))-r)/r))))
