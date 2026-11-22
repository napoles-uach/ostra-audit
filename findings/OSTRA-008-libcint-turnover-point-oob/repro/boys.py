"""Accuracy of libcint gamma_inc_like (Boys function F_m(t)) against mpmath."""
import ctypes, numpy as np, mpmath as mp, glob, os, pyscf
mp.mp.dps=50
L=ctypes.CDLL(glob.glob(os.path.dirname(pyscf.__file__)+'/lib/deps/lib/libcint.so*')[0])
L.gamma_inc_like.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.c_double,ctypes.c_int]
def cint(t,m):
    f=(ctypes.c_double*(m+1))(); L.gamma_inc_like(f,t,m); return np.array(f[:])
def ref(t,m):
    t=mp.mpf(t)
    return [0.5*t**(-(i+mp.mpf(0.5)))*mp.gammainc(i+mp.mpf(0.5),0,t) for i in range(m+1)]
M=40
tt=np.logspace(-8,2.3,300)
rel=np.zeros((len(tt),M+1))
for k,t in enumerate(tt):
    f=cint(float(t),M); r=ref(float(t),M)
    rel[k]=[abs(float((mp.mpf(f[i])-r[i])/r[i])) for i in range(M+1)]
for m in (0,1,2,4,8,12,16,20,24,30,35,40):
    i=np.argmax(rel[:,m]); print("m=%2d  max rel err %.2e at t=%.4g"%(m,rel[i,m],tt[i]))
np.save("boys_rel.npy",rel); np.save("boys_t.npy",tt)
