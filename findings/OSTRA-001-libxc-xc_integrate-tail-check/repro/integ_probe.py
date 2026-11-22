import ctypes, glob, os, numpy as np, pyscf
from scipy.special import k0, k1, modstruve
from scipy.integrate import quad
L = ctypes.CDLL(glob.glob(os.path.join(os.path.dirname(pyscf.__file__), "lib","deps","lib*","libxc.so*"))[0])
L.xc_bessel_K0.restype = ctypes.c_double; L.xc_bessel_K0.argtypes = [ctypes.c_double]
# R-style vectorised integrand: void f(double *x, int n, void *ex)
FN = ctypes.CFUNCTYPE(None, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_void_p)
def mk(power):
    def f(x, n, ex):
        for i in range(n):
            t = x[i]; x[i] = 2*L.xc_bessel_K0(t)*(t**power)
    return FN(f)
f1, f2 = mk(0), mk(1)
L.xc_integrate.restype = ctypes.c_double
L.xc_integrate.argtypes = [FN, ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
exact1 = lambda R: 2*(np.pi*R/2)*(k0(R)*modstruve(-1,R) + k1(R)*modstruve(0,R))   # 2*int_0^R K0
exact2 = lambda R: 2*(1 - R*k1(R))                                                 # 2*int_0^R t K0
print(" R      xc_integrate(2K0)   exact          xc_integrate(2tK0)  exact")
for R in (10, 50, 70, 76, 77, 77.5, 78, 80, 100, 300):
    a1 = L.xc_integrate(f1, None, 0.0, R); a2 = L.xc_integrate(f2, None, 0.0, R)
    print(f"{R:6.1f}  {a1:+.12f}   {exact1(R):+.12f}   {a2:+.12f}   {exact2(R):+.12f}")
