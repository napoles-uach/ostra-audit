# Explorer 1, 2026-11-23. Patched libcint v6.2.1 (table extended to m = 71), Linux x86_64 (80-bit long double).
# Max relative error of F_m just above t0(m) (T in [1.0001, 1.3] * t0) vs 50-digit quadrature.
import ctypes, numpy as np, mpmath as mp
mp.mp.dps = 50
lib = ctypes.CDLL('libcint.so')   # patched build
lib.gamma_inc_like.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
lib.fmt_erfc_like.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.lgamma_inc_like.argtypes = [ctypes.c_void_p, ctypes.c_longdouble, ctypes.c_int]
def t0(m): return 0.5*float((mp.fac2(2*m-1)/mp.sqrt(2*m-1))**(mp.mpf(1)/(m-1)))
def ref(m, T, lo=0): return mp.quad(lambda u: u**(2*m)*mp.exp(-T*u*u), [lo, 1])
for m in (39, 41, 45, 47, 55, 63, 71):
    e1 = e2 = e3 = 0
    for T in np.linspace(t0(m)*1.0001, t0(m)*1.3, 15):
        f = np.zeros(m+1); lib.gamma_inc_like(f.ctypes.data, T, m); e1 = max(e1, abs(f[m]-float(ref(m, T)))/float(ref(m, T)))
        f = np.zeros(m+1); lib.fmt_erfc_like(f.ctypes.data, T, 0.2, m); e2 = max(e2, abs(f[m]-float(ref(m, T, .2)))/float(ref(m, T, .2)))
        g = np.zeros(m+1, dtype=np.longdouble); lib.lgamma_inc_like(g.ctypes.data, T, m)
        e3 = max(e3, float(abs(mp.mpf(str(g[m]))-ref(m, T))/ref(m, T)))
    print(f"m={m:2d} t0={t0(m):6.2f} double gamma {e1:.1e}  double erfc(.2) {e2:.1e}  long double gamma {e3:.1e}")
