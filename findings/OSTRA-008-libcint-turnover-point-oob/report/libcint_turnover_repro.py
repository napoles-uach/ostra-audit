"""Out-of-bounds read of TURNOVER_POINT[m] in libcint src/fmt.c for m >= 40.

Usage:  python3 libcint_turnover_repro.py /path/to/libcint.so
Needs:  numpy, mpmath.  Calls the exported long-double Boys-like function
fmt_lerfc_like(f, t, lower, m), which the SR Rys-root code calls with
m = 2*nroots - 1 (rys_wheeler.c, lnaive_jacobi_moments), and compares
F_m(t) = int_lower^1 u^(2m) exp(-t u^2) du with a 40-digit quadrature.
"""
import sys, ctypes, mpmath as mp
mp.mp.dps = 40
L = ctypes.CDLL(sys.argv[1])
LD = ctypes.c_longdouble
L.fmt_lerfc_like.argtypes = [ctypes.POINTER(LD), LD, LD, ctypes.c_int]
lower = 0.2
print(" m    t     rel. error of F_m")
for m in (37, 39, 41, 43, 45, 47):
    for t in (0.5, 5.0, 30.0):
        f = (LD * (m + 1))()
        L.fmt_lerfc_like(f, LD(t), LD(lower), m)
        ref = mp.quad(lambda u: u**(2*m) * mp.exp(-t*u*u), [lower, 1])
        err = abs((mp.mpf(float(f[m])) - ref) / ref)
        print("%2d  %5.1f   %.1e" % (m, t, float(err)))
