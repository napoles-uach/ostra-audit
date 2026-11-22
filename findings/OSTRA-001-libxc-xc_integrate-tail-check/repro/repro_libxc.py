"""Minimal reproductions of the findings from the libxc 7.1.2 probe.

Requires: pip install pyscf  (bundles libxc 7.1.2; this script calls the raw libxc C API via ctypes)
"""
import ctypes, glob, os
import numpy as np
import pyscf

L = ctypes.CDLL(glob.glob(os.path.join(os.path.dirname(pyscf.__file__), "lib", "deps", "lib*", "libxc.so*"))[0])
L.xc_func_alloc.restype = ctypes.c_void_p
L.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
L.xc_version_string.restype = ctypes.c_char_p
P = lambda a: a.ctypes.data_as(ctypes.c_void_p)


def func(fid, nspin=1):
    p = L.xc_func_alloc(); assert L.xc_func_init(p, fid, nspin) == 0; return ctypes.c_void_p(p)


print("libxc", L.xc_version_string().decode())

# 1) LDA_X_1D_SOFT (id 21): exchange energy collapses (and turns positive) for beta*kF > ~38.6
p = func(21)
rho = np.array([10.0, 20.0, 24.5, 24.7, 26.0, 30.0, 100.0]); n = len(rho)
zk, vrho = np.zeros(n), np.zeros(n)
L.xc_lda_exc_vxc(p, ctypes.c_size_t(n), P(rho), P(zk), P(vrho))
print("\n[1] lda_x_1d_soft (beta=1): expected eps_x -> -1/(2 beta) = -0.5 at high density")
for r, z, v in zip(rho, zk, vrho):
    print(f"  rho={r:7.2f}  eps_x={z:+.6e}  v_x={v:+.6e}")

# 2) MGGA_X_PBE_GX (id 576): v_tau jumps at alpha = 1 (kink is in the published form, Loos 2017 Eq. 28)
p = func(576)
r0, s0 = 1.0, 0.5
tauW = s0 / (8 * r0); tauU = 0.3 * (3 * np.pi**2) ** (2 / 3) * r0 ** (5 / 3)
al = np.array([0.99, 0.999, 1.001, 1.01]); n = len(al)
rho, sig, lap, tau = np.full(n, r0), np.full(n, s0), np.zeros(n), tauW + al * tauU
out = [np.zeros(n) for _ in range(5)]
L.xc_mgga_exc_vxc(p, ctypes.c_size_t(n), P(rho), P(sig), P(lap), P(tau), *map(P, out))
print("\n[2] mgga_x_pbe_gx: v_tau across alpha=1 (one-sided dF/dalpha ratio from the paper: ~6.1)")
for a, vt in zip(al, out[4]):
    print(f"  alpha={a:.3f}  v_tau={vt:+.6e}")

# 3) HYB_GGA_XC_HSE06 (id 428): ~2% jitter in v_sigma for s < 0.05 (PBE: ~1e-11)
for fid, nm in ((428, "hse06"), (101, "gga_x_pbe")):
    p = func(fid)
    r0 = 1.0; kf = (3 * np.pi**2 * r0) ** (1 / 3)
    s = np.linspace(0.010, 0.012, 6); n = len(s)
    rho, sig = np.full(n, r0), (2 * kf * r0 * s) ** 2
    zk, vr, vs = np.zeros(n), np.zeros(n), np.zeros(n)
    L.xc_gga_exc_vxc(p, ctypes.c_size_t(n), P(rho), P(sig), P(zk), P(vr), P(vs))
    print(f"\n[3] {nm}: v_sigma along s in [0.010, 0.012] (should vary smoothly)")
    print("  ", np.array2string(vs, precision=8))
