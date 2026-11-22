#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidato #4:
Singularidad analitica (1 - |zeta|)^(-1/3) de TPSS Ec. (34) y corte duro en `tpss_C0`
en 1 - |zeta| = 1e-12 (maple/tpss_c.mpl:34-36), afectando a mgga_c_tpss (231) y hyb_mgga_xc_tpssh (457).
"""
import ctypes, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
LIB_PATH = os.environ.get("LIBXC_PATH", str(ROOT / "build-libxc" / "libxc.dylib"))
L = ctypes.CDLL(LIB_PATH)

L.xc_version_string.restype = ctypes.c_char_p
L.xc_func_alloc.restype = ctypes.c_void_p
L.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
L.xc_func_end.argtypes = [ctypes.c_void_p]
L.xc_func_free.argtypes = [ctypes.c_void_p]
L.xc_func_set_dens_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_func_set_sigma_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_func_set_zeta_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_func_set_tau_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_mgga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 9


def eval_tpss_pol(fid, omz, zeta_thr=1e-20):
    rtot = 0.2
    ra, rb = 0.5 * rtot * (2.0 - omz), 0.5 * rtot * omz
    g = 0.1  # gradientes proporcionales => grad(zeta) = 0 exacto, pero d(xi^2)/d(sigma_bb) = O(1) != 0
    saa, sab, sbb = (ra / rtot) ** 2 * g, (ra * rb / rtot**2) * g, (rb / rtot) ** 2 * g
    ta = saa / (8.0 * ra) + 0.2 * (ra ** (5.0 / 3.0))
    tb = sbb / (8.0 * rb) + 0.2 * (rb ** (5.0 / 3.0))
    p = L.xc_func_alloc()
    assert L.xc_func_init(p, fid, 2) == 0
    L.xc_func_set_dens_threshold(p, 1e-30)
    L.xc_func_set_sigma_threshold(p, 1e-40)
    L.xc_func_set_tau_threshold(p, 1e-40)
    L.xc_func_set_zeta_threshold(p, zeta_thr)
    rho = np.array([[ra, rb]]); sig = np.array([[saa, sab, sbb]])
    lap = np.array([[0.0, 0.0]]); tau = np.array([[ta, tb]])
    zk = np.zeros(1); vr = np.zeros((1, 2)); vs = np.zeros((1, 3)); vl = np.zeros((1, 2)); vt = np.zeros((1, 2))
    L.xc_mgga_exc_vxc(p, 1, rho.ctypes.data, sig.ctypes.data, lap.ctypes.data, tau.ctypes.data,
                      zk.ctypes.data, vr.ctypes.data, vs.ctypes.data, vl.ctypes.data, vt.ctypes.data)
    L.xc_func_end(p); L.xc_func_free(p)
    return zk[0], vr[0], vs[0]


if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")
    for label, fid in [("mgga_c_tpss (id=231)", 231), ("hyb_mgga_xc_tpssh (id=457)", 457)]:
        print(f"=== {label} alrededor de 1 - |zeta| = 1e-12 (rho_tot=0.2, grad_zeta=0) ===")
        for omz in (1e-12 - 1e-14, 1e-12 - 1e-15, 1e-12 + 1e-15, 1e-12 + 1e-14, 1e-9, 1e-6):
            zk, vr, vs = eval_tpss_pol(fid, omz)
            print(f"  1-|z|={omz:.6e}: zk={zk:.15e}, vrho_b={vr[1]:10.4f}, vsigma_bb={vs[2]:12.6f}")
        print()
