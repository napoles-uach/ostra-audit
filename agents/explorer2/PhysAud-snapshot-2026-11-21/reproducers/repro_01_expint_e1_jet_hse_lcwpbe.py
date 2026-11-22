#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidato #1:
Cambio de rama asintotico prematuro en `xc_E1_scaled_jet` (src/expint_e1.c:140).

Afecta en prioridad alta a: HSE06 (428), HSE03 (427), HSE12 (479), LC-wPBE (478), gga_x_wpbeh (524).
Ejecucion: LIBXC_PATH=/ruta/a/libxc.dylib python3 repro_01_expint_e1_jet_hse_lcwpbe.py
"""
import ctypes, math, os
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
L.xc_gga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 5

PI = math.pi
X2S = 1.0 / (2.0 * (6.0 * PI**2) ** (1.0 / 3.0))


def eval_gga(fid, rho_arr, sig_arr):
    p = L.xc_func_alloc()
    assert L.xc_func_init(p, fid, 1) == 0
    L.xc_func_set_dens_threshold(p, 1e-30)
    L.xc_func_set_sigma_threshold(p, 1e-40)
    rho = np.ascontiguousarray(rho_arr, dtype=np.float64)
    sig = np.ascontiguousarray(sig_arr, dtype=np.float64)
    n = len(rho)
    zk, vrho, vsig = np.zeros(n), np.zeros(n), np.zeros(n)
    L.xc_gga_exc_vxc(p, n, rho.ctypes.data, sig.ctypes.data, zk.ctypes.data, vrho.ctypes.data, vsig.ctypes.data)
    L.xc_func_end(p); L.xc_func_free(p)
    return zk, vrho, vsig


def rho_sigma_for_aux5(aux5, omega, s=0.2):
    Ha1, Ha2, Ha3, Ha4, Ha5 = 0.00979681, 0.0410834, 0.187440, 0.00120824, 0.0347188
    Hs = (Ha1 * s**2 + Ha2 * s**4) / (1.0 + Ha3 * s**4 + Ha4 * s**5 + Ha5 * s**6)
    aux4 = aux5 * 4.0 * 1.0161144 / 9.0
    w = math.sqrt((aux4 - s**2 * Hs) / 1.455915450052607)
    rho = ((omega / w) ** 3) / (3.0 * PI**2)
    sig = 4.0 * ((s / X2S) * (rho / 2.0) ** (4.0 / 3.0)) ** 2
    return rho, sig


if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")
    if hasattr(L, "xc_E1_scaled_d1"):
        L.xc_E1_scaled_d1.argtypes = [ctypes.c_double]
        L.xc_E1_scaled_d1.restype = ctypes.c_double
        for d in (1e-9, 1e-14):
            dm, dp = L.xc_E1_scaled_d1(18.0 - d), L.xc_E1_scaled_d1(18.0 + d)
            print(f"xc_E1_scaled_d1(18 +- {d:.0e}): d1(-)={dm:.15e}, d1(+)={dp:.15e}, salto_rel={abs(dp-dm)/abs(dm):.3e}")
        print()

    for name, fid, omega in [("HSE06 (id=428)", 428, 0.11), ("LC-wPBE (id=478)", 478, 0.40)]:
        print(f"=== {name}, omega={omega}, s=0.2 ===")
        for d in (1e-9, 1e-13):
            rm, sm = rho_sigma_for_aux5(18.0 - d, omega)
            rp, sp = rho_sigma_for_aux5(18.0 + d, omega)
            zkm, vrm, _ = eval_gga(fid, [rm], [sm])
            zkp, vrp, _ = eval_gga(fid, [rp], [sp])
            print(f"  aux5 = 18 +- {d:.0e} (rho={rm:.3e}): salto_rel(zk)={abs(zkp[0]-zkm[0])/abs(zkm[0]):.3e}, "
                  f"salto_rel(vrho)={abs(vrp[0]-vrm[0])/abs(vrm[0]):.3e}")
        r0, s0 = rho_sigma_for_aux5(18.5, omega)
        for h_rel in (1e-4, 1e-5, 1e-6):
            h = h_rel * r0
            pts = np.array([r0 - 2*h, r0 - h, r0, r0 + h, r0 + 2*h])
            zk, vrho, _ = eval_gga(fid, pts, np.full(5, s0))
            E = pts * zk
            fd4 = (E[0] - 8*E[1] + 8*E[3] - E[4]) / (12.0 * h)
            print(f"  Consistencia FD4 en aux5=18.5 (h_rel={h_rel:.0e}): |vrho - FD4|/|FD4| = {abs(vrho[2]-fd4)/abs(fd4):.3e}")
