#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidato #2:
Discontinuidades en gga_x_wpbeh / familia HSE y LC-wPBE:
  (a) s_scaling_2 en s = 1.0  (maple/scaling.mpl:17)
  (b) s_scaling_2 en s = 15.0 (maple/scaling.mpl:18)
  (c) eb1 / term1 en w = omega/k_F = 14.0 (maple/gga_exc/gga_x_wpbeh.mpl:66,294) + pico puntual en w == 14.0
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


def sigma_from_s(rho, s):
    return 4.0 * ((s / X2S) * (rho / 2.0) ** (4.0 / 3.0)) ** 2


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


if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")
    for name, fid, omega in [("HSE06 (id=428)", 428, 0.11), ("LC-wPBE (id=478)", 478, 0.40), ("gga_x_wpbeh (id=524)", 524, 0.0)]:
        print(f"=== {name} ===")
        for sc, rho0 in [(1.0, 0.08), (15.0, 1e-4)]:
            for d in (1e-9, 1e-14):
                zkm, _, vsm = eval_gga(fid, [rho0], [sigma_from_s(rho0, sc - d)])
                zkp, _, vsp = eval_gga(fid, [rho0], [sigma_from_s(rho0, sc + d)])
                print(f"  s_c={sc:4.1f} +- {d:.0e} (rho={rho0:.2g}): salto_rel(zk)={abs(zkp[0]-zkm[0])/abs(zkm[0]):.3e}, "
                      f"salto_rel(vsigma)={abs(vsp[0]-vsm[0])/max(abs(vsm[0]),abs(vsp[0])):.3e}")
        if omega > 0.0:
            for d in (1e-9, 1e-13):
                rm = ((omega / (14.0 - d)) ** 3) / (3.0 * PI**2)
                rp = ((omega / (14.0 + d)) ** 3) / (3.0 * PI**2)
                zkm, vrm, _ = eval_gga(fid, [rm], [sigma_from_s(rm, 0.5)])
                zkp, vrp, _ = eval_gga(fid, [rp], [sigma_from_s(rp, 0.5)])
                print(f"  w_c=14.0 +- {d:.0e} (s=0.5): salto_rel(zk)={abs(zkp[0]-zkm[0])/abs(zkm[0]):.3e}, "
                      f"salto_rel(vrho)={abs(vrp[0]-vrm[0])/abs(vrm[0]):.3e}")
            r = ((omega / 14.0) ** 3) / (3.0 * PI**2)
            for _ in range(8):
                r = float(np.nextafter(r, -np.inf))
            pts = []
            for _ in range(17):
                pts.append(r); r = float(np.nextafter(r, np.inf))
            zk_u, _, _ = eval_gga(fid, pts, [sigma_from_s(x, 0.5) for x in pts])
            spike = (np.max(np.abs(zk_u)) - np.min(np.abs(zk_u))) / np.min(np.abs(zk_u))
            print(f"  Pico puntual en barrido de 17 ULPs alrededor de w == 14.0: salto_rel(zk) = {spike:.3e}")
        print()
