#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidato #3:
Corte prematuro `pc07_thr = 1/40` en `pc07_fab` (maple/mgga_exc/mgga_k_pc07.mpl:24-27)
reutilizado sin escalar con 1/b en r2SCAN-L (718, 719), SCAN-L (700, 702), revSCAN-L (701, 703) y PC07-opt (634).
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
L.xc_func_set_tau_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_mgga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 9

PI = math.pi
X2S = 1.0 / (2.0 * (6.0 * PI**2) ** (1.0 / 3.0))


def q_for_pc07_z(target_z):
    q = (target_z - 1.0) * 9.0 / 20.0
    for _ in range(40):
        def zfn(qq):
            dlt = 8.0 * qq * qq / 81.0
            return (1.0 + 20.0 * qq / 9.0 + dlt) / math.sqrt(1.0 + dlt * dlt)
        q -= (zfn(q) - target_z) / ((zfn(q + 1e-8) - zfn(q - 1e-8)) / 2e-8)
    return q


def eval_at_z(fid, z_val, rho=0.1):
    lapl = 2.0 * (q_for_pc07_z(z_val) / (X2S**2)) * ((rho / 2.0) ** (5.0 / 3.0))
    p = L.xc_func_alloc()
    assert L.xc_func_init(p, fid, 1) == 0
    L.xc_func_set_dens_threshold(p, 1e-30)
    L.xc_func_set_sigma_threshold(p, 1e-40)
    L.xc_func_set_tau_threshold(p, 1e-40)
    ra = np.array([rho]); sa = np.array([0.0]); la = np.array([lapl]); ta = np.array([0.1])
    zk, vr, vs, vl, vt = [np.zeros(1) for _ in range(5)]
    L.xc_mgga_exc_vxc(p, 1, ra.ctypes.data, sa.ctypes.data, la.ctypes.data, ta.ctypes.data,
                      zk.ctypes.data, vr.ctypes.data, vs.ctypes.data, vl.ctypes.data, vt.ctypes.data)
    L.xc_func_end(p); L.xc_func_free(p)
    return zk[0], vl[0]


if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")
    for a, b, tag in [(0.5389, 3.0, "PC07 (b=3)"), (1.784720, 0.258304, "PC07-opt / r2SCAN-L (b=0.258304)")]:
        val = ((1.0 + math.exp(-40.0 / 39.0)) ** b) / ((1.0 + math.exp(40.0 - 40.0 / 39.0)) ** b)
        print(f"pc07_fab0(z_lo = {a/40:.6f}) en {tag}: {val:.3e}")
    print()

    for label, fid, zc in [
        ("CONTROL mgga_k_pc07 (id=543, b=3.0)", 543, 0.5389 / 40.0),
        ("mgga_k_pc07_opt     (id=634, b=0.258)", 634, 1.784720 / 40.0),
        ("mgga_x_r2scanl      (id=718, b=0.258)", 718, 1.784720 / 40.0),
        ("mgga_c_r2scanl      (id=719, b=0.258)", 719, 1.784720 / 40.0),
        ("mgga_x_scanl        (id=700, b=0.258)", 700, 1.784720 / 40.0),
        ("mgga_c_scanl        (id=702, b=0.258)", 702, 1.784720 / 40.0),
    ]:
        print(f"=== {label} en z_lo = {zc:.6f} ===")
        for d in (1e-9, 1e-13):
            zkm, vlm = eval_at_z(fid, zc - d)
            zkp, vlp = eval_at_z(fid, zc + d)
            rj_zk = abs(zkp - zkm) / max(abs(zkm), abs(zkp), 1e-30)
            rj_vl = abs(vlp - vlm) / max(abs(vlm), abs(vlp), 1e-30) if max(abs(vlm), abs(vlp)) > 1e-25 else 0.0
            print(f"  d={d:.0e}: salto_rel(zk)={rj_zk:.3e}, vlapl(-)={vlm:.4e}, vlapl(+)={vlp:.4e}, salto_rel(vlapl)={rj_vl:.3e}")
