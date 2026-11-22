#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidatos secundarios #6, #7 y #8:
  - Candidato #6: mgga_c_lak (id=328, maple/mgga_exc/mgga_c_lak.mpl:61,64) — m_max(rs*a, eps) corrompe alpha < 0
  - Candidato #7: lda_c_pk09 (id=554, maple/lda_exc/lda_c_pk09.mpl:25,44) — desborde y NaN cerca del polo rs ~ 22.741825
  - Candidato #8: gga_x_lb   (id=160, maple/gga_vxc/gga_x_lb.mpl:31) — salto 0.321% en vrho en xs = 300
                  mgga_x_br89_explicit (id=586, maple/mgga_exc/mgga_x_br89.mpl:36) — banda muerta |Q| < 5e-13
"""
import ctypes, math, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
LIB_PATH = os.environ.get("LIBXC_PATH", str(ROOT / "build-libxc" / "libxc.dylib"))
L = ctypes.CDLL(LIB_PATH)

L.xc_version_string.restype = ctypes.c_char_p
L.xc_functional_get_number.argtypes = [ctypes.c_char_p]
L.xc_func_alloc.restype = ctypes.c_void_p
L.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
L.xc_func_end.argtypes = [ctypes.c_void_p]
L.xc_func_free.argtypes = [ctypes.c_void_p]
L.xc_func_set_dens_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
L.xc_func_set_fhc_enforcement.argtypes = [ctypes.c_void_p, ctypes.c_int]
L.xc_lda_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 3
L.xc_gga_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 4
L.xc_mgga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 9

PI = math.pi
RS_FACTOR = (3.0 / (4.0 * PI)) ** (1.0 / 3.0)
K_C = 0.3 * (6.0 * PI**2) ** (2.0 / 3.0)

if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")

    # 1. LAK_C (mgga_c_lak) en alpha = +-1e-9 y alpha = -0.2 (con FHC desactivado)
    fid_lak = L.xc_functional_get_number(b"mgga_c_lak")
    p = L.xc_func_alloc(); L.xc_func_init(p, fid_lak, 1); L.xc_func_set_fhc_enforcement(p, 0)
    print("=== Candidato #6: mgga_c_lak (maple/mgga_exc/mgga_c_lak.mpl:61,64) ===")
    rho0, sig0 = 0.1, 0.02
    tau_w = sig0 / (8.0 * rho0)
    tau_u = 2.0 * K_C * ((rho0 / 2.0) ** (5.0 / 3.0))
    for a in (-0.2, -1e-9, -1e-14, +1e-14, +1e-9, +0.2):
        ra = np.array([rho0]); sa = np.array([sig0]); la = np.array([0.0]); ta = np.array([tau_w + a * tau_u])
        zk, vr, vs, vl, vt = [np.zeros(1) for _ in range(5)]
        L.xc_mgga_exc_vxc(p, 1, ra.ctypes.data, sa.ctypes.data, la.ctypes.data, ta.ctypes.data,
                          zk.ctypes.data, vr.ctypes.data, vs.ctypes.data, vl.ctypes.data, vt.ctypes.data)
        print(f"  alpha={a:+.1e}: zk={zk[0]:.15e}, vtau={vt[0]:+.8e}")
    L.xc_func_end(p); L.xc_func_free(p)

    # 2. PK09 (lda_c_pk09) cerca del polo real rs = 22.739412132480812 (zeta = 0.4)
    fid_pk = L.xc_functional_get_number(b"lda_c_pk09")
    p = L.xc_func_alloc(); L.xc_func_init(p, fid_pk, 2)
    print("\n=== Candidato #7: lda_c_pk09 (maple/lda_exc/lda_c_pk09.mpl:25,44; zeta=0.4) ===")
    rs_pole = 22.739412132480812
    for drs in (-1e-2, -1e-4, -1e-6, 0.0, +1e-6, +1e-4, +1e-2):
        rs = rs_pole + drs
        rtot = (RS_FACTOR / rs) ** 3
        rho = np.array([[0.7 * rtot, 0.3 * rtot]])
        zk = np.zeros(1); vr = np.zeros((1, 2))
        L.xc_lda_exc_vxc(p, 1, rho.ctypes.data, zk.ctypes.data, vr.ctypes.data)
        print(f"  rs = rs_polo {drs:+.0e} ({rs:.8f}): zk={zk[0]:12.5e}, vrho=[{vr[0,0]:12.4e}, {vr[0,1]:12.4e}]")
    L.xc_func_end(p); L.xc_func_free(p)

    # 3. LB94 (gga_x_lb) en xs = 300 +- delta
    fid_lb = L.xc_functional_get_number(b"gga_x_lb")
    p = L.xc_func_alloc(); L.xc_func_init(p, fid_lb, 1); L.xc_func_set_dens_threshold(p, 1e-30)
    print("\n=== Candidato #8a: gga_x_lb (maple/gga_vxc/gga_x_lb.mpl:31) en xs = 300 ===")
    r_lb = 1e-4
    for d in (1e-9, 1e-13):
        vrs = []
        for sgn in (-1.0, +1.0):
            xs = 300.0 * (1.0 + sgn * d)
            sig = 4.0 * (xs * (r_lb / 2.0) ** (4.0 / 3.0)) ** 2
            ra = np.array([r_lb]); sa = np.array([sig]); vr = np.zeros(1); vs = np.zeros(1)
            L.xc_gga_vxc(p, 1, ra.ctypes.data, sa.ctypes.data, vr.ctypes.data, vs.ctypes.data)
            vrs.append(vr[0])
        print(f"  xs = 300*(1 +- {d:.0e}): vrho(-)={vrs[0]:.12e}, vrho(+)={vrs[1]:.12e}, salto_rel={abs(vrs[1]-vrs[0])/abs(vrs[0]):.3e}")
    L.xc_func_end(p); L.xc_func_free(p)

    # 4. BR89 / BR89_explicit en la banda muerta |Q| < 5e-13 (maple/mgga_exc/mgga_x_br89.mpl:36)
    fid_br = L.xc_functional_get_number(b"mgga_x_br89_explicit")
    p = L.xc_func_alloc(); L.xc_func_init(p, fid_br, 1)
    print("\n=== Candidato #8b: mgga_x_br89_explicit (maple/mgga_exc/mgga_x_br89.mpl:36) en Q = 0 y |Q| = 5e-13 ===")
    rho_b = 0.1; xs_b = 0.5; sig_b = 4.0 * (xs_b * (rho_b / 2.0) ** (4.0 / 3.0)) ** 2
    t_s = 1.2 * K_C; tau_b = 2.0 * t_s * ((rho_b / 2.0) ** (5.0 / 3.0))
    for q_val in (-1e-12, -4.9e-13, -1e-15, +1e-15, +4.9e-13, +1e-12):
        u_s = 6.0 * q_val + 4.0 * 0.8 * t_s - 0.5 * 0.8 * (xs_b**2)
        lap_b = 2.0 * u_s * ((rho_b / 2.0) ** (5.0 / 3.0))
        ra = np.array([rho_b]); sa = np.array([sig_b]); la = np.array([lap_b]); ta = np.array([tau_b])
        zk, vr, vs, vl, vt = [np.zeros(1) for _ in range(5)]
        L.xc_mgga_exc_vxc(p, 1, ra.ctypes.data, sa.ctypes.data, la.ctypes.data, ta.ctypes.data,
                          zk.ctypes.data, vr.ctypes.data, vs.ctypes.data, vl.ctypes.data, vt.ctypes.data)
        print(f"  Q={q_val:+.2e}: zk={zk[0]:.15e}, vlapl={vl[0]:+.8e}")
    L.xc_func_end(p); L.xc_func_free(p)
