#!/usr/bin/env python3
"""Reproductor minimo autocontenido — Candidato #5:
Violacion de contrato de API en `xc_func_set_fhc_enforcement` (src/functionals.c:511-518):
no propaga el apagado/encendido de XC_FLAGS_ENFORCE_FHC a `p->func_aux[ii]`, convirtiendose
en un no-op silencioso en los 38 meta-GGAs mixtos/hibridos (TPSSh, r2SCAN0, SCAN0, SCAN-rVV10, PW6B95...).
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
L.xc_func_set_fhc_enforcement.argtypes = [ctypes.c_void_p, ctypes.c_int]
L.xc_mgga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.c_void_p] * 9


def eval_mgga(fid, fhc_on):
    # Punto con tau = 0.02 < tau_W = sigma/(8*rho) = 0.10 / 0.96 = 0.10417 (alpha < 0, tipico en esferas PAW)
    rho = np.array([0.12]); sig = np.array([0.10]); lap = np.array([0.0]); tau = np.array([0.02])
    zk, vr, vs, vl, vt = [np.zeros(1) for _ in range(5)]
    p = L.xc_func_alloc()
    assert L.xc_func_init(p, fid, 1) == 0
    L.xc_func_set_fhc_enforcement(p, 1 if fhc_on else 0)
    L.xc_mgga_exc_vxc(p, 1, rho.ctypes.data, sig.ctypes.data, lap.ctypes.data, tau.ctypes.data,
                      zk.ctypes.data, vr.ctypes.data, vs.ctypes.data, vl.ctypes.data, vt.ctypes.data)
    L.xc_func_end(p); L.xc_func_free(p)
    return zk[0], vt[0]


if __name__ == "__main__":
    print(f"libxc version={L.xc_version_string().decode()} ({LIB_PATH})\n")
    pairs = [
        ("mgga_x_r2scan (puro, id=497)", 497),
        ("hyb_mgga_xc_r2scan0 (mixto, id=660)", 660),
        ("mgga_x_tpss (puro, id=202)", 202),
        ("hyb_mgga_xc_tpssh (mixto, id=457)", 457),
        ("mgga_x_scan (puro, id=263)", 263),
        ("hyb_mgga_x_scan0 (mixto, id=264)", 264),
    ]
    for label, fid in pairs:
        zk1, vt1 = eval_mgga(fid, True)
        zk0, vt0 = eval_mgga(fid, False)
        print(f"{label:38s}: |dzk|={abs(zk1-zk0):.3e}, |dvtau|={abs(vt1-vt0):.3e}  (vtau FHC=1: {vt1:+.5e}, FHC=0: {vt0:+.5e})")
