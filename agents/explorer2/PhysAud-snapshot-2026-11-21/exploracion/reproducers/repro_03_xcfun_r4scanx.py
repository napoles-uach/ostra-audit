#!/usr/bin/env python3
"""Minimal reproducer 03: XCFun 2.1.1 R4SCANX (id=54) uses wrong damping parameters
dp4 = 0.232, da4 = 0.232 instead of published dp4 = 0.802, da4 = 0.178
(Furness et al., J. Chem. Phys. 156, 034109 (2022)).
Proved to 3.7e-15 against Libxc mgga_x_r4scan with xc_func_set_ext_params.
Impact: +0.46 to +0.93 kcal/mol silent error on H2O, N2, CH4 (def2-TZVP).
"""
import ctypes, os
import numpy as np
from pyscf import gto, dft
from pyscf.dft import xcfun, libxc

LIBXC_PATH = os.environ.get("LIBXC_PATH", "/Users/napolesd/Documents/Antigravity/PhysAud/build-libxc/libxc.dylib")
lib = ctypes.CDLL(LIBXC_PATH)
lib.xc_func_alloc.restype = ctypes.c_void_p
lib.xc_func_init.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
lib.xc_func_init.restype = ctypes.c_int
lib.xc_func_end.argtypes = [ctypes.c_void_p]
lib.xc_func_free.argtypes = [ctypes.c_void_p]
lib.xc_functional_get_number.argtypes = [ctypes.c_char_p]
lib.xc_functional_get_number.restype = ctypes.c_int
lib.xc_func_get_info.argtypes = [ctypes.c_void_p]
lib.xc_func_get_info.restype = ctypes.c_void_p
lib.xc_func_info_get_n_ext_params.argtypes = [ctypes.c_void_p]
lib.xc_func_info_get_n_ext_params.restype = ctypes.c_int
lib.xc_func_info_get_ext_params_name.argtypes = [ctypes.c_void_p, ctypes.c_int]
lib.xc_func_info_get_ext_params_name.restype = ctypes.c_char_p
lib.xc_func_info_get_ext_params_default_value.argtypes = [ctypes.c_void_p, ctypes.c_int]
lib.xc_func_info_get_ext_params_default_value.restype = ctypes.c_double
lib.xc_func_set_ext_params.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
lib.xc_mgga_exc_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t] + [ctypes.POINTER(ctypes.c_double)] * 9

rho = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 0.1, 0.1, 0.3], dtype=np.float64)
s   = np.array([0.01, 0.1,  0.2, 0.4, 0.8, 1.5, 3.0, 0.3, 0.6, 0.25], dtype=np.float64)
alp = np.array([0.1,  0.5,  0.8, 1.0, 1.2, 1.8, 3.0, 0.01, 2.5, 0.9], dtype=np.float64)
kf = (3 * np.pi**2 * rho)**(1/3)
sig = (2 * kf * rho * s)**2
tau = sig / (8 * rho) + alp * 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
lapl = np.zeros_like(rho)
pt = np.vstack([rho, np.sqrt(sig), np.zeros_like(rho), np.zeros_like(rho), lapl, tau])

ex_x, vx_x, _, _ = xcfun.eval_xc("R4SCANX", pt, spin=0, deriv=1)
ex_l, vx_l, _, _ = libxc.eval_xc("MGGA_X_R4SCAN", pt, spin=0, deriv=1)

p = lib.xc_func_alloc()
assert lib.xc_func_init(p, lib.xc_functional_get_number(b"mgga_x_r4scan"), 1) == 0
info = lib.xc_func_get_info(p)
n_ext = lib.xc_func_info_get_n_ext_params(info)
params = [lib.xc_func_info_get_ext_params_default_value(info, i) for i in range(n_ext)]
for i in range(n_ext):
    nm = lib.xc_func_info_get_ext_params_name(info, i).decode()
    if nm in ("_dp4", "_da4"):
        params[i] = 0.232
lib.xc_func_set_ext_params(p, (ctypes.c_double * n_ext)(*params))

N = len(rho)
exc_m = np.zeros(N); vr_m = np.zeros(N); vs_m = np.zeros(N); vl_m = np.zeros(N); vt_m = np.zeros(N)
dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
lib.xc_mgga_exc_vxc(p, N, dp(rho), dp(sig), dp(lapl), dp(tau), dp(exc_m), dp(vr_m), dp(vs_m), dp(vl_m), dp(vt_m))
lib.xc_func_end(p); lib.xc_func_free(p)

print("=== 1. Exact 15-digit parameter proof ===")
print(f"  Default Libxc (dp4=0.802, da4=0.178) vs XCFun: max |dexc| = {np.max(np.abs(ex_x - ex_l)):.3e}")
print(f"  Libxc (dp4=0.232, da4=0.232) vs XCFun:")
print(f"    max |dexc|    = {np.max(np.abs(ex_x - exc_m)):.3e}")
print(f"    max |dvrho|   = {np.max(np.abs(vx_x[0] - vr_m)):.3e}")
print(f"    max |dvsigma| = {np.max(np.abs(vx_x[1] - vs_m)):.3e}")
print(f"    max |dvtau|   = {np.max(np.abs(vx_x[3] - vt_m)):.3e}")

print("\n=== 2. Observable SCF impact (r4SCAN, def2-TZVP, grid level 4) ===")
for name, atom in [
    ("H2O", "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"),
    ("N2",  "N 0 0 0; N 0 0 1.0977"),
    ("CH4", "C 0 0 0; H 0.629 0.629 0.629; H -0.629 -0.629 0.629; H -0.629 0.629 -0.629; H 0.629 -0.629 -0.629"),
]:
    mol = gto.M(atom=atom, basis="def2-tzvp", verbose=0)
    mf_l = dft.RKS(mol); mf_l.xc = "MGGA_X_R4SCAN,MGGA_C_R2SCAN"; mf_l.grids.level = 4; mf_l.conv_tol = 1e-10; el = mf_l.kernel()
    mf_x = dft.RKS(mol); mf_x._numint.libxc = dft.xcfun; mf_x.xc = "R4SCANX,R2SCANC"; mf_x.grids.level = 4; mf_x.conv_tol = 1e-10; ex = mf_x.kernel()
    print(f"  {name:4s}: E(Libxc)={el:.8f}, E(XCFun)={ex:.8f}, diff={ex-el:+.6f} Eh ({(ex-el)*627.5095:+.4f} kcal/mol)")
