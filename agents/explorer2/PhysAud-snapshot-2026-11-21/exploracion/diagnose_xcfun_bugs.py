import os, ctypes
import numpy as np
from pyscf.dft import xcfun, libxc

LIBXC_PATH = "/Users/napolesd/Documents/Antigravity/PhysAud/build-libxc/libxc.dylib"
lib = ctypes.CDLL(LIBXC_PATH)

def make_pt(rho, s, alpha):
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
    tau_w = sig / (8 * rho)
    tau = tau_w + alpha * tau_unif
    return np.array([[rho], [np.sqrt(sig)], [0.0], [0.0], [0.0], [tau]])

# =====================================================================
# 1. EXACT PROOF: XCFun R4SCANX == Libxc MGGA_X_R4SCAN(dp4=0.232, da4=0.232)
# =====================================================================
print("=== 1. XCFun R4SCANX vs Libxc mgga_x_r4scan(ext_params) ===")
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
lib.xc_mgga_exc_vxc.argtypes = [
    ctypes.c_void_p, ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]

fid = lib.xc_functional_get_number(b"mgga_x_r4scan")
p = lib.xc_func_alloc()
assert lib.xc_func_init(p, fid, 1) == 0
info = lib.xc_func_get_info(p)
n_ext = lib.xc_func_info_get_n_ext_params(info)
defaults = []
names = []
for i in range(n_ext):
    nm = lib.xc_func_info_get_ext_params_name(info, i).decode()
    dv = lib.xc_func_info_get_ext_params_default_value(info, i)
    names.append(nm)
    defaults.append(dv)
    print(f"  ext_param[{i}] {nm} = {dv}")

rho_pts = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 0.1, 0.1, 0.3], dtype=np.float64)
s_pts   = np.array([0.01, 0.1,  0.2, 0.4, 0.8, 1.5, 3.0, 0.3, 0.6, 0.25], dtype=np.float64)
a_pts   = np.array([0.1,  0.5,  0.8, 1.0, 1.2, 1.8, 3.0, 0.01, 2.5, 0.9], dtype=np.float64)
kf = (3 * np.pi**2 * rho_pts)**(1/3)
sig = (2 * kf * rho_pts * s_pts)**2
tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho_pts**(5/3)
tau_w = sig / (8 * rho_pts)
tau = tau_w + a_pts * tau_unif
lapl = np.zeros_like(rho_pts)
pt_pyscf = np.vstack([rho_pts, np.sqrt(sig), np.zeros_like(rho_pts), np.zeros_like(rho_pts), lapl, tau])

ex_xcfun, vx_xcfun, _, _ = xcfun.eval_xc("R4SCANX", pt_pyscf, spin=0, deriv=1)
ex_libxc_def, vx_libxc_def, _, _ = libxc.eval_xc("MGGA_X_R4SCAN", pt_pyscf, spin=0, deriv=1)
print(f"  Default Libxc vs XCFun max |dexc| = {np.max(np.abs(ex_xcfun - ex_libxc_def)):.3e}")

# Try setting dp4=0.232, da4=0.232 (check names!)
params = list(defaults)
for i, nm in enumerate(names):
    if "dp4" in nm.lower() or "da4" in nm.lower() or "d_p4" in nm.lower() or "d_a4" in nm.lower():
        params[i] = 0.232
arr = (ctypes.c_double * n_ext)(*params)
lib.xc_func_set_ext_params(p, arr)

N = len(rho_pts)
exc_m = np.zeros(N); vr_m = np.zeros(N); vs_m = np.zeros(N); vl_m = np.zeros(N); vt_m = np.zeros(N)
dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
lib.xc_mgga_exc_vxc(p, N, dp(rho_pts), dp(sig), dp(lapl), dp(tau),
                    dp(exc_m), dp(vr_m), dp(vs_m), dp(vl_m), dp(vt_m))
print("  Libxc(dp4=0.232, da4=0.232) vs XCFun R4SCANX:")
print(f"    max |dexc|    = {np.max(np.abs(ex_xcfun - exc_m)):.3e}")
print(f"    max |dvrho|   = {np.max(np.abs(vx_xcfun[0] - vr_m)):.3e}")
print(f"    max |dvsigma| = {np.max(np.abs(vx_xcfun[1] - vs_m)):.3e}")
print(f"    max |dvtau|   = {np.max(np.abs(vx_xcfun[3] - vt_m)):.3e}")
lib.xc_func_end(p); lib.xc_func_free(p)

# =====================================================================
# 2. XCFun REVTPSSX vs Libxc MGGA_X_REVTPSS
# =====================================================================
print("\n=== 2. XCFun REVTPSSX vs Libxc MGGA_X_REVTPSS ===")
print(f"{'rho':>6} {'s':>6} {'alpha':>6} | {'ex_xcfun':>13} {'ex_libxc':>13} {'ratio':>10} | {'Fx_xcfun':>10} {'Fx_libxc':>10}")
for rho, s, alpha in [
    (0.1, 0.0, 1.0),
    (0.1, 0.0, 0.5),
    (0.1, 0.0, 2.0),
    (0.1, 0.5, 1.0),
    (0.1, 0.5, 0.0),
    (0.1, 0.5, 2.0),
    (0.1, 1.0, 1.0),
    (0.1, 2.0, 1.0),
    (1.0, 0.5, 1.0),
    (1.0, 0.0, 1.0),
]:
    pt = make_pt(rho, s, alpha)
    e_xc = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0]
    e_c  = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ex_x = e_xc - e_c
    ex_l = libxc.eval_xc("MGGA_X_REVTPSS", pt, spin=0, deriv=0)[0][0]
    ex_lda = -0.75 * (3 / np.pi)**(1/3) * rho**(1/3)
    print(f"{rho:6.2f} {s:6.2f} {alpha:6.2f} | {ex_x:13.8f} {ex_l:13.8f} {ex_x/ex_l:10.6f} | {ex_x/ex_lda:10.6f} {ex_l/ex_lda:10.6f}")

# =====================================================================
# 3. XCFun REVTPSSC vs Libxc MGGA_C_REVTPSS
# =====================================================================
print("\n=== 3. XCFun REVTPSSC vs Libxc MGGA_C_REVTPSS ===")
print(f"{'rho':>6} {'s':>6} {'alpha':>6} | {'ec_xcfun':>13} {'ec_libxc':>13} {'diff':>12} | {'ec_tpss_l':>13}")
for rho, s, alpha in [
    (0.1, 0.0, 1.0),
    (0.1, 0.0, 0.0),
    (0.1, 0.0, 2.0),
    (0.1, 0.5, 1.0),
    (0.1, 0.5, 0.0),
    (0.1, 0.5, 1.2),
    (0.1, 1.0, 1.0),
    (1.0, 0.5, 1.0),
    (10.0, 0.5, 1.0),
    (0.001, 0.5, 1.0),
]:
    pt = make_pt(rho, s, alpha)
    ec_x = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ec_l = libxc.eval_xc("MGGA_C_REVTPSS", pt, spin=0, deriv=0)[0][0]
    ec_t = libxc.eval_xc("MGGA_C_TPSS", pt, spin=0, deriv=0)[0][0]
    print(f"{rho:6.3f} {s:6.2f} {alpha:6.2f} | {ec_x:13.8f} {ec_l:13.8f} {ec_x-ec_l:+12.4e} | {ec_t:13.8f}")

# =====================================================================
# 4. XCFun SCANX & SCANC edge cases (s -> 0 and alpha -> 1)
# =====================================================================
print("\n=== 4. XCFun SCANX & SCANC edge cases ===")
for s in [0.0, 1e-12, 1e-8, 1e-6, 1e-4, 1e-2]:
    for alpha in [0.5, 1.0, 1.0 + 1e-10, 1.2]:
        pt = make_pt(0.1, s, alpha)
        ex, vx, _, _ = xcfun.eval_xc("SCANX", pt, spin=0, deriv=1)
        ec, vc, _, _ = xcfun.eval_xc("SCANC", pt, spin=0, deriv=1)
        exl, vxl, _, _ = libxc.eval_xc("MGGA_X_SCAN", pt, spin=0, deriv=1)
        ecl, vcl, _, _ = libxc.eval_xc("MGGA_C_SCAN", pt, spin=0, deriv=1)
        bad_x = not np.isfinite(ex[0]) or not np.isfinite(vx[1][0]) or abs(vx[1][0]-vxl[1][0])/(abs(vxl[1][0])+1e-12) > 1e-4
        bad_c = not np.isfinite(ec[0]) or not np.isfinite(vc[1][0]) or abs(vc[1][0]-vcl[1][0])/(abs(vcl[1][0])+1e-12) > 1e-4
        if bad_x or bad_c:
            print(f"  s={s:7.1e} a={alpha:<12.10g} | SCANX: ex={ex[0]:+.4e} vs={vx[1][0]:+.4e} (L:{vxl[1][0]:+.4e}) | SCANC: ec={ec[0]:+.4e} vs={vc[1][0]:+.4e} (L:{vcl[1][0]:+.4e})")
