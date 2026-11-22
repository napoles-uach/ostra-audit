#!/usr/bin/env python3
"""Comprehensive boundary continuity tests for branch switches across libxc (maple/ and src/).

Prioritized by real-world usage:
  Tier 1: PBE, PBE0, HSE06 (and LC-wPBE/wpbeh), SCAN, r2SCAN, TPSS, B3LYP, wB97X
  Tier 2: Shared helpers in maple/util.mpl, maple/attenuation.mpl, maple/scaling.mpl, src/expint_e1.c
  Tier 3: Functional-specific piecewise / series / clamp switches across LDA, GGA, MGGA

For each branch switch at x_c, evaluates functional outputs (zk and first derivatives
vrho, vsigma, vtau, vlapl) at x_c +- delta_wide (~1e-9) and x_c +- delta_tight (~1e-14),
recording relative and absolute jumps. A true C0/C1 discontinuity exhibits a jump that
does NOT shrink by 10^5 between delta_wide and delta_tight.
"""
import csv
import ctypes
import math
import os
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
if "LIBXC_PATH" not in os.environ:
    os.environ["LIBXC_PATH"] = str(ROOT / "build-libxc" / "libxc.dylib")

import lxc  # noqa: E402

# Bind additional threshold / FHC controls and vxc/fxc entrypoints on lxc.L
lxc.L.xc_func_set_sigma_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
lxc.L.xc_func_set_zeta_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
lxc.L.xc_func_set_tau_threshold.argtypes = [ctypes.c_void_p, ctypes.c_double]
lxc.L.xc_func_set_fhc_enforcement.argtypes = [ctypes.c_void_p, ctypes.c_int]
lxc.L.xc_lda_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p]
lxc.L.xc_gga_vxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
lxc.L.xc_mgga_vxc.argtypes = [
    ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
    ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
]
lxc.L.xc_lda_fxc.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p]


def eval_func(f, x, include_fxc=False):
    if f.flags & lxc.FLAG_EXC:
        out = f.exc_vxc(x)
    else:
        rho = np.ascontiguousarray(x["rho"])
        n = rho.shape[0]
        ns = f.nspin
        nsig = 1 if ns == 1 else 3
        zk = np.zeros(n)
        vrho = np.zeros((n, ns))
        out = {"zk": zk, "vrho": vrho}
        if f.family == "LDA":
            lxc.L.xc_lda_vxc(ctypes.c_void_p(f.p), ctypes.c_size_t(n), lxc._p(rho), lxc._p(vrho))
        elif f.family == "GGA":
            sig = np.ascontiguousarray(x["sigma"])
            vs = np.zeros((n, nsig))
            lxc.L.xc_gga_vxc(ctypes.c_void_p(f.p), ctypes.c_size_t(n), lxc._p(rho), lxc._p(sig), lxc._p(vrho), lxc._p(vs))
            out["vsigma"] = vs
        elif f.family == "MGGA":
            sig = np.ascontiguousarray(x["sigma"])
            lap = np.ascontiguousarray(x["lapl"])
            tau = np.ascontiguousarray(x["tau"])
            vs = np.zeros((n, nsig))
            vl = np.zeros((n, ns))
            vt = np.zeros((n, ns))
            lxc.L.xc_mgga_vxc(
                ctypes.c_void_p(f.p), ctypes.c_size_t(n), lxc._p(rho), lxc._p(sig), lxc._p(lap), lxc._p(tau),
                lxc._p(vrho), lxc._p(vs), lxc._p(vl), lxc._p(vt),
            )
            out.update(vsigma=vs, vlapl=vl, vtau=vt)
    if include_fxc and f.family == "LDA" and (f.flags & lxc.FLAG_FXC):
        rho = np.ascontiguousarray(x["rho"])
        n = rho.shape[0]
        n2 = 1 if f.nspin == 1 else 3
        v2rho2 = np.zeros((n, n2))
        lxc.L.xc_lda_fxc(ctypes.c_void_p(f.p), ctypes.c_size_t(n), lxc._p(rho), lxc._p(v2rho2))
        out["v2rho2"] = v2rho2
    return out

# Exact constants matching maple/util.mpl and src/util.h
PI = math.pi
EPS = float(np.finfo(np.float64).eps)
RS_FACTOR = (3.0 / (4.0 * PI)) ** (1.0 / 3.0)
X2S = 1.0 / (2.0 * (6.0 * PI**2) ** (1.0 / 3.0))
XT2S = 1.0 / (2.0 * (3.0 * PI**2) ** (1.0 / 3.0))
K_FACTOR_C = 0.3 * (6.0 * PI**2) ** (2.0 / 3.0)
X_FACTOR_C = (3.0 / 8.0) * (3.0 / PI) ** (1.0 / 3.0) * (4.0 ** (2.0 / 3.0))
M_C = 137.0359996287515


def rho_from_rs(rs, dim=3):
    if dim == 2:
        return (1.0 / math.sqrt(PI) / rs) ** 2
    return (RS_FACTOR / rs) ** 3


def unpol_sigma_from_s(rho, s):
    """Total sigma for unpolarized 3D density such that gga_s(xs0) = XT2S * xt = s."""
    # For z=0: n_spin = rho/2, sigma_spin = sigma/4, xs0 = sqrt(sigma/4)/(rho/2)^(4/3)
    # gga_s(xs0) = X2S * xs0 = s  => sigma = 4 * ((s / X2S) * (rho / 2)^(4/3))^2
    return 4.0 * ((s / X2S) * (rho / 2.0) ** (4.0 / 3.0)) ** 2


def unpol_sigma_from_xs(rho, xs):
    """Total sigma for unpolarized 3D density given Libxc reduced gradient xs = |grad rho_s|/rho_s^(4/3)."""
    return 4.0 * (xs * (rho / 2.0) ** (4.0 / 3.0)) ** 2


def unpol_tau_from_alpha(rho, sigma, alpha, eta=0.0):
    """Total tau for unpolarized 3D MGGA given alpha = (t - p/(8*X2S^2))/(K_FACTOR_C + eta*p/(8*X2S^2))."""
    # tau_s = tau / 2, rho_s = rho / 2, sigma_s = sigma / 4
    # t_s = tau_s / rho_s^(5/3), xs^2 = sigma_s / rho_s^(8/3)
    # p/(8*X2S^2) = xs^2 / 8 => tau_W_s = sigma_s / (8 * rho_s)
    rho_s = rho / 2.0
    sig_s = sigma / 4.0
    xs2 = sig_s / (rho_s ** (8.0 / 3.0))
    t_s = xs2 / 8.0 + alpha * (K_FACTOR_C + eta * xs2 / 8.0)
    return 2.0 * t_s * (rho_s ** (5.0 / 3.0))


def unpol_lapl_from_u(rho, u):
    """Total lapl for unpolarized 3D MGGA given reduced u_s = lapl_s / rho_s^(5/3)."""
    rho_s = rho / 2.0
    return 2.0 * u * (rho_s ** (5.0 / 3.0))


def rel_jump(a, b, floor=1e-30):
    denom = max(abs(a), abs(b), floor)
    return float(abs(b - a) / denom)


def eval_pair(fid, nspin, make_pt, xc, relax_thresholds=True, disable_fhc=False, exact_deltas=None, include_fxc=False):
    f = lxc.Func(fid, nspin)
    if relax_thresholds:
        lxc.L.xc_func_set_dens_threshold(f.p, 1e-30)
        lxc.L.xc_func_set_sigma_threshold(f.p, 1e-40)
        lxc.L.xc_func_set_tau_threshold(f.p, -1e6 if disable_fhc else 1e-40)
    if disable_fhc:
        lxc.L.xc_func_set_fhc_enforcement(f.p, 0)

    if exact_deltas is not None:
        d_wide, d_tight = exact_deltas
    elif 1e-3 <= abs(xc) <= 20.0:
        d_wide, d_tight = 1e-9, 1e-13
    elif xc == 0.0:
        d_wide, d_tight = 1e-9, 1e-14
    else:
        d_wide = abs(xc) * 1e-9
        ulp = abs(np.nextafter(xc, np.inf) - xc)
        d_tight = max(abs(xc) * 1e-13, 8.0 * ulp)

    res = {}
    for tag, d in [("wide", d_wide), ("tight", d_tight)]:
        xm = xc - d
        xp = xc + d
        if xm == xc:
            xm = float(np.nextafter(xc, -np.inf))
        if xp == xc:
            xp = float(np.nextafter(xc, np.inf))
        om = eval_func(f, make_pt(xm), include_fxc=include_fxc)
        op = eval_func(f, make_pt(xp), include_fxc=include_fxc)

        zk_m, zk_p = float(om["zk"][0]), float(op["zk"][0])
        zk_rj = rel_jump(zk_m, zk_p) if (f.flags & lxc.FLAG_EXC) else 0.0
        zk_aj = abs(zk_p - zk_m) if (f.flags & lxc.FLAG_EXC) else 0.0

        # Maximum relative jump across all active derivative channels
        max_der_rj = 0.0
        max_der_aj = 0.0
        worst_chan = "vrho[0]"
        worst_m, worst_p = 0.0, 0.0
        chans = ("vrho", "vsigma", "vtau", "vlapl", "v2rho2") if include_fxc else ("vrho", "vsigma", "vtau", "vlapl")
        for k in chans:
            if k not in om:
                continue
            am = om[k][0]
            ap = op[k][0]
            for j in range(len(am)):
                # Avoid false alarms on channels that are identically ~0 (< 1e-25) on both sides
                if max(abs(am[j]), abs(ap[j])) < 1e-25:
                    continue
                rj = rel_jump(am[j], ap[j])
                aj = abs(ap[j] - am[j])
                if rj >= max_der_rj:
                    max_der_rj = rj
                    max_der_aj = aj
                    worst_chan = f"{k}[{j}]"
                    worst_m, worst_p = float(am[j]), float(ap[j])

        res[tag] = {
            "delta": d,
            "zk_m": zk_m, "zk_p": zk_p, "zk_rel_jump": zk_rj, "zk_abs_jump": zk_aj,
            "der_chan": worst_chan, "der_m": worst_m, "der_p": worst_p,
            "der_rel_jump": max_der_rj, "der_abs_jump": max_der_aj,
        }
    f.close()
    return res


def build_test_cases():
    cases = []

    # =========================================================================
    # TIER 1: HIGH-USAGE PRIORITY FUNCTIONALS
    # (PBE, PBE0, HSE06, SCAN, r2SCAN, TPSS, B3LYP, wB97X)
    # =========================================================================

    # 1. PBE / PBE0: gga_c_pbe (130) & hyb_gga_xc_pbeh (406) -> m_saturate(y=1) in fH (reached via f_pbe_grad_p / mgga_c_m08;
    #    note gga_c_pbe itself uses pbe_H_g = 1/(1+A*f1) without m_saturate, while mgga_c_m08 (233) uses f_pbe_grad_p!)
    #    Let's test BOTH gga_c_pbe around A*f1 = 1 AND mgga_c_m08 around A*f1 = 1!
    def _pbe_A(rs):
        # At z=0, mphi=1, f_pw(rs=1, 0)
        f_pw = lxc.Func(12, 1)  # lda_c_pw
        zk = float(f_pw.exc_vxc({"rho": np.array([[rho_from_rs(rs)]])})["zk"][0])
        f_pw.close()
        gamma = (1.0 - math.log(2.0)) / PI**2
        beta = 0.06672455060314922
        return beta / (gamma * math.expm1(-zk / gamma))

    A_rs1 = _pbe_A(1.0)
    # Solve A*(t2 + A*t2^2) = y => A^2 t2^2 + A t2 - y = 0 => A*t2 = (-1 + sqrt(1 + 4*y))/2
    def _pt_pbe_msat(y, rs=1.0, mgga=False):
        rho = rho_from_rs(rs)
        At2 = (-1.0 + math.sqrt(1.0 + 4.0 * y)) / 2.0
        t2 = At2 / A_rs1
        # tt2 = xt2 / (16 * 2^(2/3) * rs) => xt2 = t2 * 16 * 2^(2/3) * rs
        xt2 = t2 * 16.0 * (2.0 ** (2.0 / 3.0)) * rs
        sigma = xt2 * (rho ** (8.0 / 3.0))
        d = {"rho": np.array([[rho]]), "sigma": np.array([[sigma]])}
        if mgga:
            d["lapl"] = np.array([[0.0]])
            d["tau"] = np.array([[unpol_tau_from_alpha(rho, sigma, 1.2)]])
        return d

    cases.append(("1_alta", "PBE0 (hyb_gga_xc_pbeh)", 406, 1, "maple/util.mpl:163", "m_saturate (control en PBE0)", "y = A*(t^2+A*t^4)", 1.0,
                  lambda y: _pt_pbe_msat(y, mgga=False), True, False, None))
    cases.append(("2_infra", "M08-HX_C (mgga_c_m08 via gga_c_pbe)", 235, 1, "maple/gga_exc/gga_c_pbe.mpl:73", "m_saturate(A*f1)", "y = A*(t^2+A*t^4)", 1.0,
                  lambda y: _pt_pbe_msat(y, mgga=True), True, False, None))

    # 2. HSE06 (428), LC-wPBE (478), gga_x_wpbeh (524): ALL 6 BRANCH SWITCHES!
    rho_hse = 0.08  # typical valence density
    # 2a. s_scaling_2 at s = 1 (maple/scaling.mpl:17)
    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/scaling.mpl:17", "s_scaling_2 (corte inferior)", "s", 1.0,
                  lambda s: {"rho": np.array([[rho_hse]]), "sigma": np.array([[unpol_sigma_from_s(rho_hse, s)]])}, True, False, None))
    cases.append(("1_alta", "LC-wPBE (hyb_gga_xc_lc_wpbe)", 478, 1, "maple/scaling.mpl:17", "s_scaling_2 (corte inferior)", "s", 1.0,
                  lambda s: {"rho": np.array([[rho_hse]]), "sigma": np.array([[unpol_sigma_from_s(rho_hse, s)]])}, True, False, None))
    cases.append(("1_alta", "wPBEh_X (gga_x_wpbeh, w=0)", 524, 1, "maple/scaling.mpl:17", "s_scaling_2 (corte inferior)", "s", 1.0,
                  lambda s: {"rho": np.array([[rho_hse]]), "sigma": np.array([[unpol_sigma_from_s(rho_hse, s)]])}, True, False, None))

    # 2b. s_scaling_2 at s = 15 (maple/scaling.mpl:18)
    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/scaling.mpl:18", "s_scaling_2 (corte superior)", "s", 15.0,
                  lambda s: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_s(1e-4, s)]])}, True, False, None))
    cases.append(("1_alta", "LC-wPBE (hyb_gga_xc_lc_wpbe)", 478, 1, "maple/scaling.mpl:18", "s_scaling_2 (corte superior)", "s", 15.0,
                  lambda s: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_s(1e-4, s)]])}, True, False, None))

    # 2c. wpbeh_EG at s = EGscut = 0.08 (shadowed by wpbeh_scut = 0.13!)
    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/gga_exc/gga_x_wpbeh.mpl:121", "wpbeh_EG (sombra bajo 0.13)", "s", 0.08,
                  lambda s: {"rho": np.array([[rho_hse]]), "sigma": np.array([[unpol_sigma_from_s(rho_hse, s)]])}, True, False, None))

    # 2d. term2345_smooth + f_wpbeh_quad at s = wpbeh_scut = 0.13
    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/gga_exc/gga_x_wpbeh.mpl:208", "term2345_smooth", "s", 0.13,
                  lambda s: {"rho": np.array([[rho_hse]]), "sigma": np.array([[unpol_sigma_from_s(rho_hse, s)]])}, True, False, None))

    # 2e. t10 at aux4(w, s) = wpbeh_D = 0.57786348 (fix s=0.2, vary w via rho for HSE06 omega=0.11)
    def _pt_hse_aux4(aux4_val, omega=0.11, s=0.2):
        Ha1, Ha2, Ha3, Ha4, Ha5 = 0.00979681, 0.0410834, 0.187440, 0.00120824, 0.0347188
        Hs = (Ha1*s**2 + Ha2*s**4)/(1.0 + Ha3*s**4 + Ha4*s**5 + Ha5*s**6)
        eb1 = 1.455915450052607
        w = math.sqrt((aux4_val - s**2 * Hs) / eb1)
        kF = omega / w
        rho = (kF**3) / (3.0 * PI**2)
        return {"rho": np.array([[rho]]), "sigma": np.array([[unpol_sigma_from_s(rho, s)]])}

    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/gga_exc/gga_x_wpbeh.mpl:232", "t10 (log vs log1p)", "aux4", 0.57786348,
                  lambda a4: _pt_hse_aux4(a4, 0.11, 0.2), True, False, None))

    # 2f. xc_E1_scaled_jet at aux5(w, s) = 18.0 (src/expint_e1.c:140) in HSE06 and LC-wPBE!
    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "src/expint_e1.c:140", "xc_E1_scaled_jet (z=18)", "aux5", 18.0,
                  lambda a5: _pt_hse_aux4(a5 * 4.0 * 1.0161144 / 9.0, 0.11, 0.2), True, False, None))
    cases.append(("1_alta", "LC-wPBE (hyb_gga_xc_lc_wpbe)", 478, 1, "src/expint_e1.c:140", "xc_E1_scaled_jet (z=18)", "aux5", 18.0,
                  lambda a5: _pt_hse_aux4(a5 * 4.0 * 1.0161144 / 9.0, 0.40, 0.2), True, False, None))

    # 2g. term1 / eb1 at w = nu = wcutoff = 14.0 (maple/gga_exc/gga_x_wpbeh.mpl:66,294)
    def _pt_hse_w(w, omega=0.11, s=0.5):
        kF = omega / w
        rho = (kF**3) / (3.0 * PI**2)
        return {"rho": np.array([[rho]]), "sigma": np.array([[unpol_sigma_from_s(rho, s)]])}

    cases.append(("1_alta", "HSE06 (hyb_gga_xc_hse06)", 428, 1, "maple/gga_exc/gga_x_wpbeh.mpl:294", "term1 / eb1 (wcutoff=14)", "w = omega/k_F", 14.0,
                  lambda w: _pt_hse_w(w, 0.11, 0.5), True, False, None))
    cases.append(("1_alta", "LC-wPBE (hyb_gga_xc_lc_wpbe)", 478, 1, "maple/gga_exc/gga_x_wpbeh.mpl:294", "term1 / eb1 (wcutoff=14)", "w = omega/k_F", 14.0,
                  lambda w: _pt_hse_w(w, 0.40, 0.5), True, False, None))

    # 3. SCAN: mgga_x_scan (263) & mgga_c_scan (267)
    rho_scan = 0.1
    sig_scan = unpol_sigma_from_s(rho_scan, 0.4)
    # 3a. scan_f_alpha_left_cutoff (c1x=0.667 for X, c1c=0.64 for C)
    c1x, c2x, dx = 0.667, 0.8, 1.24
    cut_left_x = -math.log(EPS) / (-math.log(EPS) + c1x)
    cut_right_x = (-math.log(EPS / dx) + c2x) / (-math.log(EPS / dx))
    for label, acut, loc in [
        ("scan_f_alpha_left", cut_left_x, "maple/mgga_exc/mgga_x_scan.mpl:28"),
        ("scan_f_alpha (a=1)", 1.0, "maple/mgga_exc/mgga_x_scan.mpl:33"),
        ("scan_f_alpha_right", cut_right_x, "maple/mgga_exc/mgga_x_scan.mpl:32"),
    ]:
        cases.append(("1_alta", "SCAN_X (mgga_x_scan)", 263, 1, loc, label, "alpha", acut,
                      lambda a: {"rho": np.array([[rho_scan]]), "sigma": np.array([[sig_scan]]),
                                 "lapl": np.array([[0.0]]), "tau": np.array([[unpol_tau_from_alpha(rho_scan, sig_scan, a)]])},
                      True, False, None))
    # 3b. scan_gx (m_one_minus_recexp) at p^(1/4)/4.9479 = -1/log(EPS) => s = (4.9479 / (-log(EPS)))^2
    s_cut_scangx = (4.9479 / (-math.log(EPS))) ** 2
    cases.append(("1_alta", "SCAN_X (mgga_x_scan)", 263, 1, "maple/mgga_exc/mgga_x_scan.mpl:63", "scan_gx / m_one_minus_recexp", "s", s_cut_scangx,
                  lambda s: {"rho": np.array([[rho_scan]]), "sigma": np.array([[unpol_sigma_from_s(rho_scan, s)]]),
                             "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(rho_scan, unpol_sigma_from_s(rho_scan, s), 0.5)]])},
                  True, False, None))

    # 4. r2SCAN: mgga_x_r2scan (497) & mgga_c_r2scan (498) at alpha_reg = 0 and alpha_reg = 2.5
    eta_r2 = 0.001
    for fid, fn in [(497, "r2SCAN_X (mgga_x_r2scan)"), (498, "r2SCAN_C (mgga_c_r2scan)")]:
        # alpha_reg = 0 requires disable_fhc=True so alpha_reg < 0 is not clamped by FHC!
        cases.append(("1_alta", fn, fid, 1, "maple/scan_f_alpha.mpl:29", "r2scan_f_alpha (a=0)", "alpha_reg", 0.0,
                      lambda a: {"rho": np.array([[rho_scan]]), "sigma": np.array([[sig_scan]]),
                                 "lapl": np.array([[0.0]]), "tau": np.array([[unpol_tau_from_alpha(rho_scan, sig_scan, a, eta=eta_r2)]])},
                      True, True, None))
        cases.append(("1_alta", fn, fid, 1, "maple/scan_f_alpha.mpl:30", "r2scan_f_alpha (a=2.5)", "alpha_reg", 2.5,
                      lambda a: {"rho": np.array([[rho_scan]]), "sigma": np.array([[sig_scan]]),
                                 "lapl": np.array([[0.0]]), "tau": np.array([[unpol_tau_from_alpha(rho_scan, sig_scan, a, eta=eta_r2)]])},
                      True, False, None))

    # 5. TPSS: mgga_c_tpss (231) & hyb_mgga_xc_tpssh (457)
    # 5a. tpss_C0 at 1 - |z| = 1e-12 WITH uniform polarization (grad zeta = 0) vs WITH non-zero |grad zeta|!
    def _pt_tpss_c0(omz, nonprop_grad=True):
        # 1 - z = omz => z = 1 - omz, rho_tot = 0.2 => rho_a = 0.1*(2 - omz), rho_b = 0.1*omz
        rho_tot = 0.2
        ra = 0.5 * rho_tot * (2.0 - omz)
        rb = 0.5 * rho_tot * omz
        if nonprop_grad:
            # Physical nodal-like scaling: sigma_bb ~ c * rho_b (so |grad psi_b| finite at node!)
            saa = 0.01
            sbb = 0.05 * rb  # vanishes linearly with rb, obeying |grad rb|^2 = 4 rb |grad sqrt(rb)|^2!
            sab = 0.0
        else:
            # Proportional gradients: grad zeta = 0 identically
            g = 0.1
            saa = (ra / rho_tot) ** 2 * g
            sbb = (rb / rho_tot) ** 2 * g
            sab = (ra * rb / rho_tot**2) * g
        # Ensure tau > tau_W in each channel (and tau_tot > tau_W_tot so aux < 1, e.g. aux ~ 0.5)
        ta = saa / (8.0 * ra) + 0.05 * K_FACTOR_C * (ra ** (5.0 / 3.0))
        tb = sbb / (8.0 * rb) + 0.05 * K_FACTOR_C * (rb ** (5.0 / 3.0))
        return {
            "rho": np.array([[ra, rb]]),
            "sigma": np.array([[saa, sab, sbb]]),
            "lapl": np.array([[0.0, 0.0]]),
            "tau": np.array([[ta, tb]]),
        }

    cases.append(("1_alta", "TPSS_C (mgga_c_tpss, grad_zeta!=0)", 231, 2, "maple/tpss_c.mpl:34", "tpss_C0 (1-|z|=1e-12, grad_z!=0)", "1 - |zeta|", 1e-12,
                  lambda d: _pt_tpss_c0(d, nonprop_grad=True), True, False, (1e-14, 1e-15)))
    cases.append(("1_alta", "TPSS_C (mgga_c_tpss, grad_zeta=0)", 231, 2, "maple/tpss_c.mpl:34", "tpss_C0 (1-|z|=1e-12, grad_z=0)", "1 - |zeta|", 1e-12,
                  lambda d: _pt_tpss_c0(d, nonprop_grad=False), True, False, (1e-14, 1e-15)))
    cases.append(("1_alta", "TPSSh (hyb_mgga_xc_tpssh, grad_z!=0)", 457, 2, "maple/tpss_c.mpl:34", "tpss_C0 (1-|z|=1e-12, grad_z!=0)", "1 - |zeta|", 1e-12,
                  lambda d: _pt_tpss_c0(d, nonprop_grad=True), True, False, (1e-14, 1e-15)))

    # 5b. tpss_aux = m_min(xt^2 / (8*t_total), 1) at tau = tau_W (ratio = 1.0, disable_fhc=True)
    def _pt_tpss_aux(ratio):
        rho = 0.1
        sig = 0.02
        tau_w = sig / (8.0 * rho)
        tau = tau_w / ratio
        return {"rho": np.array([[rho]]), "sigma": np.array([[sig]]), "lapl": np.array([[0.0]]), "tau": np.array([[tau]])}

    cases.append(("1_alta", "TPSS_C (mgga_c_tpss)", 231, 1, "maple/tpss_c.mpl:40", "tpss_aux m_min(..., 1)", "tau_W / tau", 1.0,
                  _pt_tpss_aux, True, True, None))

    # 6. B3LYP (402) / gga_x_b88 (106): x_asinh_cx(1, p) at p = x^2 = 1/16 (x = 0.25)
    cases.append(("1_alta", "B3LYP (hyb_gga_xc_b3lyp)", 402, 1, "maple/util.mpl:192", "x_asinh_cx (p = 1/16)", "p = xs^2", 1.0 / 16.0,
                  lambda p: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, math.sqrt(p))]])},
                  True, False, None))

    # 7. wB97X (464, omega=0.3): attenuation_erf(a) at a = 1/16, 1/2, 16777216
    # In lda_x_erf.mpl: a = a_cnst * rs / (1+z)^(1/3), a_cnst = (4/(9*pi))^(1/3) * omega / 2; for z=0, opz=1 => rs = a / a_cnst
    a_cnst_wb97x = ((4.0 / (9.0 * PI)) ** (1.0 / 3.0)) * 0.3 / 2.0
    for acut, lbl, loc in [
        (1.0 / 16.0, "attenuation_erf (sr_cutoff=1/16)", "maple/enforce_smooth_cutoffs.mpl:24"),
        (0.5, "attenuation_erf (mid_cutoff=1/2)", "maple/enforce_smooth_cutoffs.mpl:21"),
        (16777216.0, "attenuation_erf (lr_cutoff=2^24)", "maple/enforce_smooth_cutoffs.mpl:25"),
    ]:
        cases.append(("1_alta", "wB97X (hyb_gga_xc_wb97x)", 464, 1, loc, lbl, "a = omega/(2*k_F)", acut,
                      lambda a: {"rho": np.array([[rho_from_rs(a / a_cnst_wb97x)]]),
                                 "sigma": np.array([[unpol_sigma_from_s(rho_from_rs(a / a_cnst_wb97x), 0.3)]])},
                      True, False, None))

    # =========================================================================
    # TIER 2: SHARED INFRASTRUCTURE HELPERS IN maple/util.mpl & attenuation.mpl
    # =========================================================================

    # 8. e1_scaled_G(mu=100) in gga_c_ft97 (114), plus xc_E1_scaled_jet(mu=18) in gga_c_ft97 (114)!
    # In gga_c_ft97.mpl: for z=0, rsa = rs * 2^(1/3), ga2 = 0 (at sigma -> tiny), mu_ba(rsa, 0) = C1 * rsa / kssp0(rsa)^2
    # Let's solve mu_ba(rsa, 0) = target_mu numerically with Newton!
    C0_ft = (1.0 - math.log(2.0)) / (2.0 * PI**2)
    C1_ft = 4.0 * C0_ft / 3.0
    def _rsa_for_ft97_mu(target_mu):
        rsa = target_mu * (1.291551074 - 0.349064173)**2 / C1_ft
        for _ in range(30):
            k = 1.291551074 + 0.349064173 * math.expm1(-0.08327588 * (rsa ** 0.8))
            dk = -0.349064173 * 0.08327588 * 0.8 * (rsa ** -0.2) * math.exp(-0.08327588 * (rsa ** 0.8))
            mu = C1_ft * rsa / (k * k)
            dmu = C1_ft / (k * k) - 2.0 * C1_ft * rsa * dk / (k ** 3)
            rsa -= (mu - target_mu) / dmu
        return rsa / (2.0 ** (1.0 / 3.0))

    cases.append(("2_infra", "FT97_C (gga_c_ft97)", 114, 1, "maple/util.mpl:41", "e1_scaled_G (mu=100)", "mu_ba", 100.0,
                  lambda mu: {"rho": np.array([[rho_from_rs(_rsa_for_ft97_mu(mu))]]), "sigma": np.array([[0.0]])},
                  True, False, None))
    cases.append(("2_infra", "FT97_C (gga_c_ft97)", 114, 1, "src/expint_e1.c:140", "xc_E1_scaled_jet (mu=18)", "mu_ba", 18.0,
                  lambda mu: {"rho": np.array([[rho_from_rs(_rsa_for_ft97_mu(mu))]]), "sigma": np.array([[0.0]])},
                  True, False, None))

    # 9. m_saturate2(y=1) in gga_x_rge2 (142): y = (10/81)*t/0.8040 = 1, where t = X2S^2 * p => s^2 = 0.8040 / (10/81)
    s_rge2 = math.sqrt(0.8040 / (10.0 / 81.0))
    cases.append(("2_infra", "RGE2_X (gga_x_rge2)", 142, 1, "maple/util.mpl:175", "m_saturate2 (y=1)", "s", s_rge2,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    by_name = {lxc.name(i): i for i in lxc.all_ids()}

    # 10. m_expm1_over(y=1/16) in hyb_gga_x_pbe_erf_laikov & gga_x_vmt_pbe
    laikov_gamma = 2.14992827678013751815
    cases.append(("2_infra", "Laikov_SR_X (hyb_gga_x_pbe_erf_laikov)", by_name["hyb_gga_x_pbe_erf_laikov"], 1,
                  "maple/util.mpl:203", "m_expm1_over (gamma*s^2=1/16)", "s", math.sqrt(1.0 / (16.0 * laikov_gamma)),
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    # 11. m_bound(|y|=1) in hyb_mgga_xc_b98: b98_gamma_ss * q = 1, q = (K - t + x^2/8 + u/4)/K => at x=0, u=0: t = K*(1 - q)
    def _pt_b98_y(y_val, gamma=1.6):
        rho = 0.1
        q = y_val / gamma
        t_s = K_FACTOR_C * (1.0 - q)
        tau = 2.0 * t_s * ((rho / 2.0) ** (5.0 / 3.0))
        return {"rho": np.array([[rho]]), "sigma": np.array([[0.0]]), "lapl": np.array([[0.0]]), "tau": np.array([[tau]])}

    cases.append(("2_infra", "B98 (hyb_mgga_xc_b98)", by_name["hyb_mgga_xc_b98"], 1,
                  "maple/util.mpl:219", "m_bound (gamma_ss*q = +1)", "gamma*q", 1.0,
                  _pt_b98_y, True, True, None))

    # 12. mgga_series_w at 2*K/(K + t) = 1/5 => t = 9*K_FACTOR_C
    def _pt_m06l_w(v):
        rho = 0.1
        sig = unpol_sigma_from_s(rho, 0.5)
        t_s = K_FACTOR_C * (2.0 / v - 1.0)
        tau = 2.0 * t_s * ((rho / 2.0) ** (5.0 / 3.0))
        return {"rho": np.array([[rho]]), "sigma": np.array([[sig]]), "lapl": np.array([[0.0]]), "tau": np.array([[tau]])}

    cases.append(("2_infra", "M06-L_X (mgga_x_m06_l)", by_name["mgga_x_m06_l"], 1, "maple/util.mpl:636", "mgga_series_w (1+w = 1/5)", "1 + w", 0.2,
                  _pt_m06l_w, True, False, None))

    # 12b. attenuation_yukawa in lda_x_yukawa (omega=0.3): a = a_cnst * rs at a = 1/16 and a = 1/2
    for acut, lbl, loc in [
        (1.0 / 16.0, "attenuation_yukawa (sr=1/16)", "maple/attenuation.mpl:134"),
        (0.5, "attenuation_yukawa (mid=1/2)", "maple/attenuation.mpl:134"),
    ]:
        cases.append(("2_infra", "Yukawa_LDA_X (lda_x_yukawa)", by_name["lda_x_yukawa"], 1, loc, lbl, "a", acut,
                      lambda a: {"rho": np.array([[rho_from_rs(a / a_cnst_wb97x)]])},
                      True, False, None))

    # 12c. attenuation_erf_f2 / f3 in hyb_mgga_x_pjs18 (omega=0.33? check at s=0 so tm_f0(0)=1, a = a_cnst * rs!)
    # Get exact cam_omega of hyb_mgga_x_pjs18 via C API:
    lxc.L.xc_hyb_cam_coef.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    _fp = lxc.Func(by_name["hyb_mgga_x_pjs18"], 1)
    _om, _al, _be = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    lxc.L.xc_hyb_cam_coef(ctypes.c_void_p(_fp.p), ctypes.byref(_om), ctypes.byref(_al), ctypes.byref(_be))
    _fp.close()
    a_cnst_pjs18 = ((4.0 / (9.0 * PI)) ** (1.0 / 3.0)) * float(_om.value) / 2.0
    for acut, lbl, loc in [
        (1.0 / 16.0, "attenuation_erf_f2/f3 (sr=1/16)", "maple/attenuation.mpl:67"),
        (0.25, "attenuation_erf_f2/f3 (mid=1/4)", "maple/attenuation.mpl:81"),
    ]:
        def _pt_pjs18(a_val):
            # At s = 0.2, tm_f0(p): let's use s = 0.0 so p = 0 => tm_f0(0) = 1.0, js18_H(0)=0, js18_G(0, t) != 0 when t != K!
            # Wait: to activate BOTH f2 and f3, let s = 0.2 -> compute exact tm_f0(0.04) from mgga_x_tm.mpl?
            # Simpler: at s=0, f2 is active with t = 1.5*K! And at a = 0.25, f2 and f3 both switch at mid_cutoff = 1/4!
            rs = a_val / a_cnst_pjs18
            rho = rho_from_rs(rs)
            tau = 2.0 * (1.5 * K_FACTOR_C) * ((rho / 2.0) ** (5.0 / 3.0))
            return {"rho": np.array([[rho]]), "sigma": np.array([[0.0]]), "lapl": np.array([[0.0]]), "tau": np.array([[tau]])}
        cases.append(("2_infra", "PJS18_X (hyb_mgga_x_pjs18)", by_name["hyb_mgga_x_pjs18"], 1, loc, lbl, "A", acut,
                      _pt_pjs18, True, False, None))

    # =========================================================================
    # TIER 3: ALL FUNCTIONAL-SPECIFIC PIECEWISE / SERIES / CLAMP SWITCHES
    # =========================================================================

    # 13. LB94 (gga_x_lb, ID 160): p < 300^2 (xs = 300) in maple/gga_vxc/gga_x_lb.mpl:31
    cases.append(("3_sec", "LB94 (gga_x_lb)", by_name["gga_x_lb"], 1, "maple/gga_vxc/gga_x_lb.mpl:31", "lb_f0 (asintota p < 300^2)", "xs", 300.0,
                  lambda xs: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_xs(1e-4, xs)]])},
                  True, False, None))

    # 14. BR89 (mgga_x_br89, ID 206; mgga_x_b00; mgga_x_br89_explicit, ID 586):
    #     br89_Q = (u - 4*gamma*t + gamma*p/2)/6, gamma = 0.8 (for br89).
    #     Test BOTH Q = 0 and Q = +5e-13 and Q = -5e-13!
    def _pt_br89_Q(Q_val, gamma=0.8):
        rho = 0.1
        xs = 0.5
        sig = unpol_sigma_from_xs(rho, xs)
        t_s = 1.2 * K_FACTOR_C
        tau = 2.0 * t_s * ((rho / 2.0) ** (5.0 / 3.0))
        u_s = 6.0 * Q_val + 4.0 * gamma * t_s - 0.5 * gamma * (xs ** 2)
        lapl = unpol_lapl_from_u(rho, u_s)
        return {"rho": np.array([[rho]]), "sigma": np.array([[sig]]), "lapl": np.array([[lapl]]), "tau": np.array([[tau]])}

    for fn_key, gamma_val in [("mgga_x_br89", 0.8), ("mgga_x_br89_explicit", 0.8)]:
        fid_br = by_name[fn_key]
        cases.append(("3_sec", fn_key, fid_br, 1, "maple/mgga_exc/mgga_x_br89.mpl:37", "br89_cQ (salto en Q=0)", "Q", 0.0,
                      lambda q, g=gamma_val: _pt_br89_Q(q, g), True, False, (1e-14, 1e-16)))
        cases.append(("3_sec", fn_key, fid_br, 1, "maple/mgga_exc/mgga_x_br89.mpl:36", "br89_cQ (banda |Q|=5e-13)", "Q", 5.0e-13,
                      lambda q, g=gamma_val: _pt_br89_Q(q, g), True, False, (1e-14, 1e-16)))

    # 15. Perdew-Zunger 1981 (lda_c_pz, ID 9): rs = 1.0 (maple/lda_exc/lda_c_pz.mpl:35)
    cases.append(("3_sec", "PZ81 (lda_c_pz)", by_name["lda_c_pz"], 1, "maple/lda_exc/lda_c_pz.mpl:35", "ec (rs = 1)", "rs", 1.0,
                  lambda rs: {"rho": np.array([[rho_from_rs(rs)]])}, True, False, None))

    # 16. Gordon-Kim 1972 (lda_c_gk72): rs = 0.7 and rs = 10.0 (maple/lda_exc/lda_c_gk72.mpl:22)
    for rcut in (0.7, 10.0):
        cases.append(("3_sec", "GK72 (lda_c_gk72)", by_name["lda_c_gk72"], 1, "maple/lda_exc/lda_c_gk72.mpl:22", f"f (rs = {rcut})", "rs", rcut,
                      lambda rs: {"rho": np.array([[rho_from_rs(rs)]])}, True, False, None))

    # 17. Hedin-Lundqvist (lda_c_hl, r=21): rs / 21 = 4 => rs = 84.0 (maple/lda_exc/lda_c_hl.mpl:46)
    cases.append(("3_sec", "HL (lda_c_hl)", by_name["lda_c_hl"], 1, "maple/lda_exc/lda_c_hl.mpl:46", "hl_f0 (rs/r = 4)", "rs", 84.0,
                  lambda rs: {"rho": np.array([[rho_from_rs(rs)]])}, True, False, None))

    # 18. Paziani et al. 2006 (lda_c_pw_erf / lda_c_pmgb06): 0.7524 * rs = 1 (maple/lda_exc/lda_c_pw_erf.mpl:83)
    cases.append(("3_sec", "PMGB06 (lda_c_pmgb06)", by_name["lda_c_pmgb06"], 1, "maple/lda_exc/lda_c_pw_erf.mpl:83", "lda_c_pw_erf_g0c (0.7524*rs=1)", "rs", 1.0 / 0.7524,
                  lambda rs: {"rho": np.array([[rho_from_rs(rs)]])}, True, False, None))

    # 19. Proynov-Kong 2009 (lda_c_pk09): pole screen |f_r_den(rs)| <= 1e-12 near rs ~ 22.7418!
    # Find exact root of b_24_i polynomial:
    b24 = [-109.74263493216910, 16.2663129444242415, 54.4034331373908366, -25.154009904187990, 1.0]
    rs_pole = 22.74
    for _ in range(30):
        val = sum(b24[i] * rs_pole**i for i in range(5))
        der = sum(i * b24[i] * rs_pole**(i-1) for i in range(1, 5))
        rs_pole -= val / der
    # Boundary of |f_r_den(rs)| = 1e-12 is at rs_cut = rs_pole + 1e-12 / |der_pole|
    der_pole = sum(i * b24[i] * rs_pole**(i-1) for i in range(1, 5))
    rs_pk09_cut = rs_pole + 1e-12 / abs(der_pole)
    cases.append(("3_sec", "PK09 (lda_c_pk09, polo f_r)", by_name["lda_c_pk09"], 2, "maple/lda_exc/lda_c_pk09.mpl:25", "f_r (|den|<=1e-12 cerca rs=22.74)", "rs", rs_pk09_cut,
                  lambda rs: {"rho": np.array([[0.7 * rho_from_rs(rs), 0.3 * rho_from_rs(rs)]])}, True, False, (1e-15, 1e-16)))

    # 20. 2D AMGB (lda_c_2d_amgb): |zeta| = 0.5 (maple/lda_exc/lda_c_2d_amgb.mpl:48)
    cases.append(("3_sec", "2D-AMGB (lda_c_2d_amgb)", by_name["lda_c_2d_amgb"], 2, "maple/lda_exc/lda_c_2d_amgb.mpl:48", "amgb_ex6_residual (|z|=0.5)", "zeta", 0.5,
                  lambda z: {"rho": np.array([[0.05 * (1.0 + z), 0.05 * (1.0 - z)]])}, True, False, None))

    # 21. Relativistic LDA exchange (lda_x_rel): beta(rs) = (9*pi/4)^(1/3)/(rs*M_C) = 0.5
    rs_rel_cut = ((9.0 * PI / 4.0) ** (1.0 / 3.0)) / (0.5 * M_C)
    cases.append(("3_sec", "LDA_X_REL (lda_x_rel)", by_name["lda_x_rel"], 1, "maple/lda_exc/lda_x_rel.mpl:34", "lda_x_rel_bracket (b=0.5)", "rs", rs_rel_cut,
                  lambda rs: {"rho": np.array([[rho_from_rs(rs)]])}, True, False, None))

    # 22. revTCA (gga_c_revtca): msinc at revtca_aa * s / rs = EPS^(1/4) (polarized z=0.5 so z^4 != 0)
    revtca_aa = PI * ((9.0 * PI / 4.0) ** (1.0 / 3.0))
    sinc_cut = EPS ** 0.25
    def _pt_revtca(x):
        rs = 1.0
        rho = rho_from_rs(rs)
        s_tot = x * rs / revtca_aa
        xt = s_tot / XT2S
        sig_tot = (xt * (rho ** (4.0 / 3.0))) ** 2
        z = 0.5
        ra, rb = 0.5 * rho * (1 + z), 0.5 * rho * (1 - z)
        return {"rho": np.array([[ra, rb]]), "sigma": np.array([[0.25*sig_tot, 0.25*sig_tot, 0.25*sig_tot]])}

    cases.append(("3_sec", "revTCA (gga_c_revtca)", by_name["gga_c_revtca"], 2, "maple/gga_exc/gga_c_revtca.mpl:24", "msinc (x=eps^(1/4))", "x = aa*s/rs", sinc_cut,
                  _pt_revtca, True, False, None))

    # 23. W94 (gga_c_w94): 1 - |z| = 0.5 => z = 0.5 (maple/gga_exc/gga_c_w94.mpl:30)
    cases.append(("3_sec", "W94 (gga_c_w94)", by_name["gga_c_w94"], 2, "maple/gga_exc/gga_c_w94.mpl:30", "w94_log_az (1-|z|=0.5)", "zeta", 0.5,
                  lambda z: {"rho": np.array([[0.05 * (1 + z), 0.05 * (1 - z)]]), "sigma": np.array([[0.01, 0.005, 0.01]])},
                  True, False, None))

    # 24. zPBEint (gga_c_zpbeint, alpha=2.4): zpbeint_u = alpha * t^3 * log(mphi(z)) = -0.5
    def _pt_zpbeint(u_val):
        # Note u_val ~ -0.5 (negative)
        z = 0.6
        rs = 1.0
        rho = rho_from_rs(rs)
        mphi = 0.5 * ((1 + z)**(2.0/3.0) + (1 - z)**(2.0/3.0))
        alpha_p = 2.4
        t = (u_val / (alpha_p * math.log(mphi))) ** (1.0 / 3.0)
        xt2 = (t ** 2) * 16.0 * (2.0 ** (2.0 / 3.0)) * (mphi ** 2) * rs
        sig_tot = xt2 * (rho ** (8.0 / 3.0))
        ra, rb = 0.5 * rho * (1 + z), 0.5 * rho * (1 - z)
        return {"rho": np.array([[ra, rb]]), "sigma": np.array([[0.25*sig_tot, 0.25*sig_tot, 0.25*sig_tot]])}

    cases.append(("3_sec", "zPBEint (gga_c_zpbeint)", by_name["gga_c_zpbeint"], 2, "maple/gga_exc/gga_c_zpbeint.mpl:56", "f (zpbeint_u = -0.5)", "u", -0.5,
                  _pt_zpbeint, True, False, None))

    # 25. LKT kinetic (gga_k_lkt, a=1.3): a * s = 0.25 => s = 0.25 / 1.3
    cases.append(("3_sec", "LKT (gga_k_lkt)", by_name["gga_k_lkt"], 1, "maple/gga_exc/gga_k_lkt.mpl:53", "lkt_sech (a*s = 1/4)", "s", 0.25 / 1.3,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    # 26. Meyer kinetic (gga_k_meyer): w = (s/6)^2 at w = 1/16 (s=1.5), w = 1 (s=6), and p = meyer_p_cut (s = 6/sqrt(1e-3))
    for scut, lbl, loc in [
        (1.5, "meyer_one_minus_feta_over_w (w=1/16)", "maple/gga_exc/gga_k_meyer.mpl:43"),
        (6.0, "meyer_lambda (w=1)", "maple/gga_exc/gga_k_meyer.mpl:66"),
        (6.0 / math.sqrt(1e-3), "meyer_f (v=1e-3)", "maple/gga_exc/gga_k_meyer.mpl:93"),
    ]:
        cases.append(("3_sec", "Meyer (gga_k_meyer)", by_name["gga_k_meyer"], 1, loc, lbl, "s", scut,
                      lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                      True, False, None))

    # 27. BPCCAC (gga_x_bpccac): xs = 19 (maple/gga_exc/gga_x_bpccac.mpl:40)
    cases.append(("3_sec", "BPCCAC (gga_x_bpccac)", by_name["gga_x_bpccac"], 1, "maple/gga_exc/gga_x_bpccac.mpl:40", "bpccac_f (xs=19)", "xs", 19.0,
                  lambda xs: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, xs)]])},
                  True, False, None))

    # 28. Chachiyo exchange (gga_x_chachiyo): q = chachiyo_cx * xs = 1/16
    chachiyo_cx = (2.0 / 9.0) * ((PI / 3.0) ** (1.0 / 3.0)) * (2.0 ** (-1.0 / 3.0))
    cases.append(("3_sec", "Chachiyo_X (gga_x_chachiyo)", by_name["gga_x_chachiyo"], 1, "maple/gga_exc/gga_x_chachiyo.mpl:54", "chachiyo_R (q=1/16)", "xs", 1.0 / (16.0 * chachiyo_cx),
                  lambda xs: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, xs)]])},
                  True, False, None))

    # 29. FD-LB94 (gga_x_fd_lb94): X = xs / 2^(1/3) = 0.25 => xs = 0.25 * 2^(1/3)
    cases.append(("3_sec", "FD-LB94 (gga_x_fd_lb94)", by_name["gga_x_fd_lb94"], 1, "maple/gga_exc/gga_x_fd_lb94.mpl:64", "fd_f (X=1/4)", "xs", 0.25 * (2.0 ** (1.0 / 3.0)),
                  lambda xs: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, xs)]])},
                  True, False, None))

    # 30. GG99 (gga_x_gg99): p = gg99_b^(2/3) => xs = (4*sqrt(3)*pi^3)^(1/3)
    xs_gg99 = (4.0 * math.sqrt(3.0) * (PI ** 3)) ** (1.0 / 3.0)
    cases.append(("3_sec", "GG99 (gga_x_gg99)", by_name["gga_x_gg99"], 1, "maple/gga_exc/gga_x_gg99.mpl:74", "gg99_r2 (p = b^(2/3))", "xs", xs_gg99,
                  lambda xs: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, xs)]])},
                  True, False, None))

    # 31. HTBS (gga_x_htbs): s = 0.6 and s = 2.6
    for scut in (0.6, 2.6):
        cases.append(("3_sec", "HTBS (gga_x_htbs)", by_name["gga_x_htbs"], 1, "maple/gga_exc/gga_x_htbs.mpl:40", f"htbs_f0 (s={scut})", "s", scut,
                      lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                      True, False, None))

    # 32. LV-rPW86 (gga_x_lv_rpw86): s = 100
    cases.append(("3_sec", "LV-rPW86 (gga_x_lv_rpw86)", by_name["gga_x_lv_rpw86"], 1, "maple/gga_exc/gga_x_lv_rpw86.mpl:33", "enforce_smooth_lr (s=100)", "s", 100.0,
                  lambda s: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_s(1e-4, s)]])},
                  True, False, None))

    # 33. NCAP (gga_x_ncap): s = 0.25 (p = s^2 = 1/16)
    cases.append(("3_sec", "NCAP (gga_x_ncap)", by_name["gga_x_ncap"], 1, "maple/gga_exc/gga_x_ncap.mpl:34", "ncap_tanh_asinh (s^2=1/16)", "s", 0.25,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    # 34. PBEpow (gga_x_pbepow) & PBEtrans (gga_x_pbetrans): s = 1.0
    cases.append(("3_sec", "PBEpow (gga_x_pbepow)", by_name["gga_x_pbepow"], 1, "maple/gga_exc/gga_x_pbepow.mpl:35", "pbepow_f0 (s^2=1)", "s", 1.0,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))
    cases.append(("3_sec", "PBEtrans (gga_x_pbetrans)", by_name["gga_x_pbetrans"], 1, "maple/gga_exc/gga_x_pbetrans.mpl:50", "pbetrans_f0 (s^2=1)", "s", 1.0,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    # 35. PBE-GX (mgga_x_gx): alpha = 1.0 (known C0 / C1-kink by design!)
    cases.append(("3_sec", "PBE-GX (mgga_x_gx)", by_name["mgga_x_gx"], 1, "maple/mgga_exc/mgga_x_gx.mpl:37", "gx_f_a (alpha=1, diseno)", "alpha", 1.0,
                  lambda a: {"rho": np.array([[0.1]]), "sigma": np.array([[sig_scan]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(0.1, sig_scan, a)]])},
                  True, False, None))

    # 36. regTPSS (mgga_x_regtpss): |alpha| = 0.5 and 3*s^2 = 1 => s = 1/sqrt(3)
    cases.append(("3_sec", "regTPSS_X (mgga_x_regtpss)", by_name["mgga_x_regtpss"], 1, "maple/mgga_exc/mgga_x_regtpss.mpl:37", "reg_one_minus_f_a (|a|=0.5)", "alpha", 0.5,
                  lambda a: {"rho": np.array([[0.1]]), "sigma": np.array([[sig_scan]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(0.1, sig_scan, a)]])},
                  True, False, None))
    cases.append(("3_sec", "regTPSS_X (mgga_x_regtpss)", by_name["mgga_x_regtpss"], 1, "maple/mgga_exc/mgga_x_regtpss.mpl:67", "reg_one_minus_A (3*s^2=1)", "s", 1.0 / math.sqrt(3.0),
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(0.1, unpol_sigma_from_s(0.1, s), 0.3)]])},
                  True, False, None))

    # 37. LAK correlation (mgga_c_lak): alpha = 0 (maple/mgga_exc/mgga_c_lak.mpl:61) and lak_beta_t(rs) = 0 (line 79)
    cases.append(("3_sec", "LAK_C (mgga_c_lak)", by_name["mgga_c_lak"], 1, "maple/mgga_exc/mgga_c_lak.mpl:61", "lak_fc (alpha=0)", "alpha", 0.0,
                  lambda a: {"rho": np.array([[0.1]]), "sigma": np.array([[sig_scan]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(0.1, sig_scan, a)]])},
                  True, True, None))
    # Find root of lak_mu_sc(rs) = 0 => lak_Cma(rs)*(1 + 6*lak_mu_a) = 1 + 6*lak_mu_ax
    lak_h0x = 1.174
    lak_mu_ax = -(97.0 + 3.0*lak_h0x + math.sqrt((3.0*lak_h0x)**2 + 74166.0*lak_h0x - 64175.0)) / 1200.0
    lak_mu_a = -lak_mu_ax / 2.0
    lak_Cs0 = -16.0 * PI * ((3.0 * PI**2) ** (1.0 / 3.0)) * 2.568 / (3000.0 * (10.0 / 81.0))
    target_cma = (1.0 + 6.0 * lak_mu_ax) / (1.0 + 6.0 * lak_mu_a)
    # Check if lak_Cma(rs) crosses target_cma for rs > 0
    # lak_Cs0 < 0, target_cma: 1 + 6*(-0.246) = -0.476, 1 + 3*(0.246) = 1.738 => target_cma < 0!
    rs_lak = 5.0
    for _ in range(50):
        def _cma(r):
            return lak_Cs0 * (1.0 + 0.1 * r**0.65) / ((1.0 + 0.065 * r**0.9) * (1.0 + 0.03 * r**1.2))
        f0 = _cma(rs_lak) - target_cma
        df = (_cma(rs_lak + 1e-7) - _cma(rs_lak - 1e-7)) / 2e-7
        rs_lak -= f0 / df
    if rs_lak > 0:
        cases.append(("3_sec", "LAK_C (mgga_c_lak)", by_name["mgga_c_lak"], 1, "maple/mgga_exc/mgga_c_lak.mpl:79", "lak_H1 (lak_beta_t(rs)=0)", "rs", rs_lak,
                      lambda rs: {"rho": np.array([[rho_from_rs(rs)]]), "sigma": np.array([[unpol_sigma_from_s(rho_from_rs(rs), 0.5)]]),
                                  "lapl": np.array([[0.0]]),
                                  "tau": np.array([[unpol_tau_from_alpha(rho_from_rs(rs), unpol_sigma_from_s(rho_from_rs(rs), 0.5), 0.5)]])},
                      True, False, None))

    # 38. FT98 (mgga_x_ft98): q3 = xs^4 - u^2 - b2 at q3 = eps^(1/4)
    # b2 in mgga_x_ft98.c: let's read or test mgga_x_ms (ms_amax = 1e6)!
    cases.append(("3_sec", "MS0 (mgga_x_ms0)", by_name["mgga_x_ms0"], 1, "maple/mgga_exc/mgga_x_ms.mpl:29", "ms_fa clamp (alpha=1e6)", "alpha", 1e6,
                  lambda a: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_s(1e-4, 0.2)]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(1e-4, unpol_sigma_from_s(1e-4, 0.2), a)]])},
                  True, False, None))

    # 39. mBEEF (mgga_x_mbeef): enforce_smooth_lr at alpha = 1e4
    cases.append(("3_sec", "mBEEF (mgga_x_mbeef)", by_name["mgga_x_mbeef"], 1, "maple/mbeef.mpl:38", "mbeef_one_minus_xj (alpha=1e4)", "alpha", 1e4,
                  lambda a: {"rho": np.array([[1e-4]]), "sigma": np.array([[unpol_sigma_from_s(1e-4, 0.2)]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(1e-4, unpol_sigma_from_s(1e-4, 0.2), a)]])},
                  True, False, None))

    # 40. CASE21 (hyb_gga_xc_case21): B-spline interior knot u = 3/7 in xbspline (u = gammax*s^2 / (1 + gammax*s^2))
    gammax_c21 = (0.2195149727645171 / 0.8040)
    u_knot = 3.0 / 7.0  # k=3, Nsp=10 => qmin=-3/7, qmax=10/7, dq=1/7 => knots at m/7 for m=-3..10!
    s_c21 = math.sqrt((u_knot / (1.0 - u_knot)) / gammax_c21)
    cases.append(("3_sec", "CASE21 (hyb_gga_xc_case21)", by_name["hyb_gga_xc_case21"], 1, "src/hyb_gga_xc_case21.c:53", "xc_bspline (nudo u=3/7)", "s", s_c21,
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]])},
                  True, False, None))

    # 41. PC07 / PC07opt / SCAN-L / r2SCAN-L (maple/mgga_exc/mgga_k_pc07.mpl:27):
    #     At p = 0 (sigma = 0): f_W = 0, Delta = 8*q^2/81, B = 1 + 20*q/9, GE4 = 1 + 20*q/9 + 8*q^2/81, S = sqrt(1 + Delta^2).
    #     Thus pc07_z(0, q) = (Delta + B)^2 / (S * GE4) = GE4 / sqrt(1 + (8*q^2/81)^2)!
    #     Solve pc07_z(0, q) = target_z with Newton to machine precision!
    def _q_for_pc07_z(target_z):
        q = (target_z - 1.0) * 9.0 / 20.0
        for _ in range(40):
            def _zfn(qq):
                dlt = 8.0 * (qq ** 2) / 81.0
                ge4 = 1.0 + 20.0 * qq / 9.0 + dlt
                return ge4 / math.sqrt(1.0 + dlt * dlt)
            f0 = _zfn(q) - target_z
            df = (_zfn(q + 1e-8) - _zfn(q - 1e-8)) / 2e-8
            q -= f0 / df
        return q

    def _pt_pc07_z(target_z):
        rho = 0.1
        q = _q_for_pc07_z(target_z)
        u_s = q / (X2S ** 2)
        return {"rho": np.array([[rho]]), "sigma": np.array([[0.0]]), "lapl": np.array([[unpol_lapl_from_u(rho, u_s)]]), "tau": np.array([[0.1]])}

    # Test original PC07 (a=0.5389, b=3 -> zlo=0.5389/40, zhi=39*0.5389/40) AND PC07_OPT / SCAN-L / r2SCAN-L (a=1.784720, b=0.258304)!
    cases.append(("3_sec", "PC07 (mgga_k_pc07, b=3)", by_name["mgga_k_pc07"], 1, "maple/mgga_exc/mgga_k_pc07.mpl:27",
                  "pc07_fab (zlo = a/40, b=3)", "z_pc07", 0.5389 / 40.0, _pt_pc07_z, True, False, None))
    for prio_l, fn_l, label_l in [
        ("1_alta", "mgga_x_r2scanl", "r2SCAN-L_X (mgga_x_r2scanl)"),
        ("1_alta", "mgga_c_r2scanl", "r2SCAN-L_C (mgga_c_r2scanl)"),
        ("1_alta", "mgga_x_scanl", "SCAN-L_X (mgga_x_scanl)"),
        ("1_alta", "mgga_c_scanl", "SCAN-L_C (mgga_c_scanl)"),
        ("3_sec", "mgga_k_pc07_opt", "PC07-opt (mgga_k_pc07_opt)"),
    ]:
        cases.append((prio_l, label_l, by_name[fn_l], 1, "maple/mgga_exc/mgga_k_pc07.mpl:27",
                      "pc07_fab (zlo = a/40, b=0.258)", "z_pc07", 1.784720 / 40.0, _pt_pc07_z, True, False, None))
    cases.append(("1_alta", "r2SCAN-L_X (mgga_x_r2scanl)", by_name["mgga_x_r2scanl"], 1, "maple/mgga_exc/mgga_k_pc07.mpl:27",
                  "pc07_fab (zhi = 39*a/40)", "z_pc07", 39.0 * 1.784720 / 40.0, _pt_pc07_z, True, False, None))

    # 42. TB09 (mgga_x_tb09, vxc-only): br89_cQ at Q=0 and Q=5e-13 (gamma=0.8)
    cases.append(("3_sec", "TB09 (mgga_x_tb09)", by_name["mgga_x_tb09"], 1, "maple/mgga_exc/mgga_x_br89.mpl:37",
                  "br89_cQ (salto en Q=0)", "Q", 0.0, lambda q: _pt_br89_Q(q, 0.8), True, False, (1e-14, 1e-16)))
    cases.append(("3_sec", "TB09 (mgga_x_tb09)", by_name["mgga_x_tb09"], 1, "maple/mgga_exc/mgga_x_br89.mpl:36",
                  "br89_cQ (banda |Q|=5e-13)", "Q", 5.0e-13, lambda q: _pt_br89_Q(q, 0.8), True, False, (1e-14, 1e-16)))

    # 43. GG99 inner series switches: gg99_asinh_sq (v = 1/16) and gg99_f0 (w = r^2 = 1/16)
    gg99_b = 4.0 * math.sqrt(3.0) * (PI ** 3)
    gg99_a = 1.0 / ((6.0 ** 0.5) * PI)  # check exact value via root finding on v(p) = 1/16? Wait: read gg99_a!
    # Let's check gg99_a: in gga_x_gg99, at small p, r ~ asinh(sqrt(v)) -> let's also test FT98, EDMGGA, EEL, CSK1, LAK_X, RPP, PBE_ERF_GWS!

    # 44. FT98 (mgga_x_ft98): q3 = xs^4 - u^2 - b2 (b2 = 66.364138) at q3 = +eps^(1/4) (fix u = 0 => xs = (b2 + q3)^(1/4))
    ft98_b2 = 66.364138
    cases.append(("3_sec", "FT98 (mgga_x_ft98)", by_name["mgga_x_ft98"], 1, "maple/mgga_exc/mgga_x_ft98.mpl:60",
                  "ft98_q20 (|q3| = eps^(1/4))", "q3", EPS ** 0.25,
                  lambda q3: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_xs(0.1, (ft98_b2 + q3) ** 0.25)]]),
                              "lapl": np.array([[0.0]]), "tau": np.array([[0.2]])},
                  True, False, None))

    # 45. EDMGGA (mgga_x_edmgga): Qb = b98_q(xs^2, u, t) = -1e4 / 0.704 (at xs=0, u=0 => t = K*(1 - Qb))
    def _pt_edmgga(Qb):
        rho = 1e-3
        t_s = K_FACTOR_C * (1.0 - Qb)
        tau = 2.0 * t_s * ((rho / 2.0) ** (5.0 / 3.0))
        return {"rho": np.array([[rho]]), "sigma": np.array([[0.0]]), "lapl": np.array([[0.0]]), "tau": np.array([[tau]])}

    cases.append(("3_sec", "EDMGGA (mgga_x_edmgga)", by_name["mgga_x_edmgga"], 1, "maple/mgga_exc/mgga_x_edmgga.mpl:50",
                  "edmgga_x (Qb = -1e4/a)", "Qb", -1.0e4 / 0.704, _pt_edmgga, True, False, None))

    # 46. EEL (mgga_x_eel): eel_G(s^2) at s^2 = x0 = 0.14118847627290322299 => s = sqrt(x0)
    eel_x0 = 0.14118847627290322299
    cases.append(("3_sec", "EEL (mgga_x_eel)", by_name["mgga_x_eel"], 1, "maple/mgga_exc/mgga_x_eel.mpl:24",
                  "eel_G (s^2 = x0)", "s", math.sqrt(eel_x0),
                  lambda s: {"rho": np.array([[0.1]]), "sigma": np.array([[unpol_sigma_from_s(0.1, s)]]),
                             "lapl": np.array([[0.0]]), "tau": np.array([[unpol_tau_from_alpha(0.1, unpol_sigma_from_s(0.1, s), 1.0)]])},
                  True, False, None))

    # 47. LAK_X (mgga_x_lak): lak_fx at alpha = EPS (disable_fhc=True!)
    cases.append(("3_sec", "LAK_X (mgga_x_lak)", by_name["mgga_x_lak"], 1, "maple/mgga_exc/mgga_x_lak.mpl:46",
                  "lak_fx (alpha = eps)", "alpha", EPS,
                  lambda a: {"rho": np.array([[0.1]]), "sigma": np.array([[sig_scan]]), "lapl": np.array([[0.0]]),
                             "tau": np.array([[unpol_tau_from_alpha(0.1, sig_scan, a)]])},
                  True, True, (1e-13, 4e-16)))

    # 48. CSK1 (mgga_k_csk1, a=1): csk_zI at z = -(-log(EPS))^(-1) (at p=0: z = 20*q/9 => q = 9*z/20)
    z_csk_cut = -((-math.log(EPS)) ** (-1.0))
    cases.append(("3_sec", "CSK1 (mgga_k_csk1)", by_name["mgga_k_csk1"], 1, "maple/mgga_exc/mgga_k_csk.mpl:40",
                  "csk_zI (z = -1/(-log(eps)))", "z_csk", z_csk_cut,
                  lambda zc: {"rho": np.array([[0.1]]), "sigma": np.array([[0.0]]),
                              "lapl": np.array([[unpol_lapl_from_u(0.1, (9.0 * zc / 20.0) / (X2S ** 2))]]), "tau": np.array([[0.1]])},
                  True, False, None))

    # 49. PBE_ERF_GWS (gga_x_pbe_erf_gws): pbe_x_erf_gws_num_smooth at mu_t = nu/2 = 1/16 and 1/2
    _fgws = lxc.Func(by_name["gga_x_pbe_erf_gws"], 1)
    _om_gws = ctypes.c_double()
    lxc.L.xc_hyb_cam_coef(ctypes.c_void_p(_fgws.p), ctypes.byref(_om_gws), ctypes.byref(_al), ctypes.byref(_be))
    _fgws.close()
    a_cnst_gws = ((4.0 / (9.0 * PI)) ** (1.0 / 3.0)) * float(_om_gws.value) / 2.0
    for acut, lbl in [(1.0 / 16.0, "gws_num_smooth (sr=1/16)"), (0.5, "gws_num_smooth (mid=1/2)")]:
        cases.append(("3_sec", "PBE-erf-GWS (gga_x_pbe_erf_gws)", by_name["gga_x_pbe_erf_gws"], 1,
                      "maple/gga_exc/gga_x_pbe_erf_gws.mpl:64", lbl, "mu_t", acut,
                      lambda a: {"rho": np.array([[rho_from_rs(a / a_cnst_gws)]]),
                                 "sigma": np.array([[unpol_sigma_from_s(rho_from_rs(a / a_cnst_gws), 0.3)]])},
                      True, False, None))

    # 50. MBRXC_BG (mgga_x_mbrxc_bg): mbrxc_v at x(Q) = sqrt(EPS) (requires negative tau via disable_fhc=True)
    x_mbr_cut = math.sqrt(EPS)
    rhs_mbr = ((1.0 + x_mbr_cut) ** (5.0 / 3.0)) * math.exp(-2.0 * x_mbr_cut / 3.0) / (x_mbr_cut - 3.0)
    Q_mbr_cut = ((32.0 * PI) ** (2.0 / 3.0)) / (6.0 * rhs_mbr)
    cases.append(("3_sec", "MBRXC-BG (mgga_x_mbrxc_bg)", by_name["mgga_x_mbrxc_bg"], 1,
                  "maple/mgga_exc/mgga_x_mbrxc_bg.mpl:56", "mbrxc_v (x = sqrt(eps))", "Q_mbr", Q_mbr_cut,
                  lambda qm: {"rho": np.array([[0.1]]), "sigma": np.array([[0.0]]), "lapl": np.array([[0.0]]),
                              "tau": np.array([[2.0 * ((qm + K_FACTOR_C) / (2.0 * 0.074746)) * ((0.05) ** (5.0 / 3.0))]])},
                  True, True, None))

    # 51. RPP (mgga_k_rpp): rpp_ap at x = x0 = 0.819411 and x = 0 (vary q via lapl at p = 0)
    # At p = 0: rpp_x(0, q) = 1 + 20*q/9 + rpp_bqq * q^2 * exp(-(q/rpp_c2)^2)
    rscan_fx = [1.0, -0.667, -0.4445555, -0.663086601049, 1.451297044490, -0.887998041597, 0.234528941479, -0.023185843322]
    rpp_dfx = sum(i * rscan_fx[7 - i + 1] for i in range(1, 8))  # wait: in Maple 1-indexed rscan_fx[8-i] for i=1..7 -> bisect directly on mgga_k_rpp!
    _frpp = lxc.Func(by_name["mgga_k_rpp"], 1)
    # At p = 0, zk(rho=0.1, lapl) = K_FACTOR_C * (0.05)^(5/3)*2/0.1 * rpp_ap(rpp_x(0, q)):
    _k0_rpp = float(_frpp.exc_vxc({"rho": np.array([[0.1]]), "sigma": np.array([[0.0]]), "lapl": np.array([[0.0]]), "tau": np.array([[0.1]])})["zk"][0])
    # For q = 0, rpp_x(0, 0) = 1.0 => rpp_ap(1.0) = 1.0 => _k0_rpp is exact Thomas-Fermi zk!
    def _bisect_rpp_q(target_ap):
        lo, hi = -0.8, 0.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            ap = float(_frpp.exc_vxc({"rho": np.array([[0.1]]), "sigma": np.array([[0.0]]),
                                      "lapl": np.array([[unpol_lapl_from_u(0.1, mid / (X2S**2))]]), "tau": np.array([[0.1]])})["zk"][0]) / _k0_rpp
            if ap < target_ap:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    q_rpp_x0 = _bisect_rpp_q(0.819411)
    q_rpp_0 = _bisect_rpp_q(1e-30)
    _frpp.close()
    for qcut, lbl in [(q_rpp_x0, "rpp_ap (x = x0 = 0.819411)"), (q_rpp_0, "rpp_ap (x = 0)")]:
        cases.append(("3_sec", "RPP (mgga_k_rpp)", by_name["mgga_k_rpp"], 1,
                      "maple/mgga_exc/mgga_k_rpp.mpl:65", lbl, "q_rpp", qcut,
                      lambda qq: {"rho": np.array([[0.1]]), "sigma": np.array([[0.0]]),
                                  "lapl": np.array([[unpol_lapl_from_u(0.1, qq / (X2S**2))]]), "tau": np.array([[0.1]])},
                      True, False, None))

    return cases


def main():
    cases = build_test_cases()
    rows = []
    hdr = f"{'Prio':7s} {'Functional':32s} {'Cut / Helper':34s} {'xc':>11s} {'zk(1e-9)':>10s} {'zk(1e-14)':>10s} {'der(1e-9)':>10s} {'der(1e-14)':>11s} {'Chan':9s}"
    print(hdr)
    print("-" * len(hdr))
    for prio, fname, fid, nspin, loc, helper, var, xc, make_pt, relax, no_fhc, edeltas in cases:
        res = eval_pair(fid, nspin, make_pt, xc, relax_thresholds=relax, disable_fhc=no_fhc, exact_deltas=edeltas)
        rw, rt = res["wide"], res["tight"]
        # Flag persistent jump if tight jump > 1e-9 and does not decay >= 100x from wide jump
        persistent = (rt["zk_rel_jump"] > 1e-9 and rt["zk_rel_jump"] > 0.1 * rw["zk_rel_jump"]) or (
            rt["der_rel_jump"] > 1e-8 and rt["der_rel_jump"] > 0.1 * rw["der_rel_jump"]
        )
        mark = " <-- PERSISTENT" if persistent else ""
        print(f"{prio:7s} {fname[:32]:32s} {helper[:34]:34s} {xc:11.4g} {rw['zk_rel_jump']:10.2e} {rt['zk_rel_jump']:10.2e} "
              f"{rw['der_rel_jump']:10.2e} {rt['der_rel_jump']:11.2e} {rt['der_chan']:9s}{mark}")
        rows.append({
            "prioridad": prio,
            "funcional": fname,
            "func_id": fid,
            "nspin": nspin,
            "archivo_linea": loc,
            "cambio_de_rama": helper,
            "variable": var,
            "x_c": xc,
            "delta_wide": rw["delta"],
            "delta_tight": rt["delta"],
            "zk_rel_jump_wide": rw["zk_rel_jump"],
            "zk_rel_jump_tight": rt["zk_rel_jump"],
            "zk_abs_jump_tight": rt["zk_abs_jump"],
            "der_canal_peor": rt["der_chan"],
            "der_rel_jump_wide": rw["der_rel_jump"],
            "der_rel_jump_tight": rt["der_rel_jump"],
            "der_abs_jump_tight": rt["der_abs_jump"],
            "zk_izq_tight": rt["zk_m"],
            "zk_der_tight": rt["zk_p"],
            "der_izq_tight": rt["der_m"],
            "der_der_tight": rt["der_p"],
            "salto_persistente": persistent,
        })

    out_csv = ROOT / "tests" / "results_branch_boundaries.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {len(rows)} boundary test results to {out_csv}")


if __name__ == "__main__":
    main()
