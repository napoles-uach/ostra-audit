#!/usr/bin/env python3
"""Boundary continuity tests for all C-level special functions in libxc (src/*.c, src/*.h).

Measures value f(x) and first derivative f'(x) (exact jet where available, and centered FD
away from the boundary on each branch) at x_c +- 1e-9 and x_c +- 1e-14 (or scaled by |x_c|
when |x_c| != O(1)).
"""
import csv
import ctypes
import math
import os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DYLIB = os.environ.get("LIBXC_PATH", str(ROOT / "build-libxc" / "libxc.dylib"))
L = ctypes.CDLL(DYLIB)

# Bind C special functions exported in libxc.dylib
DBL = ctypes.c_double
INT = ctypes.c_int

for fn in [
    "xc_x_minus_log1p", "xc_atan_tail3", "xc_log1p_tail3",
    "xc_bessel_I0_scaled", "xc_bessel_I0",
    "xc_bessel_I1_scaled", "xc_bessel_I1",
    "xc_bessel_K0_scaled", "xc_bessel_K0",
    "xc_bessel_K1_scaled", "xc_bessel_K1",
    "xc_erfcx", "xc_erfcx_d1",
    "xc_E1_scaled_d1",
    "LambertW", "xc_lambertw_d1",
    "xc_dilogarithm",
    "xc_mgga_x_br89_get_x", "xc_mgga_x_mbrxc_get_x",
]:
    f = getattr(L, fn)
    f.argtypes = [DBL]
    f.restype = DBL

L.xc_exp_tail.argtypes = [INT, DBL]
L.xc_exp_tail.restype = DBL
L.xc_erf_tail.argtypes = [INT, DBL]
L.xc_erf_tail.restype = DBL
L.xc_expint_e1_impl.argtypes = [DBL, INT]
L.xc_expint_e1_impl.restype = DBL


def rel_jump(a, b, ref_scale=0.0):
    denom = max(abs(a), abs(b), abs(ref_scale), 1e-300)
    return abs(b - a) / denom


def one_sided_deriv(f, x, side, xc):
    """4-point one-sided finite difference strictly on one side of x_c (never crossing x_c or 0)."""
    scale = max(abs(xc), 1.0) if xc != 0.0 else 1.0
    h = 1e-6 * scale
    if xc != 0.0:
        # Never step across 0 (important when xc ~ 2*sqrt(eps) and f uses fabs(x))
        h = min(h, 0.05 * abs(xc))
    s = 1.0 if side > 0 else -1.0
    f0 = f(x)
    f1 = f(x + s * h)
    f2 = f(x + 2.0 * s * h)
    f3 = f(x + 3.0 * s * h)
    # 3rd-order forward/backward difference: (-11 f0 + 18 f1 - 9 f2 + 2 f3) / (6 s h)
    return (-11.0 * f0 + 18.0 * f1 - 9.0 * f2 + 2.0 * f3) / (6.0 * s * h)


def probe_c_cut(name, file_line, f, xc, df=None):
    # Use exact 1e-9 and 1e-14 when 1e-3 <= |xc| <= 20, else proportional to |xc|
    if 1e-3 <= abs(xc) <= 20.0:
        d_wide = 1e-9
        d_tight = 1e-14
    elif xc == 0.0:
        d_wide = 1e-9
        d_tight = 1e-14
    else:
        d_wide = abs(xc) * 1e-9
        d_tight = max(abs(xc) * 1e-14, abs(np.nextafter(xc, np.inf) - xc) * 4.0)

    results = {}
    for tag, d in [("1e-9", d_wide), ("1e-14", d_tight)]:
        xm = xc - d
        xp = xc + d
        if xm == xc:
            xm = float(np.nextafter(xc, -np.inf))
        if xp == xc:
            xp = float(np.nextafter(xc, np.inf))
        vm = f(xm)
        vp = f(xp)
        if df is not None:
            dm = df(xm)
            dp = df(xp)
        else:
            dm = one_sided_deriv(f, xm, -1, xc)
            dp = one_sided_deriv(f, xp, +1, xc)
        # Avoid zero-crossing metric division artifact when f(xc) = 0 and f'(xc) != 0 (e.g. Li2(0) = 0, W(0) = 0)
        val_scale = max(abs(dm), abs(dp)) if min(abs(vm), abs(vp)) < 1e-6 * max(abs(dm), abs(dp)) else 0.0
        der_scale = max(abs(vm), abs(vp), 1.0) if max(abs(dm), abs(dp)) < 1e-6 else 0.0
        results[tag] = {
            "xm": xm, "xp": xp,
            "vm": vm, "vp": vp, "val_rel_jump": rel_jump(vm, vp, val_scale), "val_abs_jump": abs(vp - vm),
            "dm": dm, "dp": dp, "der_rel_jump": rel_jump(dm, dp, der_scale), "der_abs_jump": abs(dp - dm),
        }
    return {
        "name": name,
        "file_line": file_line,
        "xc": xc,
        **results,
    }


def main():
    eps = float(np.finfo(np.float64).eps)
    sq_eps = math.sqrt(eps)
    cbrt_eps = eps ** (1.0 / 3.0)

    cases = [
        ("xc_x_minus_log1p(+0.25)", "src/util_inline.h:15", L.xc_x_minus_log1p, 0.25, lambda x: x / (1.0 + x)),
        ("xc_x_minus_log1p(-0.25)", "src/util_inline.h:15", L.xc_x_minus_log1p, -0.25, lambda x: x / (1.0 + x)),
        ("xc_atan_tail3(+0.6)", "src/util_inline.h:37", L.xc_atan_tail3, 0.6, lambda x: x**4 / (1.0 + x*x)),
        ("xc_log1p_tail3(+0.5)", "src/util_inline.h:54", L.xc_log1p_tail3, 0.5, lambda x: -(x**3) / (1.0 + x)),
        ("xc_exp_tail(n=3, x=+4)", "src/util_inline.h:77", lambda x: L.xc_exp_tail(3, x), 4.0, lambda x: L.xc_exp_tail(2, x)),
        ("xc_exp_tail(n=3, x=-4)", "src/util_inline.h:77", lambda x: L.xc_exp_tail(3, x), -4.0, lambda x: L.xc_exp_tail(2, x)),
        ("xc_erf_tail(n=3, z=+2)", "src/util_inline.h:102", lambda z: L.xc_erf_tail(3, z), 2.0, lambda z: (2.0 / math.sqrt(math.pi)) * L.xc_exp_tail(3, -z*z)),
        ("xc_bessel_I0_scaled(2*sqrt(eps))", "src/bessel.c:54", L.xc_bessel_I0_scaled, 2.0 * sq_eps, None),
        ("xc_bessel_I0_scaled(3.0)", "src/bessel.c:56", L.xc_bessel_I0_scaled, 3.0, None),
        ("xc_bessel_I0_scaled(8.0)", "src/bessel.c:58", L.xc_bessel_I0_scaled, 8.0, None),
        ("xc_bessel_I0(2*sqrt(eps))", "src/bessel.c:78", L.xc_bessel_I0, 2.0 * sq_eps, None),
        ("xc_bessel_I0(3.0)", "src/bessel.c:80", L.xc_bessel_I0, 3.0, None),
        ("xc_bessel_I1_scaled(2*sqrt(2*eps))", "src/bessel.c:133", L.xc_bessel_I1_scaled, 2.0 * math.sqrt(2.0) * sq_eps, None),
        ("xc_bessel_I1_scaled(3.0)", "src/bessel.c:135", L.xc_bessel_I1_scaled, 3.0, None),
        ("xc_bessel_I1_scaled(8.0)", "src/bessel.c:138", L.xc_bessel_I1_scaled, 8.0, None),
        ("xc_bessel_I1(2*sqrt(2*eps))", "src/bessel.c:163", L.xc_bessel_I1, 2.0 * math.sqrt(2.0) * sq_eps, None),
        ("xc_bessel_I1(3.0)", "src/bessel.c:165", L.xc_bessel_I1, 3.0, None),
        ("xc_bessel_K0_scaled(2.0)", "src/bessel.c:205", L.xc_bessel_K0_scaled, 2.0, None),
        ("xc_bessel_K0_scaled(8.0)", "src/bessel.c:207", L.xc_bessel_K0_scaled, 8.0, None),
        ("xc_bessel_K0(2.0)", "src/bessel.c:224", L.xc_bessel_K0, 2.0, lambda x: -L.xc_bessel_K1(x)),
        ("xc_bessel_K1_scaled(2.0)", "src/bessel.c:263", L.xc_bessel_K1_scaled, 2.0, None),
        ("xc_bessel_K1_scaled(8.0)", "src/bessel.c:266", L.xc_bessel_K1_scaled, 8.0, None),
        ("xc_bessel_K1(2.0)", "src/bessel.c:289", L.xc_bessel_K1, 2.0, None),
        ("xc_expint_e1_scaled(-10)", "src/expint_e1.c:100", lambda x: L.xc_expint_e1_impl(x, 1), -10.0, None),
        ("xc_expint_e1_scaled(-4)", "src/expint_e1.c:103", lambda x: L.xc_expint_e1_impl(x, 1), -4.0, None),
        ("xc_expint_e1_scaled(-1)", "src/expint_e1.c:106", lambda x: L.xc_expint_e1_impl(x, 1), -1.0, None),
        ("xc_expint_e1_scaled(+1)", "src/expint_e1.c:113", lambda x: L.xc_expint_e1_impl(x, 1), 1.0, L.xc_E1_scaled_d1),
        ("xc_expint_e1_scaled(+4)", "src/expint_e1.c:116", lambda x: L.xc_expint_e1_impl(x, 1), 4.0, L.xc_E1_scaled_d1),
        ("xc_E1_scaled_jet(+18)", "src/expint_e1.c:140", lambda x: L.xc_expint_e1_impl(x, 1), 18.0, L.xc_E1_scaled_d1),
        ("xc_erfcx(-6.1)", "src/faddeeva.c:503", L.xc_erfcx, -6.1, L.xc_erfcx_d1),
        ("xc_erfcx(0.0)", "src/faddeeva.c:491", L.xc_erfcx, 0.0, L.xc_erfcx_d1),
        ("xc_erfcx(bin_50: x=4.0)", "src/faddeeva.c:81", L.xc_erfcx, 4.0, L.xc_erfcx_d1),
        ("xc_erfcx_jet(8.0)", "src/faddeeva.c:518", L.xc_erfcx, 8.0, L.xc_erfcx_d1),
        ("xc_erfcx(50.0)", "src/faddeeva.c:492", L.xc_erfcx, 50.0, L.xc_erfcx_d1),
        ("xc_erfcx(5e7)", "src/faddeeva.c:494", L.xc_erfcx, 5e7, L.xc_erfcx_d1),
        ("LambertW(+cbrt(eps))", "src/special_functions.c:45", L.LambertW, cbrt_eps, L.xc_lambertw_d1),
        ("LambertW(-cbrt(eps))", "src/special_functions.c:45", L.LambertW, -cbrt_eps, L.xc_lambertw_d1),
        ("LambertW(-0.31408624)", "src/special_functions.c:49", L.LambertW, -0.3140862435046707, L.xc_lambertw_d1),
        ("LambertW(+1.14987649)", "src/special_functions.c:53", L.LambertW, 1.149876485041417, L.xc_lambertw_d1),
        ("xc_dilogarithm(-1.0)", "src/special_functions.c:220", L.xc_dilogarithm, -1.0, None),
        ("xc_dilogarithm(0.0)", "src/special_functions.c:217", L.xc_dilogarithm, 0.0, None),
        ("xc_dilogarithm(0.5)", "src/special_functions.c:204", L.xc_dilogarithm, 0.5, None),
        ("xc_dilogarithm(1.0)", "src/special_functions.c:199", L.xc_dilogarithm, 1.0, None),
        ("xc_dilogarithm(2.0)", "src/special_functions.c:193", L.xc_dilogarithm, 2.0, None),
        ("br89_get_x(+2*eps)", "src/mgga_x_br89_get_x_inline.h:32", L.xc_mgga_x_br89_get_x, 2.0 * eps, None),
        ("mbrxc_get_x(+2*eps)", "src/mgga_x_mbrxc_get_x_inline.h:24", L.xc_mgga_x_mbrxc_get_x, 2.0 * eps, None),
    ]

    rows = []
    print(f"{'Function & Cut':36s} {'Loc':32s} {'val_jump(1e-9)':>14s} {'val_jump(1e-14)':>15s} {'der_jump(1e-9)':>14s} {'der_jump(1e-14)':>15s}")
    print("-" * 132)
    for name, loc, fn, xc, dfn in cases:
        res = probe_c_cut(name, loc, fn, xc, dfn)
        r9 = res["1e-9"]
        r14 = res["1e-14"]
        print(f"{name:36s} {loc:32s} {r9['val_rel_jump']:14.3e} {r14['val_rel_jump']:15.3e} {r9['der_rel_jump']:14.3e} {r14['der_rel_jump']:15.3e}")
        rows.append({
            "helper": name,
            "file_line": loc,
            "xc": xc,
            "val_rel_jump_1e9": r9["val_rel_jump"],
            "val_rel_jump_1e14": r14["val_rel_jump"],
            "der_rel_jump_1e9": r9["der_rel_jump"],
            "der_rel_jump_1e14": r14["der_rel_jump"],
            "vm_1e14": r14["vm"],
            "vp_1e14": r14["vp"],
            "dm_1e14": r14["dm"],
            "dp_1e14": r14["dp"],
        })

    out_csv = ROOT / "tests" / "results_c_helpers_boundaries.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {out_csv}")


if __name__ == "__main__":
    main()
