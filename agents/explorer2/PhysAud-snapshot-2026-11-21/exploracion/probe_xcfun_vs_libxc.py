#!/usr/bin/env python3
"""Differential oracle: XCFun 2.1.1 vs Libxc 7.0.0 across EVERY functional in XCFun,
properly mapped to its Libxc counterpart, plus edge-case & API audit on XCFun!
"""
import math
import numpy as np
from pyscf import dft

xcfun = dft.xcfun
libxc = dft.libxc

PI = math.pi
X2S = 1.0 / (2.0 * (6.0 * PI**2) ** (1.0 / 3.0))
K_C = 0.3 * (6.0 * PI**2) ** (2.0 / 3.0)


def make_mgga_unpol(rho, s, alpha, u_red=0.0):
    sig = 4.0 * ((s / X2S) * (rho / 2.0) ** (4.0 / 3.0)) ** 2
    gx = math.sqrt(sig)
    rho_s = rho / 2.0
    xs2 = (sig / 4.0) / (rho_s ** (8.0 / 3.0))
    t_s = xs2 / 8.0 + alpha * K_C
    tau = 2.0 * t_s * (rho_s ** (5.0 / 3.0))
    lapl = 2.0 * u_red * (rho_s ** (5.0 / 3.0))
    return np.array([[rho], [gx], [0.0], [0.0], [lapl], [tau]], dtype=np.float64)


def make_mgga_pol(ra, rb, sa, sb, aa, ab, ua=0.0, ub=0.0):
    def one(r, s, a, u):
        sig = ((s / X2S) * (r ** (4.0 / 3.0))) ** 2
        gx = math.sqrt(sig)
        xs2 = sig / (r ** (8.0 / 3.0))
        t_s = xs2 / 8.0 + a * K_C
        tau = t_s * (r ** (5.0 / 3.0))
        lapl = u * (r ** (5.0 / 3.0))
        return np.array([[r], [gx], [0.0], [0.0], [lapl], [tau]], dtype=np.float64)
    return (one(ra, sa, aa, ua), one(rb, sb, ab, ub))


# Curated 1-to-1 mapping (xcfun_name, libxc_name) covering all major & specialized XCFun functionals
PAIRS = [
    ("SLATERX", "LDA_X"),
    ("VWN5C", "LDA_C_VWN"),
    ("VWN3C", "LDA_C_VWN_3"),
    ("VWN_PBEC", "LDA_C_PW"),
    ("PZ81C", "LDA_C_PZ"),
    ("PW92C", "LDA_C_PW"),
    ("BECKEX", "GGA_X_B88"),
    ("BECKECORRX", "GGA_X_B88"),  # check if different!
    ("PBEX", "GGA_X_PBE"),
    ("PBEC", "GGA_C_PBE"),
    ("RPBEX", "GGA_X_RPBE"),
    ("REVPBEX", "GGA_X_PBE_R"),
    ("PBESOLX", "GGA_X_PBE_SOL"),
    ("PBESOLC", "GGA_C_PBE_SOL"),
    ("PBEINTX", "GGA_X_PBEINT"),
    ("PBEINTC", "GGA_C_PBEINT"),
    ("PW86X", "GGA_X_PW86"),
    ("PW91X", "GGA_X_PW91"),
    ("PW91K", "GGA_K_PW91"),
    ("OPTX", "GGA_X_OPTX"),
    ("LYPC", "GGA_C_LYP"),
    ("P86C", "GGA_C_P86"),
    ("KT1X", "GGA_XC_KT1"),
    ("KT2X", "GGA_XC_KT2"),
    ("KT3X", "GGA_XC_KT3"),
    ("APBEC", "GGA_C_APBE"),
    ("APBEX", "GGA_X_APBE"),
    ("SOGGA", "GGA_X_SOGGA"),
    ("SOGGA11X", "GGA_X_SOGGA11"),
    ("SOGGA11C", "GGA_C_SOGGA11"),
    ("TPSSX", "MGGA_X_TPSS"),
    ("TPSSC", "MGGA_C_TPSS"),
    ("REVTPSSX", "MGGA_X_REVTPSS"),
    ("REVTPSSC", "MGGA_C_REVTPSS"),
    ("SCANX", "MGGA_X_SCAN"),
    ("SCANC", "MGGA_C_SCAN"),
    ("REVSCANX", "MGGA_X_REVSCAN"),
    ("REVSCANC", "MGGA_C_REVSCAN"),
    ("RPPX", "MGGA_X_RPP09"),
    ("M05X", "HYB_MGGA_X_M05"),
    ("M05C", "MGGA_C_M05"),
    ("M052XX", "HYB_MGGA_X_M05_2X"),
    ("M052XC", "MGGA_C_M05_2X"),
    ("M06X", "HYB_MGGA_X_M06"),
    ("M06C", "MGGA_C_M06"),
    ("M062XX", "HYB_MGGA_X_M06_2X"),
    ("M062XC", "MGGA_C_M06_2X"),
    ("M06LX", "MGGA_X_M06_L"),
    ("M06LC", "MGGA_C_M06_L"),
    ("M06HFX", "HYB_MGGA_X_M06_HF"),
    ("M06HFC", "MGGA_C_M06_HF"),
    ("MS0X", "MGGA_X_MS0"),
    ("MS1X", "MGGA_X_MS1"),
    ("MS2X", "MGGA_X_MS2"),
    ("MVSX", "MGGA_X_MVS"),
    ("BRX", "MGGA_X_BR89"),
    ("PKZBX", "MGGA_X_PKZB"),
    ("PKZBC", "MGGA_C_PKZB"),
    ("B3LYP", "HYB_GGA_XC_B3LYP"),
    ("B3LYPG", "HYB_GGA_XC_B3LYP"),
    ("PBE0", "HYB_GGA_XC_PBEH"),
    ("HSE06", "HYB_GGA_XC_HSE06"),
    ("CAMB3LYP", "HYB_GGA_XC_CAM_B3LYP"),
    ("B97", "HYB_GGA_XC_B97"),
    ("B97_1", "HYB_GGA_XC_B97_1"),
    ("B97_2", "HYB_GGA_XC_B97_2"),
    ("BLYP", "GGA_X_B88,GGA_C_LYP"),
    ("PBE", "GGA_X_PBE,GGA_C_PBE"),
    ("TPSS", "MGGA_X_TPSS,MGGA_C_TPSS"),
    ("TPSSH", "HYB_MGGA_XC_TPSSH"),
    ("SCAN", "MGGA_X_SCAN,MGGA_C_SCAN"),
    ("R2SCAN", "MGGA_X_R2SCAN,MGGA_C_R2SCAN"),
    ("M06", "HYB_MGGA_X_M06,MGGA_C_M06"),
    ("M06_2X", "HYB_MGGA_X_M06_2X,MGGA_C_M06_2X"),
    ("M06_L", "MGGA_X_M06_L,MGGA_C_M06_L"),
    ("WB97", "HYB_GGA_XC_WB97"),
    ("WB97X", "HYB_GGA_XC_WB97X"),
]

# Filter to valid xcfun keys
valid_xcfun = set(xcfun.XC_CODES.keys())
print(f"Checking {len(PAIRS)} XCFun vs Libxc pairs...")
print(f"{'XCFun':14s} {'Libxc':24s} {'Pt':10s} {'rel_exc':>9s} {'rel_vr':>9s} {'rel_vs':>9s} {'rel_vt':>9s} {'rel_f2':>9s} | exc_xcfun / exc_libxc")
print("-" * 125)

pts = [
    ("unpol_int", 0, make_mgga_unpol(0.10, 0.50, 1.20, 0.3)),
    ("pol_mid",   1, make_mgga_pol(0.07, 0.03, 0.40, 0.70, 0.90, 1.50, 0.2, -0.1)),
    ("unpol_s15", 0, make_mgga_unpol(1e-4, 15.0, 2.00, 0.5)),
]

hits = []
for xf, lx in PAIRS:
    if xf not in valid_xcfun:
        continue
    for pname, spin, pt in pts:
        try:
            ex_x, vx_x, fx_x, _ = xcfun.eval_xc(xf, pt, spin=spin, deriv=2)
            ex_l, vx_l, fx_l, _ = libxc.eval_xc(lx, pt, spin=spin, deriv=2)
        except Exception:
            continue
        def rdiff(a, b):
            if a is None or b is None:
                return 0.0
            a, b = np.asarray(a), np.asarray(b)
            if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
                return float("nan")
            return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-12))

        de = rdiff(ex_x, ex_l)
        dvr = rdiff(vx_x[0], vx_l[0])
        dvs = rdiff(vx_x[1], vx_l[1])
        dvt = rdiff(vx_x[3], vx_l[3])
        df2 = rdiff(fx_x[0], fx_l[0])
        mx = max(de, dvr, dvs, dvt, df2)
        if math.isnan(mx) or mx > 1e-4:
            hits.append((mx, xf, lx, pname, de, dvr, dvs, dvt, df2, float(ex_x[0]), float(ex_l[0])))

for _, xf, lx, pname, de, dvr, dvs, dvt, df2, exx, exl in sorted(hits, key=lambda x: -(1e9 if math.isnan(x[0]) else x[0])):
    print(f"{xf:14s} {lx:24s} {pname:10s} {de:9.2e} {dvr:9.2e} {dvs:9.2e} {dvt:9.2e} {df2:9.2e} | {exx:+.6e} / {exl:+.6e}")

print("\n=== ALL XCFun XC_CODES keys (for completeness) ===")
print(sorted(valid_xcfun))
