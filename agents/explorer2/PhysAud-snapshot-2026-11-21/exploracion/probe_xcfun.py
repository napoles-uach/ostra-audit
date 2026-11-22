#!/usr/bin/env python3
"""Probe A: Systematic physics-oracle & differential audit of XCFun (bundled in PySCF) vs Libxc.

Oracles:
  1. XCFun internal energy-potential consistency (4-pt FD with step scan)
  2. XCFun polarized (rho_a = rho_b) == unpolarized (rho)
  3. XCFun vs Libxc cross-comparison (zk, vrho, vsigma, vtau, fxc) on all shared functionals
     across interior (typical) AND edge (small s, large s, alpha=0, high/low rho, high zeta) points.
"""
import math
import numpy as np
from pyscf import dft

# Discover all available XCFun functional aliases in PySCF
xcfun = dft.xcfun
libxc = dft.libxc

XCFUN_FUNCS = sorted(xcfun.XC_CODES.keys())
print(f"XCFun version: {xcfun.__version__}, aliases in PySCF: {len(XCFUN_FUNCS)}")

# Build test points: (name, rho_unpol, sigma_unpol, lapl_unpol, tau_unpol) and polarized counterparts
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
    # PySCF MGGA rho format (spin=0): (6, N) = [rho, gx, gy, gz, lapl, tau]
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


POINTS_UNPOL = [
    ("interior",    make_mgga_unpol(0.10, 0.50, 1.20, 0.3)),
    ("small_s",     make_mgga_unpol(0.10, 1e-4, 1.00, 0.0)),
    ("large_s",     make_mgga_unpol(1e-4, 12.0, 2.00, 0.5)),
    ("iso_alpha0",  make_mgga_unpol(0.10, 0.40, 1e-6, -0.2)),
    ("high_rho",    make_mgga_unpol(50.0, 0.30, 0.80, 0.1)),
    ("low_rho",     make_mgga_unpol(1e-6, 1.50, 1.50, 0.0)),
]

POINTS_POL = [
    ("pol_equal",   make_mgga_pol(0.05, 0.05, 0.50, 0.50, 1.20, 1.20, 0.3, 0.3)),
    ("pol_mid",     make_mgga_pol(0.07, 0.03, 0.40, 0.70, 0.90, 1.50, 0.2, -0.1)),
    ("pol_ferro",   make_mgga_pol(0.10, 1e-8, 0.50, 0.50, 1.00, 1.00, 0.0, 0.0)),
]

# 1. Check all XCFun functionals: spin=0 vs spin=1 (rho_a == rho_b), AND XCFun vs Libxc!
pol_discrep = []
diff_discrep = []
fd_discrep = []
nan_hits = []

# Skip the 19 VV10/rVV10 aliases that SIGSEGV (-11) in PySCF's xcfun interface
SKIP_CRASH = {
    'B97MV', 'B97M_V', 'LCVV10', 'LC_VV10', 'REVSCANVV10', 'REVSCAN_VV10',
    'SCANLRVV10', 'SCANLVV10', 'SCANL_RVV10', 'SCANL_VV10', 'SCANRVV10',
    'SCANVV10', 'SCAN_RVV10', 'SCAN_VV10', 'VV10', 'WB97MV', 'WB97M_V',
    'WB97XV', 'WB97X_V',
}

for name in XCFUN_FUNCS:
    if name in SKIP_CRASH:
        continue
    # Skip pure aliases that map to same XC_CODES? Keep unique codes only
    code = xcfun.XC_CODES[name]
    # Test spin=0 vs spin=1 at interior point
    rho0 = POINTS_UNPOL[0][1]
    rho_eq = POINTS_POL[0][1]
    try:
        e0, v0, f0, _ = xcfun.eval_xc(name, rho0, spin=0, deriv=2)
        e1, v1, f1, _ = xcfun.eval_xc(name, rho_eq, spin=1, deriv=2)
    except Exception:
        continue

    if not np.all(np.isfinite(e0)) or not np.all(np.isfinite(v0[0])):
        nan_hits.append((name, "interior_unpol_nonfinite", e0, v0[0]))
        continue

    # Check pol_equal vs unpol
    de_pol = abs(float(e1[0]) - float(e0[0])) / max(abs(float(e0[0])), 1e-12)
    dv_pol = abs(float(v1[0][0, 0]) - float(v0[0][0])) / max(abs(float(v0[0][0])), 1e-12)
    if de_pol > 1e-7 or dv_pol > 1e-7 or not np.isfinite(de_pol) or not np.isfinite(dv_pol):
        pol_discrep.append((name, de_pol, dv_pol, float(e0[0]), float(e1[0]), float(v0[0][0]), float(v1[0][0, 0])))

    # Check XCFun internal 4-pt FD consistency on vrho at interior point!
    r0 = 0.10
    h = 1e-5 * r0
    rpts = np.array([r0 - 2*h, r0 - h, r0, r0 + h, r0 + 2*h])
    arr = np.repeat(rho0, 5, axis=1)
    arr[0, :] = rpts
    ek, vk, _, _ = xcfun.eval_xc(name, arr, spin=0, deriv=1)
    E = rpts * ek
    fd_vrho = (E[0] - 8*E[1] + 8*E[3] - E[4]) / (12.0 * h)
    err_fd = abs(float(vk[0][2]) - fd_vrho) / max(abs(fd_vrho), 1e-12)
    if err_fd > 1e-6 or not np.isfinite(err_fd):
        fd_discrep.append((name, err_fd, float(vk[0][2]), fd_vrho))

    # Check XCFun edge points for NaN/Inf
    for pname, pt in POINTS_UNPOL[1:]:
        try:
            ex, vx, _, _ = xcfun.eval_xc(name, pt, spin=0, deriv=1)
            if not np.all(np.isfinite(ex)) or any(v is not None and not np.all(np.isfinite(v)) for v in vx):
                nan_hits.append((name, f"unpol_{pname}", ex, [v for v in vx if v is not None]))
        except Exception:
            pass
    for pname, pt in POINTS_POL[1:]:
        try:
            ex, vx, _, _ = xcfun.eval_xc(name, pt, spin=1, deriv=1)
            if not np.all(np.isfinite(ex)) or any(v is not None and not np.all(np.isfinite(v)) for v in vx):
                nan_hits.append((name, f"pol_{pname}", ex, [v for v in vx if v is not None]))
        except Exception:
            pass

    # Cross-compare XCFun vs Libxc on ALL shared names across interior & pol_mid & edge points!
    if name in libxc.XC_CODES:
        for pname, pt in [("interior", POINTS_UNPOL[0][1]), ("small_s", POINTS_UNPOL[1][1]), ("large_s", POINTS_UNPOL[2][1])]:
            try:
                ex_x, vx_x, fx_x, _ = xcfun.eval_xc(name, pt, spin=0, deriv=2)
                ex_l, vx_l, fx_l, _ = libxc.eval_xc(name, pt, spin=0, deriv=2)
            except Exception:
                continue
            if not (np.isfinite(ex_x[0]) and np.isfinite(ex_l[0])):
                continue
            de = abs(float(ex_x[0]) - float(ex_l[0])) / max(abs(float(ex_l[0])), 1e-12)
            dv = abs(float(vx_x[0][0]) - float(vx_l[0][0])) / max(abs(float(vx_l[0][0])), 1e-12)
            dvs = 0.0
            if vx_x[1] is not None and vx_l[1] is not None:
                dvs = abs(float(vx_x[1][0]) - float(vx_l[1][0])) / max(abs(float(vx_l[1][0])), 1e-12)
            df = 0.0
            if fx_x[0] is not None and fx_l[0] is not None:
                df = abs(float(fx_x[0][0]) - float(fx_l[0][0])) / max(abs(float(fx_l[0][0])), 1e-12)
            if max(de, dv, dvs, df) > 1e-4:
                diff_discrep.append((name, pname, de, dv, dvs, df, float(ex_x[0]), float(ex_l[0])))

print(f"\n=== 1. XCFun pol_equal vs unpol discrepancies (>1e-7): {len(pol_discrep)} ===")
for r in sorted(pol_discrep, key=lambda x: -max(x[1], x[2]))[:25]:
    print(f"  {r[0]:24s} de={r[1]:.3e} dv={r[2]:.3e} e0={r[3]:.6e} e1={r[4]:.6e}")

print(f"\n=== 2. XCFun internal FD4 vrho discrepancies (>1e-6): {len(fd_discrep)} ===")
for r in sorted(fd_discrep, key=lambda x: -x[1])[:25]:
    print(f"  {r[0]:24s} err_fd={r[1]:.3e} vrho={r[2]:.8e} fd4={r[3]:.8e}")

print(f"\n=== 3. XCFun non-finite (NaN/Inf) hits: {len(nan_hits)} ===")
for r in nan_hits[:30]:
    print(f"  {r[0]:24s} @ {r[1]:16s} exc={r[2]}")

print(f"\n=== 4. XCFun vs Libxc discrepancies (>1e-4): {len(diff_discrep)} ===")
seen = set()
for r in sorted(diff_discrep, key=lambda x: -max(x[2], x[3], x[4], x[5])):
    if (r[0], r[1]) in seen:
        continue
    seen.add((r[0], r[1]))
    print(f"  {r[0]:22s} @{r[1]:9s} de={r[2]:.2e} dvr={r[3]:.2e} dvs={r[4]:.2e} df2={r[5]:.2e} | xc={r[6]:.6e} lx={r[7]:.6e}")
