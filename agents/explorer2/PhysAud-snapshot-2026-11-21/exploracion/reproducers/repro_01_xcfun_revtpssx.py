#!/usr/bin/env python3
"""Minimal reproducer 01: XCFun 2.1.1 REVTPSSX (id=44) triple bug.
1) Registered as XC_GGA (type=1) instead of XC_MGGA (type=2) -> NaN when called alone.
2) Missing square in denominator: (1 + sqrt(e)*p) instead of (1 + sqrt(e)*p)^2.
3) Dimensionful (p * sigma)^2 instead of p^2 inside sqrt(0.5*(9/25*z^2 + p^2)),
   violating uniform coordinate scaling, creating a pole at x = -kappa = -0.804,
   and causing -3.786 Eh (-2,376 kcal/mol) error in E_x on H2O (SCF diverges to +13.84 Eh).
"""
import numpy as np
from pyscf import gto, dft
from pyscf.dft import xcfun, libxc

def make_pt(rho, s, alpha):
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
    tau = sig / (8 * rho) + alpha * tau_unif
    return np.array([[rho], [np.sqrt(sig)], [0.0], [0.0], [0.0], [tau]])

print("=== 1. Registration bug ===")
print(f"XCFun xc_type('REVTPSSX') = {xcfun.xc_type('REVTPSSX')} (expected MGGA)")
pt_gga = make_pt(0.1, 0.5, 1.0)[:4]
ex_alone, v_alone, _, _ = xcfun.eval_xc("REVTPSSX", pt_gga, spin=0, deriv=1)
print(f"eval_xc('REVTPSSX') alone: exc={ex_alone[0]}, vrho={v_alone[0][0]}, vsigma={v_alone[1][0]}")

print("\n=== 2. Exact reverse-engineered formula proof (15-digit match) ===")
ex_lda = lambda r: -0.75 * (3 / np.pi)**(1/3) * r**(1/3)
kap, b, c, e, mu, mu_ge = 0.804, 0.40, 2.35204, 2.1677, 0.14, 10/81
max_err = 0.0
for rho, s, alpha in [(0.1, 0.5, 1.0), (1.0, 0.5, 1.0), (3.0, 0.5, 1.0), (3.5, 0.5, 1.0)]:
    pt = make_pt(rho, s, alpha)
    # Pairing with REVTPSSC forces MGGA evaluation; subtract REVTPSSC to isolate REVTPSSX
    ex_obs = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0] - xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ex_ref = libxc.eval_xc("MGGA_X_REVTPSS", pt, spin=0, deriv=0)[0][0]
    p = s**2
    z = (5/3 * p) / (5/3 * p + alpha)
    qb = 9/20 * (alpha - 1) / np.sqrt(1 + b * alpha * (alpha - 1)) + 2*p/3
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    num = (
        (mu_ge + c * z**3 / (1 + z**2)**2) * p
        + 146/2025 * qb**2
        - 73/405 * qb * np.sqrt(0.5 * (9/25 * z**2 + (p * sig)**2))
        + mu_ge**2 / kap * p**2
        + 2 * np.sqrt(e) * mu_ge * 9/25 * z**2
        + e * mu * p**3
    )
    x_bug = num / (1.0 + np.sqrt(e) * p)
    ex_pred = ex_lda(rho) * (1.0 + kap - kap / (1.0 + x_bug / kap))
    max_err = max(max_err, abs(ex_obs - ex_pred))
    print(f"  rho={rho:3.1f} s={s} a={alpha} | Fx_xcfun={ex_obs/ex_lda(rho):+.6f} (Libxc={ex_ref/ex_lda(rho):+.6f}) | |obs-pred|={abs(ex_obs-ex_pred):.2e}")
print(f"Max |XCFun - buggy_formula| = {max_err:.3e}")

print("\n=== 3. Impact on H2O (def2-TZVP, grid level 4) ===")
mol = gto.M(atom="O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587", basis="def2-tzvp", verbose=0)
mf_l = dft.RKS(mol); mf_l.xc = "MGGA_X_REVTPSS,MGGA_C_REVTPSS"; mf_l.grids.level = 4; mf_l.conv_tol = 1e-10
e_l = mf_l.kernel()
dm = mf_l.make_rdm1()
ao = dft.numint.eval_ao(mol, mf_l.grids.coords, deriv=2)
rho_g = dft.numint.eval_rho(mol, ao, dm, xctype="MGGA", with_lapl=True)
wt = mf_l.grids.weights
ex_xc = xcfun.eval_xc("REVTPSSX,REVTPSSC", rho_g, spin=0, deriv=0)[0]
ex_c  = xcfun.eval_xc("REVTPSSC", rho_g, spin=0, deriv=0)[0]
el_x  = libxc.eval_xc("MGGA_X_REVTPSS", rho_g, spin=0, deriv=0)[0]
Ex_xcfun = float(np.sum(wt * rho_g[0] * (ex_xc - ex_c)))
Ex_libxc = float(np.sum(wt * rho_g[0] * el_x))
print(f"  E_x(XCFun) = {Ex_xcfun:.6f} Eh, E_x(Libxc) = {Ex_libxc:.6f} Eh, diff = {Ex_xcfun - Ex_libxc:+.6f} Eh ({(Ex_xcfun - Ex_libxc)*627.5095:+.1f} kcal/mol)")
mf_x = dft.RKS(mol); mf_x._numint.libxc = dft.xcfun; mf_x.xc = "REVTPSSX,REVTPSSC"; mf_x.grids.level = 4; mf_x.max_cycle = 15
e_x = mf_x.kernel()
print(f"  SCF H2O: E(Libxc) = {e_l:.6f} Eh, E(XCFun, 15 cyc) = {e_x:.6f} Eh (converged={mf_x.converged})")
