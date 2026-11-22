#!/usr/bin/env python3
"""Minimal reproducer 04: XCFun 2.1.1 SCANX/SCANC edge bugs (alpha=1 and s=0)
and TFK (Thomas-Fermi kinetic) spin-polarization bug.
"""
import numpy as np
from pyscf.dft import xcfun, libxc

def make_pt(rho, s, alpha):
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
    tau = sig / (8 * rho) + alpha * tau_unif
    return np.array([[rho], [np.sqrt(sig)], [0.0], [0.0], [0.0], [tau]])

print("=== 1. XCFun SCANX & SCANC at alpha = 1.0 and s = 0.0 (rho = 0.1) ===")
for s, alpha in [(0.0, 1.0), (1e-4, 1.0), (1e-2, 1.0), (0.0, 1.2), (1e-4, 1.2), (1e-8, 1.0 + 1e-10)]:
    pt = make_pt(0.1, s, alpha)
    ex, vx, _, _ = xcfun.eval_xc("SCANX", pt, spin=0, deriv=1)
    ec, vc, _, _ = xcfun.eval_xc("SCANC", pt, spin=0, deriv=1)
    exl, vxl, _, _ = libxc.eval_xc("MGGA_X_SCAN", pt, spin=0, deriv=1)
    ecl, vcl, _, _ = libxc.eval_xc("MGGA_C_SCAN", pt, spin=0, deriv=1)
    print(f"  s={s:7.1e} a={alpha:<12.10g} | SCANX: ex={ex[0]:+.4e} vs={vx[1][0]:+.4e} (L:{vxl[1][0]:+.4e}) | SCANC: ec={ec[0]:+.4e} vs={vc[1][0]:+.4e} (L:{vcl[1][0]:+.4e})")

print("\n=== 2. XCFun TFK (Thomas-Fermi kinetic) ignores spin polarization ===")
CF = 0.3 * (3 * np.pi**2)**(2/3)
for zeta in [0.0, 0.4, 0.8, 1.0]:
    rho = 0.2
    ra = 0.5 * rho * (1 + zeta)
    rb = 0.5 * rho * (1 - zeta)
    pa = np.array([[ra]])
    pb = np.array([[rb]])
    t_xcfun = xcfun.eval_xc("TFK", (pa, pb), spin=1, deriv=0)[0][0]
    t_libxc = libxc.eval_xc("LDA_K_TF", (pa, pb), spin=1, deriv=0)[0][0]
    t_exact = CF * rho**(2/3) * 0.5 * ((1 + zeta)**(5/3) + (1 - zeta)**(5/3))
    print(f"  zeta={zeta:3.1f}: XCFun_TFK={t_xcfun:.8f}, Libxc_LDA_K_TF={t_libxc:.8f}, exact={t_exact:.8f}, rel_err={(t_xcfun-t_exact)/t_exact:+.2%}")
