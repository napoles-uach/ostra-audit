#!/usr/bin/env python3
"""Minimal reproducer 02: XCFun 2.1.1 REVTPSSC (id=43) silent z^2 vs z^3 exponent bug + (s=0, alpha=0) NaN.
In Perdew et al., PRL 103, 026403 (2009), Eq. (24) is:
    eps_c = eps_c^PKZB * (1 + d * eps_c^PKZB * z^3),  d = 2.8,  z = tau_W / tau.
XCFun 2.1.1 evaluates z^2 instead of z^3, passing UEG (z=0), 1-electron/He (z=1), and internal FD4,
while silently shifting molecular correlation energies by +1.5 to +2.1 kcal/mol!
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

print("=== 1. Proof of z^2 vs z^3 exponent bug across (rho, s, alpha) ===")
d = 2.8
max_z3 = max_z2 = 0.0
for rho in [0.001, 0.1, 1.0, 5.0]:
    for s in [0.1, 0.5, 1.5]:
        for alpha in [0.01, 0.5, 1.0, 2.0]:
            pt = make_pt(rho, s, alpha)
            z = (5/3 * s**2) / (5/3 * s**2 + alpha)
            ec_x = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
            ec_l = libxc.eval_xc("MGGA_C_REVTPSS", pt, spin=0, deriv=0)[0][0]
            a_q = d * z**3
            y = (-1.0 + np.sqrt(1.0 + 4.0 * a_q * ec_l)) / (2.0 * a_q) if a_q > 1e-15 else ec_l
            ec_z2 = y * (1.0 + d * y * z**2)
            max_z3 = max(max_z3, abs(ec_x - ec_l))
            max_z2 = max(max_z2, abs(ec_x - ec_z2))
print(f"  max |XCFun - Libxc (exact z^3)| = {max_z3:.4e}")
print(f"  max |XCFun - Libxc_with_z^2|    = {max_z2:.4e}  (remaining ~1e-7 is beta_a=0.066725 vs 0.06672455)")

pt_zero = make_pt(0.1, 0.0, 0.0)
print(f"  Edge (s=0, alpha=0): XCFun={xcfun.eval_xc('REVTPSSC', pt_zero, spin=0, deriv=0)[0][0]}, Libxc={libxc.eval_xc('MGGA_C_REVTPSS', pt_zero, spin=0, deriv=0)[0][0]:.8f}")

print("\n=== 2. Observable SCF impact (TPSSx + revTPSSc, def2-TZVP, grid level 4) ===")
for name, atom in [
    ("H2O", "O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587"),
    ("N2",  "N 0 0 0; N 0 0 1.0977"),
    ("CH4", "C 0 0 0; H 0.629 0.629 0.629; H -0.629 -0.629 0.629; H -0.629 0.629 -0.629; H 0.629 -0.629 -0.629"),
]:
    mol = gto.M(atom=atom, basis="def2-tzvp", verbose=0)
    mf_l = dft.RKS(mol); mf_l.xc = "MGGA_X_TPSS,MGGA_C_REVTPSS"; mf_l.grids.level = 4; mf_l.conv_tol = 1e-10; el = mf_l.kernel()
    mf_x = dft.RKS(mol); mf_x._numint.libxc = dft.xcfun; mf_x.xc = "TPSSX,REVTPSSC"; mf_x.grids.level = 4; mf_x.conv_tol = 1e-10; ex = mf_x.kernel()
    print(f"  {name:4s}: E(Libxc)={el:.8f}, E(XCFun)={ex:.8f}, diff={ex-el:+.6f} Eh ({(ex-el)*627.5095:+.4f} kcal/mol)")
