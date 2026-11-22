import numpy as np
from pyscf.dft import xcfun, libxc

def make_pt(rho, s, alpha):
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
    tau_w = sig / (8 * rho)
    tau = tau_w + alpha * tau_unif
    return np.array([[rho], [np.sqrt(sig)], [0.0], [0.0], [0.0], [tau]])

# =====================================================================
# 1. VERIFY EXACT ROOT CAUSE OF XCFun REVTPSSC: z^2 vs z^3 in Eq. (24)!
# =====================================================================
print("=== 1. Verify XCFun REVTPSSC == revTPSSc with z^2 instead of z^3 in Eq. (24) ===")
# From Libxc: ec_revtpss = eps_PKZB * (1 + d * eps_PKZB * z^3), d = 2.8
# Solve quadratic for eps_PKZB from Libxc ec_revtpss:
#   d * z^3 * eps_PKZB^2 + eps_PKZB - ec_libxc = 0
# Then predict ec_buggy = eps_PKZB * (1 + d * eps_PKZB * z^2) !
d_param = 2.8
max_diff_z3 = 0.0
max_diff_z2 = 0.0
for rho in [0.001, 0.05, 0.1, 0.5, 1.0, 5.0]:
    for s in [0.1, 0.5, 1.0, 2.0]:
        for alpha in [0.01, 0.5, 1.0, 2.0, 10.0]:
            pt = make_pt(rho, s, alpha)
            z = (5/3 * s**2) / (5/3 * s**2 + alpha)
            ec_x = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
            ec_l = libxc.eval_xc("MGGA_C_REVTPSS", pt, spin=0, deriv=0)[0][0]
            # Solve d*z^3 * y^2 + y - ec_l = 0 (take root near ec_l)
            a_q = d_param * z**3
            if a_q > 1e-15:
                y = (-1.0 + np.sqrt(1.0 + 4.0 * a_q * ec_l)) / (2.0 * a_q)
            else:
                y = ec_l
            ec_pred_z2 = y * (1.0 + d_param * y * z**2)
            max_diff_z3 = max(max_diff_z3, abs(ec_x - ec_l))
            max_diff_z2 = max(max_diff_z2, abs(ec_x - ec_pred_z2))
print(f"  Across 120 (rho, s, alpha) points:")
print(f"    max |XCFun - Libxc (z^3)|          = {max_diff_z3:.4e}")
print(f"    max |XCFun - Libxc_with_z^2_bug|   = {max_diff_z2:.4e}")
# Check if the remaining ~1e-7 is purely beta_a = 0.066725 vs 0.0667245506!
pt = make_pt(1.0, 0.5, 1.0)
z = (5/3 * 0.25) / (5/3 * 0.25 + 1.0)
ec_x = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
ec_l = libxc.eval_xc("MGGA_C_REVTPSS", pt, spin=0, deriv=0)[0][0]
y = (-1.0 + np.sqrt(1.0 + 4.0 * d_param * z**3 * ec_l)) / (2.0 * d_param * z**3)
print(f"  Single point (rho=1, s=0.5, a=1): XCFun={ec_x:.10f}, Libxc(z^3)={ec_l:.10f}, Libxc(z^2)={y*(1+d_param*y*z**2):.10f}, residual={ec_x - y*(1+d_param*y*z**2):+.3e}")

# Also test spin-polarized!
print("  Spin-polarized check (rho_a=0.12, rho_b=0.05, grad_a=0.1, grad_b=0.03, tau_a=0.08, tau_b=0.04):")
ra, rb = 0.12, 0.05
ga, gb = 0.10, 0.03
ta, tb = 0.08, 0.04
pt_a = np.array([[ra], [ga], [0.0], [0.0], [0.0], [ta]])
pt_b = np.array([[rb], [gb], [0.0], [0.0], [0.0], [tb]])
ec_xp = xcfun.eval_xc("REVTPSSC", (pt_a, pt_b), spin=1, deriv=0)[0][0]
ec_lp = libxc.eval_xc("MGGA_C_REVTPSS", (pt_a, pt_b), spin=1, deriv=0)[0][0]
z_pol = ((ga + gb)**2 / (8 * (ra + rb))) / (ta + tb)
yp = (-1.0 + np.sqrt(1.0 + 4.0 * d_param * z_pol**3 * ec_lp)) / (2.0 * d_param * z_pol**3)
print(f"    XCFun={ec_xp:.10f}, Libxc(z^3)={ec_lp:.10f}, Libxc(z^2)={yp*(1+d_param*yp*z_pol**2):.10f}, residual={ec_xp - yp*(1+d_param*yp*z_pol**2):+.3e}")

# =====================================================================
# 2. INSPECT POLE & FORMULA BUG IN XCFun REVTPSSX!
# =====================================================================
print("\n=== 2. Scan XCFun REVTPSSX pole in rho (s=0.5, alpha=1.0) ===")
ex_lda = lambda r: -0.75 * (3 / np.pi)**(1/3) * r**(1/3)
for rho in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]:
    pt = make_pt(rho, 0.5, 1.0)
    ex_x = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0] - xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    Fx = ex_x / ex_lda(rho)
    # Invert Fx = 1 + kappa - kappa / (1 + x/kappa) => x = kappa * (Fx - 1) / (1 + kappa - Fx)
    kap = 0.804
    x_eff = kap * (Fx - 1.0) / (1.0 + kap - Fx)
    print(f"  rho={rho:4.1f}: Fx={Fx:+12.6f}, x_eff={x_eff:+14.6f}")

print("\n=== 2b. Check XCFun REVTPSSX at small rho (rho=1e-4) vs Libxc TPSS / revTPSS ===")
# Compare x_eff(xcfun) vs x_libxc_revtpss across (s, alpha) at rho=1e-4!
for s, alpha in [(0.1, 1.0), (0.5, 1.0), (1.0, 1.0), (0.5, 0.0), (0.5, 2.0), (0.0, 0.5), (0.0, 2.0)]:
    rho = 1e-4
    pt = make_pt(rho, s, alpha)
    ex_x = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0] - xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ex_l = libxc.eval_xc("MGGA_X_REVTPSS", pt, spin=0, deriv=0)[0][0]
    ex_t = libxc.eval_xc("MGGA_X_TPSS", pt, spin=0, deriv=0)[0][0]
    Fx_x = ex_x / ex_lda(rho); Fx_l = ex_l / ex_lda(rho); Fx_t = ex_t / ex_lda(rho)
    kap = 0.804
    xx = kap * (Fx_x - 1) / (1 + kap - Fx_x)
    xl = kap * (Fx_l - 1) / (1 + kap - Fx_l)
    print(f"  s={s:3.1f} a={alpha:3.1f} | Fx_xcfun={Fx_x:.8f} Fx_revtpss={Fx_l:.8f} Fx_tpss={Fx_t:.8f} | x_xcfun={xx:.8f} x_libxc={xl:.8f} diff_x={xx-xl:+.6e}")
