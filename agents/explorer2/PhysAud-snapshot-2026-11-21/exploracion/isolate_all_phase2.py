import subprocess, sys
import numpy as np
from pyscf import gto, dft
from pyscf.dft import xcfun, libxc
from pyscf.pbc import gto as pgto, scf as pscf, dft as pdft

def make_pt(rho, s, alpha):
    kf = (3 * np.pi**2 * rho)**(1/3)
    sig = (2 * kf * rho * s)**2
    tau_unif = 0.3 * (3 * np.pi**2)**(2/3) * rho**(5/3)
    tau_w = sig / (8 * rho)
    tau = tau_w + alpha * tau_unif
    return np.array([[rho], [np.sqrt(sig)], [0.0], [0.0], [0.0], [tau]])

# =====================================================================
# 1. ISOLATE XCFun REVTPSSX scaling bug!
# =====================================================================
print("=== 1. XCFun REVTPSSX scaling of (Fx - 1) vs rho and s ===")
# At alpha=1, s=0.5, vary rho:
ex_lda_fn = lambda r: -0.75 * (3 / np.pi)**(1/3) * r**(1/3)
for rho in [0.01, 0.1, 0.5, 1.0, 2.0, 10.0]:
    pt = make_pt(rho, 0.5, 1.0)
    ex_x = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0] - xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ex_l = libxc.eval_xc("MGGA_X_REVTPSS", pt, spin=0, deriv=0)[0][0]
    Fx_x = ex_x / ex_lda_fn(rho)
    Fx_l = ex_l / ex_lda_fn(rho)
    # Also test what s_eff in Libxc gives the same Fx_x!
    # Binary search s_eff in [0, 20]
    lo, hi = 0.0, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        em = libxc.eval_xc("MGGA_X_REVTPSS", make_pt(rho, mid, 1.0), spin=0, deriv=0)[0][0] / ex_lda_fn(rho)
        if em < Fx_x: lo = mid
        else: hi = mid
    s_eff = 0.5 * (lo + hi)
    print(f"  rho={rho:6.2f} s=0.50 a=1.0 | Fx_x={Fx_x:.8f} Fx_l={Fx_l:.8f} | s_eff={s_eff:.8f} | s_eff/s={s_eff/0.5:.8f} | (s_eff/s)*rho^(1/3)={(s_eff/0.5)*rho**(1/3):.8f}")

# Check at s=0.01 (small s expansion: Fx = 1 + mu_eff * s^2)
print("  Small-s (s=0.01, alpha=1.0) effective coefficient (Fx-1)/s^2:")
for rho in [0.01, 0.1, 1.0, 10.0]:
    s = 0.01
    pt = make_pt(rho, s, 1.0)
    ex_x = xcfun.eval_xc("REVTPSSX,REVTPSSC", pt, spin=0, deriv=0)[0][0] - xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
    ex_l = libxc.eval_xc("MGGA_X_REVTPSS", pt, spin=0, deriv=0)[0][0]
    c_x = (ex_x / ex_lda_fn(rho) - 1.0) / s**2
    c_l = (ex_l / ex_lda_fn(rho) - 1.0) / s**2
    print(f"    rho={rho:6.2f}: c_xcfun={c_x:.8f}, c_libxc={c_l:.8f}, ratio={c_x/c_l:.8f}, ratio*rho^(2/3)={(c_x/c_l)*rho**(2/3):.8f}, ratio*rho^(1/3)={(c_x/c_l)*rho**(1/3):.8f}")

# =====================================================================
# 2. ISOLATE XCFun REVTPSSC difference vs Libxc (check GGA_C_REGTPSS vs PBE!)
# =====================================================================
print("\n=== 2. XCFun REVTPSSC vs Libxc MGGA_C_REVTPSS & TPSS ===")
# Check spin-polarized zeta dependence too!
for rho in [0.1, 1.0]:
    for s in [0.5]:
        for alpha in [0.0, 0.5, 1.0, 2.0, 10.0]:
            pt = make_pt(rho, s, alpha)
            ec_x = xcfun.eval_xc("REVTPSSC", pt, spin=0, deriv=0)[0][0]
            ec_l = libxc.eval_xc("MGGA_C_REVTPSS", pt, spin=0, deriv=0)[0][0]
            ec_regpbe = libxc.eval_xc("GGA_C_REGTPSS", pt[:4], spin=0, deriv=0)[0][0]
            # Wait: z = tau_w / tau = (5/3 * s^2) / (5/3 * s^2 + alpha)
            z = (5/3 * s**2) / (5/3 * s**2 + alpha)
            print(f"  rho={rho:4.1f} s={s} a={alpha:4.1f} (z={z:.4f}) | xcfun={ec_x:.9f} libxc={ec_l:.9f} diff={ec_x-ec_l:+.4e} | regtpss_gga={ec_regpbe:.9f}")

# =====================================================================
# 3. PySCF PBC Skewed Cell: RHF (no XC grid!) vs RKS (UniformGrids vs BeckeGrids)
# =====================================================================
print("\n=== 3. PySCF PBC Skewed Cell: RHF vs RKS (GDF) ===")
a0 = np.array([[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]])
Us = [
    ("orig",  np.eye(3)),
    ("skew1", np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])),
    ("skew2", np.array([[1, 2, 1], [0, 1, 1], [0, 0, 1]])),
]
e_rhf0 = e_becke0 = None
for lbl, U in Us:
    cell = pgto.M(atom="Si 0 0 0; Si 1.3575 1.3575 1.3575", a=U @ a0, basis="gth-szv", pseudo="gth-pbe", verbose=0)
    e_rhf = pscf.RHF(cell).density_fit().set(conv_tol=1e-10).kernel()
    mf_b = pdft.RKS(cell).density_fit().set(xc="pbe", conv_tol=1e-10)
    mf_b.grids = pdft.gen_grid.BeckeGrids(cell)
    mf_b.grids.level = 3
    e_becke = mf_b.kernel()
    if lbl == "orig":
        e_rhf0, e_becke0 = e_rhf, e_becke
    print(f"  {lbl:6s}: dE_RHF(GDF)={e_rhf-e_rhf0:+.3e} | dE_RKS(GDF,Becke3)={e_becke-e_becke0:+.3e}", flush=True)

# =====================================================================
# 4. libcint int2e_erfc impact on real cc-pV5Z / cc-pV6Z basis with HSE06 / LC-wPBE!
# =====================================================================
print("\n=== 4. libcint int2e_erfc (SR-ERI) on real Ne cc-pV5Z / cc-pV6Z ===")
for basis in ["cc-pv5z", "cc-pv6z"]:
    mol = gto.M(atom="Ne 0 0 0", basis=basis, verbose=0)
    lmax = max(mol.bas_angular(i) for i in range(mol.nbas))
    # Pick only the highest-l shell to test fast!
    hi_sh = [i for i in range(mol.nbas) if mol.bas_angular(i) == lmax]
    print(f"  Ne {basis}: lmax={lmax}, nbas={mol.nbas}, max-l shells={hi_sh}")
    # Test partition of unity on (hi, hi | hi, hi) for realistic omegas!
    shls = (hi_sh[0], hi_sh[0], hi_sh[0], hi_sh[0])
    eri_full = mol.intor_by_shell("int2e", shls)
    for w in [0.11, 0.2, 0.33, 0.4, 1.0, 2.0, 5.0, 10.0, 20.0]:
        with mol.with_range_coulomb(w):
            eri_lr = mol.intor_by_shell("int2e", shls)
        with mol.with_range_coulomb(-w):
            eri_sr = mol.intor_by_shell("int2e", shls)
        err = np.max(np.abs((eri_lr + eri_sr) - eri_full))
        if err > 1e-10:
            print(f"    w={w:5.2f}: max|LR+SR-Full|={err:.3e}, max|Full|={np.max(np.abs(eri_full)):.3e}, max|SR|={np.max(np.abs(eri_sr)):.3e}")
        else:
            print(f"    w={w:5.2f}: CLEAN ({err:.2e})")
