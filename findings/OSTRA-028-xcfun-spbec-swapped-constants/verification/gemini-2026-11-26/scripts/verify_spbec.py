#!/usr/bin/env python3
"""
Independent verification (Phase 1, blind) of XCFun 'SPBEC' (src/functionals/spbec.cpp)
against first-principles implementation of Swart, Sola & Bickelhaupt (JCP 131, 094103, 2009)
and libxc GGA_C_SPBE.
"""
import numpy as np
from pyscf import gto, dft

# ============================================================================
# Independent from-scratch implementations of LDA correlation (PW92 and VWN5)
# and sPBE gradient correction H_sPBE
# ============================================================================

def pw92_G(rs, A, a1, b1, b2, b3, b4):
    sq = np.sqrt(rs)
    Q1 = 2.0 * A * (b1 * sq + b2 * rs + b3 * rs * sq + b4 * rs * rs)
    return -2.0 * A * (1.0 + a1 * rs) * np.log(1.0 + 1.0 / Q1)

def pw92_eps(na, nb):
    """Perdew-Wang 1992 LDA correlation energy per electron (exact modified constants)."""
    n = na + nb
    zeta = (na - nb) / n
    rs = (3.0 / (4.0 * np.pi * n)) ** (1.0 / 3.0)
    # PW92 Table I / Perdew-Wang 1992 (same high-precision values used in PBE/libxc/XCFun)
    ec0 = pw92_G(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
    ec1 = pw92_G(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
    ac  = pw92_G(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)
    f_zeta = ((1.0 + zeta)**(4.0/3.0) + (1.0 - zeta)**(4.0/3.0) - 2.0) / (2.0**(4.0/3.0) - 2.0)
    fpp0 = 4.0 / (9.0 * (2.0**(1.0/3.0) - 1.0))
    z4 = zeta**4
    return ec0 - ac * f_zeta / fpp0 * (1.0 - z4) + (ec1 - ec0) * f_zeta * z4

def vwn_Q(x, A, x0, b, c):
    X = x * x + b * x + c
    X0 = x0 * x0 + b * x0 + c
    Q = np.sqrt(4.0 * c - b * b)
    return A * (
        np.log(x * x / X)
        + 2.0 * b / Q * np.arctan(Q / (2.0 * x + b))
        - b * x0 / X0 * (
            np.log((x - x0)**2 / X)
            + 2.0 * (b + 2.0 * x0) / Q * np.arctan(Q / (2.0 * x + b))
        )
    )

def vwn5_eps(na, nb):
    """Vosko-Wilk-Nusair VWN5 LDA correlation energy per electron."""
    n = na + nb
    zeta = (na - nb) / n
    rs = (3.0 / (4.0 * np.pi * n)) ** (1.0 / 3.0)
    x = np.sqrt(rs)
    # VWN5 parameters (Table 5 of Vosko, Wilk, Nusair 1980, Can. J. Phys. 58, 1200)
    ec0 = vwn_Q(x, 0.0310907, -0.10498, 3.72744, 12.9352)
    ec1 = vwn_Q(x, 0.01554535, -0.32500, 7.06042, 18.0578)
    ac  = vwn_Q(x, -1.0 / (6.0 * np.pi * np.pi), -0.0047584, 1.13107, 13.0045)
    f_zeta = ((1.0 + zeta)**(4.0/3.0) + (1.0 - zeta)**(4.0/3.0) - 2.0) / (2.0**(4.0/3.0) - 2.0)
    fpp0 = 4.0 / (9.0 * (2.0**(1.0/3.0) - 1.0))
    z4 = zeta**4
    return ec0 + ac * f_zeta / fpp0 * (1.0 - z4) + (ec1 - ec0) * f_zeta * z4

def spbe_eps(na, nb, gnn, lda="pw92", gamma=None, beta=None, buggy_expm1=False):
    """
    sPBE correlation energy per electron:
      eps_c^sPBE = eps_c^LDA + gamma * phi^3 * ln(1 + (beta/gamma) * t^2 / (1 + A * t^2))
      A = (beta/gamma) / expm1(-eps_c^LDA / (gamma * phi^3))
    """
    if gamma is None:
        gamma = (1.0 - np.log(2.0)) / (np.pi**2)  # 0.031090690869654895
    if beta is None:
        beta = 0.06672455060314922
    bg = beta / gamma

    n = na + nb
    zeta = (na - nb) / n
    phi = 0.5 * ((1.0 + zeta)**(2.0/3.0) + (1.0 - zeta)**(2.0/3.0))
    phi3 = phi**3
    t2 = ((np.pi / 3.0)**(1.0/3.0) / 16.0) * gnn / (phi * phi * n**(7.0/3.0))

    eps_lda = pw92_eps(na, nb) if lda == "pw92" else vwn5_eps(na, nb)
    denom = (bg * phi3) if buggy_expm1 else (gamma * phi3)
    A = bg / np.expm1(-eps_lda / denom)
    H = gamma * phi3 * np.log1p(bg * t2 / (1.0 + A * t2))
    return eps_lda + H

# Exact PBE constants vs rounded constants in paper
GAMMA_EXACT = (1.0 - np.log(2.0)) / (np.pi**2)
BETA_EXACT = 0.06672455060314922
GAMMA_PAPER = 0.031091
BETA_PAPER = 0.066725

print("=" * 80)
print("1. MOLECULAR DENSITY: H2O (RKS PBE/cc-pVTZ)")
print("=" * 80)
mol = gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis="cc-pvtz", verbose=0)
mf = dft.RKS(mol)
mf.xc = "PBE"
mf.kernel()
dm = mf.make_rdm1()
grids = mf.grids

ni_xcfun = dft.numint.NumInt(); ni_xcfun.libxc = dft.xcfun
ni_libxc = dft.numint.NumInt(); ni_libxc.libxc = dft.libxc

_, E_xcfun_spbec, _ = ni_xcfun.nr_rks(mol, grids, ",SPBEC", dm)
_, E_libxc_spbec, _ = ni_libxc.nr_rks(mol, grids, ",GGA_C_SPBE", dm)
_, E_xcfun_pbec, _  = ni_xcfun.nr_rks(mol, grids, ",PBEC", dm)
_, E_xcfun_vpbec, _ = ni_xcfun.nr_rks(mol, grids, ",VWN_PBEC", dm)
_, E_libxc_pbec, _  = ni_libxc.nr_rks(mol, grids, ",GGA_C_PBE", dm)

# Evaluate our independent sPBE directly on the same H2O grid
ao = ni_libxc.eval_ao(mol, grids.coords, deriv=1)
rho_g = ni_libxc.eval_rho(mol, ao, dm, xctype="GGA")
n_g = rho_g[0]
gnn_g = rho_g[1]**2 + rho_g[2]**2 + rho_g[3]**2
mask = n_g > 1e-12
na_g = 0.5 * n_g[mask]
nb_g = 0.5 * n_g[mask]
gnn_m = gnn_g[mask]
w_m = grids.weights[mask]

E_ours_buggy = np.sum(w_m * n_g[mask] * spbe_eps(na_g, nb_g, gnn_m, lda="vwn5", gamma=BETA_PAPER, beta=GAMMA_PAPER, buggy_expm1=True))
E_ours_vwn5_paper = np.sum(w_m * n_g[mask] * spbe_eps(na_g, nb_g, gnn_m, lda="vwn5", gamma=GAMMA_PAPER, beta=BETA_PAPER))
E_ours_vwn5_exact = np.sum(w_m * n_g[mask] * spbe_eps(na_g, nb_g, gnn_m, lda="vwn5", gamma=GAMMA_EXACT, beta=BETA_EXACT))
E_ours_pw92_paper = np.sum(w_m * n_g[mask] * spbe_eps(na_g, nb_g, gnn_m, lda="pw92", gamma=GAMMA_PAPER, beta=BETA_PAPER))
E_ours_pw92_exact = np.sum(w_m * n_g[mask] * spbe_eps(na_g, nb_g, gnn_m, lda="pw92", gamma=GAMMA_EXACT, beta=BETA_EXACT))

print(f"XCFun 'SPBEC' (buggy as-is)               : {E_xcfun_spbec:+.12f} Eh")
print(f"Our reproduction of buggy spbec.cpp       : {E_ours_buggy:+.12f} Eh  (diff vs XCFun = {E_ours_buggy - E_xcfun_spbec:+.3e})")
print(f"libxc 'GGA_C_SPBE'                        : {E_libxc_spbec:+.12f} Eh")
print(f"Our sPBE (PW92, exact PBE consts)         : {E_ours_pw92_exact:+.12f} Eh  (diff vs libxc = {E_ours_pw92_exact - E_libxc_spbec:+.3e})")
print(f"Our sPBE (PW92, rounded 0.031091/0.066725): {E_ours_pw92_paper:+.12f} Eh  (diff vs libxc = {E_ours_pw92_paper - E_libxc_spbec:+.3e})")
print(f"Our sPBE (VWN5, exact PBE consts)         : {E_ours_vwn5_exact:+.12f} Eh  (diff vs libxc = {E_ours_vwn5_exact - E_libxc_spbec:+.3e})")
print(f"Our sPBE (VWN5, rounded 0.031091/0.066725): {E_ours_vwn5_paper:+.12f} Eh  (diff vs libxc = {E_ours_vwn5_paper - E_libxc_spbec:+.3e})")
print(f"--- Error magnitude on H2O ---")
print(f"XCFun SPBEC - libxc GGA_C_SPBE            : {E_xcfun_spbec - E_libxc_spbec:+.6f} Eh ({(E_xcfun_spbec - E_libxc_spbec)*627.509474:+.2f} kcal/mol, rel={(E_xcfun_spbec - E_libxc_spbec)/abs(E_libxc_spbec)*100:+.2f}%)")
print(f"XCFun SPBEC - sPBE(VWN5,fixed)            : {E_xcfun_spbec - E_ours_vwn5_paper:+.6f} Eh ({(E_xcfun_spbec - E_ours_vwn5_paper)*627.509474:+.2f} kcal/mol)")
print(f"Ref PBE diff (VWN_PBEC - PBEC)            : {E_xcfun_vpbec - E_xcfun_pbec:+.6f} Eh ({(E_xcfun_vpbec - E_xcfun_pbec)*627.509474:+.2f} kcal/mol)")

print("\n" + "=" * 80)
print("2. SPIN-POLARIZED SINGLE POINT TEST + 4-PT FD DERIVATIVES")
print("=" * 80)
na_p, nb_p = 0.35, 0.15
ga_p = np.array([0.4, 0.1, 0.0])
gb_p = np.array([-0.1, 0.3, 0.2])
gnn_p = float(np.sum((ga_p + gb_p)**2))
n_p = na_p + nb_p
rho_pol = (
    np.array([[na_p], [ga_p[0]], [ga_p[1]], [ga_p[2]]]),
    np.array([[nb_p], [gb_p[0]], [gb_p[1]], [gb_p[2]]]),
)
exc_xcf, vxc_xcf, _, _ = dft.xcfun.eval_xc(",SPBEC", rho_pol, spin=1, deriv=1)
exc_lib, vxc_lib, _, _ = dft.libxc.eval_xc(",GGA_C_SPBE", rho_pol, spin=1, deriv=1)

eps_ours_buggy = float(spbe_eps(na_p, nb_p, gnn_p, lda="vwn5", gamma=BETA_PAPER, beta=GAMMA_PAPER, buggy_expm1=True))
eps_ours_pw92  = float(spbe_eps(na_p, nb_p, gnn_p, lda="pw92", gamma=GAMMA_EXACT, beta=BETA_EXACT))
eps_ours_vwn5  = float(spbe_eps(na_p, nb_p, gnn_p, lda="vwn5", gamma=GAMMA_PAPER, beta=BETA_PAPER))

print(f"Point: na={na_p}, nb={nb_p}, gnn={gnn_p:.4f}")
print(f"  XCFun 'SPBEC' eps_c                     : {exc_xcf[0]:+.15f}")
print(f"  Our buggy model eps_c                   : {eps_ours_buggy:+.15f}  (diff = {eps_ours_buggy - exc_xcf[0]:+.3e})")
print(f"  libxc 'GGA_C_SPBE' eps_c                : {exc_lib[0]:+.15f}")
print(f"  Our sPBE(PW92, exact) eps_c             : {eps_ours_pw92:+.15f}  (diff vs libxc = {eps_ours_pw92 - exc_lib[0]:+.3e})")
print(f"  Our sPBE(VWN5, paper) eps_c             : {eps_ours_vwn5:+.15f}  (diff vs libxc = {eps_ours_vwn5 - exc_lib[0]:+.3e})")
print(f"  Error XCFun SPBEC vs libxc GGA_C_SPBE   : {exc_xcf[0] - exc_lib[0]:+.6e} (rel = {(exc_xcf[0]-exc_lib[0])/abs(exc_lib[0])*100:+.2f}%)")
print(f"  vrho XCFun SPBEC                        : {vxc_xcf[0][0]}")
print(f"  vrho libxc GGA_C_SPBE                   : {vxc_lib[0][0]}")
print(f"  vsigma XCFun SPBEC                      : {vxc_xcf[1][0]}")
print(f"  vsigma libxc GGA_C_SPBE                 : {vxc_lib[1][0]}")

print("\n" + "=" * 80)
print("3. PHYSICAL EXACT LIMIT CHECK: rapidly varying limit (t -> inf, eps_c -> 0)")
print("   Exact PBE/sPBE constraint: lim_{t->inf} (eps_LDA + H) = 0")
print("=" * 80)
for grad_mag in [1.0, 10.0, 100.0, 1e4]:
    rho_lt = np.array([[0.1], [grad_mag], [0.0], [0.0]])
    e_x, _, _, _ = dft.xcfun.eval_xc(",SPBEC", rho_lt, spin=0, deriv=0)
    e_l, _, _, _ = dft.libxc.eval_xc(",GGA_C_SPBE", rho_lt, spin=0, deriv=0)
    print(f"  n=0.1, |grad n|={grad_mag:7.1f}: XCFun SPBEC eps_c = {e_x[0]:+.10e}, libxc GGA_C_SPBE eps_c = {e_l[0]:+.10e}")
