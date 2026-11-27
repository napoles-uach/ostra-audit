#!/usr/bin/env python3
"""
Independent verification (Phase 1, blind) of XCFun 'TW' (src/functionals/tw.cpp)
against exact analytical / physical oracles and libxc GGA_K_VW.
"""
import numpy as np
from pyscf import gto, dft

print("=" * 78)
print("1. EXACT ANALYTICAL ORACLE: 1-electron H atom (exact 1s orbital)")
print("   rho_a(r) = (1/pi)*exp(-2r), rho_b = 0, |grad rho_a| = 2*rho_a")
print("   Exact T_W = int |grad rho|^2 / (8*rho) d^3r = 0.5 Eh")
print("   Analytic integral of buggy XCFun TW = (1/8)*(4*rho^2)^2/rho = 2/(27*pi^2)")
print(f"   2 / (27 * pi^2) = {2.0 / (27.0 * np.pi**2):.15f} Eh")
print("=" * 78)

# Build a PySCF numerical grid centered on H atom (level 5 for high accuracy)
mol_h = gto.M(atom="H 0 0 0", basis="cc-pvqz", spin=1, verbose=0)
grids = dft.gen_grid.Grids(mol_h)
grids.level = 5
grids.build()
coords = grids.coords
weights = grids.weights
r = np.linalg.norm(coords, axis=1)

# Exact 1s density and gradient (spin-polarized: alpha=exact 1s, beta=0)
rho_a = (1.0 / np.pi) * np.exp(-2.0 * r)
# grad rho_a = -2 * rho_a * r_vec / r
r_safe = np.where(r > 0, r, 1.0)
grad_a = (-2.0 * rho_a / r_safe)[:, None] * coords
grad_a[r == 0] = 0.0
rho_b = np.zeros_like(rho_a)
grad_b = np.zeros_like(grad_a)

rho_spin_exact = (
    np.vstack([rho_a, grad_a[:, 0], grad_a[:, 1], grad_a[:, 2]]),
    np.vstack([rho_b, grad_b[:, 0], grad_b[:, 1], grad_b[:, 2]]),
)

# Direct formula evaluations
gaa = np.sum(grad_a**2, axis=1)
tw_exact_grid = np.sum(weights * gaa / (8.0 * rho_a))
tw_buggy_grid = np.sum(weights * (1.0 / 8.0) * (gaa**2) / rho_a)

exc_xcfun_tw, _, _, _ = dft.xcfun.eval_xc("TW", rho_spin_exact, spin=1, deriv=1)
exc_xcfun_vwk, _, _, _ = dft.xcfun.eval_xc("VWK", rho_spin_exact, spin=1, deriv=1)
exc_libxc_vw, _, _, _ = dft.libxc.eval_xc("GGA_K_VW,", rho_spin_exact, spin=1, deriv=1)

# Note: eval_xc returns energy per particle eps (so E = sum(w * n * eps))
E_xcfun_tw = np.sum(weights * rho_a * exc_xcfun_tw)
E_xcfun_vwk = np.sum(weights * rho_a * exc_xcfun_vwk)
E_libxc_vw = np.sum(weights * rho_a * exc_libxc_vw)

print(f"Grid integral exact N = int rho d3r       : {np.sum(weights * rho_a):.15f}")
print(f"Grid integral exact T_W                   : {tw_exact_grid:.15f} Eh (err = {tw_exact_grid - 0.5:+.3e})")
print(f"Grid integral buggy formula (gaa^2)/(8n)  : {tw_buggy_grid:.15f} Eh (diff vs 2/(27pi^2) = {tw_buggy_grid - 2/(27*np.pi**2):+.3e})")
print(f"XCFun 'TW'  (exact H 1s, spin=1)          : {E_xcfun_tw:.15f} Eh (err vs 0.5 = {E_xcfun_tw - 0.5:+.6e}, ratio = {E_xcfun_tw / 0.5:.6e})")
print(f"XCFun 'VWK' (exact H 1s, spin=1)          : {E_xcfun_vwk:.15f} Eh (err vs 0.5 = {E_xcfun_vwk - 0.5:+.3e})")
print(f"libxc 'GGA_K_VW' (exact H 1s, spin=1)     : {E_libxc_vw:.15f} Eh (err vs 0.5 = {E_libxc_vw - 0.5:+.3e})")

print("\n" + "=" * 78)
print("2. PYSCF UHF H ATOM (cc-pVQZ) via mf._numint.libxc = dft.xcfun vs dft.libxc")
print("=" * 78)
mf_h = dft.UKS(mol_h)
mf_h.xc = "HF"
mf_h.kernel()
dm_h = mf_h.make_rdm1()
T_kin_h = np.einsum("ij,ji->", mol_h.intor("int1e_kin"), dm_h[0] + dm_h[1])

ni_xcfun = dft.numint.NumInt()
ni_xcfun.libxc = dft.xcfun
ni_libxc = dft.numint.NumInt()
ni_libxc.libxc = dft.libxc

_, E_pyscf_tw, _ = ni_xcfun.nr_uks(mol_h, grids, "TW", dm_h)
_, E_pyscf_vwk, _ = ni_xcfun.nr_uks(mol_h, grids, "VWK", dm_h)
_, E_pyscf_libxc, _ = ni_libxc.nr_uks(mol_h, grids, "GGA_K_VW,", dm_h)

print(f"Exact 1e kinetic energy <T> (UHF/cc-pVQZ) : {T_kin_h:.12f} Eh")
print(f"PySCF + XCFun 'TW'                        : {E_pyscf_tw:.12f} Eh (diff = {E_pyscf_tw - T_kin_h:+.6e}, ratio = {E_pyscf_tw / T_kin_h:.6e})")
print(f"PySCF + XCFun 'VWK'                       : {E_pyscf_vwk:.12f} Eh (diff = {E_pyscf_vwk - T_kin_h:+.3e})")
print(f"PySCF + libxc 'GGA_K_VW'                  : {E_pyscf_libxc:.12f} Eh (diff = {E_pyscf_libxc - T_kin_h:+.3e})")

print("\n" + "=" * 78)
print("3. MOLECULAR DENSITY: H2O (RKS PBE/cc-pVTZ, closed-shell unpolarized)")
print("=" * 78)
mol_w = gto.M(
    atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587",
    basis="cc-pvtz",
    verbose=0,
)
mf_w = dft.RKS(mol_w)
mf_w.xc = "PBE"
mf_w.kernel()
dm_w = mf_w.make_rdm1()
grids_w = mf_w.grids

_, E_w_tw, _ = ni_xcfun.nr_rks(mol_w, grids_w, "TW", dm_w)
_, E_w_vwk, _ = ni_xcfun.nr_rks(mol_w, grids_w, "VWK", dm_w)
_, E_w_libxc, _ = ni_libxc.nr_rks(mol_w, grids_w, "GGA_K_VW,", dm_w)

print(f"H2O PySCF + XCFun 'TW'                    : {E_w_tw:.12e} Eh")
print(f"H2O PySCF + XCFun 'VWK'                   : {E_w_vwk:.12f} Eh")
print(f"H2O PySCF + libxc 'GGA_K_VW'              : {E_w_libxc:.12f} Eh")
print(f"H2O ratio (XCFun TW / libxc GGA_K_VW)     : {E_w_tw / E_w_libxc:.6e} (orders of magnitude!)")
print(f"H2O diff  (XCFun VWK - libxc GGA_K_VW)    : {E_w_vwk - E_w_libxc:+.3e} Eh")

print("\n" + "=" * 78)
print("4. POINTWISE TEST: UNPOLARIZED AND PARTIALLY SPIN-POLARIZED")
print("   Distinguishing total-density |grad n|^2/(8n) vs spin-scaled (gaa/na + gbb/nb)/8")
print("=" * 78)

# Point 1: unpolarized (na=nb=0.25, ga=gb=[0.3, 0.4, 0.0])
# gaa=gbb=gab=0.25, gnn = 1.0, n = 0.5
rho_unpol = np.array([[0.5], [0.6], [0.8], [0.0]])
exc_u_tw, vxc_u_tw, _, _ = dft.xcfun.eval_xc("TW", rho_unpol, spin=0, deriv=1)
exc_u_vwk, vxc_u_vwk, _, _ = dft.xcfun.eval_xc("VWK", rho_unpol, spin=0, deriv=1)
exc_u_lib, vxc_u_lib, _, _ = dft.libxc.eval_xc("GGA_K_VW,", rho_unpol, spin=0, deriv=1)

n_u = 0.5
gnn_u = 1.0
e_u_exact = gnn_u / (8.0 * n_u)  # energy density t_W = 0.25
e_u_buggy = (1.0 / 8.0) * (0.25 + 0.25)**2 / n_u  # = 0.0625
print(f"Unpolarized point (n=0.5, |grad n|^2=1.0):")
print(f"  Exact t_W = |grad n|^2/(8n)             : {e_u_exact:.15f}")
print(f"  Buggy (gaa+gbb)^2/(8n)                  : {e_u_buggy:.15f}")
print(f"  XCFun 'TW'  n*eps                       : {n_u * exc_u_tw[0]:.15f}")
print(f"  XCFun 'VWK' n*eps                       : {n_u * exc_u_vwk[0]:.15f}")
print(f"  libxc 'GGA_K_VW' n*eps                  : {n_u * exc_u_lib[0]:.15f}")
print(f"  (Note: even removing pow(...,2), (gaa+gbb)/(8n) = {(0.25+0.25)/(8*n_u):.15f}, STILL off by 2x!)")

# Point 2: partially spin-polarized with non-collinear gradients
# na=0.35, nb=0.15, ga=[0.4, 0.1, 0.0], gb=[-0.1, 0.3, 0.2]
ga = np.array([0.4, 0.1, 0.0])
gb = np.array([-0.1, 0.3, 0.2])
na, nb = 0.35, 0.15
gaa_p = float(np.dot(ga, ga))
gab_p = float(np.dot(ga, gb))
gbb_p = float(np.dot(gb, gb))
gnn_p = gaa_p + 2.0 * gab_p + gbb_p
n_p = na + nb
rho_pol = (
    np.array([[na], [ga[0]], [ga[1]], [ga[2]]]),
    np.array([[nb], [gb[0]], [gb[1]], [gb[2]]]),
)
exc_p_tw, _, _, _ = dft.xcfun.eval_xc("TW", rho_pol, spin=1, deriv=1)
exc_p_vwk, _, _, _ = dft.xcfun.eval_xc("VWK", rho_pol, spin=1, deriv=1)
exc_p_lib, _, _, _ = dft.libxc.eval_xc("GGA_K_VW,", rho_pol, spin=1, deriv=1)

e_p_total = gnn_p / (8.0 * n_p)
e_p_spin_scaled = (gaa_p / na + gbb_p / nb) / 8.0
e_p_buggy = (1.0 / 8.0) * (gaa_p + gbb_p)**2 / n_p
print(f"\nPartially polarized point (na={na}, nb={nb}, gaa={gaa_p:.3f}, gab={gab_p:.3f}, gbb={gbb_p:.3f}, gnn={gnn_p:.3f}):")
print(f"  Total-density |grad n|^2 / (8n)         : {e_p_total:.15f}")
print(f"  Spin-scaled (gaa/na + gbb/nb) / 8       : {e_p_spin_scaled:.15f}")
print(f"  Buggy (gaa+gbb)^2 / (8n)                : {e_p_buggy:.15f}")
print(f"  XCFun 'TW'       n*eps                  : {n_p * exc_p_tw[0]:.15f}")
print(f"  XCFun 'VWK'      n*eps                  : {n_p * exc_p_vwk[0]:.15f}")
print(f"  libxc 'GGA_K_VW' n*eps                  : {n_p * exc_p_lib[0]:.15f}")
print(f"  Diff XCFun VWK vs libxc GGA_K_VW        : {n_p * (exc_p_vwk[0] - exc_p_lib[0]):+.3e}")
print(f"  Diff Total-density vs libxc GGA_K_VW    : {e_p_total - n_p * exc_p_lib[0]:+.6e}")
