#!/usr/bin/env python3
"""Probe B: Physics-oracle audit of PySCF PBC (periodic boundary conditions):
  1. Unimodular cell equivalence (skewed lattice a' = U a, det U = 1) for Ewald, FFTDF, GDF, RSDF, MDF!
  2. Supercell invariance (2x1x1 supercell energy == 2 * primitive cell energy).
  3. Fractional coordinate translation invariance (r_i -> r_i + t).
  4. PBC analytical nuclear gradients vs 4-pt finite differences.
"""
import numpy as np
from pyscf.pbc import gto as pgto, scf as pscf, dft as pdft, df as pdf


def test_ewald_and_scf_skewed_cells():
    print("=== 1. PySCF PBC Unimodular Cell Equivalence (a' = U @ a, det(U)=1) ===")
    a0 = np.array([[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]])  # FCC Si primitive (Angstrom)
    # Unimodular matrices with increasing skewness (det = +1)
    Us = [
        ("orig", np.eye(3)),
        ("skew1", np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])),
        ("skew2", np.array([[1, 2, 1], [0, 1, 1], [0, 0, 1]])),
        ("skew3", np.array([[1, 3, 2], [0, 1, 2], [0, 0, 1]])),
        ("skew4", np.array([[1, 4, 3], [0, 1, 3], [0, 0, 1]])),
    ]
    e_ewald_ref = None
    e_gdf_ref = None
    e_rsdf_ref = None
    e_fft_ref = None
    for lbl, U in Us:
        a = U @ a0
        cell = pgto.M(
            atom="Si 0 0 0; Si 1.3575 1.3575 1.3575",
            a=a,
            basis="gth-szv",
            pseudo="gth-pbe",
            verbose=0,
        )
        ew = cell.energy_nuc()
        if e_ewald_ref is None:
            e_ewald_ref = ew
        # Test FFTDF, GDF, RSDF PBE at Gamma point
        mf_fft = pdft.RKS(cell); mf_fft.xc = "pbe"; mf_fft.conv_tol = 1e-10; e_fft = mf_fft.kernel()
        mf_gdf = pdft.RKS(cell).density_fit(); mf_gdf.xc = "pbe"; mf_gdf.conv_tol = 1e-10; e_gdf = mf_gdf.kernel()
        mf_rsdf = pdft.RKS(cell).rs_density_fit(); mf_rsdf.xc = "pbe"; mf_rsdf.conv_tol = 1e-10; e_rsdf = mf_rsdf.kernel()
        if e_fft_ref is None:
            e_fft_ref, e_gdf_ref, e_rsdf_ref = e_fft, e_gdf, e_rsdf
        print(f"  {lbl:6s}: dE_ewald={ew - e_ewald_ref:+.3e} | dE_FFT={e_fft - e_fft_ref:+.3e} | "
              f"dE_GDF={e_gdf - e_gdf_ref:+.3e} | dE_RSDF={e_rsdf - e_rsdf_ref:+.3e}")


def test_pbc_translation_and_supercell():
    print("\n=== 2. PySCF PBC Translation & Supercell Invariance (LiH / Si, GDF & RSDF & FFTDF) ===")
    a = np.diag([4.0, 4.0, 4.0])
    c0 = pgto.M(atom="Li 0 0 0; H 2.0 2.0 2.0", a=a, basis="gth-szv", pseudo="gth-pbe", verbose=0)
    # Fractional shift (0.17, 0.31, 0.43)
    sh = np.array([0.17, 0.31, 0.43]) * 4.0
    c_sh = pgto.M(atom=f"Li {sh[0]} {sh[1]} {sh[2]}; H {2+sh[0]} {2+sh[1]} {2+sh[2]}", a=a, basis="gth-szv", pseudo="gth-pbe", verbose=0)
    # 2x1x1 supercell
    a2 = np.diag([8.0, 4.0, 4.0])
    c2 = pgto.M(
        atom="Li 0 0 0; H 2 2 2; Li 4 0 0; H 6 2 2",
        a=a2, basis="gth-szv", pseudo="gth-pbe", verbose=0,
    )
    for tag, mk in [
        ("FFTDF-PBE", lambda c, k: pdft.KRKS(c, kpts=c.make_kpts(k)).set(xc="pbe", conv_tol=1e-10)),
        ("GDF-PBE",   lambda c, k: pdft.KRKS(c, kpts=c.make_kpts(k)).density_fit().set(xc="pbe", conv_tol=1e-10)),
        ("RSDF-PBE",  lambda c, k: pdft.KRKS(c, kpts=c.make_kpts(k)).rs_density_fit().set(xc="pbe", conv_tol=1e-10)),
        ("GDF-RHF",   lambda c, k: pscf.KRHF(c, kpts=c.make_kpts(k)).density_fit().set(conv_tol=1e-10)),
    ]:
        e0 = mk(c0, [2, 1, 1]).kernel()
        e_sh = mk(c_sh, [2, 1, 1]).kernel()
        e_sc = mk(c2, [1, 1, 1]).kernel()
        print(f"  {tag:10s}: E0={e0:.10f} | shift_err={e_sh - e0:+.3e} | supercell(E2-2*E0)={e_sc - 2*e0:+.3e}")


def test_pbc_gradients():
    print("\n=== 3. PySCF PBC Analytical Nuclear Gradients vs 4-pt FD ===")
    a = np.diag([4.2, 4.2, 4.2])
    x0 = 1.9
    h = 1e-3  # Bohr
    for tag, method in [
        ("GDF-RHF", lambda c: pscf.RHF(c).density_fit().set(conv_tol=1e-11)),
        ("GDF-PBE", lambda c: pdft.RKS(c).density_fit().set(xc="pbe", conv_tol=1e-11)),
        ("GDF-r2SCAN", lambda c: pdft.RKS(c).density_fit().set(xc="r2scan", conv_tol=1e-11)),
        ("GDF-HSE06", lambda c: pdft.RKS(c).density_fit().set(xc="hse06", conv_tol=1e-11)),
    ]:
        def run_e(x_bohr):
            c = pgto.M(atom=f"Li 0 0 0; H {x_bohr:.12f} 0.3 0.2", a=a, basis="gth-szv", pseudo="gth-pbe", unit="Bohr", verbose=0)
            mf = method(c)
            return mf.kernel(), mf
        try:
            e0, mf0 = run_e(x0)
            g_ana = mf0.nuc_grad_method().kernel()[1, 0]
            em2, _ = run_e(x0 - 2*h)
            em1, _ = run_e(x0 - h)
            ep1, _ = run_e(x0 + h)
            ep2, _ = run_e(x0 + 2*h)
            g_fd = (em2 - 8*em1 + 8*ep1 - ep2) / (12.0 * h)
            print(f"  {tag:12s}: g_ana={g_ana:+.10f}, g_fd4={g_fd:+.10f}, |diff|={abs(g_ana-g_fd):.3e}")
        except Exception as exc:
            print(f"  {tag:12s}: error ({exc})")


if __name__ == "__main__":
    test_ewald_and_scf_skewed_cells()
    test_pbc_translation_and_supercell()
    test_pbc_gradients()
