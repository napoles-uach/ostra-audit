#!/usr/bin/env python3
"""Probe C: Physics-oracle audit of PySCF advanced DFT/Solvent/TDDFT/Non-collinear GKS engines:
  1. Non-collinear GKS ('mcol') SO(3) spin-quantization-axis rotational invariance!
  2. VV10 non-local correlation (nlc='vv10', wB97M-V, r2SCAN+VV10):
     (a) polarized (rho_a == rho_b) == unpolarized,
     (b) analytical nuclear gradient vs 4-pt FD!
  3. Implicit solvation (ddCOSMO, PCM: IEF-PCM, C-PCM, COSMO, SS(V)PE):
     analytical nuclear gradient vs 4-pt FD & rotational invariance!
  4. TDDFT / TDA (RKS/UKS with MGGA / RSH / XCFun) & analytical TDDFT gradient vs 4-pt FD!
"""
import math
import numpy as np
from pyscf import gto, dft, grad, solvent, tdscf


def test_mcol_spin_rotation_invariance():
    print("=== 1. PySCF Non-Collinear GKS (collinear='mcol') Spin-Rotation Invariance ===")
    # Open-shell O2 triplet (spin=2) or OH radical (spin=1): rotate spin density matrix by SU(2)
    # U(theta, phi) = exp(-i theta/2 sigma_y) on spin block; E_xc MUST be invariant to 1e-12!
    mol = gto.M(atom="O 0 0 0; H 0 0 0.97", basis="def2-svp", spin=1, verbose=0)
    mf_u = dft.UKS(mol); mf_u.xc = "pbe"; mf_u.grids.level = 4; mf_u.conv_tol = 1e-11; mf_u.kernel()
    dma, dmb = mf_u.make_rdm1()
    nao = mol.nao
    # Build GKS spin-diagonal DM along z: [[dma, 0], [0, dmb]]
    dm_z = np.zeros((2*nao, 2*nao), dtype=np.complex128)
    dm_z[:nao, :nao] = dma
    dm_z[nao:, nao:] = dmb

    for xc in ["lda,vwn", "pbe", "r2scan", "tpss", "b3lyp"]:
        mf_g = dft.GKS(mol)
        mf_g.xc = xc
        mf_g.collinear = "mcol"
        mf_g.grids.level = 4
        mf_g.grids.build()
        e_vals = []
        for theta, phi in [(0.0, 0.0), (0.7, 0.0), (math.pi/2, 0.0), (1.1, 0.8), (math.pi, 0.0)]:
            # SU(2) spin rotation matrix
            c, s = math.cos(theta/2), math.sin(theta/2)
            u = np.array([[c, -s*np.exp(-1j*phi)], [s*np.exp(1j*phi), c]])
            U = np.kron(u, np.eye(nao))
            dm_rot = U @ dm_z @ U.conj().T
            vxc = mf_g.get_veff(mol, dm_rot)
            e_vals.append(float(np.real(vxc.exc)))
        spread = max(e_vals) - min(e_vals)
        print(f"  mcol {xc:10s}: E_xc(z)={e_vals[0]:.10f}, max_spin_rot_spread={spread:.3e} Eh  (vals={[round(x, 8) for x in e_vals]})")


def test_vv10_oracles():
    print("\n=== 2. PySCF VV10 (nlc) Oracles: RKS vs UKS(equal spin) & Analytical Gradient vs 4-pt FD ===")
    mol0 = gto.M(atom="He 0 0 0; Ne 0 0 2.8", basis="def2-svp", verbose=0, unit="Bohr")
    for xc, nlc in [("r2scan", "vv10"), ("wb97m_v", "vv10"), ("b97m_v", "vv10")]:
        mf_r = dft.RKS(mol0); mf_r.xc = xc; mf_r.nlc = nlc; mf_r.grids.level = 4; mf_r.nlcgrids.level = 3; mf_r.conv_tol = 1e-11
        e_r = mf_r.kernel()
        mf_u = dft.UKS(mol0); mf_u.xc = xc; mf_u.nlc = nlc; mf_u.grids.level = 4; mf_u.nlcgrids.level = 3; mf_u.conv_tol = 1e-11
        e_u = mf_u.kernel(dm0=(mf_r.make_rdm1()/2, mf_r.make_rdm1()/2))
        g_ana = mf_r.nuc_grad_method().set(grid_response=True).kernel()[1, 2]
        h = 1e-3
        def eval_e(z):
            m = gto.M(atom=f"He 0 0 0; Ne 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
            mf = dft.RKS(m); mf.xc = xc; mf.nlc = nlc; mf.grids.level = 4; mf.nlcgrids.level = 3; mf.conv_tol = 1e-11
            return mf.kernel()
        g_fd = (eval_e(2.8 - 2*h) - 8*eval_e(2.8 - h) + 8*eval_e(2.8 + h) - eval_e(2.8 + 2*h)) / (12.0 * h)
        print(f"  {xc:8s}+{nlc}: RKS-UKS={e_r - e_u:+.3e} | g_ana={g_ana:+.9f}, g_fd4={g_fd:+.9f}, |diff|={abs(g_ana - g_fd):.3e}")


def test_solvent_gradients():
    print("\n=== 3. PySCF Implicit Solvent (ddCOSMO, IEF-PCM, C-PCM, COSMO, SS(V)PE) Gradients vs 4-pt FD ===")
    h = 1e-3  # Bohr
    z0 = 1.85
    for tag, setup in [
        ("ddCOSMO-RHF", lambda m: m.RHF().ddCOSMO().set(conv_tol=1e-11)),
        ("IEFPCM-RHF",  lambda m: m.RHF().PCM().set(conv_tol=1e-11)),
        ("CPCM-RHF",    lambda m: _set_pcm(m.RHF().PCM(), "C-PCM")),
        ("COSMO-RHF",   lambda m: _set_pcm(m.RHF().PCM(), "COSMO")),
        ("SSVPE-RHF",   lambda m: _set_pcm(m.RHF().PCM(), "SS(V)PE")),
        ("IEFPCM-PBE",  lambda m: m.RKS().set(xc="pbe").PCM().set(conv_tol=1e-11)),
    ]:
        def run_e(z):
            m = gto.M(atom=f"H 0 0 0; F 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
            mf = setup(m)
            return mf.kernel(), mf
        try:
            e0, mf0 = run_e(z0)
            g_ana = mf0.nuc_grad_method().kernel()[1, 2]
            g_fd = (run_e(z0 - 2*h)[0] - 8*run_e(z0 - h)[0] + 8*run_e(z0 + h)[0] - run_e(z0 + 2*h)[0]) / (12.0 * h)
            print(f"  {tag:13s}: g_ana={g_ana:+.9f}, g_fd4={g_fd:+.9f}, |diff|={abs(g_ana - g_fd):.3e}")
        except Exception as exc:
            print(f"  {tag:13s}: error ({exc})")


def _set_pcm(mf, method):
    mf.with_solvent.method = method
    mf.conv_tol = 1e-11
    return mf


def test_tddft_gradients():
    print("\n=== 4. PySCF TDDFT / TDA Gradients (RSH / MGGA / XCFun) vs 4-pt FD ===")
    h = 1e-3  # Bohr
    z0 = 2.10
    for tag, xc, use_tda in [
        ("TDA-PBE",    "pbe",    True),
        ("TDDFT-PBE",  "pbe",    False),
        ("TDDFT-HSE06","hse06",  False),
        ("TDDFT-wB97X","wb97x",  False),
        ("TDDFT-r2SCAN","r2scan",False),
    ]:
        def run_e(z):
            m = gto.M(atom=f"C 0 0 0; O 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
            mf = dft.RKS(m); mf.xc = xc; mf.grids.level = 4; mf.conv_tol = 1e-11
            e0 = mf.kernel()
            td = mf.TDA() if use_tda else mf.TDDFT()
            td.nstates = 3; td.conv_tol = 1e-9
            ex = td.kernel()[0][0]
            return e0 + ex, td
        try:
            _, td0 = run_e(z0)
            g_ana = td0.nuc_grad_method().set(grid_response=True).kernel(state=1)[1, 2]
            g_fd = (run_e(z0 - 2*h)[0] - 8*run_e(z0 - h)[0] + 8*run_e(z0 + h)[0] - run_e(z0 + 2*h)[0]) / (12.0 * h)
            print(f"  {tag:14s}: g_ana={g_ana:+.9f}, g_fd4={g_fd:+.9f}, |diff|={abs(g_ana - g_fd):.3e}")
        except Exception as exc:
            print(f"  {tag:14s}: error ({exc})")


if __name__ == "__main__":
    # test_mcol_spin_rotation_invariance()  # requires optional mcfun package
    test_vv10_oracles()
    test_solvent_gradients()
    test_tddft_gradients()
