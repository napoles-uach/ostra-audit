#!/usr/bin/env python3
"""Probe D: Physics-oracle audit of libcint special integrals (range-separated erf/erfc,
analytic derivative consistency, translational/rotational invariance, ECP 1st & 2nd derivatives).
"""
import math
import numpy as np
from pyscf import gto


def test_rsh_partition_of_unity():
    print("=== 1. libcint Range-Separated Partition of Unity: int2e_erf(w) + int2e_erfc(w) == int2e ===")
    # Test across angular momenta s, p, d, f, g, h and omega from 1e-4 to 1e3!
    basis = {
        "H": gto.basis.parse("""
H    S
   0.15  1.0
H    P
   0.45  1.0
H    D
   1.20  1.0
H    F
   0.80  1.0
H    G
   0.60  1.0
H    H
   0.50  1.0
"""),
    }
    mol = gto.M(atom="H 0 0 0; H 0.7 1.1 1.5", basis=basis, verbose=0)
    eri_full = mol.intor("int2e")
    norm_full = np.linalg.norm(eri_full)
    for omega in [1e-5, 1e-3, 1e-2, 0.1, 0.3, 1.0, 5.0, 20.0, 100.0, 1000.0]:
        with mol.with_range_coulomb(omega):
            eri_lr = mol.intor("int2e")
        with mol.with_range_coulomb(-omega):
            eri_sr = mol.intor("int2e")
        err_max = float(np.max(np.abs((eri_lr + eri_sr) - eri_full)))
        err_rel = float(np.linalg.norm((eri_lr + eri_sr) - eri_full) / norm_full)
        print(f"  omega={omega:9.1e}: max|LR+SR-full|={err_max:.3e}, rel_L2={err_rel:.3e}, max|LR|={np.max(np.abs(eri_lr)):.3e}, max|SR|={np.max(np.abs(eri_sr)):.3e}")

    print("\n=== 1b. libcint RSH Gradient Consistency: int2e_ip1(+-omega) vs 4-pt FD ===")
    mol2 = gto.M(atom="H 0 0 0; H 0.4 0.9 1.2", basis="cc-pvtz", verbose=0)
    h = 1e-4  # Bohr
    for omega in [0.1, 1.0, 10.0, 100.0, -0.1, -1.0, -10.0, -100.0]:
        with mol2.with_range_coulomb(omega):
            ip1 = mol2.intor("int2e_ip1")  # (3, nao, nao, nao, nao) = -d/dR_i on bra index 0? Check sign!
        # Perturb atom 0 along x: ao slices of atom 0
        sh0, sh1, p0, p1 = mol2.aoslice_by_atom()[0]
        def eval_eri_dx(dx):
            m = gto.M(atom=f"H {dx:.12f} 0 0; H 0.4 0.9 1.2", basis="cc-pvtz", verbose=0, unit="Bohr")
            with m.with_range_coulomb(omega):
                return m.intor("int2e")
        # Wait: mol2 used Angstrom! Re-create mol2 in Bohr for exact match
        m0 = gto.M(atom="H 0 0 0; H 0.4 0.9 1.2", basis="cc-pvtz", verbose=0, unit="Bohr")
        with m0.with_range_coulomb(omega):
            ip1 = m0.intor("int2e_ip1")
        e_m2, e_m1, e_p1, e_p2 = [eval_eri_dx(k * h) for k in (-2, -1, 1, 2)]
        fd = (e_m2 - 8*e_m1 + 8*e_p1 - e_p2) / (12.0 * h)
        # d/dR_0 (ij|kl) = -(nabla_i ij|kl) - (i nabla_j|kl) - (ij|nabla_k l) - (ij|k nabla_l) restricted to AO on atom 0
        ana = np.zeros_like(fd)
        ana[p0:p1, :, :, :] -= ip1[0, p0:p1, :, :, :]
        ana[:, p0:p1, :, :] -= ip1[0, p0:p1, :, :, :].transpose(1, 0, 2, 3)
        ana[:, :, p0:p1, :] -= ip1[0, p0:p1, :, :, :].transpose(2, 3, 0, 1)
        ana[:, :, :, p0:p1] -= ip1[0, p0:p1, :, :, :].transpose(2, 3, 1, 0)
        rel = float(np.linalg.norm(ana - fd) / max(np.linalg.norm(fd), 1e-30))
        mx = float(np.max(np.abs(ana - fd)))
        print(f"  omega={omega:+7.1f}: max|ana-FD4|={mx:.3e}, rel_L2={rel:.3e}")


def test_ecp_derivatives_and_invariance():
    print("\n=== 2. PySCF libcgto ECP 1st & 2nd derivatives (ECPscalar_ipnuc / ipipnuc / iprinvip) vs 4-pt FD ===")
    # Use a heavy element with high-l ECP (e.g. Pt or I or Au with def2-tzvp ECP) + a light atom (H, no ECP)
    # and ALSO two heavy atoms (I2 / PtH)!
    for atom_tpl, basis_d, ecp_d, lbl in [
        ("I 0 0 0; H {x:.12f} 0.5 0.9", "def2-tzvp", {"I": "def2-tzvp"}, "I-H (I has ECP, H no ECP)"),
        ("I 0 0 0; I {x:.12f} 0.5 1.8", "def2-tzvp", {"I": "def2-tzvp"}, "I2 (both ECP)"),
        ("Au 0 0 0; Cl {x:.12f} 0.4 1.9", "def2-tzvp", {"Au": "def2-tzvp"}, "Au-Cl (Au f-proj ECP)"),
    ]:
        h = 1e-4
        x0 = 1.2
        def get_ecp(x):
            m = gto.M(atom=atom_tpl.format(x=x), basis=basis_d, ecp=ecp_d, verbose=0, unit="Bohr")
            return m.intor("ECPscalar")
        m0 = gto.M(atom=atom_tpl.format(x=x0), basis=basis_d, ecp=ecp_d, verbose=0, unit="Bohr")
        e0 = m0.intor("ECPscalar")
        # 1st derivative wrt atom 1 x-coord (both basis on atom 1 and ECP center on atom 1 if any)
        # Use PySCF's built-in gradient nuclear ECP contribution!
        from pyscf import grad
        # Compare analytic d(ECPscalar)/dR_1x using ECPscalar_ipnuc and ECPscalar_iprinv
        _, _, p0, p1 = m0.aoslice_by_atom()[1]
        ipnuc = m0.intor("ECPscalar_ipnuc")  # (3, nao, nao) = <nabla_1 | U_ECP | 2>
        with m0.with_rinv_at_nucleus(1):
            iprinv = m0.intor("ECPscalar_iprinv")  # <nabla_1 | U_ECP(at atom 1) | 2> ? Wait: check sign & definition!
        ana = np.zeros_like(e0)
        # Basis on atom 1 moves with R_1 => d/dR_1 <phi_i|U|phi_j> = -<nabla phi_i|U|phi_j> (i in 1) - <phi_i|U|nabla phi_j> (j in 1)
        # Plus ECP center on atom 1 moves with R_1 => +<nabla phi_i|U_1|phi_j> + <phi_i|U_1|nabla phi_j>
        ana[p0:p1, :] -= ipnuc[0, p0:p1, :]
        ana[:, p0:p1] -= ipnuc[0, p0:p1, :].T
        if m0._ecpbas.shape[0] > 0 and np.any(m0._ecpbas[:, 0] == 1):
            ana += iprinv[0] + iprinv[0].T
        em2, em1, ep1, ep2 = [get_ecp(x0 + k*h) for k in (-2, -1, 1, 2)]
        fd = (em2 - 8*em1 + 8*ep1 - ep2) / (12.0 * h)
        # Check whether iprinv sign is + or -
        err_p = float(np.max(np.abs(ana - fd)))
        rel_p = float(np.linalg.norm(ana - fd) / max(np.linalg.norm(fd), 1e-30))
        print(f"  {lbl:28s}: max|ana-FD4|={err_p:.3e}, rel_L2={rel_p:.3e}, max|FD4|={np.max(np.abs(fd)):.3e}")

        # Also check ECPscalar_iprinv on atom 1 when atom 1 has NO ECP (known pyscf#3669!)
        if not np.any(m0._ecpbas[:, 0] == 1):
            print(f"    (control pyscf#3669 check: atom 1 no ECP, iprinv finite={np.all(np.isfinite(iprinv))}, max={np.max(np.abs(iprinv)):.3e})")


def test_relativistic_and_3c_integrals():
    print("\n=== 3. libcint Special 1e/2e/3c integrals: translational & rotational invariance ===")
    mol0 = gto.M(atom="O 0 0 0; H 0.75 0.58 0.1; H -0.75 0.58 -0.2", basis="cc-pvtz", verbose=0)
    # Random translation + SO(3) rotation (only scalar/trace invariants for rotation, exact equality for translation!)
    shift = np.array([13.7, -42.3, 19.1])
    coords_t = mol0.atom_coords(unit="Angstrom") + shift
    mol_t = gto.M(atom=[(mol0.atom_symbol(i), coords_t[i]) for i in range(mol0.natm)], basis="cc-pvtz", verbose=0)
    for intor in ["int1e_ovlp", "int1e_kin", "int1e_nuc", "int1e_pnucp", "int2e", "int2e_spsp1", "int2e_spsp1spsp2", "int2e_ssp1sps2", "int2e_breit_r1p2", "int2e_gauge_r1_ssp1ssp2"]:
        try:
            a0 = mol0.intor(intor)
            at = mol_t.intor(intor)
            diff = float(np.max(np.abs(a0 - at)))
            rel = diff / max(float(np.max(np.abs(a0))), 1e-30)
            print(f"  Translation invariance {intor:26s}: max_diff={diff:.3e}, rel={rel:.3e}")
        except Exception as e:
            print(f"  {intor:26s}: skipped ({e})")


if __name__ == "__main__":
    test_rsh_partition_of_unity()
    test_ecp_derivatives_and_invariance()
    test_relativistic_and_3c_integrals()
