import numpy as np
from pyscf import gto, dft
from pyscf.pbc import gto as pgto, dft as pdft

# =====================================================================
# 1. Two-center Ne2 cc-pV5Z SR-ERI test (l=5 H shells on two centers!)
# =====================================================================
mol = gto.M(atom="Ne 0 0 0; Ne 0 0 2.5", basis="cc-pv5z", verbose=0, unit="Bohr")
h_sh = [i for i in range(mol.nbas) if mol.bas_angular(i) == 5]
print(f"=== 1. Ne2 cc-pV5Z (two-center l=5 shells {h_sh}) ===")
for shls in [(h_sh[0], h_sh[0], h_sh[1], h_sh[1]), (h_sh[0], h_sh[1], h_sh[0], h_sh[1])]:
    eri_full = mol.intor_by_shell("int2e", shls)
    print(f"  shls={shls}, max|Full|={np.max(np.abs(eri_full)):.3e}")
    for w in [0.11, 0.2, 0.33, 0.4, 1.0, 2.0, 5.0, 10.0, 20.0]:
        with mol.with_range_coulomb(w):
            eri_lr = mol.intor_by_shell("int2e", shls)
        with mol.with_range_coulomb(-w):
            eri_sr = mol.intor_by_shell("int2e", shls)
        err = np.max(np.abs((eri_lr + eri_sr) - eri_full))
        if err > 1e-10:
            print(f"    w={w:5.2f}: ERR max|LR+SR-Full|={err:.3e}, max|LR|={np.max(np.abs(eri_lr)):.3e}, max|SR|={np.max(np.abs(eri_sr)):.3e}")
        else:
            print(f"    w={w:5.2f}: clean ({err:.2e})")

# =====================================================================
# 2. PySCF PBC Skewed Cell RKS grid/Exc diagnostic
# =====================================================================
print("\n=== 2. PySCF PBC Skewed Cell RKS grid/Exc diagnostic ===")
a0 = np.array([[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]])
U2 = np.array([[1, 2, 1], [0, 1, 1], [0, 0, 1]])
for lbl, a in [("orig", a0), ("skew2", U2 @ a0)]:
    cell = pgto.M(atom="Si 0 0 0; Si 1.3575 1.3575 1.3575", a=a, basis="gth-szv", pseudo="gth-pbe", verbose=0)
    mf = pdft.RKS(cell).density_fit().set(xc="pbe", conv_tol=1e-10)
    e = mf.kernel()
    s = mf.scf_summary
    print(f"  {lbl} (default): mesh={cell.mesh}, grids={type(mf.grids).__name__}, ngrids={len(mf.grids.weights)}, E={e:.8f}, Exc={s['exc']:.8f}, Ecoul={s['coul']:.8f}, E1={s['e1']:.8f}")
    mf2 = pdft.RKS(cell).density_fit().set(xc="pbe", conv_tol=1e-10)
    mf2.grids = pdft.gen_grid.BeckeGrids(cell)
    mf2.grids.level = 5
    e2 = mf2.kernel()
    s2 = mf2.scf_summary
    print(f"  {lbl} (Becke5) : ngrids={len(mf2.grids.weights)}, E={e2:.8f}, Exc={s2['exc']:.8f}, Ecoul={s2['coul']:.8f}, E1={s2['e1']:.8f}")

# =====================================================================
# 3. TDDFT-HSE06 & LC-wPBE gradient mismatch: grid level 4, 6, 8 + TDA!
# =====================================================================
print("\n=== 3. TDDFT-HSE06 & LC-wPBE gradient vs grid level (CO, def2-svp) ===")
z0 = 2.10
h = 1e-3
for xc in ["hse06", "lc_wpbe", "cam_b3lyp"]:
    for lvl in [4, 6, 8]:
        def run_td(z):
            m = gto.M(atom=f"C 0 0 0; O 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
            mf = dft.RKS(m); mf.xc = xc; mf.grids.level = lvl; mf.conv_tol = 1e-12
            e0 = mf.kernel()
            td = mf.TDDFT(); td.nstates = 3; td.conv_tol = 1e-11
            ex = td.kernel()[0][0]
            return e0 + ex, td
        _, td0 = run_td(z0)
        g_ana = td0.nuc_grad_method().set(grid_response=True).kernel(state=1)[1, 2]
        g_fd = (run_td(z0-2*h)[0] - 8*run_td(z0-h)[0] + 8*run_td(z0+h)[0] - run_td(z0+2*h)[0]) / (12*h)
        print(f"  {xc:10s} grid={lvl}: g_ana={g_ana:+.9f}, g_fd4={g_fd:+.9f}, |diff|={abs(g_ana-g_fd):.3e}", flush=True)
