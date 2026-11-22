import numpy as np
from pyscf import gto, dft

print("=== Checking RKS vs TDDFT analytical gradient vs FD4 (CO, def2-svp, grid level 6) ===")
z0 = 2.10
h = 1e-3

for tag, xc in [
    ("PBE", "pbe"),
    ("PBE0", "pbe0"),
    ("CAM-B3LYP", "cam_b3lyp"),
    ("LC-wPBE", "lc_wpbe"),
    ("HSE06", "hse06"),
]:
    def run_gs(z):
        m = gto.M(atom=f"C 0 0 0; O 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
        mf = dft.RKS(m); mf.xc = xc; mf.grids.level = 6; mf.conv_tol = 1e-12
        e0 = mf.kernel()
        return e0, mf

    e0, mf0 = run_gs(z0)
    g_gs_ana = mf0.nuc_grad_method().set(grid_response=True).kernel()[1, 2]
    g_gs_fd = (run_gs(z0-2*h)[0] - 8*run_gs(z0-h)[0] + 8*run_gs(z0+h)[0] - run_gs(z0+2*h)[0]) / (12*h)

    def run_td(z):
        e, mf = run_gs(z)
        td = mf.TDDFT(); td.nstates = 3; td.conv_tol = 1e-11
        ex = td.kernel()[0][0]
        return e + ex, td

    _, td0 = run_td(z0)
    g_td_ana = td0.nuc_grad_method().set(grid_response=True).kernel(state=1)[1, 2]
    g_td_fd = (run_td(z0-2*h)[0] - 8*run_td(z0-h)[0] + 8*run_td(z0+h)[0] - run_td(z0+2*h)[0]) / (12*h)

    print(f"  {tag:10s} | GS |ana-FD4|={abs(g_gs_ana-g_gs_fd):.3e} | TD |ana-FD4|={abs(g_td_ana-g_td_fd):.3e} (ana={g_td_ana:+.8f}, fd={g_td_fd:+.8f})", flush=True)
