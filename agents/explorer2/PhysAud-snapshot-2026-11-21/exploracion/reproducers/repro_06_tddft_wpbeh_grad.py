#!/usr/bin/env python3
"""Minimal reproducer 06: Observable impact of Libxc s_scaling_2 discontinuity at s=1
on PySCF TDDFT excited-state analytical gradients (HSE06, LC-wPBE, GGA_X_WPBEH) vs
smooth controls (CAM-B3LYP, GGA_X_HJS_PBE, PBE).
On CO (def2-SVP), CAM-B3LYP converges to 2.8e-9 Eh/Bohr at grid level 8, while
HSE06 (3.4e-5) and LC-wPBE (3.9e-5) stall and oscillate 12,000x-14,000x worse!
"""
import numpy as np
from pyscf import gto, dft

z0 = 2.10
h = 1e-3

def check_td_grad(xc, grid_level=6):
    def run_td(z):
        m = gto.M(atom=f"C 0 0 0; O 0 0 {z:.12f}", basis="def2-svp", verbose=0, unit="Bohr")
        mf = dft.RKS(m); mf.xc = xc; mf.grids.level = grid_level; mf.conv_tol = 1e-12
        e0 = mf.kernel()
        td = mf.TDDFT(); td.nstates = 3; td.conv_tol = 1e-11
        ex = td.kernel()[0][0]
        return e0 + ex, td
    _, td0 = run_td(z0)
    g_ana = td0.nuc_grad_method().set(grid_response=True).kernel(state=1)[1, 2]
    g_fd4 = (run_td(z0-2*h)[0] - 8*run_td(z0-h)[0] + 8*run_td(z0+h)[0] - run_td(z0+2*h)[0]) / (12*h)
    return g_ana, g_fd4, abs(g_ana - g_fd4)

print("=== TDDFT analytical vs FD4 gradient (CO, def2-SVP, h=1e-3 Bohr) ===")
for xc, lvl in [
    ("cam_b3lyp", 6),
    ("GGA_X_PBE,GGA_C_PBE", 6),
    ("GGA_X_HJS_PBE,GGA_C_PBE", 6),
    ("GGA_X_WPBEH,GGA_C_PBE", 6),
    ("hse06", 6),
    ("lc_wpbe", 6),
]:
    ga, gf, d = check_td_grad(xc, lvl)
    print(f"  {xc:25s} (grid={lvl}): g_ana={ga:+.9f}, g_fd4={gf:+.9f}, |diff|={d:.3e}", flush=True)
