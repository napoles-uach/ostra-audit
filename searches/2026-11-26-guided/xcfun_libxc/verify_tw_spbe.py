# Independent checks for two XCFun-vs-libxc discrepancies found in the guided differential search.
import numpy as np
from pyscf import gto, dft
from pyscf.dft import xcfun, libxc
# --- TW (von Weizsacker kinetic energy): exact for a one-orbital density, T_W = T_s. H atom: T = 0.5 Eh.
mol = gto.M(atom="H 0 0 0", spin=1, basis="aug-cc-pv5z", verbose=0)
mf = dft.UKS(mol); mf.xc = "hf"; mf.kernel()
g = dft.gen_grid.Grids(mol); g.level = 6; g.build()
ni = dft.numint.NumInt(); ao = ni.eval_ao(mol, g.coords, deriv=1)
dm = mf.make_rdm1(); rho = np.array([ni.eval_rho(mol, ao, d, xctype="GGA") for d in dm])
n = rho[0, 0] + rho[1, 0]
T_exact = mf.mol.intor("int1e_kin").__mul__(dm[0]).sum()
print(f"H atom: <T> of the orbital = {T_exact:.6f} Eh (exact 0.5)")
for mod, name in ((libxc, "GGA_K_VW"), (xcfun, "TW")):
    e = mod.eval_xc(name, rho, spin=1, deriv=0)[0]
    print(f"   {mod.__name__.split('.')[-1]:6s} {name:9s} T_W = {np.sum(g.weights * n * e):.6e} Eh")
# --- SPBEC vs GGA_C_SPBE: rebuild sPBE correlation from the paper's formula (Swart, Sola, Bickelhaupt 2009),
#     PBE constants beta = 0.066725, gamma = (1 - ln 2)/pi^2, H = gamma phi^3 ln[1 + (beta/gamma) t^2 / (1 + A t^2)],
#     A = (beta/gamma) / (exp(-eps/(gamma phi^3)) - 1), with VWN5 eps (as XCFun) and PW92 eps (as PBE).
mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469", basis="def2-svp", verbose=0)
mf = dft.RKS(mol); mf.xc = "pbe"; mf.kernel()
g = mf.grids; ao = ni.eval_ao(mol, g.coords, deriv=1); rho = ni.eval_rho(mol, ao, mf.make_rdm1(), xctype="GGA")
n = rho[0]; s2 = (rho[1:4] ** 2).sum(0); m = n > 1e-10
beta, gamma = 0.066725, (1 - np.log(2)) / np.pi ** 2
kF = (3 * np.pi ** 2 * n[m]) ** (1 / 3); ks = np.sqrt(4 * kF / np.pi); t2 = s2[m] / (2 * ks * n[m]) ** 2
for lda in ("LDA_C_VWN", "LDA_C_PW"):
    eps = libxc.eval_xc(lda, n[m], spin=0, deriv=0)[0]
    A = (beta / gamma) / np.expm1(-eps / gamma)
    H = gamma * np.log(1 + (beta / gamma) * t2 / (1 + A * t2))
    print(f"H2O sPBE correlation from the paper formula ({lda} eps): {np.sum(g.weights[m] * n[m] * (eps + H)):.6f} Eh")
for mod, name in ((libxc, "GGA_C_SPBE"), (xcfun, "SPBEC"), (libxc, "GGA_C_PBE")):
    e = mod.eval_xc(name, rho, spin=0, deriv=0)[0]
    print(f"   {mod.__name__.split('.')[-1]:6s} {name:11s} E_c = {np.sum(g.weights * n * e):.6f} Eh")
