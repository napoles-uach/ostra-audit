"""Does the high-l ERI precision loss in libcint change HF energies with a real basis?
Oracle: total energy is invariant under a rigid rotation of the molecule.
Step 1: SCF in orientation A -> D_A, E_A.
Step 2: rotate D_A exactly to orientation B (Wigner matrices of real spherical harmonics),
        evaluate E_B = E[D_B] with integrals computed in orientation B. dE = E_B - E_A
        measures the effect of integral errors on the energy, weighted by a converged density.
Repeated with the i (l=6) shells removed from the basis to attribute the error.
"""
import sys, time
import numpy as np
from pyscf import gto, scf, lib
from pyscf.symm.sph import real_sph_vec
from scipy.spatial.transform import Rotation

lib.num_threads(2)
basis_name = sys.argv[1] if len(sys.argv) > 1 else "cc-pv6z"
drop_l = int(sys.argv[2]) if len(sys.argv) > 2 else -1
geom_name = sys.argv[3] if len(sys.argv) > 3 else "N2"

geoms = {
    "N2": [["N", (0, 0, 0)], ["N", (0, 0, 2.074)]],
    "N2s": [["N", (0, 0, 0)], ["N", (0, 0, 4.5)]],
    "CO2": [["C", (0, 0, 0)], ["O", (0, 0, 2.196)], ["O", (0, 0, -2.196)]],
    "H2O": [["O", (0, 0, 0)], ["H", (0, 1.431, 1.108)], ["H", (0, -1.431, 1.108)]],
}

def build(coords, symbols):
    b = gto.basis.load(basis_name, symbols[0]) if False else None
    bas = {}
    for s in set(symbols):
        shells = gto.basis.load(basis_name, s)
        if drop_l >= 0:
            shells = [sh for sh in shells if sh[0] < drop_l]
        bas[s] = shells
    m = gto.Mole(atom=[[s, c] for s, c in zip(symbols, coords)], basis=bas, unit="bohr", verbose=0, max_memory=5500)
    m.build()
    return m

symbols = [a[0] for a in geoms[geom_name]]
cA = np.array([a[1] for a in geoms[geom_name]], float)
cA = cA - cA.mean(0) + np.array([0.11, -0.07, 0.05])
Rm = Rotation.random(random_state=2026).as_matrix()
cB = cA @ Rm.T
mA, mB = build(cA, symbols), build(cB, symbols)
print(f"{geom_name} {basis_name} drop_l={drop_l}: nao={mA.nao}, max l={max(mA.bas_angular(i) for i in range(mA.nbas))}", flush=True)

# AO rotation matrix: phi_B(r) = phi_A(R^T (r - ...)) ; build per shell from real spherical harmonics
pts = np.random.default_rng(0).normal(size=(400, 3)); pts /= np.linalg.norm(pts, axis=1)[:, None]
Dl = {}
for l in range(7):
    Y = np.array(real_sph_vec(pts, l, True)[l]).T
    Yr = np.array(real_sph_vec(pts @ Rm.T, l, True)[l]).T
    Dl[l] = np.linalg.lstsq(Yr, Y, rcond=None)[0]
U = np.zeros((mA.nao, mA.nao))
off = 0
for ib in range(mA.nbas):
    l, nc = mA.bas_angular(ib), mA.bas_nctr(ib)
    d = 2 * l + 1
    for c in range(nc):
        U[off:off + d, off:off + d] = Dl[l]; off += d
SA, SB = mA.intor("int1e_ovlp"), mB.intor("int1e_ovlp")
for cand in (U, U.T):
    if np.abs(cand.T @ SA @ cand - SB).max() < 1e-10:
        U = cand; break
print(f"  overlap check |U^T S_A U - S_B| = {np.abs(U.T @ SA @ U - SB).max():.2e}", flush=True)

t = time.time()
mfA = scf.RHF(mA); mfA.conv_tol = 1e-11; mfA.direct_scf_tol = 1e-14; mfA.max_cycle = 60
EA = mfA.kernel()
print(f"  SCF A: E = {EA:.12f}  converged={mfA.converged}  ({time.time() - t:.0f} s)", flush=True)
DA = mfA.make_rdm1()
DB = U.T @ DA @ U
mfB = scf.RHF(mB); mfB.direct_scf_tol = 1e-14
t = time.time()
EB = mfB.energy_tot(DB)
hA = mfA.get_hcore(); hB = mfB.get_hcore()
J_A, K_A = mfA.get_jk(mA, DA); J_B, K_B = mfB.get_jk(mB, DB)
e1A, e1B = np.einsum("ij,ji", hA, DA), np.einsum("ij,ji", hB, DB)
eJA, eJB = 0.5 * np.einsum("ij,ji", J_A, DA), 0.5 * np.einsum("ij,ji", J_B, DB)
eKA, eKB = -0.25 * np.einsum("ij,ji", K_A, DA), -0.25 * np.einsum("ij,ji", K_B, DB)
print(f"  E[D_B] in orientation B = {EB:.12f}   ({time.time() - t:.0f} s)")
print(f"  dE total = {EB - EA:+.3e} Eh   (one-electron {e1B - e1A:+.3e}, Coulomb {eJB - eJA:+.3e}, exchange {eKB - eKA:+.3e})")
