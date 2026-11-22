import numpy as np, sys, io, contextlib
sys.path.insert(0, ".")
from dftd4.interface import DampingParam as D4Param, DispersionModel as D4Model
from dftd3.interface import RationalDampingParam as D3BJ, DispersionModel as D3Model
src = open("probe_disp.py").read()
ns = {}
exec(src[: src.index("# ---------------- molecular oracles")], ns)
fd_grad, fd_strain, rel = ns["fd_grad"], ns["fd_strain"], ns["rel"]
P3 = np.array([True] * 3)

def d4(cut=None, width=0.0):
    def f(n, p, lat=None):
        m = D4Model(n, p, lattice=lat, periodic=P3)
        if cut: m.set_realspace_cutoff(cut[0], cut[1], cut[2], width, width)
        return m.get_dispersion(D4Param(method="pbe0"), grad=True)
    return f

def d3(cut=None, width=0.0):
    def f(n, p, lat=None):
        m = D3Model(n, p, lattice=lat, periodic=P3)
        if cut: m.set_realspace_cutoff(cut[0], cut[1], cut[2], width, width)
        return m.get_dispersion(D3BJ(method="b3lyp"), grad=True)
    return f

print("== A. absolute net force in high-symmetry crystals (relative metric was dividing by ~0)")
for name, (n, p, lat) in (("Si", ns["si_diamond"]()), ("NaCl", ns["nacl"]())):
    for mname, f in (("D4", d4()), ("D3BJ", d3())):
        G = f(n, p, lat)["gradient"]
        print(f"  {name} {mname}: max|G| = {np.abs(G).max():.2e} Eh/bohr, |sum G| = {np.abs(G.sum(0)).max():.2e}")

print("== C. D4 partial charges and energy: cell choice and supercell (D3 has no charges)")
T = np.array([[1, 1, 0], [0, 1, 0], [-1, 0, 1]])
for name, (n, p, lat) in (("NaCl", ns["nacl"]()), ("triclinic", ns["triclinic_mixed"]())):
    q0 = D4Model(n, p, lattice=lat, periodic=P3).get_properties()["partial charges"]
    q1 = D4Model(n, p, lattice=T @ lat, periodic=P3).get_properties()["partial charges"]
    shifts = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)]) @ lat
    nsc = np.tile(n, 8); psc = (p[None] + shifts[:, None]).reshape(-1, 3)
    qs = D4Model(nsc, psc, lattice=2 * lat, periodic=P3).get_properties()["partial charges"].reshape(8, -1)
    print(f"  {name}: charges {np.round(q0, 5)}")
    print(f"     equivalent cell: max|dq| = {np.abs(q1 - q0).max():.2e}   supercell: max|dq| = {np.abs(qs - q0[None]).max():.2e}")
    cn0 = D4Model(n, p, lattice=lat, periodic=P3).get_properties()["coordination numbers"]
    cn1 = D4Model(n, p, lattice=T @ lat, periodic=P3).get_properties()["coordination numbers"]
    print(f"     coordination numbers, equivalent cell: max|dCN| = {np.abs(cn1 - cn0).max():.2e}")
