"""Physics oracles for dftd4 (D4) and simple-dftd3 (D3-BJ / D3-zero): molecules and periodic systems."""
import numpy as np
from scipy.spatial.transform import Rotation
from dftd4.interface import DampingParam as D4Param, DispersionModel as D4Model
from dftd3.interface import RationalDampingParam as D3BJ, ZeroDampingParam as D3Zero, DispersionModel as D3Model

rng = np.random.default_rng(5)

# ---------------- systems (bohr) ----------------
def benzene_dimer():
    ang = np.arange(6) * np.pi / 3
    c = np.c_[2.64 * np.cos(ang), 2.64 * np.sin(ang), np.zeros(6)]
    h = np.c_[4.68 * np.cos(ang), 4.68 * np.sin(ang), np.zeros(6)]
    mono = np.r_[c, h]; num = np.r_[[6] * 6, [1] * 6]
    pos = np.r_[mono, mono + [1.2, 0.3, 6.8]]
    return np.r_[num, num], pos

def random_cluster(n=24):
    elems = np.array([1, 6, 7, 8, 9, 15, 16, 17, 35, 26, 53, 78])
    num = rng.choice(elems, n)
    pos = []
    while len(pos) < n:
        p = rng.uniform(-7, 7, 3)
        if all(np.linalg.norm(p - q) > 2.4 for q in pos):
            pos.append(p)
    return num, np.array(pos)

def si_diamond():
    a = 10.26
    lat = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    return np.array([14, 14]), np.array([[0, 0, 0], [0.25, 0.25, 0.25]]) @ lat, lat

def nacl():
    a = 10.66
    lat = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    return np.array([11, 17]), np.array([[0, 0, 0], [0.5, 0.5, 0.5]]) @ lat, lat

def triclinic_mixed():
    lat = np.array([[9.1, 0.0, 0.0], [2.3, 8.4, 0.0], [-1.7, 1.9, 10.2]])
    frac = rng.uniform(0, 1, (6, 3)); num = np.array([6, 8, 1, 1, 16, 7])
    return num, frac @ lat, lat

# ---------------- wrappers ----------------
METHODS = {
    "D4/PBE0": lambda n, p, lat=None: D4Model(n, p, lattice=lat, periodic=None if lat is None else np.array([True] * 3)).get_dispersion(D4Param(method="pbe0"), grad=True),
    "D3BJ/B3LYP": lambda n, p, lat=None: D3Model(n, p, lattice=lat, periodic=None if lat is None else np.array([True] * 3)).get_dispersion(D3BJ(method="b3lyp"), grad=True),
    "D3zero/PBE": lambda n, p, lat=None: D3Model(n, p, lattice=lat, periodic=None if lat is None else np.array([True] * 3)).get_dispersion(D3Zero(method="pbe"), grad=True),
}

def rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)

def line(tag, val, tol, note=""):
    flag = "  <-- CHECK" if val > tol else ""
    print(f"  {tag:52s} {val:9.2e}{flag} {note}", flush=True)

def fd_grad(f, n, p, lat, h):
    g = np.zeros_like(p)
    for i in range(len(p)):
        for k in range(3):
            e = []
            for s in (-2, -1, 1, 2):
                q = p.copy(); q[i, k] += s * h
                e.append(f(n, q, lat)["energy"])
            g[i, k] = (e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h)
    return g

def fd_strain(f, n, p, lat, h):
    s = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            e = []
            for t in (-2, -1, 1, 2):
                eps = np.eye(3); eps[i, j] += t * h
                e.append(f(n, p @ eps.T, lat @ eps.T)["energy"])
            s[i, j] = (e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h)
    return s

# ---------------- molecular oracles ----------------
for sysname, (num, pos) in (("benzene dimer", benzene_dimer()), ("random 24-atom cluster", random_cluster())):
    for mname, f in METHODS.items():
        print(f"--- {sysname} / {mname}")
        r0 = f(num, pos); E0, G0 = r0["energy"], r0["gradient"]
        line("energy (Eh)", abs(E0), np.inf, f"E = {E0:.10e}")
        line("translation invariance (shift 1e3 bohr)", abs(f(num, pos + [700, -400, 500])["energy"] - E0) / abs(E0), 1e-12)
        Rm = Rotation.random(random_state=3).as_matrix()
        r1 = f(num, pos @ Rm.T)
        line("rotation: energy", abs(r1["energy"] - E0) / abs(E0), 1e-12)
        line("rotation: gradient covariance", rel(r1["gradient"], G0 @ Rm.T), 1e-10)
        perm = rng.permutation(len(num))
        r2 = f(num[perm], pos[perm])
        line("atom permutation: energy", abs(r2["energy"] - E0) / abs(E0), 1e-12)
        line("atom permutation: gradient", rel(r2["gradient"], G0[perm]), 1e-10)
        line("sum of forces (should be 0)", np.max(np.abs(G0.sum(0))) / np.max(np.abs(G0)), 1e-12)
        line("torque (should be 0)", np.max(np.abs(np.cross(pos, G0).sum(0))) / np.max(np.abs(np.cross(pos, G0))), 1e-10)
        best = min(rel(G0, fd_grad(f, num, pos, None, h)) for h in (1e-3, 1e-4))
        line("gradient vs finite differences (best of h=1e-3,1e-4)", best, 1e-7)
    # dissociation: fragments far apart -> additivity
    if sysname == "benzene dimer":
        for mname, f in METHODS.items():
            half = len(num) // 2
            far = pos.copy(); far[half:] += [0, 0, 400.0]
            eA = f(num[:half], pos[:half])["energy"]; eB = f(num[half:], pos[half:])["energy"]
            line(f"{mname}: dissociation E(A..B @400 bohr) - E(A) - E(B)", abs(f(num, far)["energy"] - eA - eB) / abs(eA), 1e-10)

# ---------------- periodic oracles ----------------
for sysname, (num, pos, lat) in (("Si diamond", si_diamond()), ("NaCl rocksalt", nacl()), ("triclinic mixed", triclinic_mixed())):
    for mname, f in METHODS.items():
        print(f"--- {sysname} (periodic) / {mname}")
        r0 = f(num, pos, lat); E0, G0, V0 = r0["energy"], r0["gradient"], r0["virial"]
        line("energy per cell (Eh)", abs(E0), np.inf, f"E = {E0:.10e}")
        # translation of all atoms, and wrapping an atom by a lattice vector
        line("translation of all atoms", abs(f(num, pos + [1.1, -0.7, 0.3], lat)["energy"] - E0) / abs(E0), 1e-12)
        p2 = pos.copy(); p2[0] += lat[0] - 2 * lat[2]
        line("atom shifted by a lattice vector", abs(f(num, p2, lat)["energy"] - E0) / abs(E0), 1e-12)
        # equivalent cell choice: a' = a + b, b' = b, c' = c - a (unimodular)
        T = np.array([[1, 1, 0], [0, 1, 0], [-1, 0, 1]])
        line("equivalent cell (unimodular transform)", abs(f(num, pos, T @ lat)["energy"] - E0) / abs(E0), 1e-10)
        # supercell 2x2x2 -> 8 E
        shifts = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)]) @ lat
        nsc = np.tile(num, 8); psc = (pos[None] + shifts[:, None]).reshape(-1, 3)
        rsc = f(nsc, psc, 2 * lat)
        line("supercell 2x2x2: E_sc / 8 vs E_cell", abs(rsc["energy"] / 8 - E0) / abs(E0), 1e-9)
        Rm = Rotation.random(random_state=9).as_matrix()
        line("rotation of cell + atoms", abs(f(num, pos @ Rm.T, lat @ Rm.T)["energy"] - E0) / abs(E0), 1e-12)
        line("sum of forces (should be 0)", np.max(np.abs(G0.sum(0))) / max(np.max(np.abs(G0)), 1e-300), 1e-10)
        if len(num) > 2 or sysname == "triclinic mixed":
            best = min(rel(G0, fd_grad(f, num, pos, lat, h)) for h in (1e-3, 1e-4))
            line("gradient vs finite differences", best, 1e-7)
        S = min((fd_strain(f, num, pos, lat, h) for h in (1e-3, 1e-4)), key=lambda s: rel(V0, s))
        cands = {"virial = dE/deps": V0, "virial = -dE/deps": -V0, "virial^T = dE/deps": V0.T}
        best_name = min(cands, key=lambda k: rel(cands[k], S))
        line(f"virial vs strain derivative ({best_name})", rel(cands[best_name], S), 1e-7)
