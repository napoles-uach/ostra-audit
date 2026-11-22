"""libcint: derivative integrals vs finite differences; high-l ERIs vs rotation invariants."""
import numpy as np
from pyscf import gto
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(7)
L = "SPDFGHIKLM"

def mk(geom, shells, cart=False):
    """shells: list (per atom) of list of (l, exp)"""
    atoms = [[f"X{i}", c] for i, c in enumerate(geom)]
    basis = {f"X{i}": gto.basis.parse("".join(f"X {L[l]}\n {e:.17g} 1.0\n" for l, e in sh)) for i, sh in enumerate(shells)}
    m = gto.Mole(atom=atoms, basis=basis, unit="bohr", verbose=0, cart=cart)
    m.build(dump_input=False, parse_arg=False)
    return m

def rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)

def fd(fun, geom, atom, k, h=1e-4):
    vals = []
    for s in (-2, -1, 1, 2):
        g = geom.copy(); g[atom, k] += s * h
        vals.append(fun(g))
    return (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12 * h)

print("== first-derivative integrals vs 4-point finite differences ==")
for lmax in (2, 4, 6):
    geom = rng.normal(size=(2, 3)) * 1.2
    shells = [[(l, 0.7 + 0.3 * l) for l in range(lmax + 1)], [(l, 1.1) for l in range(lmax + 1)]]
    m = mk(geom, shells)
    n = m.nao; off = m.aoslice_by_atom()
    # d/dR_A <i|j> = -<d i|j> for i on A, -<i|d j> for j on A
    ip = m.intor("int1e_ipovlp")  # <nabla i|j>, (3,n,n)
    ipk = m.intor("int1e_ipkin")
    for name, ipm, f in (("overlap", ip, lambda g: mk(g, shells).intor("int1e_ovlp")),
                         ("kinetic", ipk, lambda g: mk(g, shells).intor("int1e_kin"))):
        worst = 0
        for A in range(2):
            p0, p1 = off[A, 2], off[A, 3]
            for k in range(3):
                ana = np.zeros((n, n))
                ana[p0:p1, :] -= ipm[k, p0:p1, :]
                ana[:, p0:p1] -= ipm[k, p0:p1, :].T
                worst = max(worst, rel(ana, fd(f, geom, A, k)))
        print(f"  l<={lmax}  d{name}/dR      max rel dev {worst:9.2e}")
    if n <= 60:
        ip1 = m.intor("int2e_ip1", aosym="s1").reshape(3, n, n, n, n)  # (nabla i j|k l)
        worst = 0
        for A in range(2):
            p0, p1 = off[A, 2], off[A, 3]
            for k in range(3):
                ana = np.zeros((n, n, n, n))
                g = ip1[k]
                ana[p0:p1] -= g[p0:p1]
                ana[:, p0:p1] -= g[p0:p1].transpose(1, 0, 2, 3)
                ana[:, :, p0:p1] -= g[p0:p1].transpose(2, 3, 0, 1)
                ana[:, :, :, p0:p1] -= g[p0:p1].transpose(2, 3, 1, 0)
                num = fd(lambda gg: mk(gg, shells).intor("int2e", aosym="s1").reshape(n, n, n, n), geom, A, k)
                worst = max(worst, rel(ana, num))
        print(f"  l<={lmax}  dERI/dR (nao={n})  max rel dev {worst:9.2e}")

print("== high-l ERIs: rotation invariants, several separations ==")
for lmax in (4, 5, 6):
    for R in (0.0, 0.5, 3.0, 20.0, 200.0):
        geom = np.array([[0, 0, 0], [0, 0, R if R > 0 else 1e-12]]) + np.array([0.2, -0.1, 0.3])
        shells = [[(lmax, 0.9)], [(lmax, 1.7)]]
        m = mk(geom, shells)
        Rm = Rotation.random(random_state=lmax * 10 + int(R)).as_matrix()
        m2 = mk(geom @ Rm.T, shells)
        n = m.nao
        e1 = m.intor("int2e", aosym="s1"); e2 = m2.intor("int2e", aosym="s1")
        f1, f2 = np.linalg.norm(e1), np.linalg.norm(e2)
        v1 = np.linalg.eigvalsh(e1.reshape(n * n, n * n)); v2 = np.linalg.eigvalsh(e2.reshape(n * n, n * n))
        print(f"  l={lmax} R={R:6.1f}  |ERI|_F rel diff {abs(f1 - f2) / f1:9.2e}   pair-matrix eigvals rel diff {rel(v2, v1):9.2e}   min eig/max {v1[0] / v1[-1]:9.2e}")
