"""libcint oracles without references: symmetry, invariance, consistency."""
import numpy as np
from pyscf import gto
from scipy.spatial.transform import Rotation
from scipy.linalg import eigh as geigh

rng = np.random.default_rng(1)

def rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)

def mol_with_l(lmax, exps, geom, cart=False):
    """atoms X (ghost) at geom, each with one shell per l<=lmax per exponent"""
    atoms = [[f"X{i}", c] for i, c in enumerate(geom)]
    basis = {}
    for i in range(len(geom)):
        txt = "".join(f"X {'SPDFGHIKLM'[l]}\n {e:.17g} 1.0\n" for l in range(lmax + 1) for e in exps)
        basis[f"X{i}"] = gto.basis.parse(txt)
    m = gto.Mole(atom=atoms, basis=basis, unit="bohr", verbose=0, cart=cart)
    m.build(dump_input=False, parse_arg=False)
    return m

results = []
def report(tag, val, note=""):
    results.append((tag, val, note))
    print(f"{tag:58s} {val:9.2e} {note}", flush=True)

geom3 = rng.normal(size=(3, 3)) * 1.5
for lmax, exps in ((4, [0.9]), (6, [1.3]), (3, [3e4, 0.05])):
    lab = f"l<={lmax} exps={exps}"
    m = mol_with_l(lmax, exps, geom3)
    S, T = m.intor("int1e_ovlp"), m.intor("int1e_kin")
    # (a) hermiticity
    report(f"[{lab}] S, T symmetric", max(rel(S, S.T), rel(T, T.T)))
    # (b) ERI 8-fold symmetry (small system only)
    if m.nao <= 80:
        eri = m.intor("int2e", aosym="s1")
        n = m.nao; e4 = eri.reshape(n, n, n, n)
        sym = max(rel(e4, e4.transpose(1, 0, 2, 3)), rel(e4, e4.transpose(0, 1, 3, 2)), rel(e4, e4.transpose(2, 3, 0, 1)))
        report(f"[{lab}] ERI 8-fold permutational symmetry (nao={n})", sym)
        # (c) ERI as a Gram matrix (ij|kl) over pairs must be positive semidefinite
        ev = np.linalg.eigvalsh(eri.reshape(n * n, n * n))
        report(f"[{lab}] ERI pair-matrix min eig / max eig", ev[0] / ev[-1], "(>= -eps expected)")
    # (d) translation invariance, small and huge shift
    for shift in (3.7, 1e3, 1e5):
        d = rng.normal(size=3); d *= shift / np.linalg.norm(d)
        m2 = mol_with_l(lmax, exps, geom3 + d)
        report(f"[{lab}] translation by {shift:g} bohr: S,T", max(rel(m2.intor("int1e_ovlp"), S), rel(m2.intor("int1e_kin"), T)))
    # (e) rotation invariance of basis-independent invariants
    Rm = Rotation.random(random_state=3).as_matrix()
    m3 = mol_with_l(lmax, exps, geom3 @ Rm.T)
    S3, T3 = m3.intor("int1e_ovlp"), m3.intor("int1e_kin")
    report(f"[{lab}] rotation: eigvals of S", rel(np.linalg.eigvalsh(S3), np.linalg.eigvalsh(S)))
    report(f"[{lab}] rotation: eigvals of S^-1/2 T S^-1/2", rel(geigh(T3, S3, eigvals_only=True), geigh(T, S, eigvals_only=True)))
    # (f) cartesian vs spherical consistency
    mc = mol_with_l(lmax, exps, geom3, cart=True)
    c2s = mc.cart2sph_coeff()
    Sc = mc.intor("int1e_ovlp")
    report(f"[{lab}] cart->sph transform of S vs sph S", rel(c2s.T @ Sc @ c2s, S))
    # (g) erf + erfc = full Coulomb (nuclear attraction and ERI)
    if m.nao <= 80:
        w = 0.37
        with m.with_range_coulomb(w):
            lr = m.intor("int2e", aosym="s1")
        with m.with_range_coulomb(-w):
            sr = m.intor("int2e", aosym="s1")
        report(f"[{lab}] ERI: erf(w) + erfc(w) = full", rel(lr + sr, eri))

# (h) dipole origin shift: <i|r-R|j> = <i|r|j> - R S
m = mol_with_l(4, [0.8], geom3)
R0 = np.array([0.3, -1.2, 2.5])
r0 = m.intor("int1e_r")
with m.with_common_orig(R0):
    rR = m.intor("int1e_r")
report("dipole origin shift <r-R> = <r> - R S (l<=4)", rel(rR, r0 - R0[:, None, None] * m.intor("int1e_ovlp")[None]))
# (i) kinetic = -1/2 <i|nabla^2|j> ; <i|nabla^2|j> via int1e_ipipovlp trace? use p^2 identity: T = 1/2 sum_k <d_k i|d_k j>
pp = m.intor("int1e_ipovlpip")  # <nabla i| . |nabla j> components (9, n, n)
Tid = 0.5 * (pp[0] + pp[4] + pp[8])
report("T = 1/2 sum_k <d_k i|d_k j> (l<=4)", rel(Tid, m.intor("int1e_kin")))
