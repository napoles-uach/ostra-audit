"""libcint oracles with closed-form references: single s-type Gaussians."""
import numpy as np
from math import erf, sqrt, pi
from pyscf import gto

def mol_s(centers, exps):
    """one uncontracted normalized s function per (center, exponent), dummy ghost atoms"""
    atom = [["ghost-H", c] for c in centers]
    basis = {}
    atoms = []
    for i, (c, a) in enumerate(zip(centers, exps)):
        lab = f"X{i}"
        atoms.append([lab, c])
        basis[lab] = [[0, [a, 1.0]]]
    m = gto.M(atom=atoms, basis=basis, unit="bohr", spin=None, charge=0, verbose=0)
    return m

def ghost(atoms_basis):
    pass

worst = {}
def rec(tag, got, ref):
    err = abs(got - ref) / max(abs(ref), 1e-300)
    if err > worst.get(tag, (0,))[0]:
        worst[tag] = (err, got, ref)

# use gto.Mole with ghost atoms (charge 0) so nuclear attraction can be controlled
def build(spec):
    atoms, basis = [], {}
    for i, (c, a) in enumerate(spec):
        lab = f"GHOST-H{i}" if False else f"X-{i}"
        atoms.append([f"X{i}", c])
        basis[f"X{i}"] = gto.basis.parse(f"X S\n {a:.17g} 1.0")
    m = gto.Mole()
    m.atom = atoms; m.basis = basis; m.unit = "bohr"; m.verbose = 0
    m.charge = 0; m.spin = 0
    # ghost atoms: zero nuclear charge
    m.build(dump_input=False, parse_arg=False)
    return m

rng = np.random.default_rng(0)
exps = 10 ** np.linspace(-6, 9, 61)
# (1) self-overlap of a normalized s function = 1 ; kinetic <s|T|s> = 3a/2
for a in exps:
    m = build([((0, 0, 0), a)])
    S = m.intor("int1e_ovlp")[0, 0]
    T = m.intor("int1e_kin")[0, 0]
    rec("S_self", S, 1.0)
    rec("T_self", T, 1.5 * a)
# (2) overlap of two s at distance R: (2 sqrt(ab)/(a+b))^{3/2} exp(-ab R^2/(a+b))
for _ in range(4000):
    a, b = 10 ** rng.uniform(-5, 7, 2)
    R = 10 ** rng.uniform(-6, 3)
    mu = a * b / (a + b)
    if mu * R * R > 600:
        continue
    m = build([((0, 0, 0), a), ((0, 0, R), b)])
    S = m.intor("int1e_ovlp")[0, 1]
    ref = (2 * sqrt(a * b) / (a + b)) ** 1.5 * np.exp(-mu * R * R)
    if ref > 1e-250:
        rec("S_pair", S, ref)
# (3) Coulomb of two normalized s densities |phi_a|^2 |phi_b|^2 at distance R:
#     (aa|bb) = erf(sqrt(g) R)/R, g = (2a)(2b)/(2a+2b); R->0: 2 sqrt(g/pi)
for _ in range(4000):
    a, b = 10 ** rng.uniform(-5, 7, 2)
    R = 0.0 if rng.random() < 0.05 else 10 ** rng.uniform(-8, 4)
    g = (2 * a) * (2 * b) / (2 * a + 2 * b)
    m = build([((0, 0, 0), a), ((0, 0, R if R > 0 else 0.0), b)]) if R > 0 else build([((0, 0, 0), a), ((0, 0, 1e-300), b)])
    eri = m.intor("int2e", aosym="s1").reshape(2, 2, 2, 2)[0, 0, 1, 1]
    x = sqrt(g) * R
    ref = 2 * sqrt(g / pi) * (1 - x * x / 3) if x < 1e-6 else erf(x) / R
    rec("ERI_ssss", eri, ref)
# (4) nuclear attraction of normalized s at distance R from a unit point charge:
#     <s|1/|r-C||s> = erf(sqrt(2a) R)/R
from pyscf.gto import moleintor
for _ in range(3000):
    a = 10 ** rng.uniform(-5, 8)
    R = 10 ** rng.uniform(-8, 4)
    m = build([((0, 0, 0), a)])
    with m.with_rinv_origin((0, 0, R)):
        v = m.intor("int1e_rinv")[0, 0]
    x = sqrt(2 * a) * R
    ref = 2 * sqrt(2 * a / pi) * (1 - x * x / 3) if x < 1e-6 else erf(x) / R
    rec("rinv", v, ref)

for k, (e, g, r) in sorted(worst.items()):
    print(f"{k:10s} max rel err {e:9.2e}   got {g:.17g}  ref {r:.17g}")
