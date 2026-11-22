"""Oracles for libcint ECP integrals (scalar part), via PySCF."""
import numpy as np
from pyscf import gto
from scipy.linalg import eigh
from scipy.spatial.transform import Rotation

L = "SPDFGHIK"
rng = np.random.default_rng(42)

def rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)

def custom_ecp(lmax_proj, nelec=10):
    """semi-local ECP: local term + projectors l = 0..lmax_proj, two Gaussians each, r^0 and r^2 powers"""
    blocks = [[-1, [[], [], [[3.1, -2.5], [0.9, 0.4]], [], [], [], []]]]
    for l in range(lmax_proj + 1):
        blocks.append([l, [[], [], [[5.0 + l, 30.0 / (l + 1)], [1.3 + 0.2 * l, -4.0]], [], [], [], []]])
    return [nelec, blocks]

def mol_custom(geom, lbas, ecp_lmax, cart=False, exps=(0.4, 2.5)):
    """atom 0 carries the ECP (as 'Kr', Z=36 -> 26 electrons left), other atoms are ghosts with basis"""
    atoms = [["Kr", geom[0]]] + [[f"X{i}", c] for i, c in enumerate(geom[1:], 1)]
    btxt = "".join(f"X {L[l]}\n {e:.17g} 1.0\n" for l in range(lbas + 1) for e in exps)
    basis = {"Kr": gto.basis.parse(btxt.replace("X ", "Kr ")), **{f"X{i}": gto.basis.parse(btxt) for i in range(1, len(geom))}}
    m = gto.Mole(atom=atoms, basis=basis, ecp={"Kr": custom_ecp(ecp_lmax)}, unit="bohr", verbose=0, cart=cart, spin=None)
    m.spin = m.nelectron % 2
    m.build(dump_input=False, parse_arg=False)
    return m

def ecp_deriv(m):
    """analytic d V_ecp / d R_A for every atom A, exactly as pyscf.grad.rhf assembles it"""
    ip = m.intor("ECPscalar_ipnuc", comp=3)
    out = []
    for A in range(m.natm):
        p0, p1 = m.aoslice_by_atom()[A, 2:]
        ecp_atoms = set(m._ecpbas[:, gto.ATOM_OF])
        if A in ecp_atoms:
            with m.with_rinv_at_nucleus(A):
                v = m.intor("ECPscalar_iprinv", comp=3)
        else:
            v = np.zeros((3, m.nao, m.nao))
        v[:, p0:p1] -= ip[:, p0:p1]
        out.append(v + v.transpose(0, 2, 1))
    return out

def fd(builder, geom, A, k, h=1e-4):
    vals = []
    for s in (-2, -1, 1, 2):
        g = geom.copy(); g[A, k] += s * h
        vals.append(builder(g).intor("ECPscalar"))
    return (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12 * h)

def run_case(tag, geom, lbas, ecp_lmax):
    b = lambda g, cart=False: mol_custom(g, lbas, ecp_lmax, cart)
    m = b(geom)
    V, S = m.intor("ECPscalar"), m.intor("int1e_ovlp")
    print(f"--- {tag}: basis l<={lbas}, ECP projectors l<={ecp_lmax}, nao={m.nao}")
    print(f"  hermiticity                              {rel(V, V.T):9.2e}")
    for shift in (2.3, 1e3):
        d = rng.normal(size=3); d *= shift / np.linalg.norm(d)
        print(f"  translation {shift:7g} bohr                 {rel(b(geom + d).intor('ECPscalar'), V):9.2e}")
    Rm = Rotation.random(random_state=7).as_matrix()
    m2 = b(geom @ Rm.T)
    ev1 = eigh(V, S, eigvals_only=True); ev2 = eigh(m2.intor("ECPscalar"), m2.intor("int1e_ovlp"), eigvals_only=True)
    print(f"  rotation: generalized eigenvalues        {rel(ev2, ev1):9.2e}")
    mc = b(geom, cart=True); c2s = mc.cart2sph_coeff()
    print(f"  cartesian -> spherical vs spherical      {rel(c2s.T @ mc.intor('ECPscalar') @ c2s, V):9.2e}")
    der = ecp_deriv(m)
    worst = 0; where = None
    for A in range(m.natm):
        for k in range(3):
            num = fd(b, geom, A, k)
            r = rel(der[A][k], num)
            if r > worst: worst, where = r, (A, "xyz"[k])
    print(f"  derivative vs finite differences         {worst:9.2e}  (worst at atom {where})")
    # translational sum rule: sum_A dV/dR_A = 0
    tot = sum(der)
    print(f"  sum over atoms of dV/dR_A (should be 0)  {np.max(np.abs(tot)) / np.max(np.abs(der[0])):9.2e}")
    return m

g2 = np.array([[0.1, -0.2, 0.05], [1.3, 0.9, 2.1]])
g3 = np.array([[0.0, 0.0, 0.0], [1.8, 0.4, -0.7], [-0.6, 2.2, 1.1]])
run_case("2 centres", g2, 2, 2)
run_case("2 centres, high l", g2, 4, 4)
run_case("3 centres (basis on two ghosts + ECP atom)", g3, 3, 3)
run_case("3 centres, high l", g3, 5, 4)

