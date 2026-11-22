"""EEQ charges (multicharge, used by dftd4) depend on the choice of an equivalent unit cell."""
import numpy as np
import dftd4
from dftd4.interface import DispersionModel, DampingParam

periodic = np.array([True, True, True])
a = 10.66  # NaCl rocksalt, bohr
prim = 0.5 * a * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
numbers = np.array([11, 17])
positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) @ prim

def run(lattice):
    model = DispersionModel(numbers, positions, lattice=lattice, periodic=periodic)
    q = model.get_properties()["partial charges"]
    e = model.get_dispersion(DampingParam(method="pbe0"), grad=False)["energy"]
    return q[0], e

print("dftd4", dftd4.__version__)
q_ref, e_ref = run(prim)
print(f"primitive cell        q(Na) = {q_ref:.8f}   E_disp = {e_ref:.10e} Eh")
for k in range(1, 6):
    # a' = a + k b : same lattice, same crystal, only the cell description changes
    lattice = np.array([[1, k, 0], [0, 1, 0], [0, 0, 1]]) @ prim
    q, e = run(lattice)
    print(f"a' = a + {k} b          q(Na) = {q:.8f}   E_disp = {e:.10e} Eh   "
          f"dq = {q - q_ref:+.2e}   dE/E = {(e - e_ref) / abs(e_ref):+.2e}")

# Verification: an LLL-reduced basis of the same lattice restores the reference charge.
def lll(B, delta=0.75):
    B = B.astype(float).copy(); n = len(B); k = 1
    def gs(B):
        Bs = np.zeros_like(B); mu = np.zeros((n, n))
        for i in range(n):
            Bs[i] = B[i]
            for j in range(i):
                mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j]); Bs[i] -= mu[i, j] * Bs[j]
        return Bs, mu
    Bs, mu = gs(B)
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k, j])
            if q:
                B[k] -= q * B[j]; Bs, mu = gs(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]):
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]; Bs, mu = gs(B); k = max(k - 1, 1)
    return B

for k in range(1, 6):
    lattice = np.array([[1, k, 0], [0, 1, 0], [0, 0, 1]]) @ prim
    q, e = run(lll(lattice))
    print(f"a' = a + {k} b, LLL-reduced first: q(Na) = {q:.8f}   dq = {q - q_ref:+.2e}")
