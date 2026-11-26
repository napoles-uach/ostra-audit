# Open candidate (explorer 1, 2026-11-19): triclinic FD deviations in periodic D3. Check analytic gradient and
# virial against central finite differences on equivalent cells (primitive and sheared a' = a + 2b).
import sys, numpy as np
sys.path.insert(0, "/home/claude/impact2")
from cells import MATS, variants
from dftd3.interface import RationalDampingParam, DispersionModel
AU = 1 / 0.52917721067
P = RationalDampingParam(method="pbe0")
def calc(num, xyz, lat, grad=True):
    return DispersionModel(num, xyz, lattice=lat, periodic=np.array([True] * 3)).get_dispersion(P, grad=grad)
for name in ("NaCl (rocksalt)", "MgO (rocksalt)", "ZnO (wurtzite)"):
    for vname, at in variants(MATS[name][0]).items():
        if vname not in ("primitive", "a' = a + 2b", "conventional cubic", "hexagonal, gamma = 60"): continue
        num, xyz, lat = at.numbers, at.positions * AU, np.array(at.cell) * AU
        r = calc(num, xyz, lat); g, sig = r["gradient"], r["virial"]
        h = 1e-4; gfd = np.zeros_like(g)
        for i in range(len(num)):
            for k in range(3):
                d = np.zeros_like(xyz); d[i, k] = h
                gfd[i, k] = (calc(num, xyz + d, lat, False)["energy"] - calc(num, xyz - d, lat, False)["energy"]) / (2 * h)
        sfd = np.zeros((3, 3)); e = 1e-5
        for a in range(3):
            for b in range(3):
                eps = np.eye(3); eps[a, b] += e; ep = calc(num, xyz @ eps.T, lat @ eps.T, False)["energy"]
                eps = np.eye(3); eps[a, b] -= e; em = calc(num, xyz @ eps.T, lat @ eps.T, False)["energy"]
                sfd[a, b] = (ep - em) / (2 * e)
        print(f"{name:16s} {vname:22s} max|g-gFD| = {np.abs(g - gfd).max():.1e} (|g| {np.abs(g).max():.1e})   "
              f"max|sigma-sigmaFD| = {np.abs(sig - sfd).max():.1e} (|sigma| {np.abs(sig).max():.1e})", flush=True)
