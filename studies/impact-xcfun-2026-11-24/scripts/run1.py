# One species, one functional, one backend. Prints a JSON line.
import sys, json, time
from pyscf import gto, dft
from mols import MOLS, ATOMS
sp, xc, be, basis = sys.argv[1:5]
if sp in MOLS: atom, spin = MOLS[sp]
else: atom, spin = f"{sp} 0 0 0", ATOMS[sp]
mol = gto.M(atom=atom, spin=spin, basis=basis, verbose=0)
mf = dft.UKS(mol) if spin else dft.RKS(mol)
mf.xc = xc; mf.grids.level = 5; mf.conv_tol = 1e-10; mf.max_cycle = 200
if be != "libxc": mf._numint.libxc = dft.xcfun
t = time.time(); out = {"sp": sp, "xc": xc, "be": be}
try:
    e = mf.kernel()
    if not mf.converged:
        mf = mf.newton(); e = mf.kernel()
    out.update(E=e, conv=bool(mf.converged))
    import os, numpy
    if os.environ.get("SAVEDM"):
        numpy.save(f"dm/{sp}__{xc.replace(',','_')}__{be}.npy", numpy.asarray(mf.make_rdm1()))
except Exception as ex:
    out.update(E=None, err=f"{type(ex).__name__}: {str(ex)[:100]}")
out["t"] = round(time.time() - t, 1)
print(json.dumps(out))
