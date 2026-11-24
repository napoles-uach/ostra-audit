# Non-self-consistent XC energy of a stored density with the XCFun library loaded in this process (LD_PRELOAD).
import sys, json, numpy
from pyscf import gto, dft
from mols import MOLS, ATOMS
dmfile, xc, tag = sys.argv[1:4]
sp = dmfile.split("/")[-1].split("__")[0]
if sp in MOLS: atom, spin = MOLS[sp]
else: atom, spin = f"{sp} 0 0 0", ATOMS[sp]
mol = gto.M(atom=atom, spin=spin, basis=sys.argv[4], verbose=0)
dm = numpy.load(dmfile)
g = dft.gen_grid.Grids(mol); g.level = 5; g.build()
ni = dft.numint.NumInt(); ni.libxc = dft.xcfun
if spin: n, exc, v = ni.nr_uks(mol, g, xc, dm)
else:    n, exc, v = ni.nr_rks(mol, g, xc, dm)
print(json.dumps({"sp": sp, "xc": xc, "tag": tag, "Exc": exc}))
