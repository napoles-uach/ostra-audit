import sys, numpy
from pyscf import gto, dft
mol = gto.M(atom="H 0 0 0", spin=1, basis="def2-tzvp", verbose=0)
mf = dft.UKS(mol); mf.xc = "mgga_x_revtpss,mgga_c_revtpss"; mf.grids.level = 5; mf.kernel()
dm = mf.make_rdm1()
ni = dft.numint.NumInt(); ni.libxc = dft.xcfun
for xc in ("revtpssx,", ",revtpssc", "tpssx,", ",tpssc"):
    try: print(xc, ni.nr_uks(mol, mf.grids, xc, dm)[1])
    except Exception as e: print(xc, type(e).__name__, str(e)[:60])
