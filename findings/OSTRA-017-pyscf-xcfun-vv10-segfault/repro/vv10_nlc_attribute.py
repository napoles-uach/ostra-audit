# The supported route (mf.nlc) works with XCFun; only xc strings that name a VV10 functional crash.
import sys
from pyscf import gto, dft
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
xc, nlc, be = sys.argv[1], sys.argv[2], sys.argv[3]
mf = dft.RKS(mol); mf.xc = xc
if nlc != "-": mf.nlc = nlc
if be == "xcfun": mf._numint.libxc = dft.xcfun
print(xc, nlc, be, "E=%.8f" % mf.kernel())
