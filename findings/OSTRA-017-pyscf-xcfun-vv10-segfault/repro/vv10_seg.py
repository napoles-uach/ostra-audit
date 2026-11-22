import sys
from pyscf import gto, dft
mol=gto.M(atom="H 0 0 0; H 0 0 0.74",basis="sto-3g",verbose=0)
mf=dft.RKS(mol); mf.xc=sys.argv[1]
if sys.argv[2]=="xcfun": mf._numint.libxc=dft.xcfun
print(sys.argv[1],sys.argv[2],"E=%.8f"%mf.kernel())
