import sys
from pyscf import gto
lb,lc=int(sys.argv[1]),int(sys.argv[2])
ecp={"Na":[10,[[-1,[[],[[1.0,2.0]]]]]+[[l,[[],[],[[1.5,3.0]]]] for l in range(lc+1)]]}
mol=gto.M(atom="Na 0 0 0; H 0 0 1.5",basis={"Na":[[lb,(1.2,1.0)]],"H":[[0,(1.0,1.0)]]},ecp=ecp,spin=None,charge=0,verbose=0)
v=mol.intor(sys.argv[3] if len(sys.argv)>3 else "ECPscalar_ipnuc"); print("lb=%d lc=%d ok, max|v|=%.3e"%(lb,lc,abs(v).max()))
