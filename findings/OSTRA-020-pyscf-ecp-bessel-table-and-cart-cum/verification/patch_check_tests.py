import sys, numpy as np
from pyscf import gto
print("pyscf from", __import__("pyscf").__file__.split("/b020/")[-1][:20])
def syn(lb, lc, op):
    ecp={"Na":[10,[[-1,[[],[[1.0,2.0]]]]]+[[l,[[],[],[[1.5,3.0]]]] for l in range(lc+1)]]}
    mol=gto.M(atom="Na 0 0 0; H 0 0 1.5",basis={"Na":[[lb,(1.2,1.0)]],"H":[[0,(1.0,1.0)]]},ecp=ecp,verbose=0)
    return mol.intor(op)
case=sys.argv[1]
if case=="agh":
    mol=gto.M(atom="Ag 0 0 0; H 0 0 1.62",basis={"Ag":"cc-pv5z-pp","H":"sto-3g"},ecp={"Ag":"def2-svp"},verbose=0)
    for op in ("ECPscalar_ipnuc","ECPscalar_iprinv"): v=mol.intor(op); print(op,"ok",v.shape)
elif case=="syn":
    v=syn(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]); print("ok",np.abs(v).max()); np.save(sys.argv[5],v)
