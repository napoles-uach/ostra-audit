import sys
from pyscf import gto
b,e,op=sys.argv[1],sys.argv[2],sys.argv[3]
mol=gto.M(atom="Ag 0 0 0; H 0 0 1.62",basis={"Ag":b,"H":"sto-3g"},ecp={"Ag":e},verbose=0)
ls=sorted(set(mol._bas[:,1])); lc=sorted(set(mol._ecpbas[:,1]))
v=mol.intor(op); print("basis l:",ls,"ecp channels:",lc,op,"ok",v.shape)
