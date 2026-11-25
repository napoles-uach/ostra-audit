import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import dftd3; from dftd3 import library
from dftd3.interface import GeometricCounterpoise
assert dftd3.__file__.startswith(sys.argv[1]), dftd3.__file__
B = 1/0.52917721
out = {}
def gcp(num, xyz, method, basis=None):
    kw = dict(method=method) if basis is None else dict(method=method, basis=basis)
    return GeometricCounterpoise(np.array(num), np.array(xyz)*B, **kw)
# 019
for num, m, b in (([36,8],"hf","minis"),([54,8],"hf","minis"),([87,8],"r2scan3c",None),([88,8],"r2scan3c",None),([89,8],"r2scan3c",None),([19,8],"r2scan3c",None),([20,8],"r2scan3c",None),([90,8],"r2scan3c",None)):
    r = gcp(num, [[0,0,0],[0,0,2.2]], m, b).get_counterpoise(grad=True)
    out[f"019 {num} {m}"] = [float(r["energy"]), float(np.abs(r["gradient"]).max())]
# 025
w = [[-1.551007,-0.114520,0],[-1.934259,0.762503,0],[-0.599677,0.040712,0],[1.350625,0.111469,0],[1.680398,-0.373741,-0.758561],[1.680398,-0.373741,0.758561]]
g = gcp([8,1,1,8,1,1], w, "hf3c"); full = g.get_counterpoise(grad=True)
e = np.array(0.0); gr = np.zeros((6,3)); library.get_counterpoise(g._mol, g._gcp, e, gr, None)
out["025 max|full grad|"] = float(np.abs(full["gradient"]).max()); out["025 max|grad, sigma NULL - full|"] = float(np.abs(gr-full["gradient"]).max())
# regression: common methods on a small organic molecule + water dimer
mol = ([6,8,1,1], [[0,0,0],[0,0,1.205],[0,0.943,-0.587],[0,-0.943,-0.587]])
for m, b in (("hf3c",None),("pbeh3c",None),("r2scan3c",None),("b973c",None),("hf","def2svp"),("b3lyp","631gd")):
    for name, (num, xyz) in (("H2CO", mol), ("H2O2", ([8,1,1,8,1,1], w))):
        try:
            r = gcp(num, xyz, m, b).get_counterpoise(grad=True)
            out[f"reg {name} {m} {b}"] = [float(r["energy"])] + r["gradient"].ravel().tolist() + r["virial"].ravel().tolist()
        except Exception as ex: out[f"reg {name} {m} {b}"] = str(ex)
json.dump(out, open(sys.argv[2], "w"))
