import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from dftd3.interface import GeometricCounterpoise
for z in (94, 95, 100, 103):
    r = GeometricCounterpoise(np.array([z, 8]), np.array([[0,0,0],[0,0,4.16]]), method="r2scan3c").get_counterpoise(grad=True)
    print(z, r["energy"], np.abs(r["gradient"]).max())
