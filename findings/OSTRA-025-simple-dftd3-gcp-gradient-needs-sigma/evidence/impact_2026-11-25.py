# Explorer 1, 2026-11-25. dftd3 1.7.0 Python package. Same C API call with and without sigma.
import numpy as np
from dftd3 import library
from dftd3.interface import GeometricCounterpoise
num = np.array([8, 1, 1, 8, 1, 1])  # water dimer, Angstrom -> bohr
xyz = np.array([[-1.551007, -0.114520, 0.0], [-1.934259, 0.762503, 0.0], [-0.599677, 0.040712, 0.0],
                [1.350625, 0.111469, 0.0], [1.680398, -0.373741, -0.758561], [1.680398, -0.373741, 0.758561]]) / 0.52917721
g = GeometricCounterpoise(num, xyz, method="hf3c")
full = g.get_counterpoise(grad=True)
e = np.array(0.0); grad = np.zeros((6, 3))
library.get_counterpoise(g._mol, g._gcp, e, grad, None)   # dftd3_get_counterpoise(..., gradient, sigma = NULL)
print("max |grad| with sigma  :", np.abs(full["gradient"]).max())   # 0.0215 Eh/bohr
print("max |grad| sigma = NULL:", np.abs(grad).max())               # 0.0
