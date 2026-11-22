"""D4 (dftd4 Python API) for M-O dimers: EEQ charges for Z = 112, 118 (beyond multicharge max_elem = 103)."""
import numpy as np
from dftd4.interface import DampingParam, DispersionModel
for Z, lab in ((80, "Hg"), (112, "Cn"), (118, "Og")):
    m = DispersionModel(np.array([Z, 8]), np.array([[0, 0, 0], [0, 0, 3.6]]))
    q = m.get_properties()["partial charges"]
    e = m.get_dispersion(DampingParam(method="pbe"), grad=False)["energy"]
    print(lab, "charges", np.round(q, 4), "E_disp %.6e" % e)
