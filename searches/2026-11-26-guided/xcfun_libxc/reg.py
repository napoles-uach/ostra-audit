import numpy as np, sys
from pyscf.dft import xcfun
rng = np.random.default_rng(1); N = 2000
rho = np.vstack([rng.uniform(1e-4, 2, N), rng.normal(0, .5, (3, N)), rng.normal(0, 1, N), rng.uniform(0.01, 3, N)])
vals = []
for xc in ("PBE", "B3LYP", "TPSS", "SCAN", "M06L", "PW91", "LDA", "BP86", "TF"):
    t = xcfun.xc_type(xc); r = rho if t == "MGGA" else (rho[:4] if t == "GGA" else rho[0])
    exc, vxc, fxc, _ = xcfun.eval_xc(xc, r, deriv=2)
    for a in (exc,) + tuple(vxc) + tuple(fxc):
        if a is not None: vals.append(np.ravel(a))
np.save(sys.argv[1], np.concatenate(vals))
