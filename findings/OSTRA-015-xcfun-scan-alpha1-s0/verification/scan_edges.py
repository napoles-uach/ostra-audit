"""XCFun vs libxc SCAN at alpha = 1 and s = 0, and FD check of v_sigma at s = 0 (PySCF backends)."""
import numpy as np, warnings; warnings.filterwarnings("ignore")
from pyscf.dft import xcfun, libxc
rho = 0.1; tu = 0.3 * (3*np.pi**2)**(2/3) * rho**(5/3)
for a in (1.0, 1.2):
    R = np.array([[rho], [0.0], [0.0], [0.0], [0.0], [a * tu]])
    for code, lx in (("SCANX", "MGGA_X_SCAN"), ("SCANC", "MGGA_C_SCAN")):
        e1, v1 = xcfun.eval_xc(code, R, spin=0, deriv=1)[:2]; e2, v2 = libxc.eval_xc(lx, R, spin=0, deriv=1)[:2]
        print("s=0 alpha=%.1f %s exc xcfun=%s libxc=%s | vsigma xcfun=%s libxc=%s" % (a, code, e1, e2, v1[1], v2[1]))
tau = 1.2 * tu
for mod, c in ((xcfun, "SCANX"), (libxc, "MGGA_X_SCAN")):
    ex = lambda sg: (mod.eval_xc(c, np.array([[rho], [np.sqrt(sg/3)]*1, [np.sqrt(sg/3)], [np.sqrt(sg/3)], [0.0], [tau]], dtype=float).reshape(6, 1), spin=0, deriv=0)[0] * rho)[0]
    print(mod.__name__, "FD vsigma at s=0:", [(ex(d) - ex(0.0)) / d for d in (1e-6, 1e-8)])
