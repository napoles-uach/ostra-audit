"""xc_func_set_fhc_enforcement(p, 0) at a point with tau < tau_W: effective for pure TPSS-x, no-op for TPSSh/TPSS0."""
import numpy as np, ctypes, lxc
L = lxc.L; L.xc_func_set_fhc_enforcement.argtypes = [ctypes.c_void_p, ctypes.c_int]
r = np.array([[0.1]]); sig = np.array([[0.05]]); tau = 0.5 * sig / (8 * r)
for fid in (202, 457, 396):
    out = []
    for on in (1, 0):
        f = lxc.Func(fid, 1); L.xc_func_set_fhc_enforcement(ctypes.c_void_p(f.p), on)
        out.append(f.exc_vxc({"rho": r, "sigma": sig, "tau": tau, "lapl": np.zeros_like(r)})["zk"][0]); f.close()
    print("%-22s zk(FHC on)=%.12f zk(FHC off)=%.12f changed: %s" % (lxc.name(fid), out[0], out[1], out[0] != out[1]))
