"""Energy-potential consistency of LC-wPBE (vrho vs 4-point FD of rho*zk at fixed sigma, best of 3 steps).
Usage: LIBXC_PATH=/path/libxc.so python3 lcwpbe_consistency.py   (needs probes/libxc/lxc.py on PYTHONPATH)"""
import numpy as np, lxc
for fid in (478, 428, 524, 101):
    f = lxc.Func(fid, 1)
    rho = np.logspace(-6, 3, 181); s = np.array([0.05, 0.2, 0.5, 1.5, 3, 6])
    R, S = [a.ravel() for a in np.meshgrid(rho, s, indexing="ij")]
    ev = lambda r, sg: f.exc_vxc({"rho": r[:, None], "sigma": sg[:, None]})
    sig = (2 * (3 * np.pi**2 * R)**(1/3) * R * S)**2; v = ev(R, sig)["vrho"][:, 0]
    best = np.full(len(R), np.inf)
    for h in (1e-3, 1e-4, 1e-5):
        d = h * R; E = lambda r: r * ev(r, sig)["zk"]
        fd = (-E(R + 2*d) + 8*E(R + d) - 8*E(R - d) + E(R - 2*d)) / (12*d)
        best = np.minimum(best, np.abs(fd - v) / (np.abs(v) + 1e-30))
    i = np.argmax(best); print(lxc.name(fid), "max %.1e at rho=%.3g s=%g" % (best[i], R[i], S[i]))
