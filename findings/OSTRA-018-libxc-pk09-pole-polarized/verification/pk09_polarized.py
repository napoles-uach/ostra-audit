"""lda_c_pk09 around r_s = 22.739412132480812, unpolarized-equivalent (zeta = 0) and zeta = 0.4. Needs probes/libxc/lxc.py."""
import numpy as np, lxc
fid = [i for i in lxc.all_ids() if lxc.name(i) == "lda_c_pk09"][0]
f = lxc.Func(fid, 2); rs0 = 22.739412132480812
rs = rs0 + np.array([-1e-4, -1e-6, 0, 1e-6, 1e-4]); rho = 3 / (4 * np.pi * rs**3)
for z in (0.0, 0.4):
    o = f.exc_vxc({"rho": np.stack([rho * (1 + z) / 2, rho * (1 - z) / 2], 1)})
    print("zeta=%.1f zk=%s vrho_up=%s" % (z, np.round(o["zk"], 6), o["vrho"][:, 0]))
rr = np.linspace(rs0 - 0.05, rs0 + 0.05, 20001); rho = 3 / (4 * np.pi * rr**3)
o = f.exc_vxc({"rho": np.stack([rho * 0.7, rho * 0.3], 1)}); bad = ~np.isfinite(o["zk"])
print("zeta=0.4: non-finite points in +-0.05:", bad.sum(), (rr[bad].min() - rs0, rr[bad].max() - rs0) if bad.any() else "")
